# -*- coding: utf-8 -*-
# ============================================================================
# 04_ar1_h2_sensitivity.R
# Sensitivity analysis: Broad-sense heritability (H2) re-estimation under
# an AR(1) working residual covariance structure using nlme vs. primary lme4.
# ============================================================================

suppressPackageStartupMessages({
  library(nlme)
  library(dplyr)
})

# Define relative paths for reproducibility and data availability
out_dir <- file.path("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(out_dir, "plot_year_wide.csv")
h2_file <- file.path(out_dir, "heritability_final_by_cohort.csv")

d <- read.csv(input_file, stringsAsFactors = FALSE)
h2_lme4 <- read.csv(h2_file, stringsAsFactors = FALSE)

d <- d %>% mutate(
  year_f   = factor(calendar_year),
  year_num = calendar_year - 2001,
  plot     = factor(plot),
  genotype = factor(genotype_ID)
)

TRAITS <- c("SummerPH", "FW", "DW")
COHORTS <- c("2001-established", "2002-established")

# Fits lme model with blocked random effects and optional AR(1) autocorrelation
fit_one <- function(dd, trait, use_ar1) {
  f <- as.formula(paste0(trait, " ~ year_f"))
  rnd <- list(
    genotype = pdBlocked(list(pdIdent(~1), pdIdent(~year_f - 1))),
    plot     = ~1
  )
  ctrl <- lmeControl(
    opt = "optim",
    msMaxIter = 800,
    msMaxEval = 800,
    msVerbose = FALSE,
    returnObject = TRUE
  )
  args <- list(f, random = rnd, data = dd, method = "REML", control = ctrl)
  
  if (use_ar1) {
    args$correlation <- corAR1(form = ~ year_num | genotype / plot)
  }
  do.call(lme, args)
}

# Extracts variance components and calculates H2 from lme object
h2_from_lme <- function(m, dd, trait) {
  vc <- as.matrix(VarCorr(m))
  rn <- rownames(vc)
  int_pos <- which(rn == "(Intercept)")
  
  V_G   <- as.numeric(vc[int_pos[1], "Variance"])   # Genotype intercept
  V_P   <- as.numeric(vc[int_pos[2], "Variance"])   # Plot intercept
  gy_rn <- grep("^year_f", rn, value = TRUE)
  V_GY  <- if (length(gy_rn)) as.numeric(vc[gy_rn[1], "Variance"]) else NA_real_
  V_e   <- as.numeric(vc["Residual", "Variance"])
  
  if (is.na(V_G))  V_G  <- 0
  if (is.na(V_GY)) V_GY <- 0
  
  ny <- length(unique(dd$calendar_year))
  r  <- dd %>% group_by(genotype, year_f) %>% summarise(n = n(), .groups = "drop") %>% pull(n)
  r  <- 1 / mean(1 / r)
  H2 <- V_G / (V_G + V_GY / ny + V_e / (ny * r))
  
  list(V_G = V_G, V_GY = V_GY, V_plot = V_P, V_e = V_e, ny = ny, r = r, H2 = H2)
}

# Extracts AR(1) autocorrelation parameter (rho)
rho_of <- function(m) {
  if (is.null(m$modelStruct$corStruct)) return(NA_real_)
  tryCatch(
    as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE))[1],
    error = function(e) NA_real_
  )
}

res <- list()
for (coh in COHORTS) {
  for (trait in TRAITS) {
    dd <- d %>% filter(cohort == coh, !is.na(.data[[trait]]))
    if (nrow(dd) < 30) next
    cat("Fitting model:", coh, trait, "(n =", nrow(dd), ")\n")

    m_id  <- tryCatch(fit_one(dd, trait, FALSE), error = function(e) NULL)
    m_ar1 <- tryCatch(fit_one(dd, trait, TRUE),  error = function(e) NULL)

    s_id  <- if (!is.null(m_id))  h2_from_lme(m_id,  dd, trait) else NULL
    s_ar1 <- if (!is.null(m_ar1)) h2_from_lme(m_ar1, dd, trait) else NULL

    ref <- h2_lme4 %>% filter(cohort == coh, trait == !!trait) %>% pull(H2)
    ref <- if (length(ref)) ref[1] else NA_real_

    res[[paste(coh, trait)]] <- data.frame(
      cohort          = coh,
      trait           = trait,
      H2_lme4_primary = ref,
      H2_nlme_ID      = if (is.null(s_id))  NA_real_ else s_id$H2,
      H2_nlme_AR1     = if (is.null(s_ar1)) NA_real_ else s_ar1$H2,
      rho_AR1         = if (is.null(m_ar1)) NA_real_ else rho_of(m_ar1),
      AIC_ID          = if (is.null(m_id))  NA_real_ else AIC(m_id),
      AIC_AR1         = if (is.null(m_ar1)) NA_real_ else AIC(m_ar1),
      V_G_AR1         = if (is.null(s_ar1)) NA_real_ else s_ar1$V_G,
      V_GY_AR1        = if (is.null(s_ar1)) NA_real_ else s_ar1$V_GY,
      V_e_AR1         = if (is.null(s_ar1)) NA_real_ else s_ar1$V_e,
      stringsAsFactors = FALSE
    )
  }
}

res_tab <- do.call(rbind, res)
res_tab$delta_nlme_vs_lme4   <- res_tab$H2_nlme_ID  - res_tab$H2_lme4_primary
res_tab$delta_AR1_vs_ID      <- res_tab$H2_nlme_AR1 - res_tab$H2_nlme_ID
res_tab$delta_AR1_vs_primary <- res_tab$H2_nlme_AR1 - res_tab$H2_lme4_primary
rownames(res_tab) <- NULL

output_file <- file.path(out_dir, "h2_ar1_sensitivity.csv")
write.csv(res_tab, output_file, row.names = FALSE)

cat("\n=== H2 Sensitivity Summary (Primary lme4 vs. nlme ID vs. nlme AR1) ===\n")
print(res_tab, digits = 3)

cat(
  "\nMax |delta H2| (AR1 vs Primary lme4):",
  max(abs(res_tab$delta_AR1_vs_primary), na.rm = TRUE),
  "\n"
)
cat(
  "Max |delta H2| (AR1 vs nlme ID):",
  max(abs(res_tab$delta_AR1_vs_ID), na.rm = TRUE),
  "\n"
)
cat(
  "Max |delta H2| (nlme ID vs Primary lme4, structural equivalence check):",
  max(abs(res_tab$delta_nlme_vs_lme4), na.rm = TRUE),
  "\n"
)
cat("\nResults written to:", output_file, "\n")
