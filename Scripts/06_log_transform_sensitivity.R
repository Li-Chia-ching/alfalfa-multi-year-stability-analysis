# -*- coding: utf-8 -*-
# ============================================================================
# 06_log_transform_sensitivity.R
# Sensitivity analysis: Log-transformation log(x + c) evaluation for FW/DW
# (and SummerPH for comparison) to assess impact on:
#   1. Residual diagnostics (normality, homoscedasticity, outliers)
#   2. Broad-sense heritability (H2) estimates
#   3. Genotypic BLUP rank consistency (Spearman correlation & top-5 overlap)
#   4. Model singularity/boundary estimates
# Transformation offset: c = 0.5 * minimum positive value for each trait.
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
})

# Define relative paths for reproducibility and data availability
out_dir <- file.path("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(out_dir, "plot_year_wide.csv")
d <- read.csv(input_file, stringsAsFactors = FALSE)
d <- d %>% mutate(
  year_f   = factor(calendar_year),
  plot     = factor(plot),
  genotype = factor(genotype_ID)
)

TRAITS  <- c("SummerPH", "FW", "DW")
COHORTS <- c("2001-established", "2002-established")

# Calculate minimum positive value per trait for log-transformation offset
offset_of <- function(v) {
  pos <- v[!is.na(v) & v > 0]
  if (length(pos) == 0) return(NA_real_)
  0.5 * min(pos)
}

offsets <- sapply(TRAITS, function(t) offset_of(d[[t]]))
cat("Calculated log transformation offsets (0.5 * min positive value):\n")
print(offsets)

# Apply log transformation: log(x + offset)
for (t in TRAITS) {
  d[[paste0("log_", t)]] <- log(d[[t]] + offsets[t])
}

fit_primary <- function(dd, trait) {
  f <- as.formula(paste0(trait, " ~ year_f + (1|genotype) + (1|genotype:year_f) + (1|plot)"))
  lmer(f, data = dd, REML = TRUE, control = lmerControl(calc.derivs = FALSE))
}

# Brown-Forsythe test (median-based Levene's test) for heteroscedasticity
hetero_test <- function(resid, fitted, ngrp = 4) {
  g <- cut(
    fitted,
    breaks = quantile(fitted, probs = seq(0, 1, length.out = ngrp + 1)),
    include.lowest = TRUE,
    labels = FALSE
  )
  g <- factor(g)
  if (nlevels(g) < 2) return(c(stat = NA_real_, p = NA_real_))
  
  z <- abs(resid - tapply(resid, g, median)[as.integer(g)])
  aov_fit <- tryCatch(anova(lm(z ~ g)), error = function(e) NULL)
  if (is.null(aov_fit)) return(c(stat = NA_real_, p = NA_real_))
  c(stat = unname(aov_fit[1, "F value"]), p = unname(aov_fit[1, "Pr(>F)"]))
}

results <- list()
for (coh in COHORTS) {
  for (t in TRAITS) {
    log_trait <- paste0("log_", t)
    dd_orig   <- d %>% filter(cohort == coh, !is.na(.data[[t]]))
    dd_log    <- d %>% filter(cohort == coh, !is.na(.data[[log_trait]]))

    m_orig <- tryCatch(fit_primary(dd_orig, t), error = function(e) NULL)
    m_log  <- tryCatch(fit_primary(dd_log,  log_trait), error = function(e) NULL)

    # Calculate broad-sense heritability (H2)
    h2 <- function(m, trait_col, dd) {
      vc  <- as.data.frame(VarCorr(m))
      VG  <- vc$vcov[vc$grp == "genotype"][1]
      VGY <- vc$vcov[vc$grp == "genotype:year_f"][1]
      VE  <- sigma(m)^2
      if (is.na(VG))  VG  <- 0
      if (is.na(VGY)) VGY <- 0
      ny  <- length(unique(dd$calendar_year))
      r   <- dd %>% group_by(genotype, year_f) %>% summarise(n = n(), .groups = "drop") %>% pull(n)
      r   <- 1 / mean(1 / r)
      VG / (VG + VGY / ny + VE / (ny * r))
    }

    # Standardized residual diagnostics
    diag_of <- function(m, dd) {
      r  <- residuals(m, type = "pearson", scaled = TRUE)
      f  <- fitted(m)
      n  <- length(r)
      sw <- tryCatch(shapiro.test(r)$p.value, error = function(e) NA_real_)
      ht <- hetero_test(r, f)
      c(n = n, shapiro_p = sw, bf_p = unname(ht["p"]), pct_z3 = 100 * mean(abs(r) > 3))
    }

    # BLUP rank consistency (Spearman correlation & top-5 overlap)
    rank_agree <- function(m_orig, m_log, dd_orig, dd_log, t, lt) {
      b1 <- ranef(m_orig)$genotype[[1]] + fixef(m_orig)["(Intercept)"]
      b2 <- ranef(m_log)$genotype[[1]] + fixef(m_log)["(Intercept)"]
      df <- data.frame(
        g     = rownames(ranef(m_orig)$genotype),
        orig  = b1,
        logv  = b2,
        stringsAsFactors = FALSE
      )
      # Back-transform log BLUPs to original scale for comparative evaluation
      df$log_back <- exp(df$logv) - offsets[t]
      rho  <- cor.test(df$orig, df$log_back, method = "spearman")$estimate
      top1 <- df$g[order(df$orig, decreasing = TRUE)[1:5]]
      top2 <- df$g[order(df$log_back, decreasing = TRUE)[1:5]]
      c(spearman = as.numeric(rho), top5_overlap = length(intersect(top1, top2)))
    }

    h2_orig   <- if (!is.null(m_orig)) h2(m_orig, t, dd_orig) else NA_real_
    h2_log    <- if (!is.null(m_log))  h2(m_log,  log_trait, dd_log) else NA_real_
    diag_orig <- if (!is.null(m_orig)) diag_of(m_orig, dd_orig) else rep(NA_real_, 4)
    diag_log  <- if (!is.null(m_log))  diag_of(m_log,  dd_log)  else rep(NA_real_, 4)
    agree     <- if (!is.null(m_orig) && !is.null(m_log)) {
      rank_agree(m_orig, m_log, dd_orig, dd_log, t, log_trait)
    } else {
      rep(NA_real_, 2)
    }

    results[[paste(coh, t)]] <- data.frame(
      cohort         = coh,
      trait          = t,
      log_offset     = offsets[t],
      H2_orig        = h2_orig,
      H2_log         = h2_log,
      H2_delta       = h2_log - h2_orig,
      singular_orig  = if (is.null(m_orig)) NA else isSingular(m_orig),
      singular_log   = if (is.null(m_log))  NA else isSingular(m_log),
      n_obs          = as.integer(diag_orig["n"]),
      shapiro_p_orig = diag_orig["shapiro_p"],
      shapiro_p_log  = diag_log["shapiro_p"],
      bf_p_orig      = diag_orig["bf_p"],
      bf_p_log       = diag_log["bf_p"],
      pct_z3_orig    = diag_orig["pct_z3"],
      pct_z3_log     = diag_log["pct_z3"],
      rank_spearman  = agree["spearman"],
      top5_overlap   = agree["top5_overlap"],
      stringsAsFactors = FALSE
    )
  }
}

res_tab <- do.call(rbind, results)
rownames(res_tab) <- NULL

output_file <- file.path(out_dir, "log_transform_sensitivity.csv")
write.csv(res_tab, output_file, row.names = FALSE)

cat("\n=== Log-Transformation Sensitivity Summary ===\n")
print(res_tab, digits = 3)
cat("\nResults written to:", output_file, "\n")
