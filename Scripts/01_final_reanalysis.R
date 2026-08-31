# -*- coding: utf-8 -*-
# ============================================================================
# 01_final_reanalysis.R
# Statistical reanalysis for the alfalfa phenotype dataset.
# Input: plot_year_wide.csv (aggregated plot x year dataset)
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(nlme)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Define relative paths for reproducibility and data availability
input_file <- file.path("output", "plot_year_wide.csv")
out_dir <- file.path("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load aggregated dataset
d <- read.csv(input_file, stringsAsFactors = FALSE)
d <- d %>% mutate(
  year_f = factor(calendar_year),
  year_num = calendar_year - 2001,
  plot = factor(plot),
  genotype = factor(genotype_ID)
)

TRAITS <- c("SummerPH", "FW", "DW")
cat("Plot x Year rows:", nrow(d), "| Unique genotypes:", n_distinct(d$genotype), "\n")
cat("Cohort distribution:\n")
print(table(d$cohort))

# ============================================================================
# 1. Covariance Structure Comparison (nlme)
# ============================================================================
gfun <- function(m, nm) {
  if (is.null(m)) return(data.frame(model = nm, AIC = NA, BIC = NA, logLik = NA, rho = NA))
  rho <- tryCatch(as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE)),
                  error = function(e) NA_real_)
  if (length(rho) == 0) rho <- NA_real_
  data.frame(model = nm, AIC = AIC(m), BIC = BIC(m),
             logLik = as.numeric(logLik(m)), rho = rho[1])
}

cov_res <- list()
for (coh in c("2001-established", "2002-established")) {
  dd <- d %>% filter(cohort == coh)
  for (trait in TRAITS) {
    ddt <- dd %>% filter(!is.na(.data[[trait]]))
    if (nrow(ddt) < 20) next
    f <- as.formula(paste0(trait, " ~ year_f"))
    ctrl <- lmeControl(opt = "optim", msMaxIter = 500)
    m_id  <- tryCatch(lme(f, random = ~1|genotype/plot, data = ddt, method = "REML", control = ctrl), error = function(e) NULL)
    m_cs  <- tryCatch(lme(f, random = ~1|genotype/plot, data = ddt, method = "REML",
                          correlation = corCompSymm(form = ~1|genotype/plot), control = ctrl), error = function(e) NULL)
    m_ar1 <- tryCatch(lme(f, random = ~1|genotype/plot, data = ddt, method = "REML",
                          correlation = corAR1(form = ~year_num|genotype/plot), control = ctrl), error = function(e) NULL)
    cov_res[[paste(coh, trait)]] <- cbind(cohort = coh, trait = trait,
                                          rbind(gfun(m_id, "ID"), gfun(m_cs, "CS"), gfun(m_ar1, "AR1")))
  }
}
cov_tab <- do.call(rbind, cov_res)
rownames(cov_tab) <- NULL
write.csv(cov_tab, file.path(out_dir, "covariance_model_comparison_final.csv"), row.names = FALSE)

cat("\n=== Covariance Model Comparison ===\n")
print(cov_tab, digits = 2)

# ============================================================================
# 2. Main Effects & Variance Components Estimation (lme4)
# ============================================================================
varcomp <- list(); h2 <- list(); blups <- list()
for (coh in c("2001-established", "2002-established")) {
  dd <- d %>% filter(cohort == coh)
  for (trait in TRAITS) {
    ddt <- dd %>% filter(!is.na(.data[[trait]]))
    if (nrow(ddt) < 20) next
    f <- as.formula(paste0(trait, " ~ year_f + (1|genotype) + (1|genotype:year_f) + (1|plot)"))
    m <- tryCatch(lmer(f, data = ddt, REML = TRUE, control = lmerControl(calc.derivs = FALSE)), error = function(e) NULL)
    if (is.null(m)) { cat("  Model fitting failed:", coh, trait, "\n"); next }
    
    vc <- as.data.frame(VarCorr(m))
    VG  <- vc$vcov[vc$grp == "genotype"][1]
    VGY <- vc$vcov[vc$grp == "genotype:year_f"][1]
    VP  <- vc$vcov[vc$grp == "plot"][1]
    VE  <- sigma(m)^2
    if (is.na(VG)) VG <- 0; if (is.na(VGY)) VGY <- 0; if (is.na(VP)) VP <- 0
    
    varcomp[[paste(coh, trait)]] <- data.frame(cohort = coh, trait = trait,
                                               V_G = VG, V_GY = VGY, V_plot = VP, V_e = VE,
                                               singular = isSingular(m), stringsAsFactors = FALSE)
    ny <- length(unique(ddt$calendar_year))
    r  <- ddt %>% group_by(genotype, year_f) %>% summarise(n = n(), .groups = "drop") %>% pull(n)
    r  <- 1 / mean(1 / r)
    H2 <- VG / (VG + VGY / ny + VE / (ny * r))
    h2[[paste(coh, trait)]] <- data.frame(cohort = coh, trait = trait, H2 = H2, ny = ny, r_harmonic = r)
    
    rr <- ranef(m)$genotype
    blups[[paste(coh, trait)]] <- data.frame(
      cohort = coh, trait = trait, genotype = rownames(rr),
      BLUP = fixef(m)["(Intercept)"] + rr[[1]], stringsAsFactors = FALSE)
  }
}

varcomp_tab <- do.call(rbind, varcomp); rownames(varcomp_tab) <- NULL
h2_tab      <- do.call(rbind, h2);      rownames(h2_tab) <- NULL
blup_tab    <- do.call(rbind, blups);   rownames(blup_tab) <- NULL

write.csv(varcomp_tab, file.path(out_dir, "variance_components_final_by_cohort.csv"), row.names = FALSE)
write.csv(h2_tab,      file.path(out_dir, "heritability_final_by_cohort.csv"), row.names = FALSE)
write.csv(blup_tab,    file.path(out_dir, "BLUP_final_by_cohort.csv"), row.names = FALSE)

cat("\n=== Variance Components ===\n"); print(varcomp_tab, digits = 3)
cat("\n=== Broad-Sense Heritability (H2) ===\n"); print(h2_tab, digits = 3)

for (coh in c("2001-established", "2002-established")) {
  cat("\n=== Top 5 BLUPs for ", coh, " ===\n", sep = "")
  for (trait in TRAITS) {
    top <- blup_tab %>% filter(cohort == coh, trait == !!trait) %>% arrange(desc(BLUP)) %>% head(5)
    cat("  ", trait, ":", paste(sprintf("%s(%.2f)", top$genotype, top$BLUP), collapse = ", "), "\n")
  }
}

# ============================================================================
# 3. Variance Boundary Diagnostics for 2002 Fresh Weight (FW)
# ============================================================================
cat("\n=== Boundary Diagnosis: 2002 Fresh Weight V_G ===\n")
dd2002fw <- d %>% filter(cohort == "2002-established", !is.na(FW))
m2002fw <- lmer(FW ~ year_f + (1|genotype) + (1|genotype:year_f) + (1|plot),
                data = dd2002fw, REML = TRUE, control = lmerControl(calc.derivs = FALSE))
vc2002 <- as.data.frame(VarCorr(m2002fw))
cat("V_G:", vc2002$vcov[vc2002$grp == "genotype"][1], " singular:", isSingular(m2002fw), "\n")

ci <- tryCatch(confint(m2002fw, oldNames = FALSE, quiet = TRUE), error = function(e) NULL)
if (!is.null(ci)) {
  write.csv(as.data.frame(ci), file.path(out_dir, "2002_FW_variance_boundary_final.csv"))
  cat("\n=== 2002 FW Variance Component Profile CIs ===\n"); print(ci)
}

cat("\n=== All output files successfully written to:", out_dir, "===\n")
