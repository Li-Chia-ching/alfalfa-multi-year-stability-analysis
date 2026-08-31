# -*- coding: utf-8 -*-
# ============================================================================
# 05_residual_diagnostics.R
# Diagnostic evaluation of primary linear mixed-effects models (lme4).
# Evaluates homoscedasticity and normality of residuals per cohort x trait model:
#   1. Residual-vs-fitted scatterplots and normal Q-Q plots (PNG format).
#   2. Quantitative diagnostic summary (Shapiro-Wilk test, Brown-Forsythe 
#      heteroscedasticity test, and extreme standardized residual counts).
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
})

# Define relative paths for reproducibility and data availability
out_dir <- file.path("output")
fig_dir <- file.path(out_dir, "diagnostics")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(out_dir, "plot_year_wide.csv")
d <- read.csv(input_file, stringsAsFactors = FALSE)
d <- d %>% mutate(
  year_f   = factor(calendar_year),
  plot     = factor(plot),
  genotype = factor(genotype_ID)
)

TRAITS  <- c("SummerPH", "FW", "DW")
COHORTS <- c("2001-established", "2002-established")

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

summary_rows <- list()
plots <- list()

for (coh in COHORTS) {
  for (trait in TRAITS) {
    dd <- d %>% filter(cohort == coh, !is.na(.data[[trait]]))
    if (nrow(dd) < 30) next
    
    m <- tryCatch(fit_primary(dd, trait), error = function(e) NULL)
    if (is.null(m)) {
      cat("Model fitting failed:", coh, trait, "\n")
      next
    }

    fitted_v <- fitted(m)
    resid_r  <- residuals(m, type = "response")
    resid_v  <- residuals(m, type = "pearson", scaled = TRUE)

    # Normality evaluation (Shapiro-Wilk)
    n_use <- min(length(resid_v), 5000)
    set.seed(20260831)
    samp  <- sample(seq_along(resid_v), n_use)
    sw    <- tryCatch(shapiro.test(resid_v[samp]), error = function(e) NULL)

    # Heteroscedasticity evaluation (Brown-Forsythe)
    ht <- hetero_test(resid_r, fitted_v)

    # Count of extreme standardized residuals
    ext <- sum(abs(resid_v) > 3)

    # Skewness and kurtosis
    n_obs <- length(resid_r)
    m2 <- mean((resid_r - mean(resid_r))^2)
    skew <- mean((resid_r - mean(resid_r))^3) / m2^1.5
    kurt <- mean((resid_r - mean(resid_r))^4) / m2^2 - 3

    summary_rows[[paste(coh, trait)]] <- data.frame(
      cohort            = coh,
      trait             = trait,
      n                 = nrow(dd),
      shapiro_W         = if (is.null(sw)) NA_real_ else unname(sw$statistic),
      shapiro_p         = if (is.null(sw)) NA_real_ else sw$p.value,
      bf_F              = unname(ht["stat"]),
      bf_p              = unname(ht["p"]),
      resid_skew        = skew,
      resid_excess_kurt = kurt,
      n_stdresid_gt3    = ext,
      pct_stdresid_gt3  = 100 * ext / nrow(dd),
      stringsAsFactors  = FALSE
    )

    plots[[paste(coh, trait)]] <- list(
      dd = dd, trait = trait, coh = coh,
      fitted = fitted_v, resid = resid_v,
      rawresid = resid_r
    )
    
    cat(sprintf("%-18s %-9s n=%3d  Shapiro p=%.3g  BF p=%.3g  |z|>3: %d\n",
                coh, trait, nrow(dd),
                if (is.null(sw)) NA else sw$p.value,
                unname(ht["p"]), ext))
  }
}

sum_tab <- do.call(rbind, summary_rows)
rownames(sum_tab) <- NULL
output_summary_file <- file.path(out_dir, "residual_diagnostics_summary.csv")
write.csv(sum_tab, output_summary_file, row.names = FALSE)

# ----------------------------------------------------------------------------
# Graphic Output: Residual Plots
# ----------------------------------------------------------------------------
npanel <- length(plots)
nc <- 2; nr <- ceiling(npanel / nc)

png(file.path(fig_dir, "resid_vs_fitted.png"), width = 1400, height = 460 * nr, res = 160)
par(mfrow = c(nr, nc), mar = c(4.2, 4.2, 2.6, 0.8), oma = c(0, 0, 2, 0))
for (p in plots) {
  plot(p$fitted, p$resid, pch = 19, cex = 0.45, col = "#0072B2",
       xlab = "Fitted values", ylab = "Pearson residual",
       main = paste0(p$coh, " | ", p$trait))
  abline(h = 0, col = "#D55E00", lwd = 1.2, lty = 2)
  lines(lowess(p$fitted, p$resid), col = "#E69F00", lwd = 1.4)
}
mtext("Residuals vs Fitted — Primary lme4 Model", outer = TRUE, line = 0.4, cex = 1.1)
dev.off()

png(file.path(fig_dir, "qq_residuals.png"), width = 1400, height = 460 * nr, res = 160)
par(mfrow = c(nr, nc), mar = c(4.2, 4.2, 2.6, 0.8), oma = c(0, 0, 2, 0))
for (p in plots) {
  qqnorm(p$resid, pch = 19, cex = 0.45, col = "#0072B2",
         xlab = "Theoretical quantiles", ylab = "Sample quantiles",
         main = paste0(p$coh, " | ", p$trait))
  qqline(p$resid, col = "#D55E00", lwd = 1.2)
}
mtext("Normal Q-Q — Pearson Residuals, Primary lme4 Model", outer = TRUE, line = 0.4, cex = 1.1)
dev.off()

cat("\n=== Residual Diagnostics Summary ===\n")
print(sum_tab, digits = 3)
cat("\nDiagnostic metrics written to:", output_summary_file, "\n")
cat("Diagnostic plots generated:\n  -", file.path(fig_dir, "resid_vs_fitted.png"),
    "\n  -", file.path(fig_dir, "qq_residuals.png"), "\n")
