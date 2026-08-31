# -*- coding: utf-8 -*-
# ============================================================================
# 02_stability.R
# Stability and sensitivity analyses: GxY interactions, single-site AMMI
# temporal stability, stand-age sensitivity, and cross-cohort comparisons.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Define relative paths for reproducibility and data availability
input_file <- file.path("output", "plot_year_wide.csv")
out_dir <- file.path("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(input_file, stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------
# 1. Genotype x Year (GxY) Raw Means
# ----------------------------------------------------------------------------
gxy <- d %>%
  group_by(cohort, genotype_ID, calendar_year) %>%
  summarise(
    FW = mean(FW, na.rm = TRUE),
    DW = mean(DW, na.rm = TRUE),
    SummerPH = mean(SummerPH, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(gxy, file.path(out_dir, "genotype_year_interaction_final.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 2. Single-Site Temporal AMMI Analysis (2001 Cohort, 4-Year Balanced Data)
# ----------------------------------------------------------------------------
ammi_pc <- list()
for (trait in c("SummerPH", "FW", "DW")) {
  M <- gxy %>%
    filter(cohort == "2001-established", !is.na(.data[[trait]])) %>%
    select(genotype_ID, calendar_year, all_of(trait)) %>%
    pivot_wider(names_from = calendar_year, values_from = all_of(trait))
  
  complete <- M[complete.cases(M), ]
  if (nrow(complete) < 2) next
  
  X <- as.matrix(complete[, -1])
  rownames(X) <- complete$genotype_ID
  
  gr <- rowMeans(X)
  cg <- colMeans(X)
  gm <- mean(X)
  
  I <- X - gr - rep(cg, each = nrow(X)) + gm
  s <- svd(I)
  ev <- s$d^2 / sum(s$d^2) * 100
  
  ammi_pc[[trait]] <- data.frame(
    trait = trait,
    n_geno = nrow(X),
    PC1 = ev[1],
    PC2 = ev[2],
    cum = ev[1] + ev[2],
    stringsAsFactors = FALSE
  )
}

ammi_tab <- do.call(rbind, ammi_pc)
rownames(ammi_tab) <- NULL
write.csv(ammi_tab, file.path(out_dir, "AMMI_temporal_stability_final.csv"), row.names = FALSE)

cat("=== AMMI Temporal Stability (2001 Cohort Complete Genotypes) ===\n")
print(ammi_tab, digits = 2)

# ----------------------------------------------------------------------------
# 3. Stand Age Sensitivity Analysis
# ----------------------------------------------------------------------------
d$stand_age <- d$calendar_year - d$establishment_year
same_age <- d %>%
  filter(stand_age %in% 1:3) %>%
  group_by(cohort, stand_age) %>%
  summarise(
    FW = mean(FW, na.rm = TRUE),
    DW = mean(DW, na.rm = TRUE),
    SummerPH = mean(SummerPH, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(same_age, file.path(out_dir, "same_stand_age_sensitivity_final.csv"), row.names = FALSE)

cat("\n=== Stand Age Sensitivity Analysis ===\n")
print(as.data.frame(same_age), digits = 3)

# ----------------------------------------------------------------------------
# 4. Cross-Cohort Descriptive Comparison (L33 & Algonquin)
# ----------------------------------------------------------------------------
cross <- list()
for (g in c("L33", "Algonquin")) {
  sub <- d %>%
    filter(genotype_ID == g) %>%
    group_by(cohort, stand_age, calendar_year) %>%
    summarise(
      FW = mean(FW, na.rm = TRUE),
      DW = mean(DW, na.rm = TRUE),
      SummerPH = mean(SummerPH, na.rm = TRUE),
      .groups = "drop"
    )
  cross[[g]] <- cbind(genotype = g, sub)
}

cross_tab <- do.call(rbind, cross)
rownames(cross_tab) <- NULL
write.csv(cross_tab, file.path(out_dir, "L33_Algonquin_cross_cohort_final.csv"), row.names = FALSE)

cat("\n=== Cross-Cohort Comparison (L33 / Algonquin) ===\n")
print(as.data.frame(cross_tab), digits = 2)

cat("\n=== Stability Analysis Completed ===\n")
