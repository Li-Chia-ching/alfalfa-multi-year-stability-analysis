# -*- coding: utf-8 -*-
# ============================================================================
# 03_revision_checks.R
# Revision verification and sensitivity checks:
#   1) Within-cohort BLUP Pearson correlation matrix across traits.
#   2) Single-site GGE temporal stability analysis (2001 cohort).
#   3) Cohort-level summary statistics of genotype means.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Define relative paths for reproducibility and data availability
out_dir <- file.path("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 1. Within-Cohort BLUP Pearson Correlations
# ----------------------------------------------------------------------------
blup_file <- file.path(out_dir, "BLUP_final_by_cohort.csv")
blup <- read.csv(blup_file, stringsAsFactors = FALSE)

traits <- c("SummerPH", "FW", "DW")
pairs <- list(c("SummerPH", "FW"), c("SummerPH", "DW"), c("FW", "DW"))
cor_rows <- list()

for (coh in c("2001-established", "2002-established")) {
  for (pr in pairs) {
    t1 <- pr[1]; t2 <- pr[2]
    x <- blup %>% filter(cohort == coh, trait == t1) %>% select(genotype, BLUP) %>% rename(x = BLUP)
    y <- blup %>% filter(cohort == coh, trait == t2) %>% select(genotype, BLUP) %>% rename(y = BLUP)
    m <- merge(x, y, by = "genotype")
    if (nrow(m) < 3) next
    
    ct <- cor.test(m$x, m$y)
    cor_rows[[paste(coh, t1, t2)]] <- data.frame(
      cohort = coh, trait1 = t1, trait2 = t2, n = nrow(m),
      r = as.numeric(ct$estimate), p = ct$p.value, stringsAsFactors = FALSE
    )
  }
}

cor_tab <- do.call(rbind, cor_rows)
rownames(cor_tab) <- NULL
write.csv(cor_tab, file.path(out_dir, "BLUP_correlation_by_cohort.csv"), row.names = FALSE)

cat("=== Within-Cohort BLUP Pearson Correlations ===\n")
print(cor_tab, digits = 4)

# ----------------------------------------------------------------------------
# 2. Single-Site Temporal GGE Stability Analysis (2001 Cohort)
# ----------------------------------------------------------------------------
input_file <- file.path(out_dir, "plot_year_wide.csv")
d <- read.csv(input_file, stringsAsFactors = FALSE)

gxy <- d %>%
  group_by(cohort, genotype_ID, calendar_year) %>%
  summarise(
    SummerPH = mean(SummerPH, na.rm = TRUE),
    FW = mean(FW, na.rm = TRUE),
    DW = mean(DW, na.rm = TRUE),
    .groups = "drop"
  )

gge_rows <- list()
for (trait in c("SummerPH", "FW", "DW")) {
  M <- gxy %>%
    filter(cohort == "2001-established", !is.na(.data[[trait]])) %>%
    select(genotype_ID, calendar_year, all_of(trait)) %>%
    pivot_wider(names_from = calendar_year, values_from = all_of(trait))
  
  complete <- M[complete.cases(M), ]
  X <- as.matrix(complete[, -1])
  rownames(X) <- complete$genotype_ID
  
  # Environment (year) centering: G + GxY
  Xc <- sweep(X, 2, colMeans(X))
  s <- svd(Xc)
  ev <- s$d^2 / sum(s$d^2) * 100
  
  gge_rows[[trait]] <- data.frame(
    trait = trait,
    n_geno = nrow(X),
    PC1 = ev[1],
    PC2 = ev[2],
    cum = ev[1] + ev[2],
    stringsAsFactors = FALSE
  )
}

gge_tab <- do.call(rbind, gge_rows)
rownames(gge_tab) <- NULL
write.csv(gge_tab, file.path(out_dir, "GGE_temporal_stability_final.csv"), row.names = FALSE)

cat("\n=== Single-Site GGE Temporal Stability (2001 Cohort Complete Genotypes) ===\n")
print(gge_tab, digits = 2)

# ----------------------------------------------------------------------------
# 3. Cohort-Level Genotype Mean Descriptive Summary
# ----------------------------------------------------------------------------
gmean <- d %>%
  group_by(cohort, genotype_ID) %>%
  summarise(
    SummerPH = mean(SummerPH, na.rm = TRUE),
    FW = mean(FW, na.rm = TRUE),
    DW = mean(DW, na.rm = TRUE),
    n_year = n_distinct(calendar_year),
    .groups = "drop"
  )

cohort_summary <- gmean %>%
  group_by(cohort) %>%
  summarise(
    n_geno = n(),
    FW_mean = mean(FW),
    FW_median = median(FW),
    DW_mean = mean(DW),
    DW_median = median(DW),
    SummerPH_mean = mean(SummerPH),
    SummerPH_median = median(SummerPH),
    .groups = "drop"
  )

write.csv(cohort_summary, file.path(out_dir, "cohort_genotype_mean_summary.csv"), row.names = FALSE)

cat("\n=== Cohort Genotype-Level Summary Statistics ===\n")
print(as.data.frame(cohort_summary), digits = 3)

cat("\n=== Verification and Sensitivity Checks Completed ===\n")
