# 04_estimate_blup.R
# Best Linear Unbiased Prediction (BLUP) of genotypic values from the LMM.
# BLUP(g) = population intercept + random Genotype effect.
# Also builds a wide BLUP profile matrix (used later for clustering).
source("00_config.R")
models <- readRDS(file.path(INT, "models.rds"))

blup_list <- list()
for (t in TRAITS) {
  m   <- models[[t]]
  re  <- ranef(m, drop = TRUE)[["Genotype"]]      # named vector, Genotype effects
  int <- fixef(m)[1]                              # population mean (Year baseline)
  val <- sort(int + re, decreasing = TRUE)
  blup_list[[t]] <- data.frame(Trait = t,
                               Genotype = names(val),
                               BLUP = as.numeric(val),
                               stringsAsFactors = FALSE)
}
blup_all  <- bind_rows(blup_list)
blup_wide <- blup_all %>% pivot_wider(names_from = Trait, values_from = BLUP)
write.csv(blup_all,  file.path(TAB, "BLUP_genotypic_values.csv"), row.names = FALSE)
saveRDS(blup_wide,   file.path(INT, "blup_wide.rds"))

cat("=== Top 5 genotypes by BLUP per trait ===\n")
for (t in TRAITS) {
  cat(sprintf("\n-- %s (%s) --\n", t, TRAIT_UNIT[t]))
  print(blup_all %>% filter(Trait == t) %>% slice_max(BLUP, n = 5) %>% mutate(BLUP = round(BLUP, 3)), row.names = FALSE)
}
cat("\nBLUP table -> 04_Results/Tables/BLUP_genotypic_values.csv\n")
