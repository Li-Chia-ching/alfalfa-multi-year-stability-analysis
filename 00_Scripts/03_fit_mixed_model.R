# 03_fit_mixed_model.R
# Fit a linear mixed model (LMM) for each core trait:
#   y = mu + Year(fixed) + Genotype(random) + Genotype:Year(random) + Rep:Year(random) + e
# Single location, Year = environment. Rep nested in Year = incomplete block.
# Outputs: ANOVA (Type III, lmerTest), variance components, residual diagnostics.
source("00_config.R")
core <- readRDS(file.path(INT, "core_rep.rds"))   # experimental unit = Genotype x Rep x Year

fit_one <- function(trait) {
  d <- core %>% filter(!is.na(.data[[trait]])) %>%
    mutate(Year  = factor(Year),
           Genotype = factor(Genotype),
           RepYr = interaction(Rep, Year))   # block within year
  f <- as.formula(paste0(trait, " ~ Year + (1|Genotype) + (1|Genotype:Year) + (1|RepYr)"))
  lmer(f, data = d, REML = TRUE, control = lmerControl(calc.derivs = FALSE))
}

models <- setNames(lapply(TRAITS, fit_one), TRAITS)
saveRDS(models, file.path(INT, "models.rds"))

# --- ANOVA (Type III) + variance components --------------------------------
anova_all <- bind_rows(lapply(TRAITS, function(t) {
  a <- as.data.frame(anova(models[[t]])); a$term <- rownames(a); a$trait <- t; a
}))
varcomp_all <- bind_rows(lapply(TRAITS, function(t) {
  vc <- as.data.frame(VarCorr(models[[t]])); vc$trait <- t; vc
}))
write.csv(anova_all,  file.path(STAT, "ANOVA_mixed_model.csv"),  row.names = FALSE)
write.csv(varcomp_all, file.path(STAT, "Variance_Components.csv"), row.names = FALSE)

cat("=== ANOVA (Type III) ===\n"); print(anova_all[, c("trait","term","Sum Sq","Mean Sq","NumDF","DenDF","F value","Pr(>F)")], row.names = FALSE)
cat("\n=== Variance components ===\n"); print(varcomp_all[, c("trait","grp","var1","vcov","sdcor")], row.names = FALSE)

# --- Residual diagnostics (saved to Supplementary) -------------------------
for (t in TRAITS) {
  m <- models[[t]]
  rs <- data.frame(fitted = fitted(m), resid = resid(m))
  p <- ggplot(rs, aes(sample = resid)) + stat_qq(distribution = qnorm) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
    labs(title = paste("QQ plot of residuals -", TRAIT_LABEL[t]), x = "Theoretical quantiles", y = "Sample quantiles") +
    FIG_THEME
  ggsave(file.path(SUPP, paste0("FigS_resid_qq_", t, ".png")), p, width = 12, height = 10, units = "cm", dpi = 300)
  p2 <- ggplot(rs, aes(x = fitted, y = resid)) + geom_point(size = 0.6, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, colour = "red") +
    labs(title = paste("Residuals vs fitted -", TRAIT_LABEL[t]), x = "Fitted", y = "Residual") + FIG_THEME
  ggsave(file.path(SUPP, paste0("FigS_resid_fit_", t, ".png")), p2, width = 12, height = 10, units = "cm", dpi = 300)
}
cat("\nResidual diagnostics written to 04_Results/Supplementary/\n")
