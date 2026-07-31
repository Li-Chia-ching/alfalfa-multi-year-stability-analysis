# 05_heritability.R
# Broad-sense heritability (H2) on a genotype-mean basis across tested
# environments (years), from the LMM variance components:
#   H2 = sigma2_G / (sigma2_G + sigma2_GE/n_env + sigma2_e/(n_env*n_rep))
# sigma2_e = residual variance; n_env = #years with data; n_rep = mean reps/cell.
source("00_config.R")
models <- readRDS(file.path(INT, "models.rds"))
core   <- readRDS(file.path(INT, "core_rep.rds"))   # rep-level experimental unit

get_vc <- function(m, grp) {
  v <- VarCorr(m)
  if (grp %in% names(v)) as.numeric(v[[grp]]) else 0
}

h2_rows <- list()
for (t in TRAITS) {
  m   <- models[[t]]
  sg2 <- get_vc(m, "Genotype")
  sge2<- get_vc(m, "Genotype:Year")
  se2 <- sigma(m)^2                                  # residual (plot error)
  d   <- core %>% filter(!is.na(.data[[t]]))
  n_env <- n_distinct(d$Year)
  n_rep <- mean(d %>% group_by(Genotype, Year) %>% summarise(k = n(), .groups = "drop") %>% pull(k))
  H2 <- sg2 / (sg2 + sge2 / n_env + se2 / (n_env * n_rep))
  h2_rows[[t]] <- data.frame(Trait = t, Unit = TRAIT_UNIT[t],
                             sigma2_G = sg2, sigma2_GE = sge2, sigma2_e = se2,
                             n_env = n_env, n_rep_mean = n_rep, H2 = H2,
                             stringsAsFactors = FALSE)
}
h2_all <- bind_rows(h2_rows)
h2_all$H2 <- round(h2_all$H2, 4)
write.csv(h2_all, file.path(TAB, "Heritability_Summary.csv"), row.names = FALSE)

cat("=== Broad-sense heritability ===\n")
print(h2_all, row.names = FALSE)
cat("\nHeritability -> 04_Results/Tables/Heritability_Summary.csv\n")
