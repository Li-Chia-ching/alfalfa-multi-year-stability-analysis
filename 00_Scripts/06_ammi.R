# 06_ammi.R
# Additive Main Effects and Multiplicative Interaction (AMMI) analysis.
# For each trait: build complete genotype x environment (Year) matrix,
# remove additive main effects, SVD the interaction -> IPCAs.
# No external package; base SVD only (Reproducibility First).
source("00_config.R")
gxe_list <- readRDS(file.path(INT, "gxe_list.rds"))

ammi_one <- function(M) {
  yrs <- names(M)[names(M) != "Genotype"]
  X <- as.matrix(M[, yrs]); rownames(X) <- M$Genotype
  gr <- rowMeans(X, na.rm = TRUE); cg <- colMeans(X, na.rm = TRUE); gm <- mean(X, na.rm = TRUE)
  I  <- X - gr - rep(cg, each = nrow(X)) + gm          # interaction (additive removed)
  s  <- svd(I)
  list(I = I, U = s$u, V = s$v, D = s$d, years = yrs, geno = M$Genotype)
}

ammi_res <- list()
for (t in TRAITS) {
  M    <- gxe_list[[t]]
  yrs  <- names(M)[names(M) != "Genotype"]
  comp <- M[complete.cases(M[, yrs]), ]               # genotypes with all envs
  a    <- ammi_one(comp)
  a$n_used <- nrow(comp); a$n_total <- nrow(M)
  ammi_res[[t]] <- a
}
saveRDS(ammi_res, file.path(INT, "ammi_res.rds"))

# --- Export IPCA tables + variance explained -------------------------------
ammi_sum <- bind_rows(lapply(TRAITS, function(t) {
  a <- ammi_res[[t]]
  data.frame(Trait = t, PC = 1:length(a$D),
             Eigenvalue = a$D^2,
             CumVarProp = cumsum(a$D^2) / sum(a$D^2),
             GenoUsed = a$n_used, GenoTotal = a$n_total,
             stringsAsFactors = FALSE)
}))
write.csv(ammi_sum, file.path(TAB, "AMMI_variance_explained.csv"), row.names = FALSE)

for (t in TRAITS) {
  a <- ammi_res[[t]]; k <- min(2, length(a$D))
  gs <- data.frame(Genotype = a$geno, a$U[, 1:k, drop = FALSE])
  colnames(gs)[-1] <- paste0("IPCA", 1:k)
  es <- data.frame(Environment = a$years, a$V[, 1:k, drop = FALSE])
  colnames(es)[-1] <- paste0("IPCA", 1:k)
  write.csv(gs, file.path(TAB, paste0("AMMI_genotype_IPCA_", t, ".csv")), row.names = FALSE)
  write.csv(es, file.path(TAB, paste0("AMMI_env_IPCA_", t, ".csv")),   row.names = FALSE)
}
cat("=== AMMI variance explained (cumulative) ===\n")
print(ammi_sum, row.names = FALSE)
cat("\nAMMI tables -> 04_Results/Tables/\n")
