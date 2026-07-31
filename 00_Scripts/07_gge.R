# 07_gge.R
# GGE (Genotype + Genotype-by-Environment) biplot analysis.
# Environment-centred matrix -> SVD; symmetric scaling for biplot.
# Adds the average-environment (AE) axis to rank genotypes by
# mean performance (projection) and stability (perpendicular distance).
source("00_config.R")
gxe_list <- readRDS(file.path(INT, "gxe_list.rds"))

gge_one <- function(M) {
  yrs <- names(M)[names(M) != "Genotype"]
  X <- as.matrix(M[, yrs]); rownames(X) <- M$Genotype
  Xc <- scale(X, center = colMeans(X, na.rm = TRUE), scale = FALSE)  # env-centred
  s  <- svd(Xc)
  D <- s$d; U <- s$u; V <- s$v
  G <- U %*% diag(sqrt(D)); E <- V %*% diag(sqrt(D))   # genotype / environment scores
  G <- G[, 1:2, drop = FALSE]; E <- E[, 1:2, drop = FALSE]   # keep 2 PCs for biplot
  rownames(G) <- rownames(X); rownames(E) <- yrs
  list(geno_coord = G, env_coord = E, D = D, yrs = yrs, geno = rownames(X))
}

gge_res <- list()
for (t in TRAITS) {
  M    <- gxe_list[[t]]
  yrs  <- names(M)[names(M) != "Genotype"]
  comp <- M[complete.cases(M[, yrs]), ]
  g    <- gge_one(comp)
  g$n_used <- nrow(comp); g$n_total <- nrow(M)
  ae  <- colMeans(g$env_coord); ae <- ae / sqrt(sum(ae^2))   # unit AE vector
  proj <- as.vector(g$geno_coord %*% ae)                    # mean performance
  perp <- sqrt(rowSums(g$geno_coord^2) - proj^2)            # stability (distance)
  g$rank <- data.frame(Genotype = g$geno, PC1 = g$geno_coord[, 1],
                       PC2 = g$geno_coord[, 2], MeanPerf = proj,
                       Stability = perp, stringsAsFactors = FALSE)
  gge_res[[t]] <- g
}
saveRDS(gge_res, file.path(INT, "gge_res.rds"))

for (t in TRAITS) {
  g <- gge_res[[t]]
  write.csv(g$rank, file.path(TAB, paste0("GGE_genotype_coords_", t, ".csv")), row.names = FALSE)
}
cat("=== GGE: top 5 genotypes by mean performance (AE projection) ===\n")
for (t in TRAITS) {
  cat(sprintf("\n-- %s (n_used=%d/%d) --\n", t, gge_res[[t]]$n_used, gge_res[[t]]$n_total))
  print(gge_res[[t]]$rank %>% arrange(desc(MeanPerf)) %>% slice_head(n = 5) %>%
          mutate(PC1 = round(PC1,2), PC2 = round(PC2,2), MeanPerf = round(MeanPerf,2),
                 Stability = round(Stability,2)), row.names = FALSE)
}
cat("\nGGE coordinates -> 04_Results/Tables/\n")
