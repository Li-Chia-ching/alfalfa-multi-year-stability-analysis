# 08_cluster.R
# Cluster genotypes by their multi-trait performance profile (BLUP values).
# Standardise each trait (z-score), Euclidean distance, Ward linkage.
# Optimal k chosen by maximising average silhouette width (k = 2..6).
source("00_config.R")
blup_wide <- readRDS(file.path(INT, "blup_wide.rds"))

mat <- blup_wide %>% select(all_of(TRAITS)) %>% as.matrix()
zs  <- scale(mat)                                  # standardise columns
d   <- dist(zs, method = "euclidean")
hc  <- hclust(d, method = "ward.D2")

# optimal k by silhouette
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(cluster)
k_grid <- 2:6
sil <- sapply(k_grid, function(k) {
  cl <- cutree(hc, k)
  if (length(unique(cl)) < 2) -1 else mean(silhouette(cl, d)[, 3])
})
k_best <- k_grid[which.max(sil)]
cl <- cutree(hc, k_best)

clust_df <- blup_wide %>% mutate(Cluster = cl)
centres <- clust_df %>% group_by(Cluster) %>%
  summarise(across(all_of(TRAITS), mean), n = n(), .groups = "drop")

write.csv(clust_df, file.path(TAB, "Cluster_membership.csv"), row.names = FALSE)
write.csv(centres,  file.path(TAB, "Cluster_centres.csv"),   row.names = FALSE)

cat("Silhouette by k:", paste(k_grid, round(sil,3), sep=":", collapse="  "), "\n")
cat("Selected k =", k_best, "\n")
cat("\n=== Cluster centres (BLUP scale) ===\n")
print(centres, row.names = FALSE)
cat("\nCluster membership -> 04_Results/Tables/Cluster_membership.csv\n")
saveRDS(list(hc = hc, k = k_best, clust_df = clust_df, centres = centres, zs = zs),
        file.path(INT, "cluster_res.rds"))
