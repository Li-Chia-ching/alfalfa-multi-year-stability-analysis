# 08_cluster.R
# Objective: Cluster genotypes based on their multi-trait BLUP values (performance profile analysis).
# Steps: Trait standardization (z-score) -> Euclidean distance -> Ward's hierarchical clustering.
# The optimal number of clusters (k) is determined by maximizing the average silhouette width (k = 2~6).
# Output: Cluster membership table, cluster centers table, and cluster result object (including dendrogram, etc.).

# Load configuration (includes paths, trait lists, etc.)
source("00_config.R")

# Read wide-format BLUP data (one genotype per row, one trait per column)
blup_wide <- readRDS(file.path(INT, "blup_wide.rds"))

# ---------- CRITICAL FIX: Set genotype names as matrix row names ----------
# The original code only extracted trait columns, losing genotype identifiers, resulting in numeric labels in plots.
# Now, select "Genotype" and all trait columns first, then use column_to_rownames to set genotype names as row names.
library(dplyr)
library(tibble)   # Provides column_to_rownames

mat <- blup_wide %>%
  select(Genotype, all_of(TRAITS)) %>%    # Keep genotype column
  column_to_rownames("Genotype") %>%      # Convert Genotype column to row names
  as.matrix()                             # Convert to numeric matrix (row names are now genotype names)

# Standardize each trait (mean = 0, SD = 1) to avoid scaling effects
zs <- scale(mat)

# Calculate Euclidean distance between genotypes
d <- dist(zs, method = "euclidean")

# Perform hierarchical clustering using Ward's method (minimizing within-cluster variance)
hc <- hclust(d, method = "ward.D2")

# ---------- Select optimal k via silhouette width ----------
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(cluster)

k_grid <- 2:6
sil <- sapply(k_grid, function(k) {
  cl <- cutree(hc, k)
  if (length(unique(cl)) < 2) {
    -1   # If there is only one cluster, silhouette width is meaningless
  } else {
    mean(silhouette(cl, d)[, 3])   # Get the average silhouette width
  }
})

k_best <- k_grid[which.max(sil)]
cl <- cutree(hc, k_best)   # Final cluster assignment

# ---------- Generate result tables ----------
# Merge cluster results back into the original dataframe (keeping the genotype column)
clust_df <- blup_wide %>%
  mutate(Cluster = cl)

# Calculate the center for each cluster (mean for each trait) and intra-cluster sample size
centres <- clust_df %>%
  group_by(Cluster) %>%
  summarise(
    across(all_of(TRAITS), mean),
    n = n(),
    .groups = "drop"
  )

# Write CSV tables
write.csv(clust_df, file.path(TAB, "Cluster_membership.csv"), row.names = FALSE)
write.csv(centres,  file.path(TAB, "Cluster_centres.csv"),   row.names = FALSE)

# ---------- Console output summary ----------
cat("Silhouette by k:", paste(k_grid, round(sil, 3), sep = ":", collapse = "  "), "\n")
cat("Selected k =", k_best, "\n")
cat("\n=== Cluster centres (BLUP scale) ===\n")
print(centres, row.names = FALSE)
cat("\nCluster membership -> 04_Results/Tables/Cluster_membership.csv\n")

# Save the cluster result object for downstream plotting (dendrogram, cluster assignment, etc.)
saveRDS(
  list(
    hc = hc,           # Hierarchical clustering object (leaf labels are now genotype names)
    k = k_best,
    clust_df = clust_df,
    centres = centres,
    zs = zs
  ),
  file.path(INT, "cluster_res.rds")
)