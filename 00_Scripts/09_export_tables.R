# 09_export_tables.R
# Compile a master genotype summary (BLUP x traits + cluster + GGE stability)
# and a shortlist of recommended top cultivars (high-yield, stable cluster).
source("00_config.R")
blup_wide <- readRDS(file.path(INT, "blup_wide.rds"))
clust     <- read.csv(file.path(TAB, "Cluster_membership.csv"))
gge_fw    <- read.csv(file.path(TAB, "GGE_genotype_coords_FW.csv"))

master <- blup_wide %>%
  left_join(clust %>% select(Genotype, Cluster), by = "Genotype") %>%
  left_join(gge_fw %>% select(Genotype, MeanPerf, Stability), by = "Genotype")
write.csv(master, file.path(TAB, "Genotype_Summary.csv"), row.names = FALSE)

# Recommended shortlist: highest fresh-weight BLUP within the top cluster
rec <- master %>% filter(Cluster == 1) %>% arrange(desc(FW)) %>% slice_head(n = 10) %>%
  select(Genotype, SummerPH, FW, DW, Cluster, MeanPerf, Stability)
write.csv(rec, file.path(TAB, "Recommended_Top_Cultivars.csv"), row.names = FALSE)

cat("Master genotype summary rows:", nrow(master), "\n")
cat("Recommended top cultivars (Cluster 1, by FW):\n")
print(rec, row.names = FALSE)
cat("\nTables -> 04_Results/Tables/\n")
