# 10_make_figures.R
# Generate all manuscript figures with scaled-up, highly readable typography.
source("00_config.R")
library(ggdendro)

gxe_list   <- readRDS(file.path(INT, "gxe_list.rds"))
gge_res    <- readRDS(file.path(INT, "gge_res.rds"))
ammi_res   <- readRDS(file.path(INT, "ammi_res.rds"))
cluster_res<- readRDS(file.path(INT, "cluster_res.rds"))
blup_wide  <- readRDS(file.path(INT, "blup_wide.rds"))
h2         <- read.csv(file.path(TAB, "Heritability_Summary.csv"))
core       <- readRDS(file.path(INT, "core.rds"))

t <- "FW"; yrs <- names(gxe_list[[t]])[names(gxe_list[[t]]) != "Genotype"]
top <- gxe_list[[t]] %>% mutate(m = rowMeans(select(., all_of(yrs)), na.rm = TRUE)) %>%
  arrange(desc(m)) %>% slice_head(n = 10) %>% pull(Genotype)

# ---- Fig 1: Genotype x environment interaction -----------------------------
df1 <- gxe_list[[t]] %>% filter(Genotype %in% top) %>%
  pivot_longer(all_of(yrs), names_to = "Year", values_to = "value") %>%
  mutate(Year = factor(Year))

f1 <- ggplot(df1, aes(x = Year, y = value, group = Genotype, colour = Genotype)) +
  geom_line(linewidth = 1.2, alpha = 0.75) + 
  geom_point(size = 3.2, alpha = 0.9) +
  labs(y = "Fresh weight (g plot-1)", x = "Year (environment)") + 
  FIG_THEME + 
  scale_colour_viridis_d(option = "plasma", end = 0.9)

export_fig_with_data(f1, df1, "Fig1_FW_genotype_environment", w = FIG_W, h = FIG_H)

# ---- Fig 2: GGE biplot -----------------------------------------------------
g <- gge_res[["FW"]]
G <- as.data.frame(g$geno_coord); colnames(G) <- c("PC1","PC2"); G$Genotype <- g$geno
E <- as.data.frame(g$env_coord);  colnames(E) <- c("PC1","PC2"); E$Environment <- g$yrs
score <- as.matrix(g$geno_coord); env <- as.matrix(g$env_coord)
winners <- apply(env, 1, function(e) which.max(score %*% e))
win_idx <- unique(winners)
ang <- atan2(G$PC2[win_idx], G$PC1[win_idx]); ord <- win_idx[order(ang)]
poly <- G[ord, c("PC1","PC2")]
geno_ang <- atan2(G$PC2, G$PC1); win_ang <- atan2(G$PC2[win_idx], G$PC1[win_idx])
sector <- sapply(geno_ang, function(a) {
  d <- abs(((a - win_ang + pi) %% (2*pi)) - pi); win_idx[which.min(d)]
})
G$sector <- factor(sector)

f2 <- ggplot() +
  geom_polygon(data = poly, aes(x = PC1, y = PC2), fill = NA, colour = "black", linetype = 2, linewidth = 0.8, alpha = 0.5) +
  geom_segment(data = E, aes(x = 0, y = 0, xend = PC1*1.12, yend = PC2*1.12),
               arrow = arrow(length = unit(0.12,"cm")), colour = "#D55E00", linewidth = 1, alpha = 0.8) +
  geom_text(data = E, aes(x = PC1*1.18, y = PC2*1.18, label = Environment), colour = "#D55E00", size = 4.5, fontface = "bold") +
  geom_point(data = G, aes(x = PC1, y = PC2, colour = sector), size = 2.5, alpha = 0.75) +
  geom_text_repel(data = G %>% filter(Genotype %in% top), aes(x = PC1, y = PC2, label = Genotype),
                  size = 4, fontface = "bold", segment.color = 'grey50', alpha = 0.9) +
  labs(x = "PC1 (environment-centred)", y = "PC2") + 
  FIG_THEME + PAL + guides(colour = "none")

export_fig_with_data(f2, G, "Fig2_GGE_biplot_FW", w = FIG_W, h = FIG_H)

# ---- Fig 3: AMMI biplot ----------------------------------------------------
a <- ammi_res[["FW"]]
GA <- as.data.frame(a$U[, 1:2]); colnames(GA) <- c("IPCA1","IPCA2"); GA$Genotype <- a$geno
EA <- as.data.frame(a$V[, 1:2]); colnames(EA) <- c("IPCA1","IPCA2"); EA$Environment <- a$years

f3 <- ggplot() +
  geom_segment(data = EA, aes(x = 0, y = 0, xend = IPCA1, yend = IPCA2),
               arrow = arrow(length = unit(0.1,"cm")), colour = "#D55E00", linewidth = 1, alpha = 0.8) +
  geom_text(data = EA, aes(x = IPCA1*1.1, y = IPCA2*1.1, label = Environment), colour = "#D55E00", size = 4.5, fontface = "bold") +
  geom_point(data = GA, aes(x = IPCA1, y = IPCA2), size = 2.5, colour = "#0072B2", alpha = 0.7) +
  geom_text_repel(data = GA %>% filter(Genotype %in% top), aes(x = IPCA1, y = IPCA2, label = Genotype),
                  size = 4, fontface = "bold", segment.color = 'grey50') +
  labs(x = "IPCA1", y = "IPCA2") + 
  FIG_THEME

export_fig_with_data(f3, GA, "Fig3_AMMI_biplot_FW", w = FIG_W, h = FIG_H)

# ---- Fig 4: Circular cluster dendrogram (Enlarged Canvas & Font) -----------
dd <- dendro_data(cluster_res$hc, k = cluster_res$k)
leaf_cl <- setNames(cluster_res$clust_df$Cluster, cluster_res$clust_df$Genotype)
n_labels <- nrow(dd$labels)

label_data <- dd$labels %>%
  mutate(
    cluster = factor(leaf_cl[label]),
    raw_angle = 90 - (360 * x / n_labels),
    # 调整标签朝向：使文字始终可读
    hjust = ifelse(raw_angle < -90 | raw_angle > 90, 1, 0),
    final_angle = ifelse(raw_angle < -90 | raw_angle > 90, raw_angle + 180, raw_angle)
  )

f4 <- ggplot() +
  geom_segment(data = segment(dd), aes(x = x, y = y, xend = xend, yend = yend), 
               colour = "grey50", alpha = 0.6, linewidth = 0.6) +
  # 最普通的 geom_text：固定位置，使用角度和 hjust
  geom_text(data = label_data, 
            aes(x = x, y = -1.0,          # 将标签放在更外圈（原始为 -0.5）
                label = label, colour = cluster,
                angle = final_angle, hjust = hjust),
            size = 3.8, fontface = "bold") +
  coord_polar(theta = "x") +
  scale_y_reverse(expand = expansion(mult = c(0.2, 0.60))) +   # 底部留白更多（原始 0.15, 0.45）
  labs(x = NULL, y = NULL) + 
  FIG_THEME + PAL + guides(colour = "none") +
  theme(
    axis.line = element_blank(), axis.text = element_blank(), 
    axis.ticks = element_blank(), axis.title = element_blank(), 
    panel.grid = element_blank()
  )

# 画布尺寸略为增大，为标签提供物理空间
export_fig_with_data(f4, label_data, "Fig4_cluster_dendrogram", w = 28, h = 28)

# ---- Fig 5: Heritability ---------------------------------------------------
h2p <- h2 %>% mutate(Trait = factor(Trait, levels = TRAITS))

f5 <- ggplot(h2p, aes(x = Trait, y = H2, fill = Trait)) +
  geom_bar(stat = "identity", width = 0.4, alpha = 0.85) +
  geom_text(aes(label = formatC(round(H2, 3), format = "f", digits = 3)), vjust = -0.8, size = 4.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(y = "H2", x = "Trait") +
  FIG_THEME + PAL + guides(fill = "none")

export_fig_with_data(f5, h2p, "Fig5_heritability", w = FIG_W, h = FIG_H)

# ---- Fig S1: Biological Survival Rate --------------------------------------
survival_data <- core %>%
  left_join(cluster_res$clust_df %>% select(Genotype, Cluster), by = "Genotype") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster, Year) %>%
  summarise(
    Total_Plots = n(),
    Survived = sum(SummerPH > 0, na.rm = TRUE),
    Survival_Rate = Survived / Total_Plots,
    .groups = "drop"
  ) %>%
  mutate(Cluster = factor(Cluster))

fS1 <- ggplot(survival_data, aes(x = factor(Year), y = Survival_Rate, fill = Cluster)) +
  geom_col(position = position_dodge(width = 0.55), width = 0.45, alpha = 0.85) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Trial Year", y = "Survival Rate") + 
  FIG_THEME + PAL + 
  theme(legend.position = "bottom", legend.text = element_text(size = 14))

export_fig_with_data(fS1, survival_data, "FigS1_data_coverage_survival", w = 20, h = 14)

# ---- Fig S2: Trait correlation ---------------------------------------------
fS2 <- ggplot(blup_wide, aes(x = SummerPH, y = FW, colour = factor(cluster_res$clust_df$Cluster))) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = "Summer plant height BLUP (cm)", y = "Fresh weight BLUP (g plot-1)") +
  FIG_THEME + PAL + guides(colour = "none")

export_fig_with_data(fS2, blup_wide, "FigS2_trait_correlation", w = FIG_W, h = FIG_H)

cat("All figures exported with scaled typography and optimized dimensions.\n")