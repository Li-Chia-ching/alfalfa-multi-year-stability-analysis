# ============================================================================
# make_figures_revised.R — manuscript figures
#
# Main figures (Fig 1-4) are based on ALL 86 genotypes via a linear mixed
# model (unbalanced, REML). Supplementary figures (Fig S1-S3) are based on the
# 43 complete-genotype subset (AMMI / GGE require a balanced genotype x year
# matrix). This matches the manuscript structure: 86-genotype mixed model for
# heritability / BLUP / ranking, 43-genotype AMMI/GGE as exploratory appendix.
#
#   Fig. 1   Variance partitioning + broad-sense heritability (86 genotypes)
#   Fig. 2   Pearson correlation heatmap of trait BLUPs (86 genotypes)
#   Fig. 3   BLUP ranking with 95% CI (86 genotypes)
#   Fig. 4   2001 vs 2002 establishment-cohort comparison (all genotypes)
#   Fig. S1  AMMI biplot (43 complete genotypes)
#   Fig. S2  GGE biplot (43 complete genotypes)
#   Fig. S3  Mean performance vs. ASV / GGE distance (43 complete genotypes)
#
# H2 uses the genotype-mean basis formula, identical to analysis_86.R:
#   H2 = sigma_G^2 / [sigma_G^2 + sigma_GxY^2/Y + sigma_e^2/(Y*r)]
# where Y = number of years and r = harmonic mean of replications per
# genotype x year cell (read from the data, not hard-coded). The Year x Rep
# block component is a nuisance term and is excluded from the H2 denominator.
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

ROOT <- "C:/Users/lijia/Documents/R Workplace/Sino-Australian_Alfalfa_Project"
RAW  <- file.path(ROOT, "01_Raw_Phenotype_Data")
FIGD <- file.path(ROOT, "04_Results", "Manuscript_Figures")
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("SummerHeight", "FreshWeight", "DryWeight")
TRAIT_LABEL <- c(
  SummerHeight = "Summer plant height (cm)",
  FreshWeight  = "Fresh weight (kg 5 m\u00B2)",
  DryWeight    = "Dry weight (kg 5 m\u00B2)"
)
TRAIT_SHORT <- c(
  SummerHeight = "Summer height",
  FreshWeight  = "Fresh weight",
  DryWeight    = "Dry weight"
)
Y <- 4L

# Consistent Arial sans-serif. The previous "odd" look came from (a) no pinned
# family (device default) and (b) literal-bracket legend labels. We pin Arial
# and render via cairo for uniform anti-aliasing and Unicode coverage.
FONT <- "Arial"

# Okabe-Ito colourblind palette.
cb <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

save_fig <- function(p, name, w, h) {
  # PNG via ragg (auto-selected by ggsave; high-quality anti-aliasing + Arial).
  ggsave(file.path(FIGD, paste0(name, ".png")), p,
         width = w, height = h, units = "in", dpi = 300, bg = "white")
  # PDF via cairo for embedded Arial glyphs.
  ggsave(file.path(FIGD, paste0(name, ".pdf")), p,
         width = w, height = h, units = "in", bg = "white",
         device = grDevices::cairo_pdf)
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------
read_year <- function(y) {
  f <- file.path(RAW, paste0("Sino-Australian Alfalfa Project Data - ", y, "data.csv"))
  d <- read.csv(f, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  d$Year <- y
  d
}

raw <- do.call(rbind, lapply(2002:2005, read_year))
num <- function(x) suppressWarnings(as.numeric(x))
raw$Genotype <- trimws(raw$Line)
raw$Rep <- as.integer(raw$Rep)
raw$Plot <- num(raw$Plot)
raw$SummerHeight <- num(raw$Mean_Summer)
raw$FreshWeight  <- num(raw$FW)
raw$DryWeight    <- num(raw$DW)

# 2002 FW/DW conditional swap. A temporary copy of FW is saved BEFORE the
# assignment so the two columns cannot collapse to the same value.
swap <- raw$Year == 2002 &
  !is.na(raw$FreshWeight) & !is.na(raw$DryWeight) &
  raw$DryWeight > raw$FreshWeight
fw_old <- raw$FreshWeight[swap]
raw$FreshWeight[swap] <- raw$DryWeight[swap]
raw$DryWeight[swap]   <- fw_old

# Establishment cohort is a PLOT-level attribute (Plot <= 129 = 2001 block,
# > 129 = 2002 block). Genotype L33 has plots in BOTH blocks (6 plots).
raw$Estab <- ifelse(raw$Plot <= 129, 2001L, 2002L)

# Experimental unit = Genotype x Rep x Year. Estab is NOT part of this grouping
# so L33 collapses to its three replications (matching analysis_86.R). The
# plot-level Estab is kept in `raw` and used only by Fig. 4.
core_rep <- raw %>%
  group_by(Year, Genotype, Rep) %>%
  summarise(across(all_of(TRAITS), mean_or_na), .groups = "drop")

# 43 genotypes with complete data for all three traits in all four years.
complete43 <- core_rep %>%
  group_by(Genotype, Year) %>%
  summarise(across(all_of(TRAITS), ~ any(!is.na(.))), .groups = "drop") %>%
  group_by(Genotype) %>%
  summarise(n_complete = sum(SummerHeight & FreshWeight & DryWeight), .groups = "drop") %>%
  filter(n_complete == Y) %>%
  pull(Genotype)

if (length(complete43) != 43L) {
  stop(sprintf("Expected 43 complete genotypes, found %d.", length(complete43)))
}

# ---------------------------------------------------------------------------
# Main analysis: ALL 86 genotypes (unbalanced, REML)
# ---------------------------------------------------------------------------
d86 <- core_rep %>%
  mutate(
    Genotype = factor(Genotype),
    Year = factor(Year, levels = as.character(2002:2005)),
    Rep = factor(Rep)
  )

N_YEAR <- nlevels(d86$Year)

get_vc <- function(m, grp) {
  vc <- as.data.frame(VarCorr(m))
  v <- vc$vcov[vc$grp == grp]
  if (length(v)) v[1] else 0
}

harm_mean_rep <- function(dat, t) {
  k <- dat %>%
    filter(!is.na(.data[[t]])) %>%
    group_by(Genotype, Year) %>%
    summarise(n = n(), .groups = "drop") %>%
    pull(n)
  1 / mean(1 / k)
}

fit_models <- setNames(lapply(TRAITS, function(t) {
  d <- d86 %>% filter(!is.na(.data[[t]]))
  f <- as.formula(
    paste0(t, " ~ Year + (1|Genotype) + (1|Genotype:Year) + (1|Year:Rep)")
  )
  lmer(f, data = d, REML = TRUE, control = lmerControl(calc.derivs = FALSE))
}), TRAITS)

vc_tab <- do.call(rbind, lapply(TRAITS, function(t) {
  m <- fit_models[[t]]
  vg  <- get_vc(m, "Genotype")
  vgy <- get_vc(m, "Genotype:Year")
  vry <- get_vc(m, "Year:Rep")
  ve  <- sigma(m)^2
  r   <- harm_mean_rep(d86, t)
  h2  <- vg / (vg + vgy / N_YEAR + ve / (N_YEAR * r))
  data.frame(
    Trait = t, V_G = vg, V_GY = vgy, V_RY = vry, V_e = ve,
    r_harmonic = r, H2 = h2, stringsAsFactors = FALSE
  )
}))

blup_list <- lapply(TRAITS, function(t) {
  m <- fit_models[[t]]
  rr <- ranef(m, condVar = TRUE)
  se <- sqrt(attr(rr$Genotype, "postVar")[1, 1, ])
  data.frame(
    Genotype = rownames(rr$Genotype),
    BLUP = as.numeric(fixef(m)[1] + rr$Genotype[[1]]),
    SE = se,
    stringsAsFactors = FALSE
  )
})
names(blup_list) <- TRAITS

# ---------------------------------------------------------------------------
# Supplementary analysis: AMMI / GGE on the 43 complete genotypes
# ---------------------------------------------------------------------------
d43 <- core_rep %>%
  filter(Genotype %in% complete43) %>%
  mutate(
    Genotype = factor(Genotype),
    Year = factor(Year, levels = as.character(2002:2005)),
    Rep = factor(Rep)
  )

gxe <- function(t) {
  d43 %>%
    filter(!is.na(.data[[t]])) %>%
    group_by(Genotype, Year) %>%
    summarise(mu = mean(.data[[t]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Year, values_from = mu) %>%
    arrange(Genotype)
}

# SVD-based AMMI. "scores" = U*D (IPCA scores used for ASV);
# "biplot_*" = U*sqrt(D), V*sqrt(D) (symmetric biplot scaling).
ammi <- lapply(TRAITS, function(t) {
  M <- gxe(t)
  yrs <- setdiff(names(M), "Genotype")
  X <- as.matrix(M[, yrs])
  rownames(X) <- M$Genotype
  gr <- rowMeans(X); cg <- colMeans(X); gm <- mean(X)
  I <- X - gr - rep(cg, each = nrow(X)) + gm
  s <- svd(I)
  nd <- length(s$d)
  list(
    U = s$u, V = s$v, D = s$d,
    scores = s$u %*% diag(s$d, nrow = nd, ncol = nd),
    biplot_geno = s$u %*% diag(sqrt(s$d), nrow = nd, ncol = nd),
    biplot_env  = s$v %*% diag(sqrt(s$d), nrow = nd, ncol = nd),
    geno = M$Genotype, yrs = yrs
  )
})
names(ammi) <- TRAITS
ammi_pct <- sapply(TRAITS, function(t) {
  d <- ammi[[t]]$D^2; d / sum(d) * 100
})

# SVD-based GGE: environment-centred genotype means.
gge <- lapply(TRAITS, function(t) {
  M <- gxe(t)
  yrs <- setdiff(names(M), "Genotype")
  X <- as.matrix(M[, yrs])
  rownames(X) <- M$Genotype
  Xc <- scale(X, center = colMeans(X), scale = FALSE)
  s <- svd(Xc)
  nd <- length(s$d)
  list(
    G = s$u %*% diag(sqrt(s$d), nrow = nd, ncol = nd),
    E = s$v %*% diag(sqrt(s$d), nrow = nd, ncol = nd),
    D = s$d, geno = M$Genotype, yrs = yrs, M = M
  )
})
names(gge) <- TRAITS
gge_pct <- sapply(TRAITS, function(t) {
  d <- gge[[t]]$D^2; d / sum(d) * 100
})

# ---------------------------------------------------------------------------
# Fig. 1 — Variance partitioning + broad-sense heritability (86 genotypes)
# ---------------------------------------------------------------------------
vc_long <- vc_tab %>%
  mutate(Trait = factor(Trait, levels = TRAITS)) %>%
  pivot_longer(c(V_G, V_GY, V_RY, V_e), names_to = "Component", values_to = "value") %>%
  group_by(Trait) %>%
  mutate(prop = value / sum(value) * 100) %>%
  ungroup() %>%
  mutate(
    Component = factor(
      Component,
      levels = c("V_G", "V_GY", "V_RY", "V_e")
    )
  )

p1a <- ggplot(vc_long, aes(x = Trait, y = prop, fill = Component)) +
  geom_col(width = 0.58) +
  scale_fill_manual(
    values = c(cb[5], cb[2], cb[1], cb[8]),
    name = "Variance component",
    labels = parse(text = c("V[G]", "V[G %*% Y]", "V[Y %*% Rep]", "V[e]"))
  ) +
  scale_x_discrete(labels = TRAIT_SHORT) +
  labs(x = NULL, y = "Proportion of total variance (%)") +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(legend.position = "top")

p1b <- ggplot(vc_tab, aes(x = Trait, y = H2)) +
  geom_col(width = 0.58, fill = cb[5]) +
  geom_text(aes(label = sprintf("%.3f", H2)), vjust = -0.35, size = 3.4) +
  scale_x_discrete(labels = TRAIT_SHORT) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.08)),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(x = NULL, y = expression("Broad-sense heritability (" * H^2 * ")")) +
  theme_classic(base_size = 11, base_family = FONT)

p1 <- p1a + p1b + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
save_fig(p1, "Fig1_variance_heritability", 10, 4.5)

# ---------------------------------------------------------------------------
# Fig. 2 — Pearson correlation heatmap of trait BLUPs (86 genotypes)
# ---------------------------------------------------------------------------
blup_wide <- do.call(rbind, lapply(TRAITS, function(t) {
  x <- blup_list[[t]]; x$Trait <- t; x
})) %>%
  select(Genotype, Trait, BLUP) %>%
  pivot_wider(names_from = Trait, values_from = BLUP)

cor_mat <- cor(blup_wide[, TRAITS], use = "pairwise.complete.obs", method = "pearson")
cor_df <- as.data.frame(as.table(cor_mat)) %>%
  rename(Trait1 = Var1, Trait2 = Var2, r = Freq) %>%
  mutate(
    Trait1 = factor(Trait1, levels = rev(TRAITS)),
    Trait2 = factor(Trait2, levels = TRAITS)
  )

p2 <- ggplot(cor_df, aes(x = Trait2, y = Trait1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(
    low = cb[6], mid = "white", high = cb[5],
    midpoint = 0, limits = c(-1, 1), name = "Pearson r"
  ) +
  scale_x_discrete(labels = TRAIT_SHORT) +
  scale_y_discrete(labels = TRAIT_SHORT) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.line = element_blank()
  )
save_fig(p2, "Fig2_blup_correlation_heatmap", 5.8, 5.0)

# ---------------------------------------------------------------------------
# Fig. 3 — BLUP ranking + 95% CI (86 genotypes)
# ---------------------------------------------------------------------------
blup_rank <- do.call(rbind, lapply(TRAITS, function(t) {
  b <- blup_list[[t]] %>% arrange(desc(BLUP))
  b$Trait <- t
  b$rank <- seq_len(nrow(b))
  b$ci_lo <- b$BLUP - 1.96 * b$SE
  b$ci_hi <- b$BLUP + 1.96 * b$SE
  b
}))
blup_rank$Trait <- factor(blup_rank$Trait, levels = TRAITS)
blup_rank$top5 <- blup_rank$rank <= 5L

p3 <- ggplot(blup_rank, aes(x = rank, y = BLUP, colour = top5)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.3, linewidth = 0.3, alpha = 0.7) +
  geom_point(size = 0.9) +
  geom_text_repel(
    data = subset(blup_rank, top5),
    aes(label = Genotype), size = 2.6, colour = cb[6],
    nudge_x = 1.2, direction = "y", segment.size = 0.2, max.overlaps = 20
  ) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 1,
             labeller = labeller(Trait = TRAIT_LABEL)) +
  scale_colour_manual(values = c(`FALSE` = cb[8], `TRUE` = cb[6]), guide = "none") +
  labs(x = "Rank (sorted by BLUP)", y = "BLUP genotypic value") +
  theme_classic(base_size = 10, base_family = FONT) +
  theme(strip.text = element_text(face = "bold", size = 9))
save_fig(p3, "Fig3_blup_ranking", 7, 8)

# ---------------------------------------------------------------------------
# Fig. 4 — Establishment-cohort comparison (2001 vs 2002, all genotypes)
# ---------------------------------------------------------------------------
# Cohort is a plot-level attribute; L33 has plots in both cohorts, so it
# contributes one genotype-level mean to each. Aggregate to the experimental
# unit (Genotype x Rep x Year) first, then to a genotype x cohort mean over the
# 2003-2005 fair-comparison window.
plot_rep <- raw %>%
  group_by(Year, Genotype, Rep, Estab) %>%
  summarise(across(all_of(TRAITS), mean_or_na), .groups = "drop")

estab_means <- plot_rep %>%
  filter(Year >= 2003, Year <= 2005) %>%
  group_by(Genotype, Estab) %>%
  summarise(across(all_of(TRAITS), mean_or_na), .groups = "drop") %>%
  mutate(
    Cohort = factor(
      Estab, levels = c(2001, 2002),
      labels = c("Established in 2001", "Established in 2002")
    )
  )

cohort_pvals <- sapply(TRAITS, function(t) {
  a <- estab_means[[t]][estab_means$Estab == 2001]
  b <- estab_means[[t]][estab_means$Estab == 2002]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 2L || length(b) < 2L) return(NA_real_)
  t.test(a, b, var.equal = FALSE)$p.value
})

cohort_sig <- ifelse(
  is.na(cohort_pvals), "NA",
  ifelse(cohort_pvals < 0.001, "***",
         ifelse(cohort_pvals < 0.01, "**",
                ifelse(cohort_pvals < 0.05, "*", "ns")))
)

cohort_long <- estab_means %>%
  select(Genotype, Cohort, all_of(TRAITS)) %>%
  pivot_longer(all_of(TRAITS), names_to = "Trait", values_to = "value") %>%
  mutate(Trait = factor(Trait, levels = TRAITS))

cohort_y <- sapply(TRAITS, function(t) {
  x <- cohort_long$value[cohort_long$Trait == t]
  max(x, na.rm = TRUE) * 1.08
})

sig_df <- data.frame(
  Trait = factor(TRAITS, levels = TRAITS),
  y = as.numeric(cohort_y),
  label = unname(cohort_sig),
  stringsAsFactors = FALSE
)

p4 <- ggplot(cohort_long, aes(x = Cohort, y = value, fill = Cohort)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.4) +
  facet_wrap(~ Trait, scales = "free_y", labeller = labeller(Trait = TRAIT_LABEL)) +
  geom_text(data = sig_df, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c(cb[5], cb[6]), guide = "none") +
  labs(x = NULL, y = "Genotypic mean (2003\u20132005)") +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(strip.text = element_text(face = "bold"))
save_fig(p4, "Fig4_establishment_cohorts", 10, 4)

# ---------------------------------------------------------------------------
# Fig. S1 — AMMI biplot (43 complete genotypes)
# ---------------------------------------------------------------------------
ammi_plot_data <- do.call(rbind, lapply(TRAITS, function(t) {
  a <- ammi[[t]]
  g <- data.frame(Type = "Genotype", Label = a$geno,
                  x = a$biplot_geno[, 1], y = a$biplot_geno[, 2])
  e <- data.frame(Type = "Year", Label = a$yrs,
                  x = a$biplot_env[, 1], y = a$biplot_env[, 2])
  rbind(g, e) %>% mutate(Trait = t, PC1pct = ammi_pct[1, t], PC2pct = ammi_pct[2, t])
}))
ammi_plot_data$Trait <- factor(ammi_plot_data$Trait, levels = TRAITS)

pS1 <- ggplot(ammi_plot_data, aes(x = x, y = y, colour = Type, shape = Type)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(size = 1.2) +
  geom_text_repel(aes(label = Label), size = 2.4, max.overlaps = 50, show.legend = FALSE) +
  facet_wrap(~ Trait, scales = "free", ncol = 3,
             labeller = labeller(Trait = setNames(
               sprintf("%s\nIPCA1 %.1f%%  IPCA2 %.1f%%",
                       TRAIT_SHORT[TRAITS], ammi_pct[1, ], ammi_pct[2, ]),
               TRAITS))) +
  scale_colour_manual(values = c(cb[5], cb[6])) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "IPCA1 score", y = "IPCA2 score", colour = NULL, shape = NULL) +
  theme_classic(base_size = 10, base_family = FONT) +
  theme(legend.position = "top")
save_fig(pS1, "FigS1_ammi_biplot", 12, 4.2)

# ---------------------------------------------------------------------------
# Fig. S2 — GGE biplot (43 complete genotypes)
# ---------------------------------------------------------------------------
gge_plot <- do.call(rbind, lapply(TRAITS, function(t) {
  g <- gge[[t]]
  geo <- data.frame(Type = "Genotype", Label = g$geno,
                    x = g$G[, 1], y = g$G[, 2], Trait = t)
  env <- data.frame(Type = "Year", Label = g$yrs,
                    x = g$E[, 1], y = g$E[, 2], Trait = t)
  rbind(geo, env)
}))
gge_plot$Trait <- factor(gge_plot$Trait, levels = TRAITS)

hull_lines <- do.call(rbind, lapply(TRAITS, function(t) {
  g <- gge[[t]]
  pts <- g$G[, 1:2, drop = FALSE]
  h <- chull(pts)
  hp <- pts[h, , drop = FALSE]
  segs <- lapply(seq_along(h), function(i) {
    A <- hp[i, ]; B <- hp[if (i == nrow(hp)) 1 else i + 1, ]
    d <- B - A; denom <- sum(d * d)
    F <- if (denom == 0) A else A + (-sum(A * d) / denom) * d
    data.frame(x = 0, y = 0, xend = F[1], yend = F[2], Trait = t)
  })
  do.call(rbind, segs)
}))
hull_lines$Trait <- factor(hull_lines$Trait, levels = TRAITS)

pS2 <- ggplot() +
  geom_segment(data = hull_lines,
               aes(x = x, y = y, xend = xend, yend = yend),
               linetype = 3, colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(data = gge_plot, aes(x = x, y = y, colour = Type, shape = Type), size = 1.2) +
  geom_text_repel(data = gge_plot,
                  aes(x = x, y = y, label = Label, colour = Type),
                  size = 2.3, max.overlaps = 50, show.legend = FALSE) +
  facet_wrap(~ Trait, scales = "free", ncol = 3,
             labeller = labeller(Trait = setNames(
               sprintf("%s\nPC1 %.1f%%  PC2 %.1f%%",
                       TRAIT_SHORT[TRAITS], gge_pct[1, ], gge_pct[2, ]),
               TRAITS))) +
  scale_colour_manual(values = c(cb[5], cb[6])) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "GGE PC1", y = "GGE PC2", colour = NULL, shape = NULL) +
  theme_classic(base_size = 10, base_family = FONT) +
  theme(legend.position = "top")
save_fig(pS2, "FigS2_gge_biplot", 12, 4.2)

# ---------------------------------------------------------------------------
# Fig. S3 — Mean performance vs. ASV and GGE distance (43 complete genotypes)
# ---------------------------------------------------------------------------
stability_list <- lapply(TRAITS, function(t) {
  a <- ammi[[t]]
  ipca1 <- a$scores[, 1]; ipca2 <- a$scores[, 2]
  ratio <- (a$D[1]^2) / (a$D[2]^2)
  asv <- sqrt((ratio * ipca1)^2 + ipca2^2)

  g <- gge[[t]]
  ideal_x <- max(g$G[, 1], na.rm = TRUE)
  gge_dist <- sqrt((g$G[, 1] - ideal_x)^2 + g$G[, 2]^2)

  M <- gxe(t)
  means <- rowMeans(M[, setdiff(names(M), "Genotype"), drop = FALSE])

  data.frame(
    Genotype = M$Genotype, Mean = means,
    ASV = as.numeric(asv), GGE_Distance = as.numeric(gge_dist),
    Trait = t, stringsAsFactors = FALSE
  )
})

stability_df <- do.call(rbind, stability_list) %>%
  mutate(Trait = factor(Trait, levels = TRAITS)) %>%
  pivot_longer(c(ASV, GGE_Distance), names_to = "Metric", values_to = "StabilityIndex") %>%
  mutate(
    Metric = factor(
      Metric, levels = c("ASV", "GGE_Distance"),
      labels = c("AMMI Stability Value (ASV)", "GGE distance to ideal genotype")
    )
  )

pS3 <- ggplot(stability_df, aes(x = Mean, y = StabilityIndex)) +
  geom_point(size = 1.7, colour = cb[5]) +
  geom_text_repel(aes(label = Genotype), size = 2.1, max.overlaps = 20, colour = "grey30") +
  facet_grid(Trait ~ Metric, scales = "free",
             labeller = labeller(Trait = TRAIT_SHORT)) +
  labs(x = "Mean genotype performance", y = "Stability index") +
  theme_classic(base_size = 10, base_family = FONT) +
  theme(strip.text = element_text(face = "bold"),
        panel.spacing = grid::unit(0.9, "lines"))
save_fig(pS3, "FigS3_mean_stability", 11, 8)

# ---------------------------------------------------------------------------
# Console checks
# ---------------------------------------------------------------------------
cat("Figures saved to:", FIGD, "\n")
cat("Genotypes in main analysis:", nlevels(d86$Genotype), "\n")
cat("Complete genotypes (AMMI/GGE subset):", length(complete43), "\n")
cat("Years:", N_YEAR,
    " Replications/year (harmonic mean):",
    paste(sprintf("%.3f", vc_tab$r_harmonic), collapse = ", "), "\n")
cat("H2 (86 genotypes):",
    paste(vc_tab$Trait, sprintf("%.4f", vc_tab$H2), collapse = "; "), "\n")
cat("Cohort comparison p-values:",
    paste(names(cohort_pvals), format(cohort_pvals, digits = 5), collapse = "; "), "\n")
cat("Generated: Fig1, Fig2, Fig3, Fig4 (86-genotype main) + FigS1, FigS2, FigS3 (43-genotype supplement).\n")
