# ============================================================================
# make_figures.R — Plan A: Figure generation for 43 complete genotypes
# ----------------------------------------------------------------------------
# Generates Fig1–Fig7 in PNG (300 dpi) and PDF formats.
# All labels and legends are in English, with corrected units.
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# ----------------------------------------------------------------------------
# Project directories
# ----------------------------------------------------------------------------
# All paths are relative to the project root.
# Run this script with the project root as the working directory.
ROOT <- "."

RAW <- file.path(
  ROOT,
  "01_Raw_Phenotype_Data"
)

FIGD <- file.path(
  ROOT,
  "04_Results",
  "Genotype43",
  "Figures"
)

dir.create(
  FIGD,
  showWarnings = FALSE,
  recursive = TRUE
)

TRAITS <- c(
  "SummerHeight",
  "FreshWeight",
  "DryWeight"
)

TRAIT_LABEL <- c(
  SummerHeight = "Summer plant height (cm)",
  FreshWeight  = "Fresh weight (kg 5 m²)",
  DryWeight    = "Dry weight (kg 5 m²)"
)

TRAIT_SHORT <- c(
  SummerHeight = "Summer height",
  FreshWeight  = "Fresh weight",
  DryWeight    = "Dry weight"
)

Y <- 4

# Colorblind-friendly palette (Okabe–Ito)
cb <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7",
  "#999999"
)

# ----------------------------------------------------------------------------
# Figure export function
# ----------------------------------------------------------------------------
save_fig <- function(p, name, w, h) {

  ggsave(
    file.path(
      FIGD,
      paste0(name, ".png")
    ),
    p,
    width = w,
    height = h,
    units = "in",
    dpi = 300
  )

  ggsave(
    file.path(
      FIGD,
      paste0(name, ".pdf")
    ),
    p,
    width = w,
    height = h,
    units = "in"
  )
}

# ----------------------------------------------------------------------------
# Data preparation
# Same data-processing logic as analysis_main.R
# ----------------------------------------------------------------------------

read_year <- function(y) {

  f <- file.path(
    RAW,
    paste0(
      "Sino-Australian Alfalfa Project Data - ",
      y,
      "data.csv"
    )
  )

  d <- read.csv(
    f,
    fileEncoding = "UTF-8-BOM",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  d$Year <- y

  d
}

raw <- do.call(
  rbind,
  lapply(
    2002:2005,
    read_year
  )
)

num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

raw$Genotype <- trimws(raw$Line)
raw$Rep      <- as.integer(raw$Rep)

raw$SummerHeight <- num(raw$Mean_Summer)
raw$FreshWeight  <- num(raw$FW)
raw$DryWeight    <- num(raw$DW)

# ----------------------------------------------------------------------------
# Conditional FW/DW swap in 2002
# Only rows with DW > FW are swapped.
# ----------------------------------------------------------------------------

swap <- raw$Year == 2002 &
        !is.na(raw$FreshWeight) &
        !is.na(raw$DryWeight) &
        raw$DryWeight > raw$FreshWeight

fw_old <- raw$FreshWeight[swap]

raw$FreshWeight[swap] <- raw$DryWeight[swap]
raw$DryWeight[swap]   <- fw_old

# ----------------------------------------------------------------------------
# Establishment year
# Plots <= 129 are assigned to 2001; all others to 2002.
# ----------------------------------------------------------------------------

raw$Estab <- ifelse(
  raw$Plot <= 129,
  2001,
  2002
)

# ----------------------------------------------------------------------------
# Aggregate to genotype × year × replicate × establishment-year level
# ----------------------------------------------------------------------------

core_rep <- raw %>%
  group_by(
    Year,
    Genotype,
    Rep,
    Estab
  ) %>%
  summarise(
    across(
      all_of(TRAITS),
      ~ if (all(is.na(.))) {
          NA_real_
        } else {
          mean(., na.rm = TRUE)
        }
    ),
    .groups = "drop"
  )

# ----------------------------------------------------------------------------
# Identify the 43 complete genotypes
# ----------------------------------------------------------------------------

complete43 <- core_rep %>%
  group_by(
    Genotype,
    Year
  ) %>%
  summarise(
    across(
      all_of(TRAITS),
      ~ any(!is.na(.))
    ),
    .groups = "drop"
  ) %>%
  group_by(Genotype) %>%
  summarise(
    n_complete = sum(
      SummerHeight &
        FreshWeight &
        DryWeight
    ),
    .groups = "drop"
  ) %>%
  filter(
    n_complete == Y
  ) %>%
  pull(Genotype)

# ----------------------------------------------------------------------------
# Missing-data pattern
# 86 genotypes × 4 years
# ----------------------------------------------------------------------------

missing_df <- core_rep %>%
  group_by(
    Genotype,
    Year
  ) %>%
  summarise(
    n_trait =
      (any(!is.na(SummerHeight))) +
      (any(!is.na(FreshWeight))) +
      (any(!is.na(DryWeight))),
    .groups = "drop"
  ) %>%
  mutate(
    Status = case_when(
      n_trait == 3 ~ "Complete",
      n_trait == 0 ~ "Absent",
      TRUE ~ "Partial"
    ),

    Status = factor(
      Status,
      levels = c(
        "Complete",
        "Partial",
        "Absent"
      )
    )
  ) %>%
  mutate(
    Genotype = factor(Genotype),

    Year = factor(
      Year,
      levels = as.character(2002:2005)
    )
  )

# ----------------------------------------------------------------------------
# Genotype ordering:
# descending number of complete years, then alphabetical order
# ----------------------------------------------------------------------------

geno_order <- missing_df %>%
  group_by(Genotype) %>%
  summarise(
    n_comp = sum(
      n_trait == 3
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_comp),
    Genotype
  ) %>%
  pull(Genotype)

missing_df$Genotype <- factor(
  missing_df$Genotype,
  levels = geno_order
)

# ----------------------------------------------------------------------------
# Genotype means for complete vs incomplete groups
# Fair comparison window: 2003–2005
# ----------------------------------------------------------------------------

group_label <- function(g) {
  ifelse(
    g %in% complete43,
    "Complete (n = 43)",
    "Incomplete (n = 43)"
  )
}

sub_means <- raw %>%
  filter(
    Year >= 2003
  ) %>%
  group_by(Genotype) %>%
  summarise(
    across(
      all_of(TRAITS),
      ~ mean(., na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Group = group_label(Genotype)
  ) %>%
  pivot_longer(
    all_of(TRAITS),
    names_to = "Trait",
    values_to = "value"
  ) %>%
  mutate(
    Trait = factor(
      Trait,
      levels = TRAITS
    ),

    Group = factor(
      Group,
      levels = c(
        "Complete (n = 43)",
        "Incomplete (n = 43)"
      )
    )
  )

# ----------------------------------------------------------------------------
# Mixed models + BLUPs + variance components + AMMI + GGE
# ----------------------------------------------------------------------------

d43 <- core_rep %>%
  filter(
    Genotype %in% complete43
  ) %>%
  mutate(
    Genotype = factor(Genotype),
    Year     = factor(Year),
    Rep      = factor(Rep)
  )

get_vc <- function(m, grp) {

  vc <- as.data.frame(
    VarCorr(m)
  )

  v <- vc$vcov[
    vc$grp == grp
  ]

  if (length(v)) {
    v[1]
  } else {
    0
  }
}

fit_models <- setNames(
  lapply(
    TRAITS,
    function(t) {

      d <- d43 %>%
        filter(
          !is.na(.data[[t]])
        )

      f <- as.formula(
        paste0(
          t,
          " ~ Year + (1|Genotype) + ",
          "(1|Genotype:Year) + (1|Year:Rep)"
        )
      )

      lmer(
        f,
        data = d,
        REML = TRUE,
        control = lmerControl(
          calc.derivs = FALSE
        )
      )
    }
  ),
  TRAITS
)

blup_list <- lapply(
  TRAITS,
  function(t) {

    m <- fit_models[[t]]

    rr <- ranef(
      m,
      condVar = TRUE
    )

    se <- sqrt(
      attr(
        rr$Genotype,
        "postVar"
      )[1, 1, ]
    )

    data.frame(
      Genotype = rownames(
        rr$Genotype
      ),

      BLUP = as.numeric(
        fixef(m)[1] +
          rr$Genotype[[1]]
      ),

      SE = se,

      stringsAsFactors = FALSE
    )
  }
)

names(blup_list) <- TRAITS

# ----------------------------------------------------------------------------
# AMMI / GGE coordinates
# Base R SVD
# ----------------------------------------------------------------------------

gxe <- function(t) {

  d43 %>%
    filter(
      !is.na(.data[[t]])
    ) %>%
    group_by(
      Genotype,
      Year
    ) %>%
    summarise(
      mu = mean(
        .data[[t]]
      ),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Year,
      values_from = mu
    ) %>%
    arrange(Genotype)
}

ammi <- lapply(
  TRAITS,
  function(t) {

    M <- gxe(t)

    yrs <- setdiff(
      names(M),
      "Genotype"
    )

    X <- as.matrix(
      M[, yrs]
    )

    rownames(X) <- M$Genotype

    gr <- rowMeans(X)
    cg <- colMeans(X)
    gm <- mean(X)

    # Remove additive genotype and environment main effects
    I <- X -
         gr -
         rep(
           cg,
           each = nrow(X)
         ) +
         gm

    s <- svd(I)

    list(
      U    = s$u,
      V    = s$v,
      D    = s$d,
      geno = M$Genotype,
      yrs  = yrs
    )
  }
)

names(ammi) <- TRAITS

ammi_pct <- sapply(
  TRAITS,
  function(t) {

    d <- ammi[[t]]$D^2

    round(
      d / sum(d) * 100,
      2
    )
  }
)

gge <- lapply(
  TRAITS,
  function(t) {

    M <- gxe(t)

    yrs <- setdiff(
      names(M),
      "Genotype"
    )

    X <- as.matrix(
      M[, yrs]
    )

    rownames(X) <- M$Genotype

    # Environment-centered GGE matrix
    Xc <- scale(
      X,
      center = colMeans(X),
      scale = FALSE
    )

    s <- svd(Xc)

    G <- s$u %*%
         diag(
           sqrt(s$d)
         )

    E <- s$v %*%
         diag(
           sqrt(s$d)
         )

    list(
      G    = G,
      E    = E,
      D    = s$d,
      geno = M$Genotype,
      yrs  = yrs
    )
  }
)

names(gge) <- TRAITS

gge_pct <- sapply(
  TRAITS,
  function(t) {

    d <- gge[[t]]$D^2

    round(
      d / sum(d) * 100,
      2
    )
  }
)

# ============================================================================
# Fig1 — Missing-data pattern heatmap (86 × 4)
# ============================================================================

p1 <- ggplot(
  missing_df,
  aes(
    x = Year,
    y = Genotype,
    fill = Status
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.1
  ) +

  scale_fill_manual(
    values = c(
      "#009E73",
      "#F0E442",
      "#D9D9D9"
    ),
    name = "Observation",
    drop = FALSE
  ) +

  scale_x_discrete(
    expand = c(0, 0)
  ) +

  labs(
    x = "Year",
    y = "Genotype",
    title = "Data availability across 86 genotypes and four years"
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    axis.text.y = element_text(
      size = 5.5
    ),

    legend.position = "top",

    plot.title = element_text(
      face = "bold",
      size = 11
    )
  )

save_fig(
  p1,
  "Fig1_missing_pattern",
  5.5,
  9.5
)

# ============================================================================
# Fig2 — Complete vs incomplete subset comparison
#         Boxplots + significance
# ============================================================================

# Welch's t-test P values for 2003–2005
pvals <- sapply(
  TRAITS,
  function(t) {

    a <- sub_means$value[
      sub_means$Trait == t &
        sub_means$Group == "Complete (n = 43)"
    ]

    b <- sub_means$value[
      sub_means$Trait == t &
        sub_means$Group == "Incomplete (n = 43)"
    ]

    t.test(
      a,
      b
    )$p.value
  }
)

sig_labels <- setNames(
  ifelse(
    pvals < 0.001,
    "***",
    ifelse(
      pvals < 0.01,
      "**",
      ifelse(
        pvals < 0.05,
        "*",
        "ns"
      )
    )
  ),
  TRAITS
)

p2 <- ggplot(
  sub_means,
  aes(
    x = Group,
    y = value,
    fill = Group
  )
) +

  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.6
  ) +

  geom_jitter(
    width = 0.12,
    size = 0.7,
    alpha = 0.4
  ) +

  facet_wrap(
    ~ Trait,
    scales = "free_y",
    labeller = labeller(
      Trait = TRAIT_LABEL
    )
  ) +

  scale_fill_manual(
    values = c(
      "#0072B2",
      "#D55E00"
    ),
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = "Trait value (genotypic mean, 2003–2005)",
    title = "Complete vs incomplete genotype subsets"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    ),

    strip.text = element_text(
      face = "bold"
    )
  )

# Add significance labels within each facet
sig_df <- data.frame(
  Trait = factor(
    TRAITS,
    levels = TRAITS
  ),

  y = sapply(
    TRAITS,
    function(t) {

      max(
        sub_means$value[
          sub_means$Trait == t
        ],
        na.rm = TRUE
      ) * 1.08
    }
  ),

  label = sig_labels[TRAITS],

  stringsAsFactors = FALSE
)

p2 <- p2 +
  geom_text(
    data = sig_df,
    aes(
      x = 1.5,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  )

save_fig(
  p2,
  "Fig2_subset_comparison",
  10,
  4
)

# ============================================================================
# Fig3 — Variance-component proportions
#         Stacked bar chart
# ============================================================================

vc_tab <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      m <- fit_models[[t]]

      data.frame(
        Trait = t,

        V_G = get_vc(
          m,
          "Genotype"
        ),

        V_GY = get_vc(
          m,
          "Genotype:Year"
        ),

        V_RY = get_vc(
          m,
          "Year:Rep"
        ),

        V_e = sigma(m)^2,

        stringsAsFactors = FALSE
      )
    }
  )
)

vc_long <- vc_tab %>%

  mutate(
    Trait = factor(
      Trait,
      levels = TRAITS
    )
  ) %>%

  pivot_longer(
    c(
      V_G,
      V_GY,
      V_RY,
      V_e
    ),
    names_to = "Component",
    values_to = "value"
  ) %>%

  group_by(Trait) %>%
  mutate(
    prop = value / sum(value) * 100
  ) %>%

  mutate(
    Component = factor(
      Component,
      levels = c(
        "V_G",
        "V_GY",
        "V_RY",
        "V_e"
      ),

      labels = c(
        "V[G] (genotype)",
        "V[GY] (G × Year)",
        "V[RY] (block)",
        "V[e] (residual)"
      )
    )
  )

p3 <- ggplot(
  vc_long,
  aes(
    x = Trait,
    y = prop,
    fill = Component
  )
) +

  geom_col(
    width = 0.55
  ) +

  scale_fill_manual(
    values = c(
      "#0072B2",
      "#56B4E9",
      "#E69F00",
      "#D9D9D9"
    ),
    name = "Variance component",

    labels = c(
      expression(V[G]),
      expression(V[GY]),
      expression(V[RY]),
      expression(V[e])
    )
  ) +

  labs(
    x = "Trait",
    y = "Proportion of total variance (%)",
    title = "Partitioning of phenotypic variance"
  ) +

  scale_x_discrete(
    labels = c(
      "Summer height",
      "Fresh weight",
      "Dry weight"
    )
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    ),

    legend.position = "right"
  )

save_fig(
  p3,
  "Fig3_variance_components",
  7.5,
  4
)

# ============================================================================
# Fig4 — BLUP ranking with 95% CI error bars
#         Top five genotypes highlighted
# ============================================================================

blup_rank <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      b <- blup_list[[t]]

      b$Trait <- t

      b %>%
        arrange(
          desc(BLUP)
        ) %>%
        mutate(
          rank = 1:n(),
          ci_lo = BLUP - 1.96 * SE,
          ci_hi = BLUP + 1.96 * SE
        )
    }
  )
)

blup_rank$Trait <- factor(
  blup_rank$Trait,
  levels = TRAITS
)

blup_rank$top5 <- blup_rank$rank <= 5

p4 <- ggplot(
  blup_rank,
  aes(
    x = rank,
    y = BLUP,
    colour = top5
  )
) +

  geom_errorbar(
    aes(
      ymin = ci_lo,
      ymax = ci_hi
    ),
    width = 0.3,
    linewidth = 0.3,
    alpha = 0.7
  ) +

  geom_point(
    size = 0.9
  ) +

  geom_text_repel(
    data = subset(
      blup_rank,
      top5
    ),

    aes(
      label = Genotype
    ),

    size = 2.6,
    colour = "#D55E00",
    nudge_x = 1.2,
    direction = "y",
    segment.size = 0.2,
    max.overlaps = 20
  ) +

  facet_wrap(
    ~ Trait,
    scales = "free_y",
    ncol = 1,
    labeller = labeller(
      Trait = TRAIT_LABEL
    )
  ) +

  scale_colour_manual(
    values = c(
      "#999999",
      "#D55E00"
    ),
    guide = "none"
  ) +

  labs(
    x = "Rank (sorted by BLUP)",
    y = "BLUP genotypic value",
    title = "Genotypic BLUPs with 95% confidence intervals"
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    ),

    strip.text = element_text(
      face = "bold",
      size = 9
    )
  )

save_fig(
  p4,
  "Fig4_blup_ranking",
  7,
  8
)

# ============================================================================
# Fig5 — AMMI biplot
#         IPCA1 vs IPCA2
# ============================================================================

ammi_plot_data <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      a <- ammi[[t]]

      g <- data.frame(
        Type = "Genotype",
        Label = a$geno,
        x = a$U[, 1],
        y = a$U[, 2]
      )

      e <- data.frame(
        Type = "Year",
        Label = a$yrs,
        x = a$V[, 1],
        y = a$V[, 2]
      )

      rbind(
        g,
        e
      ) %>%
        mutate(
          Trait = t,

          xlab = sprintf(
            "IPCA1 (%.1f%%)",
            ammi_pct[1, t]
          ),

          ylab = sprintf(
            "IPCA2 (%.1f%%)",
            ammi_pct[2, t]
          )
        )
    }
  )
)

ammi_plot_data$Trait <- factor(
  ammi_plot_data$Trait,
  levels = TRAITS
)

p5 <- ggplot(
  ammi_plot_data,
  aes(
    x = x,
    y = y,
    colour = Type,
    shape = Type
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +

  geom_point(
    size = 1.2
  ) +

  geom_text_repel(
    aes(
      label = Label
    ),
    size = 2.4,
    max.overlaps = 50,
    show.legend = FALSE
  ) +

  facet_wrap(
    ~ Trait,
    scales = "free",
    ncol = 3,

    labeller = labeller(
      Trait = setNames(
        sprintf(
          "%s\nIPCA1 %.1f%%  IPCA2 %.1f%%",
          TRAIT_SHORT[TRAITS],
          ammi_pct[1, ],
          ammi_pct[2, ]
        ),
        TRAITS
      )
    )
  ) +

  scale_colour_manual(
    values = c(
      "#0072B2",
      "#D55E00"
    )
  ) +

  scale_shape_manual(
    values = c(
      16,
      17
    )
  ) +

  labs(
    x = "IPCA1",
    y = "IPCA2",
    colour = NULL,
    shape = NULL,
    title = "AMMI biplots (genotypes vs years)"
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    ),

    legend.position = "top"
  )

save_fig(
  p5,
  "Fig5_ammi_biplot",
  12,
  4.2
)

# ============================================================================
# Fig6 — GGE biplot
#         Which-won-where
# ============================================================================

gge_plot <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      g <- gge[[t]]

      geo <- data.frame(
        Type = "Genotype",
        Label = g$geno,
        x = g$G[, 1],
        y = g$G[, 2],
        Trait = t
      )

      env <- data.frame(
        Type = "Year",
        Label = g$yrs,
        x = g$E[, 1],
        y = g$E[, 2],
        Trait = t
      )

      rbind(
        geo,
        env
      )
    }
  )
)

gge_plot$Trait <- factor(
  gge_plot$Trait,
  levels = TRAITS
)

# ----------------------------------------------------------------------------
# Convex hulls + perpendicular projections from the origin
# ----------------------------------------------------------------------------

hull_lines <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      g <- gge[[t]]

      pts <- g$G[, 1:2]

      h <- chull(
        pts
      )

      hp <- pts[h, ]

      segs <- lapply(
        seq_along(h),
        function(i) {

          A <- hp[i, ]

          B <- hp[
            if (i == nrow(hp)) {
              1
            } else {
              i + 1
            },
          ]

          d <- B - A

          tt <- -sum(
            A * d
          ) / sum(
            d * d
          )

          F <- A + tt * d

          data.frame(
            x = 0,
            y = 0,
            xend = F[1],
            yend = F[2],
            Trait = t
          )
        }
      )

      do.call(
        rbind,
        segs
      )
    }
  )
)

hull_lines$Trait <- factor(
  hull_lines$Trait,
  levels = TRAITS
)

p6 <- ggplot() +

  geom_segment(
    data = hull_lines,

    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),

    linetype = 3,
    colour = "grey50",
    linewidth = 0.4
  ) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +

  geom_point(
    data = gge_plot,

    aes(
      x = x,
      y = y,
      colour = Type,
      shape = Type
    ),

    size = 1.2
  ) +

  geom_text_repel(
    data = gge_plot,

    aes(
      x = x,
      y = y,
      label = Label,
      colour = Type
    ),

    size = 2.3,
    max.overlaps = 50,
    show.legend = FALSE
  ) +

  facet_wrap(
    ~ Trait,
    scales = "free",
    ncol = 3,

    labeller = labeller(
      Trait = setNames(
        sprintf(
          "%s\nPC1 %.1f%%  PC2 %.1f%%",
          TRAIT_SHORT[TRAITS],
          gge_pct[1, ],
          gge_pct[2, ]
        ),
        TRAITS
      )
    )
  ) +

  scale_colour_manual(
    values = c(
      "#0072B2",
      "#D55E00"
    )
  ) +

  scale_shape_manual(
    values = c(
      16,
      17
    )
  ) +

  labs(
    x = "PC1",
    y = "PC2",
    colour = NULL,
    shape = NULL,
    title = "GGE biplots (which-won-where)"
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    ),

    legend.position = "top"
  )

save_fig(
  p6,
  "Fig6_gge_biplot",
  12,
  4.2
)

# ============================================================================
# Fig7 — Summer height vs fresh weight BLUP scatter plot
# ============================================================================

sh <- blup_list[["SummerHeight"]]
fw <- blup_list[["FreshWeight"]]

sc <- data.frame(
  Genotype = sh$Genotype,
  Height   = sh$BLUP,
  FW       = fw$BLUP
)

ct <- cor.test(
  sc$Height,
  sc$FW
)

p7 <- ggplot(
  sc,
  aes(
    x = Height,
    y = FW
  )
) +

  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "#0072B2",
    fill = "#56B4E9",
    alpha = 0.25,
    linewidth = 0.8
  ) +

  geom_point(
    size = 2,
    colour = "#0072B2"
  ) +

  geom_text_repel(
    aes(
      label = Genotype
    ),
    size = 2.2,
    max.overlaps = 15
  ) +

  annotate(
    "text",
    x = min(sc$Height) + 2,
    y = max(sc$FW) - 0.5,

    label = sprintf(
      "r = %.3f\np = %.3g",
      ct$estimate,
      ct$p.value
    ),

    hjust = 0,
    size = 3.5
  ) +

  labs(
    x = "Summer height BLUP (cm)",
    y = "Fresh weight BLUP (kg 5 m²)",
    title = "Correlation between summer height and fresh weight BLUPs"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

save_fig(
  p7,
  "Fig7_blup_scatter",
  6.5,
  6
)

# ----------------------------------------------------------------------------
# Final console summary
# ----------------------------------------------------------------------------

cat(
  "All figures saved to:",
  FIGD,
  "\n"
)

cat(
  "P-values for subset comparison:",
  paste(
    names(pvals),
    round(
      pvals,
      6
    ),
    collapse = "; "
  ),
  "\n"
)
