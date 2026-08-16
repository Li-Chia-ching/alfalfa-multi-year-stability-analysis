# ============================================================================
# analysis_main.R — Plan A: Main analysis based on 43 genotypes with complete
#                    four-year records
# ----------------------------------------------------------------------------
# Data:
#   Four CSV files from 2002–2005 in 01_Raw_Phenotype_Data/
#   Only the 43 genotypes with complete records for all four years are used.
#
# Methods:
#   lme4 linear mixed models
#   Year is treated as a fixed effect.
#   Genotype, Genotype:Year, and Year:Rep are treated as random effects.
#
#   AMMI / GGE analyses use base R svd(), consistent with the 03_Analysis
#   pipeline and without dependencies on agricolae or metan.
#
# Notes:
#   1) Original columns are mapped as:
#      Line -> Genotype
#      Mean_Summer -> SummerHeight
#      FW -> FreshWeight
#      DW -> DryWeight
#      Year is extracted from the file name.
#
#   2) Text markers such as #DIV/0!, "-", death, and missing plants are
#      converted to NA by as.numeric().
#
#   3) The 2002 FW/DW columns contain anomalies. Based on the physical
#      constraint FW >= DW, rows with DW > FW are conditionally swapped.
#
#   4) Raw observations are aggregated to the
#      genotype × year × replicate experimental-unit level
#      to avoid pseudoreplication.
#
#   5) H2 is interpreted only as broad-sense heritability and is not described
#      as additive genetic variance.
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
})

# ----------------------------------------------------------------------------
# Project directories
# ----------------------------------------------------------------------------
# All paths are relative to the project working directory.
# The script should be run from the project root directory.
ROOT <- "."
RAW  <- file.path(ROOT, "01_Raw_Phenotype_Data")
OUT  <- file.path(ROOT, "04_Results", "Genotype43")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TRAITS     <- c("SummerHeight", "FreshWeight", "DryWeight")
TRAIT_UNIT <- c(
  SummerHeight = "cm",
  FreshWeight  = "kg/5m2",
  DryWeight    = "kg/5m2"
)

Y <- 4  # Number of years (environments)

# ----------------------------------------------------------------------------
# 1. Read the four annual CSV files
# ----------------------------------------------------------------------------
read_year <- function(y) {
  f <- file.path(
    RAW,
    paste0("Sino-Australian Alfalfa Project Data - ", y, "data.csv")
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

raw <- do.call(rbind, lapply(2002:2005, read_year))

# ----------------------------------------------------------------------------
# 2. Column mapping and data cleaning
# ----------------------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))

raw$Genotype     <- trimws(raw$Line)
raw$Rep          <- as.integer(raw$Rep)
raw$SummerHeight <- num(raw$Mean_Summer)
raw$FreshWeight  <- num(raw$FW)
raw$DryWeight    <- num(raw$DW)

# Conditional FW/DW swap in 2002
# Physical constraint: FreshWeight >= DryWeight
swap <- raw$Year == 2002 &
        !is.na(raw$FreshWeight) &
        !is.na(raw$DryWeight) &
        raw$DryWeight > raw$FreshWeight

fw_old <- raw$FreshWeight[swap]

raw$FreshWeight[swap] <- raw$DryWeight[swap]
raw$DryWeight[swap]   <- fw_old

cat(sprintf(
  "Number of rows with conditional FW/DW swaps in 2002 = %d\n",
  sum(swap)
))

# ----------------------------------------------------------------------------
# 3. Aggregate to genotype × year × replicate experimental units
# ----------------------------------------------------------------------------
core <- raw[, c("Year", "Genotype", "Rep", TRAITS)]

core_rep <- core %>%
  group_by(Year, Genotype, Rep) %>%
  summarise(
    across(
      all_of(TRAITS),
      ~ if (all(is.na(.))) NA_real_ else mean(., na.rm = TRUE)
    ),
    .groups = "drop"
  )

# ----------------------------------------------------------------------------
# 4. Identify the 43 complete genotypes
# ----------------------------------------------------------------------------
# A genotype is considered complete when, for each of the four years,
# all three traits have at least one non-missing observation.
# ----------------------------------------------------------------------------
complete43 <- core_rep %>%
  group_by(Genotype, Year) %>%
  summarise(
    across(all_of(TRAITS), ~ any(!is.na(.))),
    .groups = "drop"
  ) %>%
  group_by(Genotype) %>%
  summarise(
    n_complete = sum(
      SummerHeight & FreshWeight & DryWeight
    ),
    .groups = "drop"
  ) %>%
  filter(n_complete == Y) %>%
  pull(Genotype)

cat(sprintf(
  "Number of complete genotypes = %d (expected: 43)\n",
  length(complete43)
))

# Filter to the 43 complete genotypes
d43 <- core_rep %>%
  filter(Genotype %in% complete43) %>%
  mutate(
    Genotype = factor(Genotype),
    Year     = factor(Year),
    Rep      = factor(Rep)
  )

cat(sprintf(
  "Filtered data rows = %d; genotypes = %d; years = %s\n",
  nrow(d43),
  nlevels(d43$Genotype),
  paste(levels(d43$Year), collapse = ",")
))

# ----------------------------------------------------------------------------
# 5. Mixed models, variance components, and broad-sense heritability
# ----------------------------------------------------------------------------
get_vc <- function(m, grp) {
  vc  <- as.data.frame(VarCorr(m))
  val <- vc$vcov[vc$grp == grp]

  if (length(val) == 0) 0 else val[1]
}

# Harmonic mean number of replicates for each trait
harm_mean_rep <- function(trait) {
  k <- d43 %>%
    filter(!is.na(.data[[trait]])) %>%
    group_by(Genotype, Year) %>%
    summarise(n = n(), .groups = "drop") %>%
    pull(n)

  1 / mean(1 / k)
}

fit_models <- setNames(
  lapply(TRAITS, function(t) {

    d <- d43 %>%
      filter(!is.na(.data[[t]]))

    f <- as.formula(
      paste0(
        t,
        " ~ Year + (1|Genotype) + (1|Genotype:Year) + (1|Year:Rep)"
      )
    )

    lmer(
      f,
      data = d,
      REML = TRUE,
      control = lmerControl(calc.derivs = FALSE)
    )
  }),
  TRAITS
)

# Variance component table
vc_tab <- do.call(
  rbind,
  lapply(TRAITS, function(t) {

    m <- fit_models[[t]]

    data.frame(
      Trait = t,
      V_G   = get_vc(m, "Genotype"),
      V_GY  = get_vc(m, "Genotype:Year"),
      V_RY  = get_vc(m, "Year:Rep"),
      V_e   = sigma(m)^2,
      stringsAsFactors = FALSE
    )
  })
)

print(vc_tab)

# Broad-sense heritability point estimate and 95% CI
# Parametric bootstrap using bootMer(), nsim = 1000
h2_point <- function(m, t) {

  V_G  <- get_vc(m, "Genotype")
  V_GY <- get_vc(m, "Genotype:Year")
  V_e  <- sigma(m)^2
  r    <- harm_mean_rep(t)

  V_G / (
    V_G +
      V_GY / Y +
      V_e / (Y * r)
  )
}

h2_fun_factory <- function(t) {

  r <- harm_mean_rep(t)

  function(m) {

    V_G  <- get_vc(m, "Genotype")
    V_GY <- get_vc(m, "Genotype:Year")
    V_e  <- sigma(m)^2

    V_G / (
      V_G +
        V_GY / Y +
        V_e / (Y * r)
    )
  }
}

set.seed(123)

h2_res <- lapply(TRAITS, function(t) {

  m <- fit_models[[t]]
  h2 <- h2_point(m, t)

  b <- tryCatch(
    bootMer(
      m,
      h2_fun_factory(t),
      nsim = 1000,
      seed = 123,
      .progress = "none"
    ),
    error = function(e) NULL
  )

  ci <- if (!is.null(b)) {
    quantile(
      b$t,
      c(0.025, 0.975),
      na.rm = TRUE
    )
  } else {
    c(NA, NA)
  }

  data.frame(
    Trait       = t,
    H2          = h2,
    H2_lo       = ci[1],
    H2_hi       = ci[2],
    r_harmonic  = harm_mean_rep(t),
    stringsAsFactors = FALSE
  )
})

h2_tab <- do.call(rbind, h2_res)

h2_tab$H2         <- round(h2_tab$H2, 3)
h2_tab$H2_lo      <- round(h2_tab$H2_lo, 3)
h2_tab$H2_hi      <- round(h2_tab$H2_hi, 3)

print(h2_tab)

# ----------------------------------------------------------------------------
# 6. BLUP extraction, standard errors, and ranking
# ----------------------------------------------------------------------------
blup_list <- lapply(TRAITS, function(t) {

  m <- fit_models[[t]]

  rr <- ranef(
    m,
    condVar = TRUE
  )

  # Random-effect deviation (BLUP)
  g <- rr$Genotype[[1]]

  # Conditional posterior variance
  # 1 × 1 × n array
  pv <- attr(
    rr$Genotype,
    "postVar"
  )

  se <- sqrt(pv[1, 1, ])

  data.frame(
    Genotype = rownames(rr$Genotype),
    BLUP_dev = as.numeric(g),
    BLUP     = as.numeric(fixef(m)[1] + g),
    SE       = se,
    stringsAsFactors = FALSE
  )
})

names(blup_list) <- TRAITS

# Merge the three trait BLUP tables.
# The 43 genotypes are aligned by the same factor levels.
blup_out <- data.frame(
  Genotype = blup_list[[1]]$Genotype,

  SH_BLUP = blup_list[[1]]$BLUP,
  SH_SE   = blup_list[[1]]$SE,

  FW_BLUP = blup_list[[2]]$BLUP,
  FW_SE   = blup_list[[2]]$SE,

  DW_BLUP = blup_list[[3]]$BLUP,
  DW_SE   = blup_list[[3]]$SE,

  stringsAsFactors = FALSE
)

write.csv(
  blup_out,
  file.path(OUT, "BLUP_43.csv"),
  row.names = FALSE
)

# Top 10 genotypes for each trait
top10 <- lapply(TRAITS, function(t) {

  blup_list[[t]] %>%
    arrange(desc(BLUP)) %>%
    slice_head(n = 10) %>%
    mutate(
      across(
        c(BLUP, BLUP_dev, SE),
        ~ round(., 3)
      )
    )
})

names(top10) <- TRAITS

# ----------------------------------------------------------------------------
# 7. AMMI analysis (base R SVD, 43 complete genotypes only)
# ----------------------------------------------------------------------------
ammi_one <- function(t) {

  M <- d43 %>%
    filter(!is.na(.data[[t]])) %>%
    group_by(Genotype, Year) %>%
    summarise(
      mu = mean(.data[[t]]),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Year,
      values_from = mu
    ) %>%
    arrange(Genotype)

  yrs <- setdiff(
    names(M),
    "Genotype"
  )

  X <- as.matrix(M[, yrs])
  rownames(X) <- M$Genotype

  gr <- rowMeans(X)
  cg <- colMeans(X)
  gm <- mean(X)

  # Interaction matrix after removing additive genotype and environment
  # main effects
  I <- X -
       gr -
       rep(cg, each = nrow(X)) +
       gm

  s <- svd(I)

  list(
    geno  = M$Genotype,
    years = yrs,
    D     = s$d,
    U     = s$u,
    V     = s$v
  )
}

ammi_res <- setNames(
  lapply(TRAITS, ammi_one),
  TRAITS
)

ammi_tab <- do.call(
  rbind,
  lapply(TRAITS, function(t) {

    a <- ammi_res[[t]]

    ev  <- a$D^2
    tot <- sum(ev)

    data.frame(
      Trait     = t,
      PC1_pct   = round(ev[1] / tot * 100, 2),
      PC2_pct   = round(ev[2] / tot * 100, 2),
      PC1_cum   = round(ev[1] / tot * 100, 2),
      PC12_cum  = round((ev[1] + ev[2]) / tot * 100, 2),
      stringsAsFactors = FALSE
    )
  })
)

print(ammi_tab)

# Export AMMI biplot data (IPCA coordinates)
for (t in TRAITS) {

  a <- ammi_res[[t]]
  k <- min(2, length(a$D))

  gs <- data.frame(
    Genotype = a$geno,
    a$U[, 1:k, drop = FALSE]
  )

  colnames(gs)[-1] <- paste0(
    "IPCA",
    1:k
  )

  es <- data.frame(
    Year = a$years,
    a$V[, 1:k, drop = FALSE]
  )

  colnames(es)[-1] <- paste0(
    "IPCA",
    1:k
  )

  write.csv(
    gs,
    file.path(
      OUT,
      paste0(
        "AMMI_genotype_IPCA_",
        t,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    es,
    file.path(
      OUT,
      paste0(
        "AMMI_env_IPCA_",
        t,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

# ----------------------------------------------------------------------------
# 8. GGE analysis
#    Base R SVD with environment centering and symmetric scaling
# ----------------------------------------------------------------------------
gge_one <- function(t) {

  M <- d43 %>%
    filter(!is.na(.data[[t]])) %>%
    group_by(Genotype, Year) %>%
    summarise(
      mu = mean(.data[[t]]),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Year,
      values_from = mu
    ) %>%
    arrange(Genotype)

  yrs <- setdiff(
    names(M),
    "Genotype"
  )

  X <- as.matrix(M[, yrs])
  rownames(X) <- M$Genotype

  # Environment centering
  Xc <- scale(
    X,
    center = colMeans(X),
    scale = FALSE
  )

  s <- svd(Xc)

  list(
    D     = s$d,
    U     = s$u,
    V     = s$v,
    geno  = M$Genotype,
    years = yrs
  )
}

gge_res <- setNames(
  lapply(TRAITS, gge_one),
  TRAITS
)

gge_tab <- do.call(
  rbind,
  lapply(TRAITS, function(t) {

    a <- gge_res[[t]]

    ev  <- a$D^2
    tot <- sum(ev)

    data.frame(
      Trait     = t,
      PC1_pct   = round(ev[1] / tot * 100, 2),
      PC2_pct   = round(ev[2] / tot * 100, 2),
      PC12_cum  = round((ev[1] + ev[2]) / tot * 100, 2),
      stringsAsFactors = FALSE
    )
  })
)

print(gge_tab)

# ----------------------------------------------------------------------------
# 9. Correlation between SummerHeight and FreshWeight BLUPs
# ----------------------------------------------------------------------------
sh <- blup_list[["SummerHeight"]]$BLUP
fw <- blup_list[["FreshWeight"]]$BLUP

ct <- cor.test(
  sh,
  fw,
  method = "pearson"
)

cor_res <- data.frame(
  r  = round(ct$estimate, 4),
  p  = signif(ct$p.value, 4),
  lo = round(ct$conf.int[1], 4),
  hi = round(ct$conf.int[2], 4),
  n  = length(sh),
  stringsAsFactors = FALSE
)

print(cor_res)

# ----------------------------------------------------------------------------
# 10. Write results_summary.md
# ----------------------------------------------------------------------------
lines <- c(

  "# Analysis Results Summary — Plan A (43 Complete Genotypes)",

  "",

  sprintf(
    "**Generated on:** %s",
    Sys.Date()
  ),

  "**Method:** lme4 linear mixed model `trait ~ Year + (1|Genotype) + (1|Genotype:Year) + (1|Year:Rep)`",

  "**Data:** Four CSV files from 2002–2005; only the 43 genotypes with complete records for all four years and all three traits (`complete43`) were used.",

  "",

  "---",

  "",

  "## 1. Variance Components",

  "",

  "| Trait | V_G | V_GY | V_RY (Year:Rep) | V_e |",
  "|---|---|---|---|---|",

  sprintf(
    "| SummerHeight | %.3f | %.3f | %.3f | %.3f |",
    vc_tab$V_G[1],
    vc_tab$V_GY[1],
    vc_tab$V_RY[1],
    vc_tab$V_e[1]
  ),

  sprintf(
    "| FreshWeight | %.3f | %.3f | %.3f | %.3f |",
    vc_tab$V_G[2],
    vc_tab$V_GY[2],
    vc_tab$V_RY[2],
    vc_tab$V_e[2]
  ),

  sprintf(
    "| DryWeight | %.3f | %.3f | %.3f | %.3f |",
    vc_tab$V_G[3],
    vc_tab$V_GY[3],
    vc_tab$V_RY[3],
    vc_tab$V_e[3]
  ),

  "",

  "## 2. Broad-Sense Heritability H²",
  "",
  "95% confidence intervals were obtained using 1,000 parametric bootstrap replicates.",

  "",

  "> H² = V_G / (V_G + V_GY/Y + V_e/(Y·r)), where Y = 4 and r = harmonic mean number of replicates.",

  "",

  "| Trait | H² | 95% CI Lower | 95% CI Upper | r (Harmonic Mean) |",
  "|---|---|---|---|---|",

  sprintf(
    "| SummerHeight | %.3f | %.3f | %.3f | %.3f |",
    h2_tab$H2[1],
    h2_tab$H2_lo[1],
    h2_tab$H2_hi[1],
    h2_tab$r_harmonic[1]
  ),

  sprintf(
    "| FreshWeight | %.3f | %.3f | %.3f | %.3f |",
    h2_tab$H2[2],
    h2_tab$H2_lo[2],
    h2_tab$H2_hi[2],
    h2_tab$r_harmonic[2]
  ),

  sprintf(
    "| DryWeight | %.3f | %.3f | %.3f | %.3f |",
    h2_tab$H2[3],
    h2_tab$H2_lo[3],
    h2_tab$H2_hi[3],
    h2_tab$r_harmonic[3]
  ),

  "",

  "## 3. BLUP Ranking (Top 10 Genotypes)",
  "",
  "Genotypic value = model intercept + genotype random-effect deviation.",
  ""
)

for (t in TRAITS) {

  lines <- c(
    lines,

    sprintf(
      "### %s (Unit: %s)",
      t,
      TRAIT_UNIT[t]
    ),

    "",

    "| Rank | Genotype | BLUP | SE |",
    "|---|---|---|---|"
  )

  for (i in 1:10) {

    lines <- c(
      lines,

      sprintf(
        "| %d | %s | %.3f | %.3f |",
        i,
        top10[[t]]$Genotype[i],
        top10[[t]]$BLUP[i],
        top10[[t]]$SE[i]
      )
    )
  }

  lines <- c(
    lines,
    ""
  )
}

lines <- c(
  lines,

  "## 4. AMMI Variance Explained (%)",
  "",

  "| Trait | IPCA1 (%) | IPCA2 (%) | IPCA1 Cumulative (%) | IPCA1 + IPCA2 Cumulative (%) |",
  "|---|---|---|---|---|",

  sprintf(
    "| SummerHeight | %.2f | %.2f | %.2f | %.2f |",
    ammi_tab$PC1_pct[1],
    ammi_tab$PC2_pct[1],
    ammi_tab$PC1_cum[1],
    ammi_tab$PC12_cum[1]
  ),

  sprintf(
    "| FreshWeight | %.2f | %.2f | %.2f | %.2f |",
    ammi_tab$PC1_pct[2],
    ammi_tab$PC2_pct[2],
    ammi_tab$PC1_cum[2],
    ammi_tab$PC12_cum[2]
  ),

  sprintf(
    "| DryWeight | %.2f | %.2f | %.2f | %.2f |",
    ammi_tab$PC1_pct[3],
    ammi_tab$PC2_pct[3],
    ammi_tab$PC1_cum[3],
    ammi_tab$PC12_cum[3]
  ),

  "",

  "## 5. GGE Variance Explained (%)",
  "",

  "| Trait | PC1 (%) | PC2 (%) | PC1 + PC2 Cumulative (%) |",
  "|---|---|---|---|",

  sprintf(
    "| SummerHeight | %.2f | %.2f | %.2f |",
    gge_tab$PC1_pct[1],
    gge_tab$PC2_pct[1],
    gge_tab$PC12_cum[1]
  ),

  sprintf(
    "| FreshWeight | %.2f | %.2f | %.2f |",
    gge_tab$PC1_pct[2],
    gge_tab$PC2_pct[2],
    gge_tab$PC12_cum[2]
  ),

  sprintf(
    "| DryWeight | %.2f | %.2f | %.2f |",
    gge_tab$PC1_pct[3],
    gge_tab$PC2_pct[3],
    gge_tab$PC12_cum[3]
  ),

  "",

  "## 6. Correlation Between SummerHeight and FreshWeight BLUPs (Pearson)",
  "",

  sprintf(
    "- r = %.4f",
    cor_res$r
  ),

  sprintf(
    "- p = %.4g",
    cor_res$p
  ),

  sprintf(
    "- 95%% CI = [%.4f, %.4f]",
    cor_res$lo,
    cor_res$hi
  ),

  sprintf(
    "- n = %d",
    cor_res$n
  ),

  "",

  "---",

  "",

  "## Notes",
  "",

  "- H² is interpreted only as **broad-sense heritability** on a genotype-mean basis and does not represent additive genetic variance.",

  "- The 2002 FW/DW anomaly was handled using a **conditional swap** based on the physical constraint FW ≥ DW; only rows with DW > FW were swapped.",

  "- Raw observations were aggregated to the **genotype × year × replicate** level to avoid pseudoreplication.",

  "- AMMI/GGE analyses were based on only **four environments (n_env = 4)**; therefore, the stability interpretation of the principal axes should be considered exploratory.",

  ""
)

writeLines(
  lines,
  file.path(ROOT, "results_summary.md")
)

cat(
  "\nAnalysis results written to results_summary.md and 04_Results/Genotype43/\n"
)
