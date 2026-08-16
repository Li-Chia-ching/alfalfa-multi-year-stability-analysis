# ============================================================================
# analysis_86.R — Mixed-model analysis using all 86 genotypes
# ----------------------------------------------------------------------------
# Purpose:
#   Core analysis for Plan B using the full set of 86 genotypes.
#
# Data:
#   Four annual CSV files from 2002–2005 in 01_Raw_Phenotype_Data/
#
# Methods:
#   lme4 linear mixed models
#   Year is treated as a fixed effect.
#   Genotype, Genotype:Year, and Year:Rep are treated as random effects.
#
# Traits:
#   SummerHeight
#   FreshWeight
#   DryWeight
#
# Notes:
#   1) Original columns are mapped as:
#      Line -> Genotype
#      Mean_Summer -> SummerHeight
#      FW -> FreshWeight
#      DW -> DryWeight
#
#   2) Year is extracted from the input file name.
#
#   3) Text markers and non-numeric entries are converted to NA through
#      as.numeric() with warnings suppressed.
#
#   4) The 2002 FW/DW columns contain anomalies. Based on the physical
#      constraint FW >= DW, rows with DW > FW are conditionally swapped.
#
#   5) Raw observations are aggregated to the genotype × year × replicate
#      level to avoid pseudoreplication.
#
#   6) Broad-sense heritability is estimated on a genotype-mean basis.
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
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

OUT <- file.path(
  ROOT,
  "04_Results",
  "Genotype86"
)

dir.create(
  OUT,
  showWarnings = FALSE,
  recursive = TRUE
)

TRAITS <- c(
  "SummerHeight",
  "FreshWeight",
  "DryWeight"
)

Y <- 4  # Number of years (environments)

# ----------------------------------------------------------------------------
# 1. Read annual data
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

# ----------------------------------------------------------------------------
# 2. Column mapping and data cleaning
# ----------------------------------------------------------------------------
num <- function(x) {
  suppressWarnings(as.numeric(x))
}

raw$Genotype <- trimws(raw$Line)
raw$Rep      <- as.integer(raw$Rep)

raw$SummerHeight <- num(raw$Mean_Summer)
raw$FreshWeight  <- num(raw$FW)
raw$DryWeight    <- num(raw$DW)

# ----------------------------------------------------------------------------
# 3. Conditional FW/DW swap for 2002
# ----------------------------------------------------------------------------
# Physical constraint:
#   FreshWeight >= DryWeight
#
# Only rows where DW > FW are swapped.
# ----------------------------------------------------------------------------
swap <- raw$Year == 2002 &
        !is.na(raw$FreshWeight) &
        !is.na(raw$DryWeight) &
        raw$DryWeight > raw$FreshWeight

fw_old <- raw$FreshWeight[swap]

raw$FreshWeight[swap] <- raw$DryWeight[swap]
raw$DryWeight[swap]   <- fw_old

cat(
  "Number of rows with conditional FW/DW swaps in 2002:",
  sum(swap),
  "\n"
)

# ----------------------------------------------------------------------------
# 4. Aggregate to genotype × year × replicate experimental units
# ----------------------------------------------------------------------------
core_rep <- raw %>%
  group_by(
    Year,
    Genotype,
    Rep
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

cat(
  "Number of genotype × year × replicate rows:",
  nrow(core_rep),
  "\n"
)

# ----------------------------------------------------------------------------
# 5. Prepare data for mixed-model analysis
# ----------------------------------------------------------------------------
d86 <- core_rep %>%
  mutate(
    Genotype = factor(Genotype),
    Year     = factor(Year),
    Rep      = factor(Rep)
  )

# ----------------------------------------------------------------------------
# 6. Extract variance components
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# 7. Fit mixed models
# ----------------------------------------------------------------------------
fit_models <- setNames(
  lapply(
    TRAITS,
    function(t) {

      d <- d86 %>%
        filter(
          !is.na(.data[[t]])
        )

      formula_t <- as.formula(
        paste0(
          t,
          " ~ Year + (1|Genotype) + ",
          "(1|Genotype:Year) + (1|Year:Rep)"
        )
      )

      lmer(
        formula_t,
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

# ----------------------------------------------------------------------------
# 8. Harmonic mean number of replicates
# ----------------------------------------------------------------------------
harm_mean_rep <- function(t) {

  k <- d86 %>%
    filter(
      !is.na(.data[[t]])
    ) %>%
    group_by(
      Genotype,
      Year
    ) %>%
    summarise(
      n = n(),
      .groups = "drop"
    ) %>%
    pull(n)

  1 / mean(1 / k)
}

# ----------------------------------------------------------------------------
# 9. Broad-sense heritability function
# ----------------------------------------------------------------------------
# H² = V_G / (V_G + V_GY/Y + V_e/(Y × r))
#
# where:
#   V_G  = genotype variance
#   V_GY = genotype × year variance
#   V_e  = residual variance
#   Y    = number of environments
#   r    = harmonic mean number of replicates
# ----------------------------------------------------------------------------
h2_fun <- function(t) {

  r <- harm_mean_rep(t)

  function(m) {

    V_G  <- get_vc(
      m,
      "Genotype"
    )

    V_GY <- get_vc(
      m,
      "Genotype:Year"
    )

    V_e <- sigma(m)^2

    V_G / (
      V_G +
        V_GY / Y +
        V_e / (Y * r)
    )
  }
}

# ----------------------------------------------------------------------------
# 10. Variance components
# ----------------------------------------------------------------------------
vc_tab <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      m <- fit_models[[t]]

      data.frame(
        Trait = t,
        V_G   = get_vc(m, "Genotype"),
        V_GY  = get_vc(m, "Genotype:Year"),
        V_RY  = get_vc(m, "Year:Rep"),
        V_e   = sigma(m)^2,
        stringsAsFactors = FALSE
      )
    }
  )
)

# ----------------------------------------------------------------------------
# 11. Broad-sense heritability with 95% CI
# ----------------------------------------------------------------------------
# Parametric bootstrap using bootMer(), nsim = 1000
# ----------------------------------------------------------------------------
set.seed(123)

h2_tab <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      m <- fit_models[[t]]

      h2 <- h2_fun(t)(m)

      b <- tryCatch(
        bootMer(
          m,
          h2_fun(t),
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
        Trait      = t,
        H2         = h2,
        H2_lo      = ci[1],
        H2_hi      = ci[2],
        r_harmonic = harm_mean_rep(t),
        stringsAsFactors = FALSE
      )
    }
  )
)

# ----------------------------------------------------------------------------
# 12. BLUP extraction and standard errors
# ----------------------------------------------------------------------------
blup_list <- lapply(
  TRAITS,
  function(t) {

    m <- fit_models[[t]]

    rr <- ranef(
      m,
      condVar = TRUE
    )

    genotype_effect <- rr$Genotype[[1]]

    posterior_var <- attr(
      rr$Genotype,
      "postVar"
    )

    se <- sqrt(
      posterior_var[1, 1, ]
    )

    data.frame(
      Genotype = rownames(rr$Genotype),

      # Genotypic value = fixed intercept + genotype random effect
      BLUP = as.numeric(
        fixef(m)[1] + genotype_effect
      ),

      SE = se,

      stringsAsFactors = FALSE
    )
  }
)

names(blup_list) <- TRAITS

# ----------------------------------------------------------------------------
# 13. Combine BLUPs across traits
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# 14. Correlation between SummerHeight and FreshWeight BLUPs
# ----------------------------------------------------------------------------
ct <- cor.test(
  blup_list[["SummerHeight"]]$BLUP,
  blup_list[["FreshWeight"]]$BLUP,
  method = "pearson"
)

# ----------------------------------------------------------------------------
# 15. Export analysis results
# ----------------------------------------------------------------------------
write.csv(
  vc_tab,
  file.path(
    OUT,
    "variance_components_86.csv"
  ),
  row.names = FALSE
)

write.csv(
  h2_tab,
  file.path(
    OUT,
    "heritability_86.csv"
  ),
  row.names = FALSE
)

write.csv(
  blup_out,
  file.path(
    OUT,
    "BLUP_86.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 16. Print summary results to the console
# ----------------------------------------------------------------------------
cat(
  "\n=== Variance Components (86 Genotypes) ===\n"
)

print(
  vc_tab,
  row.names = FALSE
)

cat(
  "\n=== Broad-Sense Heritability H² (86 Genotypes) ===\n"
)

print(
  h2_tab,
  row.names = FALSE
)

cat(
  "\n=== SummerHeight–FreshWeight BLUP Correlation (86 Genotypes) ===\n"
)

cat(
  sprintf(
    "r = %.4f  p = %.4g  95%% CI = [%.4f, %.4f]  n = %d\n",
    ct$estimate,
    ct$p.value,
    ct$conf.int[1],
    ct$conf.int[2],
    length(
      blup_list[[1]]$Genotype
    )
  )
)

cat(
  "\n=== Top 10 BLUPs for Each Trait (86 Genotypes) ===\n"
)

for (t in TRAITS) {

  b <- blup_list[[t]] %>%
    arrange(
      desc(BLUP)
    ) %>%
    head(10)

  cat(
    sprintf(
      "\n[%s]\n",
      t
    )
  )

  print(
    b,
    row.names = FALSE
  )
}
