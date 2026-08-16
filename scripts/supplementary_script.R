# ============================================================================
# supplementary_script.R
#
# Reproducible analysis for:
# "Preliminary evaluation of Sino-Australian alfalfa germplasm using
# historical data from a single site in Gansu: heritability, BLUP, and
# genotype × year interaction"
#
# INPUT:
#   supplementary_data.csv
#   (stored in the same directory as this script)
#
# OUTPUT:
#   results/
#   (contains tables reproducing the manuscript results)
#
# HOW TO RUN:
#   Rscript supplementary_script.R
#
# The script automatically detects its own directory when invoked through
# Rscript, allowing the analysis to use reproducible relative paths without
# exposing any user-specific absolute path.
#
# ANALYSIS PIPELINE
#   1. Read anonymized raw plot-level data from 2002–2005.
#   2. Correct the 2002 fresh/dry weight column-reversal error.
#   3. Aggregate sub-samples to the genotype × year × replicate level
#      to avoid pseudoreplication.
#   4. Retain the 43 genotypes with complete records for all four years.
#   5. Fit linear mixed models:
#         trait ~ Year + (1|Genotype) + (1|Genotype:Year)
#                 + (1|Year:Rep)
#      with Year as a fixed effect and the remaining terms as random effects.
#      Models are fitted using REML.
#   6. Estimate variance components, broad-sense heritability (H²) with
#      95% confidence intervals using parametric bootstrap, BLUPs, AMMI,
#      GGE, and the correlation between SummerHeight and FreshWeight BLUPs.
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
})

# ----------------------------------------------------------------------------
# Resolve the directory containing this script
# ----------------------------------------------------------------------------
# This avoids hard-coded absolute paths and makes the analysis portable across
# computers and operating systems.
#
# When executed with:
#   Rscript supplementary_script.R
# the --file argument identifies the script location directly.
#
# When run interactively, getwd() is used as a fallback.
# ----------------------------------------------------------------------------

args <- commandArgs(
  trailingOnly = FALSE
)

file_arg <- grep(
  "--file=",
  args,
  value = TRUE
)

ROOT <- if (length(file_arg) > 0) {
  dirname(
    sub(
      "--file=",
      "",
      file_arg[1]
    )
  )
} else {
  getwd()
}

DATA_FILE <- file.path(
  ROOT,
  "supplementary_data.csv"
)

OUT <- file.path(
  ROOT,
  "results"
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

Y <- 4  # Number of years (2002–2005)

# ----------------------------------------------------------------------------
# 1. Read anonymized data
# ----------------------------------------------------------------------------

d0 <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  na.strings = c(
    "",
    "NA"
  )
)

# Expected input columns:
# genotype_code,
# year,
# rep,
# plot,
# summer_height_cm,
# fresh_weight_kg5m2,
# dry_weight_kg5m2

names(d0) <- c(
  "Genotype",
  "Year",
  "Rep",
  "Plot",
  "SummerHeight",
  "FreshWeight",
  "DryWeight"
)

d0$Genotype <- factor(
  d0$Genotype
)

d0$Year <- factor(
  d0$Year
)

d0$Rep <- factor(
  d0$Rep
)

# ----------------------------------------------------------------------------
# 2. Correct the 2002 fresh/dry weight column-reversal error
# ----------------------------------------------------------------------------
# In part of the 2002 records, dry weight exceeded fresh weight.
# Based on the physical constraint FW >= DW, these records are interpreted
# as having the fresh- and dry-weight columns reversed and are swapped back.
# ----------------------------------------------------------------------------

swap <- d0$Year == 2002 &
        !is.na(d0$FreshWeight) &
        !is.na(d0$DryWeight) &
        d0$DryWeight > d0$FreshWeight

fw_old <- d0$FreshWeight[
  swap
]

d0$FreshWeight[
  swap
] <- d0$DryWeight[
  swap
]

d0$DryWeight[
  swap
] <- fw_old

cat(
  sprintf(
    "2002 fresh/dry weight records corrected: %d\n",
    sum(swap)
  )
)

# ----------------------------------------------------------------------------
# 3. Aggregate multiple sub-sample rows to one value per
#    genotype × year × replicate
# ----------------------------------------------------------------------------

core_rep <- d0 %>%
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
          mean(
            .,
            na.rm = TRUE
          )
        }
    ),
    .groups = "drop"
  )

# ----------------------------------------------------------------------------
# 4. Retain the 43 genotypes complete in all four years
#    for all three traits
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
  group_by(
    Genotype
  ) %>%
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
  pull(
    Genotype
  )

d43 <- core_rep %>%
  filter(
    Genotype %in% complete43
  ) %>%
  mutate(
    Genotype = factor(
      Genotype
    ),

    Year = factor(
      Year
    ),

    Rep = factor(
      Rep
    )
  )

cat(
  sprintf(
    "Complete genotypes: %d; analysis rows: %d\n",
    length(complete43),
    nrow(d43)
  )
)

# ----------------------------------------------------------------------------
# 5. Linear mixed model, variance components, and heritability
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

fit_models <- setNames(
  lapply(
    TRAITS,
    function(t) {

      dd <- d43 %>%
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
        data = dd,
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
# Harmonic mean number of replicates per genotype × year cell
# ----------------------------------------------------------------------------

harm_mean_rep <- function(t) {

  k <- d43 %>%
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
    pull(
      n
    )

  1 / mean(
    1 / k
  )
}

# ----------------------------------------------------------------------------
# Broad-sense heritability
#
# H² = V_G /
#      (V_G + V_GY/Y + V_e/(Y × r))
#
# where:
#   V_G  = genotype variance
#   V_GY = genotype × year variance
#   V_e  = residual variance
#   Y    = number of years
#   r    = harmonic mean number of replicates
# ----------------------------------------------------------------------------

h2_fun <- function(t) {

  r <- harm_mean_rep(
    t
  )

  function(m) {

    V_G <- get_vc(
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
# Variance components
# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------
# Broad-sense heritability with 95% CI
# Parametric bootstrap, 1000 replicates
# ----------------------------------------------------------------------------

set.seed(
  123
)

h2_tab <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      m <- fit_models[[t]]

      h2 <- h2_fun(t)(
        m
      )

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
          c(
            0.025,
            0.975
          ),
          na.rm = TRUE
        )

      } else {

        c(
          NA,
          NA
        )
      }

      data.frame(
        Trait = t,
        H2 = h2,
        H2_lo = ci[1],
        H2_hi = ci[2],
        r_harmonic = harm_mean_rep(t),

        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  vc_tab,
  file.path(
    OUT,
    "variance_components.csv"
  ),
  row.names = FALSE
)

write.csv(
  h2_tab,
  file.path(
    OUT,
    "heritability.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 6. BLUPs
#    Genotypic value = model intercept + random genotype effect
# ----------------------------------------------------------------------------

blup_list <- lapply(
  TRAITS,
  function(t) {

    m <- fit_models[[t]]

    rr <- ranef(
      m,
      condVar = TRUE
    )

    data.frame(

      Genotype = rownames(
        rr$Genotype
      ),

      BLUP = as.numeric(
        fixef(m)[1] +
          rr$Genotype[[1]]
      ),

      SE = sqrt(
        attr(
          rr$Genotype,
          "postVar"
        )[1, 1, ]
      ),

      stringsAsFactors = FALSE
    )
  }
)

names(
  blup_list
) <- TRAITS

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
  file.path(
    OUT,
    "BLUP_all.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 7. AMMI and GGE
#    Base R SVD; no external AMMI/GGE packages required
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
    arrange(
      Genotype
    )
}

# ----------------------------------------------------------------------------
# AMMI
# ----------------------------------------------------------------------------

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

    # Remove genotype and environment main effects
    I <- X -
         gr -
         rep(
           cg,
           each = nrow(X)
         ) +
         gm

    svd(I)
  }
)

names(
  ammi
) <- TRAITS

# ----------------------------------------------------------------------------
# GGE
# ----------------------------------------------------------------------------

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

    # Environment-centered matrix
    Xc <- scale(
      X,
      center = colMeans(X),
      scale = FALSE
    )

    svd(Xc)
  }
)

names(
  gge
) <- TRAITS

# ----------------------------------------------------------------------------
# AMMI / GGE variance explained
# ----------------------------------------------------------------------------

pc_tab <- do.call(
  rbind,
  lapply(
    TRAITS,
    function(t) {

      a <- ammi[[t]]
      g <- gge[[t]]

      ev_a <- a$d^2 /
              sum(a$d^2) *
              100

      ev_g <- g$d^2 /
              sum(g$d^2) *
              100

      data.frame(

        Trait = t,

        AMMI_PC1 = round(
          ev_a[1],
          2
        ),

        AMMI_PC2 = round(
          ev_a[2],
          2
        ),

        AMMI_cum = round(
          ev_a[1] +
            ev_a[2],
          2
        ),

        GGE_PC1 = round(
          ev_g[1],
          2
        ),

        GGE_PC2 = round(
          ev_g[2],
          2
        ),

        GGE_cum = round(
          ev_g[1] +
            ev_g[2],
          2
        ),

        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  pc_tab,
  file.path(
    OUT,
    "ammi_gge_pc.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 8. Correlation between SummerHeight and FreshWeight BLUPs
# ----------------------------------------------------------------------------

ct <- cor.test(
  blup_list[["SummerHeight"]]$BLUP,
  blup_list[["FreshWeight"]]$BLUP
)

cor_tab <- data.frame(

  r = ct$estimate,

  r_lo = ct$conf.int[1],

  r_hi = ct$conf.int[2],

  p = ct$p.value,

  n = length(
    blup_list[[1]]$Genotype
  ),

  stringsAsFactors = FALSE
)

write.csv(
  cor_tab,
  file.path(
    OUT,
    "correlation.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# Print summary
# ----------------------------------------------------------------------------

cat(
  "\n=== Variance Components ===\n"
)

print(
  vc_tab,
  row.names = FALSE
)

cat(
  "\n=== Broad-Sense Heritability (95% CI) ===\n"
)

print(
  h2_tab,
  row.names = FALSE
)

cat(
  "\n=== AMMI / GGE Variance Explained (%) ===\n"
)

print(
  pc_tab,
  row.names = FALSE
)

cat(
  "\n=== Summer Height vs Fresh Weight BLUP Correlation ===\n"
)

print(
  cor_tab,
  row.names = FALSE
)

cat(
  sprintf(
    "\nAnalysis completed. Results written to: %s\n",
    OUT
  )
)
