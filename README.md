
---

# Multi-year Evaluation of Sino-Australian Alfalfa Cultivars in Gansu

This repository contains the R scripts used for the statistical analysis, stability analysis, and figure generation reported in:

**"Multi-year evaluation of Sino-Australian alfalfa cultivars in Gansu: heritability, BLUP, and GGE/AMMI stability analysis"**

The workflow evaluates multi-year phenotypic performance, genetic variation, genotype × year interaction, genotypic BLUPs, broad-sense heritability, and stability patterns of Sino-Australian alfalfa germplasm evaluated at a single site in Gansu, China.

The analysis emphasizes transparent data processing and reproducibility. Raw observations are aggregated to the **genotype × year × replicate** experimental-unit level before statistical modeling to avoid pseudoreplication. Two complementary genotype sets are analyzed: **43 genotypes with complete records across all four years and all three traits**, and **all 86 available genotypes** for a full-data mixed-model sensitivity analysis.

---

## 1. Analysis Overview

The repository currently contains three complementary analysis components.

### Plan A — Main analysis: 43 complete genotypes

The primary multi-year analysis uses the **43 genotypes with complete observations for SummerHeight, FreshWeight, and DryWeight in all four years (2002–2005)**.

The workflow includes:

1. Import of annual phenotype files.
2. Harmonization of trait names and data types.
3. Conversion of non-numeric entries to missing values.
4. Conditional correction of the 2002 fresh-weight/dry-weight column reversal.
5. Aggregation to the genotype × year × replicate level.
6. Identification of the 43 four-year complete genotypes.
7. Linear mixed-model analysis.
8. Variance-component estimation.
9. Broad-sense heritability estimation with parametric bootstrap confidence intervals.
10. Genotypic BLUP estimation and ranking.
11. AMMI analysis using base R singular value decomposition (SVD).
12. GGE analysis using base R SVD.
13. BLUP correlation analysis between summer plant height and fresh weight.
14. Generation of manuscript figures.

### Plan B — Full-data mixed-model analysis: 86 genotypes

`analysis_86.R` evaluates the full set of **86 available genotypes** without restricting the analysis to four-year complete cases.

This analysis uses the same genotype × year × replicate aggregation and the same linear mixed-model framework as Plan A. It provides:

* variance components;
* broad-sense heritability;
* BLUPs and BLUP standard errors; and
* the Pearson correlation between SummerHeight and FreshWeight BLUPs.

Plan B is intended as a **full-data sensitivity analysis** rather than a replacement for the balanced 43-genotype analysis.

### Supplementary reproducibility script

`supplementary_script.R` is designed for distribution with supplementary materials. It uses an anonymized `supplementary_data.csv` file stored in the same directory as the script and automatically resolves its own directory when executed with `Rscript`.

It reproduces the core 43-genotype analysis, including:

* variance components;
* broad-sense heritability;
* BLUPs;
* AMMI;
* GGE; and
* BLUP correlation.

No user-specific absolute paths are required.

---

## 2. Repository Structure

The repository is organized around analysis scripts, raw input data, and generated results.

```text
Sino-Australian_Alfalfa_Project/
│
├── analysis_main.R
├── analysis_86.R
├── make_figures.R
│
├── 01_Raw_Phenotype_Data/
│   ├── Sino-Australian Alfalfa Project Data - 2002data.csv
│   ├── Sino-Australian Alfalfa Project Data - 2003data.csv
│   ├── Sino-Australian Alfalfa Project Data - 2004data.csv
│   └── Sino-Australian Alfalfa Project Data - 2005data.csv
│
├── 04_Results/
│   ├── Genotype43/
│   │   ├── BLUP_43.csv
│   │   ├── AMMI_genotype_IPCA_*.csv
│   │   ├── AMMI_env_IPCA_*.csv
│   │   └── Figures/
│   │       ├── Fig1_missing_pattern.png
│   │       ├── Fig1_missing_pattern.pdf
│   │       ├── Fig2_subset_comparison.png
│   │       ├── Fig2_subset_comparison.pdf
│   │       ├── Fig3_variance_components.png
│   │       ├── Fig3_variance_components.pdf
│   │       ├── Fig4_blup_ranking.png
│   │       ├── Fig4_blup_ranking.pdf
│   │       ├── Fig5_ammi_biplot.png
│   │       ├── Fig5_ammi_biplot.pdf
│   │       ├── Fig6_gge_biplot.png
│   │       ├── Fig6_gge_biplot.pdf
│   │       ├── Fig7_blup_scatter.png
│   │       └── Fig7_blup_scatter.pdf
│   │
│   └── Genotype86/
│       ├── variance_components_86.csv
│       ├── heritability_86.csv
│       └── BLUP_86.csv
│
├── supplementary_materials/
│   ├── supplementary_script.R
│   ├── supplementary_data.csv
│   └── results/
│       ├── variance_components.csv
│       ├── heritability.csv
│       ├── BLUP_all.csv
│       ├── ammi_gge_pc.csv
│       └── correlation.csv
│
└── results_summary.md
```

Generated result directories and individual output files may not be present in a fresh clone until the corresponding scripts are executed.

---

## 3. Input Data

The main analysis expects four annual CSV files in:

```text
01_Raw_Phenotype_Data/
```

The scripts read the 2002, 2003, 2004, and 2005 files and harmonize the relevant columns as follows:

| Original column | Analysis variable |
| --------------- | ----------------- |
| `Line`          | `Genotype`        |
| `Rep`           | `Rep`             |
| `Mean_Summer`   | `SummerHeight`    |
| `FW`            | `FreshWeight`     |
| `DW`            | `DryWeight`       |

The analysis uses three traits:

* `SummerHeight`
* `FreshWeight`
* `DryWeight`

Year is extracted from the file name.

Raw observations containing non-numeric markers are converted to missing values during numeric conversion.

### Data privacy and repository distribution

The main raw phenotype files may be excluded from the public repository because of data availability restrictions.

The supplementary workflow instead uses an anonymized file:

```text
supplementary_data.csv
```

with the following expected columns:

```text
genotype_code
year
rep
plot
summer_height_cm
fresh_weight_kg5m2
dry_weight_kg5m2
```

---

## 4. Data Cleaning and Experimental-Unit Definition

### 4.1 Conditional correction of the 2002 FW/DW column reversal

For the 2002 records, rows in which:

```text
DryWeight > FreshWeight
```

are conditionally interpreted as having the fresh- and dry-weight columns reversed.

Only those rows are swapped.

This correction is based on the physical constraint that fresh weight should not be lower than dry weight.

### 4.2 Aggregation to the experimental-unit level

Where multiple sub-sample observations occur within the same:

```text
Genotype × Year × Replicate
```

combination, the observations are averaged before statistical analysis.

This aggregation prevents multiple sub-samples from being incorrectly treated as independent biological replicates and therefore avoids pseudoreplication.

### 4.3 Complete-genotype definition

For Plan A, a genotype is classified as complete when all three traits have at least one non-missing observation in **each of the four years**.

The resulting primary analysis set contains:

```text
43 genotypes × 4 years
```

---

## 5. Linear Mixed Model

For each trait, the following linear mixed model is fitted using `lme4::lmer()` with restricted maximum likelihood (REML):

```text
Trait ~ Year + (1|Genotype) + (1|Genotype:Year) + (1|Year:Rep)
```

Equivalently,

$$
Y_{ijk} = \mu + Year_i + G_j + (G \times Year)_{ij} + Rep(Year)_{ik} + e_{ijk}
$$

where:

* `Year` is treated as a **fixed environmental effect**;
* `Genotype` is treated as a **random genetic effect**;
* `Genotype × Year` is treated as a **random genotype-by-year interaction effect**;
* `Year:Rep` represents the replicate/block component within year; and
* $e$ is the residual error.

The model is fitted separately for:

* SummerHeight
* FreshWeight
* DryWeight

---

## 6. Variance Components

The mixed model provides the following variance components:

| Symbol   | Model term      | Interpretation                       |
| -------- | --------------- | ------------------------------------ |
| $V_G$    | `Genotype`      | Genotypic variance                   |
| $V_{GY}$ | `Genotype:Year` | Genotype × year interaction variance |
| $V_{RY}$ | `Year:Rep`      | Replicate/block variance within year |
| $V_e$    | Residual        | Residual variance                    |

These components are reported for each trait.

---

## 7. Broad-Sense Heritability

Broad-sense heritability is estimated on a genotype-mean basis across the four years as:

$$
H^2 =
\frac{V_G}
{V_G + \frac{V_{GY}}{Y} + \frac{V_e}{Yr}}
$$

where:

* $V_G$ is the genotype variance component;
* $V_{GY}$ is the genotype × year interaction variance component;
* $V_e$ is the residual variance;
* $Y = 4$ is the number of years/environments; and
* $r$ is the harmonic mean number of replicates per genotype-year cell.

Confidence intervals are estimated using a **parametric bootstrap with 1,000 `bootMer()` replicates**.

Throughout the analysis, $H^2$ is interpreted as **broad-sense heritability**. It is not interpreted as additive genetic variance or narrow-sense heritability.

---

## 8. BLUP Prediction

Genotypic values are estimated using best linear unbiased prediction (BLUP).

For each genotype:

$$
BLUP_g = \hat{\mu} + \hat{g}
$$

where:

* $\hat{\mu}$ is the fitted model intercept; and
* $\hat{g}$ is the predicted random genotype effect.

BLUP standard errors are obtained from the conditional variance of the random effects.

BLUPs are used for:

* genotype ranking;
* visualization of genotypic performance; and
* correlation analysis among traits.

---

## 9. AMMI Analysis

Additive Main Effects and Multiplicative Interaction (AMMI) analysis is conducted using the genotype × year mean matrix.

The interaction matrix is obtained after removing additive genotype and year main effects:

$$
I = X - \bar{X}*{g} - \bar{X}*{e} + \bar{X}
$$

Singular value decomposition is then applied using base R:

```r
svd(I)
```

No external AMMI package is required.

The resulting principal interaction coordinates are used for AMMI biplots.

Because the analysis contains only four years/environments, AMMI principal-axis interpretation should be considered **exploratory**, particularly when assessing stability.

---

## 10. GGE Analysis

GGE analysis is based on the genotype + genotype × environment component of the genotype × year matrix.

The matrix is centered by environment before singular value decomposition:

```r
Xc <- scale(X, center = colMeans(X), scale = FALSE)
svd(Xc)
```

The resulting genotype and environment coordinates are used to construct GGE biplots, including the **which-won-where** representation.

As with AMMI, the relatively small number of environments limits the strength of stability inference.

---

## 11. Correlation of BLUPs

Pearson correlation is calculated between genotype BLUPs for:

* SummerHeight
* FreshWeight

The correlation analysis uses the same genotype-level BLUPs derived from the mixed models.

The reported statistics include:

* Pearson correlation coefficient ($r$);
* P value;
* 95% confidence interval; and
* sample size.

---

## 12. Figure Generation

`make_figures.R` generates the manuscript figures for the 43-genotype analysis in both PNG and PDF formats.

The current script generates:

| Figure | Content                                                          |
| ------ | ---------------------------------------------------------------- |
| Fig. 1 | Data availability pattern across 86 genotypes and four years     |
| Fig. 2 | Complete vs incomplete genotype subset comparison                |
| Fig. 3 | Proportion of phenotypic variance attributed to model components |
| Fig. 4 | BLUP ranking with 95% confidence intervals                       |
| Fig. 5 | AMMI biplots                                                     |
| Fig. 6 | GGE which-won-where biplots                                      |
| Fig. 7 | Summer height vs fresh weight BLUP correlation                   |

Figures are exported at:

```text
04_Results/Genotype43/Figures/
```

PNG files are generated at **300 dpi**, and PDF versions are also produced.

---

## 13. Reproducible Execution

### Main 43-genotype analysis

Run from the project root:

```bash
Rscript analysis_main.R
```

The script writes the main tabular results to:

```text
04_Results/Genotype43/
```

and generates:

```text
results_summary.md
```

in the project root.

### Full 86-genotype analysis

Run:

```bash
Rscript analysis_86.R
```

Outputs are written to:

```text
04_Results/Genotype86/
```

### Figure generation

After the main 43-genotype analysis data are available, run:

```bash
Rscript make_figures.R
```

Figures are written to:

```text
04_Results/Genotype43/Figures/
```

### Supplementary reproducibility script

Place the anonymized supplementary data beside the script:

```text
supplementary_materials/
├── supplementary_script.R
└── supplementary_data.csv
```

Then run:

```bash
Rscript supplementary_script.R
```

Results are written to:

```text
supplementary_materials/results/
```

The supplementary script automatically resolves its own directory and therefore does not require a machine-specific absolute path.

---

## 14. Software Requirements

The current scripts use R packages directly rather than a broad `tidyverse` meta-package.

Recommended R version:

```text
R >= 4.3
```

Required packages for the main analyses:

```r
install.packages(c(
  "lme4",
  "dplyr",
  "tidyr"
))
```

Required packages for figure generation:

```r
install.packages(c(
  "lme4",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "scales"
))
```

AMMI and GGE decomposition are implemented with base R `svd()` and therefore do not require `agricolae`, `metan`, or other specialized stability-analysis packages.

---

## 15. Reproducibility and Path Handling

All project analysis scripts use portable relative paths.

The main project scripts are designed to be executed from the project root, while the supplementary script automatically resolves the directory in which it is located.

No user-specific Windows paths, usernames, or local workstation directories are embedded in the analysis code.

This design allows the repository to be copied to another computer or cloned from version control without editing machine-specific paths.

---

## 16. Relationship Between the 43- and 86-Genotype Analyses

The two genotype sets serve different analytical purposes.

### 43-genotype analysis

The 43-genotype dataset provides a balanced set of genotypes with complete observations across all four years and all three traits. It is used for the principal multi-year genotype × year analysis and for the manuscript figures.

### 86-genotype analysis

The 86-genotype analysis retains all available genotypes and therefore maximizes the available sample size for mixed-model estimation.

Because the completeness of observations differs among genotypes, this analysis is used as a complementary **sensitivity analysis** rather than being treated as identical to the balanced 43-genotype analysis.

Results from the two analyses should therefore be compared for consistency rather than interpreted as estimates from the same analytical population.

---

## 17. Notes on Interpretation

This repository is intended to reproduce the computational analyses underlying the manuscript.

In particular:

* broad-sense heritability should not be interpreted as narrow-sense heritability or additive genetic variance;
* BLUPs represent model-based genotype-level predictions rather than raw phenotype means;
* AMMI and GGE analyses are based on only four years/environments and should therefore be interpreted cautiously;
* the 43-genotype analysis and 86-genotype analysis represent different analytical populations; and
* the conditional 2002 FW/DW correction is an explicit data-cleaning step and should be retained when reproducing the reported results.
