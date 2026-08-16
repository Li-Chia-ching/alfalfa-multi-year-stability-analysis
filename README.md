
# Multi-year Evaluation of Sino-Australian Alfalfa Cultivars in Gansu

This repository contains the R scripts used for the statistical analysis, exploratory genotype × year stability analysis, and manuscript figure generation reported in:

**"Multi-year evaluation of Sino-Australian alfalfa cultivars in Gansu: heritability, BLUP, and GGE/AMMI stability analysis"**

The workflow evaluates multi-year phenotypic performance, genetic variation, genotype × year interaction, broad-sense heritability, genotypic BLUPs, genotype ranking, establishment-cohort differences, and exploratory stability patterns of Sino-Australian alfalfa germplasm evaluated at a single site in Gansu, China.

The current repository distinguishes two analytical populations according to their roles in the manuscript:

- **All 86 available genotypes** are used for the main mixed-model analysis and the main manuscript figures.
- **The 43 genotypes with complete observations across all four years and all three traits** are used as a balanced exploratory subset for AMMI and GGE analyses and the supplementary stability figures.

This separation is intentional. Mixed-model estimation can accommodate the unbalanced 86-genotype dataset, whereas the exploratory AMMI/GGE workflow is restricted to the balanced 43-genotype subset.

---

## 1. Analysis Overview

### 1.1 Main manuscript analysis: all 86 genotypes

The main analysis retains all **86 available genotypes** and uses a linear mixed model fitted by restricted maximum likelihood (REML).

For each trait, the model is:

```text
Trait ~ Year + (1|Genotype) + (1|Genotype:Year) + (1|Year:Rep)
```

The main analysis provides:

- genotype, genotype × year, replicate-within-year, and residual variance components;
- broad-sense heritability on a genotype-mean basis;
- genotype-level BLUPs and associated standard errors;
- BLUP-based genotype ranking;
- Pearson correlations among trait BLUPs; and
- establishment-cohort comparisons using the 2003–2005 fair-comparison window.

These analyses provide the quantitative basis for **Fig. 1–Fig. 4**.

### 1.2 Exploratory balanced-subset analysis: 43 complete genotypes

A subset of **43 genotypes** satisfies the completeness criterion of having non-missing observations for all three traits in all four years (2002–2005).

This balanced subset is used for exploratory genotype × year stability analysis:

- AMMI analysis;
- GGE biplot analysis;
- AMMI Stability Value (ASV);
- GGE distance to the ideal genotype; and
- mean-performance versus stability-index plots.

These analyses provide **Fig. S1–Fig. S3** and are intended to support exploratory interpretation rather than strong claims about stability based on only four environments.

### 1.3 Reproducibility workflow

`supplementary_script.R` is intended for distribution with supplementary materials. It uses an anonymized `supplementary_data.csv` file stored beside the script and automatically resolves its own directory when executed with `Rscript`.

The supplementary workflow is designed to reproduce the analysis components selected for public distribution without requiring user-specific absolute paths.

---

## 2. Figure System

The current manuscript figure architecture is:

| Figure      |      Dataset | Content                                               | Role                   |
| ----------- | -----------: | ----------------------------------------------------- | ---------------------- |
| **Fig. 1**  | 86 genotypes | Variance partitioning and broad-sense heritability    | Main                   |
| **Fig. 2**  | 86 genotypes | Pearson correlation heatmap of trait BLUPs            | Main                   |
| **Fig. 3**  | 86 genotypes | BLUP ranking with 95% CI; top 5 genotypes highlighted | Main                   |
| **Fig. 4**  | 86 genotypes | 2001 vs. 2002 establishment-cohort comparison         | Main                   |
| **Fig. S1** | 43 genotypes | AMMI biplot                                           | Exploratory supplement |
| **Fig. S2** | 43 genotypes | GGE biplot                                            | Exploratory supplement |
| **Fig. S3** | 43 genotypes | Mean performance vs. ASV and GGE distance             | Exploratory supplement |

The previous missing-pattern heatmap, complete-vs-incomplete subset comparison, and summer-height-versus-fresh-weight BLUP scatter plot are **not part of the current manuscript figure set**.

---

## 3. Repository Structure

A typical project layout is:

```text
Sino-Australian_Alfalfa_Project/
│
├── analysis_main.R
├── analysis_86.R
├── make_figures.R
├── supplementary_script.R
│
├── 01_Raw_Phenotype_Data/
│   ├── Sino-Australian Alfalfa Project Data - 2002data.csv
│   ├── Sino-Australian Alfalfa Project Data - 2003data.csv
│   ├── Sino-Australian Alfalfa Project Data - 2004data.csv
│   └── Sino-Australian Alfalfa Project Data - 2005data.csv
│
├── 04_Results/
│   ├── Genotype43/
│   ├── Genotype86/
│   └── Manuscript_Figures/
│       ├── Fig1_variance_heritability.png
│       ├── Fig1_variance_heritability.pdf
│       ├── Fig2_blup_correlation_heatmap.png
│       ├── Fig2_blup_correlation_heatmap.pdf
│       ├── Fig3_blup_ranking.png
│       ├── Fig3_blup_ranking.pdf
│       ├── Fig4_establishment_cohorts.png
│       ├── Fig4_establishment_cohorts.pdf
│       ├── FigS1_ammi_biplot.png
│       ├── FigS1_ammi_biplot.pdf
│       ├── FigS2_gge_biplot.png
│       ├── FigS2_gge_biplot.pdf
│       ├── FigS3_mean_stability.png
│       └── FigS3_mean_stability.pdf
│
├── supplementary_materials/
│   ├── supplementary_script.R
│   ├── supplementary_data.csv
│   └── results/
│
└── results_summary.md
```

Generated result directories and individual output files may not be present in a fresh clone until the corresponding scripts are executed.

---

## 4. Input Data

The main workflow expects four annual phenotype files in:

```text
01_Raw_Phenotype_Data/
```

The scripts read the 2002, 2003, 2004, and 2005 files and harmonize the relevant columns as follows:

| Original column | Analysis variable |
| --------------- | ----------------- |
| `Line`          | `Genotype`        |
| `Rep`           | `Rep`             |
| `Plot`          | `Plot`            |
| `Mean_Summer`   | `SummerHeight`    |
| `FW`            | `FreshWeight`     |
| `DW`            | `DryWeight`       |

The three analyzed traits are:

- `SummerHeight`
- `FreshWeight`
- `DryWeight`

Year is assigned from the annual input file.

Non-numeric entries are converted to missing values during numeric conversion.

---

## 5. Data Cleaning and Experimental-Unit Definition

### 5.1 Conditional correction of the 2002 FW/DW column reversal

For 2002 records, rows satisfying:

```text
DryWeight > FreshWeight
```

are interpreted as rows in which the fresh-weight and dry-weight columns were reversed.

Only those rows are swapped.

The current figure-generation script explicitly stores the original fresh-weight values before assignment so that the two columns cannot collapse to the same values during correction.

### 5.2 Genotype × year × replicate aggregation

The experimental unit for the mixed-model analysis is:

```text
Genotype × Year × Replicate
```

When multiple observations occur within the same experimental unit, the trait values are averaged before statistical modeling.

This aggregation avoids treating subsamples as independent biological replicates.

### 5.3 Definition of the 43 complete genotypes

A genotype is classified as complete when all three traits have at least one non-missing observation in every year from 2002 to 2005.

The resulting balanced subset contains:

```text
43 genotypes × 4 years
```

This subset is reserved for the exploratory AMMI/GGE stability analyses.

---

## 6. Linear Mixed Model

For each trait, the main 86-genotype analysis fits the following linear mixed model using `lme4::lmer()` with REML:

```text
Trait ~ Year + (1|Genotype) + (1|Genotype:Year) + (1|Year:Rep)
```

Equivalently,

$$
Y_{ijk} =
\mu +
Year_i +
G_j +
(G \times Year)_{ij} +
Rep(Year)_{ik} +
e_{ijk}
$$

where:

- `Year` is treated as a **fixed environmental effect**;
- `Genotype` is treated as a **random genetic effect**;
- `Genotype × Year` is treated as a **random interaction effect**;
- `Year:Rep` represents the replicate/block component nested within year; and
- $e$ is the residual error.

The model is fitted separately for:

- summer plant height;
- fresh weight; and
- dry weight.

---

## 7. Variance Components

The mixed model provides:

| Symbol   | Model term      | Interpretation                       |
| -------- | --------------- | ------------------------------------ |
| $V_G$    | `Genotype`      | Genotypic variance                   |
| $V_{GY}$ | `Genotype:Year` | Genotype × year interaction variance |
| $V_{RY}$ | `Year:Rep`      | Replicate/block variance within year |
| $V_e$    | Residual        | Residual variance                    |

For Fig. 1, these components are expressed as percentages of the sum of the four modeled variance components.

---

## 8. Broad-Sense Heritability

Broad-sense heritability is estimated on a genotype-mean basis across the four years:

$$
H^2 =
\frac{V_G}
{V_G + \frac{V_{GY}}{Y} + \frac{V_e}{Yr}}
$$

where:

- $V_G$ is the genotype variance component;
- $V_{GY}$ is the genotype × year interaction variance component;
- $V_e$ is the residual variance;
- $Y$ is the number of years/environments; and
- $r$ is the harmonic mean number of replicates per genotype × year cell.

In the current figure-generation script, the number of years is read from the analyzed data, and the harmonic mean of replication counts is calculated from the observed genotype × year cells rather than hard-coded.

The `Year × Rep` component is treated as a nuisance block component and is excluded from the genotype-mean heritability denominator.

Throughout the analysis, $H^2$ refers to **broad-sense heritability** and should not be interpreted as narrow-sense heritability or additive genetic variance.

---

## 9. BLUP Prediction and Ranking

Genotypic values are estimated using best linear unbiased prediction (BLUP).

For each genotype:

$$
BLUP_g = \hat{\mu} + \hat{g}
$$

where:

- $\hat{\mu}$ is the fitted model intercept; and
- $\hat{g}$ is the predicted random genotype effect.

BLUP standard errors are obtained from the conditional variance of the random genotype effects.

BLUPs are used for:

- genotype ranking;
- the Fig. 3 ranking plot;
- cross-trait BLUP correlation analysis; and
- genotype-level visualization.

Fig. 3 shows BLUP values with approximate 95% confidence intervals calculated as:

$$
BLUP \pm 1.96 \times SE
$$

The five highest-ranking genotypes for each trait are highlighted and labeled.

---

## 10. BLUP Correlation Analysis

Fig. 2 summarizes Pearson correlations among genotype-level BLUPs for:

- summer plant height;
- fresh weight; and
- dry weight.

The full three-trait correlation matrix is visualized as a heatmap.

The correlation analysis is based on the BLUPs estimated from the 86-genotype mixed models rather than on raw phenotypic observations.

---

## 11. Establishment-Cohort Comparison

Fig. 4 compares genotypic performance between two establishment cohorts:

- **Established in 2001**
- **Established in 2002**

The establishment cohort is derived from the plot identifier:

```text
Plot <= 129  →  2001
Plot > 129   →  2002
```

The comparison uses the **2003–2005** period as a fair comparison window.

The analysis first aggregates observations to the genotype × replicate × year × establishment-cohort level and then obtains genotype-level means within cohort.

An important feature of the original field structure is that genotype **L33 occurs in both establishment blocks**. The current workflow therefore retains establishment as a plot-level attribute rather than forcing L33 into a single cohort.

The cohort comparison is exploratory and is visualized as genotype-level distributions for the three agronomic traits.

---

## 12. AMMI Analysis

AMMI (Additive Main Effects and Multiplicative Interaction) analysis is conducted on the balanced 43-genotype × 4-year mean matrix.

The interaction matrix is obtained by removing genotype and year main effects:

$$
I = X - \bar{X}_{g} - \bar{X}_{e} + \bar{X}
$$

Singular value decomposition is then applied using base R:

```r
svd(I)
```

The first two interaction principal components are used for the AMMI biplots in **Fig. S1**.

Because the analysis contains only four environments, AMMI results are presented as **exploratory** evidence of genotype × year interaction patterns rather than as strong stability claims.

---

## 13. GGE Analysis

GGE analysis is based on the genotype + genotype × environment component of the balanced 43-genotype matrix.

The genotype × year matrix is centered by environment before singular value decomposition:

```r
Xc <- scale(X, center = colMeans(X), scale = FALSE)
svd(Xc)
```

The resulting genotype and environment coordinates are used to construct GGE biplots, including the **which-won-where** representation in **Fig. S2**.

As with AMMI, the small number of environments limits the strength of stability inference.

---

## 14. Stability Indices

Fig. S3 relates mean genotype performance to two exploratory stability metrics.

### 14.1 AMMI Stability Value

The AMMI Stability Value (ASV) is calculated from the first two AMMI genotype IPCA scores.

Lower ASV values indicate greater stability.

### 14.2 GGE distance to the ideal genotype

For each trait, the GGE distance is calculated as the Euclidean distance from the genotype to the ideal point defined by:

- maximum observed GGE PC1 score; and
- PC2 = 0.

Smaller distances indicate greater proximity to the ideal genotype in the GGE biplot.

These indices are used for exploratory visualization only and should not be interpreted as independent estimates of stability.

---

## 15. Figure Generation

`make_figures.R` generates the current manuscript figure set as both PNG and PDF files.

### Main figures

| Figure | File stem                       | Dataset |
| ------ | ------------------------------- | ------: |
| Fig. 1 | `Fig1_variance_heritability`    |      86 |
| Fig. 2 | `Fig2_blup_correlation_heatmap` |      86 |
| Fig. 3 | `Fig3_blup_ranking`             |      86 |
| Fig. 4 | `Fig4_establishment_cohorts`    |      86 |

### Supplementary figures

| Figure  | File stem              | Dataset |
| ------- | ---------------------- | ------: |
| Fig. S1 | `FigS1_ammi_biplot`    |      43 |
| Fig. S2 | `FigS2_gge_biplot`     |      43 |
| Fig. S3 | `FigS3_mean_stability` |      43 |

Figures are written to:

```text
04_Results/Manuscript_Figures/
```

PNG output is generated at **300 dpi**. PDF output uses the Cairo PDF device.

The current figure script uses **Arial** as the base font for consistent manuscript typography.

---

## 16. Reproducible Execution

### 16.1 Main 86-genotype analysis

Run:

```bash
Rscript analysis_86.R
```

The main 86-genotype tabular results are written to:

```text
04_Results/Genotype86/
```

### 16.2 Existing 43-genotype analysis workflow

The repository also retains `analysis_main.R` for the 43-genotype analytical workflow used elsewhere in the project. It should be treated separately from the main 86-genotype manuscript mixed-model analysis described above.

### 16.3 Manuscript figure generation

Run:

```bash
Rscript make_figures.R
```

The script reads the four annual raw phenotype files, performs the required preprocessing, fits the 86-genotype mixed models, constructs the 43-genotype balanced subset for AMMI/GGE, and writes:

```text
04_Results/Manuscript_Figures/
```

### 16.4 Supplementary reproducibility workflow

Place the anonymized supplementary data beside the supplementary script:

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

The supplementary script is designed to resolve its own directory and does not require a machine-specific absolute path.

---

## 17. Software Requirements

Recommended R version:

```text
R >= 4.3
```

### Main analytical packages

```r
install.packages(c(
  "lme4",
  "dplyr",
  "tidyr"
))
```

### Figure-generation packages

```r
install.packages(c(
  "lme4",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "patchwork"
))
```

AMMI and GGE decomposition are implemented with base R `svd()` and therefore do not require specialized stability-analysis packages such as `agricolae` or `metan`.

---

## 18. Path Handling and Portability

The main analysis and figure-generation scripts currently use the project root defined in the script configuration. For a different local installation, update the root path:

```r
ROOT <- "C:/Users/lijia/Documents/R Workplace/Sino-Australian_Alfalfa_Project"
```

The supplementary workflow is designed separately for redistribution and automatically resolves the directory in which the supplementary script is located.

When publishing or sharing the repository, user-specific absolute paths should be replaced with relative or project-root-aware paths where appropriate.

---

## 19. Relationship Between the 86- and 43-Genotype Analyses

The two genotype sets have distinct analytical purposes.

### 86-genotype analysis

The 86-genotype dataset retains all available genotypes and maximizes the available information for mixed-model estimation.

It is the basis for:

- variance components;
- broad-sense heritability;
- BLUPs;
- BLUP correlations;
- genotype ranking; and
- the main manuscript figures.

### 43-genotype analysis

The 43-genotype dataset is a balanced subset with complete records for all traits and all four years.

It is retained for:

- AMMI;
- GGE;
- ASV;
- GGE distance; and
- supplementary exploratory figures.

The two analyses should therefore not be described as estimates from the same analytical population. Their roles are complementary.

---

## 20. Interpretation Notes

This repository is intended to reproduce the computational analyses underlying the manuscript.

In particular:

- broad-sense heritability should not be interpreted as narrow-sense heritability or additive genetic variance;
- BLUPs are model-based genotype-level predictions rather than raw phenotype means;
- the main mixed-model analysis uses all 86 available genotypes;
- AMMI and GGE are restricted to the 43 complete genotypes because a balanced genotype × year matrix is used;
- stability interpretation is exploratory because only four years/environments are available;
- the establishment-cohort analysis uses the 2003–2005 fair-comparison window; and
- the conditional 2002 FW/DW correction is an explicit preprocessing step and should be retained when reproducing the reported results.

---

## 21. Version

Current release:

```text
v2.1.0
```

This release updates the repository to the manuscript-aligned figure architecture and analysis roles described above.
