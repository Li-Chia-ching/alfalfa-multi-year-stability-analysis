# Cohort-Specific Alfalfa Visualization Pipeline

This script generates seven publication-ready figures from cohort-specific longitudinal analyses of historical *Medicago sativa* L. germplasm established in 2001 and 2002.

## Outputs

| Figure | Output                            | Main purpose                                                    |
| ------ | --------------------------------- | --------------------------------------------------------------- |
| Fig. 1 | `Fig1_same_stand_age_sensitivity` | Compare cohort performance at matched stand ages                |
| Fig. 2 | `Fig2_variance_heritability`      | Visualize variance components and broad-sense heritability      |
| Fig. 3 | `Fig3_BLUP_Caterpillar_v3`        | Rank genotype BLUPs by cohort and trait                         |
| Fig. 4 | `Fig4_ammi_biplot`                | Assess temporal stability using AMMI                            |
| Fig. 5 | `Fig5_gge_biplot`                 | Visualize genotype performance and G×year interaction using GGE |
| Fig. 6 | `Fig6_mean_stability`             | Relate mean performance to AMMI stability value (ASV)           |
| Fig. 7 | `Fig7_blup_correlation_cohort`    | Compare trait correlations between cohorts                      |

Figures are exported as PDF and PNG; Figure 3 also includes TIFF and input metadata.

## Methods

The pipeline visualizes results from linear mixed-effects models, including variance components, broad-sense heritability, and BLUPs. AMMI and GGE analyses use singular value decomposition (SVD), while ASV is used to summarize temporal stability. `adjustText` is used for label collision avoidance in biplots.

## Requirements

Python 3.8+ with:

```bash
pip install numpy pandas matplotlib adjustText
```

## Input

The script expects the following CSV files in the analysis-results directory:

```text
same_stand_age_sensitivity_final.csv
variance_components_final_by_cohort.csv
heritability_final_by_cohort.csv
BLUP_final_by_cohort.csv
plot_year_wide.csv
BLUP_correlation_by_cohort.csv
```

## Usage

Specify the input and output directories:

```bash
python make_all_figures_portable.py \
  --results-dir path/to/results \
  --output-dir path/to/figures
```

All output directories are created automatically.
