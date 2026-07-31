
# Multi-year stability analysis of Sino-Australian alfalfa cultivars

This repository contains R scripts for the statistical analysis and visualization described in:

**"Multi-year evaluation of Sino-Australian alfalfa cultivars in Gansu: heritability, BLUP, and GGE/AMMI stability analysis"**

The workflow evaluates multi-environment performance, genetic variation, and stability of alfalfa cultivars using mixed models, BLUP prediction, broad-sense heritability estimation, AMMI and GGE biplot analyses.

---

## Overview

The analysis pipeline includes:

1. Raw data import and harmonization
2. Phenotypic data cleaning and genotype × environment matrix construction
3. Linear mixed model fitting
4. Genotypic BLUP estimation
5. Broad-sense heritability calculation
6. AMMI stability analysis
7. GGE biplot analysis
8. Multi-trait clustering of genotypes
9. Automated table and figure generation

---

## Repository structure

```

00_Scripts/
Main R scripts

01_Data/
User-provided raw phenotypic data
(not included due to data availability restrictions)

02_Output/
Generated tables and figures
(may updates in published paper)

03_Documentation/
Additional workflow information
(not included due to data availability restrictions)

````

---

## Reproducible workflow

### Requirements

R version ≥ 4.3

Required packages:

- tidyverse
- lme4
- lmerTest
- ggplot2
- ggpubr
- ggrepel
- cluster
- ggdendro


Install packages:

```r
install.packages(c(
"tidyverse",
"lme4",
"lmerTest",
"ggrepel",
"cluster",
"ggdendro"
))
````

---

## Running the analysis

Place the raw phenotype files in:

```
01_Data/raw/
```

Then run:

```r
source("00_Scripts/00_run_pipeline.R")
```

The pipeline automatically generates:

* mixed model statistics
* variance components
* BLUP estimates
* heritability estimates
* AMMI outputs
* GGE coordinates
* clustering results
* manuscript figures

---

## Statistical methods

### Linear mixed model

For each trait, a linear mixed model was fitted:

$$
Y = \mu + Year + Genotype + Genotype \times Year + Rep(Year) + e
$$

where:

- $Year$ was treated as a fixed environmental effect.
- $Genotype$ and $Genotype \times Year$ interaction were treated as random effects.
- $Rep(Year)$ represents replication nested within year.


---

### BLUP prediction

Genotypic values were estimated using best linear unbiased prediction (BLUP):

$$
BLUP_g = \mu + \hat{g}
$$

where:

- $\mu$ represents the population mean.
- $\hat{g}$ represents the predicted genotypic effect.


---

### Broad-sense heritability

Broad-sense heritability ($H^2$) across environments was estimated as:

$$
H^2 =
\frac{\sigma_G^2}
{\sigma_G^2+\frac{\sigma_{GE}^2}{n_e}
+\frac{\sigma_e^2}{n_e n_r}}
$$

where:

- $\sigma_G^2$ = genotypic variance component.
- $\sigma_{GE}^2$ = genotype × environment interaction variance.
- $\sigma_e^2$ = residual variance.
- $n_e$ = number of environments (years).
- $n_r$ = mean number of replications per genotype within environment.


---

### Stability analysis

Genotype × environment interaction was evaluated using:

- Additive Main Effects and Multiplicative Interaction (AMMI) analysis.
- Genotype plus Genotype-by-Environment (GGE) biplot analysis.

---

## Data availability

Raw phenotypic data are not included in this repository.Data access is available from the corresponding author upon reasonable request.
