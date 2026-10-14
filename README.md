# Research Compendium for 'Exploring time, distance and climate as drivers of technological variability in the Italian Aurignacian'

**Compendium DOI**

[![DOI](https://zenodo.org/badge/DOI/XX.XXXX/zenodo.XXXXXXX.svg)](https://doi.org/XX.XXXX/zenodo.XXXXXXX)

The content available at the provided DOI reproduces the results as documented in the publication. The files hosted at https://github.com/ArmandoFalcucci/FLEXDIST-Aurignacian represent the developmental versions and might have undergone modifications since the paper's publication.

## Maintainer of this Repository

[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--XXXX--XXXX-green.svg)](https://orcid.org/0000-0002-XXXX-XXXX) Armando Falcucci (A.Falcucci@soton.ac.uk)

## Paper in preparation

Falcucci et al., Time, space, and environment: technological variability in the Italian Aurignacian and the limits of matrix statistics. *Journal to be defined*.

## Abstract

<!-- REPLACE with the final manuscript abstract before deposit. The version
     below is the pre-revision text and no longer matches the paper. -->

The Aurignacian, long treated as the archaeological signature of early *Homo sapiens* in Europe, is usually divided into chrono-cultural phases whose spatial and temporal structure remains poorly quantified. We introduce a data-driven workflow to test how geography, time, and environment structure technological variability, applying it to six assemblages from five stratified sequences south of the Alps and across peninsular Italy dated between approximately 42,000 and 34,000 years ago. Lithic assemblages are summarized through a mixed-data dimensionality reduction technique, and the resulting technological distances are compared with least-cost path, chronological, and paleoclimatic distances using matrix correlation statistics with leave-one-out sensitivity analysis. Our approach independently recovers a Protoaurignacian–Early Aurignacian distinction consistent with previous site-based assessments, structured primarily by chronology rather than by geographic proximity. Apparent correlations between technology and environmental variables, principally annual precipitation, are shown to be driven by a single outlying assemblage and by uneven geographic sampling, collapsing once these are accounted for. We argue that matrix correlation statistics cannot support confirmatory inference about competing drivers at the sample sizes typical of regional Paleolithic studies, and we outline diagnostics that guard against over-interpretation. These results provide a template, and a cautionary framework, for wider analyses of technological change at the onset of the Upper Paleolithic.

## Dataset Attribution

The lithic datasets analysed in this research come from the renewed technological analyses of the studied sequences. Each site-based dataset is published open-access under a CC-BY-4.0 license in its own repository:

- **Grotta della Cala**: Falcucci, Kitagawa, et al. (2025). Available at: https://doi.org/10.5281/zenodo.15430431.
- **Grotta di Castelcivita**: Falcucci, Arrighi, et al. (2024). Available at: https://doi.org/10.5281/zenodo.10639552.
- **Grotta di Fumane**: Falcucci, Giusti, et al. (2024). Available at: https://doi.org/10.5281/zenodo.10965413.
- **Riparo Bombrini**: Falcucci, Bertola, et al. (2025). Available at: https://doi.org/10.5281/zenodo.15363594.
- **Grotta Paglicci**: Analyzed by A. Falcucci in this paper. Data available in this compendium.

Paleoclimatic reconstructions were extracted from the Krapp et al. (2021) global model via the `pastclim` package (Leonardi et al., 2023). The downscaled model shipped here was produced with WorldClim 2.1 as the high-resolution reference.

## Overview of Contents

- **data**: blank and core attribute datasets, the FLEXDIST-ready subsets, the derived distance matrices, and the precomputed downscaled palaeoclimate model (`Europe_bioclims_10min.nc`).
- **script**: R scripts in run order, plus the Supplementary Material `.Rmd`.
- **output**: `figures/`, `bayesian-models/` (OxCal), and `rds/` (intermediate caches passed between scripts).

## How to Reproduce

The scripts are numbered in run order and pass objects to one another through `.rds` caches in `output/rds/`, so each runs standalone once the previous one has completed.

```
script/00a_flexdist-blanks.R   # FLEXDIST of blank attributes
                               #   -> data/distance_matrix_blanks.csv
                               #   -> output/rds/04_scores_blanks.rds
script/00b_flexdist-cores.R    # FLEXDIST of core attributes
                               #   -> data/distance_matrix_cores.csv
                               #   -> output/rds/04_scores_cores.rds
script/01_distances.R          # coordinates, dates, climate and least-cost-path
                               #   predictors -> output/rds/01_distances.rds
script/02_mantel.R             # combined lithic distance, simple and partial
                               #   Mantel tests, leave-one-out
                               #   -> output/rds/02_mantel.rds ; Figures 6-7
script/03_figures_tables.R     # assemblage comparability: Table 2, Figures 2-3
                               #   (reads the raw datasets only)
script/04_null_diagnostics.R   # exact family-wise permutation test, assemblage
                               #   size check, specimen-level bootstrap
                               #   -> output/rds/04_diagnostics.rds ; Figure S5
script/supplementary-material.Rmd   # knit last; reads the 01, 02 and 04 caches
```

`03_figures_tables.R` is independent of scripts 01, 02 and 04 and can be run at any time.

Steps:

1. Download the entire repository.
2. Open the `FLEXDIST-Aurignacian.Rproj` project file in a **fresh R session**, with no workspace restored from a previous session.
3. Confirm that `data/Europe_bioclims_10min.nc` is present.
4. Run the scripts in the order above, then knit `script/supplementary-material.Rmd`.

The downscaled paleoclimate model is provided so that it need not be rebuilt. To regenerate it from scratch, set `RUN_DOWNSCALING <- TRUE` at the top of `script/01_distances.R` (slow; requires large WorldClim and Krapp2021 downloads).

### Notes on the analysis

- **Group distances.** FLEXDIST (Rathmann et al., 2023) returns distances between individual artifacts. Scripts `00a` and `00b` extend this to groups: the group distance matrix is computed once per imputation run and summarised as the cell-wise median across runs, mirroring the specimen-level construction. Group positions are the coordinatewise median of the standardized PC scores.
- **Reproducibility of the random components.** The imputation and the bootstrap draw random numbers; each script sets its seed at the top. Permutation results are exact rather than sampled, because at six assemblages all 720 relabellings are enumerated in full, and so do not depend on the seed.
- **Runtime.** `00a` is the slowest step, since the ordination is repeated across imputation runs on 2,530 blanks. `01_distances.R` is slow on first use if the climate model is regenerated.

For consistent results, the `renv` package (v.1.1.5) was used, following the procedures outlined in its vignette. All analyses were performed using R 4.5.1. Necessary packages are recorded in `renv.lock`.

## Licenses

**Code**: MIT License (https://opensource.org/licenses/MIT), copyright holders: Armando Falcucci (2026).

**Data and Intellectual Work**: Creative Commons Attribution 4.0 International License (https://creativecommons.org/licenses/by/4.0/), copyright holder: the authors (2026).
