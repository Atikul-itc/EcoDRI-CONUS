# R Analysis

R scripts that reproduce the numerical results reported in Islam et al. (2024) from the CSV outputs of the GEE pipeline.

## Scripts

| Script | Reproduces | Requires |
|---|---|---|
| `USDM_validation_kappa.R` | kappa metrics | Paired-sample CSVs from `gee/EcoDRI_USDM_Validation_Sampling.js` |
| `threshold_recalibration.R` | threshold values | Same paired-sample CSVs |
| `yield_validation.R` | crop correlations | County CSVs from `gee/EcoDRI_County_Aggregation.js` + NASS API access |
| `ecological_validation.R` | primary ecological result | Ecoregion CSVs from `gee/EcoDRI_Ecoregion_Aggregation.js` and `gee/RAP_Ecoregion_Aggregation.js` |
| `ablation_analysis.R` | weighting-scheme comparison | Paired-sample CSVs (with raw components) |

## Prerequisites

Tested on R 4.3 and later. Install dependencies once:

```r
install.packages(c("tidyverse", "irr", "caret", "nloptr", "lme4", "rnassqs"))
```

For yield validation, register for a [USDA NASS API token](https://quickstats.nass.usda.gov/api) and set the environment variable before running:

```bash
export NASS_QS_TOKEN="your-token-here"
```

## Expected inputs

Each script expects the corresponding GEE-output CSV folder to be present in the working directory. Folder names default to:

- `EcoDRI_USDM_Validation/` (paired samples)
- `EcoDRI_County_Aggregation/` (county summaries)
- `EcoDRI_Ecoregion_Aggregation/` (EcoDRI at ecoregions)
- `RAP_Ecoregion_Aggregation/` (RAP NPP at ecoregions)

Edit the `INPUT_DIR` variable at the top of each script if your folders are elsewhere.

## Outputs

Each script writes to a script-specific `outputs_*/` folder: CSV tables and a headline summary text file. These reproduce the numbers in the corresponding paper section.

## Reproducing the paper end-to-end

```bash
# In R, working directory containing all four GEE-output folders:
Rscript r_analysis/threshold_recalibration.R     # ~2 min
Rscript r_analysis/USDM_validation_kappa.R       # ~1 min
Rscript r_analysis/yield_validation.R            # ~10 min (NASS downloads)
Rscript r_analysis/ecological_validation.R       # ~1 min
Rscript r_analysis/ablation_analysis.R           # ~15 min (optimization)
```
