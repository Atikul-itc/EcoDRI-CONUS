# EcoDRI-CONUS

Google Earth Engine and R analysis code for the weekly Ecological Drought Index (EcoDRI) dataset covering the conterminous United States, 2001-2024.

**Citation:** Hoque, A. *et al.* A weekly multi-source Ecological Drought Index (EcoDRI) for the conterminous United States, 2001-2024. *Scientific Data* (in review).

**Published dataset:** [Zenodo DOI [10.5281/zenodo.22119391](https://doi.org/10.5281/zenodo.22119391)].

**Interactive viewer:** [EcoDRI Earth Engine app](https://ee-atikul.projects.earthengine.app/view/ecodri-conus).

## Repository structure

```
EcoDRI-CONUS/
├── gee/               # Earth Engine JavaScript pipeline
│   ├── EcoDRI_Weekly_Composite.js
│   ├── EcoDRI_Batch_Export.js
│   ├── EcoDRI_USDM_Validation_Sampling.js
│   ├── EcoDRI_County_Aggregation.js
│   ├── EcoDRI_Ecoregion_Aggregation.js
│   └── RAP_Ecoregion_Aggregation.js
├── r_analysis/        # R scripts reproducing the paper's numerical results
│   ├── USDM_validation_kappa.R
│   ├── threshold_recalibration.R
│   ├── yield_validation.R
│   ├── ecological_validation.R
│   ├── ablation_analysis.R
│   └── README.md
├── LICENSE
├── CITATION.cff
└── README.md
```

## GEE pipeline (Methods §2, §4)

| File | Purpose | Corresponds to |
|---|---|---|
| `gee/EcoDRI_Weekly_Composite.js`     | Generates a single weekly EcoDRI composite from MODIS, GLDAS/SMAP, and GRIDMET inputs; applies Köppen-weighted aggregation; classifies into six ordinal categories. | Methods §2.3-§2.6 |
| `gee/EcoDRI_Batch_Export.js`          | Driver script that submits weekly export tasks for the full 2001-2024 archive.                                                                                    | Methods §2.5, §2.8 |
| `gee/EcoDRI_USDM_Validation_Sampling.js` | Draws stratified pixel samples paired with USDM classifications; produces the CSVs used in the kappa and per-class analyses.                                    | Methods §2.7; §4.2 |
| `gee/EcoDRI_County_Aggregation.js`   | Aggregates weekly EcoDRI to US counties over the full growing season and crop-specific critical windows.                                                          | Methods §2.7; §4.3 |
| `gee/EcoDRI_Ecoregion_Aggregation.js` | Aggregates weekly EcoDRI to EPA Level III ecoregions (rangeland-dominated only) for pairing with RAP biomass.                                                    | Methods §2.7; §4.1 |
| `gee/RAP_Ecoregion_Aggregation.js`   | Aggregates the Rangeland Analysis Platform 16-day NPP product to matching ecoregions.                                                                             | Methods §2.7; §4.1 |

## R analysis (reproduces numerical results)

Five scripts reproduce the numbers in the paper from the CSV outputs of the GEE pipeline. See `r_analysis/README.md` for details.

## How to run

### Step 1: Weekly EcoDRI generation (GEE)

1. Sign in to the [GEE Code Editor](https://code.earthengine.google.com/).
2. Update the asset path prefix in each script from `projects/ee-atikul/...` to your own asset location.
3. Run `EcoDRI_Weekly_Composite.js` on one date to inspect the output.
4. Run `EcoDRI_Batch_Export.js` to build the full 2001-2024 archive (~720 export tasks; processes over ~24 hours).

### Step 2: Validation sampling and aggregation (GEE)

Once the weekly assets are built, run the four validation-sampling and aggregation scripts. They produce CSVs downloaded to your Google Drive.

### Step 3: Analysis (R)

Download the CSVs to a local working directory. Install R dependencies:

```r
install.packages(c("tidyverse", "irr", "caret", "nloptr", "lme4", "rnassqs"))
```

Run each R script from `r_analysis/`. Output tables and headline summaries land in `outputs_*/` folders.

## Prerequisites

- A Google Earth Engine account with sufficient asset storage (~30 GB for the weekly collection).
- Pre-computed week-of-year percentile climatology assets (referenced as `CLIM_*` at the top of `EcoDRI_Weekly_Composite.js`).
- A user-provided Köppen-Geiger classification shapefile as an asset (`Aridity_Koppen`).
- R 4.3 or later.
- A [USDA NASS API token](https://quickstats.nass.usda.gov/api) for the yield download step.

## Data availability

The weekly EcoDRI composites (2001-2024) are archived at [Zenodo DOI 10.5281/zenodo.22119391](https://zenodo.org/records/22119391) as multi-band GeoTIFFs.

## License

Code is released under the [MIT License](LICENSE). The EcoDRI dataset itself is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).


