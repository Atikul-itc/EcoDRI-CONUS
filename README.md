# EcoDRI-CONUS

Google Earth Engine code for producing the weekly Ecological Drought Index (EcoDRI) dataset for the conterminous United States (2001-2024).

**Companion paper:** Islam, A. *et al.* A weekly multi-source Ecological Drought Index (EcoDRI) for the conterminous United States, 2001-2024. *Scientific Data* (in review).

**Published dataset:** [Zenodo DOI 10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx) (to be assigned upon publication).

**Interactive viewer:** [EcoDRI Earth Engine app](https://ee-atikul.projects.earthengine.app/view/ecodri-conus).

## Repository contents

All scripts are Google Earth Engine JavaScript, intended to be pasted into the [GEE Code Editor](https://code.earthengine.google.com/) or imported as a repository. The scripts implement the computational pipeline described in Sections 2 (Methods) and 4 (Technical Validation) of the paper.

| File | Purpose | Corresponds to |
|---|---|---|
| `gee/EcoDRI_Weekly_Composite.js`     | Generates a single weekly EcoDRI composite: reads MODIS, GLDAS/SMAP, GRIDMET inputs; computes VCI, TCI, SMCI; applies Köppen-weighted aggregation; classifies into six ordinal categories. | Methods §2.3-§2.6 |
| `gee/EcoDRI_Batch_Export.js`          | Driver script that submits weekly export tasks for the full 2001-2024 study period.                                                                                                       | Methods §2.5, §2.8 |
| `gee/EcoDRI_USDM_Validation_Sampling.js` | Draws stratified pixel samples paired with contemporaneous USDM classifications; produces the CSVs underlying the kappa and per-class metrics.                                          | Methods §2.7; Results §4.2 |
| `gee/EcoDRI_County_Aggregation.js`   | Aggregates weekly EcoDRI to US counties over the full growing season and crop-specific critical windows.                                                                                  | Methods §2.7; Results §4.3 |
| `gee/EcoDRI_Ecoregion_Aggregation.js` | Aggregates weekly EcoDRI to EPA Level III ecoregions (rangeland-dominated only), for pairing with RAP herbaceous biomass.                                                                | Methods §2.7; Results §4.1 |
| `gee/RAP_Ecoregion_Aggregation.js`   | Aggregates the Rangeland Analysis Platform 16-day NPP product to matching ecoregions for the ecological validation.                                                                       | Methods §2.7; Results §4.1 |

## How to run

1. Sign in to the GEE Code Editor.
2. Create a new repository or copy the scripts into your workspace.
3. Update the asset path prefix at the top of each script (currently `projects/ee-atikul/...`) to point at your own asset location.
4. Run `EcoDRI_Weekly_Composite.js` on a single date to inspect the output. Once you are satisfied, run `EcoDRI_Batch_Export.js` to submit the full 2001-2024 archive.
5. After the weekly assets are built, run the validation-sampling and aggregation scripts to reproduce the paired CSVs used in the validation analyses.

## Prerequisites

- A GEE account with sufficient asset storage (~30 GB for the weekly EcoDRI collection).
- Pre-computed week-of-year percentile climatology assets for NDVI, LST, and soil moisture. These are used by `EcoDRI_Weekly_Composite.js` (see the `CLIM_*` asset paths near the top).
- The user-provided Köppen-Geiger classification shapefile as an asset (`Aridity_Koppen`).

## Downstream analysis (R)

The R analysis scripts that produce the specific validation numbers reported in the paper (kappa 0.44, per-crop yield correlations, per-ecoregion RAP correlations) are available on request. Contact the corresponding author.

## Data availability

The weekly EcoDRI composites (2001-2024) are archived at [Zenodo DOI 10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx) as multi-band GeoTIFFs. They are also mirrored in the Google Earth Engine asset collection `projects/ee-atikul/assets/EcoDRI_v3_2_growing_season`, accessible via the GEE catalog once the dataset receives its DOI.

## License

Code is released under the [MIT License](LICENSE). The EcoDRI dataset itself is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Citation

If you use EcoDRI in your research, please cite the companion paper and the Zenodo archive:

```
Islam, A., Buenemann, M., et al. (2024). A weekly multi-source Ecological Drought
Index (EcoDRI) for the conterminous United States, 2001-2024. Scientific Data.

Islam, A. et al. (2024). EcoDRI-CONUS: Weekly Ecological Drought Index for the
conterminous United States, 2001-2024 [Data set]. Zenodo. https://doi.org/10.xxxx/xxxxx
```

## Contact

Atikul Islam (`atikul` [at] `nmsu.edu`)
Buenemann Lab, Department of Geography and Environmental Studies
New Mexico State University
