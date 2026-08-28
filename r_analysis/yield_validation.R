# =============================================================================
# Agricultural Yield Validation
#
# Reads USDA NASS county yields, applies LOWESS detrending, joins with
# county-aggregated EcoDRI summaries, and computes correlations for corn,
# soybeans, and winter wheat over full-season and phenology-aligned
# critical windows.
#
# Reproduces the r = -0.45, -0.36, -0.11 corn/soy/wheat correlations
# reported in Section 4.3 and the Full/Critical comparison shown in Figure 4b.
# =============================================================================

library(tidyverse)
library(rnassqs)     # USDA NASS API client

# ----- CONFIGURATION ---------------------------------------------------------
COUNTY_DIR  <- "EcoDRI_County_Aggregation/"    # from GEE county aggregation
OUTPUT_DIR  <- "outputs_yield_validation/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# NASS API token: get from https://quickstats.nass.usda.gov/api
# then set as env variable NASS_QS_TOKEN before running.
nass_token <- Sys.getenv("NASS_QS_TOKEN")
if (nass_token == "") stop("Set NASS_QS_TOKEN environment variable.")
nassqs_auth(nass_token)

YEARS <- 2010:2024
CROPS <- list(
  corn         = list(commodity = "CORN",     window_label = "corn"),
  soybean      = list(commodity = "SOYBEANS", window_label = "soybean"),
  winter_wheat = list(commodity = "WHEAT",    window_label = "winter_wheat",
                       class_desc = "WINTER")
)

# ----- DOWNLOAD YIELDS FROM NASS ---------------------------------------------
download_crop <- function(commodity, class_desc = NULL) {
  params <- list(
    commodity_desc = commodity,
    statisticcat_desc = "YIELD",
    unit_desc = "BU / ACRE",
    agg_level_desc = "COUNTY",
    year__GE = min(YEARS), year__LE = max(YEARS)
  )
  if (!is.null(class_desc)) params$class_desc <- class_desc
  nassqs(params) |>
    transmute(
      GEOID = paste0(state_fips_code, county_code),
      year = as.integer(year),
      yield = as.numeric(Value)
    ) |>
    filter(!is.na(yield))
}

cat("Downloading NASS yields...\n")
yields <- imap_dfr(CROPS, function(spec, crop_name) {
  cat("  ", crop_name, "...\n")
  df <- download_crop(spec$commodity, spec$class_desc)
  df$crop <- crop_name
  df
})

# ----- LOWESS DETRENDING (per county-crop) -----------------------------------
detrend_yields <- yields |>
  group_by(GEOID, crop) |>
  filter(n() >= 8) |>
  group_modify(~ {
    fit <- lowess(.x$year, .x$yield, f = 0.75)
    trend <- approx(fit$x, fit$y, xout = .x$year, rule = 2)$y
    resid <- .x$yield - trend
    tibble(
      year = .x$year,
      yield = .x$yield,
      trend = trend,
      residual = resid,
      yield_zscore = resid / sd(resid, na.rm = TRUE)
    )
  }) |>
  ungroup()

# ----- LOAD COUNTY-AGGREGATED EcoDRI -----------------------------------------
county_files <- list.files(COUNTY_DIR, pattern = "EcoDRI_county_.*\\.csv$",
                            full.names = TRUE)
county_df <- map_dfr(county_files, ~ read_csv(.x, col_types = cols(
  GEOID = col_character(), STATEFP = col_character(), NAME = col_character(),
  year = col_integer(), window = col_character(),
  EcoDRI_mean = col_double(), EcoDRI_max = col_double(),
  EcoDRI_p90 = col_double(),
  .default = col_double()
))) |>
  mutate(GEOID = str_pad(GEOID, 5, pad = "0"))

# ----- MERGE AND CORRELATE ---------------------------------------------------
compute_crop_correlations <- function(crop_name, window_label) {
  window_df <- county_df |> filter(window == window_label)
  merged <- inner_join(
    detrend_yields |> filter(crop == crop_name),
    window_df, by = c("GEOID", "year")
  )
  if (nrow(merged) < 30) return(NULL)
  tibble(
    crop = crop_name,
    window = window_label,
    n = nrow(merged),
    pearson_r = cor(merged$yield_zscore, merged$EcoDRI_mean),
    spearman_r = cor(merged$yield_zscore, merged$EcoDRI_mean,
                     method = "spearman"),
    p_value = cor.test(merged$yield_zscore, merged$EcoDRI_mean)$p.value
  )
}

results <- expand_grid(
  crop = names(CROPS),
  window = c("full_season", "corn", "soybean", "winter_wheat")
) |>
  filter(window == "full_season" | window == crop) |>
  pmap_dfr(compute_crop_correlations)

cat("\n=== Yield correlations (Section 4.3, Figure 4b) ===\n")
print(results)
write_csv(results, file.path(OUTPUT_DIR, "yield_correlations.csv"))

cat("\nOutputs written to:", OUTPUT_DIR, "\n")
