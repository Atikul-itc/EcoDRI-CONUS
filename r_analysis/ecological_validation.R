# =============================================================================
# Ecological Validation: EcoDRI vs RAP Herbaceous Biomass
#
# Reads ecoregion-aggregated EcoDRI and RAP NPP CSVs, applies per-ecoregion
# LOWESS detrending to biomass, computes per-ecoregion correlations, and
# summarizes by EPA Level I biome.
#
# Reproduces the primary ecological validation result: median per-ecoregion
# r = -0.66 within four water-limited rangeland biomes (Section 4.1).
# =============================================================================

library(tidyverse)
library(lme4)   # for mixed-effects model

# ----- CONFIGURATION ---------------------------------------------------------
ECODRI_DIR  <- "EcoDRI_Ecoregion_Aggregation/"
RAP_DIR     <- "RAP_Ecoregion_Aggregation/"
OUTPUT_DIR  <- "outputs_ecological_validation/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

WATER_LIMITED_BIOMES <- c(
  "GREAT PLAINS",
  "NORTH AMERICAN DESERTS",
  "NORTHWESTERN FORESTED MOUNTAINS",
  "SOUTHERN SEMI-ARID HIGHLANDS"
)

# ----- LOAD ------------------------------------------------------------------
ecodri <- list.files(ECODRI_DIR, pattern = "EcoDRI_ecoregion_peak_.*\\.csv$",
                      full.names = TRUE) |>
  map_dfr(~ read_csv(.x, col_types = cols(
    na_l3code = col_character(), na_l3name = col_character(),
    na_l1name = col_character(), year = col_integer(),
    .default = col_double()
  )))

rap <- list.files(RAP_DIR, pattern = "RAP_ecoregion_peak_.*\\.csv$",
                   full.names = TRUE) |>
  map_dfr(~ read_csv(.x, col_types = cols(
    na_l3code = col_character(), na_l3name = col_character(),
    na_l1name = col_character(), year = col_integer(),
    .default = col_double()
  )))

# ----- DETREND RAP HERBACEOUS NPP (per ecoregion) ---------------------------
detrend_zscore <- function(df, col) {
  if (nrow(df) < 8) return(df |> mutate("{col}_z" := NA_real_))
  y <- df[[col]]; x <- df$year
  keep <- !is.na(y) & is.finite(y)
  if (sum(keep) < 8) return(df |> mutate("{col}_z" := NA_real_))
  sm <- lowess(x[keep], y[keep], f = 0.75)
  trend <- approx(sm$x, sm$y, xout = x, rule = 2)$y
  resid <- y - trend
  s <- sd(resid, na.rm = TRUE)
  df[[paste0(col, "_z")]] <- if (!is.na(s) && s > 0) resid / s else NA_real_
  df
}

rap_z <- rap |>
  arrange(na_l3code, year) |>
  group_by(na_l3code) |>
  group_modify(~ detrend_zscore(.x, "NPP_herbaceous")) |>
  ungroup()

# ----- MERGE -----------------------------------------------------------------
merged <- inner_join(
  ecodri |> select(na_l3code, na_l3name, na_l1name, year, rangeland_frac,
                    EcoDRI_mean),
  rap_z |> select(na_l3code, year, NPP_herbaceous, NPP_herbaceous_z),
  by = c("na_l3code", "year")
) |> drop_na(EcoDRI_mean, NPP_herbaceous_z)

cat("Merged records:", nrow(merged),
    "across", n_distinct(merged$na_l3code), "ecoregions\n")

# ----- PER-ECOREGION PEARSON CORRELATIONS ------------------------------------
per_eco <- merged |>
  group_by(na_l3code, na_l3name, na_l1name) |>
  filter(n() >= 8) |>
  summarise(n = n(),
            r_rap = cor(EcoDRI_mean, NPP_herbaceous_z),
            .groups = "drop")

write_csv(per_eco, file.path(OUTPUT_DIR, "per_ecoregion_correlations.csv"))

# ----- PER-BIOME SUMMARY (Figure 7a, Table 6) --------------------------------
per_biome <- per_eco |>
  group_by(na_l1name) |>
  summarise(n_ecoregions = n(),
            median_r = median(r_rap, na.rm = TRUE),
            q25 = quantile(r_rap, 0.25, na.rm = TRUE),
            q75 = quantile(r_rap, 0.75, na.rm = TRUE)) |>
  arrange(median_r)

cat("\n=== Per-biome correlation summary (Figure 7a, Table 6) ===\n")
print(per_biome)
write_csv(per_biome, file.path(OUTPUT_DIR, "per_biome_summary.csv"))

# ----- WATER-LIMITED CORE RESULT ---------------------------------------------
core <- merged |>
  filter(toupper(na_l1name) %in% WATER_LIMITED_BIOMES)

core_median_r <- per_eco |>
  filter(toupper(na_l1name) %in% WATER_LIMITED_BIOMES) |>
  pull(r_rap) |>
  median(na.rm = TRUE)

cat(sprintf("\nWater-limited core: median per-ecoregion r = %.3f\n",
            core_median_r))

# ----- MIXED-EFFECTS MODEL (conditional R^2) ---------------------------------
mm <- lmer(
  NPP_herbaceous_z ~ EcoDRI_mean + (1 + EcoDRI_mean | na_l3code),
  data = core, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
var_fixed <- var(predict(mm, re.form = NA), na.rm = TRUE)
var_random <- as.numeric(VarCorr(mm)$na_l3code[1, 1])
var_resid <- attr(VarCorr(mm), "sc") ^ 2
r2_marg <- var_fixed / (var_fixed + var_random + var_resid)
r2_cond <- (var_fixed + var_random) / (var_fixed + var_random + var_resid)

cat(sprintf("Mixed-effects marginal R^2:    %.3f\n", r2_marg))
cat(sprintf("Mixed-effects conditional R^2: %.3f\n", r2_cond))

# ----- WRITE HEADLINE SUMMARY ------------------------------------------------
sink(file.path(OUTPUT_DIR, "headline_result.txt"))
cat("EcoDRI Ecological Validation: Headline Result\n")
cat("=============================================\n")
cat(sprintf("Median per-ecoregion Pearson r within 4 water-limited biomes: %.3f\n",
            core_median_r))
cat(sprintf("Mixed-effects conditional R^2 (core biomes):                  %.3f\n",
            r2_cond))
cat(sprintf("Interpretation: EcoDRI explains approximately %d percent of\n",
            round(core_median_r^2 * 100)))
cat("inter-annual herbaceous biomass variance within typical water-\n")
cat("limited rangeland ecoregions.\n")
sink()

cat("\nOutputs written to:", OUTPUT_DIR, "\n")
