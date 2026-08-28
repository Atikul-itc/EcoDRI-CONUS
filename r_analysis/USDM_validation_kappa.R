# =============================================================================
# USDM Validation: Cohen's Kappa and Per-Class Metrics
#
# Reads the paired EcoDRI-USDM CSV files produced by
# gee/EcoDRI_USDM_Validation_Sampling.js and computes:
#   - Overall and per-fold kappa (unweighted, linear, quadratic)
#   - Per-class sensitivity, specificity, F1
#   - Spearman correlation between continuous EcoDRI and USDM categories
#   - Confusion matrix (column proportions)
#
# Reproduces the numbers reported in Section 4.2 and Tables 4-5 of
# Islam et al. (2024).
# =============================================================================

library(tidyverse)
library(irr)         # for kappa2()
library(caret)       # for confusionMatrix() per-class metrics

# ----- CONFIGURATION ---------------------------------------------------------
INPUT_DIR   <- "EcoDRI_USDM_Validation/"    # folder of downloaded CSVs
OUTPUT_DIR  <- "outputs_usdm_validation/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

TRAIN_YEARS <- 2010:2019
TEST_YEARS  <- 2020:2024

# Recalibrated thresholds from Section 2.6 (Table 2, "Recalibrated" column)
THRESHOLDS <- c(1.50, 1.85, 2.40, 2.90, 3.15)

# ----- LOAD ------------------------------------------------------------------
files <- list.files(INPUT_DIR, pattern = "EcoDRI_USDM_paired_.*\\.csv$",
                    full.names = TRUE)
stopifnot(length(files) > 0)

df <- map_dfr(files, ~ read_csv(.x, col_types = cols(
  date = col_character(),
  koppen = col_integer(),
  usdm_cat = col_integer(),
  drought_category = col_integer(),
  EcoDRI = col_double(),
  VCI = col_double(), TCI = col_double(),
  SMCI = col_double(), SPI = col_double()
))) |>
  mutate(date = as.Date(date),
         year = as.integer(format(date, "%Y")),
         fold = if_else(year %in% TRAIN_YEARS, "train", "test")) |>
  drop_na(EcoDRI, usdm_cat, drought_category)

cat("Loaded", nrow(df), "paired observations across",
    length(unique(df$date)), "dates.\n")
cat("Train fold:", sum(df$fold == "train"), "| Test fold:",
    sum(df$fold == "test"), "\n\n")

# ----- KAPPA HELPERS ---------------------------------------------------------
compute_kappa_metrics <- function(pred, ref) {
  pred_f <- factor(pred, levels = 0:5)
  ref_f  <- factor(ref,  levels = 0:5)
  tibble(
    n                    = length(pred),
    accuracy             = mean(pred == ref),
    kappa_unweighted     = kappa2(cbind(pred, ref))$value,
    kappa_linear         = kappa2(cbind(pred, ref), weight = "equal")$value,
    kappa_quadratic      = kappa2(cbind(pred, ref), weight = "squared")$value,
    spearman             = cor(pred, ref, method = "spearman")
  )
}

# ----- PER-FOLD METRICS ------------------------------------------------------
per_fold <- df |>
  group_by(fold) |>
  group_modify(~ compute_kappa_metrics(.x$drought_category, .x$usdm_cat)) |>
  ungroup()

cat("=== Per-fold agreement metrics (Table 4) ===\n")
print(per_fold)
write_csv(per_fold, file.path(OUTPUT_DIR, "table4_usdm_kappa_by_fold.csv"))

# ----- PER-CLASS METRICS (test fold only, per Table 5) -----------------------
test_df <- df |> filter(fold == "test")

cm_test <- confusionMatrix(
  factor(test_df$drought_category, levels = 0:5),
  factor(test_df$usdm_cat, levels = 0:5)
)

per_class <- as.data.frame(cm_test$byClass) |>
  rownames_to_column("class") |>
  mutate(class = str_replace(class, "Class: ", "")) |>
  mutate(class_label = case_when(
    class == "0" ~ "None",         class == "1" ~ "D0 Abnormal",
    class == "2" ~ "D1 Moderate",  class == "3" ~ "D2 Severe",
    class == "4" ~ "D3 Extreme",   class == "5" ~ "D4 Exceptional"
  )) |>
  select(class_label, Sensitivity, Specificity, F1)

cat("\n=== Per-class metrics on test fold (Table 5) ===\n")
print(per_class)
write_csv(per_class, file.path(OUTPUT_DIR, "table5_per_class_metrics.csv"))

# ----- CONFUSION MATRIX (column proportions; Figure 3 source) ---------------
cm_prop <- prop.table(cm_test$table, margin = 2)
write_csv(as.data.frame.matrix(cm_prop) |> rownames_to_column("prediction"),
          file.path(OUTPUT_DIR, "figure3_confusion_matrix_column_props.csv"))

cat("\nOutputs written to:", OUTPUT_DIR, "\n")
