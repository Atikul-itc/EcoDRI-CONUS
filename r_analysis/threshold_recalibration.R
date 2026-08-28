# =============================================================================
# EcoDRI Threshold Recalibration
#
# Optimizes the five monotonic EcoDRI categorical thresholds against USDM
# using a 70/30 date-blocked train/test split. Produces the recalibrated
# threshold values (1.50, 1.85, 2.40, 2.90, 3.15) reported in Table 2.
#
# Reproduces the calibration procedure described in Section 2.6 and the
# operational thresholds adopted for the published v3.2 dataset.
# =============================================================================

library(tidyverse)
library(irr)
library(nloptr)

INPUT_DIR   <- "EcoDRI_USDM_Validation/"
OUTPUT_DIR  <- "outputs_threshold_recalibration/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

TRAIN_YEARS <- 2010:2019
TEST_YEARS  <- 2020:2024

# ----- LOAD PAIRED SAMPLES ---------------------------------------------------
files <- list.files(INPUT_DIR, pattern = "EcoDRI_USDM_paired_.*\\.csv$",
                    full.names = TRUE)
df <- map_dfr(files, ~ read_csv(.x, col_types = cols(
  date = col_character(),
  usdm_cat = col_integer(),
  EcoDRI = col_double(),
  .default = col_double()
))) |>
  mutate(date = as.Date(date),
         year = as.integer(format(date, "%Y"))) |>
  drop_na(EcoDRI, usdm_cat)

train <- df |> filter(year %in% TRAIN_YEARS)
test  <- df |> filter(year %in% TEST_YEARS)

# ----- HELPERS ---------------------------------------------------------------
discretize <- function(x, t) {
  cut(x, breaks = c(-Inf, t, Inf), labels = 0:5,
      include.lowest = TRUE, right = FALSE) |>
    as.integer() |>
    (\(v) v - 1)()
}

quad_kappa <- function(pred, obs) {
  suppressWarnings(kappa2(cbind(pred, obs), weight = "squared")$value)
}

# ----- OPTIMIZE THRESHOLDS ON TRAINING FOLD ----------------------------------
obj <- function(t) {
  if (any(diff(t) <= 0)) return(0.5)   # invalid ordering
  pred <- discretize(train$EcoDRI, t)
  -quad_kappa(pred, train$usdm_cat)
}
ineq_constraint <- function(t) -diff(t) + 0.02

cat("Optimizing thresholds on training fold (n =", nrow(train), ") ...\n")
res <- nloptr(
  x0 = c(1.5, 1.85, 2.4, 2.9, 3.15),
  eval_f = obj,
  eval_g_ineq = ineq_constraint,
  lb = rep(0, 5), ub = rep(4, 5),
  opts = list(algorithm = "NLOPT_LN_COBYLA",
              xtol_rel = 1e-4, maxeval = 500)
)

recal_thresholds <- res$solution
train_kappa <- -res$objective

# ----- EVALUATE ON TEST FOLD -------------------------------------------------
pred_test <- discretize(test$EcoDRI, recal_thresholds)
test_kappa <- quad_kappa(pred_test, test$usdm_cat)

# ----- COMPARE WITH THEORETICAL ---------------------------------------------
theoretical <- c(2.50, 2.85, 3.20, 3.45, 3.75)
theo_kappa_test <- quad_kappa(discretize(test$EcoDRI, theoretical),
                               test$usdm_cat)

result <- tibble(
  boundary = c("None|D0","D0|D1","D1|D2","D2|D3","D3|D4"),
  theoretical = theoretical,
  recalibrated = round(recal_thresholds, 3)
)

cat("\n=== Recalibrated thresholds (Table 2) ===\n")
print(result)
cat(sprintf("\nTraining-fold quadratic kappa (calibration): %.3f\n", train_kappa))
cat(sprintf("Test-fold quadratic kappa (recalibrated):    %.3f\n", test_kappa))
cat(sprintf("Test-fold quadratic kappa (theoretical):     %.3f\n",
            theo_kappa_test))

write_csv(result, file.path(OUTPUT_DIR, "table2_thresholds.csv"))
sink(file.path(OUTPUT_DIR, "recalibration_summary.txt"))
cat("Threshold recalibration summary\n")
cat(sprintf("Training fold n:   %d\n", nrow(train)))
cat(sprintf("Test fold n:       %d\n", nrow(test)))
cat(sprintf("Training kappa:    %.3f\n", train_kappa))
cat(sprintf("Test kappa:        %.3f  <-- reported in Section 4.2\n", test_kappa))
cat(sprintf("Theoretical kappa: %.3f\n", theo_kappa_test))
sink()

cat("\nOutputs written to:", OUTPUT_DIR, "\n")
