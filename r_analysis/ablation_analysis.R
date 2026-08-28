# =============================================================================
# EcoDRI Weighting Scheme Ablation
#
# Compares four weighting schemes for the multi-source drought aggregation,
# each paired with its own recalibrated thresholds, and reports held-out
# quadratic-weighted Cohen's kappa on the test fold.
#
# Schemes:
#   1. Equal weights (0.25 each)
#   2. Global optimized (single weight vector for all pixels)
#   3. Theoretical Köppen (Table 1, v3.2 operational)
#   4. Optimized Köppen (per-zone under sum-to-1 non-negativity constraint)
#
# Reproduces the ablation results reported in Section 4.4 (Sensitivity
# and robustness).
# =============================================================================

library(tidyverse)
library(irr)
library(nloptr)

INPUT_DIR   <- "EcoDRI_USDM_Validation/"
OUTPUT_DIR  <- "outputs_ablation/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

TRAIN_YEARS <- 2010:2019
TEST_YEARS  <- 2020:2024

# ----- LOAD ------------------------------------------------------------------
files <- list.files(INPUT_DIR, pattern = "EcoDRI_USDM_paired_.*\\.csv$",
                    full.names = TRUE)
df <- map_dfr(files, ~ read_csv(.x, col_types = cols(
  date = col_character(),
  koppen = col_integer(),
  usdm_cat = col_integer(),
  VCI = col_double(), TCI = col_double(),
  SMCI = col_double(), SPI = col_double(),
  .default = col_double()
))) |>
  mutate(date = as.Date(date),
         year = as.integer(format(date, "%Y"))) |>
  drop_na(VCI, TCI, SMCI, SPI, usdm_cat, koppen)

train <- df |> filter(year %in% TRAIN_YEARS)
test  <- df |> filter(year %in% TEST_YEARS)

# ----- HELPERS ---------------------------------------------------------------
d_component <- function(dat) {
  tibble(
    D_VCI  = 1 - dat$VCI / 100,
    D_TCI  = 1 - dat$TCI / 100,
    D_SMCI = 1 - dat$SMCI / 100,
    D_SPI  = pmin(pmax((-dat$SPI + 3) / 6, 0), 1)
  )
}

compute_ecodri_global <- function(dat, w) {
  d <- d_component(dat)
  4 * (w[1]*d$D_VCI + w[2]*d$D_TCI + w[3]*d$D_SMCI + w[4]*d$D_SPI)
}

compute_ecodri_koppen <- function(dat, weights_by_koppen) {
  d <- d_component(dat)
  out <- numeric(nrow(dat))
  for (k in 1:5) {
    idx <- dat$koppen == k
    if (any(idx)) {
      w <- weights_by_koppen[k, ]
      out[idx] <- 4 * (w[1]*d$D_VCI[idx] + w[2]*d$D_TCI[idx] +
                       w[3]*d$D_SMCI[idx] + w[4]*d$D_SPI[idx])
    }
  }
  out
}

discretize <- function(x, t) {
  cut(x, breaks = c(-Inf, t, Inf), labels = 0:5,
      include.lowest = TRUE, right = FALSE) |>
    as.integer() |> (\(v) v - 1)()
}
quad_kappa <- function(pred, obs) {
  suppressWarnings(kappa2(cbind(pred, obs), weight = "squared")$value)
}

optimize_thresholds <- function(ecodri, usdm,
                                 init = c(1.5, 1.85, 2.4, 2.9, 3.15)) {
  obj <- function(t) {
    if (any(diff(t) <= 0)) return(0.5)
    -quad_kappa(discretize(ecodri, t), usdm)
  }
  res <- nloptr(x0 = init, eval_f = obj,
                eval_g_ineq = function(t) -diff(t) + 0.02,
                lb = rep(0, 5), ub = rep(4, 5),
                opts = list(algorithm = "NLOPT_LN_COBYLA",
                            xtol_rel = 1e-4, maxeval = 400))
  list(thresholds = res$solution, kappa = -res$objective)
}

theoretical_koppen <- function() {
  matrix(c(0.20,0.20,0.20,0.40,   # A
           0.15,0.40,0.15,0.30,   # B
           0.20,0.25,0.20,0.35,   # C
           0.25,0.25,0.20,0.30,   # D
           0.10,0.30,0.30,0.30),  # E
         nrow = 5, byrow = TRUE)
}

optimize_global_weights <- function(train_dat) {
  softmax <- function(x) exp(x) / sum(exp(x))
  obj <- function(raw) {
    w <- softmax(raw)
    ecodri <- compute_ecodri_global(train_dat, w)
    -optimize_thresholds(ecodri, train_dat$usdm_cat)$kappa
  }
  res <- nloptr(x0 = c(0,0,0,0), eval_f = obj,
                lb = rep(-3, 4), ub = rep(3, 4),
                opts = list(algorithm = "NLOPT_LN_COBYLA",
                            xtol_rel = 1e-3, maxeval = 150))
  softmax(res$solution)
}

optimize_koppen_weights <- function(train_dat, n_min = 500) {
  softmax <- function(x) exp(x) / sum(exp(x))
  weights <- theoretical_koppen()
  for (k in 1:5) {
    sub <- train_dat |> filter(koppen == k)
    if (nrow(sub) < n_min) next
    obj <- function(raw) {
      w <- softmax(raw)
      -optimize_thresholds(compute_ecodri_global(sub, w), sub$usdm_cat)$kappa
    }
    res <- nloptr(x0 = c(0,0,0,0), eval_f = obj,
                  lb = rep(-3, 4), ub = rep(3, 4),
                  opts = list(algorithm = "NLOPT_LN_COBYLA",
                              xtol_rel = 1e-3, maxeval = 120))
    weights[k, ] <- softmax(res$solution)
  }
  weights
}

# ----- RUN ALL FOUR SCHEMES --------------------------------------------------
results <- list()

# Scheme 1: Equal
w <- c(0.25, 0.25, 0.25, 0.25)
thr <- optimize_thresholds(compute_ecodri_global(train, w), train$usdm_cat)
kappa_test <- quad_kappa(
  discretize(compute_ecodri_global(test, w), thr$thresholds),
  test$usdm_cat)
results[[1]] <- list(scheme = "1_equal", kappa_train = thr$kappa,
                      kappa_test = kappa_test)

# Scheme 2: Global optimized
w <- optimize_global_weights(train)
thr <- optimize_thresholds(compute_ecodri_global(train, w), train$usdm_cat)
kappa_test <- quad_kappa(
  discretize(compute_ecodri_global(test, w), thr$thresholds),
  test$usdm_cat)
results[[2]] <- list(scheme = "2_global_optimized", kappa_train = thr$kappa,
                      kappa_test = kappa_test)

# Scheme 3: Theoretical Köppen
w <- theoretical_koppen()
thr <- optimize_thresholds(compute_ecodri_koppen(train, w), train$usdm_cat)
kappa_test <- quad_kappa(
  discretize(compute_ecodri_koppen(test, w), thr$thresholds),
  test$usdm_cat)
results[[3]] <- list(scheme = "3_theoretical_koppen", kappa_train = thr$kappa,
                      kappa_test = kappa_test)

# Scheme 4: Optimized Köppen
w <- optimize_koppen_weights(train)
thr <- optimize_thresholds(compute_ecodri_koppen(train, w), train$usdm_cat)
kappa_test <- quad_kappa(
  discretize(compute_ecodri_koppen(test, w), thr$thresholds),
  test$usdm_cat)
results[[4]] <- list(scheme = "4_optimized_koppen", kappa_train = thr$kappa,
                      kappa_test = kappa_test)

summary_tbl <- tibble(
  scheme = map_chr(results, "scheme"),
  kappa_train = map_dbl(results, "kappa_train"),
  kappa_test = map_dbl(results, "kappa_test"),
  train_test_gap = kappa_train - kappa_test
)

cat("\n=== Ablation summary (Section 4.4) ===\n")
print(summary_tbl)
write_csv(summary_tbl, file.path(OUTPUT_DIR, "ablation_summary.csv"))

cat("\nOutputs written to:", OUTPUT_DIR, "\n")
