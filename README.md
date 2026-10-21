# hypoforgeR — Discover hypotheses, validate incremental value

**hypoforgeR** is an experimental, standalone R package for interpretable feature discovery, numerical safety checks, calendar and marketing-mix transformations, and **rolling-origin Ridge tournaments**. It does not require Python.

The question is not *“Which transformed variable correlates most with sales?”* It is **“Does an economically meaningful transformation improve a realistic baseline on later observations?”**

> **Scope:** Predictive, exploratory feature analysis. Not causal inference, an automatic marketing-mix model, a guarantee of future improvement, or proof that any particular mathematical response curve is correct.

## Install

R 4.1 or later:

```r
install.packages("remotes")
remotes::install_github("<new-owner>/HypoForge-R")
library(hypoforgeR)
```

For a local checkout, run `install.packages(".", repos = NULL, type = "source")` from the package directory. Pin a Git commit for reproducibility.

## Quick start: Is distribution nonlinear?

This example **runs without private data** and illustrates the workflow, not a benchmark claim.

```r
# Run from the repository root after installing hypoforgeR.
d <- read.csv("examples/hypothetical_weekly_sales.csv")
train <- d[seq_len(floor(nrow(d) * 0.8)), , drop = FALSE]

ranked <- hf_discover(train, target = "sales", top_k = 10L)
print(ranked)

candidates <- list(
  hf_spec("distribution__square", "square", "distribution"),
  hf_spec("distribution__cube", "cube", "distribution")
)
# Use a genuine date/time column for hf_tournament(time=...).
# The hypothetical CSV intentionally uses an abstract sequential week index.
```


`hf_discover()` ranks candidates by **absolute in-sample Pearson correlation**. It is for hypothesis generation only. `hf_tournament()` compares the *same* Ridge baseline with and without each candidate on later weeks. `mean_improvement` is **baseline RMSE minus candidate RMSE** (positive is better); `positive_fold_fraction` is the share of improved folds. Preprocessing is fitted within each training fold.

**Forecast availability warning:** In the quick start, price and distribution are labeled `known_future` for illustration. For genuine forecasts, do this **only** when their future values are available at the forecast origin. Otherwise forecast those drivers separately or use an appropriate scenario design.

## Hypothetical case study: nonlinear distribution and sparse signals

**Entirely hypothetical illustration.** The file [`examples/hypothetical_weekly_sales.csv`](examples/hypothetical_weekly_sales.csv) is deliberately constructed from round values. It contains no observations about any real product, organization, market or time period. The example is intended to demonstrate a validation workflow, **not** to claim a measured improvement or a causal effect.

### Hypothesis A — distribution and seasonal demand

In a hypothetical system, distribution can increase ahead of a seasonal demand period. The association between distribution and sales might change with seasonal activation. A linear relationship is only one candidate; squared, cubic or distribution-by-season terms are hypotheses, not established effects.

### Step 1 — Reserve an untouched holdout before exploration

```r
d <- read.csv("examples/hypothetical_weekly_sales.csv")
stopifnot(!anyDuplicated(d$week_index), all(diff(d$week_index) > 0))
cutoff <- floor(nrow(d) * 0.8)
train <- d[seq_len(cutoff), , drop = FALSE]
holdout <- d[seq.int(cutoff + 1L, nrow(d)), , drop = FALSE]
```

Do not shuffle time series. Do not inspect holdout outcomes when choosing features. A future driver can be marked `known_future` only if its value is genuinely available at the forecast origin; using retrospectively observed future values is a diagnostic, not automatically an ex-ante forecast.

### Step 2 — Specify the baseline and screen training-only hypotheses

```r
baseline <- c("distribution", "season_index", "promotion_signal")
ranked <- hf_discover(train[, c("sales", baseline)], target = "sales", top_k = 10L)
candidates <- list(
  hf_spec("distribution__square", "square", "distribution"),
  hf_spec("distribution__cube", "cube", "distribution")
)
print(ranked)
```

Screening correlations are not incremental validation. Keep a defensible baseline, check collinearity, numerical safety and future availability, and avoid selecting features on holdout outcomes.

### Step 3 — Validate with rolling origins

The R package provides `hf_tournament()` for time-ordered Ridge baseline comparisons. Its preprocessing is fitted separately inside each training fold. The provided CSV uses a sequential `week_index`; for the API's `time` contract, add a synthetic date index or use a separate genuinely dated dataset. Choose `min_train`, `horizon` and `step` based on available observations, not on any previous investigation. Compare candidate and baseline on identical folds, and inspect fold-level results. No numeric tournament outcomes are asserted here.

### Hypothesis B — sparse promotion extrapolation

The hypothetical `promotion_signal` is zero in almost all initial observations and takes a much larger value late in the sequence. These are **constructed input values**, not measured outcomes. A Ridge model trained on sparse low-valued inputs can produce unstable predictions when an unseen value is far beyond the training range. Ridge shrinkage penalizes coefficients; it does **not** cap extrapolated predictions. A numerical forecast-distortion claim would require actually fitting and evaluating a specified model.

### Step 4 — Diagnose without rewriting the test

Inspect training support, maximum observed training values, standardized distances, leverage, measurement units and structural changes. A diagnostic exclusion of an extreme observation is not a replacement for a pristine holdout. If a new interaction or saturation transform is proposed after seeing test errors, validate it on new data.

### Step 8 — A beginner-friendly checklist

1. Define the target, forecast horizon, and which future inputs will really be available.
2. Sort chronologically and reserve an untouched final period.
3. Include economically relevant original drivers in the baseline.
4. Generate interpretable candidates on training data only.
5. Compare baseline and baseline-plus-candidate on identical rolling-origin folds.
6. Inspect every fold, not just average improvement.
7. Diagnose sparse promotion indicators, unusual z-scores, measurement changes, collinearity and structural shifts.
8. Propose saturation or season interactions only as **new hypotheses**, then revalidate them.
9. Make your final claim only after evaluating a locked pipeline on a genuinely untouched period.

**Bottom line:** The useful outcome is not “the cube always wins.” It is a disciplined, auditable path from **hypothesis → time-aware test → stress test → qualified conclusion**.

## API overview

| Function | Purpose |
| --- | --- |
| `hf_spec`, `hf_evaluate` | Define and compute interpretable candidate features |
| `hf_generate`, `hf_discover` | Generate and rank exploratory hypotheses |
| `hf_safety`, `hf_select` | Numerical safety and redundancy screening |
| `hf_rolling_splits`, `hf_tournament` | Time-ordered splits and Ridge baseline comparisons |
| `hf_easter`, `hf_calendar_features` | Easter and calendar features |
| `hf_scan_recurring` | Exploratory recurring patterns in **training residuals** |
| `hf_geometric_adstock`, `hf_delayed_adstock`, `hf_weibull_adstock` | Deterministic carryover transformations |
| `hf_logistic_saturation`, `hf_michaelis_menten`, `hf_power_saturation`, `hf_tanh_saturation` | Deterministic response curves |
| `hf_reach_frequency` | Reach × saturated-frequency transform |

Adstock functions require one correctly ordered time series; never concatenate unrelated panels. Calendar features are daily-date aligned; aggregate to your **actual business week** when working with weekly data. Residual phase scans are exploratory, not significance tests.

## R and Python: related, not identical

The [Python project](https://github.com/sayac007/HypoForge-Feature-Discovery) supports its own discovery and screening pipeline. The **R** ranker uses absolute Pearson correlation and includes a native `hf_tournament()`; the Python ranker and API are different, and its current package does **not** provide an equivalent native tournament. Do not interpret results as algorithmically identical.

## License and contributions

Apache-2.0. Reproducible issues, test cases and independent validation are welcome. Run `R CMD check .` before a release.
