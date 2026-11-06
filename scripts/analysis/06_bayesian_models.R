#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(brms)
  library(posterior)
})

options(mc.cores = min(4, parallel::detectCores()))

BASE <- "."

INFILE <- file.path(
  BASE,
  "data",
  "processed",
  "manuscript1_analysis_input.tsv"
)

OUT <- file.path(
  BASE,
  "outputs",
  "bayesian"
)

FITDIR <- file.path(OUT, "fits")
TABDIR <- file.path(OUT, "tables")
DIAGDIR <- file.path(OUT, "diagnostics")

dir.create(FITDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(DIAGDIR, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. DATA
# ============================================================

d <- read_tsv(INFILE, show_col_types=FALSE)

stopifnot(nrow(d) == 3227)
stopifnot(n_distinct(d$sample_id) == 3227)

d <- d %>%
  mutate(
    lineage = factor(
      lineage_model,
      levels=c("L4","L1","L2","L3","L5","L6")
    ),

    year = suppressWarnings(
      as.numeric(collection_year_clean)
    ),

    subregion = factor(
      subregion,
      levels=c(
        "East Africa",
        "Central Africa",
        "North Africa",
        "Southern Africa",
        "West Africa"
      )
    ),

    study = factor(study_accession),

    any_resistance = if_else(
      resistance_category_clean ==
        "no_resistance_predicted",
      0L, 1L
    ),

    mdr_or_higher = if_else(
      resistance_category_clean %in%
        c("MDR","pre_XDR","XDR"),
      1L, 0L
    ),

    advanced_resistance = if_else(
      resistance_category_clean %in%
        c("pre_XDR","XDR"),
      1L, 0L
    )
  )

# Full-cohort sensitivity dataset
dfull <- d %>%
  filter(
    !is.na(lineage),
    !is.na(study)
  ) %>%
  droplevels()

stopifnot(nrow(dfull) == 3227)

# Complete-case adjusted dataset
mdata <- d %>%
  filter(
    !is.na(lineage),
    !is.na(year),
    !is.na(study),
    !is.na(subregion)
  ) %>%
  droplevels()

stopifnot(nrow(mdata) == 3158)

mean_year <- mean(mdata$year)

mdata <- mdata %>%
  mutate(
    year5 = (year - mean_year) / 5
  )

cat("\n============================================\n")
cat("FINAL BAYESIAN OBJECTIVE 1 RERUN\n")
cat("============================================\n")
cat("Full cohort:", nrow(dfull), "\n")
cat("Adjusted cohort:", nrow(mdata), "\n")
cat("Excluded because of missing year:",
    nrow(dfull)-nrow(mdata), "\n")
cat("Studies:", n_distinct(mdata$study), "\n")
cat("Mean year:", mean_year, "\n")
cat("Any resistance events:",
    sum(mdata$any_resistance), "\n")
cat("MDR-or-higher events:",
    sum(mdata$mdr_or_higher), "\n")
cat("Advanced events:",
    sum(mdata$advanced_resistance), "\n")

# ============================================================
# 2. PRIORS
# ============================================================

pri_binary <- c(
  prior(normal(0,1), class="b"),
  prior(student_t(3,0,2.5), class="Intercept"),
  prior(student_t(3,0,2.5), class="sd")
)

# Advanced resistance:
# tighter regularisation because of sparse / zero-event strata.
pri_advanced <- c(
  prior(normal(0,1), class="b"),
  prior(normal(0,0.75), class="b", coef="lineageL1"),
  prior(normal(0,0.75), class="b", coef="lineageL2"),
  prior(normal(0,0.75), class="b", coef="lineageL3"),
  prior(normal(0,0.75), class="b", coef="lineageL5"),
  prior(normal(0,0.75), class="b", coef="lineageL6"),
  prior(student_t(3,0,2.5), class="Intercept")
)

CTRL_MAIN <- list(
  adapt_delta=0.99,
  max_treedepth=15
)

# ============================================================
# 3. FIT MODELS
# ============================================================

cat("\nFitting B1 primary...\n")

B1_primary <- brm(
  any_resistance ~ lineage + year5 + (1|study),
  data=mdata,
  family=bernoulli(link="logit"),
  backend="cmdstanr",
  prior=pri_binary,
  chains=4,
  iter=4000,
  warmup=2000,
  seed=20261007,
  control=CTRL_MAIN,
  refresh=200
)

saveRDS(
  B1_primary,
  file.path(FITDIR,"B1_primary_any_resistance.rds")
)

cat("\nFitting B1 full-cohort sensitivity...\n")

B1_full <- brm(
  any_resistance ~ lineage + (1|study),
  data=dfull,
  family=bernoulli(link="logit"),
  backend="cmdstanr",
  prior=pri_binary,
  chains=4,
  iter=4000,
  warmup=2000,
  seed=20261008,
  control=CTRL_MAIN,
  refresh=200
)

saveRDS(
  B1_full,
  file.path(FITDIR,"B1_full3227_any_resistance.rds")
)

cat("\nFitting B2 primary...\n")

B2_primary <- brm(
  mdr_or_higher ~ lineage + year5 + (1|study),
  data=mdata,
  family=bernoulli(link="logit"),
  backend="cmdstanr",
  prior=pri_binary,
  chains=4,
  iter=4000,
  warmup=2000,
  seed=20261009,
  control=CTRL_MAIN,
  refresh=200
)

saveRDS(
  B2_primary,
  file.path(FITDIR,"B2_primary_MDR_or_higher.rds")
)

cat("\nFitting B2 full-cohort sensitivity...\n")

B2_full <- brm(
  mdr_or_higher ~ lineage + (1|study),
  data=dfull,
  family=bernoulli(link="logit"),
  backend="cmdstanr",
  prior=pri_binary,
  chains=4,
  iter=4000,
  warmup=2000,
  seed=20261010,
  control=CTRL_MAIN,
  refresh=200
)

saveRDS(
  B2_full,
  file.path(FITDIR,"B2_full3227_MDR_or_higher.rds")
)

cat("\nFitting B3 primary...\n")

B3_primary <- brm(
  advanced_resistance ~ lineage + year5,
  data=mdata,
  family=bernoulli(link="logit"),
  backend="cmdstanr",
  prior=pri_advanced,
  chains=4,
  iter=6000,
  warmup=3000,
  seed=20261011,
  control=list(
    adapt_delta=0.995,
    max_treedepth=15
  ),
  refresh=200
)

saveRDS(
  B3_primary,
  file.path(FITDIR,"B3_primary_advanced_resistance.rds")
)

cat("\nFitting B3 geographic sensitivity...\n")

B3_geo <- brm(
  advanced_resistance ~ lineage + subregion + year5,
  data=mdata,
  family=bernoulli(link="logit"),
  backend="cmdstanr",
  prior=pri_advanced,
  chains=4,
  iter=6000,
  warmup=3000,
  seed=20261012,
  control=list(
    adapt_delta=0.995,
    max_treedepth=15
  ),
  refresh=200
)

saveRDS(
  B3_geo,
  file.path(FITDIR,"B3_geographic_sensitivity.rds")
)

# ============================================================
# 4. POSTERIOR LINEAGE OR EXTRACTION
# ============================================================

extract_lineage_OR <- function(
  fit,
  model,
  outcome,
  N,
  reference="L4"
) {

  draws <- as_draws_df(fit)

  terms <- paste0(
    "b_lineage",
    c("L1","L2","L3","L5","L6")
  )

  result <- map_dfr(
    terms,
    function(term) {

      if (!term %in% names(draws)) {
        return(NULL)
      }

      x <- draws[[term]]
      or <- exp(x)

      tibble(
        model=model,
        outcome=outcome,
        N=N,
        reference=reference,
        lineage=sub("^b_lineage","",term),

        OR_median=median(or),
        OR_mean=mean(or),

        CrI_2.5=quantile(or,0.025),
        CrI_97.5=quantile(or,0.975),

        Pr_OR_gt_1=mean(or > 1),
        Pr_OR_lt_1=mean(or < 1),

        Pr_beta_gt_0=mean(x > 0),
        Pr_beta_lt_0=mean(x < 0)
      )
    }
  )

  result
}

POST_OR <- bind_rows(

  extract_lineage_OR(
    B1_primary,
    "B1 primary: lineage + year + study RE",
    "Any genomic resistance",
    nrow(mdata)
  ),

  extract_lineage_OR(
    B1_full,
    "B1 sensitivity: full 3227 + study RE",
    "Any genomic resistance",
    nrow(dfull)
  ),

  extract_lineage_OR(
    B2_primary,
    "B2 primary: lineage + year + study RE",
    "MDR-or-higher",
    nrow(mdata)
  ),

  extract_lineage_OR(
    B2_full,
    "B2 sensitivity: full 3227 + study RE",
    "MDR-or-higher",
    nrow(dfull)
  ),

  extract_lineage_OR(
    B3_primary,
    "B3 primary: lineage + year",
    "Pre-XDR/XDR",
    nrow(mdata)
  ),

  extract_lineage_OR(
    B3_geo,
    "B3 sensitivity: lineage + subregion + year",
    "Pre-XDR/XDR",
    nrow(mdata)
  )
)

write_tsv(
  POST_OR,
  file.path(
    TABDIR,
    "Table_BAYES_FINAL_01_posterior_lineage_ORs.tsv"
  )
)

# ============================================================
# 5. MCMC DIAGNOSTICS
# ============================================================

diag_one <- function(fit, model) {

  draws <- posterior::as_draws_df(fit)

  summ <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )

  fixed <- summ %>%
    filter(
      grepl("^b_", variable) |
      grepl("^sd_", variable)
    )

  np <- suppressWarnings(
    tryCatch(
      brms::nuts_params(fit),
      error=function(e) NULL
    )
  )

  div <- NA_integer_
  max_td <- NA_integer_

  if (!is.null(np)) {

    div <- sum(
      np$Parameter == "divergent__" &
      np$Value == 1
    )

    td <- np$Value[
      np$Parameter == "treedepth__"
    ]

    if (length(td) > 0) {
      max_td <- max(td)
    }
  }

  tibble(
    model=model,
    max_Rhat=max(fixed$rhat, na.rm=TRUE),
    min_Bulk_ESS=min(fixed$ess_bulk, na.rm=TRUE),
    min_Tail_ESS=min(fixed$ess_tail, na.rm=TRUE),
    divergent_transitions=div,
    maximum_observed_treedepth=max_td
  )
}

DIAG <- bind_rows(
  diag_one(B1_primary,"B1_primary"),
  diag_one(B1_full,"B1_full3227"),
  diag_one(B2_primary,"B2_primary"),
  diag_one(B2_full,"B2_full3227"),
  diag_one(B3_primary,"B3_primary"),
  diag_one(B3_geo,"B3_geographic_sensitivity")
)

write_tsv(
  DIAG,
  file.path(
    TABDIR,
    "Table_BAYES_FINAL_02_MCMC_diagnostics.tsv"
  )
)

# ============================================================
# 6. STUDY RANDOM-EFFECT SD
# ============================================================

extract_sd <- function(fit,model) {

  d <- as_draws_df(fit)

  term <- grep(
    "^sd_study__Intercept$",
    names(d),
    value=TRUE
  )

  if (length(term) == 0) {
    return(
      tibble(
        model=model,
        SD_median=NA_real_,
        CrI_2.5=NA_real_,
        CrI_97.5=NA_real_
      )
    )
  }

  x <- d[[term[1]]]

  tibble(
    model=model,
    SD_median=median(x),
    CrI_2.5=quantile(x,0.025),
    CrI_97.5=quantile(x,0.975)
  )
}

RE_SD <- bind_rows(
  extract_sd(B1_primary,"B1_primary"),
  extract_sd(B1_full,"B1_full3227"),
  extract_sd(B2_primary,"B2_primary"),
  extract_sd(B2_full,"B2_full3227")
)

write_tsv(
  RE_SD,
  file.path(
    TABDIR,
    "Table_BAYES_FINAL_03_study_random_effect_SD.tsv"
  )
)

# ============================================================
# 7. COMPLETE MODEL SUMMARIES
# ============================================================

capture.output(
  summary(B1_primary),
  file=file.path(DIAGDIR,"B1_primary_summary.txt")
)

capture.output(
  summary(B1_full),
  file=file.path(DIAGDIR,"B1_full3227_summary.txt")
)

capture.output(
  summary(B2_primary),
  file=file.path(DIAGDIR,"B2_primary_summary.txt")
)

capture.output(
  summary(B2_full),
  file=file.path(DIAGDIR,"B2_full3227_summary.txt")
)

capture.output(
  summary(B3_primary),
  file=file.path(DIAGDIR,"B3_primary_summary.txt")
)

capture.output(
  summary(B3_geo),
  file=file.path(DIAGDIR,"B3_geographic_sensitivity_summary.txt")
)

# ============================================================
# 8. MODEL-COHORT AUDIT
# ============================================================

AUDIT <- tibble(
  dataset=c(
    "Full descriptive cohort",
    "Primary adjusted models"
  ),
  N=c(
    nrow(dfull),
    nrow(mdata)
  ),
  any_resistance_events=c(
    sum(dfull$any_resistance),
    sum(mdata$any_resistance)
  ),
  MDR_or_higher_events=c(
    sum(dfull$mdr_or_higher),
    sum(mdata$mdr_or_higher)
  ),
  advanced_events=c(
    sum(dfull$advanced_resistance),
    sum(mdata$advanced_resistance)
  ),
  studies=c(
    n_distinct(dfull$study),
    n_distinct(mdata$study)
  )
)

write_tsv(
  AUDIT,
  file.path(
    TABDIR,
    "Table_BAYES_FINAL_00_model_cohort_audit.tsv"
  )
)

cat("\n============================================\n")
cat("POSTERIOR LINEAGE ORs\n")
cat("============================================\n")
print(POST_OR, n=Inf)

cat("\n============================================\n")
cat("MCMC DIAGNOSTICS\n")
cat("============================================\n")
print(DIAG, n=Inf)

cat("\n============================================\n")
cat("STUDY RANDOM-EFFECT SDs\n")
cat("============================================\n")
print(RE_SD, n=Inf)

cat("\nBayesian rerun complete.\n")
