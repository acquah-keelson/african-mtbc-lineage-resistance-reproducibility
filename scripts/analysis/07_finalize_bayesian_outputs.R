#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(brms)
  library(rstan)
})

BASE <- file.path(
  "outputs",
  "bayesian"
)

FITDIR  <- file.path(BASE, "fits")
TABDIR  <- file.path(BASE, "tables")
DIAGDIR <- file.path(BASE, "diagnostics")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(DIAGDIR, recursive=TRUE, showWarnings=FALSE)

fit_files <- c(
  B1_primary = "B1_primary_any_resistance.rds",
  B1_full    = "B1_full3227_any_resistance.rds",
  B2_primary = "B2_primary_MDR_or_higher.rds",
  B2_full    = "B2_full3227_MDR_or_higher.rds",
  B3_primary = "B3_primary_advanced_resistance.rds",
  B3_geo     = "B3_geographic_sensitivity.rds"
)

fits <- lapply(
  file.path(FITDIR, fit_files),
  readRDS
)

names(fits) <- names(fit_files)

cat("Loaded all six Bayesian fits.\n")
cat("Fit names:", paste(names(fits), collapse=", "), "\n")

# ============================================================
# POSTERIOR LINEAGE ORs FROM brms::fixef()
# ============================================================

extract_lineage_OR <- function(fit, model, outcome, N) {

  x <- brms::fixef(
    fit,
    summary = FALSE
  )

  wanted <- c(
    "lineageL1",
    "lineageL2",
    "lineageL3",
    "lineageL5",
    "lineageL6"
  )

  available <- intersect(
    wanted,
    colnames(x)
  )

  if (length(available) == 0) {
    cat("\nAvailable fixed effects:\n")
    print(colnames(x))
    stop(paste("No lineage coefficients found for", model))
  }

  map_dfr(
    available,
    function(term) {

      beta <- x[, term]
      OR <- exp(beta)

      tibble(
        model = model,
        outcome = outcome,
        N = N,
        reference = "L4",
        lineage = sub("^lineage", "", term),

        beta_median = median(beta),

        OR_median = median(OR),

        OR_mean = mean(OR),

        CrI_2.5 =
          unname(
            quantile(
              OR,
              0.025
            )
          ),

        CrI_97.5 =
          unname(
            quantile(
              OR,
              0.975
            )
          ),

        Pr_OR_gt_1 =
          mean(
            OR > 1
          ),

        Pr_OR_lt_1 =
          mean(
            OR < 1
          )
      )
    }
  )
}

POST_OR <- bind_rows(

  extract_lineage_OR(
    fits$B1_primary,
    "B1 primary: lineage + year + study RE",
    "Any genomic resistance",
    3158
  ),

  extract_lineage_OR(
    fits$B1_full,
    "B1 sensitivity: full 3227 + study RE",
    "Any genomic resistance",
    3227
  ),

  extract_lineage_OR(
    fits$B2_primary,
    "B2 primary: lineage + year + study RE",
    "MDR-or-higher",
    3158
  ),

  extract_lineage_OR(
    fits$B2_full,
    "B2 sensitivity: full 3227 + study RE",
    "MDR-or-higher",
    3227
  ),

  extract_lineage_OR(
    fits$B3_primary,
    "B3 primary: lineage + year",
    "Pre-XDR/XDR",
    3158
  ),

  extract_lineage_OR(
    fits$B3_geo,
    "B3 sensitivity: lineage + subregion + year",
    "Pre-XDR/XDR",
    3158
  )
)

write_tsv(
  POST_OR,
  file.path(
    TABDIR,
    "Table_BAYES_FINAL_01_posterior_lineage_ORs.tsv"
  )
)

cat("Posterior lineage OR table written.\n")

# ============================================================
# MCMC DIAGNOSTICS FROM EMBEDDED stanfit
# ============================================================

diag_one <- function(fit, model) {

  sf <- fit$fit

  if (!inherits(sf, "stanfit")) {
    stop(paste(model, "does not contain a stanfit object"))
  }

  s <- rstan::summary(sf)$summary

  # Use all sampled parameters except lp__
  keep <- rownames(s) != "lp__"

  ss <- s[
    keep,
    ,
    drop = FALSE
  ]

  max_rhat <- max(
    ss[, "Rhat"],
    na.rm = TRUE
  )

  min_neff <- min(
    ss[, "n_eff"],
    na.rm = TRUE
  )

  sampler <- rstan::get_sampler_params(
    sf,
    inc_warmup = FALSE
  )

  divergent <- sum(
    vapply(
      sampler,
      function(z)
        sum(
          z[, "divergent__"],
          na.rm = TRUE
        ),
      numeric(1)
    )
  )

  max_td <- max(
    vapply(
      sampler,
      function(z)
        max(
          z[, "treedepth__"],
          na.rm = TRUE
        ),
      numeric(1)
    )
  )

  tibble(
    model = model,
    max_Rhat = max_rhat,
    min_n_eff = min_neff,
    divergent_transitions = divergent,
    max_observed_treedepth = max_td
  )
}

DIAG <- bind_rows(

  diag_one(fits$B1_primary, "B1_primary"),
  diag_one(fits$B1_full,    "B1_full3227"),
  diag_one(fits$B2_primary, "B2_primary"),
  diag_one(fits$B2_full,    "B2_full3227"),
  diag_one(fits$B3_primary, "B3_primary"),
  diag_one(fits$B3_geo,     "B3_geographic_sensitivity")
)

write_tsv(
  DIAG,
  file.path(
    TABDIR,
    "Table_BAYES_FINAL_02_MCMC_diagnostics.tsv"
  )
)

cat("MCMC diagnostics written.\n")

# ============================================================
# STUDY RANDOM-EFFECT SD USING VarCorr()
# ============================================================

extract_study_sd <- function(fit, model) {

  vc <- brms::VarCorr(
    fit,
    summary = FALSE
  )

  if (!"study" %in% names(vc)) {
    return(
      tibble(
        model = model,
        SD_median = NA_real_,
        CrI_2.5 = NA_real_,
        CrI_97.5 = NA_real_
      )
    )
  }

  sd_obj <- vc$study$sd

  # Convert whatever shape brms returns to numeric draws
  vals <- as.numeric(sd_obj)

  vals <- vals[
    is.finite(vals)
  ]

  tibble(
    model = model,

    SD_median =
      median(vals),

    CrI_2.5 =
      unname(
        quantile(
          vals,
          0.025
        )
      ),

    CrI_97.5 =
      unname(
        quantile(
          vals,
          0.975
        )
      )
  )
}

RE_SD <- bind_rows(

  extract_study_sd(
    fits$B1_primary,
    "B1_primary"
  ),

  extract_study_sd(
    fits$B1_full,
    "B1_full3227"
  ),

  extract_study_sd(
    fits$B2_primary,
    "B2_primary"
  ),

  extract_study_sd(
    fits$B2_full,
    "B2_full3227"
  )
)

write_tsv(
  RE_SD,
  file.path(
    TABDIR,
    "Table_BAYES_FINAL_03_study_random_effect_SD.tsv"
  )
)

cat("Study random-effect SD table written.\n")

# ============================================================
# SAVE BRMS SUMMARIES
# ============================================================

for (nm in names(fits)) {

  capture.output(
    summary(
      fits[[nm]]
    ),
    file =
      file.path(
        DIAGDIR,
        paste0(
          nm,
          "_summary.txt"
        )
      )
  )
}

# ============================================================
# SCREEN OUTPUT
# ============================================================

cat("\n====================================================\n")
cat("POSTERIOR LINEAGE ODDS RATIOS\n")
cat("====================================================\n")

print(
  POST_OR,
  n = Inf,
  width = Inf
)

cat("\n====================================================\n")
cat("MCMC DIAGNOSTICS\n")
cat("====================================================\n")

print(
  DIAG,
  n = Inf,
  width = Inf
)

cat("\n====================================================\n")
cat("STUDY RANDOM-EFFECT SD\n")
cat("====================================================\n")

print(
  RE_SD,
  n = Inf,
  width = Inf
)

cat("\nFINAL BAYESIAN CLOSURE COMPLETE\n")

