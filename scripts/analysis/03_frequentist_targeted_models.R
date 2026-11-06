suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lme4)
  library(brglm2)
})

BASE <- "."

INFILE <- file.path(
  BASE,
  "data",
  "processed",
  "manuscript1_analysis_input.tsv"
)

OUTDIR <- file.path(
  BASE,
  "outputs",
  "tables",
  "frequentist"
)

MODELDIR <- file.path(
  BASE,
  "outputs",
  "models",
  "frequentist"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODELDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD CORRECTED COHORT
# ============================================================

d <- read_tsv(INFILE, show_col_types = FALSE)

stopifnot(nrow(d) == 3227)
stopifnot(n_distinct(d$sample_id) == 3227)

cat("============================================================\n")
cat("OBJECTIVE 1 FINAL TARGETED RERUN\n")
cat("============================================================\n")
cat("Full corrected cohort:", nrow(d), "\n")

# ============================================================
# 2. STANDARDISE VARIABLES
# ============================================================

d <- d %>%
  mutate(
    lineage = factor(
      lineage_model,
      levels = c("L4", "L1", "L2", "L3", "L5", "L6")
    ),

    year = suppressWarnings(
      as.numeric(collection_year_clean)
    ),

    subregion = factor(
      subregion,
      levels = c(
        "East Africa",
        "Central Africa",
        "North Africa",
        "Southern Africa",
        "West Africa"
      )
    ),

    study_accession = factor(study_accession),

    any_resistance = if_else(
      resistance_category_clean ==
        "no_resistance_predicted",
      0L, 1L
    ),

    mdr_or_higher = if_else(
      resistance_category_clean %in%
        c("MDR", "pre_XDR", "XDR"),
      1L, 0L
    ),

    advanced_resistance = if_else(
      resistance_category_clean %in%
        c("pre_XDR", "XDR"),
      1L, 0L
    )
  )

# ============================================================
# 3. DESCRIPTIVE FULL COHORT — N = 3227
# ============================================================

full_lineage <- d %>%
  group_by(lineage) %>%
  summarise(
    n = n(),

    susceptible =
      sum(any_resistance == 0),

    any_resistant =
      sum(any_resistance),

    any_resistant_pct =
      100 * mean(any_resistance),

    MDR_or_higher =
      sum(mdr_or_higher),

    MDR_or_higher_pct =
      100 * mean(mdr_or_higher),

    preXDR_XDR =
      sum(advanced_resistance),

    preXDR_XDR_pct =
      100 * mean(advanced_resistance),

    .groups = "drop"
  )

write_tsv(
  full_lineage,
  file.path(
    OUTDIR,
    "Table_01_full3227_outcomes_by_lineage.tsv"
  )
)

country_lineage <- d %>%
  count(
    country_final_curated,
    lineage,
    name = "n"
  )

write_tsv(
  country_lineage,
  file.path(
    OUTDIR,
    "Table_02_country_by_lineage_3227.tsv"
  )
)

subregion_lineage <- d %>%
  count(
    subregion,
    lineage,
    name = "n"
  )

write_tsv(
  subregion_lineage,
  file.path(
    OUTDIR,
    "Table_03_subregion_by_lineage_3227.tsv"
  )
)

resistance_category_table <- d %>%
  count(
    resistance_category_clean,
    name = "n"
  ) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  arrange(desc(n))

write_tsv(
  resistance_category_table,
  file.path(
    OUTDIR,
    "Table_04_resistance_categories_3227.tsv"
  )
)

# ============================================================
# 4. FINAL COMPLETE-CASE MODEL COHORT
# ============================================================

mdata <- d %>%
  filter(
    !is.na(lineage),
    !is.na(year),
    !is.na(study_accession),
    !is.na(subregion)
  ) %>%
  droplevels()

if (nrow(mdata) != 3158) {
  stop(
    paste(
      "Expected final adjusted-model N=3158; found",
      nrow(mdata)
    )
  )
}

mean_year <- mean(mdata$year)

mdata <- mdata %>%
  mutate(
    year5 = (year - mean_year) / 5
  )

cat("\nFinal adjusted-model N:", nrow(mdata), "\n")
cat("Contributing studies:",
    n_distinct(mdata$study_accession), "\n")
cat("Mean year:", round(mean_year, 3), "\n")
cat("Any resistance:",
    sum(mdata$any_resistance), "\n")
cat("MDR-or-higher:",
    sum(mdata$mdr_or_higher), "\n")
cat("Pre-XDR/XDR:",
    sum(mdata$advanced_resistance), "\n")

# ============================================================
# 5. MODEL-FITTING HELPER WITH CONTROLLED OPTIMIZER FALLBACK
# ============================================================

fit_glmer_safe <- function(formula, data, label) {

  cat("\n------------------------------------------------------------\n")
  cat("FITTING:", label, "\n")
  cat("------------------------------------------------------------\n")

  optimizers <- c(
    "bobyqa",
    "nloptwrap"
  )

  fits <- list()

  for (opt in optimizers) {

    cat("Trying optimizer:", opt, "\n")

    fit <- tryCatch(
      {
        if (opt == "bobyqa") {
          glmer(
            formula,
            data = data,
            family = binomial(),
            control = glmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 5e5)
            )
          )
        } else {
          glmer(
            formula,
            data = data,
            family = binomial(),
            control = glmerControl(
              optimizer = "nloptwrap",
              optCtrl = list(maxeval = 5e5)
            )
          )
        }
      },
      error = function(e) {
        cat(
          "FAILED:",
          conditionMessage(e),
          "\n"
        )
        NULL
      }
    )

    if (!is.null(fit)) {

      grad <- NA_real_

      derivs <- fit@optinfo$derivs

      if (
        !is.null(derivs) &&
        !is.null(derivs$gradient)
      ) {
        grad <- max(
          abs(derivs$gradient)
        )
      }

      fits[[opt]] <- list(
        fit = fit,
        singular = isSingular(
          fit,
          tol = 1e-4
        ),
        gradient = grad,
        AIC = AIC(fit)
      )

      cat("SUCCESS\n")
      cat(
        "Singular:",
        isSingular(fit),
        "\n"
      )
      cat(
        "Max gradient:",
        grad,
        "\n"
      )
      cat(
        "AIC:",
        AIC(fit),
        "\n"
      )
    }
  }

  if (length(fits) == 0) {
    stop(
      paste(
        "No optimizer successfully fitted",
        label
      )
    )
  }

  # Prefer bobyqa when successful and acceptable.
  if (
    "bobyqa" %in% names(fits) &&
    !fits[["bobyqa"]]$singular &&
    (
      is.na(fits[["bobyqa"]]$gradient) ||
      fits[["bobyqa"]]$gradient < 0.01
    )
  ) {
    chosen <- "bobyqa"
  } else if (
    "nloptwrap" %in% names(fits)
  ) {
    chosen <- "nloptwrap"
  } else {
    chosen <- names(fits)[1]
  }

  cat(
    "Selected optimizer:",
    chosen,
    "\n"
  )

  attr(
    fits[[chosen]]$fit,
    "selected_optimizer"
  ) <- chosen

  fits[[chosen]]$fit
}

# ============================================================
# 6. CRUDE MODELS
# ============================================================

crude_any <- glm(
  any_resistance ~ lineage,
  data = mdata,
  family = binomial()
)

crude_mdr <- glm(
  mdr_or_higher ~ lineage,
  data = mdata,
  family = binomial()
)

# ============================================================
# 7. PRIMARY MODELS
#
# Outcome ~ lineage + year + study random intercept
#
# Study RE is primary because 101/104 studies occur in only
# one subregion; study and subregion are therefore largely
# nested/confounded.
# ============================================================

primary_any <- fit_glmer_safe(
  any_resistance ~
    lineage + year5 +
    (1 | study_accession),
  mdata,
  "PRIMARY ANY RESISTANCE"
)

primary_mdr <- fit_glmer_safe(
  mdr_or_higher ~
    lineage + year5 +
    (1 | study_accession),
  mdata,
  "PRIMARY MDR-OR-HIGHER"
)

# ============================================================
# 8. GEOGRAPHIC STRUCTURE
#
# A mixed-effects model containing both subregion and a
# study-accession random intercept was not retained because
# subregion was nearly nested within study:
# 101/104 contributing studies occurred in only one subregion.
#
# Geographic variation is therefore reported descriptively
# from the complete 3,227-isolate country/subregion cohort.
#
# ============================================================

# ============================================================
# 9. ADVANCED RESISTANCE
#
# Bias-reduced logistic regression retained because L5/L6
# have zero pre-XDR/XDR events.
#
# Primary: lineage + year
# Sensitivity: lineage + subregion + year
# ============================================================

advanced_primary <- glm(
  advanced_resistance ~
    lineage + year5,
  data = mdata,
  family = binomial(),
  method = brglmFit,
  type = "AS_mean"
)

advanced_sensitivity <- glm(
  advanced_resistance ~
    lineage + subregion + year5,
  data = mdata,
  family = binomial(),
  method = brglmFit,
  type = "AS_mean"
)

# ============================================================
# 10. COEFFICIENT EXTRACTION
# ============================================================

extract_or <- function(
  model,
  outcome,
  model_name
) {

  cf <- summary(model)$coefficients

  zcol <- grep(
    "^Pr\\(",
    colnames(cf),
    value = TRUE
  )[1]

  out <- data.frame(
    term = rownames(cf),
    beta = cf[, "Estimate"],
    SE = cf[, "Std. Error"],
    p_value = cf[, zcol],
    stringsAsFactors = FALSE
  )

  out %>%
    filter(grepl("^lineage", term)) %>%
    mutate(
      lineage =
        sub("^lineage", "", term),

      OR =
        exp(beta),

      conf_low =
        exp(beta - 1.96 * SE),

      conf_high =
        exp(beta + 1.96 * SE),

      outcome =
        outcome,

      model =
        model_name
    ) %>%
    select(
      outcome,
      model,
      lineage,
      beta,
      SE,
      OR,
      conf_low,
      conf_high,
      p_value
    )
}

OR_TABLE <- bind_rows(

  extract_or(
    crude_any,
    "Any genomic resistance",
    "Crude logistic"
  ),

  extract_or(
    primary_any,
    "Any genomic resistance",
    "Primary: lineage + year + study RE"
  ),

  extract_or(
    crude_mdr,
    "MDR-or-higher",
    "Crude logistic"
  ),

  extract_or(
    primary_mdr,
    "MDR-or-higher",
    "Primary: lineage + year + study RE"
  ),

  extract_or(
    advanced_primary,
    "Pre-XDR/XDR",
    "Primary bias-reduced: lineage + year"
  ),

  extract_or(
    advanced_sensitivity,
    "Pre-XDR/XDR",
    "Sensitivity bias-reduced: lineage + subregion + year"
  )
)

write_tsv(
  OR_TABLE,
  file.path(
    OUTDIR,
    "Table_05_FINAL_lineage_ORs_all_models.tsv"
  )
)

# ============================================================
# 11. MODEL DIAGNOSTICS
# ============================================================

mixed_diag <- function(
  model,
  outcome,
  specification
) {

  grad <- NA_real_

  if (
    !is.null(model@optinfo$derivs) &&
    !is.null(
      model@optinfo$derivs$gradient
    )
  ) {
    grad <- max(
      abs(
        model@optinfo$derivs$gradient
      )
    )
  }

  tibble(
    outcome = outcome,
    specification = specification,

    optimizer =
      attr(
        model,
        "selected_optimizer"
      ),

    N =
      nobs(model),

    studies =
      n_distinct(mdata$study_accession),

    singular =
      isSingular(
        model,
        tol = 1e-4
      ),

    max_gradient =
      grad,

    study_RE_SD =
      attr(
        VarCorr(model)$study_accession,
        "stddev"
      )[1],

    AIC =
      AIC(model)
  )
}

DIAG <- bind_rows(

  mixed_diag(
    primary_any,
    "Any genomic resistance",
    "Primary: lineage + year + study RE"
  ),

  mixed_diag(
    primary_mdr,
    "MDR-or-higher",
    "Primary: lineage + year + study RE"
  )
)

write_tsv(
  DIAG,
  file.path(
    OUTDIR,
    "Table_06_FINAL_mixed_model_diagnostics.tsv"
  )
)

# ============================================================
# 12. STUDY/SUBREGION IDENTIFIABILITY TABLE
# ============================================================

study_geo <- mdata %>%
  distinct(
    study_accession,
    subregion
  ) %>%
  count(
    study_accession,
    name = "n_subregions"
  )

identifiability <- tibble(
  complete_case_N = nrow(mdata),
  contributing_studies =
    n_distinct(
      mdata$study_accession
    ),
  studies_one_subregion =
    sum(
      study_geo$n_subregions == 1
    ),
  studies_multiple_subregions =
    sum(
      study_geo$n_subregions > 1
    )
)

write_tsv(
  identifiability,
  file.path(
    OUTDIR,
    "Table_07_study_subregion_identifiability.tsv"
  )
)

# ============================================================
# 13. SAVE MODELS
# ============================================================

saveRDS(
  primary_any,
  file.path(
    MODELDIR,
    "primary_any_resistance.rds"
  )
)

saveRDS(
  primary_mdr,
  file.path(
    MODELDIR,
    "primary_MDR_or_higher.rds"
  )
)

saveRDS(
  advanced_primary,
  file.path(
    MODELDIR,
    "advanced_primary_brglm.rds"
  )
)

saveRDS(
  advanced_sensitivity,
  file.path(
    MODELDIR,
    "advanced_subregion_sensitivity_brglm.rds"
  )
)

# ============================================================
# 14. FINAL AUDIT
# ============================================================

cat("\n============================================================\n")
cat("FINAL COHORT SUMMARY\n")
cat("============================================================\n")

cat(
  "Full descriptive cohort:",
  nrow(d),
  "\n"
)

cat(
  "Adjusted-model cohort:",
  nrow(mdata),
  "\n"
)

cat(
  "Excluded for missing year:",
  nrow(d) - nrow(mdata),
  "\n"
)

cat(
  "Contributing studies:",
  n_distinct(mdata$study_accession),
  "\n"
)

cat(
  "Studies restricted to one subregion:",
  sum(
    study_geo$n_subregions == 1
  ),
  "\n"
)

cat(
  "Studies spanning multiple subregions:",
  sum(
    study_geo$n_subregions > 1
  ),
  "\n"
)

cat("\n============================================================\n")
cat("PRIMARY LINEAGE RESULTS\n")
cat("============================================================\n")

print(
  OR_TABLE %>%
    filter(
      grepl("^Primary", model)
    )
)

cat("\n============================================================\n")
cat("MIXED MODEL DIAGNOSTICS\n")
cat("============================================================\n")

print(DIAG)

cat("\n============================================================\n")
cat("TARGETED OBJECTIVE 1 RERUN COMPLETE\n")
cat("============================================================\n")
