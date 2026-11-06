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
  "frequentist_core"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

d <- read_tsv(INFILE, show_col_types = FALSE)

cat("Input rows:", nrow(d), "\n")

# ------------------------------------------------------------
# Resolve columns
# ------------------------------------------------------------

lineage_col <- intersect(
  c("lineage_model", "lineage_final", "tbprof_main_lineage"),
  names(d)
)[1]

year_candidates <- intersect(
  c(
    "collection_year_clean",
    "collection_year",
    "year",
    "sampling_year",
    "collection_year_num",
    "collection_year_numeric"
  ),
  names(d)
)

if (length(year_candidates) == 0) {
  stop("No candidate year column found")
}

year_nonmissing <- sapply(
  year_candidates,
  function(x) {
    sum(!is.na(suppressWarnings(as.numeric(d[[x]]))))
  }
)

year_col <- year_candidates[
  which.max(year_nonmissing)
]

cat(
  "Selected year column:",
  year_col,
  "with",
  max(year_nonmissing),
  "numeric values\n"
)

subregion_col <- intersect(
  c("subregion", "african_subregion", "subregion_final"),
  names(d)
)[1]

stopifnot(
  !is.na(lineage_col),
  !is.na(year_col),
  !is.na(subregion_col)
)

# ------------------------------------------------------------
# Standardise lineage
# ------------------------------------------------------------

d <- d %>%
  mutate(
    lineage = case_when(
      grepl("1$", .data[[lineage_col]]) ~ "L1",
      grepl("2$", .data[[lineage_col]]) ~ "L2",
      grepl("3$", .data[[lineage_col]]) ~ "L3",
      grepl("4$", .data[[lineage_col]]) ~ "L4",
      grepl("5$", .data[[lineage_col]]) ~ "L5",
      grepl("6$", .data[[lineage_col]]) ~ "L6",
      TRUE ~ as.character(.data[[lineage_col]])
    ),
    year = suppressWarnings(
      as.numeric(.data[[year_col]])
    ),
    subregion = as.factor(.data[[subregion_col]]),
    study_accession = as.factor(study_accession)
  )

d$lineage <- factor(
  d$lineage,
  levels = c("L4", "L1", "L2", "L3", "L5", "L6")
)

# ------------------------------------------------------------
# Derive outcomes from the frozen resistance classification
# ------------------------------------------------------------

stopifnot("resistance_category_clean" %in% names(d))

d <- d %>%
  mutate(
    any_resistance = if_else(
      resistance_category_clean == "no_resistance_predicted",
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

# ------------------------------------------------------------
# Final complete-case cohort
# ------------------------------------------------------------

model_df <- d %>%
  filter(
    !is.na(lineage),
    !is.na(subregion),
    !is.na(year),
    !is.na(study_accession)
  )

cat("\nComplete-case N:", nrow(model_df), "\n")

if (nrow(model_df) != 3158) {
  stop(
    paste(
      "Expected 3158 complete cases, found",
      nrow(model_df)
    )
  )
}

# Same centering/scaling strategy as previous final model
mean_year <- mean(model_df$year)

model_df <- model_df %>%
  mutate(
    year5 = (year - mean_year) / 5
  )

cat("Mean year:", mean_year, "\n")
cat(
  "Any resistance:",
  sum(model_df$any_resistance), "/",
  nrow(model_df), "\n"
)
cat(
  "MDR-or-higher:",
  sum(model_df$mdr_or_higher), "/",
  nrow(model_df), "\n"
)
cat(
  "Advanced:",
  sum(model_df$advanced_resistance), "/",
  nrow(model_df), "\n"
)
cat(
  "Studies:",
  n_distinct(model_df$study_accession), "\n"
)

# ------------------------------------------------------------
# Crude models
# ------------------------------------------------------------

crude_any <- glm(
  any_resistance ~ lineage,
  data = model_df,
  family = binomial()
)

crude_mdr <- glm(
  mdr_or_higher ~ lineage,
  data = model_df,
  family = binomial()
)

# ------------------------------------------------------------
# Primary adjusted mixed models
# ------------------------------------------------------------

m_any <- glmer(
  any_resistance ~
    lineage + subregion + year5 +
    (1 | study_accession),
  data = model_df,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

m_mdr <- glmer(
  mdr_or_higher ~
    lineage + subregion + year5 +
    (1 | study_accession),
  data = model_df,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# ------------------------------------------------------------
# Advanced resistance
# Same strategy as previous final model:
# bias-reduced logistic regression because L5/L6 have zero events.
# ------------------------------------------------------------

m_adv <- glm(
  advanced_resistance ~
    lineage + subregion + year5,
  data = model_df,
  family = binomial(),
  method = brglmFit,
  type = "AS_mean"
)

# ------------------------------------------------------------
# Tidy lineage ORs
# ------------------------------------------------------------

tidy_or <- function(model, outcome, model_name) {

  cf <- summary(model)$coefficients

  # Both glm/brglm2 and glmer models return the four columns
  # Estimate, Std. Error, z value and Pr(>|z|)
  out <- data.frame(
    term = rownames(cf),
    estimate_log = cf[, "Estimate"],
    SE = cf[, "Std. Error"],
    p_value = cf[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  out <- out %>%
    filter(grepl("^lineage", term)) %>%
    mutate(
      OR = exp(estimate_log),
      conf_low = exp(estimate_log - 1.96 * SE),
      conf_high = exp(estimate_log + 1.96 * SE),
      lineage = sub("^lineage", "", term),
      outcome = outcome,
      model = model_name
    ) %>%
    select(
      outcome,
      model,
      lineage,
      OR,
      conf_low,
      conf_high,
      p_value
    )

  return(out)
}

res <- bind_rows(
  tidy_or(
    crude_any,
    "Any genomic resistance",
    "Crude logistic"
  ),
  tidy_or(
    m_any,
    "Any genomic resistance",
    "Adjusted mixed-effects logistic"
  ),
  tidy_or(
    crude_mdr,
    "MDR-or-higher",
    "Crude logistic"
  ),
  tidy_or(
    m_mdr,
    "MDR-or-higher",
    "Adjusted mixed-effects logistic"
  ),
  tidy_or(
    m_adv,
    "Pre-XDR/XDR",
    "Bias-reduced adjusted logistic"
  )
)

write_tsv(
  res,
  file.path(
    OUTDIR,
    "Table_FINAL_country_corrected_lineage_ORs.tsv"
  )
)

# ------------------------------------------------------------
# Outcome-by-lineage descriptive table
# ------------------------------------------------------------

outcomes <- model_df %>%
  group_by(lineage) %>%
  summarise(
    n = n(),
    any_resistant = sum(any_resistance),
    any_resistant_pct = 100 * mean(any_resistance),

    MDR_or_higher = sum(mdr_or_higher),
    MDR_or_higher_pct = 100 * mean(mdr_or_higher),

    preXDR_XDR = sum(advanced_resistance),
    preXDR_XDR_pct = 100 * mean(advanced_resistance),

    .groups = "drop"
  )

write_tsv(
  outcomes,
  file.path(
    OUTDIR,
    "Table_FINAL_country_corrected_outcomes_by_lineage.tsv"
  )
)

# ------------------------------------------------------------
# Geographic tables
# ------------------------------------------------------------

country_table <- d %>%
  count(
    country = country_final_curated,
    lineage,
    name = "n"
  )

write_tsv(
  country_table,
  file.path(
    OUTDIR,
    "Table_FINAL_country_by_lineage.tsv"
  )
)

subregion_table <- d %>%
  count(
    subregion,
    lineage,
    name = "n"
  )

write_tsv(
  subregion_table,
  file.path(
    OUTDIR,
    "Table_FINAL_subregion_by_lineage.tsv"
  )
)

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

diagnostics <- tibble(
  model = c(
    "Any resistance",
    "MDR-or-higher"
  ),
  N = c(
    nobs(m_any),
    nobs(m_mdr)
  ),
  studies = c(
    n_distinct(model_df$study_accession),
    n_distinct(model_df$study_accession)
  ),
  singular = c(
    isSingular(m_any),
    isSingular(m_mdr)
  ),
  study_RE_SD = c(
    sqrt(as.numeric(VarCorr(m_any)$study_accession)),
    sqrt(as.numeric(VarCorr(m_mdr)$study_accession))
  ),
  AIC = c(
    AIC(m_any),
    AIC(m_mdr)
  )
)

write_tsv(
  diagnostics,
  file.path(
    OUTDIR,
    "Table_FINAL_country_corrected_model_diagnostics.tsv"
  )
)

# Save model objects
saveRDS(
  m_any,
  file.path(OUTDIR, "model_any_resistance.rds")
)

saveRDS(
  m_mdr,
  file.path(OUTDIR, "model_MDR_or_higher.rds")
)

saveRDS(
  m_adv,
  file.path(OUTDIR, "model_preXDR_XDR_brglm.rds")
)

cat("\n========================================\n")
cat("FINAL OR TABLE\n")
cat("========================================\n")
print(res)

cat("\n========================================\n")
cat("DIAGNOSTICS\n")
cat("========================================\n")
print(diagnostics)

cat("\nCompleted successfully.\n")
