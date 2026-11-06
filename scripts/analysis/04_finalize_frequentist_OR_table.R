suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lme4)
})

BASE <- "."

META <- file.path(
  BASE,
  "data",
  "processed",
  "manuscript1_analysis_input.tsv"
)

OLD_OR <- file.path(
  BASE,
  "outputs",
  "tables",
  "frequentist",
  "Table_05_FINAL_lineage_ORs_all_models.tsv"
)

OUTDIR <- file.path(
  BASE,
  "outputs",
  "tables",
  "frequentist"
)

d <- read_tsv(META, show_col_types = FALSE)

stopifnot(nrow(d) == 3227)

d <- d %>%
  mutate(
    lineage = factor(
      lineage_model,
      levels = c("L4","L1","L2","L3","L5","L6")
    ),

    any_resistance = if_else(
      resistance_category_clean == "no_resistance_predicted",
      0L, 1L
    ),

    mdr_or_higher = if_else(
      resistance_category_clean %in%
        c("MDR","pre_XDR","XDR"),
      1L, 0L
    )
  )

# ------------------------------------------------------------
# Full-cohort crude models: N = 3227
# ------------------------------------------------------------

crude_any <- glm(
  any_resistance ~ lineage,
  data = d,
  family = binomial()
)

crude_mdr <- glm(
  mdr_or_higher ~ lineage,
  data = d,
  family = binomial()
)

extract_crude <- function(model, outcome) {

  cf <- summary(model)$coefficients

  out <- data.frame(
    term = rownames(cf),
    beta = cf[, "Estimate"],
    SE = cf[, "Std. Error"],
    p_value = cf[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  out %>%
    filter(grepl("^lineage", term)) %>%
    mutate(
      outcome = outcome,
      model = "Crude logistic: full N=3227",
      lineage = sub("^lineage", "", term),
      OR = exp(beta),
      conf_low = exp(beta - 1.96 * SE),
      conf_high = exp(beta + 1.96 * SE)
    ) %>%
    select(
      outcome, model, lineage,
      beta, SE, OR,
      conf_low, conf_high, p_value
    )
}

crude_final <- bind_rows(
  extract_crude(
    crude_any,
    "Any genomic resistance"
  ),
  extract_crude(
    crude_mdr,
    "MDR-or-higher"
  )
)

old <- read_tsv(
  OLD_OR,
  show_col_types = FALSE
)

# Keep only primary adjusted models.
# Drop old crude rows and the discarded geographic sensitivity.
adjusted <- old %>%
  filter(
    model == "Primary: lineage + year + study RE" |
    model == "Primary bias-reduced: lineage + year"
  )

final <- bind_rows(
  crude_final,
  adjusted
)

# ------------------------------------------------------------
# Suppress invalid L5 advanced-resistance result
# ------------------------------------------------------------

final <- final %>%
  mutate(
    estimation_status = case_when(
      outcome == "Pre-XDR/XDR" &
        lineage == "L5" ~
        "Not estimable: zero events and bias-reduced model did not converge",

      TRUE ~ "Estimated"
    ),

    beta = if_else(
      outcome == "Pre-XDR/XDR" &
        lineage == "L5",
      NA_real_, beta
    ),

    SE = if_else(
      outcome == "Pre-XDR/XDR" &
        lineage == "L5",
      NA_real_, SE
    ),

    OR = if_else(
      outcome == "Pre-XDR/XDR" &
        lineage == "L5",
      NA_real_, OR
    ),

    conf_low = if_else(
      outcome == "Pre-XDR/XDR" &
        lineage == "L5",
      NA_real_, conf_low
    ),

    conf_high = if_else(
      outcome == "Pre-XDR/XDR" &
        lineage == "L5",
      NA_real_, conf_high
    ),

    p_value = if_else(
      outcome == "Pre-XDR/XDR" &
        lineage == "L5",
      NA_real_, p_value
    )
  )

write_tsv(
  final,
  file.path(
    OUTDIR,
    "Table_05_FINAL_REPORTABLE_lineage_ORs.tsv"
  )
)

cat("============================================================\n")
cat("FINAL REPORTABLE OBJECTIVE 1 OR TABLE\n")
cat("============================================================\n")
cat("Crude cohort N: 3227\n")
cat("Adjusted cohort N: 3158\n")
cat("L5 Pre-XDR/XDR: suppressed as not estimable\n\n")

print(final, row.names = FALSE)
