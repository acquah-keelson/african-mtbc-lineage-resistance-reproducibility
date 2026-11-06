suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

BASE <- "."

F <- file.path(
  BASE,
  "data",
  "processed",
  "manuscript1_analysis_input.tsv"
)

d <- read_tsv(F, show_col_types = FALSE) %>%
  mutate(
    year = as.numeric(collection_year_clean)
  ) %>%
  filter(!is.na(year))

cat("============================================================\n")
cat("STUDY × SUBREGION IDENTIFIABILITY AUDIT\n")
cat("============================================================\n")

cat("Complete-case isolates:", nrow(d), "\n")
cat("Observed studies:", n_distinct(d$study_accession), "\n")
cat("Subregions:", n_distinct(d$subregion), "\n\n")

x <- d %>%
  distinct(study_accession, subregion) %>%
  count(study_accession, name = "n_subregions")

cat("Studies represented in exactly 1 subregion:",
    sum(x$n_subregions == 1), "\n")

cat("Studies represented in >1 subregion:",
    sum(x$n_subregions > 1), "\n")

cat("\nDistribution of number of subregions per study:\n")
print(table(x$n_subregions))

cat("\nStudies spanning >1 subregion:\n")
multi <- x %>%
  filter(n_subregions > 1)

if (nrow(multi) == 0) {
  cat("NONE\n")
} else {
  print(multi)

  print(
    d %>%
      semi_join(multi, by = "study_accession") %>%
      count(study_accession, subregion) %>%
      arrange(study_accession, subregion)
  )
}

cat("\n============================================================\n")
cat("STUDY × SUBREGION SAMPLE COUNTS\n")
cat("============================================================\n")

study_region <- d %>%
  count(study_accession, subregion, name = "n") %>%
  arrange(subregion, desc(n))

print(study_region, n = Inf)

cat("\n============================================================\n")
cat("LINEAGE × SUBREGION × STUDY SUPPORT\n")
cat("============================================================\n")

support <- d %>%
  group_by(lineage_model, subregion) %>%
  summarise(
    isolates = n(),
    studies = n_distinct(study_accession),
    .groups = "drop"
  )

print(support, n = Inf)

