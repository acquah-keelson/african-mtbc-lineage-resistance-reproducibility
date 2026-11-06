# Final manuscript analysis scripts

The scripts under `analysis/` constitute the final reproducible statistical workflow for Manuscript 1.

Historical or provenance scripts are kept separately under `archive_provenance/` and are not part of the final modelling workflow.

## Script order

### 01_study_subregion_audit.R

Audits study and African-subregion structure and verifies the extent of study-geography entanglement.

### 02_frequentist_core_models.R

Fits the primary frequentist lineage models for:

- any genomic resistance;
- MDR or higher.

It also produces the principal descriptive lineage and resistance tables.

### 03_frequentist_targeted_models.R

Fits targeted sensitivity analyses, including:

- advanced-resistance bias-reduced logistic regression;
- advanced resistance adjusted for African subregion;
- full-cohort sensitivity models where applicable.

### 04_finalize_frequentist_OR_table.R

Consolidates the final frequentist lineage odds ratios into reportable machine-readable tables.

### 05_resistance_groups_L5_missingness.py

Reconstructs and audits the eight mutually exclusive genomic resistance groups and evaluates Lineage 5 collection-year missingness.

### 06_bayesian_models.R

Fits the six validated Bayesian models:

- primary any resistance;
- full-cohort any resistance;
- primary MDR or higher;
- full-cohort MDR or higher;
- primary advanced resistance;
- advanced resistance plus geography.

This is the computationally expensive Bayesian fitting step.

### 07_finalize_bayesian_outputs.R

Loads the saved Bayesian model objects and generates:

- posterior lineage odds-ratio summaries;
- MCMC diagnostics;
- study random-effect standard-deviation summaries.

### 08_bayesian_prior_sensitivity.R

Runs or reloads Bayesian prior-sensitivity models.

For any resistance and MDR-or-higher models, the population-level coefficient prior width is varied.

For advanced-resistance models, lineage-specific prior width is varied while other fixed-effect priors retain the validated specification.

### 09_plot_prior_sensitivity.R

Generates the final Bayesian prior-sensitivity forest plots from the completed sensitivity table.

## Running the workflow

Run from the repository root.

Frequentist and audit steps:

    Rscript scripts/analysis/01_study_subregion_audit.R
    Rscript scripts/analysis/02_frequentist_core_models.R
    Rscript scripts/analysis/03_frequentist_targeted_models.R
    Rscript scripts/analysis/04_finalize_frequentist_OR_table.R
    python3 scripts/analysis/05_resistance_groups_L5_missingness.py

Bayesian steps:

    micromamba run -n obj1_bayes Rscript scripts/analysis/06_bayesian_models.R
    micromamba run -n obj1_bayes Rscript scripts/analysis/07_finalize_bayesian_outputs.R
    micromamba run -n obj1_bayes Rscript scripts/analysis/08_bayesian_prior_sensitivity.R
    micromamba run -n obj1_bayes Rscript scripts/analysis/09_plot_prior_sensitivity.R

## Frozen audit conditions

The final workflow must reproduce:

    Principal cohort                              3227
    Collection-year-adjusted cohort               3158
    Full-cohort studies                            112
    Year-resolved studies                          104
    L2 advanced-resistance events                   39
    Total advanced-resistance events                86

Selected frequentist checkpoints:

    L2 advanced, year-adjusted OR                 7.698
    L2 advanced + subregion OR                    3.214
    L6 any-resistance adjusted OR                 0.394

Selected Bayesian checkpoints:

    L6 any-resistance posterior median OR         0.429
    L5 MDR-or-higher primary OR                   3.589
    L5 MDR-or-higher full-cohort OR               2.458
    L2 advanced posterior median OR               7.189
    L2 advanced + geography OR                    3.049

Prior-sensitivity analysis must contain 60 lineage/model/prior combinations and zero divergent transitions in the validated run.

If these quantities do not reproduce, the discrepancy must be reconciled before public release.

## Model specification

Formal model equations, priors and interpretation constraints are documented in:

`docs/model_specification.md`
