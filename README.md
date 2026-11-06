# Study composition and geography shape lineage–drug resistance associations in African MTBC

Reproducibility repository for the manuscript:

**Study composition and geography shape lineage–drug resistance associations in African Mycobacterium tuberculosis complex**

This repository contains the frozen analysis inputs, statistical code, model outputs, audit files and manuscript figures used for the final lineage–drug resistance analyses.

## Study scope

The analysis uses publicly available whole-genome sequencing data from African MTBC isolates processed through a harmonised genomic workflow.

The assembled genomic collection is not a probability sample of tuberculosis cases in Africa. Lineage frequencies and lineage-resistance associations therefore describe this assembled collection and must not be interpreted as population prevalence estimates or intrinsic lineage effects.

## Frozen analysis denominators

| Quantity | N |
|---|---:|
| Sequence libraries assessed | 3,576 |
| Genomes passing sequence QC | 3,491 |
| Principal unmixed L1-L6 cohort | 3,227 |
| Countries represented | 22 |
| Study accessions in full cohort | 112 |
| Isolates with resolved collection year | 3,158 |
| Study accessions in year-resolved cohort | 104 |
| Isolates missing collection year | 69 |
| Advanced-resistance isolates | 86 |

All 69 isolates with missing collection year belong to Lineage 5.

## Repository contents

- `data/processed/` - frozen modelling and public supporting datasets
- `docs/` - statistical model specification
- `environment/` - package versions and R session information
- `figures/` - final manuscript figures
- `outputs/` - model objects, tables, diagnostics, audits and logs
- `scripts/analysis/` - final manuscript analysis workflow
- `scripts/archive_provenance/` - provenance scripts not required for the final statistical reproduction

## Analysis datasets

### manuscript1_analysis_input.tsv

`data/processed/manuscript1_analysis_input.tsv` is the frozen compact modelling dataset used directly by the analysis scripts.

It contains 3,227 isolates and eight variables. Lineage labels used for modelling are L1-L6.

### manuscript1_public_analysis_dataset.tsv

`data/processed/manuscript1_public_analysis_dataset.tsv` is the richer accession-level supporting dataset intended for public reuse and manuscript supplementary data.

It contains the same 3,227 sequencing runs and 31 variables.

The `major_lineage` field uses `lineage1-lineage6`, corresponding directly to modelling labels `L1-L6`.

The modelling and public datasets contain identical isolate sets, with zero discrepancies in study accession, country, subregion, collection year or normalized lineage.

## Principal resistance outcomes

Three binary outcomes are analysed:

1. Any genomic resistance
2. MDR or higher
3. Advanced resistance

Lineage 4 is the reference lineage in all lineage contrasts.

Detailed definitions and statistical equations are provided in `docs/model_specification.md`.

## Reproducibility workflow

Run all commands from the repository root.

Frequentist and audit analyses:

    Rscript scripts/analysis/01_study_subregion_audit.R
    Rscript scripts/analysis/02_frequentist_core_models.R
    Rscript scripts/analysis/03_frequentist_targeted_models.R
    Rscript scripts/analysis/04_finalize_frequentist_OR_table.R
    python3 scripts/analysis/05_resistance_groups_L5_missingness.py

Bayesian analyses:

    micromamba run -n obj1_bayes Rscript scripts/analysis/06_bayesian_models.R
    micromamba run -n obj1_bayes Rscript scripts/analysis/07_finalize_bayesian_outputs.R
    micromamba run -n obj1_bayes Rscript scripts/analysis/08_bayesian_prior_sensitivity.R
    micromamba run -n obj1_bayes Rscript scripts/analysis/09_plot_prior_sensitivity.R

Script 06 refits the principal Bayesian models and is computationally more expensive than the remaining scripts. Saved `.rds` model objects are excluded from the GitHub repository to avoid versioning large binary files and are included with the archived Zenodo release. Running script 06 regenerates the principal Bayesian model objects locally.

## Reproducibility checkpoints

Expected counts:

    Full cohort                     3227
    Collection-year-resolved        3158
    Full-cohort studies              112
    Year-resolved studies            104
    Any resistance                  1405
    MDR or higher                   1053
    Advanced resistance               86

Selected frequentist estimates:

    Any resistance, L6 vs L4: OR = 0.394
    MDR or higher, L1 vs L4: OR = 1.970
    Advanced resistance, L2 vs L4: OR = 7.698
    Advanced + subregion, L2 vs L4: OR = 3.214

Selected Bayesian estimates:

    Any resistance, L6 vs L4: posterior median OR = 0.429
    Advanced resistance, L2 vs L4: posterior median OR = 7.189
    Advanced + geography, L2 vs L4: posterior median OR = 3.049

## Prior sensitivity

The validated prior-sensitivity analysis reproduced all 60 lineage/model/prior combinations exactly from the frozen October 2026 reference output.

Important examples:

    L6 any resistance
    SD 0.5 -> 0.505
    SD 1.0 -> 0.429
    SD 1.5 -> 0.411

    L5 MDR or higher
    SD 0.5 -> 1.850
    SD 1.0 -> 3.589
    SD 1.5 -> 5.121

    L2 advanced resistance
    SD 0.5  -> 6.053
    SD 0.75 -> 7.189
    SD 1.5  -> 7.771

    L2 advanced resistance + geography
    SD 0.5  -> 2.690
    SD 0.75 -> 3.049
    SD 1.5  -> 3.302

The Lineage 5 MDR-or-higher estimate is prior-sensitive.

The Lineage 2 advanced-resistance association remains directionally positive across prior specifications but attenuates substantially after adjustment for African subregion.

## Interpretation constraints

Important limitations include:

- the dataset was purposively assembled from public genomic archives;
- lineage frequencies are not population prevalence estimates;
- adjustment for study does not eliminate ascertainment bias;
- study and geography are strongly entangled;
- 101 of 104 studies in the year-resolved cohort represent only one African subregion;
- patient treatment history and several patient-level determinants were unavailable;
- Lineage 5 has substantial non-random collection-year missingness;
- sparse advanced-resistance outcomes create unstable estimates for lineages with few or zero events.

## Software

Exact package and session information are stored under `environment/`.

Validated Bayesian software includes:

    R        4.5.3
    brms     2.23.0
    cmdstanr 0.9.0
    CmdStan  2.39.0

## Raw sequence data

Raw sequencing reads are not redistributed in this repository.

The underlying sequencing runs originate from public ENA/SRA resources. Run-, sample- and study-level accession identifiers are included in the public supporting dataset.

## Release strategy

1. Complete the final local reproducibility audit.
2. Version-control the repository on GitHub.
3. Tag the manuscript submission release.
4. Archive the tagged release with Zenodo.
5. Add the GitHub URL and Zenodo DOI to the manuscript Data Summary.

## Licence and citation

See `LICENSE` for reuse conditions and `CITATION.cff` for citation metadata.
