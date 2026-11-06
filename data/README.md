# Data files

This directory contains the frozen datasets and machine-readable values used for Manuscript 1.

## processed/manuscript1_analysis_input.tsv

Authoritative compact modelling input used directly by the statistical analysis scripts.

Rows: 3,227 isolates.

Columns:

- `sample_id`
- `study_accession`
- `country_final_curated`
- `subregion`
- `collection_year_clean`
- `lineage_model`
- `resistance_category_clean`
- `resistance_drugs`

`lineage_model` uses labels `L1`-`L6`.

The collection-year-adjusted analyses use the 3,158 rows with non-missing `collection_year_clean`.

## processed/manuscript1_public_analysis_dataset.tsv

Public accession-level supporting dataset.

Rows: 3,227 sequencing runs.

Columns: 31.

It contains public run, sample and study accessions; curated geographic metadata; sequence-quality summaries; lineage assignments; genomic resistance classifications; individual resistance indicators; and resistance genes and variants.

The `major_lineage` field uses `lineage1`-`lineage6`, corresponding directly to modelling labels `L1`-`L6`.

The public dataset and compact modelling dataset have been cross-audited and contain identical isolate sets, with zero discrepancies in study accession, country, subregion, collection year or normalized lineage.

## processed/data_dictionary.tsv

Variable definitions for the public supporting dataset.

## processed/Figure2_forest_plot_values.tsv

Machine-readable values underlying the principal lineage-resistance forest plot.

## Raw sequencing data

Raw FASTQ files are not redistributed here.

The underlying sequence data are publicly available through ENA/SRA. Public run, sample and study accessions are retained in `manuscript1_public_analysis_dataset.tsv`.

## Dataset scope

The repository contains genomic and study-level metadata used for this analysis.

No patient names or direct personal identifiers are included.

## Frozen denominators

Expected dataset-level counts are:

    Principal cohort                 3227
    Countries                          22
    Full-cohort study accessions      112
    Collection-year-resolved         3158
    Year-resolved study accessions    104
    Missing collection year            69

All 69 isolates with missing collection year belong to Lineage 5.

## Resistance outcomes

The main binary outcomes are:

- any genomic resistance;
- MDR or higher;
- advanced resistance.

The public dataset also contains the final mutually exclusive eight-group genomic resistance classification.

Detailed definitions are provided in `../docs/model_specification.md`.
