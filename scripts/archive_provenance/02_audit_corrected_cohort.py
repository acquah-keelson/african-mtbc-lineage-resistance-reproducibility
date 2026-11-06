#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path.home() / "projects/tb_africa_qc"

F = (
    BASE /
    "objectives/objective1_population_structure/"
    "rerun_20261004_country_corrected/metadata/"
    "objective1_core_3227_country_corrected.tsv"
)

df = pd.read_csv(F, sep="\t", dtype=str)

print("=" * 80)
print("OBJECTIVE 1 CORRECTED COHORT AUDIT")
print("=" * 80)

print("Rows:", len(df))
print("Unique samples:", df.sample_id.nunique())

lin = next(
    c for c in ["lineage_model", "lineage_final", "tbprof_main_lineage"]
    if c in df.columns
)

year = next(
    c for c in ["collection_year_clean", "collection_year", "decimal_date"]
    if c in df.columns
)

country = next(
    c for c in ["country_final_curated", "country_clean", "country"]
    if c in df.columns
)

subregion = next(
    c for c in ["subregion", "african_subregion", "subregion_final"]
    if c in df.columns
)

print("\nLineage:")
print(df[lin].value_counts())

print("\nCountry:")
print(df[country].value_counts())

print("\nSubregion:")
print(df[subregion].value_counts())

yr = pd.to_numeric(df[year], errors="coerce")

print("\nWith year:", yr.notna().sum())
print("Without year:", yr.isna().sum())

cc = (
    df[lin].notna()
    & df[country].notna()
    & df[subregion].notna()
    & yr.notna()
)

print("\nExpected geography+year complete cases:", cc.sum())

print("\nMissing-year isolates by lineage:")
print(df.loc[yr.isna(), lin].value_counts())

print("\nBackfilled isolates by study:")
print(
    df.loc[df["country_backfilled"].astype(str) == "True",
           "study_accession"]
    .value_counts()
)
