#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path.home() / "projects/tb_africa_qc"

INFILE = (
    BASE /
    "objectives/objective1_population_structure/"
    "rerun_20260919/metadata/"
    "objective1_core_L1_L6_resistance_resolved.tsv"
)

OUTDIR = (
    BASE /
    "objectives/objective1_population_structure/"
    "rerun_20261004_country_corrected/metadata"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "objective1_core_3227_country_corrected.tsv"
AUDIT = OUTDIR / "country_backfill_audit.tsv"

# Verified from ENA study titles / BioSample-publication provenance.
STUDY_COUNTRY = {
    "PRJEB27244": "Ghana",
    "PRJEB9545": "Republic of the Congo",
    "PRJEB7727": "Sierra Leone",
    "PRJEB40777": "Cameroon",
    "PRJEB41245": "Nigeria",
    "PRJNA1441484": "Senegal",
    "PRJNA233395": "Uganda",
    "PRJNA233396": "Uganda",
    "PRJNA233397": "Uganda",
    "PRJNA233398": "Uganda",
    "PRJNA233399": "Uganda",
}

COUNTRY_SUBREGION = {
    "Ghana": "West Africa",
    "Sierra Leone": "West Africa",
    "Nigeria": "West Africa",
    "Senegal": "West Africa",
    "Cameroon": "Central Africa",
    "Republic of the Congo": "Central Africa",
    "Uganda": "East Africa",
}

df = pd.read_csv(INFILE, sep="\t", dtype=str)

if len(df) != 3227:
    raise RuntimeError(f"Expected 3227 rows, found {len(df)}")

if df["sample_id"].nunique() != 3227:
    raise RuntimeError("sample_id is not unique")

# Preserve original fields
for c in ["country_final_curated", "country_clean", "country"]:
    if c in df.columns:
        df[c + "_before_backfill"] = df[c]

if "is_country_missing" in df.columns:
    df["is_country_missing_before_backfill"] = df["is_country_missing"]

if "is_africa_confirmed" in df.columns:
    df["is_africa_confirmed_before_backfill"] = df["is_africa_confirmed"]

# Identify canonical country column
country_col = next(
    (c for c in ["country_final_curated", "country_clean", "country"]
     if c in df.columns),
    None
)

if country_col is None:
    country_col = "country_final_curated"
    df[country_col] = ""

# Identify subregion column if present
subregion_col = next(
    (c for c in [
        "subregion",
        "african_subregion",
        "subregion_final",
        "who_subregion"
    ] if c in df.columns),
    None
)

if subregion_col is None:
    subregion_col = "subregion"
    df[subregion_col] = ""

df["country_backfilled"] = False
df["country_backfill_source"] = ""
df["country_backfill_basis"] = ""

audit_rows = []

for idx, row in df.iterrows():

    study = str(row.get("study_accession", "")).strip()

    current_country = str(row.get(country_col, "")).strip()

    is_missing = (
        current_country == ""
        or current_country.lower() in {"nan", "none", "na"}
    )

    if is_missing and study in STUDY_COUNTRY:

        new_country = STUDY_COUNTRY[study]
        new_subregion = COUNTRY_SUBREGION[new_country]

        df.at[idx, country_col] = new_country
        df.at[idx, subregion_col] = new_subregion

        if "country_clean" in df.columns:
            df.at[idx, "country_clean"] = new_country

        if "country_final_curated" in df.columns:
            df.at[idx, "country_final_curated"] = new_country

        if "is_country_missing" in df.columns:
            df.at[idx, "is_country_missing"] = "False"

        if "is_africa_confirmed" in df.columns:
            df.at[idx, "is_africa_confirmed"] = "True"

        if "is_confirmed_non_africa" in df.columns:
            df.at[idx, "is_confirmed_non_africa"] = "False"

        df.at[idx, "country_backfilled"] = True
        df.at[idx, "country_backfill_source"] = study
        df.at[idx, "country_backfill_basis"] = (
            "ENA study title/BioSample/publication provenance audit"
        )

        audit_rows.append({
            "sample_id": row["sample_id"],
            "study_accession": study,
            "lineage": row.get(
                "lineage_model",
                row.get("lineage_final", "")
            ),
            "old_country": current_country,
            "new_country": new_country,
            "new_subregion": new_subregion,
            "collection_year": row.get(
                "collection_year_clean",
                row.get("collection_year", "")
            ),
        })

audit = pd.DataFrame(audit_rows)

print("=" * 80)
print("COUNTRY BACKFILL AUDIT")
print("=" * 80)

print("Input rows:", len(df))
print("Backfilled:", len(audit))

if len(audit) != 74:
    raise RuntimeError(
        f"Expected 74 backfilled isolates, found {len(audit)}"
    )

print("\nBackfilled by study:")
print(
    audit["study_accession"]
    .value_counts()
    .to_string()
)

print("\nBackfilled by country:")
print(
    audit["new_country"]
    .value_counts()
    .to_string()
)

print("\nBackfilled by lineage:")
print(
    audit["lineage"]
    .value_counts()
    .to_string()
)

# Check missing country after backfill
remaining_missing = (
    df[country_col]
    .fillna("")
    .astype(str)
    .str.strip()
    .isin(["", "nan", "None", "NA"])
)

print("\nRemaining missing countries:", remaining_missing.sum())

# Year availability
year_col = next(
    (c for c in [
        "collection_year_clean",
        "collection_year",
        "decimal_date"
    ] if c in df.columns),
    None
)

if year_col:
    yr = pd.to_numeric(df[year_col], errors="coerce")
    print("With resolved year/date:", yr.notna().sum())
    print("Without resolved year/date:", yr.isna().sum())

df.to_csv(OUTFILE, sep="\t", index=False)
audit.to_csv(AUDIT, sep="\t", index=False)

print("\nWrote:")
print(OUTFILE)
print(AUDIT)
