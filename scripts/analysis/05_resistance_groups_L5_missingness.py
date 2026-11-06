#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

BASE = Path(".")

INFILE = (
    BASE
    / "data/processed/manuscript1_analysis_input.tsv"
)

OUTDIR = (
    BASE
    / "outputs/tables/sensitivity"
)

AUDITDIR = (
    BASE
    / "outputs/audits"
)

OUTDIR.mkdir(parents=True, exist_ok=True)
AUDITDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INFILE, sep="\t", dtype=str)

if len(df) != 3227:
    raise RuntimeError(f"Expected 3227 isolates, found {len(df)}")

required = [
    "sample_id",
    "lineage_model",
    "resistance_category_clean",
    "resistance_drugs",
    "collection_year_clean",
    "country_final_curated",
    "subregion",
]

missing = [x for x in required if x not in df.columns]

if missing:
    raise RuntimeError(f"Missing required columns: {missing}")


# ============================================================
# HELPERS
# ============================================================

def clean_missing(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if x.lower() in {"", "nan", "none", "na", "n/a", "null"}:
        return ""
    return x


def parse_drugs(x):
    """
    Parse resistance_drugs while retaining unique drug names.
    Existing project files have used comma/semicolon/pipe-like
    delimiters across different stages, so tolerate all here.
    """
    x = clean_missing(x)

    if not x:
        return []

    # remove simple list punctuation/quotes
    x = x.replace("[", "").replace("]", "")
    x = x.replace("'", "").replace('"', "")

    parts = re.split(r"[;,|]+", x)

    parts = [
        p.strip()
        for p in parts
        if p.strip()
    ]

    return sorted(set(parts))


def norm_token(x):
    return (
        str(x)
        .strip()
        .lower()
        .replace("_", " ")
        .replace("-", " ")
    )


def is_inh(drugs):
    d = [norm_token(x) for x in drugs]

    return any(
        x == "inh"
        or x == "isoniazid"
        or "isoniazid" in x
        for x in d
    )


def is_rif(drugs):
    d = [norm_token(x) for x in drugs]

    return any(
        x in {"rif", "rifampicin", "rifampin", "rifapentine"}
        or "rifampicin" in x
        or "rifampin" in x
        or "rifapentine" in x
        for x in d
    )


# ============================================================
# EIGHT-GROUP CLASSIFICATION
#
# Mirrors the documented hierarchy of the existing Objective 1
# mutually-exclusive resistance-group workflow:
#
# 1 gXDR-TB
# 2 gPre-XDR-TB
# 3 gMDR-TB
# 4 gRR-TB
# 5 gHR-TB
# 6 gMono-other
# 7 gPoly-resistance
# 8 gDS-TB
# ============================================================

def classify(row):

    cat = clean_missing(
        row["resistance_category_clean"]
    )

    drugs = parse_drugs(
        row["resistance_drugs"]
    )

    n = len(drugs)

    inh = is_inh(drugs)
    rif = is_rif(drugs)

    # Highest-resistance categories first
    if cat == "XDR":
        return "gXDR-TB"

    if cat == "pre_XDR":
        return "gPre-XDR-TB"

    if cat == "MDR":
        return "gMDR-TB"

    # Genomically susceptible
    if cat == "no_resistance_predicted" or n == 0:
        return "gDS-TB"

    # Rifampicin resistance without isoniazid resistance.
    # This is NOT restricted to rifampicin mono-resistance.
    if rif and not inh:
        return "gRR-TB"

    # Isoniazid-only
    if inh and not rif and n == 1:
        return "gHR-TB"

    # Exactly one drug other than INH/RIF
    if not inh and not rif and n == 1:
        return "gMono-other"

    # Multiple resistance drugs without meeting MDR/pre-XDR/XDR/RR
    return "gPoly-resistance"


df["genomic_resistance_group"] = df.apply(
    classify,
    axis=1
)


# ============================================================
# VALIDATE GROUP ASSIGNMENTS
# ============================================================

ORDER = [
    "gDS-TB",
    "gHR-TB",
    "gRR-TB",
    "gMono-other",
    "gPoly-resistance",
    "gMDR-TB",
    "gPre-XDR-TB",
    "gXDR-TB",
]

bad = sorted(
    set(df["genomic_resistance_group"])
    - set(ORDER)
)

if bad:
    raise RuntimeError(
        f"Unexpected groups: {bad}"
    )

if df["genomic_resistance_group"].isna().any():
    raise RuntimeError(
        "Some isolates were not classified"
    )


# ============================================================
# TABLE 1: EIGHT-GROUP DISTRIBUTION, N=3227
# ============================================================

gcounts = (
    df["genomic_resistance_group"]
    .value_counts()
    .reindex(ORDER, fill_value=0)
    .rename_axis("Genomic_resistance_group")
    .reset_index(name="No_isolates")
)

gcounts["Percent"] = (
    100 * gcounts["No_isolates"] / len(df)
)

gcounts.to_csv(
    OUTDIR / "Table_04A_eight_genomic_resistance_groups_3227.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# TABLE 2: LINEAGE × EIGHT GROUPS
# ============================================================

lineage_tab = pd.crosstab(
    df["lineage_model"],
    df["genomic_resistance_group"]
)

lineage_tab = lineage_tab.reindex(
    columns=ORDER,
    fill_value=0
)

lineage_tab["Total"] = lineage_tab.sum(axis=1)

lineage_tab.to_csv(
    OUTDIR / "Table_04B_lineage_by_eight_resistance_groups_3227.tsv",
    sep="\t"
)


# ============================================================
# TABLE 3: CROSSWALK BETWEEN SEVEN-LEVEL PIPELINE VARIABLE
#          AND EIGHT-GROUP PRESENTATION
# ============================================================

crosswalk = pd.crosstab(
    df["resistance_category_clean"],
    df["genomic_resistance_group"]
)

crosswalk = crosswalk.reindex(
    columns=ORDER,
    fill_value=0
)

crosswalk.to_csv(
    OUTDIR / "Table_04C_crosswalk_seven_to_eight_groups_3227.tsv",
    sep="\t"
)


# ============================================================
# TABLE 4: SUBREGION × EIGHT GROUPS
# ============================================================

subregion_tab = pd.crosstab(
    df["subregion"],
    df["genomic_resistance_group"]
)

subregion_tab = subregion_tab.reindex(
    columns=ORDER,
    fill_value=0
)

subregion_tab["Total"] = subregion_tab.sum(axis=1)

subregion_tab.to_csv(
    OUTDIR / "Table_04D_subregion_by_eight_resistance_groups_3227.tsv",
    sep="\t"
)


# ============================================================
# L5 COLLECTION-YEAR MISSINGNESS AUDIT
# ============================================================

# Binary outcomes exactly as used in the final regression work
df["any_resistance"] = (
    df["resistance_category_clean"]
    != "no_resistance_predicted"
).astype(int)

df["mdr_or_higher"] = (
    df["resistance_category_clean"]
    .isin(["MDR", "pre_XDR", "XDR"])
).astype(int)

df["advanced_resistance"] = (
    df["resistance_category_clean"]
    .isin(["pre_XDR", "XDR"])
).astype(int)

df["has_collection_year"] = (
    pd.to_numeric(
        df["collection_year_clean"],
        errors="coerce"
    )
    .notna()
)

l5 = df.loc[
    df["lineage_model"] == "L5"
].copy()

if len(l5) != 94:
    raise RuntimeError(
        f"Expected 94 L5 isolates, found {len(l5)}"
    )

l5_missingness = (
    l5
    .groupby(
        "has_collection_year",
        dropna=False
    )
    .agg(
        n=("sample_id", "size"),
        any_resistant=("any_resistance", "sum"),
        MDR_or_higher=("mdr_or_higher", "sum"),
        preXDR_XDR=("advanced_resistance", "sum"),
    )
    .reset_index()
)

l5_missingness["any_resistant_pct"] = (
    100
    * l5_missingness["any_resistant"]
    / l5_missingness["n"]
)

l5_missingness["MDR_or_higher_pct"] = (
    100
    * l5_missingness["MDR_or_higher"]
    / l5_missingness["n"]
)

l5_missingness["preXDR_XDR_pct"] = (
    100
    * l5_missingness["preXDR_XDR"]
    / l5_missingness["n"]
)

l5_missingness.to_csv(
    OUTDIR / "Table_08_L5_collection_year_missingness_resistance.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# L5 MISSING-YEAR STUDY COMPOSITION
# ============================================================

l5_missing_studies = (
    l5.loc[~l5["has_collection_year"]]
    .groupby(
        [
            "study_accession",
            "country_final_curated",
            "subregion",
        ],
        dropna=False,
    )
    .agg(
        n=("sample_id", "size"),
        any_resistant=("any_resistance", "sum"),
        MDR_or_higher=("mdr_or_higher", "sum"),
        preXDR_XDR=("advanced_resistance", "sum"),
    )
    .reset_index()
)

l5_missing_studies[
    "any_resistant_pct"
] = (
    100
    * l5_missing_studies["any_resistant"]
    / l5_missing_studies["n"]
)

l5_missing_studies[
    "MDR_or_higher_pct"
] = (
    100
    * l5_missing_studies["MDR_or_higher"]
    / l5_missing_studies["n"]
)

l5_missing_studies.to_csv(
    OUTDIR / "Table_09_L5_missing_year_by_study.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# AUDIT
# ============================================================

with open(
    AUDITDIR / "09_final_eight_group_and_L5_missingness_audit.txt",
    "w"
) as out:

    def w(x=""):
        out.write(str(x) + "\n")

    w("=" * 72)
    w("FINAL EIGHT-GROUP + L5 DATE-MISSINGNESS AUDIT")
    w("=" * 72)

    w()
    w(f"Total cohort: {len(df)}")

    w()
    w("Eight-group counts:")
    w(gcounts.to_string(index=False))

    w()
    w("Crosswalk from resistance_category_clean:")
    w(crosswalk.to_string())

    w()
    w("L5 collection-year availability:")
    w(l5_missingness.to_string(index=False))

    w()
    w("L5 missing-year studies:")
    w(l5_missing_studies.to_string(index=False))

    w()
    w("Checks:")
    w(
        f"Eight-group sum = "
        f"{gcounts['No_isolates'].sum()}"
    )
    w(
        f"Lineage table sum = "
        f"{lineage_tab['Total'].sum()}"
    )
    w(
        f"Subregion table sum = "
        f"{subregion_tab['Total'].sum()}"
    )


# ============================================================
# SCREEN OUTPUT
# ============================================================

print("=" * 72)
print("EIGHT-GROUP DISTRIBUTION — N=3227")
print("=" * 72)
print(gcounts.to_string(index=False))

print()
print("=" * 72)
print("SEVEN-LEVEL → EIGHT-GROUP CROSSWALK")
print("=" * 72)
print(crosswalk.to_string())

print()
print("=" * 72)
print("L5 COLLECTION-YEAR MISSINGNESS")
print("=" * 72)
print(l5_missingness.to_string(index=False))

print()
print("=" * 72)
print("L5 MISSING-YEAR STUDIES")
print("=" * 72)
print(l5_missing_studies.to_string(index=False))

print()
print("All eight-group isolates:",
      gcounts["No_isolates"].sum())
print("DONE")
