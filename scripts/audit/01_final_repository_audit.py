#!/usr/bin/env python3

from pathlib import Path
import csv
import sys

ROOT = Path(__file__).resolve().parents[2]

errors = []

def fail(msg):
    errors.append(msg)

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

model_path = ROOT / "data/processed/manuscript1_analysis_input.tsv"
public_path = ROOT / "data/processed/manuscript1_public_analysis_dataset.tsv"

if not model_path.exists():
    fail(f"Missing {model_path}")
if not public_path.exists():
    fail(f"Missing {public_path}")

if not errors:
    model = read_tsv(model_path)
    public = read_tsv(public_path)

    if len(model) != 3227:
        fail(f"Model input rows = {len(model)}, expected 3227")

    if len(public) != 3227:
        fail(f"Public dataset rows = {len(public)}, expected 3227")

    model_ids = {r["sample_id"] for r in model}
    public_ids = {r["run_accession"] for r in public}

    if model_ids != public_ids:
        fail("Model/public isolate ID sets differ")

    dated = sum(bool(r["collection_year_clean"].strip()) for r in model)
    if dated != 3158:
        fail(f"Dated isolates = {dated}, expected 3158")

    studies_full = len({r["study_accession"] for r in model})
    if studies_full != 112:
        fail(f"Full studies = {studies_full}, expected 112")

    studies_dated = len({
        r["study_accession"]
        for r in model
        if r["collection_year_clean"].strip()
    })
    if studies_dated != 104:
        fail(f"Dated studies = {studies_dated}, expected 104")

    lineage_expected = {
        "L1": 247,
        "L2": 257,
        "L3": 372,
        "L4": 1987,
        "L5": 94,
        "L6": 270,
    }

    lineage_counts = {}
    for r in model:
        lineage_counts[r["lineage_model"]] = (
            lineage_counts.get(r["lineage_model"], 0) + 1
        )

    if lineage_counts != lineage_expected:
        fail(
            f"Lineage counts differ: {lineage_counts}"
        )

    missing_year = [
        r for r in model
        if not r["collection_year_clean"].strip()
    ]

    if len(missing_year) != 69:
        fail(
            f"Missing collection year = {len(missing_year)}, expected 69"
        )

    if any(r["lineage_model"] != "L5" for r in missing_year):
        fail("Not all missing-year isolates are L5")

required_files = [
    "docs/model_specification.md",
    "README.md",
    "data/README.md",
    "scripts/README.md",
    "checksums/SHA256SUMS.txt",
    "outputs/tables/frequentist/Table_05_FINAL_REPORTABLE_lineage_ORs.tsv",
    "outputs/bayesian/tables/Table_BAYES_FINAL_01_posterior_lineage_ORs.tsv",
    "outputs/bayesian/tables/Table_BAYES_FINAL_02_MCMC_diagnostics.tsv",
    "outputs/bayesian/prior_sensitivity/tables/Table_OCT2026_prior_sensitivity.tsv",
]

for rel in required_files:
    if not (ROOT / rel).exists():
        fail(f"Missing required file: {rel}")

if errors:
    print("REPOSITORY AUDIT: FAIL")
    for e in errors:
        print(" -", e)
    sys.exit(1)

print("REPOSITORY AUDIT: PASS")
print("Principal cohort: 3227")
print("Year-resolved cohort: 3158")
print("Full studies: 112")
print("Year-resolved studies: 104")
print("Missing collection year: 69, all L5")
print("Lineages: L1=247 L2=257 L3=372 L4=1987 L5=94 L6=270")
