#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RNA-seq Metadata Harmonization Tool
-----------------------------------
- Reads GEO *series_matrix.txt.gz* files
- Extracts sample accessions (GSM) and titles
- Classifies treatment/control using simple keyword rules
- Performs QC (missing values, counts)
- Saves a harmonized CSV + two PNG plots

Run:
    python src/harmonize_metadata.py
"""

import os
import re
import gzip
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend so plots save without a GUI
import matplotlib.pyplot as plt


# -----------------------------
# USER: list your series_matrix files here
# -----------------------------
DATASETS = [
    "metadata/GSE27651_series_matrix.txt.gz",
    "metadata/GSE36668_series_matrix.txt.gz",
    "metadata/GSE66957_series_matrix.txt.gz",
    "metadata/GSE54388_series_matrix.txt.gz",
]

# Output folder & file
RESULTS_DIR = "results"
OUTPUT_CSV  = os.path.join(RESULTS_DIR, "harmonized_metadata_with_treatment.csv")
PLOT_STUDY = os.path.join(RESULTS_DIR, "sample_count_per_study.png")
PLOT_TREAT = os.path.join(RESULTS_DIR, "sample_count_per_treatment.png")


def ensure_dirs():
    os.makedirs(RESULTS_DIR, exist_ok=True)


def study_id_from_path(path: str) -> str:
    """
    Extracts GSE accession like 'GSE66957' from a file path.
    """
    m = re.search(r"(GSE\d+)", os.path.basename(path))
    return m.group(1) if m else "UNKNOWN_STUDY"


def split_geo_values(line: str):
    """
    GEO lines look like:
    !Sample_title    "Normal-ovarian sample1"    "Ovarian sample1" ...
    We split on tabs, drop the first token, strip quotes and whitespace.
    """
    parts = line.strip().split("\t")[1:]
    return [p.strip().strip('"') for p in parts]


def parse_series_matrix_gz(path: str):
    """
    Read only the metadata lines from a GEO series_matrix.txt.gz:
    - !Sample_title
    - !Sample_geo_accession
    Returns dict with lists (may be empty if not found).
    """
    titles = []
    accessions = []

    with gzip.open(path, "rt", encoding="latin-1", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = split_geo_values(line)
            elif line.startswith("!Sample_geo_accession"):
                accessions = split_geo_values(line)

    return {
        "titles": titles,
        "accessions": accessions
    }


def classify_treatment(title: str) -> str:
    """
    VERY SIMPLE keyword-based classification.
    Customize as you like (add rules for your studies).
    """
    t = title.lower()

    control_kw = ["normal", "control", "healthy", "non-tumor", "non tumor", "untreated", "vehicle"]
    treated_kw = ["tumor", "carcinoma", "cancer", "treated", "ovarian cancer"]

    if any(k in t for k in control_kw):
        return "control"
    if any(k in t for k in treated_kw):
        return "treated"
    return "unknown"


def load_expression_columns(path: str):
    """
    Fallback: read the series_matrix as a table (skipping metadata lines)
    and return the column names after 'ID_REF' -> these are usually GSM IDs.
    """
    df = pd.read_csv(path, sep="\t", comment="!", low_memory=False)
    cols = list(df.columns)
    # Typically first column is ID_REF; remaining are sample IDs.
    if len(cols) > 1:
        return cols[1:]
    return []


def process_one_file(path: str) -> pd.DataFrame:
    """
    For a single series_matrix file:
    - Parse titles & accessions from metadata lines
    - If missing, fallback to expression columns for sample IDs
    - Create a per-study DataFrame with sample_id, study, treatment
    """
    study = study_id_from_path(path)

    parsed = parse_series_matrix_gz(path)
    titles = parsed["titles"]
    accessions = parsed["accessions"]

    # If no accessions parsed, fallback to expression columns
    if not accessions:
        accessions = load_expression_columns(path)

    # Titles length may not match accessions length -> align carefully
    if titles and (len(titles) == len(accessions)):
        treatments = [classify_treatment(t) for t in titles]
    else:
        # we cannot align titles confidently -> mark all as unknown
        treatments = ["unknown"] * len(accessions)

    meta_df = pd.DataFrame({
        "sample_id": accessions,
        "study": study,
        "treatment": treatments
    })
    return meta_df


def qc_and_plots(combined: pd.DataFrame):
    """
    Print QC to console and save two bar plots.
    """
    print("\n--- Missing values per column ---")
    print(combined.isnull().sum())

    print("\n--- Sample count per study ---")
    print(combined.groupby("study").size())

    print("\n--- Sample count per treatment ---")
    print(combined.groupby("treatment").size())

    # plots
    # 1) per study
    counts_study = combined["study"].value_counts().sort_index()
    plt.figure(figsize=(7,4))
    counts_study.plot(kind="bar")
    plt.title("Sample count per study")
    plt.xlabel("Study")
    plt.ylabel("Number of samples")
    plt.tight_layout()
    plt.savefig(PLOT_STUDY, dpi=200)
    plt.close()

    # 2) per treatment
    counts_treat = combined["treatment"].value_counts().sort_index()
    plt.figure(figsize=(7,4))
    counts_treat.plot(kind="bar")
    plt.title("Sample count per treatment")
    plt.xlabel("Treatment")
    plt.ylabel("Number of samples")
    plt.tight_layout()
    plt.savefig(PLOT_TREAT, dpi=200)
    plt.close()

    print(f"\nSaved plots:\n- {PLOT_STUDY}\n- {PLOT_TREAT}")


def main():
    ensure_dirs()

    # sanity check: show which files we found
    print("Checking files:")
    for f in DATASETS:
        print(("FOUND  " if os.path.exists(f) else "MISSING") + f"  {f}")

    all_meta = []
    for f in DATASETS:
        if not os.path.exists(f):
            continue
        try:
            meta_df = process_one_file(f)
            all_meta.append(meta_df)
        except Exception as e:
            print(f"[WARN] failed to process {f}: {e}")

    if not all_meta:
        print("No metadata collected. Are your paths correct?")
        return

    combined = pd.concat(all_meta, ignore_index=True)

    # QC + plots
    qc_and_plots(combined)

    # Save final CSV
    combined.to_csv(OUTPUT_CSV, index=False)
    print(f"\nHarmonized metadata saved to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()

