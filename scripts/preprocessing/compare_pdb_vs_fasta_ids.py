from pathlib import Path
import os
import re
import yaml
import pandas as pd

from dotenv import load_dotenv
from Bio import SeqIO


# =========================================================
# CONFIG
# =========================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]

load_dotenv(PROJECT_DIR / ".env")

CONFIG_FILE = PROJECT_DIR / "config" / "config.yaml"

with open(CONFIG_FILE, "r") as f:
    cfg = yaml.safe_load(f)


def require_env(name):
    value = os.getenv(name)
    if not value:
        raise EnvironmentError(f"Missing environment variable: {name}")
    return Path(value).expanduser().resolve()


PDB_ROOT = require_env("PDB_ROOT")
FASTA_ROOT = require_env("FASTA_ROOT")

RESULTS_ROOT = (PROJECT_DIR / cfg["paths"]["results_root"]).resolve()
QC_DIR = RESULTS_ROOT / "qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

NANOBODY_FASTA = FASTA_ROOT / cfg["input"]["nanobody_fasta"]
ANTIGEN_FASTA = FASTA_ROOT / cfg["input"]["antigen_fasta"]

OUTPUT_FILE = QC_DIR / "pdb_vs_fasta_id_comparison.csv"


# =========================================================
# FUNCTIONS
# =========================================================

def normalize_id(x):
    """
    Normalize IDs for comparison.

    Examples:
      1bzq_1 -> 1bzq
      7OAP_2 -> 7oap
      5O04_1_Nano_26 -> 5o04
    """
    x = str(x).strip()

    # Special Nano naming pattern:
    # 5O04_1_Nano_26 -> 5O04
    if "_Nano_" in x:
        return x.split("_")[0].lower()

    # Normal mutant suffix:
    # 1bzq_1 -> 1bzq
    # 7OAP_2 -> 7oap
    if re.search(r"_[0-9]+$", x):
        return x.rsplit("_", 1)[0].lower()

    return x.lower()


def load_fasta_ids(path):
    ids = []
    for rec in SeqIO.parse(path, "fasta"):
        raw_id = rec.id.strip()
        ids.append({
            "raw_id": raw_id,
            "normalized_id": normalize_id(raw_id),
        })
    return pd.DataFrame(ids)


def load_pdb_folder_ids(pdb_root):
    rows = []
    for p in sorted(pdb_root.iterdir()):
        if p.is_dir():
            rows.append({
                "raw_id": p.name,
                "normalized_id": normalize_id(p.name),
                "pdb_dir": str(p),
            })
    return pd.DataFrame(rows)


# =========================================================
# MAIN
# =========================================================

def main():
    print("PDB root:", PDB_ROOT)
    print("Nanobody FASTA:", NANOBODY_FASTA)
    print("Antigen FASTA:", ANTIGEN_FASTA)
    print("Output:", OUTPUT_FILE)

    pdb_df = load_pdb_folder_ids(PDB_ROOT)
    nano_df = load_fasta_ids(NANOBODY_FASTA)
    antigen_df = load_fasta_ids(ANTIGEN_FASTA)

    pdb_ids = set(pdb_df["normalized_id"])
    nano_ids = set(nano_df["normalized_id"])
    antigen_ids = set(antigen_df["normalized_id"])

    all_ids = sorted(pdb_ids | nano_ids | antigen_ids)

    rows = []

    for id_ in all_ids:
        rows.append({
            "normalized_id": id_,
            "in_pdb_folders": id_ in pdb_ids,
            "in_nanobody_fasta": id_ in nano_ids,
            "in_antigen_fasta": id_ in antigen_ids,
            "pdb_raw_ids": ";".join(sorted(pdb_df.loc[pdb_df["normalized_id"] == id_, "raw_id"].unique())),
            "nanobody_raw_ids": ";".join(sorted(nano_df.loc[nano_df["normalized_id"] == id_, "raw_id"].unique())),
            "antigen_raw_ids": ";".join(sorted(antigen_df.loc[antigen_df["normalized_id"] == id_, "raw_id"].unique())),
        })

    out = pd.DataFrame(rows)

    def classify(row):
        if row["in_pdb_folders"] and row["in_nanobody_fasta"] and row["in_antigen_fasta"]:
            return "MATCHED_ALL"
        if row["in_pdb_folders"] and not row["in_nanobody_fasta"] and not row["in_antigen_fasta"]:
            return "PDB_ONLY"
        if not row["in_pdb_folders"] and row["in_nanobody_fasta"] and row["in_antigen_fasta"]:
            return "FASTA_ONLY"
        if row["in_pdb_folders"] and row["in_nanobody_fasta"] and not row["in_antigen_fasta"]:
            return "MISSING_ANTIGEN_FASTA"
        if row["in_pdb_folders"] and not row["in_nanobody_fasta"] and row["in_antigen_fasta"]:
            return "MISSING_NANOBODY_FASTA"
        return "OTHER_MISMATCH"

    out["status"] = out.apply(classify, axis=1)

    out.to_csv(OUTPUT_FILE, index=False)

    print("\nSummary:")
    print(out["status"].value_counts(dropna=False))

    print("\nCounts:")
    print("PDB normalized IDs:", len(pdb_ids))
    print("Nanobody FASTA normalized IDs:", len(nano_ids))
    print("Antigen FASTA normalized IDs:", len(antigen_ids))

    print("\nPDB-only IDs:")
    print(out.loc[out["status"] == "PDB_ONLY", ["normalized_id", "pdb_raw_ids"]].to_string(index=False))


if __name__ == "__main__":
    main()