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

OUTPUT_FILE = QC_DIR / "id_mapping_qc.csv"


# =========================================================
# ID FUNCTIONS
# =========================================================

def normalize_id(x):
    """
    Normalize PDB / FASTA IDs for matching.

    Examples:
      1bzq       -> 1bzq
      1bzq_1     -> 1bzq
      7OAP_2     -> 7oap
      5O04_1_Nano_26 -> 5o04

    Note:
      IDs are normalized to lowercase.
    """
    x = str(x).strip()

    if "_Nano_" in x:
        return x.split("_")[0].lower()

    if re.search(r"_[0-9]+$", x):
        return x.rsplit("_", 1)[0].lower()

    return x.lower()


def classify_pdb_id(pdb_id):
    """
    Classify PDB folder ID pattern.

    Normal mutant IDs:
      1bzq_1, 1bzq_2, 5m2i_3

    Special Nano IDs:
      5O04_1_Nano_26

    Reference-like / exception IDs:
      xxx_0, xxx_01, xxx_02

    No suffix:
      1bzq
    """
    pdb_id = str(pdb_id).strip()

    if "_Nano_" in pdb_id:
        return "NANO_SUFFIX"

    if re.search(r"_[0-9]+$", pdb_id):
        suffix = pdb_id.rsplit("_", 1)[1]

        if suffix.startswith("0"):
            return "REFERENCE_SUFFIX_CHECK"

        return "MUTANT_SUFFIX"

    return "NO_SUFFIX"


def load_fasta_id_map(path):
    """
    Return a dictionary:
      normalized_id -> list of raw FASTA IDs

    This allows multiple FASTA records to map to the same normalized ID.
    """
    id_map = {}

    for rec in SeqIO.parse(path, "fasta"):
        raw_id = rec.id.strip()
        norm_id = normalize_id(raw_id)

        id_map.setdefault(norm_id, [])
        id_map[norm_id].append(raw_id)

    return id_map


def find_native_pdb_file(pdb_dir):
    """
    Expected structure:
      PDB_ID/
        native/
          *.pdb

    If native/ exists, use the first .pdb inside native/.
    Otherwise, search inside the PDB folder.
    """
    native_dir = pdb_dir / "native"

    if native_dir.exists() and native_dir.is_dir():
        pdb_files = sorted(native_dir.glob("*.pdb"))
        if pdb_files:
            return pdb_files[0]

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if pdb_files:
        return pdb_files[0]

    pdb_files = sorted(pdb_dir.rglob("*.pdb"))
    if pdb_files:
        return pdb_files[0]

    return None


def scan_pdb_directories(pdb_root):
    return [p for p in sorted(pdb_root.iterdir()) if p.is_dir()]


def join_ids(ids):
    if not ids:
        return None
    return ";".join(sorted(set(ids)))


# =========================================================
# MAIN
# =========================================================

def main():
    print("Config file:", CONFIG_FILE)
    print("PDB root:", PDB_ROOT)
    print("FASTA root:", FASTA_ROOT)
    print("Nanobody FASTA:", NANOBODY_FASTA)
    print("Antigen FASTA:", ANTIGEN_FASTA)
    print("Output:", OUTPUT_FILE)

    if not PDB_ROOT.exists():
        raise FileNotFoundError(f"PDB root not found: {PDB_ROOT}")

    if not NANOBODY_FASTA.exists():
        raise FileNotFoundError(f"Nanobody FASTA not found: {NANOBODY_FASTA}")

    if not ANTIGEN_FASTA.exists():
        raise FileNotFoundError(f"Antigen FASTA not found: {ANTIGEN_FASTA}")

    nanobody_id_map = load_fasta_id_map(NANOBODY_FASTA)
    antigen_id_map = load_fasta_id_map(ANTIGEN_FASTA)

    pdb_dirs = scan_pdb_directories(PDB_ROOT)

    if not pdb_dirs:
        raise ValueError(f"No PDB directories found in: {PDB_ROOT}")

    rows = []

    for pdb_dir in pdb_dirs:
        pdb_id = pdb_dir.name
        normalized_id = normalize_id(pdb_id)
        id_case_type = classify_pdb_id(pdb_id)

        native_pdb_file = find_native_pdb_file(pdb_dir)

        nanobody_matches = nanobody_id_map.get(normalized_id, [])
        antigen_matches = antigen_id_map.get(normalized_id, [])

        pdb_file_exists = native_pdb_file is not None
        nanobody_found = len(nanobody_matches) > 0
        antigen_found = len(antigen_matches) > 0

        status_flags = []

        if not pdb_file_exists:
            status_flags.append("MISSING_PDB_FILE")

        if not nanobody_found:
            status_flags.append("MISSING_NANOBODY_FASTA")

        if not antigen_found:
            status_flags.append("MISSING_ANTIGEN_FASTA")

        if id_case_type == "REFERENCE_SUFFIX_CHECK":
            status_flags.append("REFERENCE_SUFFIX_CHECK")

        if not status_flags:
            status = "OK"
        else:
            status = ";".join(status_flags)

        rows.append({
            "pdb_id": pdb_id,
            "normalized_id": normalized_id,
            "pdb_dir": str(pdb_dir),
            "native_pdb_file": str(native_pdb_file) if native_pdb_file else None,
            "pdb_file_exists": pdb_file_exists,
            "id_case_type": id_case_type,
            "nanobody_fasta_ids": join_ids(nanobody_matches),
            "antigen_fasta_ids": join_ids(antigen_matches),
            "nanobody_fasta_found": nanobody_found,
            "antigen_fasta_found": antigen_found,
            "status": status,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_FILE, index=False)

    print("\nSaved ID mapping QC to:", OUTPUT_FILE)

    print("\nStatus counts:")
    print(df["status"].value_counts(dropna=False))

    print("\nID case counts:")
    print(df["id_case_type"].value_counts(dropna=False))

    print("\nTotal PDB folders scanned:", len(df))

    print("\nProblematic rows:")
    problem_df = df[df["status"] != "OK"]
    if problem_df.empty:
        print("None")
    else:
        print(
            problem_df[
                [
                    "pdb_id",
                    "normalized_id",
                    "id_case_type",
                    "nanobody_fasta_ids",
                    "antigen_fasta_ids",
                    "status",
                ]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()