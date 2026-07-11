from pathlib import Path
import os
import yaml
import pandas as pd

from dotenv import load_dotenv
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select, PPBuilder
from Levenshtein import distance as levenshtein_distance


# =========================================================
# PROJECT / CONFIG
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

ID_MAPPING_FILE = QC_DIR / "id_mapping_qc.csv"

NANOBODY_FASTA = FASTA_ROOT / cfg["input"]["nanobody_fasta"]

OUTPUT_DIR = RESULTS_ROOT / cfg["preprocessing"]["antigen_only_dir"]
QC_CSV = RESULTS_ROOT / cfg["preprocessing"]["qc_report"]

SIM_THRESHOLD = cfg["preprocessing"].get("nanobody_similarity_threshold", 0.80)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_ROOT.mkdir(parents=True, exist_ok=True)


# =========================================================
# FUNCTIONS
# =========================================================

def load_fasta_as_dict(path):
    seqs = {}
    for rec in SeqIO.parse(path, "fasta"):
        seqs[rec.id.strip()] = str(rec.seq)
    return seqs


def normalized_similarity(a, b):
    if not a or not b:
        return 0.0
    return 1.0 - levenshtein_distance(a, b) / max(len(a), len(b))


parser = PDBParser(QUIET=True)
ppb = PPBuilder()


def get_chain_sequence(chain):
    peptides = ppb.build_peptides(chain)
    return "".join(str(peptide.get_sequence()) for peptide in peptides)


class RemoveChainSelect(Select):
    def __init__(self, chain_to_remove):
        self.chain_to_remove = chain_to_remove

    def accept_chain(self, chain):
        return chain.id != self.chain_to_remove


def parse_single_fasta_id(fasta_ids):
    """
    id_mapping_qc.csv stores matched FASTA IDs as strings.

    Usually this is a single ID.
    If multiple IDs are present, they are separated by semicolon.
    For extraction, use the first one and record the full string in QC.
    """
    if pd.isna(fasta_ids):
        return None

    fasta_ids = str(fasta_ids).strip()

    if not fasta_ids or fasta_ids.lower() == "none":
        return None

    return fasta_ids.split(";")[0]


# =========================================================
# MAIN
# =========================================================

def main():
    print("Config file:", CONFIG_FILE)
    print("Results root:", RESULTS_ROOT)
    print("PDB root:", PDB_ROOT)
    print("FASTA root:", FASTA_ROOT)
    print("ID mapping file:", ID_MAPPING_FILE)
    print("Nanobody FASTA:", NANOBODY_FASTA)
    print("Output dir:", OUTPUT_DIR)
    print("Similarity threshold:", SIM_THRESHOLD)

    if not ID_MAPPING_FILE.exists():
        raise FileNotFoundError(
            f"ID mapping file not found: {ID_MAPPING_FILE}\n"
            "Run scripts/preprocessing/check_id_mapping.py first."
        )

    if not NANOBODY_FASTA.exists():
        raise FileNotFoundError(f"Nanobody FASTA not found: {NANOBODY_FASTA}")

    nanobody_seqs = load_fasta_as_dict(NANOBODY_FASTA)

    mapping_df = pd.read_csv(ID_MAPPING_FILE)

    required_cols = [
        "pdb_id",
        "normalized_id",
        "native_pdb_file",
        "nanobody_fasta_ids",
        "status",
    ]

    missing_cols = [c for c in required_cols if c not in mapping_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in ID mapping file: {missing_cols}")

    matched_df = mapping_df[mapping_df["status"] == "OK"].copy()

    if matched_df.empty:
        raise ValueError("No rows with status == OK found in id_mapping_qc.csv")

    print(f"\nMatched complexes to process: {len(matched_df)}")

    rows = []

    for _, row in matched_df.iterrows():
        pdb_id = row["pdb_id"]
        normalized_id = row["normalized_id"]

        pdb_path = Path(row["native_pdb_file"])
        fasta_ids_raw = row["nanobody_fasta_ids"]
        fasta_key = parse_single_fasta_id(fasta_ids_raw)

        if fasta_key is None:
            rows.append({
                "pdb_id": pdb_id,
                "normalized_id": normalized_id,
                "pdb_file": str(pdb_path),
                "fasta_ids_raw": fasta_ids_raw,
                "fasta_key_used": None,
                "nanobody_chain": None,
                "nanobody_similarity": None,
                "status": "MISSING_NANOBODY_FASTA_KEY_IN_MAPPING",
            })
            continue

        if fasta_key not in nanobody_seqs:
            rows.append({
                "pdb_id": pdb_id,
                "normalized_id": normalized_id,
                "pdb_file": str(pdb_path),
                "fasta_ids_raw": fasta_ids_raw,
                "fasta_key_used": fasta_key,
                "nanobody_chain": None,
                "nanobody_similarity": None,
                "status": "FASTA_KEY_NOT_FOUND_IN_NANOBODY_FASTA",
            })
            continue

        if not pdb_path.exists():
            rows.append({
                "pdb_id": pdb_id,
                "normalized_id": normalized_id,
                "pdb_file": str(pdb_path),
                "fasta_ids_raw": fasta_ids_raw,
                "fasta_key_used": fasta_key,
                "nanobody_chain": None,
                "nanobody_similarity": None,
                "status": "MISSING_PDB_FILE",
            })
            continue

        nb_seq = nanobody_seqs[fasta_key]

        try:
            structure = parser.get_structure(pdb_id, pdb_path)
            model = structure[0]

            chain_results = []

            for chain in model:
                chain_seq = get_chain_sequence(chain)
                sim = normalized_similarity(chain_seq, nb_seq)

                chain_results.append({
                    "chain_id": chain.id,
                    "chain_length": len(chain_seq),
                    "similarity_to_nanobody": sim,
                })

            if not chain_results:
                rows.append({
                    "pdb_id": pdb_id,
                    "normalized_id": normalized_id,
                    "pdb_file": str(pdb_path),
                    "fasta_ids_raw": fasta_ids_raw,
                    "fasta_key_used": fasta_key,
                    "nanobody_chain": None,
                    "nanobody_similarity": None,
                    "status": "NO_CHAINS_FOUND",
                })
                continue

            chain_results = sorted(
                chain_results,
                key=lambda x: x["similarity_to_nanobody"],
                reverse=True,
            )

            best = chain_results[0]
            nanobody_chain = best["chain_id"]
            nanobody_similarity = best["similarity_to_nanobody"]

            status = (
                "OK"
                if nanobody_similarity >= SIM_THRESHOLD
                else "LOW_SIMILARITY_CHECK"
            )

            output_file = OUTPUT_DIR / f"{pdb_id}_antigen_only.pdb"

            io = PDBIO()
            io.set_structure(structure)
            io.save(str(output_file), RemoveChainSelect(nanobody_chain))

            rows.append({
                "pdb_id": pdb_id,
                "normalized_id": normalized_id,
                "pdb_file": str(pdb_path),
                "fasta_ids_raw": fasta_ids_raw,
                "fasta_key_used": fasta_key,
                "nanobody_chain": nanobody_chain,
                "nanobody_similarity": nanobody_similarity,
                "nanobody_chain_length": best["chain_length"],
                "all_chains": ";".join(
                    f"{x['chain_id']}:{x['chain_length']}:{x['similarity_to_nanobody']:.3f}"
                    for x in chain_results
                ),
                "output_file": str(output_file),
                "status": status,
            })

        except Exception as e:
            rows.append({
                "pdb_id": pdb_id,
                "normalized_id": normalized_id,
                "pdb_file": str(pdb_path),
                "fasta_ids_raw": fasta_ids_raw,
                "fasta_key_used": fasta_key,
                "nanobody_chain": None,
                "nanobody_similarity": None,
                "status": "ERROR",
                "error": str(e),
            })

    qc_df = pd.DataFrame(rows)
    qc_df.to_csv(QC_CSV, index=False)

    print("\nSaved antigen-only PDBs to:", OUTPUT_DIR)
    print("Saved QC report to:", QC_CSV)

    print("\nStatus counts:")
    print(qc_df["status"].value_counts(dropna=False))

    print("\nTotal processed rows:", len(qc_df))


if __name__ == "__main__":
    main()