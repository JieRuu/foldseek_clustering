from pathlib import Path
import pandas as pd

from Bio.PDB import PDBParser, PPBuilder


PROJECT_DIR = Path(__file__).resolve().parents[2]

ANTIGEN_ONLY_DIR = PROJECT_DIR / "results" / "antigen_only_pdb"
OUTPUT_DIR = PROJECT_DIR / "results" / "qc"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "multichain_antigen_qc.tsv"

parser = PDBParser(QUIET=True)
ppb = PPBuilder()


def get_chain_sequence(chain):
    peptides = ppb.build_peptides(chain)
    return "".join(str(peptide.get_sequence()) for peptide in peptides)


def main():
    pdb_files = sorted(ANTIGEN_ONLY_DIR.glob("*.pdb"))

    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in: {ANTIGEN_ONLY_DIR}")

    rows = []

    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem

        try:
            structure = parser.get_structure(pdb_id, pdb_file)
            model = structure[0]

            chain_rows = []

            for chain in model:
                chain_seq = get_chain_sequence(chain)

                chain_rows.append({
                    "chain_id": chain.id,
                    "chain_length": len(chain_seq),
                })

            n_chains = len(chain_rows)

            rows.append({
                "pdb_id": pdb_id,
                "pdb_file": str(pdb_file),
                "n_chains": n_chains,
                "chains": ";".join(
                    f"{x['chain_id']}:{x['chain_length']}"
                    for x in chain_rows
                ),
                "status": "MULTICHAIN_CHECK" if n_chains > 1 else "OK",
            })

        except Exception as e:
            rows.append({
                "pdb_id": pdb_id,
                "pdb_file": str(pdb_file),
                "n_chains": None,
                "chains": None,
                "status": "ERROR",
                "error": str(e),
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print("Saved multichain QC to:")
    print(OUTPUT_FILE)

    print("\nStatus counts:")
    print(df["status"].value_counts(dropna=False))

    print("\nMultichain antigen-only PDBs:")
    multi = df[df["status"] == "MULTICHAIN_CHECK"]

    if multi.empty:
        print("None")
    else:
        print(
            multi[
                ["pdb_id", "n_chains", "chains"]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()