from pathlib import Path
import re
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[2]

GROUP_MAPPING_FILE = (
    PROJECT_DIR
    / "results"
    / "curation"
    / "final_group_mapping_structure_level.tsv"
)

# 如果你目前用的是 interim 檔，就改成這行：
# GROUP_MAPPING_FILE = PROJECT_DIR / "results" / "curation" / "final_group_mapping_structure_level_interim.tsv"

NANOBODY_FASTA_FILE = (
    PROJECT_DIR
    / ".."
    / "Data"
    / "Fasta"
    / "nanobodies_wt.fasta"
)

OUTPUT_DIR = PROJECT_DIR / "results" / "curation"
QC_DIR = PROJECT_DIR / "results" / "qc"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

MASTER_TABLE_OUT = OUTPUT_DIR / "nanobody_group_master_table.tsv"
QC_OUT = QC_DIR / "nanobody_mapping_qc.tsv"

MISSING_SEQUENCE = "MISSING_NANOBODY_SEQUENCE"


def clean_id(x):
    x = str(x).strip()
    x = x.replace(".pdb", "")
    return x


def structure_id_to_candidate_fasta_ids(structure_id):
    """
    Convert structure-level antigen ID to possible nanobody FASTA IDs.

    Examples:
    7Z1A_1_antigen_only -> 7Z1A_1, 7Z1A
    1RJC_antigen_only -> 1RJC, 1RJC_1
    2x89_2_antigen_only -> 2x89_2, 2x89, 2x89_1
    7OAP_2_antigen_only -> 7OAP_2, 7OAP, 7OAP_1
    """
    sid = clean_id(structure_id)

    if sid.endswith("_antigen_only"):
        base = sid.replace("_antigen_only", "")
    else:
        base = sid

    candidates = []
    candidates.append(base)

    # If ID ends with numeric suffix, also try base without suffix.
    # Example: 7OAP_02 -> 7OAP
    m = re.match(r"^(.+)_(\d+)$", base)
    if m:
        root = m.group(1)
        num = m.group(2)

        candidates.append(root)
        candidates.append(root + "_1")
        candidates.append(root + "_" + num.zfill(2))

    # If ID has no numeric suffix, also try _1.
    if not re.search(r"_\d+$", base):
        candidates.append(base + "_1")

    return list(dict.fromkeys(candidates))


def parse_fasta(path):
    records = {}
    current_id = None
    seq_chunks = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(seq_chunks)

                current_id = line[1:].split()[0].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if current_id is not None:
            records[current_id] = "".join(seq_chunks)

    return records


def main():
    group_df = pd.read_csv(GROUP_MAPPING_FILE, sep="\t")

    # Support both old and new column names
    if "pdb_id" in group_df.columns:
        structure_col = "pdb_id"
    elif "structure_id" in group_df.columns:
        structure_col = "structure_id"
    else:
        raise ValueError("Could not find structure ID column: expected 'pdb_id' or 'structure_id'.")

    if "final_group" not in group_df.columns:
        if "normalized_label" in group_df.columns:
            group_df["final_group"] = group_df["normalized_label"]
        else:
            raise ValueError("Could not find final group column: expected 'final_group' or 'normalized_label'.")

    if "final_group_id" not in group_df.columns:
        group_to_id = {
            group: f"final_group_{i + 1:03d}"
            for i, group in enumerate(sorted(group_df["final_group"].dropna().unique()))
        }
        group_df["final_group_id"] = group_df["final_group"].map(group_to_id)

    if "foldseek_cluster" not in group_df.columns:
        if "cluster_representative" in group_df.columns:
            group_df["foldseek_cluster"] = group_df["cluster_representative"]
        else:
            group_df["foldseek_cluster"] = ""

    if "label_source" not in group_df.columns:
        group_df["label_source"] = "manual"

    group_df[structure_col] = group_df[structure_col].map(clean_id)

    fasta_records = parse_fasta(NANOBODY_FASTA_FILE)

    # Case-insensitive lookup
    fasta_lookup = {
        fid.lower(): fid
        for fid in fasta_records.keys()
    }

    rows = []
    qc_rows = []

    for _, row in group_df.iterrows():
        structure_id = row[structure_col]
        candidates = structure_id_to_candidate_fasta_ids(structure_id)

        matched_fasta_id = None

        for cand in candidates:
            key = cand.lower()
            if key in fasta_lookup:
                matched_fasta_id = fasta_lookup[key]
                break

        if matched_fasta_id is None:
            nanobody_sequence = MISSING_SEQUENCE
            match_status = "missing"
        else:
            nanobody_sequence = fasta_records[matched_fasta_id]
            match_status = "matched"

        rows.append(
            {
                "final_group_id": row["final_group_id"],
                "final_group": row["final_group"],
                "foldseek_cluster": row["foldseek_cluster"],
                "structure_id": structure_id,
                "nanobody_fasta_id": matched_fasta_id if matched_fasta_id else "",
                "nanobody_sequence": nanobody_sequence,
                "label_source": row["label_source"],
                "nanobody_match_status": match_status,
                "candidate_fasta_ids": ";".join(candidates),
            }
        )

        qc_rows.append(
            {
                "structure_id": structure_id,
                "candidate_fasta_ids": ";".join(candidates),
                "matched_fasta_id": matched_fasta_id if matched_fasta_id else "",
                "match_status": match_status,
            }
        )

    master_df = pd.DataFrame(rows).drop_duplicates()
    qc_df = pd.DataFrame(qc_rows).drop_duplicates()

    master_df.to_csv(MASTER_TABLE_OUT, sep="\t", index=False)
    qc_df.to_csv(QC_OUT, sep="\t", index=False)

    print("Input structures:", group_df[structure_col].nunique())
    print("FASTA records:", len(fasta_records))
    print("Matched structures:", (qc_df["match_status"] == "matched").sum())
    print("Missing structures:", (qc_df["match_status"] == "missing").sum())

    print("\nOutput files:")
    print(MASTER_TABLE_OUT)
    print(QC_OUT)

    missing = qc_df[qc_df["match_status"] == "missing"]

    if len(missing) > 0:
        print("\nMissing nanobody sequence structures:")
        for sid in sorted(missing["structure_id"].unique()):
            print(sid)
    else:
        print("\nOK: All structures have matched nanobody sequences.")


if __name__ == "__main__":
    main()