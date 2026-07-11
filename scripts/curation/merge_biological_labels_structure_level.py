#!/usr/bin/env python3

import re
from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

LABEL_TABLE = PROJECT_ROOT / "metadata" / "antigen_clusters_t0p70_cluster_labels.tsv"
SYNONYM_TABLE = PROJECT_ROOT / "metadata" / "label_synonyms.tsv"

OUTDIR = PROJECT_ROOT / "results" / "curation"
QC_DIR = PROJECT_ROOT / "results" / "qc"

OUTDIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)


FINAL_MAPPING_OUT = OUTDIR / "final_group_mapping_structure_level.tsv"
MERGE_SUMMARY_OUT = OUTDIR / "merge_summary_structure_level.tsv"
QC_OUT = QC_DIR / "cluster_label_qc_structure_level.tsv"


def to_structure_id(x: str) -> str:
    """
    Collapse Foldseek chain-level antigen IDs to structure-level IDs.

    Example:
    6ehg_antigen_only_A -> 6ehg_antigen_only
    8CI1_1_antigen_only_B -> 8CI1_1_antigen_only
    """
    x = str(x).strip()
    m = re.match(r"^(.+_antigen_only)_[A-Za-z0-9]+$", x)
    if m:
        return m.group(1)
    return x


def main():
    labels = pd.read_csv(LABEL_TABLE, sep="\t", header=None)

    if labels.shape[1] < 3:
        raise ValueError("Label table must have at least 3 columns: cluster, id, label")

    labels = labels.iloc[:, :3]
    labels.columns = ["foldseek_cluster", "foldseek_id", "manual_label"]

    labels["foldseek_cluster"] = labels["foldseek_cluster"].astype(str).str.strip()
    labels["foldseek_id"] = labels["foldseek_id"].astype(str).str.strip()
    labels["manual_label"] = labels["manual_label"].astype(str).str.strip()

    labels["structure_id"] = labels["foldseek_id"].map(to_structure_id)

    chain_level_rows = labels[
        labels["foldseek_id"] != labels["structure_id"]
    ].copy()

    print("Input label rows:", len(labels))
    print("Unique raw Foldseek IDs:", labels["foldseek_id"].nunique())
    print("Unique structure-level IDs:", labels["structure_id"].nunique())
    print("Chain-level rows collapsed:", len(chain_level_rows))

    synonyms = pd.read_csv(SYNONYM_TABLE, sep="\t", header=None)

    if synonyms.shape[1] < 2:
        raise ValueError("Synonym table must have at least 2 columns: original_label, normalized_label")

    synonyms = synonyms.iloc[:, :2]
    synonyms.columns = ["original_label", "normalized_label"]

    synonyms["original_label"] = synonyms["original_label"].astype(str).str.strip()
    synonyms["normalized_label"] = synonyms["normalized_label"].astype(str).str.strip()

    synonym_map = dict(zip(synonyms["original_label"], synonyms["normalized_label"]))

    labels["normalized_label"] = labels["manual_label"].map(
        lambda x: synonym_map.get(x, x)
    )

    # Collapse duplicated chain-level rows after structure normalization.
    final_mapping = (
        labels[
            [
                "foldseek_cluster",
                "structure_id",
                "manual_label",
                "normalized_label",
            ]
        ]
        .drop_duplicates()
        .sort_values(["normalized_label", "structure_id", "foldseek_cluster"])
    )

    final_mapping.to_csv(FINAL_MAPPING_OUT, sep="\t", index=False)

    merge_summary = (
        final_mapping.groupby("normalized_label")
        .agg(
            n_structures=("structure_id", "nunique"),
            n_foldseek_clusters=("foldseek_cluster", "nunique"),
            foldseek_clusters=(
                "foldseek_cluster",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            structure_ids=(
                "structure_id",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            original_labels=(
                "manual_label",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
        )
        .reset_index()
        .sort_values(["n_structures", "normalized_label"], ascending=[False, True])
    )

    merge_summary.to_csv(MERGE_SUMMARY_OUT, sep="\t", index=False)

    # QC: after collapsing chain IDs, check if one structure has conflicting labels.
    qc = (
        labels.groupby("structure_id")
        .agg(
            n_raw_foldseek_ids=("foldseek_id", "nunique"),
            raw_foldseek_ids=(
                "foldseek_id",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            n_manual_labels=("manual_label", "nunique"),
            manual_labels=(
                "manual_label",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            n_normalized_labels=("normalized_label", "nunique"),
            normalized_labels=(
                "normalized_label",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            foldseek_clusters=(
                "foldseek_cluster",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
        )
        .reset_index()
        .sort_values(["n_normalized_labels", "n_raw_foldseek_ids", "structure_id"])
    )

    qc.to_csv(QC_OUT, sep="\t", index=False)

    conflicts = qc[qc["n_normalized_labels"] > 1]

    print("Final mapping rows:", len(final_mapping))
    print("Final unique structures:", final_mapping["structure_id"].nunique())
    print("Final biological groups:", merge_summary["normalized_label"].nunique())
    print("Label conflicts after normalization:", len(conflicts))

    print("\nOutput files:")
    print(FINAL_MAPPING_OUT)
    print(MERGE_SUMMARY_OUT)
    print(QC_OUT)

    if len(conflicts) > 0:
        print("\nWARNING: Some structures have conflicting normalized labels:")
        print(conflicts.to_string(index=False))
    else:
        print("\nOK: No structure has conflicting normalized labels.")


if __name__ == "__main__":
    main()