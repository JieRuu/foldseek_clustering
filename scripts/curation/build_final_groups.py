from pathlib import Path
import re
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[2]

FOLDSEEK_FILE = (
    PROJECT_DIR
    / "results"
    / "foldseek"
    / "structure_level"
    / "cluster_t0.70_structure_level.tsv"
)

LABEL_FILE = (
    PROJECT_DIR
    / "metadata"
    / "antigen_clusters_t0p70_cluster_labels.tsv"
)

SYNONYM_FILE = (
    PROJECT_DIR
    / "metadata"
    / "label_synonyms.tsv"
)

OUTPUT_DIR = PROJECT_DIR / "results" / "curation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FINAL_MAPPING_FILE = OUTPUT_DIR / "final_group_mapping_structure_level.tsv"
MERGE_SUMMARY_FILE = OUTPUT_DIR / "merge_summary_structure_level.tsv"
CLUSTER_LABEL_QC_FILE = OUTPUT_DIR / "cluster_label_qc_structure_level.tsv"

PENDING_LABEL = "MISSING_LABEL_PENDING_REVIEW"


def clean_id(x):
    x = str(x).strip()
    x = x.replace(".pdb", "")
    return x


def to_structure_id(x):
    """
    Collapse Foldseek chain-level IDs back to structure-level IDs.

    Example:
    6ehg_antigen_only_A -> 6ehg_antigen_only
    8CI1_1_antigen_only_B -> 8CI1_1_antigen_only
    """
    x = clean_id(x)
    m = re.match(r"^(.+_antigen_only)_[A-Za-z0-9]+$", x)
    if m:
        return m.group(1)
    return x


def load_foldseek_clusters(path):
    df = pd.read_csv(path, sep="\t")

    expected_cols = {"cluster_representative", "structure_id"}

    if expected_cols.issubset(set(df.columns)):
        df = df[["cluster_representative", "structure_id"]].copy()
        df = df.rename(
            columns={
                "cluster_representative": "foldseek_cluster",
                "structure_id": "pdb_id",
            }
        )
    else:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["foldseek_cluster", "pdb_id"],
        )

    df["foldseek_cluster"] = df["foldseek_cluster"].map(to_structure_id)
    df["pdb_id"] = df["pdb_id"].map(to_structure_id)

    df = (
        df.drop_duplicates()
        .sort_values(["foldseek_cluster", "pdb_id"])
        .reset_index(drop=True)
    )

    return df


def load_labels(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["original_foldseek_cluster", "pdb_id", "manual_label"],
    )

    df["original_foldseek_cluster"] = df["original_foldseek_cluster"].map(to_structure_id)
    df["pdb_id"] = df["pdb_id"].map(to_structure_id)
    df["manual_label"] = df["manual_label"].astype(str).str.strip()

    df = df[
        ["original_foldseek_cluster", "pdb_id", "manual_label"]
    ].drop_duplicates()

    return df


def load_synonyms(path):
    if not path.exists():
        return {}

    df = pd.read_csv(path, sep="\t")

    if df.empty:
        return {}

    if not {"original_label", "normalized_label"}.issubset(set(df.columns)):
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["original_label", "normalized_label"],
        )

    df["original_label"] = df["original_label"].astype(str).str.strip()
    df["normalized_label"] = df["normalized_label"].astype(str).str.strip()

    return dict(zip(df["original_label"], df["normalized_label"]))


def main():
    foldseek_df = load_foldseek_clusters(FOLDSEEK_FILE)
    label_df = load_labels(LABEL_FILE)
    synonym_map = load_synonyms(SYNONYM_FILE)

    merged = foldseek_df.merge(
        label_df[["pdb_id", "manual_label"]],
        on="pdb_id",
        how="left",
    )

    merged["label_source"] = "manual"
    merged.loc[merged["manual_label"].isna(), "label_source"] = "pending_review"

    merged["manual_label"] = merged["manual_label"].fillna(PENDING_LABEL)

    merged["normalized_label"] = merged["manual_label"].replace(synonym_map)
    merged["final_group"] = merged["normalized_label"]

    merged = (
        merged[
            [
                "foldseek_cluster",
                "pdb_id",
                "manual_label",
                "normalized_label",
                "final_group",
                "label_source",
            ]
        ]
        .drop_duplicates()
        .sort_values(["final_group", "foldseek_cluster", "pdb_id"])
        .reset_index(drop=True)
    )

    group_to_id = {
        group: f"final_group_{i + 1:03d}"
        for i, group in enumerate(sorted(merged["final_group"].unique()))
    }

    merged["final_group_id"] = merged["final_group"].map(group_to_id)

    merged = merged[
        [
            "final_group_id",
            "final_group",
            "foldseek_cluster",
            "pdb_id",
            "manual_label",
            "normalized_label",
            "label_source",
        ]
    ]

    merged.to_csv(FINAL_MAPPING_FILE, sep="\t", index=False)

    summary = (
        merged.groupby(["final_group_id", "final_group"], dropna=False)
        .agg(
            n_structures=("pdb_id", "nunique"),
            n_foldseek_clusters=("foldseek_cluster", "nunique"),
            n_pending_review=(
                "label_source",
                lambda x: (x == "pending_review").sum(),
            ),
            foldseek_clusters=(
                "foldseek_cluster",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            structure_ids=(
                "pdb_id",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            original_labels=(
                "manual_label",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
        )
        .reset_index()
        .sort_values(["n_structures", "final_group"], ascending=[False, True])
    )

    summary.to_csv(MERGE_SUMMARY_FILE, sep="\t", index=False)

    qc = (
        merged.groupby("pdb_id")
        .agg(
            n_final_groups=("final_group", "nunique"),
            final_groups=(
                "final_group",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            n_manual_labels=("manual_label", "nunique"),
            manual_labels=(
                "manual_label",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            foldseek_clusters=(
                "foldseek_cluster",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            label_sources=(
                "label_source",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
        )
        .reset_index()
        .sort_values(["n_final_groups", "pdb_id"])
    )

    qc.to_csv(CLUSTER_LABEL_QC_FILE, sep="\t", index=False)

    conflicts = qc[qc["n_final_groups"] > 1]
    pending = merged[merged["label_source"] == "pending_review"]

    # =====================================================
    # Excel output
    # =====================================================

    EXCEL_FILE = (
        OUTPUT_DIR
        / "final_group_mapping_structure_level.xlsx"
    )

    with pd.ExcelWriter(EXCEL_FILE, engine="openpyxl") as writer:

        excel_mapping = merged[
            [
                "final_group_id",
                "final_group",
                "foldseek_cluster",
                "pdb_id",
                "manual_label",
                "label_source",
            ]
        ].copy()

        excel_mapping.columns = [
            "Final Group ID",
            "Biological Group",
            "Foldseek Cluster",
            "Structure ID",
            "Manual Label",
            "Label Source",
        ]

        excel_mapping.to_excel(
            writer,
            sheet_name="Final_Group_Mapping",
            index=False,
        )

        summary.to_excel(
            writer,
            sheet_name="Merge_Summary",
            index=False,
        )

        pending.to_excel(
            writer,
            sheet_name="Pending_Review",
            index=False,
        )

    print("Input Foldseek structure-level rows:", len(foldseek_df))
    print("Input Foldseek unique structures:", foldseek_df["pdb_id"].nunique())
    print("Input manual label rows:", len(label_df))
    print("Input manual label unique structures:", label_df["pdb_id"].nunique())

    print("\nFinal mapping rows:", len(merged))
    print("Final unique structures:", merged["pdb_id"].nunique())
    print("Final groups including pending:", merged["final_group"].nunique())
    print("Pending-review structures:", pending["pdb_id"].nunique())
    print("Label conflicts:", len(conflicts))

    print("\nOutput files:")
    print(FINAL_MAPPING_FILE)
    print(MERGE_SUMMARY_FILE)
    print(CLUSTER_LABEL_QC_FILE)
    print(EXCEL_FILE)

    if len(conflicts) > 0:
        print("\nWARNING: Some structures map to multiple final groups:")
        print(conflicts.to_string(index=False))
    else:
        print("\nOK: No structure maps to multiple final groups.")

    if len(pending) > 0:
        print("\nPending-review structure IDs:")
        for pdb_id in sorted(pending["pdb_id"].unique()):
            print(pdb_id)


if __name__ == "__main__":
    main()