from pathlib import Path
import numpy as np
import pandas as pd
import Levenshtein
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


PROJECT_DIR = Path(__file__).resolve().parents[2]

MASTER_TABLE_FILE = (
    PROJECT_DIR
    / "results"
    / "curation"
    / "nanobody_group_master_table.tsv"
)

OUTPUT_DIR = PROJECT_DIR / "results" / "benchmark"
QC_DIR = PROJECT_DIR / "results" / "qc"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

# Original / detailed outputs
REPRESENTATIVES_OUT = OUTPUT_DIR / "nanobody_benchmark_representatives.tsv"
SUBCLUSTERS_OUT = OUTPUT_DIR / "nanobody_hac_subclusters.tsv"
QC_OUT = QC_DIR / "nanobody_hac_qc.tsv"

# Clean benchmark outputs
BENCHMARK_REPRESENTATIVES_CSV = OUTPUT_DIR / "benchmark_representatives.csv"
BENCHMARK_CLUSTER_MEMBERS_CSV = OUTPUT_DIR / "benchmark_cluster_members.csv"
BENCHMARK_SUMMARY_CSV = OUTPUT_DIR / "benchmark_summary.csv"
BENCHMARK_REPRESENTATIVES_TXT = OUTPUT_DIR / "benchmark_representatives.txt"

MUTATION_THRESHOLD = 5
PENDING_LABEL = "MISSING_LABEL_PENDING_REVIEW"
MISSING_SEQUENCE = "MISSING_NANOBODY_SEQUENCE"


def compute_distance_matrix(sequences):
    n = len(sequences)
    dist_matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            dist = Levenshtein.distance(sequences[i], sequences[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    return dist_matrix


def select_medoid(indices, dist_matrix):
    """
    Select the real input data point with the smallest average Levenshtein
    distance to all other members in the same subcluster.

    This is a medoid, not a centroid.
    Therefore the representative is always a real structure / nanobody
    from the original dataset.
    """
    if len(indices) == 1:
        return indices[0], 0.0

    sub_matrix = dist_matrix[np.ix_(indices, indices)]
    avg_distances = sub_matrix.mean(axis=1)

    local_medoid_idx = int(np.argmin(avg_distances))
    medoid_idx = indices[local_medoid_idx]
    medoid_avg_distance = float(avg_distances[local_medoid_idx])

    return medoid_idx, medoid_avg_distance


def process_group(group_df):
    group_df = group_df.reset_index(drop=True)
    n = len(group_df)

    # -----------------------------------------------------
    # Workflow checkpoint:
    #
    # The original workflow document contains:
    #
    # if n_seqs <= 2:
    #     final_benchmark_set.append(cluster_ids[0])
    #     continue
    #
    # That checkpoint was intended as a small-group safeguard.
    # In this production implementation, we keep the checkpoint
    # as QC metadata, but we still run HAC/medoid selection for
    # n = 2 because two nanobodies can be highly divergent.
    #
    # This prevents incorrectly collapsing two biologically distinct
    # nanobodies into one representative.
    # -----------------------------------------------------

    small_group_checkpoint = n <= 2

    if n == 1:
        row = group_df.iloc[0]
        subcluster_id = f"{row['final_group_id']}_nb_subcluster_001"

        representatives = [
            {
                "final_group_id": row["final_group_id"],
                "final_group": row["final_group"],
                "nanobody_subcluster_id": subcluster_id,
                "representative_structure_id": row["structure_id"],
                "representative_nanobody_fasta_id": row["nanobody_fasta_id"],
                "representative_sequence": row["nanobody_sequence"],
                "subcluster_size": 1,
                "medoid_avg_distance": 0.0,
                "representative_is_real_data_point": True,
                "selection_method": "medoid",
                "selection_reason": "single_member_group",
                "small_group_checkpoint_n_le_2": small_group_checkpoint,
            }
        ]

        subcluster_rows = [
            {
                "final_group_id": row["final_group_id"],
                "final_group": row["final_group"],
                "nanobody_subcluster_id": subcluster_id,
                "structure_id": row["structure_id"],
                "nanobody_fasta_id": row["nanobody_fasta_id"],
                "nanobody_sequence": row["nanobody_sequence"],
                "is_medoid": True,
                "cluster_size": 1,
            }
        ]

        qc = {
            "final_group_id": row["final_group_id"],
            "final_group": row["final_group"],
            "n_structures": 1,
            "n_nanobody_subclusters": 1,
            "max_pairwise_distance": 0.0,
            "mean_pairwise_distance": 0.0,
            "all_representatives_are_real_data_points": True,
            "small_group_checkpoint_n_le_2": small_group_checkpoint,
            "status": "single_member_group",
        }

        return representatives, subcluster_rows, qc

    sequences = group_df["nanobody_sequence"].tolist()
    dist_matrix = compute_distance_matrix(sequences)

    condensed_dist = squareform(dist_matrix)

    # Same HAC logic as the workflow document:
    # average-linkage HAC on precomputed Levenshtein distances.
    Z = linkage(condensed_dist, method="average")

    subcluster_labels = fcluster(
        Z,
        t=MUTATION_THRESHOLD,
        criterion="distance",
    )

    group_df["nanobody_subcluster_numeric"] = subcluster_labels

    representatives = []
    subcluster_rows = []

    for sub_idx, sub_label in enumerate(sorted(set(subcluster_labels)), start=1):
        indices = [
            i
            for i, label in enumerate(subcluster_labels)
            if label == sub_label
        ]

        medoid_idx, medoid_avg_distance = select_medoid(indices, dist_matrix)
        medoid_row = group_df.iloc[medoid_idx]

        subcluster_id = (
            f"{medoid_row['final_group_id']}"
            f"_nb_subcluster_{sub_idx:03d}"
        )

        representative_is_real = (
            medoid_row["structure_id"] in set(group_df["structure_id"])
            and medoid_row["nanobody_fasta_id"] in set(group_df["nanobody_fasta_id"])
        )

        representatives.append(
            {
                "final_group_id": medoid_row["final_group_id"],
                "final_group": medoid_row["final_group"],
                "nanobody_subcluster_id": subcluster_id,
                "representative_structure_id": medoid_row["structure_id"],
                "representative_nanobody_fasta_id": medoid_row["nanobody_fasta_id"],
                "representative_sequence": medoid_row["nanobody_sequence"],
                "subcluster_size": len(indices),
                "medoid_avg_distance": medoid_avg_distance,
                "representative_is_real_data_point": representative_is_real,
                "selection_method": "medoid",
                "selection_reason": "min_average_levenshtein_distance_within_subcluster",
                "small_group_checkpoint_n_le_2": small_group_checkpoint,
            }
        )

        for idx in indices:
            row = group_df.iloc[idx]

            subcluster_rows.append(
                {
                    "final_group_id": row["final_group_id"],
                    "final_group": row["final_group"],
                    "nanobody_subcluster_id": subcluster_id,
                    "structure_id": row["structure_id"],
                    "nanobody_fasta_id": row["nanobody_fasta_id"],
                    "nanobody_sequence": row["nanobody_sequence"],
                    "is_medoid": idx == medoid_idx,
                    "cluster_size": len(indices),
                }
            )

    upper = dist_matrix[np.triu_indices(n, k=1)]

    qc = {
        "final_group_id": group_df.iloc[0]["final_group_id"],
        "final_group": group_df.iloc[0]["final_group"],
        "n_structures": n,
        "n_nanobody_subclusters": len(set(subcluster_labels)),
        "max_pairwise_distance": float(upper.max()) if len(upper) else 0.0,
        "mean_pairwise_distance": float(upper.mean()) if len(upper) else 0.0,
        "all_representatives_are_real_data_points": all(
            r["representative_is_real_data_point"] for r in representatives
        ),
        "small_group_checkpoint_n_le_2": small_group_checkpoint,
        "status": "clustered",
    }

    return representatives, subcluster_rows, qc


def main():
    df = pd.read_csv(MASTER_TABLE_FILE, sep="\t")

    required_cols = {
        "final_group_id",
        "final_group",
        "structure_id",
        "nanobody_fasta_id",
        "nanobody_sequence",
    }

    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    df = df[df["nanobody_sequence"] != MISSING_SEQUENCE].copy()

    # Production benchmark:
    # exclude structures whose antigen labels are still pending.
    if "label_source" in df.columns:
        df = df[df["label_source"] != "pending_review"].copy()
    else:
        df = df[df["final_group"] != PENDING_LABEL].copy()

    representatives_all = []
    subclusters_all = []
    qc_all = []

    for _, group_df in df.groupby(["final_group_id", "final_group"], sort=True):
        representatives, subclusters, qc = process_group(group_df)
        representatives_all.extend(representatives)
        subclusters_all.extend(subclusters)
        qc_all.append(qc)

    representatives_df = pd.DataFrame(representatives_all)
    subclusters_df = pd.DataFrame(subclusters_all)
    qc_df = pd.DataFrame(qc_all)

    representatives_df.to_csv(REPRESENTATIVES_OUT, sep="\t", index=False)
    subclusters_df.to_csv(SUBCLUSTERS_OUT, sep="\t", index=False)
    qc_df.to_csv(QC_OUT, sep="\t", index=False)

    benchmark_representatives = representatives_df[
        [
            "final_group_id",
            "final_group",
            "nanobody_subcluster_id",
            "representative_structure_id",
            "representative_nanobody_fasta_id",
            "subcluster_size",
            "medoid_avg_distance",
            "representative_is_real_data_point",
            "selection_method",
            "small_group_checkpoint_n_le_2",
        ]
    ].copy()

    benchmark_representatives.to_csv(
        BENCHMARK_REPRESENTATIVES_CSV,
        index=False,
    )

    benchmark_cluster_members = subclusters_df[
        [
            "final_group_id",
            "final_group",
            "nanobody_subcluster_id",
            "structure_id",
            "nanobody_fasta_id",
            "is_medoid",
            "cluster_size",
        ]
    ].copy()

    benchmark_cluster_members.to_csv(
        BENCHMARK_CLUSTER_MEMBERS_CSV,
        index=False,
    )

    benchmark_summary = (
        subclusters_df.groupby(["final_group_id", "final_group"])
        .agg(
            n_structures=("structure_id", "nunique"),
            n_nanobody_subclusters=("nanobody_subcluster_id", "nunique"),
            n_representatives=("is_medoid", "sum"),
        )
        .reset_index()
        .sort_values(
            ["n_structures", "n_nanobody_subclusters", "final_group"],
            ascending=[False, False, True],
        )
    )

    benchmark_summary.to_csv(
        BENCHMARK_SUMMARY_CSV,
        index=False,
    )

    with open(BENCHMARK_REPRESENTATIVES_TXT, "w") as f:
        for sid in representatives_df["representative_structure_id"]:
            f.write(f"{sid}\n")

    n_input_structures = df["structure_id"].nunique()
    n_input_groups = df["final_group"].nunique()
    n_representatives = len(representatives_df)
    n_subclusters = subclusters_df["nanobody_subcluster_id"].nunique()

    all_real = representatives_df[
        "representative_is_real_data_point"
    ].all()

    small_group_qc = qc_df[qc_df["small_group_checkpoint_n_le_2"] == True]

    print("Mutation threshold:", MUTATION_THRESHOLD)
    print("Input structures with nanobody sequences:", n_input_structures)
    print("Input biological groups:", n_input_groups)
    print("Selected representatives:", n_representatives)
    print("Nanobody subclusters:", n_subclusters)
    print("All representatives are real data points:", all_real)

    print("\nSmall-group checkpoint:")
    print("Groups with <=2 structures:", len(small_group_qc))
    if len(small_group_qc) > 0:
        print(
            small_group_qc[
                [
                    "final_group_id",
                    "final_group",
                    "n_structures",
                    "n_nanobody_subclusters",
                    "max_pairwise_distance",
                    "status",
                ]
            ].to_string(index=False)
        )

    print("\nOutput files:")
    print(REPRESENTATIVES_OUT)
    print(SUBCLUSTERS_OUT)
    print(QC_OUT)
    print(BENCHMARK_REPRESENTATIVES_CSV)
    print(BENCHMARK_CLUSTER_MEMBERS_CSV)
    print(BENCHMARK_SUMMARY_CSV)
    print(BENCHMARK_REPRESENTATIVES_TXT)

    print("\nTop groups by number of nanobody subclusters:")
    print(
        benchmark_summary.sort_values(
            ["n_nanobody_subclusters", "n_structures"],
            ascending=[False, False],
        )
        .head(20)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()