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

OUTPUT_DIR = PROJECT_DIR / "results" / "qc"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SENSITIVITY_OUT = OUTPUT_DIR / "nanobody_threshold_sensitivity.tsv"
GROUP_SENSITIVITY_OUT = OUTPUT_DIR / "nanobody_threshold_sensitivity_by_group.tsv"

THRESHOLDS = [1, 2, 3, 5, 10]

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


def count_subclusters(group_df, threshold):
    n = len(group_df)

    if n <= 1:
        return 1

    sequences = group_df["nanobody_sequence"].tolist()
    dist_matrix = compute_distance_matrix(sequences)
    condensed_dist = squareform(dist_matrix)

    Z = linkage(condensed_dist, method="average")

    labels = fcluster(
        Z,
        t=threshold,
        criterion="distance",
    )

    return len(set(labels))


def main():
    df = pd.read_csv(MASTER_TABLE_FILE, sep="\t")

    df = df[df["nanobody_sequence"] != MISSING_SEQUENCE].copy()

    # Exclude pending antigen-label structures from production sensitivity analysis.
    if "label_source" in df.columns:
        df = df[df["label_source"] != "pending_review"].copy()
    else:
        df = df[df["final_group"] != PENDING_LABEL].copy()

    total_structures = df["structure_id"].nunique()
    total_groups = df["final_group"].nunique()

    summary_rows = []
    group_rows = []

    for threshold in THRESHOLDS:
        n_representatives = 0

        for (final_group_id, final_group), group_df in df.groupby(
            ["final_group_id", "final_group"],
            sort=True,
        ):
            n_subclusters = count_subclusters(group_df, threshold)
            n_representatives += n_subclusters

            group_rows.append(
                {
                    "threshold": threshold,
                    "final_group_id": final_group_id,
                    "final_group": final_group,
                    "n_structures": group_df["structure_id"].nunique(),
                    "n_nanobody_subclusters": n_subclusters,
                }
            )

        reduction_percent = (
            (total_structures - n_representatives)
            / total_structures
            * 100
        )

        summary_rows.append(
            {
                "threshold": threshold,
                "n_input_structures": total_structures,
                "n_biological_groups": total_groups,
                "n_representatives": n_representatives,
                "reduction_percent": round(reduction_percent, 2),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    group_df = pd.DataFrame(group_rows)

    summary_df.to_csv(SENSITIVITY_OUT, sep="\t", index=False)
    group_df.to_csv(GROUP_SENSITIVITY_OUT, sep="\t", index=False)

    print("Threshold sensitivity analysis complete.")
    print("Input structures:", total_structures)
    print("Biological groups:", total_groups)

    print("\nSummary:")
    print(summary_df.to_string(index=False))

    print("\nCoronavirus Spike sensitivity:")
    spike = group_df[group_df["final_group"] == "Coronavirus Spike"]
    if len(spike) > 0:
        print(
            spike[
                [
                    "threshold",
                    "n_structures",
                    "n_nanobody_subclusters",
                ]
            ].to_string(index=False)
        )
    else:
        print("Coronavirus Spike not found.")

    print("\nOutput files:")
    print(SENSITIVITY_OUT)
    print(GROUP_SENSITIVITY_OUT)


if __name__ == "__main__":
    main()