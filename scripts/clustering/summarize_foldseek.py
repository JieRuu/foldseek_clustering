from pathlib import Path
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[2]
FOLDSEEK_DIR = PROJECT_DIR / "results" / "foldseek"

THRESHOLDS = ["0.70", "0.75", "0.80", "0.85", "0.90"]

OUTPUT_FILE = FOLDSEEK_DIR / "threshold_cluster_summary.tsv"


def summarize_cluster_file(path, threshold):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["cluster_rep", "member"],
    )

    cluster_sizes = (
        df.groupby("cluster_rep")
        .size()
        .reset_index(name="cluster_size")
        .sort_values("cluster_size", ascending=False)
    )

    n_clusters = len(cluster_sizes)
    largest_cluster = cluster_sizes["cluster_size"].max()
    n_singletons = (cluster_sizes["cluster_size"] == 1).sum()
    n_clusters_ge_5 = (cluster_sizes["cluster_size"] >= 5).sum()
    n_clusters_ge_6 = (cluster_sizes["cluster_size"] >= 6).sum()

    cluster_size_file = FOLDSEEK_DIR / f"cluster_size_t{threshold}.tsv"
    cluster_sizes.to_csv(cluster_size_file, sep="\t", index=False)

    return {
        "threshold": threshold,
        "n_members": df["member"].nunique(),
        "n_clusters": n_clusters,
        "largest_cluster": largest_cluster,
        "n_clusters_ge_5": n_clusters_ge_5,
        "n_clusters_ge_6": n_clusters_ge_6,
        "n_singletons": n_singletons,
        "cluster_size_file": str(cluster_size_file),
    }


def main():
    rows = []

    for threshold in THRESHOLDS:
        cluster_file = FOLDSEEK_DIR / f"cluster_t{threshold}.tsv"

        if not cluster_file.exists():
            print(f"Skipping missing file: {cluster_file}")
            continue

        print(f"Summarizing: {cluster_file}")
        rows.append(summarize_cluster_file(cluster_file, threshold))

    if not rows:
        raise FileNotFoundError("No Foldseek cluster TSV files found.")

    summary = pd.DataFrame(rows)
    summary.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print("\nSaved summary to:")
    print(OUTPUT_FILE)

    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()