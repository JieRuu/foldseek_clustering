#!/usr/bin/env python3
"""
Sweep nanobody representative selection across mutation thresholds 5..50 (step 5).

For each threshold, runs the same HAC + medoid logic as
select_nanobody_representatives.py and writes results into a per-threshold
subfolder, so nothing gets overwritten. Also writes a combined summary
across all thresholds.

Usage:
    python scripts/curation/sweep_thresholds.py

Reads:  results/curation/nanobody_group_master_table.tsv   (the 32-group / supercluster version)
Writes: results/benchmark/threshold_sweep/t{THRESHOLD}/...
        results/benchmark/threshold_sweep/sweep_overview.tsv
"""

from pathlib import Path
import numpy as np
import pandas as pd
import Levenshtein
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


PROJECT_DIR = Path(__file__).resolve().parents[2]

MASTER_TABLE_FILE = PROJECT_DIR / "results" / "curation" / "nanobody_group_master_table.tsv"
SWEEP_DIR = PROJECT_DIR / "results" / "benchmark" / "threshold_sweep"
SWEEP_DIR.mkdir(parents=True, exist_ok=True)

THRESHOLDS = list(range(5, 51, 5))   # 5,10,15,...,50

PENDING_LABEL = "MISSING_LABEL_PENDING_REVIEW"
MISSING_SEQUENCE = "MISSING_NANOBODY_SEQUENCE"


def compute_distance_matrix(sequences):
    n = len(sequences)
    dm = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = Levenshtein.distance(sequences[i], sequences[j])
            dm[i, j] = d
            dm[j, i] = d
    return dm


def select_medoid(indices, dist_matrix):
    if len(indices) == 1:
        return indices[0], 0.0
    sub = dist_matrix[np.ix_(indices, indices)]
    avg = sub.mean(axis=1)
    local = int(np.argmin(avg))
    return indices[local], float(avg[local])


def process_group(group_df, threshold):
    group_df = group_df.reset_index(drop=True)
    n = len(group_df)

    if n == 1:
        row = group_df.iloc[0]
        sid = f"{row['final_group_id']}_nb_subcluster_001"
        rep = [{
            "final_group_id": row["final_group_id"], "final_group": row["final_group"],
            "nanobody_subcluster_id": sid,
            "representative_structure_id": row["structure_id"],
            "representative_nanobody_fasta_id": row["nanobody_fasta_id"],
            "subcluster_size": 1, "medoid_avg_distance": 0.0,
        }]
        mem = [{
            "final_group_id": row["final_group_id"], "final_group": row["final_group"],
            "nanobody_subcluster_id": sid, "structure_id": row["structure_id"],
            "nanobody_fasta_id": row["nanobody_fasta_id"], "is_medoid": True, "cluster_size": 1,
        }]
        return rep, mem

    seqs = group_df["nanobody_sequence"].tolist()
    dm = compute_distance_matrix(seqs)
    Z = linkage(squareform(dm), method="average")
    labels = fcluster(Z, t=threshold, criterion="distance")

    reps, mems = [], []
    for sub_idx, lab in enumerate(sorted(set(labels)), start=1):
        idxs = [i for i, l in enumerate(labels) if l == lab]
        m_idx, m_dist = select_medoid(idxs, dm)
        mrow = group_df.iloc[m_idx]
        sid = f"{mrow['final_group_id']}_nb_subcluster_{sub_idx:03d}"
        reps.append({
            "final_group_id": mrow["final_group_id"], "final_group": mrow["final_group"],
            "nanobody_subcluster_id": sid,
            "representative_structure_id": mrow["structure_id"],
            "representative_nanobody_fasta_id": mrow["nanobody_fasta_id"],
            "subcluster_size": len(idxs), "medoid_avg_distance": m_dist,
        })
        for i in idxs:
            r = group_df.iloc[i]
            mems.append({
                "final_group_id": r["final_group_id"], "final_group": r["final_group"],
                "nanobody_subcluster_id": sid, "structure_id": r["structure_id"],
                "nanobody_fasta_id": r["nanobody_fasta_id"],
                "is_medoid": (i == m_idx), "cluster_size": len(idxs),
            })
    return reps, mems


def main():
    df = pd.read_csv(MASTER_TABLE_FILE, sep="\t")
    df = df[df["nanobody_sequence"] != MISSING_SEQUENCE].copy()
    if "label_source" in df.columns:
        df = df[df["label_source"] != "pending_review"].copy()
    else:
        df = df[df["final_group"] != PENDING_LABEL].copy()

    total_structures = df["structure_id"].nunique()
    overview_rows = []

    for thr in THRESHOLDS:
        reps_all, mems_all = [], []
        for _, g in df.groupby(["final_group_id", "final_group"], sort=True):
            reps, mems = process_group(g, thr)
            reps_all.extend(reps)
            mems_all.extend(mems)

        reps_df = pd.DataFrame(reps_all)
        mems_df = pd.DataFrame(mems_all)

        out_dir = SWEEP_DIR / f"t{thr}"
        out_dir.mkdir(parents=True, exist_ok=True)
        reps_df.to_csv(out_dir / f"representatives_t{thr}.tsv", sep="\t", index=False)
        mems_df.to_csv(out_dir / f"cluster_members_t{thr}.tsv", sep="\t", index=False)

        n_reps = len(reps_df)
        overview_rows.append({
            "threshold": thr,
            "n_input_structures": total_structures,
            "n_representatives": n_reps,
            "reduction_percent": round((total_structures - n_reps) / total_structures * 100, 2),
        })
        print(f"threshold {thr:>2}: {n_reps} representatives  "
              f"({(total_structures - n_reps) / total_structures * 100:.1f}% reduction)  -> {out_dir}")

    overview = pd.DataFrame(overview_rows)
    overview.to_csv(SWEEP_DIR / "sweep_overview.tsv", sep="\t", index=False)

    print("\nOverview:")
    print(overview.to_string(index=False))
    print("\nAll outputs under:", SWEEP_DIR)


if __name__ == "__main__":
    main()