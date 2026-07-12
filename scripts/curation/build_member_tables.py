#!/usr/bin/env python3
"""
Build the full nanobody member table for each mutation threshold (5..50).

For every threshold, joins the per-threshold cluster membership
(results/benchmark/threshold_sweep/t{N}/cluster_members_t{N}.tsv)
with the master table's nanobody sequences, producing one row per nanobody.

Output columns (per meeting notes: PDB_ID / Nanobody_seq / antigen label / nanobody HAC label):
    PDB_ID                       - antigen structure ID
    nanobody_id                  - nanobody FASTA ID
    nanobody_seq                 - nanobody sequence
    antigen_clustering_label     - antigen group / supercluster
    nanobody_HAC_label           - nanobody HAC subcluster ID
    is_medoid                    - True if this nanobody is the subcluster medoid
    cluster_size                 - number of members in the subcluster

Usage:
    python scripts/curation/build_member_tables.py

Writes: results/benchmark/threshold_sweep/t{N}/member_table_t{N}.tsv   (one per threshold)
        results/benchmark/threshold_sweep/member_tables_all_thresholds.tsv  (combined, with a `threshold` column)
        results/benchmark/threshold_sweep/member_tables_all_thresholds.xlsx (one sheet per threshold: t5, t10, ... t50)
        results/benchmark/threshold_sweep/summary_all_thresholds.xlsx (Overview sheet + per-threshold group summary sheets)
"""

from pathlib import Path
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[2]
MASTER_TABLE = PROJECT_DIR / "results" / "curation" / "nanobody_group_master_table.tsv"
SWEEP_DIR = PROJECT_DIR / "results" / "benchmark" / "threshold_sweep"

THRESHOLDS = list(range(5, 51, 5))


def main():
    master = pd.read_csv(MASTER_TABLE, sep="\t")

    # sequence lookup by structure_id (each structure has one nanobody)
    seq_lookup = master.set_index("structure_id")["nanobody_sequence"].to_dict()

    combined = []
    sheets = {}          # threshold -> member dataframe (for the member xlsx)
    summary_sheets = {}  # threshold -> per-group summary dataframe
    overview_rows = []   # one row per threshold for the Overview sheet

    for thr in THRESHOLDS:
        members_file = SWEEP_DIR / f"t{thr}" / f"cluster_members_t{thr}.tsv"
        if not members_file.exists():
            print(f"[skip] missing: {members_file}")
            continue

        m = pd.read_csv(members_file, sep="\t")

        # attach the nanobody sequence (members file may already carry fasta_id but not sequence)
        m["nanobody_seq"] = m["structure_id"].map(seq_lookup)

        out = pd.DataFrame({
            "PDB_ID": m["structure_id"],
            "nanobody_id": m["nanobody_fasta_id"],
            "nanobody_seq": m["nanobody_seq"],
            "antigen_clustering_label": m["final_group"],
            "nanobody_HAC_label": m["nanobody_subcluster_id"],
            "is_medoid": m["is_medoid"],
            "cluster_size": m["cluster_size"],
        }).sort_values(["antigen_clustering_label", "nanobody_HAC_label", "PDB_ID"])

        out_file = SWEEP_DIR / f"t{thr}" / f"member_table_t{thr}.tsv"
        out.to_csv(out_file, sep="\t", index=False)

        sheets[thr] = out

        c = out.copy()
        c.insert(0, "threshold", thr)
        combined.append(c)

        # ---- per-group summary for this threshold ----
        grp = (
            m.groupby("final_group")
            .agg(
                n_structures=("structure_id", "nunique"),
                n_subclusters=("nanobody_subcluster_id", "nunique"),
                n_medoids=("is_medoid", "sum"),
            )
            .reset_index()
            .rename(columns={"final_group": "antigen_clustering_label"})
            .sort_values(["n_structures", "antigen_clustering_label"], ascending=[False, True])
        )
        summary_sheets[thr] = grp

        n_structures = m["structure_id"].nunique()
        n_reps = m["nanobody_subcluster_id"].nunique()
        overview_rows.append({
            "threshold": thr,
            "n_input_structures": n_structures,
            "n_biological_groups": m["final_group"].nunique(),
            "n_representatives": n_reps,
            "reduction_percent": round((n_structures - n_reps) / n_structures * 100, 2),
        })

        n_missing = out["nanobody_seq"].isna().sum()
        print(f"threshold {thr:>2}: {len(out)} nanobodies, "
              f"{out['nanobody_HAC_label'].nunique()} subclusters"
              + (f"  [WARN {n_missing} missing seq]" if n_missing else "")
              + f"  -> {out_file}")

    if combined:
        alldf = pd.concat(combined, ignore_index=True)
        all_file = SWEEP_DIR / "member_tables_all_thresholds.tsv"
        alldf.to_csv(all_file, sep="\t", index=False)
        print(f"\nCombined ({len(alldf)} rows across {len(combined)} thresholds) -> {all_file}")

        # member workbook: one sheet per threshold (t5, t10, ... t50)
        xlsx_file = SWEEP_DIR / "member_tables_all_thresholds.xlsx"
        with pd.ExcelWriter(xlsx_file, engine="openpyxl") as writer:
            for thr, df in sheets.items():
                df.to_excel(writer, sheet_name=f"t{thr}", index=False)
        print(f"Member workbook ({len(sheets)} sheets) -> {xlsx_file}")

        # summary workbook: Overview sheet + one per-group sheet per threshold
        overview_df = pd.DataFrame(overview_rows)
        summary_xlsx = SWEEP_DIR / "summary_all_thresholds.xlsx"
        with pd.ExcelWriter(summary_xlsx, engine="openpyxl") as writer:
            overview_df.to_excel(writer, sheet_name="Overview", index=False)
            for thr, df in summary_sheets.items():
                df.to_excel(writer, sheet_name=f"t{thr}", index=False)
        print(f"Summary workbook (Overview + {len(summary_sheets)} sheets) -> {summary_xlsx}")
    else:
        print("\nNo member files found. Run sweep_thresholds.py first.")


if __name__ == "__main__":
    main()