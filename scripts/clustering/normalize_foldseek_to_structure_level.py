#!/usr/bin/env python3

import argparse
import re
import pandas as pd
from pathlib import Path


def to_structure_id(x: str) -> str:
    """
    Convert Foldseek chain-level IDs back to structure-level IDs.

    Example:
    6ehg_antigen_only_A -> 6ehg_antigen_only
    5tok_1_antigen_only_C -> 5tok_1_antigen_only
    """
    m = re.match(r"^(.+_antigen_only)_[A-Za-z0-9]+$", x)
    if m:
        return m.group(1)
    return x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
        help="Foldseek cluster TSV, e.g. results/foldseek/cluster_t0.70.tsv",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for normalized structure-level files",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(
        input_path,
        sep="\t",
        header=None,
        names=["representative", "member"],
    )

    df["representative_structure"] = df["representative"].map(to_structure_id)
    df["member_structure"] = df["member"].map(to_structure_id)

    # Full mapping: original Foldseek ID -> normalized structure ID
    mapping = (
        pd.concat(
            [
                df[["representative", "representative_structure"]]
                .rename(columns={
                    "representative": "foldseek_id",
                    "representative_structure": "structure_id",
                }),
                df[["member", "member_structure"]]
                .rename(columns={
                    "member": "foldseek_id",
                    "member_structure": "structure_id",
                }),
            ],
            ignore_index=True,
        )
        .drop_duplicates()
        .sort_values(["structure_id", "foldseek_id"])
    )

    mapping.to_csv(
        outdir / "foldseek_chain_to_structure_mapping.tsv",
        sep="\t",
        index=False,
    )

    # Normalized cluster table
    normalized = (
        df[["representative_structure", "member_structure"]]
        .drop_duplicates()
        .rename(columns={
            "representative_structure": "cluster_representative",
            "member_structure": "structure_id",
        })
        .sort_values(["cluster_representative", "structure_id"])
    )

    normalized.to_csv(
        outdir / "cluster_t0.70_structure_level.tsv",
        sep="\t",
        index=False,
    )

    # QC: check whether chains from the same structure were assigned to multiple clusters
    qc = (
        df.groupby("member_structure")
        .agg(
            n_foldseek_entries=("member", "nunique"),
            n_structure_level_clusters=("representative_structure", "nunique"),
            foldseek_entries=("member", lambda x: ";".join(sorted(set(x)))),
            structure_level_clusters=(
                "representative_structure",
                lambda x: ";".join(sorted(set(x))),
            ),
        )
        .reset_index()
        .rename(columns={"member_structure": "structure_id"})
        .sort_values(["n_structure_level_clusters", "n_foldseek_entries", "structure_id"])
    )

    qc.to_csv(
        outdir / "chain_split_qc_structure_level.tsv",
        sep="\t",
        index=False,
    )

    conflicts = qc[qc["n_structure_level_clusters"] > 1]

    print("Input Foldseek rows:", len(df))
    print("Unique Foldseek IDs:", mapping["foldseek_id"].nunique())
    print("Unique structure-level IDs:", mapping["structure_id"].nunique())
    print("Normalized cluster rows:", len(normalized))
    print("Structure-level clusters:", normalized["cluster_representative"].nunique())
    print("Chain-split conflicts:", len(conflicts))

    if len(conflicts) > 0:
        print("\nWARNING: Some chain-level entries from the same structure map to different clusters.")
        print(conflicts.to_string(index=False))
    else:
        print("\nOK: No chain-split structure maps to multiple structure-level clusters.")


if __name__ == "__main__":
    main()