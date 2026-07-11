# Nanobody–Antigen Structural Clustering Pipeline

A two-stage clustering pipeline that builds a diverse, biologically representative
benchmark set of nanobody–antigen complexes for evaluating computational
structural-biology tools (e.g. binding-affinity predictors, AlphaFold Initial Guess).

The pipeline (1) clusters antigens by **3D structure** using Foldseek, merges those
clusters into biological groups using manual antigen labels, then (2) clusters the
nanobodies within each group by **sequence** (Levenshtein distance + hierarchical
agglomerative clustering) and selects a medoid representative per subcluster.

---

## Why two clustering strategies?

| Stage | What is clustered | Method | Rationale |
|-------|-------------------|--------|-----------|
| 1 | Antigens | Structural (Foldseek, TM-align) | Antigens define the binding target space; 3D fold and surface topology matter more than sequence identity. |
| 2 | Nanobodies | Sequence (Levenshtein + HAC) | The VHH backbone is nearly identical across nanobodies; functional diversity lives in the CDR sequence, so edit distance captures paratope diversity best. |

HAC is used (not k-means) because edit distances give a pairwise distance matrix,
not Cartesian coordinates — HAC operates natively on a distance matrix and lets you
cut the dendrogram at a biologically meaningful mutation threshold.

---

## Pipeline overview

```
raw PDB complexes
      │  preprocessing
      ▼
antigen-only PDBs  ──────────────┐
      │  clustering              │ qc (parallel)
      ▼                          ▼
Foldseek clusters          multichain / id checks
      │  normalize to structure level
      ▼
structure-level clusters
      │  curation: merge manual antigen labels
      ▼
biological groups (43)
      │  [optional] merge into superclusters (32)
      ▼
nanobody HAC + medoid per group
      ▼
benchmark representative set
```

---

## Repository layout

```
foldseek_clustering/
├── Makefile                 # runs the pipeline in dependency order (see below)
├── config/                  # config.yaml (paths, thresholds)
├── .env / .env.example      # PDB_ROOT, FASTA_ROOT (not committed)
├── environment.yml          # conda environment
├── metadata/                # manual label tables, synonym tables
├── scripts/
│   ├── preprocessing/       # ID matching, antigen-only PDB extraction
│   ├── qc/                  # quality-control checks
│   ├── clustering/          # Foldseek run + structure-level normalization
│   └── curation/            # label merge, master table, representative selection
└── results/                 # all generated outputs (git-ignored)
```

---

## Quick start

```bash
conda env create -f environment.yml
conda activate nanobody-af2
cp .env.example .env          # then fill in PDB_ROOT and FASTA_ROOT

make help                     # list all pipeline targets
make all                      # run the full pipeline end-to-end
```

Run individual stages (each auto-runs missing upstream steps):

```bash
make preprocessing   # ID mapping + antigen-only PDB extraction
make clustering      # Foldseek + structure-level normalization
make curation        # merge labels + build nanobody master table
make benchmark       # select nanobody representatives (final output)
make diagnostics     # all non-critical-path QC / summary scripts
```

---

## Scripts

### preprocessing/
| Script | Input | Output |
|--------|-------|--------|
| `compare_pdb_vs_fasta_ids.py` | PDB folders, nanobody/antigen FASTA | `qc/pdb_vs_fasta_id_comparison.csv` (diagnostic) |
| `check_id_mapping.py` | PDB folders, FASTA | `qc/id_mapping_qc.csv` |
| `extract_antigen_only_pdb.py` | `id_mapping_qc.csv`, nanobody FASTA | `antigen_only_pdb/*.pdb` + QC report |

Nanobody chains are identified by normalized Levenshtein similarity to the reference
nanobody sequence, then removed to leave antigen-only structures.

### qc/
| Script | Input | Output |
|--------|-------|--------|
| `check_multichain_antigens.py` | `antigen_only_pdb/*.pdb` | `qc/multichain_antigen_qc.tsv` |

### clustering/
| Script | Input | Output |
|--------|-------|--------|
| `run_foldseek.sh` | `antigen_only_pdb/*.pdb` | `foldseek/cluster_t{0.70..0.90}.tsv` |
| `normalize_foldseek_to_structure_level.py` | `cluster_t0.70.tsv` | `foldseek/structure_level/cluster_t0.70_structure_level.tsv` + chain→structure map |
| `summarize_foldseek.py` | all `cluster_t*.tsv` | `foldseek/threshold_cluster_summary.tsv` (diagnostic) |

Foldseek runs with `--alignment-type 1` (TM-align) and `-c 0.9` coverage across five
min-seq-id thresholds. Multi-chain antigens get split into chain-level entries during
DB creation; `normalize_foldseek_to_structure_level.py` collapses them back to
structure-level IDs.

### curation/
| Script | Input | Output |
|--------|-------|--------|
| `build_final_groups.py` | structure-level clusters, `metadata/` label + synonym tables | `curation/final_group_mapping_structure_level.tsv` (+ xlsx, summary, QC) |
| `build_nanobody_group_master_table.py` | final group mapping, nanobody FASTA | `curation/nanobody_group_master_table.tsv` |
| `select_nanobody_representatives.py` | master table | `benchmark/*` (representatives, subclusters, summary) |
| `mutant_threshold_sensitivity_analysis.py` | master table | `qc/nanobody_threshold_sensitivity*.tsv` (diagnostic) |

`build_final_groups.py` uses the Foldseek structure clusters as the backbone and
merges manual antigen labels onto them; unlabelled structures are flagged
`pending_review` rather than force-assigned. Representative selection runs
average-linkage HAC on the Levenshtein distance matrix per group and picks the
medoid (a real structure, never a synthetic centroid) of each subcluster.

> The mutation threshold (default 5) is adjustable — lower it to be more sensitive
> to CDR point mutations, raise it to keep only distinct nanobody lineages.

---

## Configuration

Paths and thresholds live in `config/config.yaml` and `.env`:
- `.env`: `PDB_ROOT`, `FASTA_ROOT` (absolute paths to raw data, not committed)
- `config.yaml`: results root, input FASTA names, preprocessing similarity threshold

---

## Dependencies

- `foldseek`
- Python: `pandas`, `numpy`, `scipy`, `Levenshtein`, `biopython`, `pyyaml`, `python-dotenv`, `openpyxl`

See `environment.yml` for the full pinned environment.

---

## Status

Active research project. The antigen-grouping and nanobody-representative stages are
complete; downstream negative-set construction and tool benchmarking are in progress.
