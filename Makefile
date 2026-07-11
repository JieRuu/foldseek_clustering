# =====================================================================
# Nanobody / Antigen structural-clustering pipeline
#
# Dependency graph reconstructed from each script's actual input/output
# paths (not filename guesses). Run `make help` to see targets.
#
# Stages:
#   preprocessing -> qc(parallel) -> clustering -> curation -> benchmark
#
# NOTE on the label file:
#   metadata/antigen_clusters_t0p70_cluster_labels.tsv is the MANUAL label
#   table (the merged output from the label-merge task). It is an input to
#   the curation stage, not something produced by these scripts. Update it
#   by hand / by the merge task before running `make curation`.
# =====================================================================

PY      := python
BASH    := bash

SCRIPTS := scripts
RESULTS := results
META    := metadata

# ---- key files (used as make targets / prerequisites) ----------------
ID_MAP          := $(RESULTS)/qc/id_mapping_qc.csv
ANTIGEN_DIR     := $(RESULTS)/antigen_only_pdb
ANTIGEN_STAMP   := $(ANTIGEN_DIR)/.extract.done
FOLDSEEK_STAMP  := $(RESULTS)/foldseek/.foldseek.done
STRUCT_LEVEL    := $(RESULTS)/foldseek/structure_level/cluster_t0.70_structure_level.tsv
LABELS          := $(META)/antigen_clusters_t0p70_cluster_labels.tsv
FINAL_MAPPING   := $(RESULTS)/curation/final_group_mapping_structure_level.tsv
MASTER_TABLE    := $(RESULTS)/curation/nanobody_group_master_table.tsv
BENCHMARK_STAMP := $(RESULTS)/benchmark/.benchmark.done

.PHONY: all help preprocessing qc clustering curation benchmark diagnostics clean
.DEFAULT_GOAL := help

# =====================================================================
# High-level targets
# =====================================================================
all: benchmark            ## Run the full pipeline end-to-end

preprocessing: $(ANTIGEN_STAMP)                 ## ID mapping + antigen-only PDB extraction
clustering:    $(STRUCT_LEVEL)                  ## Foldseek + structure-level normalization
curation:      $(MASTER_TABLE)                  ## Merge labels + build nanobody master table
benchmark:     $(BENCHMARK_STAMP)               ## Select nanobody representatives (final)

# =====================================================================
# 1. PREPROCESSING
# =====================================================================

# 1a. ID mapping QC  ->  id_mapping_qc.csv   (consumed by extraction)
#     Reads PDB_ROOT / FASTA_ROOT from .env + config/config.yaml
$(ID_MAP):
	$(PY) $(SCRIPTS)/preprocessing/check_id_mapping.py

# 1b. Extract antigen-only PDBs  (needs status==OK rows from id_mapping_qc)
#     Stamp file marks completion because output is a directory of *.pdb
$(ANTIGEN_STAMP): $(ID_MAP)
	$(PY) $(SCRIPTS)/preprocessing/extract_antigen_only_pdb.py
	@touch $@

# =====================================================================
# 2. QC (parallel diagnostics — not on the critical path)
# =====================================================================

# Multi-chain antigen check — reads the extracted antigen-only PDBs.
.PHONY: qc-multichain
qc-multichain: $(ANTIGEN_STAMP)
	$(PY) $(SCRIPTS)/qc/check_multichain_antigens.py

# PDB-vs-FASTA ID comparison — pure diagnostic, depends only on inputs.
.PHONY: qc-idcompare
qc-idcompare:
	$(PY) $(SCRIPTS)/preprocessing/compare_pdb_vs_fasta_ids.py

qc: qc-multichain qc-idcompare                  ## Run all QC diagnostics

# =====================================================================
# 3. CLUSTERING
# =====================================================================

# 3a. Foldseek clustering over thresholds 0.70..0.90
#     Reads antigen-only PDBs; wipes+writes results/foldseek/cluster_t*.tsv
$(FOLDSEEK_STAMP): $(ANTIGEN_STAMP)
	$(BASH) $(SCRIPTS)/clustering/run_foldseek.sh
	@touch $@

# 3b. Normalize the 0.70 chain-level clusters back to structure level.
#     --outdir MUST be results/foldseek/structure_level so build_final_groups
#     can find cluster_t0.70_structure_level.tsv (matches your history cmd 688).
$(STRUCT_LEVEL): $(FOLDSEEK_STAMP)
	$(PY) $(SCRIPTS)/clustering/normalize_foldseek_to_structure_level.py \
		--input $(RESULTS)/foldseek/cluster_t0.70.tsv \
		--outdir $(RESULTS)/foldseek/structure_level

# 3c. Threshold sensitivity summary — diagnostic, reads all 5 cluster_t*.tsv
.PHONY: clustering-summary
clustering-summary: $(FOLDSEEK_STAMP)
	$(PY) $(SCRIPTS)/clustering/summarize_foldseek.py

# =====================================================================
# 4. CURATION
# =====================================================================

# 4a. Build final group mapping (MAIN branch).
#     Uses Foldseek structure-level clusters as the backbone and merges the
#     manual labels onto them; unlabelled structures -> label_source=pending_review.
#     Produces final_group_id + label_source that the benchmark step requires.
#
#     ALTERNATIVE (not wired in): merge_biological_labels_structure_level.py
#     produces the same output filename but WITHOUT final_group_id/label_source
#     and without reading the Foldseek structure-level table. Swap in only if
#     you deliberately want the label-driven (non-Foldseek-backbone) grouping.
$(FINAL_MAPPING): $(STRUCT_LEVEL) | $(LABELS)
	$(PY) $(SCRIPTS)/curation/build_final_groups.py

# The manual label table is a hand-maintained INPUT (output of the label-merge
# task), not a pipeline product. Treat it as order-only (| prerequisite) and
# fail loudly with a clear message if it's missing, instead of make trying to
# "build" it.
$(LABELS):
	@echo "ERROR: manual label file missing: $@"; \
	 echo "  This is a hand-maintained input (the merged cluster-label table)."; \
	 echo "  Create/update it before running 'make curation'."; \
	 exit 1

# 4b. Build the nanobody group master table (attaches nanobody sequences).
$(MASTER_TABLE): $(FINAL_MAPPING)
	$(PY) $(SCRIPTS)/curation/build_nanobody_group_master_table.py

# =====================================================================
# 5. BENCHMARK (final deliverable)
# =====================================================================

# 5a. Select nanobody representatives via HAC + medoid (mutation threshold 5).
$(BENCHMARK_STAMP): $(MASTER_TABLE)
	$(PY) $(SCRIPTS)/curation/select_nanobody_representatives.py
	@touch $@

# 5b. Mutant-threshold sensitivity summary (counts only, thresholds 1..10).
.PHONY: benchmark-sensitivity
benchmark-sensitivity: $(MASTER_TABLE)
	$(PY) $(SCRIPTS)/curation/mutant_threshold_sensitivity_analysis.py

# 5c. Threshold sweep (thresholds 5..50, step 5) — writes medoid representatives
#     AND cluster members per threshold into results/benchmark/threshold_sweep/.
#     Use this when you need the actual members at each threshold, not just counts.
.PHONY: threshold-sweep
threshold-sweep: $(MASTER_TABLE)                ## Medoid + members for thresholds 5..50
	$(PY) $(SCRIPTS)/curation/sweep_thresholds.py

# =====================================================================
# Aggregate diagnostics + housekeeping
# =====================================================================
diagnostics: qc clustering-summary benchmark-sensitivity  ## Run every non-critical-path diagnostic

clean:                                          ## Remove derived results (keeps metadata/ + scripts/)
	rm -rf $(RESULTS)/antigen_only_pdb $(RESULTS)/foldseek \
	       $(RESULTS)/curation $(RESULTS)/benchmark
	rm -f  $(RESULTS)/qc/*.csv $(RESULTS)/qc/*.tsv

help:                                           ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| sort \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2}'
