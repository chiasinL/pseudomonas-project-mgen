"""
Generate anonymised example data for testing the pipeline.

Reads real input files from data/ and writes anonymised versions to example_data/.
Genome names are replaced with generic labels, gene IDs with sequential identifiers,
and the chromosome name with a generic placeholder. Large files are subsampled to
keep example_data/ lightweight while preserving all biologically relevant rows
(i.e. all exclusive gene clusters are retained).

Usage:
    python make_dummy_data.py
"""

import pandas as pd
import numpy as np

DATA_DIR = "data"
OUT_DIR = "example_data"
RANDOM_SEED = 42
N_EXTRA_GCS = 500  # random non-exclusive GCs to include alongside the real ones

GENOME_MAP = {
    "P_protegens_PBL3": "strain_pos_ref",
    "P_fluo_513":        "strain_pos_2",
    "P_fluo_1550":       "strain_pos_3",
    "P_protegens_3295":  "strain_pos_4",
    "P_fluo_1098":       "strain_neg_1",
    "P_fluo_4488":       "strain_neg_2",
    "P_putida_4512":     "strain_neg_3",
}
CHROM_MAP = {"NZ_CP051673.1": "chromosome_1"}


def anonymise_genome_names(s):
    for real, dummy in GENOME_MAP.items():
        s = s.replace(real, dummy)
    return s


# --- genomes_metadata.tsv ---
meta = pd.read_csv(f"{DATA_DIR}/genomes_metadata.tsv", sep="\t")
meta["genome"] = meta["genome"].map(GENOME_MAP)
meta.to_csv(f"{OUT_DIR}/genomes_metadata.tsv", sep="\t", index=False)
print("genomes_metadata.tsv done")

# --- gc_occurrence.tsv ---
# Identify exclusive GCs (present in all positive, absent in all negative) to always retain
gc = pd.read_csv(f"{DATA_DIR}/gc_occurrence.tsv", sep="\t", index_col=0)
positive = meta.loc[meta["activity"] == "positive", "genome"].tolist()
negative = meta.loc[meta["activity"] == "negative", "genome"].tolist()

# Map column names back to originals for filtering (gc still has real names at this point)
real_positive = [k for k, v in GENOME_MAP.items() if v in positive]
real_negative  = [k for k, v in GENOME_MAP.items() if v in negative]

exclusive_mask = (gc[real_positive] > 0).all(axis=1) & (gc[real_negative] == 0).all(axis=1)
exclusive_gcs = gc[exclusive_mask]

rng = np.random.default_rng(RANDOM_SEED)
non_exclusive = gc[~exclusive_mask]
sample_idx = rng.choice(len(non_exclusive), size=min(N_EXTRA_GCS, len(non_exclusive)), replace=False)
gc_sample = pd.concat([exclusive_gcs, non_exclusive.iloc[sample_idx]])
gc_sample.columns = [GENOME_MAP.get(c, c) for c in gc_sample.columns]
gc_sample.index.name = "gene_cluster_id"
gc_sample.to_csv(f"{OUT_DIR}/gc_occurrence.tsv", sep="\t")
print(f"gc_occurrence.tsv done ({len(exclusive_gcs)} exclusive + {min(N_EXTRA_GCS, len(non_exclusive))} random GCs)")

# --- pangenome_gene_clusters_summary.tsv ---
# Keep only rows for sampled GCs; anonymise genome_name
kept_gc_ids = set(gc_sample.index)
pan = pd.read_csv(f"{DATA_DIR}/pangenome_gene_clusters_summary.tsv", sep="\t")
pan = pan[pan["gene_cluster_id"].isin(kept_gc_ids)].copy()
pan["genome_name"] = pan["genome_name"].map(GENOME_MAP)
pan.to_csv(f"{OUT_DIR}/pangenome_gene_clusters_summary.tsv", sep="\t", index=False)
print(f"pangenome_gene_clusters_summary.tsv done ({len(pan)} rows)")

# --- reference_genome_annotations.tsv ---
# Anonymise gene IDs (index) and any occurrences of the chromosome name
gtf = pd.read_csv(f"{DATA_DIR}/reference_genome_annotations.tsv", sep="\t", index_col=0)
gene_id_map = {gid: f"gene_{i+1:05d}" for i, gid in enumerate(gtf.index)}
gtf.index = gtf.index.map(gene_id_map)
gtf.index.name = "gene_id"
# description may contain gene IDs or chromosome references — leave functional text intact
gtf.to_csv(f"{OUT_DIR}/reference_genome_annotations.tsv", sep="\t")
print(f"reference_genome_annotations.tsv done ({len(gtf)} genes)")

# --- reference_prodigal_annotations.tsv ---
# gene_callers_id is already numeric (anonymous); no strain names present
prodigal = pd.read_csv(f"{DATA_DIR}/reference_prodigal_annotations.tsv", sep="\t")
prodigal.to_csv(f"{OUT_DIR}/reference_prodigal_annotations.tsv", sep="\t", index=False)
print(f"reference_prodigal_annotations.tsv done (copied as-is, {len(prodigal)} rows)")

# --- mauve_backbone.backbone ---
# Contains only numeric coordinates; no names to anonymise — copy as-is
import shutil
shutil.copy(f"{DATA_DIR}/mauve_backbone.backbone", f"{OUT_DIR}/mauve_backbone.backbone")
print("mauve_backbone.backbone done (copied as-is)")

print("\nAll example_data files written to:", OUT_DIR)
