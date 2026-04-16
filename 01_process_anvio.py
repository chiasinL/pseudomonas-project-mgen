"""
Process anvi'o pangenome output to identify gene clusters (GCs) exclusive to
positive-group genomes.

Inputs  (data/ or example_data/):
    genomes_metadata.tsv              — genome names, activity labels, backbone positions
    gc_occurrence.tsv                 — GC × genome count matrix from anvi'o
    pangenome_gene_clusters_summary.tsv — per-gene GC membership and sequences

Output:
    results/gc_exclusive_to_positive_annotated.tsv — gene-level details for exclusive GCs
                                                      in the reference genome
"""

import pandas as pd
import argparse

# ── argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--data", default="data", help="Input data directory (default: data)")
parser.add_argument("--results", default="results", help="Output directory (default: results)")
args = parser.parse_args()

DATA = args.data
RESULTS = args.results


# ── load genome metadata ──────────────────────────────────────────────────────
def load_metadata(path):
    df = pd.read_csv(path, sep="\t")
    positive = df.loc[df["activity"] == "positive", "genome"].tolist()
    negative = df.loc[df["activity"] == "negative", "genome"].tolist()
    # first positive genome in backbone order (position 0) is the reference
    reference = df.loc[df["backbone_position"] == 0, "genome"].iloc[0]
    return reference, positive, negative


reference, positive_genomes, negative_genomes = load_metadata(f"{DATA}/genomes_metadata.tsv")

# ── load GC occurrence matrix ─────────────────────────────────────────────────
gc_df = pd.read_csv(f"{DATA}/gc_occurrence.tsv", sep="\t", index_col=0)
gc_df.index.name = "gene_cluster_id"

# ── find exclusive GCs ────────────────────────────────────────────────────────
# A GC is exclusive to the positive group when every positive genome has at
# least one gene in it and every negative genome has none.
exclusive_mask = (
    (gc_df[positive_genomes] > 0).all(axis=1) &
    (gc_df[negative_genomes] == 0).all(axis=1)
)
exclusive_gcs = gc_df[exclusive_mask]

assert len(exclusive_gcs) > 0, (
    "No exclusive gene clusters found. Check that genome names in "
    f"{DATA}/gc_occurrence.tsv match those in {DATA}/genomes_metadata.tsv."
)
print(f"Exclusive GCs found: {len(exclusive_gcs)}")

# ── annotate with gene-level details for the reference genome ─────────────────
summary_df = pd.read_csv(f"{DATA}/pangenome_gene_clusters_summary.tsv", sep="\t")
ref_summary = summary_df[summary_df["genome_name"] == reference]

annotated = (
    exclusive_gcs
    .reset_index()
    .merge(ref_summary.drop(columns="num_genes_in_gene_cluster"), on="gene_cluster_id")
)
# Number of reference-genome genes per GC, taken directly from the occurrence matrix
ref_gc_counts = gc_df[reference].rename("No_of_genes")
annotated = annotated.join(ref_gc_counts, on="gene_cluster_id")

assert len(annotated) > 0, (
    f"Merge produced no rows. Check that '{reference}' exists in "
    f"{DATA}/pangenome_gene_clusters_summary.tsv under 'genome_name'."
)

# ── save ──────────────────────────────────────────────────────────────────────
out_path = f"{RESULTS}/gc_exclusive_to_positive_annotated.tsv"
annotated.to_csv(out_path, sep="\t", index=False)
print(f"Output written to {out_path} ({len(annotated)} rows)")
