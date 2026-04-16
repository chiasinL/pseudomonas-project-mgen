"""
Merge anvi'o and progressiveMauve results into a single Excel report.

Each gene is assigned a status:
    "both"  — detected by both methods
    "mauve" — detected by progressiveMauve only
    "anvio" — detected by anvi'o only

Inputs  (results/):
    results/mauve_gene_level_with_anvio.tsv  — output of 03_integrate_results.py
    results/gc_exclusive_with_mauve_flag.tsv — output of 03_integrate_results.py

Output:
    results/merged_results.xlsx
"""

import argparse
import pandas as pd

# ── argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--results", default="results", help="Output directory (default: results)")
args = parser.parse_args()

RESULTS = args.results

# ── load inputs ───────────────────────────────────────────────────────────────
mauve_df = pd.read_csv(f"{RESULTS}/mauve_gene_level_with_anvio.tsv", sep="\t")
anvio_df = pd.read_csv(f"{RESULTS}/gc_exclusive_with_mauve_flag.tsv", sep="\t")

# ── prepare Mauve table ───────────────────────────────────────────────────────
mauve_cols = [
    "Gene_ID", "Loci", "gene_level_start", "gene_level_end",
    "GTF_annotation", "Gene_Symbol", "Gene_Synonym",
    "COG_accession", "COG_function", "KO_accession", "KEGG_function",
    "In_anvio_GCs", "gene_callers_id", "gene_cluster_id",
]
mauve_out = (
    mauve_df[mauve_cols]
    .rename(columns={
        "gene_level_start": "Start_Coordinate",
        "gene_level_end":   "End_Coordinate",
    })
    .assign(Status=lambda df: df["In_anvio_GCs"].eq("yes").map({True: "both", False: "mauve"}))
    .drop(columns="In_anvio_GCs")
)
# Prodigal-only genes (no GTF entry): use prodigal_{gene_callers_id} as the identifier
prodigal_only = mauve_out["Gene_ID"].isna() & mauve_out["gene_callers_id"].notna()
mauve_out.loc[prodigal_only, "Gene_ID"] = "prodigal_" + mauve_out.loc[prodigal_only, "gene_callers_id"].astype(str)
# Drop rows with no gene identity at all (truly intergenic LCBs with no Prodigal match)
mauve_out = mauve_out[mauve_out["Gene_ID"].notna()]

# ── prepare anvi'o-only table (genes not found by Mauve) ──────────────────────
anvio_cols = [
    "gene_cluster_id", "gene_callers_id", "No_of_genes",
    "chrom_start", "chrom_end",
    "COG_accession", "COG_function", "KO_accession", "KEGG_function",
    "In_mauve",
]
anvio_out = (
    anvio_df[anvio_cols]
    .rename(columns={"chrom_start": "Start_Coordinate", "chrom_end": "End_Coordinate"})
    .assign(
        Status=lambda df: df["In_mauve"].eq("yes").map({True: "both", False: "anvio"}),
        gene_callers_id=lambda df: df["gene_callers_id"].astype(str),
    )
    .drop(columns="In_mauve")
    # keep only genes not already represented in the Mauve table
    .query("Status == 'anvio'")
)

# ── merge and export ──────────────────────────────────────────────────────────
merged = pd.concat([mauve_out, anvio_out], ignore_index=True)

assert len(merged) > 0, "Merged table is empty — check that input files are not empty."
print(f"Merged rows: {len(merged)}")
print(merged["Status"].value_counts().to_string())

out_path = f"{RESULTS}/merged_results.xlsx"
merged.to_excel(out_path, index=False, engine="openpyxl")
print(f"Output written to {out_path}")
