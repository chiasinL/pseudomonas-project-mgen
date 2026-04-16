"""
Cross-reference anvi'o and progressiveMauve results using coordinate overlap.

For each Mauve gene-level region, checks whether it overlaps with any gene in the
exclusive anvi'o GC set (by genomic coordinates). Adds COG/KO functional annotations
from the prodigal gene calls to all Mauve genes, not just those in exclusive GCs.

Inputs  (results/ + data/ or example_data/):
    results/gc_exclusive_to_positive_annotated.tsv  — output of 01_process_anvio.py
    results/mauve_gene_level.tsv                    — output of 02_process_mauve.py
    data/reference_prodigal_annotations.tsv         — prodigal gene calls with coordinates
                                                      and COG/KO annotations

Outputs:
    results/mauve_gene_level_with_anvio.tsv   — Mauve genes annotated with COG/KO and
                                                anvi'o GC membership where applicable
    results/gc_exclusive_with_mauve_flag.tsv  — anvi'o exclusive GC genes with a flag
                                                indicating overlap with Mauve results
"""

import argparse
import pandas as pd

# ── argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--data", default="data", help="Input data directory (default: data)")
parser.add_argument("--results", default="results", help="Output directory (default: results)")
args = parser.parse_args()

DATA = args.data
RESULTS = args.results


def overlaps(a_start, a_end, b_start, b_end):
    """Return True if intervals [a_start, a_end] and [b_start, b_end] overlap."""
    return (b_start <= a_start <= b_end or
            b_start <= a_end   <= b_end or
            a_start <= b_start <= a_end or
            a_start <= b_end   <= a_end)


# ── load inputs ───────────────────────────────────────────────────────────────
anvio_df   = pd.read_csv(f"{RESULTS}/gc_exclusive_to_positive_annotated.tsv", sep="\t")
mauve_df   = pd.read_csv(f"{RESULTS}/mauve_gene_level.tsv", sep="\t")
prodigal_df = pd.read_csv(f"{DATA}/reference_prodigal_annotations.tsv", sep="\t")

# subset prodigal to genes in exclusive GCs only (for coordinate lookup)
exclusive_gene_ids = set(anvio_df["gene_callers_id"].unique())
prodigal_exclusive = (
    prodigal_df[prodigal_df["gene_callers_id"].isin(exclusive_gene_ids)]
    .set_index("gene_callers_id")
    .to_dict(orient="index")
)

# add gene_cluster_id and No_of_genes to prodigal lookup from anvio output
gc_lookup = (
    anvio_df[["gene_callers_id", "gene_cluster_id", "No_of_genes"]]
    .drop_duplicates("gene_callers_id")
    .set_index("gene_callers_id")
    .to_dict(orient="index")
)
for gid, info in prodigal_exclusive.items():
    info.update(gc_lookup.get(gid, {}))

# full prodigal dict for COG/KO annotation of all Mauve genes
prodigal_all = prodigal_df.set_index("gene_callers_id").to_dict(orient="index")


# ── add COG/KO annotations to Mauve genes via coordinate overlap ──────────────
# For each Mauve gene region, find any prodigal gene whose coordinates overlap.
# Take the first match for annotation (regions typically map to one gene).
ann_cols = ["COG_accession", "COG_function", "KO_accession", "KEGG_function"]
for col in ann_cols:
    mauve_df[col] = None

for idx, row in mauve_df.iterrows():
    if pd.isna(row["gene_level_start"]) or pd.isna(row["gene_level_end"]):
        continue
    # signed coordinates: negative strand uses negative values — normalise to positive range
    r_start = min(abs(row["gene_level_start"]), abs(row["gene_level_end"]))
    r_end   = max(abs(row["gene_level_start"]), abs(row["gene_level_end"]))
    for gid, ginfo in prodigal_all.items():
        if overlaps(r_start, r_end, ginfo["chrom_start"], ginfo["chrom_end"]):
            for col in ann_cols:
                mauve_df.at[idx, col] = ginfo.get(col)
            break


# ── flag Mauve genes that overlap with anvi'o exclusive GC genes ──────────────
mauve_df["In_anvio_GCs"] = None
mauve_df["gene_callers_id"] = None
mauve_df["gene_cluster_id"] = None

anvio_mapped_gene_ids = set()

for idx, row in mauve_df.iterrows():
    if pd.isna(row["gene_level_start"]) or pd.isna(row["gene_level_end"]):
        continue
    r_start = min(abs(row["gene_level_start"]), abs(row["gene_level_end"]))
    r_end   = max(abs(row["gene_level_start"]), abs(row["gene_level_end"]))

    matched_ids, matched_gcs = [], []
    for gid, ginfo in prodigal_exclusive.items():
        if overlaps(r_start, r_end, ginfo["chrom_start"], ginfo["chrom_end"]):
            matched_ids.append(str(gid))
            matched_gcs.append(str(ginfo.get("gene_cluster_id", "")))
            anvio_mapped_gene_ids.add(str(gid))

    if matched_ids:
        mauve_df.at[idx, "In_anvio_GCs"]    = "yes"
        mauve_df.at[idx, "gene_callers_id"] = ", ".join(matched_ids)
        mauve_df.at[idx, "gene_cluster_id"] = ", ".join(matched_gcs)

# For Mauve rows with no GTF gene ID, look up the Prodigal gene covering that
# locus (if any) and record its gene_callers_id. These are genes predicted by
# Prodigal but not annotated by PGAP, so they have no GTF entry.
no_gtf_mask = mauve_df["Gene_ID"].isna() & mauve_df["gene_callers_id"].isna()
for idx, row in mauve_df[no_gtf_mask].iterrows():
    if pd.isna(row["gene_level_start"]) or pd.isna(row["gene_level_end"]):
        continue
    r_start = min(abs(row["gene_level_start"]), abs(row["gene_level_end"]))
    r_end   = max(abs(row["gene_level_start"]), abs(row["gene_level_end"]))
    for gid, ginfo in prodigal_all.items():
        if overlaps(r_start, r_end, ginfo["chrom_start"], ginfo["chrom_end"]):
            mauve_df.at[idx, "gene_callers_id"] = str(gid)
            break

assert mauve_df["In_anvio_GCs"].notna().any(), (
    "No Mauve regions overlapped with anvi'o exclusive GC genes. "
    "Check that both methods used the same reference genome coordinates."
)

out_mauve = f"{RESULTS}/mauve_gene_level_with_anvio.tsv"
mauve_df.to_csv(out_mauve, sep="\t", index=False)
print(f"Output written to {out_mauve} ({len(mauve_df)} rows)")


# ── add In_mauve flag to anvi'o exclusive GC genes ───────────────────────────
anvio_out = anvio_df[["gene_cluster_id", "gene_callers_id", "No_of_genes",
                       "aa_sequence"]].copy()

# attach coordinates and COG/KO from prodigal
coord_ann = prodigal_df.rename(columns={"gene_callers_id": "gene_callers_id"})
anvio_out = anvio_out.merge(
    prodigal_df[["gene_callers_id"] + ann_cols + ["chrom_start", "chrom_end"]],
    on="gene_callers_id",
    how="left"
)

anvio_out["In_mauve"] = anvio_out["gene_callers_id"].astype(str).isin(anvio_mapped_gene_ids).map({True: "yes"})

out_anvio = f"{RESULTS}/gc_exclusive_with_mauve_flag.tsv"
anvio_out.to_csv(out_anvio, sep="\t", index=False)
print(f"Output written to {out_anvio} ({len(anvio_out)} rows)")
