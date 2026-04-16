"""
Process progressiveMauve backbone output to identify conserved genomic regions
exclusive to the positive-group genomes, then aggregate those regions to gene level.

Inputs  (data/ or example_data/):
    genomes_metadata.tsv          — genome names, activity labels, backbone positions
    mauve_backbone.backbone       — progressiveMauve locally collinear block (LCB) output
    reference_genome_annotations.tsv — gene coordinates and descriptions for the reference genome

Output:
    results/mauve_gene_level.tsv  — gene-level LCBs exclusive to the positive group,
                                    with reference genome annotations

Notes:
    - The backbone file encodes each genome as a pair of columns (left end, right end).
      backbone_position in genomes_metadata.tsv maps each genome to its column pair:
      position 0 → columns 0–1 (reference), position N → columns 2N–(2N+1).
    - Signed coordinates indicate strand: negative values mean the LCB is on the
      reverse strand in that genome.
    - Requires Python >= 3.10 (structural pattern matching is not used here, but
      progressiveMauve itself requires that environment).
"""

import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

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
    reference = df.loc[df["backbone_position"] == 0, "genome"].iloc[0]
    # map backbone_position → genome name (excludes reference at position 0)
    col_to_genome = (
        df.dropna(subset=["backbone_position"])
        .astype({"backbone_position": int})
        .set_index("backbone_position")["genome"]
        .to_dict()
    )
    return reference, positive, negative, col_to_genome


reference, positive_genomes, negative_genomes, col_to_genome = load_metadata(
    f"{DATA}/genomes_metadata.tsv"
)
# number of positive genomes excluding the reference (used as the conservation threshold)
n_other_positive = len(positive_genomes) - 1


# ── part 1: parse backbone → LCB regions present in the reference genome ──────
gtf_df = pd.read_csv(f"{DATA}/reference_genome_annotations.tsv", sep="\t", index_col=0)
gtf_dict = gtf_df.to_dict(orient="index")

backbone_rows = []
with open(f"{DATA}/mauve_backbone.backbone") as fh:
    next(fh)  # skip header
    for line in fh:
        fields = line.strip().split("\t")
        ref_left, ref_right = int(fields[0]), int(fields[1])

        # skip LCBs not present in the reference genome
        if ref_left == 0 and ref_right == 0:
            continue

        pos_genomes, neg_genomes = [], []
        for pos, genome in col_to_genome.items():
            if pos == 0:
                continue  # reference already handled
            col = pos * 2  # backbone_position N → column pair (2N, 2N+1)
            left, right = int(fields[col]), int(fields[col + 1])
            if left != 0 or right != 0:
                if genome in positive_genomes:
                    pos_genomes.append(genome)
                else:
                    neg_genomes.append(genome)

        locus_id = f"ref:{abs(ref_left)}-{abs(ref_right)}"
        backbone_rows.append({
            "Locus":                    locus_id,
            "No_Overlapped_Positive":   len(pos_genomes),
            "Overlapped_Positive":      ",".join(pos_genomes),
            "No_Overlapped_Negative":   len(neg_genomes),
            "Overlapped_Negative":      ",".join(neg_genomes),
            "Start_Coordinate":         ref_left,
            "End_Coordinate":           ref_right,
        })

backbone_df = pd.DataFrame(backbone_rows)
assert len(backbone_df) > 0, "No LCBs found in backbone file for the reference genome."
print(f"LCBs present in reference genome: {len(backbone_df)}")


# ── part 2: annotate LCBs with reference genome gene models ──────────────────
annotated_rows = {}
for _, row in backbone_df.iterrows():
    locus = row["Locus"]
    r_start = abs(row["Start_Coordinate"])
    r_end   = abs(row["End_Coordinate"])
    overlapping_genes = [
        gid for gid, g in gtf_dict.items()
        if (g["chrom_start"] <= r_start <= g["chrom_end"] or
            g["chrom_start"] <= r_end   <= g["chrom_end"] or
            r_start <= g["chrom_start"] <= r_end or
            r_start <= g["chrom_end"]   <= r_end)
    ]
    base = row.to_dict()
    if len(overlapping_genes) == 0:
        annotated_rows[locus] = base
    elif len(overlapping_genes) == 1:
        gid = overlapping_genes[0]
        annotated_rows[locus] = {**base,
            "GTF_annotation": gtf_dict[gid].get("description"),
            "Gene_ID":        gid,
            "Gene_Symbol":    gtf_dict[gid].get("gene"),
            "Gene_Synonym":   gtf_dict[gid].get("gene_synonym"),
        }
    else:
        # one backbone row spans multiple genes → create one row per gene
        for i, gid in enumerate(overlapping_genes, start=1):
            new_locus = f"{locus}.{i}"
            annotated_rows[new_locus] = {**base,
                "Locus":          new_locus,
                "GTF_annotation": gtf_dict[gid].get("description"),
                "Gene_ID":        gid,
                "Gene_Symbol":    gtf_dict[gid].get("gene"),
                "Gene_Synonym":   gtf_dict[gid].get("gene_synonym"),
            }

annotated_df = pd.DataFrame.from_dict(annotated_rows, orient="index")


# ── part 3: filter to LCBs exclusive to the positive group ───────────────────
# Keep LCBs present in all other positive genomes and absent from all negative genomes.
# n_other_positive excludes the reference itself.
exclusive_df = annotated_df[
    (annotated_df["No_Overlapped_Positive"] == n_other_positive) &
    (annotated_df["No_Overlapped_Negative"] == 0)
].copy()

assert len(exclusive_df) > 0, (
    "No exclusive LCBs found. Check genome names in metadata match the backbone column order."
)
print(f"LCBs exclusive to positive group: {len(exclusive_df)}")


# ── part 4: aggregate LCB regions to gene level ───────────────────────────────
# Multiple backbone rows may map to the same gene; merge them into one row per gene,
# taking the outermost coordinates (clamped to the gene's annotated boundaries).
to_aggregate = (
    exclusive_df[exclusive_df["Gene_ID"].notna()]
    .sort_values(["Gene_ID", "Start_Coordinate"])
)

gene_rows = {}
for gene_id, grp in to_aggregate.groupby("Gene_ID"):
    gene = gtf_dict.get(gene_id, {})
    starts = grp["Start_Coordinate"].values
    ends   = grp["End_Coordinate"].values

    # signed coordinates: negative strand uses negative values.
    # For negative strand, abs(starts) may not be ordered the same as starts —
    # take min/max of absolute values explicitly to get the correct outer boundary.
    is_neg = starts[0] < 0
    if is_neg:
        min_abs_start = min(abs(s) for s in starts)
        max_abs_end   = max(abs(e) for e in ends)
        gene_start = -(max(min_abs_start, gene.get("chrom_start", min_abs_start)))
        gene_end   = -(min(max_abs_end,   gene.get("chrom_end",   max_abs_end)))
    else:
        gene_start = max(starts.min(), gene.get("chrom_start", starts.min()))
        gene_end   = min(ends.max(),   gene.get("chrom_end",   ends.max()))

    first = grp.iloc[0]
    gene_rows[gene_id] = {
        "Gene_ID":                gene_id,
        "Loci":                   ",".join(grp["Locus"].tolist()),
        "gene_level_start":       gene_start,
        "gene_level_end":         gene_end,
        "GTF_annotation":         first.get("GTF_annotation"),
        "Gene_Symbol":            first.get("Gene_Symbol"),
        "Gene_Synonym":           first.get("Gene_Synonym"),
        "chrom_start":            gene.get("chrom_start"),
        "chrom_end":              gene.get("chrom_end"),
    }

# include rows where no gene ID was assigned (intergenic LCBs)
no_gene_df = exclusive_df[exclusive_df["Gene_ID"].isna()][
    ["Locus", "Start_Coordinate", "End_Coordinate",
     "GTF_annotation", "Gene_ID", "Gene_Symbol", "Gene_Synonym"]
].rename(columns={"Locus": "Loci",
                  "Start_Coordinate": "gene_level_start",
                  "End_Coordinate": "gene_level_end"})
no_gene_df["chrom_start"] = np.nan
no_gene_df["chrom_end"]   = np.nan

gene_level_df = pd.concat(
    [pd.DataFrame.from_dict(gene_rows, orient="index"), no_gene_df],
    ignore_index=True
)

# ── save ──────────────────────────────────────────────────────────────────────
out_path = f"{RESULTS}/mauve_gene_level.tsv"
gene_level_df.to_csv(out_path, sep="\t", index=False)
print(f"Output written to {out_path} ({len(gene_level_df)} rows)")
