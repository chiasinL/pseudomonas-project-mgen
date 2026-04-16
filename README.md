# Comparative Genomics Pipeline

This repository contains the analysis code accompanying the manuscript:

> *[Citation to be added upon publication]*

This repository covers the hard-to-reproduce post-processing steps for two independent comparative genomics methods — [**anvi'o pangenomics**](https://merenlab.org/2016/11/08/pangenomics-v2/) and
[**progressiveMauve whole-genome alignment**](https://darlinglab.org/mauve/user-guide/progressivemauve.html) — used to identify genes and
genomic regions present in all positive-group genomes but absent from all
negative-group genomes. Results from both methods are cross-referenced and
merged into a unified report.

> **Note:** This code is provided for reproducibility purposes only. It is not
> actively maintained and issues or pull requests are not monitored.

---

## Requirements

- Python ≥ 3.10
- pandas ≥ 2.0
- numpy ≥ 1.24
- openpyxl ≥ 3.1

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Input files

Place input files in the `data/` directory (real data) or use `example_data/`
to test with the provided anonymised dataset.

| File | Description |
|---|---|
| `genomes_metadata.tsv` | Genome names, activity labels (`positive`/`negative`), and backbone column positions |
| `gc_occurrence.tsv`* | Gene cluster × genome count matrix from anvi'o pangenome analysis |
| `pangenome_gene_clusters_summary.tsv` | Per-gene anvi'o gene cluster membership and amino acid sequences |
| `mauve_backbone.backbone` | progressiveMauve backbone output (locally collinear blocks) |
| `reference_genome_annotations.tsv` | Gene coordinates and descriptions for the reference genome (GTF-derived) |
| `reference_prodigal_annotations.tsv` | Prodigal gene calls for the reference genome with COG/KO annotations |

*`gc_occurrence.tsv` is generated following the instructions [here](https://groups.google.com/g/anvio/c/w6MB__H0jZc).


### genomes_metadata.tsv format

```
genome          activity    backbone_position
strain_pos_ref  positive    0
strain_pos_2    positive    1
...
strain_neg_1    negative    4
```

- `backbone_position` maps each genome to its column pair in the Mauve backbone
  file (position 0 = reference genome, columns 0–1; position N = columns 2N–(2N+1)).
- The genome with `backbone_position == 0` is used as the reference throughout
  the pipeline.

---

## How to run

Run scripts in order. All scripts accept `--data` and `--results` arguments to
specify input and output directories (defaults: `data/` and `results/`).

```bash
# Using real data
python 01_process_anvio.py
python 02_process_mauve.py
python 03_integrate_results.py
python 04_merge_results.py

# Using example data
python 01_process_anvio.py --data example_data
python 02_process_mauve.py --data example_data
python 03_integrate_results.py --data example_data
python 04_merge_results.py
```

---

## Output files

All output files are written to `results/`.

| Script | Output file | Description |
|---|---|---|
| 01 | `gc_exclusive_to_positive_annotated.tsv` | Gene-level details for gene clusters exclusive to the positive group (reference genome) |
| 02 | `mauve_gene_level.tsv` | Genes in LCBs exclusive to the positive group, with reference genome annotations |
| 03 | `mauve_gene_level_with_anvio.tsv` | Mauve genes annotated with COG/KO and anvi'o GC membership flag |
| 03 | `gc_exclusive_with_mauve_flag.tsv` | anvi'o exclusive GC genes with Mauve overlap flag, coordinates, and COG/KO annotations |
| 04 | `merged_results.xlsx` | Unified report combining both methods; each gene labelled `both`, `mauve`, or `anvio` |

### Expected row counts (example_data run)

| Output file | Expected rows |
|---|---|
| `gc_exclusive_to_positive_annotated.tsv` | 98 |
| `mauve_gene_level.tsv` | 223 |
| `mauve_gene_level_with_anvio.tsv` | 223 |
| `gc_exclusive_with_mauve_flag.tsv` | 98 |
| `merged_results.xlsx` | 188 (87 mauve, 74 both, 27 anvio) |

---

## Example data

The `example_data/` directory contains anonymised versions of the real input
files for testing. Genome names have been replaced with generic labels
(`strain_pos_ref`, `strain_pos_2`, …, `strain_neg_1`, …), and gene identifiers
with sequential labels (`gene_00001`, …). Large files are subsampled to include
all biologically relevant rows (the exclusive gene clusters) plus 500 randomly
selected additional rows.

To regenerate `example_data/` from `data/`:

```bash
python make_dummy_data.py
```

---

## AI assistance

Code in this repository was developed with assistance from Claude Code (Anthropic, claude-sonnet-4-6), used for code refactoring, documentation, and repository organization.