# Input data dictionary

The four TSV inputs belong in this directory when authorized for the relevant use. They are **not included in this code bundle** and are intended for public deposit upon publication, subject to approval. The scripts read tab-separated files with one header row. Researchers who have the input files can run the module locally.

## Count matrices: `MEC1_counts.tsv`, `HG3_counts.tsv`

The exact same headers are used independently in each cell line:

| Column | Meaning |
| --- | --- |
| `element_id` | Unique oligonucleotide/construct ID (text) |
| `DNA_Rep1`–`DNA_Rep4` | Nonnegative integer DNA counts across four replicates |
| `RNA_Rep1`–`RNA_Rep4` | Nonnegative integer RNA counts across four replicates |

Example (illustrative row provided during input design, **not an independent test dataset**):

```tsv
element_id	DNA_Rep1	DNA_Rep2	DNA_Rep3	DNA_Rep4	RNA_Rep1	RNA_Rep2	RNA_Rep3	RNA_Rep4
rs10005419_2_ref_C	472	446	390	323	381	491	505	534
```

The count files must contain **pre-QC counts** so that the `input > 0` and per-replicate DNA thresholds can be applied by `01_quality_control.R`.

## Shared input-library counts: `input_counts.tsv`

```tsv
element_id	input
rs10005419_2_ref_C	2771
```

Every `element_id` present in either cell-line count matrix must be present here; additional IDs in the input file are harmless but unnecessary. IDs must use the same naming scheme across files. The script checks unique IDs and nonnegative integer counts.

## Annotations: `sequence_annotations.tsv`

Required columns generated from the design workbook:

```tsv
element_id	category	variant_id	variant_class	chromosome	start	end	tested_allele
rs10005419_2_ref_C	Candidates	rs10005419	SNP	chr4	15767132	15767132	C
```

`variant_id` contains the tested rsID. `tested_allele` is the recorded nucleotide or indel allele (e.g. `C`, `A`, `len3bp`). The allelic analysis groups constructs by a **`variant_core` parsed from `element_id`**; its allele suffix indicates reference (`_ref_`) or alternate (`_alt_`, `_alt1_`, `_alt2_`, etc.). Constructs with naming conventions that cannot be parsed are reported, not used in allele tests.

Oligonucleotide IDs must stay consistent between the annotation file and both count matrices. The two cell lines need not contain identical sets of retained oligonucleotides after QC.

## Required provenance information before public release

Complete the following in this file before sharing the data:

- [ ] Genome assembly used for genomic coordinates (for example hg19 or hg38) — **not yet recorded**.
- [ ] Coordinate convention (0-based/1-based, inclusive/half-open) — **not yet recorded**.
- [ ] Source of oligonucleotide design and mapping from original Twist IDs to standardized IDs.
- [ ] Source and processing steps that generated DNA, RNA, and shared input-library counts, including any upstream filters.
- [ ] Sample/replicate description and library design.
- [ ] Permission for public redistribution of all supplied inputs and annotations.
- [ ] Appropriate data accession(s) or contact procedure if any raw/read-level data will be available elsewhere.

Note: this repository begins at the **processed count-matrix** stage; it does not recreate alignment/count generation from FASTQ files.
