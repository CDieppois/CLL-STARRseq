# STARR-seq analysis

R scripts for quality control, enhancer activity quantification and allele-specific regulatory activity in the CLL-derived MEC-1 and HG-3 cell lines. This module starts from **processed oligonucleotide count matrices**, not raw sequencing reads. The two cell lines are analyzed independently; allelic-effect results are subsequently annotated with enhancer activity in either cell line.

This is one analysis module within the broader `CLL-STARRseq` repository. Its paths, methods and software-environment report are **independent of the QTL and transcription-factor workflows**.

## Directory structure

```text
STARR-seq_analysis/
├── README.md
├── .gitignore
├── data/
│   ├── README.md
│   ├── MEC1_counts.tsv             # supplied separately
│   ├── HG3_counts.tsv              # supplied separately
│   ├── input_counts.tsv            # supplied separately
│   └── sequence_annotations.tsv    # supplied separately
├── scripts/
│   ├── 01_quality_control.R
│   ├── 02_enhancer_activity.R
│   ├── 03_allelic_effects.R
│   ├── 04_annotate_allelic_effects.R
│   └── 05_record_environment.R
├── docs/
│   ├── methods.md
│   ├── validation.md
│   └── release_checklist.md
└── results/                      # generated locally; not committed
```

**Input-data availability:** The four input TSV datasets are not distributed in this code bundle. Processed STARR-seq count matrices and associated annotations are intended to be deposited upon publication of the associated study, subject to approval. The exact input fields and examples are described in the [data dictionary](data/README.md).

## Requirements

R with the CRAN packages `dplyr`, `readr`, `tidyr`, `tibble` and Bioconductor's `MPRAnalyze`. MPRAnalyze **1.26.0** was reported as used in the original analysis; installing the current Bioconductor version does not guarantee numerical equivalence.

```r
install.packages(c("dplyr", "readr", "tidyr", "tibble"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("MPRAnalyze")
```

## Running the module

**Working directory is important.** From the top level of `CLL-STARRseq`, first enter this folder:

```bash
cd STARR-seq_analysis
Rscript scripts/01_quality_control.R
Rscript scripts/02_enhancer_activity.R
Rscript scripts/03_allelic_effects.R
Rscript scripts/04_annotate_allelic_effects.R
Rscript scripts/05_record_environment.R
```

In RStudio, set the working directory to the `STARR-seq_analysis/` folder before running any script. Script paths such as `data/MEC1_counts.tsv` and `results/qc/` are relative to **this directory**, not the top-level repository directory. Output directories are created automatically.

## Inputs

| File in `data/` | Required content |
| --- | --- |
| `MEC1_counts.tsv` | `element_id`; `DNA_Rep1`–`DNA_Rep4`; `RNA_Rep1`–`RNA_Rep4` |
| `HG3_counts.tsv` | Same columns, with HG-3 counts |
| `input_counts.tsv` | `element_id`, `input` (shared input-library counts) |
| `sequence_annotations.tsv` | `element_id`, `category`, `variant_id`, `variant_class`, `chromosome`, `start`, `end`, `tested_allele` |

Input matrices should contain counts **before** this module's QC filtering. Oligonucleotide identifiers must be consistent among files; missing shared input counts cause the QC script to stop rather than silently filling in zeros.

## Analysis steps

| Script | Method | Main output in `results/` |
| --- | --- | --- |
| `01_quality_control.R` | Retain elements with input count > 0 and DNA counts ≥ 5 in **every** replicate; no RNA minimum | `qc/MEC1_counts_filtered.tsv`, `qc/HG3_counts_filtered.tsv`, filtering reports |
| `02_enhancer_activity.R` | MPRAnalyze quantification; MAD-based p-values and Benjamini–Hochberg adjustment | `enhancer_activity/{CELL}_enhancer_activity.tsv`, exclusion reports |
| `03_allelic_effects.R` | Reference-versus-alternative MPRAnalyze tests for eligible variant/construct groups | `allelic_effects/{CELL}_allelic_effects.tsv`, eligibility reports |
| `04_annotate_allelic_effects.R` | Annotate allele contrasts with enhancer activity in MEC-1, HG-3 or either | `integrated/{CELL}_allelic_effects_annotated.tsv` |
| `05_record_environment.R` | Record the installed R, Bioconductor and package versions | `docs/software_environment.txt` |

Here `{CELL}` is `MEC1` or `HG3`. Additional outputs include summary and exclusion tables; see [methods](docs/methods.md).

### Statistical settings retained from the original analysis

- **Enhancer activity:** `estimateDepthFactors(lib.factor = "replicate", which.lib = "both")`; `analyzeQuantification(dnaDesign = ~ replicate, rnaDesign = ~ 1)`; `testEmpirical()` without designated negative controls. An oligonucleotide is called active when BH-adjusted `pval.mad` ≤ 0.05.
- **Allelic effects:** `analyzeComparative(dnaDesign = ~ replicate, rnaDesign = ~ replicate + allele, reducedDesign = ~ replicate)` and `testLrt()`. A contrast is called significant when `fdr` ≤ 0.05 and `abs(logFC)` ≥ `log2(1.5)`.
- **Enhancer annotation:** A variant construct group is active in a cell line if at least one tested allele sequence is called active. A group not assessed in a cell line is distinguished from a tested but inactive group. Enhancer activity and significant allelic effects may occur in different cell lines.

`variant_core` is a specific oligonucleotide comparison context (e.g. `rs10005419_2`). `variant_id` is the annotated variant (e.g. `rs10005419`). Use `variant_id` to count unique rsIDs across cell lines, not the contrast identifier.

## Reproducibility and status

The analysis scripts record sequences excluded during QC and absent from MPRAnalyze results; allelic tests also report eligibility and unrecognized construct IDs. Validate **numeric results**, not only successful execution or total significant counts, as described in [validation](docs/validation.md).

The MPRAnalyze pipeline was run by the project author; this reorganized module has **not been independently re-executed in a clean environment**. The statistical steps of scripts 01–04 were retained; only their header instructions and surrounding documentation were revised for this folder structure. On the original analysis machine, run `05_record_environment.R` and review the resulting `docs/software_environment.txt` before public release.

Software citation, full project-level code/data availability, acknowledgements and licensing are handled in the [main repository README](../README.md). The module's [release checklist](docs/release_checklist.md) lists remaining validation and data-provenance checks.
