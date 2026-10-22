# QTL analysis

This module contains R scripts for expression QTL (eQTL), DNA methylation QTL (mQTL), H3K27ac QTL (hQTL), and chromatin accessibility QTL (caQTL) analyses supporting the CLL STARR-seq study. It is one self-contained analysis folder within the larger `CLL-STARRseq` repository.

The workflow starts from an **analysis-ready, allele-resolved variant/genotype table**, previously processed molecular phenotype matrices, and TAD assignments. It does **not** perform genotype extraction, imputation, upstream data acquisition, or figure generation. It stops at annotated association results.

**Validation status:** This is a refactoring draft of the original QTL R Markdown notebooks. It has *not* been rerun on the controlled-access datasets or verified numerically against the original associations. Do not describe it as a validated reproducibility package until those comparisons have been completed.

## Directory structure

```text
QTL_analysis/
├── README.md
├── .gitignore
├── scripts/
│   ├── 01_prepare_phenotypes.R
│   ├── 02_run_PEER.R
│   ├── 03_prepare_MatrixEQTL_inputs.R
│   ├── 04_run_MatrixEQTL.R
│   ├── 05_annotate_QTL_results.R
│   └── utils.R
├── docs/
│   ├── data_access.md
│   ├── input_formats.md
│   └── validation.md
├── example_data/
│   ├── QTL_variant_genotypes_example.tsv
│   ├── gene_expression_example.tsv
│   ├── DNA_methylation_example.tsv
│   ├── H3K27ac_example.tsv
│   └── ATACseq_example.tsv
└── results/                       # Created by scripts; excluded from Git
```

Only **synthetic format examples** are included. The real sample-level genotype data and published molecular-phenotype source matrices are **not** included. See [data sources and access](docs/data_access.md) and [input specifications](docs/input_formats.md).

## Data requirements

Prepare `QTL_variant_genotypes.tsv` locally with unique `variant_id`, genome genotype `genotype_ref`/`genotype_alt`, the `tested_ref`/`tested_alt` alleles, TAD coordinates in hg38 *and* the hg19 TAD coordinates used for gene-expression probes, and diploid sample genotype columns (`CLL1`, `CLL2`, etc.). The real data are controlled access; the [synthetic example](example_data/QTL_variant_genotypes_example.tsv) is **only an illustration**, not a substitute for authorized inputs.

Molecular phenotype files referenced in the original analyses:

| Analysis | Source molecular matrix | PEER factors |
| --- | --- | ---: |
| eQTL | `single_probes_all_cases_RMA_expression.txt` and platform annotation `GPL13667-15572.txt` | 65 |
| mQTL | `cll_methylation_full_matrix_idibaps_norm_Beta_grch38.tsv.gz` | 45 |
| hQTL | `H3K27ac_consensus_sample_count_corrected.txt` | 10 |
| caQTL | `atac_consensus_sample_count_corrected.txt` | 10 |

The corresponding molecular datasets originate from [Speedy et al., *Nature Communications* (2019)](https://doi.org/10.1038/s41467-019-11582-2). The genotype inputs are derived from controlled-access studies deposited at EGA. Detailed links, original accessions, and data reuse restrictions are in [data_access.md](docs/data_access.md). The scripts **do not download the data**.

The scripts can also accept already-standardized local phenotype files, for example `gene_expression.tsv` or `ATACseq.tsv`, in the same input directory. See [input formats](docs/input_formats.md).

## Requirements

An R environment with `data.table`, `MatrixEQTL`, and the R `peer` package is required. The published methylation `.gz` matrix may additionally require `R.utils` with `data.table::fread()`.

```r
install.packages(c("data.table", "MatrixEQTL", "R.utils"))
```

Install a compatible version of `peer` separately for your operating system (the CRAN command above does not install it). Record R, PEER and MatrixEQTL versions from the original analysis environment; these are not yet version-locked in this module.

## Local input directory

Keep authorized input datasets **outside this GitHub repository** and point `QTL_INPUT_DIR` to the directory containing `QTL_variant_genotypes.tsv` and your molecular matrices. There is no `restricted/` folder within the repository.

Windows (RStudio session):

```r
Sys.setenv(QTL_INPUT_DIR = "D:/CLL_QTL_data")  # Change to your actual local folder
```

Linux terminal:

```bash
export QTL_INPUT_DIR=/path/to/CLL_QTL_data
```

`QTL_INPUT_DIR` must be set for each new R session or shell that runs a script. Do **not** upload the real input folder, PEER outputs or intermediate genotype matrices.

## Run the pipeline

**Run these commands from inside `QTL_analysis/`, not from the main `CLL-STARRseq` repository root.** For example, from the main repository directory:

```bash
cd QTL_analysis
Rscript scripts/01_prepare_phenotypes.R all
Rscript scripts/02_run_PEER.R all
Rscript scripts/03_prepare_MatrixEQTL_inputs.R all
Rscript scripts/04_run_MatrixEQTL.R all
Rscript scripts/05_annotate_QTL_results.R all
```

Replace `all` with `eQTL`, `mQTL`, `hQTL`, or `caQTL` to run one phenotype. From RStudio, set the working directory to `QTL_analysis/` and source the script(s) in the same order. The scripts read `scripts/utils.R` relative to that folder, and generate `results/<QTL type>/` there.

## Analysis steps

| Step | Purpose | Principal output |
| --- | --- | --- |
| `01_prepare_phenotypes.R` | Match samples, select TAD-associated features, apply phenotype-specific filtering, create PEER input | `selected_variant_feature_pairs.tsv`, `PEER_matrix_<phenotype>.txt` |
| `02_run_PEER.R` | Estimate latent factors; save residuals and inverse-normal transformed residuals | `PEER_<phenotype>_K<N>_residuals_INT.tsv` |
| `03_prepare_MatrixEQTL_inputs.R` | Orient dosage to the tested alternate allele; align sample order; prepare selected variant–feature pairs | `MatrixEQTL_genotype_dosage.tsv`, `MatrixEQTL_phenotypes.tsv`, `MatrixEQTL_selected_pairs.tsv` |
| `04_run_MatrixEQTL.R` | Run additive linear association testing | `MatrixEQTL_all_results.tsv` |
| `05_annotate_QTL_results.R` | Restrict to the selected TAD-related pairs; compute BH FDR for those pairs; attach variant/feature annotations | `<QTL type>_annotated_results.tsv`, `<QTL type>_significant_results.tsv` |

### Statistical conventions

- PEER factors: **65 gene expression, 45 DNA methylation, 10 H3K27ac, 10 ATAC-seq**.
- Gene expression GC-RMA values: threshold > 4.5 in enough patients within at least one genotype category; ATAC-seq and H3K27ac signals: > 0 with the same group-coverage criterion. Methylation uses beta-to-M-value conversion.
- Genomic selection uses hg19 expression-probe/TAD coordinates and hg38 methylation/H3K27ac/ATAC/TAD coordinates; the provided analysis-ready genotype table must supply compatible TAD intervals.
- Genotypes are encoded as dosage of the **STARR-seq tested alternative allele**, not necessarily genome ALT.
- Matrix eQTL applies the additive linear model. BH FDR is recomputed over the selected variant–feature pairs after Matrix eQTL testing, matching the original analytical design.

See [input formats](docs/input_formats.md) and [validation](docs/validation.md) for detailed definitions and checks. PEER stochastic initialization/seed, sample matching, effect signs and significance must be numerically compared with original outputs before public claims of equivalence.

## Reuse and availability

This module contains no redistributed EGA genotypes or original processed patient-level molecular matrices. The software license, overall project citation and publication status will be documented in the [main project README](../README.md) when the full repository is ready.
