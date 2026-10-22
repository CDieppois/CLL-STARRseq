# QTL input formats and source file mapping

The five QTL scripts begin with an **analysis-ready genotype file** and published molecular matrices. None of the source data are committed to this repository. `01_prepare_phenotypes.R` can read the original filenames below directly, or previously standardized phenotype TSVs if available locally.

## Variant/genotype input (local only)

`QTL_variant_genotypes.tsv`: tab-separated; one unique `variant_id` per row, metadata and patient genotype columns. Important columns:

| Name | Meaning |
|---|---|
| `variant_id` | rsID or unique association-variant identifier |
| `genotype_ref`, `genotype_alt` | Genotype encoding reference allele and comma-separated ALT alleles |
| `tested_ref`, `tested_alt` | STARR-seq-tested alleles; dosage refers to `tested_alt` |
| `tad_chromosome`, `tad_start_hg38`, `tad_end_hg38` | hg38 TAD interval for mQTL, caQTL and hQTL |
| `tad_chromosome_hg19`, `tad_start_hg19`, `tad_end_hg19` | **Pre-lifted hg19 TAD interval for eQTL** |
| `CLL1`, `CLL2`, ... | Per-sample diploid genotype calls (`0/0`, `0/1`, `1/1`, `0/2`, `2/2`, missing `./.`) |

Additional annotations such as `chromosome`, `start_hg38`, `end_hg38`, `tad_id`, `variant_class` or `category` are allowed. Do **not** include real genotype calls in the public example.

This table is already allele-resolved, genotype-merged and TAD-annotated. **Genotype extraction/imputation and liftOver are not implemented by this pipeline.** Validate the hg19 TAD coordinates against your original liftOver results.

## Molecular matrices: original filenames and formats

| QTL | Original publication filename (local) | Important source fields | Coordinate build |
|---|---|---|---|
| eQTL | `single_probes_all_cases_RMA_expression.txt` | `Probe_Set_ID`, `CLL#` GC-RMA values | hg19 alignment from GPL13667 |
| eQTL annotation | `GPL13667-15572.txt` | `ID`, `Gene Symbol`, `Alignments` | hg19 probe alignment in original notebooks |
| caQTL | `atac_consensus_sample_count_corrected.txt` | `chr`, `start`, `end`, `CLL#` accessibility signals | hg38 |
| hQTL | `H3K27ac_consensus_sample_count_corrected.txt` | `chr`, `start`, `end`, `CLL#` peak signals | hg38 |
| mQTL | `cll_methylation_full_matrix_idibaps_norm_Beta_grch38.tsv.gz` | `chr`, `start`, `end`, `name`, numeric patient columns (`1`, `2`, ...), beta values | hg38 |

The code matches patient IDs with columns named `CLL1`, `CLL2`, etc. It renames numeric methylation patient columns to `CLL#`. The gzip methylation input requires `R.utils` for `data.table::fread` to decompress it.

For gene expression, the code reads GPL13667, matches `Probe_Set_ID` to GEO `ID`, excludes probes with missing/invalid gene symbols, and extracts the first alignment position. These are source-specific transformations from the original notebook, **not** general Affymetrix annotation methods.

The script prioritizes these already-standardized local TSVs, if present in the same `QTL_INPUT_DIR`: `gene_expression.tsv`, `DNA_methylation.tsv`, `ATACseq.tsv`, `H3K27ac.tsv`. Each standardized matrix must use `feature_id`, `chromosome`, `start`, `end`, optional `gene_symbol` for eQTL, and samples `CLL#`. If these files are absent, the script attempts the original publication filenames above.

**Do not mix genome builds.** Original eQTL probes and their TADs use hg19; all other assays and corresponding TADs use hg38. Check the paper and original input coordinates if underlying datasets change.

## Local directory setup

Place `QTL_variant_genotypes.tsv` alongside the five published input files in a local directory, e.g. `D:/CLL_QTL_data`. This local directory is not tracked by GitHub.

On RStudio:

```r
Sys.setenv(QTL_INPUT_DIR = 'D:/CLL_QTL_data')
source('scripts/01_prepare_phenotypes.R')
```

On Linux shell:

```bash
export QTL_INPUT_DIR=/path/to/CLL_QTL_data
Rscript scripts/01_prepare_phenotypes.R eQTL
```

From stage 02 onward, pipeline-generated intermediates are read from `results/`.

See [data_access.md](data_access.md) for provenance and links, and the synthetic examples under `example_data/` for column formats.
