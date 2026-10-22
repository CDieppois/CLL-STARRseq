# Validation before public release

## 1. Rerun the pipeline from its true inputs

If starting at the main repository root, first run `cd STARR-seq_analysis`. The commands below and all relative paths in the R scripts assume that the current working directory is this module folder.

Check that the four files listed in [data/README.md](../data/README.md) exist and contain pre-QC counts. From the `STARR-seq_analysis/` module directory, execute:

```bash
Rscript scripts/01_quality_control.R
Rscript scripts/02_enhancer_activity.R
Rscript scripts/03_allelic_effects.R
Rscript scripts/04_annotate_allelic_effects.R
Rscript scripts/05_record_environment.R
```

A successful exit is necessary but **not sufficient** for reproducibility.

## 2. Verify the output lineage

Inspect `results/qc/qc_summary.tsv`, `results/enhancer_activity/enhancer_activity_summary.tsv`, `results/allelic_effects/allelic_effects_summary.tsv`, and `results/integrated/integration_summary.tsv` for both cell lines. Check the exclusion files and `unparsed_element_ids.tsv` for unexplained losses. In particular, verify the `*_alt_*` and `*_alt1_*` naming patterns in the eligibility report, including indels.

## 3. Count unique rsIDs and significant sequences

The following commands work in an R console with the working directory set to the `STARR-seq_analysis/` module directory:

```r
library(dplyr)
library(readr)

cell_lines <- c("MEC1", "HG3")
read_cell_results <- function(dir, suffix) {
  bind_rows(lapply(cell_lines, function(cell) {
    read_tsv(file.path("results", dir, paste0(cell, "_", suffix, ".tsv")),
             show_col_types = FALSE) %>% mutate(cell_line = cell)
  }))
}
enhancer <- read_cell_results("enhancer_activity", "enhancer_activity")
allelic <- read_cell_results("integrated", "allelic_effects_annotated")

annotation <- read_tsv("data/sequence_annotations.tsv", show_col_types = FALSE) %>%
  select(element_id, variant_id)
enhancer <- enhancer %>% left_join(annotation, by = "element_id", relationship = "many-to-one")

cat("Unique active enhancer sequences (either cell line): ",
    n_distinct(enhancer$element_id[enhancer$active_call %in% TRUE]), "\n")
cat("Unique rsIDs represented by active sequences: ",
    n_distinct(enhancer$variant_id[enhancer$active_call %in% TRUE], na.rm = TRUE), "\n")
cat("Significant allelic contrasts (cell-line tests combined): ",
    sum(allelic$allelic_effect_call %in% TRUE), "\n")
cat("Unique rsIDs with significant allelic effects (either cell line): ",
    n_distinct(allelic$variant_id[allelic$allelic_effect_call %in% TRUE], na.rm = TRUE), "\n")
cat("Unique rsIDs with allelic effects AND enhancer activity in either cell line: ",
    n_distinct(allelic$variant_id[allelic$allelic_effect_call %in% TRUE &
      allelic$sig_enhancer %in% TRUE], na.rm = TRUE), "\n")
```

An rsID can be observed in both cell lines and multiple contrast groups. **Do not sum cell-line-specific distinct-rsID totals** to obtain the unique combined total.

## 4. Compare against the original analysis

On the original and new results, compare per cell line:

- Number of pre-QC oligonucleotides, retained elements, and elements excluded by each QC criterion.
- Number of enhancer activity tests and active sequences.
- Per-sequence `alpha`, `pval.mad` and adjusted `padj` values, joining on `element_id`.
- Number of prepared/returned ref-vs-alt contrasts and per-contrast `logFC`, `fdr`, and significance call, joining on `contrast_id`.
- Unique significant rsIDs and unique rsIDs with both properties across MEC-1 and HG-3.

Record any intentional difference from earlier exploratory Excel results and explain the cause. Counts alone can match while individual values differ; compare numeric outputs, not only summaries.

## 5. Reproduce in a clean environment

After confirming the results in your currently working R installation, record exact R and package versions with `05_record_environment.R`. For more rigorous environment reproducibility, use `renv` to create a lockfile from that actual installation, then test `renv::restore()` in a new environment. Do not replace a validated environment with the latest package versions without checking the outputs again.

Until these comparisons and package-version checks have been completed, the repository is **runnable but not independently validated**.
