# Methods implemented in the analysis scripts

This document describes the STARR-seq code as currently implemented. All relative paths below are resolved from the `STARR-seq_analysis/` directory. It should be checked against the final manuscript's Methods section before publication.

## 01. Quality control

`01_quality_control.R` reads the MEC-1 and HG-3 count matrices and joins each to the shared input library on `element_id`. It stops on missing input counts, duplicate IDs, absent expected columns, or nonnegative-integer count violations. It retains an oligonucleotide only if its shared input-library count is strictly positive and all four DNA replicate counts are >= 5. The RNA counts are **not** subject to a minimum filter.

Outputs in `results/qc/`:

- `MEC1_counts_filtered.tsv`, `HG3_counts_filtered.tsv`: retained oligonucleotides with DNA and RNA columns; the joined `input` field is not retained in these matrices.
- `MEC1_filter_report.tsv`, `HG3_filter_report.tsv`: input count, DNA counts, exclusion flags and reason for each oligonucleotide.
- `qc_summary.tsv`: total, excluded by zero input, excluded by low DNA, overlap, and retained elements. Counts for failure categories **overlap** by design.

## 02. Enhancer activity

`02_enhancer_activity.R` analyzes the filtered sequences independently in each cell line. It constructs paired DNA and RNA count matrices with replicate annotations (four biological/technical replicates as supplied), estimates sequencing depth factors with `estimateDepthFactors(..., lib.factor = "replicate", which.lib = "both")`, then fits a quantification model with `dnaDesign = ~ replicate` and `rnaDesign = ~ 1`.

After `testEmpirical()` it uses **`pval.mad`**, not negative-control empirical p-values, and adjusts these p-values using Benjamini–Hochberg. `active_call = TRUE` denotes adjusted p-value <= 0.05. An unadjusted `nominal_activity_call` is also provided. The activity estimate `alpha` is extracted with `getAlpha()` and matched by `element_id`.

Outputs in `results/enhancer_activity/`:

- `{CELL}_enhancer_activity.tsv`: per-element statistical results and calls.
- `{CELL}_mpranalyze_excluded.tsv`: counts for elements missing from final MPRAnalyze test results, with diagnostics on zero counts and whether IDs remain in the object. These flags are **diagnostics**, not authoritative package-provided causes of exclusion.
- `enhancer_activity_summary.tsv`: input/tested/active/excluded counts.

## 03. Allelic effects

`03_allelic_effects.R` matches the filtered counts to `sequence_annotations.tsv`. It extracts a construct-group identifier (`variant_core`) and a ref/alt label from `element_id` using the pattern `^(.*)_(ref|alt[0-9]*)_(.+)$`, which accepts `*_alt_*` as well as `*_alt1_*`, `*_alt2_*`, etc. The same rsID can occur in multiple construct groups.

For each construct group, testing requires exactly one annotated reference and at least one alternative construct surviving QC, no duplicate allele labels, and a consistent `variant_id`. Missing/ineligible groups are reported and do not enter the model. Multiallelic groups yield separate reference-versus-alternative contrasts.

For each contrast, the same four replicates yield eight sample columns in paired ref/alt order. MPRAnalyze runs `estimateDepthFactors(..., lib.factor = "replicate", which.lib = "both")`, then `analyzeComparative(dnaDesign = ~ replicate, rnaDesign = ~ replicate + allele, reducedDesign = ~ replicate)`, followed by `testLrt()`. Significant allelic effects satisfy `fdr <= 0.05` **and** `abs(logFC) >= log2(1.5)`.

Outputs in `results/allelic_effects/`:

- `{CELL}_allelic_effects.tsv`: all tested allele contrasts, effect sizes and statistical calls.
- `{CELL}_prepared_contrasts.tsv`: mapping between reference and alternate oligonucleotide IDs for every prepared comparison.
- `{CELL}_variant_eligibility.tsv`: status/reason of each annotated construct group after QC.
- `{CELL}_mpranalyze_excluded.tsv`: prepared contrasts lacking a final LRT result.
- `unparsed_element_ids.tsv`: annotation IDs that do not match the ref/alt naming convention.
- `allelic_effects_summary.tsv`: prepared, tested and significant contrasts per cell line.

## 04. Integration of allelic effects and enhancer activity

`04_annotate_allelic_effects.R` merges tested enhancer sequences with the allele-group identifier and records a group as active in a cell line when at least one tested allele had `active_call = TRUE`. If no allele in that group was tested in that cell line, the call is **not assessed**, represented by `NA`, rather than FALSE.

The integrated allelic-effect results retain both cell-line enhancer statuses (`enhancer_in_MEC1`, `enhancer_in_HG3`), a combined annotation (`sig_enhancer`) and a readable `enhancer_in` category. These annotations **do not recompute or alter** the allelic-effect test statistics.

Outputs in `results/integrated/`:

- `{CELL}_allelic_effects_annotated.tsv`: allelic effects plus cross-cell-line enhancer annotations.
- `variant_enhancer_activity.tsv`: per-construct-group enhancer tests and calls.
- `integration_summary.tsv`: counts of contrasts by effect and enhancer status.

## Interpretation and limitations

- Enhancer activity in *either* cell line and a significant allelic effect in *either* cell line are not required to occur in the same cell line.
- `variant_core` identifies the specific comparison group. `variant_id` identifies the underlying annotated variant/rsID. Use distinct rsIDs for the overall unique-variant count.
- The reported analysis begins with processed counts; it does not reproduce primary read alignment or oligonucleotide counting.
- A successful run is not alone proof that the output reproduces the original manuscript. See [validation.md](validation.md).
