#!/usr/bin/env Rscript
# MPRAnalyze enhancer activity analysis for MEC-1 and HG-3.
# Run from the STARR-seq_analysis/ directory: Rscript scripts/02_enhancer_activity.R

suppressPackageStartupMessages({
  library(MPRAnalyze)
  library(dplyr)
  library(readr)
  library(tibble)
})

input_dir <- "results/qc"
output_dir <- "results/enhancer_activity"
cell_lines <- c("MEC1", "HG3")
dna_columns <- paste0("DNA_Rep", 1:4)
rna_columns <- paste0("RNA_Rep", 1:4)
fdr_threshold <- 0.05

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analyze_cell_line <- function(cell_line) {
  input_file <- file.path(input_dir, paste0(cell_line, "_counts_filtered.tsv"))
  if (!file.exists(input_file)) stop("Missing QC output: ", input_file)

  counts <- read_tsv(input_file, show_col_types = FALSE, progress = FALSE)
  required_columns <- c("element_id", dna_columns, rna_columns)
  missing_columns <- setdiff(required_columns, names(counts))
  if (length(missing_columns)) stop(cell_line, ": missing columns: ", paste(missing_columns, collapse = ", "))
  counts$element_id <- as.character(counts$element_id)
  if (anyNA(counts$element_id) || any(counts$element_id == "") || anyDuplicated(counts$element_id)) {
    stop(cell_line, ": invalid or duplicated element IDs")
  }
  for (column in c(dna_columns, rna_columns)) {
    values <- counts[[column]]
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) || any(values < 0) || any(values != floor(values))) {
      stop(cell_line, ": invalid nonnegative integer counts in ", column)
    }
  }
  if (!nrow(counts)) stop(cell_line, ": no sequences remaining after QC")

  dna_counts <- as.matrix(counts[, dna_columns])
  rna_counts <- as.matrix(counts[, rna_columns])
  rownames(dna_counts) <- counts$element_id
  rownames(rna_counts) <- counts$element_id
  dna_annotation <- data.frame(condition = "NS", replicate = factor(1:4), row.names = dna_columns)
  rna_annotation <- data.frame(condition = "NS", replicate = factor(1:4), row.names = rna_columns)
  element_annotation <- data.frame(element_id = counts$element_id, row.names = counts$element_id)

  mpra_object <- MpraObject(dnaCounts = dna_counts, rnaCounts = rna_counts,
                            dnaAnnot = dna_annotation, rnaAnnot = rna_annotation,
                            rowAnnot = element_annotation)
  mpra_object <- estimateDepthFactors(mpra_object, lib.factor = "replicate", which.lib = "both")
  mpra_object <- analyzeQuantification(mpra_object, dnaDesign = ~ replicate, rnaDesign = ~ 1)

  activity_tests <- as.data.frame(testEmpirical(mpra_object))
  activity_tests$element_id <- rownames(activity_tests)
  activity_estimates <- as.data.frame(getAlpha(mpra_object, by.factor = NULL))
  activity_estimates$element_id <- rownames(activity_estimates)
  if (!"alpha" %in% names(activity_estimates)) stop(cell_line, ": getAlpha() did not return an alpha column")
  if (!"pval.mad" %in% names(activity_tests)) stop(cell_line, ": testEmpirical() did not return pval.mad")
  if (anyDuplicated(activity_tests$element_id) || anyDuplicated(activity_estimates$element_id)) {
    stop(cell_line, ": duplicated IDs in MPRAnalyze results")
  }
  if (!setequal(activity_tests$element_id, activity_estimates$element_id)) {
    stop(cell_line, ": activity tests and alpha estimates contain different element IDs")
  }
  enhancer_results <- activity_tests %>%
    as_tibble() %>%
    left_join(as_tibble(activity_estimates[, c("element_id", "alpha")]), by = "element_id", relationship = "one-to-one") %>%
    mutate(padj = p.adjust(pval.mad, method = "BH"),
           nominal_activity_call = !is.na(pval.mad) & pval.mad <= 0.05,
           active_call = !is.na(padj) & padj <= fdr_threshold) %>%
    relocate(element_id, alpha, pval.mad, padj, nominal_activity_call, active_call)
  if (anyNA(enhancer_results$alpha)) stop(cell_line, ": missing alpha estimates after matching by ID")

  retained_ids <- rownames(dnaCounts(mpra_object))
  tested_ids <- enhancer_results$element_id
  excluded_ids <- setdiff(counts$element_id, tested_ids)
  # Report both IDs absent from the fitted object and those absent from final test results.
  excluded_sequences <- counts %>%
    filter(element_id %in% excluded_ids) %>%
    mutate(all_DNA_zero = rowSums(across(all_of(dna_columns))) == 0,
           all_RNA_zero = rowSums(across(all_of(rna_columns))) == 0,
           excluded_from_object = !(element_id %in% retained_ids),
           excluded_from_results = !(element_id %in% tested_ids))
  if (any(!tested_ids %in% counts$element_id)) stop(cell_line, ": test results contain unknown element IDs")
  write_tsv(enhancer_results, file.path(output_dir, paste0(cell_line, "_enhancer_activity.tsv")))
  write_tsv(excluded_sequences, file.path(output_dir, paste0(cell_line, "_mpranalyze_excluded.tsv")))
  message(cell_line, ": ", nrow(enhancer_results), " tested, ", sum(enhancer_results$active_call),
          " active (FDR <= ", fdr_threshold, "), ", length(excluded_ids), " excluded internally")
  tibble(cell_line = cell_line, input_elements = nrow(counts), tested_elements = nrow(enhancer_results),
         active_elements = sum(enhancer_results$active_call), mpranalyze_excluded = length(excluded_ids))
}

summary <- bind_rows(lapply(cell_lines, analyze_cell_line))
write_tsv(summary, file.path(output_dir, "enhancer_activity_summary.tsv"))
print(summary)
