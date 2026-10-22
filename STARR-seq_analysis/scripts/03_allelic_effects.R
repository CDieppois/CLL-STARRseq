#!/usr/bin/env Rscript
# MPRAnalyze allelic-effect analysis for MEC-1 and HG-3.
# Run from the STARR-seq_analysis/ directory: Rscript scripts/03_allelic_effects.R

suppressPackageStartupMessages({
  library(MPRAnalyze)
  library(dplyr)
  library(readr)
  library(tibble)
})

input_dir <- "results/qc"
annotation_file <- "data/sequence_annotations.tsv"
output_dir <- "results/allelic_effects"
cell_lines <- c("MEC1", "HG3")
dna_columns <- paste0("DNA_Rep", 1:4)
rna_columns <- paste0("RNA_Rep", 1:4)
fdr_threshold <- 0.05
min_fold_change <- 1.5

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The suffix of each element ID identifies its allele. The part preceding
# _ref_, _alt_, or _altN_ identifies the specific oligonucleotide/variant context.
# Example: rs10005419_2_ref_C -> variant_core=rs10005419_2, allele=ref.
# Accept both unnumbered alt IDs (e.g. _alt_len2bp) and numbered alt IDs (_alt1_A).
extract_allele_metadata <- function(element_ids) {
  matches <- regmatches(element_ids, regexec("^(.*)_(ref|alt[0-9]*)_(.+)$", element_ids, perl = TRUE))
  tibble(
    element_id = element_ids,
    variant_core = vapply(matches, function(x) if (length(x) == 4) x[2] else NA_character_, character(1)),
    allele_label = vapply(matches, function(x) if (length(x) == 4) x[3] else NA_character_, character(1))
  )
}

if (!file.exists(annotation_file)) stop("Missing annotation file: ", annotation_file)
annotations <- read_tsv(annotation_file, col_types = cols(.default = col_character()), show_col_types = FALSE)
required_annotations <- c("element_id", "variant_id", "tested_allele")
missing_annotations <- setdiff(required_annotations, names(annotations))
if (length(missing_annotations)) stop("Missing annotation columns: ", paste(missing_annotations, collapse = ", "))
annotations$element_id <- trimws(annotations$element_id)
if (anyNA(annotations$element_id) || any(annotations$element_id == "") || anyDuplicated(annotations$element_id)) {
  stop("sequence_annotations.tsv contains missing or duplicate element_id values")
}

variant_annotations <- annotations %>%
  select(all_of(required_annotations)) %>%
  left_join(extract_allele_metadata(annotations$element_id), by = "element_id", relationship = "one-to-one") %>%
  filter(!is.na(variant_core))
if (!nrow(variant_annotations)) stop("No reference/alternative IDs recognized in the annotation file")
# Unrecognized IDs may include controls or non-variant sequences; list them
# explicitly so that variant constructs with unexpected names are not overlooked.
unparsed_annotations <- annotations %>% filter(!element_id %in% variant_annotations$element_id)
write_tsv(unparsed_annotations, file.path(output_dir, "unparsed_element_ids.tsv"))
message("Annotations without a recognized ref/alt identifier: ", nrow(unparsed_annotations),
        " (see unparsed_element_ids.tsv)")
if (anyNA(variant_annotations$variant_id) || any(trimws(variant_annotations$variant_id) == "") ||
    anyNA(variant_annotations$tested_allele) || any(trimws(variant_annotations$tested_allele) == "")) {
  stop("Missing variant_id or tested_allele for one or more allele constructs")
}

analyze_cell_line <- function(cell_line) {
  input_file <- file.path(input_dir, paste0(cell_line, "_counts_filtered.tsv"))
  if (!file.exists(input_file)) stop("Missing QC output: ", input_file)
  counts <- read_tsv(input_file, show_col_types = FALSE, progress = FALSE)
  required_counts <- c("element_id", dna_columns, rna_columns)
  missing_counts <- setdiff(required_counts, names(counts))
  if (length(missing_counts)) stop(cell_line, ": missing count columns: ", paste(missing_counts, collapse = ", "))
  counts$element_id <- as.character(counts$element_id)
  if (anyNA(counts$element_id) || any(counts$element_id == "") || anyDuplicated(counts$element_id)) {
    stop(cell_line, ": missing or duplicate element IDs")
  }
  for (column in c(dna_columns, rna_columns)) {
    values <- counts[[column]]
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
        any(values < 0) || any(values != floor(values))) {
      stop(cell_line, ": invalid nonnegative integer counts in ", column)
    }
  }
  unannotated <- setdiff(counts$element_id, annotations$element_id)
  if (length(unannotated)) stop(cell_line, ": ", length(unannotated), " count elements have no annotations (e.g. ",
                                 paste(head(unannotated, 5), collapse = ", "), ")")

  # Eligibility is assessed independently in each cell line after QC.
  # We report missing ref/alt alleles rather than duplicating measurements.
  annotated <- variant_annotations %>% mutate(passed_qc = element_id %in% counts$element_id)
  variant_status <- annotated %>%
    group_by(variant_core) %>%
    summarise(
      variant_id = first(variant_id),
      annotated_alleles = n(),
      qc_alleles = sum(passed_qc),
      reference_in_annotation = sum(allele_label == "ref"),
      reference_after_qc = sum(allele_label == "ref" & passed_qc),
      alternatives_after_qc = sum(allele_label != "ref" & passed_qc),
      duplicated_allele_labels = anyDuplicated(allele_label) > 0,
      inconsistent_variant_ids = n_distinct(variant_id) != 1,
      .groups = "drop"
    ) %>%
    mutate(status = case_when(
      inconsistent_variant_ids ~ "inconsistent_variant_ids",
      reference_in_annotation != 1 ~ "invalid_reference_annotation",
      duplicated_allele_labels ~ "duplicated_allele_labels",
      reference_after_qc == 0 ~ "reference_failed_qc",
      alternatives_after_qc == 0 ~ "no_alternative_after_qc",
      TRUE ~ "eligible"
    )) %>%
    arrange(variant_core)
  write_tsv(variant_status, file.path(output_dir, paste0(cell_line, "_variant_eligibility.tsv")))

  eligible_cores <- variant_status$variant_core[variant_status$status == "eligible"]
  if (!length(eligible_cores)) stop(cell_line, ": no valid reference/alternative comparisons after QC")
  eligible <- annotated %>% filter(variant_core %in% eligible_cores, passed_qc)
  reference <- eligible %>%
    filter(allele_label == "ref") %>%
    transmute(variant_core, variant_id, reference_element_id = element_id,
              reference_allele = tested_allele)
  alternatives <- eligible %>%
    filter(allele_label != "ref") %>%
    transmute(variant_core, alternative_label = allele_label,
              alternative_element_id = element_id, alternative_allele = tested_allele)
  contrasts <- alternatives %>%
    inner_join(reference, by = "variant_core", relationship = "many-to-one") %>%
    mutate(contrast_id = paste0(variant_core, "_ref_vs_", alternative_label)) %>%
    select(contrast_id, variant_core, variant_id, alternative_label,
           reference_element_id, alternative_element_id, reference_allele, alternative_allele) %>%
    arrange(contrast_id)
  if (!nrow(contrasts) || anyDuplicated(contrasts$contrast_id)) {
    stop(cell_line, ": empty or duplicated allele contrasts")
  }
  write_tsv(contrasts, file.path(output_dir, paste0(cell_line, "_prepared_contrasts.tsv")))

  # One MPRAnalyze row per ref-vs-alt contrast; each row has 8 libraries
  # (rep1_ref, rep1_alt, ..., rep4_ref, rep4_alt), preserving original design.
  reference_indices <- match(contrasts$reference_element_id, counts$element_id)
  alternative_indices <- match(contrasts$alternative_element_id, counts$element_id)
  if (anyNA(reference_indices) || anyNA(alternative_indices)) stop(cell_line, ": unmatched allele count IDs")
  replicate_numbers <- seq_along(dna_columns)
  sample_names <- as.vector(rbind(paste0("rep", replicate_numbers, "_ref"),
                                   paste0("rep", replicate_numbers, "_alt")))
  dna_counts <- matrix(NA_real_, nrow = nrow(contrasts), ncol = length(sample_names),
                       dimnames = list(contrasts$contrast_id, sample_names))
  rna_counts <- dna_counts
  for (i in replicate_numbers) {
    dna_counts[, 2 * i - 1] <- counts[[dna_columns[i]]][reference_indices]
    dna_counts[, 2 * i] <- counts[[dna_columns[i]]][alternative_indices]
    rna_counts[, 2 * i - 1] <- counts[[rna_columns[i]]][reference_indices]
    rna_counts[, 2 * i] <- counts[[rna_columns[i]]][alternative_indices]
  }
  sample_annotation <- data.frame(
    replicate = factor(rep(paste0("rep", replicate_numbers), each = 2)),
    allele = factor(rep(c("ref", "alt"), length(replicate_numbers)), levels = c("ref", "alt")),
    row.names = sample_names
  )

  mpra_object <- MpraObject(dnaCounts = dna_counts, rnaCounts = rna_counts,
                            colAnnot = sample_annotation)
  mpra_object <- estimateDepthFactors(mpra_object, lib.factor = "replicate", which.lib = "both")
  mpra_object <- analyzeComparative(mpra_object, dnaDesign = ~ replicate,
                                    rnaDesign = ~ replicate + allele,
                                    reducedDesign = ~ replicate)
  test_results <- as.data.frame(testLrt(mpra_object))
  test_results$contrast_id <- rownames(test_results)
  if (!all(c("fdr", "logFC") %in% names(test_results))) {
    stop(cell_line, ": MPRAnalyze output lacks fdr or logFC")
  }
  if (anyNA(test_results$contrast_id) || anyDuplicated(test_results$contrast_id) ||
      any(!test_results$contrast_id %in% contrasts$contrast_id)) {
    stop(cell_line, ": unexpected or duplicated IDs returned by testLrt()")
  }

  allelic_results <- contrasts %>%
    inner_join(as_tibble(test_results), by = "contrast_id", relationship = "one-to-one") %>%
    mutate(
      fold_change_alt_vs_ref = 2^logFC,
      passes_fdr = !is.na(fdr) & fdr <= fdr_threshold,
      passes_fold_change = !is.na(logFC) & abs(logFC) >= log2(min_fold_change),
      allelic_effect_call = passes_fdr & passes_fold_change
    )

  # Diagnostic list: allele contrasts that were prepared but lack an LRT result.
  fitted_ids <- rownames(dnaCounts(mpra_object))
  excluded <- contrasts %>%
    filter(!contrast_id %in% test_results$contrast_id) %>%
    mutate(
      excluded_from_object = !contrast_id %in% fitted_ids,
      all_DNA_zero = rowSums(dna_counts[contrast_id, , drop = FALSE]) == 0,
      all_RNA_zero = rowSums(rna_counts[contrast_id, , drop = FALSE]) == 0
    )
  write_tsv(excluded, file.path(output_dir, paste0(cell_line, "_mpranalyze_excluded.tsv")))
  write_tsv(allelic_results, file.path(output_dir, paste0(cell_line, "_allelic_effects.tsv")))

  message(cell_line, ": ", nrow(contrasts), " prepared, ", nrow(allelic_results),
          " tested, ", sum(allelic_results$allelic_effect_call), " significant (FDR <= ",
          fdr_threshold, ", |FC| >= ", min_fold_change, ")")
  tibble(cell_line = cell_line, annotated_variant_groups = nrow(variant_status),
         eligible_variant_groups = length(eligible_cores), prepared_contrasts = nrow(contrasts),
         tested_contrasts = nrow(allelic_results),
         mpranalyze_excluded = nrow(excluded),
         significant_allelic_effects = sum(allelic_results$allelic_effect_call))
}

summary <- bind_rows(lapply(cell_lines, analyze_cell_line))
write_tsv(summary, file.path(output_dir, "allelic_effects_summary.tsv"))
print(summary)
