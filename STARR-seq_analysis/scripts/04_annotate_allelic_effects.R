#!/usr/bin/env Rscript
# Annotate allelic-effect results with enhancer activity observed in MEC-1 and HG-3.
# Run from the STARR-seq_analysis/ directory: Rscript scripts/04_annotate_allelic_effects.R
# This script does not recompute or filter the MPRAnalyze statistical results.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

annotation_file <- "data/sequence_annotations.tsv"
enhancer_dir <- "results/enhancer_activity"
allelic_dir <- "results/allelic_effects"
output_dir <- "results/integrated"
cell_lines <- c("MEC1", "HG3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use the same ID convention as 03_allelic_effects.R:
# rs10005419_2_ref_C and rs10005419_2_alt1_A -> rs10005419_2
# Recognize both _alt_ and _altN_ constructs, consistent with 03_allelic_effects.R.
extract_variant_core <- function(element_ids) {
  matches <- regmatches(element_ids, regexec("^(.*)_(ref|alt[0-9]*)_(.+)$", element_ids, perl = TRUE))
  vapply(matches, function(x) if (length(x) == 4) x[2] else NA_character_, character(1))
}

check_unique_ids <- function(data, id_column, description) {
  if (!id_column %in% names(data)) stop(description, ": missing column ", id_column)
  ids <- as.character(data[[id_column]])
  if (anyNA(ids) || any(trimws(ids) == "") || anyDuplicated(ids)) {
    stop(description, ": missing, blank or duplicated ", id_column, " values")
  }
}

parse_active_calls <- function(x, description) {
  values <- toupper(trimws(as.character(x)))
  parsed <- rep(NA, length(values))
  parsed[values %in% c("TRUE", "T", "1")] <- TRUE
  parsed[values %in% c("FALSE", "F", "0")] <- FALSE
  if (anyNA(parsed)) stop(description, ": active_call must be TRUE or FALSE, without missing values")
  parsed
}

if (!file.exists(annotation_file)) stop("Missing annotation file: ", annotation_file)
annotations <- read_tsv(annotation_file, col_types = cols(.default = col_character()), show_col_types = FALSE)
check_unique_ids(annotations, "element_id", "Annotations")

id_map <- annotations %>%
  transmute(element_id, variant_core = extract_variant_core(element_id)) %>%
  filter(!is.na(variant_core))
if (!nrow(id_map)) stop("No valid reference/alternative constructs found in annotations")

summarize_enhancer <- function(cell_line) {
  path <- file.path(enhancer_dir, paste0(cell_line, "_enhancer_activity.tsv"))
  if (!file.exists(path)) stop("Missing enhancer activity results: ", path)
  enhancer <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  check_unique_ids(enhancer, "element_id", paste(cell_line, "enhancer results"))
  if (!"active_call" %in% names(enhancer)) stop(cell_line, ": missing active_call column")
  enhancer$active_call <- parse_active_calls(enhancer$active_call, cell_line)
  enhancer <- enhancer %>% select(element_id, active_call)
  unexpected_ids <- setdiff(enhancer$element_id, annotations$element_id)
  if (length(unexpected_ids)) stop(cell_line, ": enhancer results contain IDs absent from sequence annotations")

  # An NA enhancer call indicates that none of the group's alleles was tested
  # (e.g. MPRAnalyze excluded them); it must not be interpreted as inactive.
  enhancer %>%
    inner_join(id_map, by = "element_id", relationship = "many-to-one") %>%
    group_by(variant_core) %>%
    summarise(n_tested = n(), n_active = sum(active_call), .groups = "drop") %>%
    rename(!!paste0("n_tested_", cell_line) := n_tested,
           !!paste0("n_active_", cell_line) := n_active)
}

mec1 <- summarize_enhancer("MEC1")
hg3 <- summarize_enhancer("HG3")
variant_summary <- id_map %>%
  distinct(variant_core) %>%
  left_join(mec1, by = "variant_core", relationship = "one-to-one") %>%
  left_join(hg3, by = "variant_core", relationship = "one-to-one") %>%
  mutate(
    enhancer_in_MEC1 = if_else(is.na(n_tested_MEC1), NA, n_active_MEC1 > 0),
    enhancer_in_HG3 = if_else(is.na(n_tested_HG3), NA, n_active_HG3 > 0),
    sig_enhancer = case_when(
      enhancer_in_MEC1 %in% TRUE | enhancer_in_HG3 %in% TRUE ~ TRUE,
      enhancer_in_MEC1 %in% FALSE & enhancer_in_HG3 %in% FALSE ~ FALSE,
      TRUE ~ NA
    ),
    enhancer_in = case_when(
      enhancer_in_MEC1 %in% TRUE & enhancer_in_HG3 %in% TRUE ~ "MEC-1 and HG-3",
      enhancer_in_MEC1 %in% TRUE & enhancer_in_HG3 %in% FALSE ~ "MEC-1",
      enhancer_in_MEC1 %in% FALSE & enhancer_in_HG3 %in% TRUE ~ "HG-3",
      enhancer_in_MEC1 %in% FALSE & enhancer_in_HG3 %in% FALSE ~ "None",
      enhancer_in_MEC1 %in% TRUE & is.na(enhancer_in_HG3) ~ "MEC-1 (HG-3 not assessed)",
      enhancer_in_HG3 %in% TRUE & is.na(enhancer_in_MEC1) ~ "HG-3 (MEC-1 not assessed)",
      TRUE ~ "Incomplete assessment"
    )
  ) %>%
  arrange(variant_core)
write_tsv(variant_summary, file.path(output_dir, "variant_enhancer_activity.tsv"))

annotate_cell_line <- function(cell_line) {
  path <- file.path(allelic_dir, paste0(cell_line, "_allelic_effects.tsv"))
  if (!file.exists(path)) stop("Missing allelic effect results: ", path)
  allelic_results <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  check_unique_ids(allelic_results, "contrast_id", paste(cell_line, "allelic results"))
  required <- c("variant_core", "allelic_effect_call")
  missing <- setdiff(required, names(allelic_results))
  if (length(missing)) stop(cell_line, ": missing allelic result columns: ", paste(missing, collapse = ", "))
  if (anyNA(allelic_results$variant_core) || any(!allelic_results$variant_core %in% variant_summary$variant_core)) {
    stop(cell_line, ": allelic results contain unknown variant groups")
  }
  allelic_results$allelic_effect_call <- parse_active_calls(allelic_results$allelic_effect_call, cell_line)
  annotated <- allelic_results %>%
    left_join(variant_summary, by = "variant_core", relationship = "many-to-one") %>%
    relocate(sig_enhancer, enhancer_in, enhancer_in_MEC1, enhancer_in_HG3,
             .after = variant_core)
  if (nrow(annotated) != nrow(allelic_results)) stop(cell_line, ": annotation changed the number of contrasts")
  output_path <- file.path(output_dir, paste0(cell_line, "_allelic_effects_annotated.tsv"))
  write_tsv(annotated, output_path)
  message(cell_line, ": annotated ", nrow(annotated), " allele contrasts; ",
          sum(annotated$allelic_effect_call), " significant allelic effects")
  tibble(cell_line = cell_line, tested_contrasts = nrow(annotated),
         significant_allelic_effects = sum(annotated$allelic_effect_call),
         enhancer_active_contrasts = sum(annotated$sig_enhancer %in% TRUE),
         significant_allelic_effects_in_active_enhancers =
           sum(annotated$allelic_effect_call & annotated$sig_enhancer %in% TRUE),
         contrasts_with_incomplete_enhancer_assessment = sum(is.na(annotated$sig_enhancer)))
}

integration_summary <- bind_rows(lapply(cell_lines, annotate_cell_line))
write_tsv(integration_summary, file.path(output_dir, "integration_summary.tsv"))
print(integration_summary)
