#!/usr/bin/env Rscript
# Quality control of STARR-seq oligonucleotide count matrices.
# Run from the STARR-seq_analysis/ directory: Rscript scripts/01_quality_control.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

input_dir <- "data"
output_dir <- "results/qc"
cell_lines <- c(MEC1 = "MEC1_counts.tsv", HG3 = "HG3_counts.tsv")
dna_columns <- paste0("DNA_Rep", 1:4)
rna_columns <- paste0("RNA_Rep", 1:4)
min_dna_count <- 5

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_valid_counts <- function(path, required_columns, label) {
  if (!file.exists(path)) stop("Missing file: ", path)
  tbl <- read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()), progress = FALSE)
  if (anyDuplicated(names(tbl))) stop(label, ": duplicate column names")
  missing_columns <- setdiff(required_columns, names(tbl))
  if (length(missing_columns)) stop(label, ": missing columns: ", paste(missing_columns, collapse = ", "))
  tbl$element_id <- trimws(tbl$element_id)
  if (anyNA(tbl$element_id) || any(tbl$element_id == "")) stop(label, ": missing element_id values")
  if (anyDuplicated(tbl$element_id)) stop(label, ": duplicate element_id values")
  value_columns <- setdiff(required_columns, "element_id")
  for (name in value_columns) {
    raw <- trimws(tbl[[name]])
    values <- suppressWarnings(as.numeric(raw))
    if (anyNA(values) || any(!is.finite(values)) || any(values < 0) || any(values != floor(values))) {
      stop(label, ": invalid or missing non-negative integer counts in column ", name)
    }
    tbl[[name]] <- values
  }
  tbl %>% select(all_of(required_columns))
}

input_counts <- read_valid_counts(file.path(input_dir, "input_counts.tsv"), c("element_id", "input"), "Input counts")

process_cell_line <- function(cell_line, filename) {
  counts <- read_valid_counts(file.path(input_dir, filename), c("element_id", dna_columns, rna_columns), cell_line)
  joined <- counts %>% left_join(input_counts, by = "element_id", relationship = "many-to-one")
  if (anyNA(joined$input)) {
    missing_ids <- joined$element_id[is.na(joined$input)]
    stop(cell_line, ": missing input counts for ", length(missing_ids), " elements (e.g. ", paste(head(missing_ids, 5), collapse = ", "), ")")
  }
  qc <- joined %>% mutate(
    fail_input = input == 0,
    fail_dna = if_any(all_of(dna_columns), ~ .x < min_dna_count),
    retained = !fail_input & !fail_dna,
    exclusion_reason = case_when(
      fail_input & fail_dna ~ "input_zero_and_low_DNA",
      fail_input ~ "input_zero",
      fail_dna ~ "low_DNA",
      TRUE ~ "retained"
    )
  )
  kept <- qc %>% filter(retained) %>% select(element_id, all_of(dna_columns), all_of(rna_columns))
  write_tsv(kept, file.path(output_dir, paste0(cell_line, "_counts_filtered.tsv")))
  write_tsv(qc %>% select(element_id, input, all_of(dna_columns), fail_input, fail_dna, retained, exclusion_reason),
            file.path(output_dir, paste0(cell_line, "_filter_report.tsv")))
  summary <- tibble(
    cell_line = cell_line, total_elements = nrow(qc),
    input_zero = sum(qc$fail_input), low_DNA_any_rep = sum(qc$fail_dna),
    both_failures = sum(qc$fail_input & qc$fail_dna),
    excluded_unique = sum(!qc$retained), retained = sum(qc$retained)
  )
  message(cell_line, ": ", summary$retained, "/", summary$total_elements, " retained")
  summary
}

summaries <- bind_rows(lapply(names(cell_lines), function(name) process_cell_line(name, cell_lines[[name]])))
write_tsv(summaries, file.path(output_dir, "qc_summary.tsv"))
print(summaries)
