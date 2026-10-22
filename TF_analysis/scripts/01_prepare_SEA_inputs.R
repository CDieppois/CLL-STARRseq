#!/usr/bin/env Rscript
# Run from TF_analysis/; generates common active/inactive FASTA sets.
suppressPackageStartupMessages({library(dplyr); library(readr); library(Biostrings)})

seq_file <- Sys.getenv('TF_SEQUENCE_FILE', unset = 'inputs/oligo_sequences.tsv')
activity_dir <- '../STARR-seq_analysis/results/enhancer_activity'
out_dir <- 'results/SEA/inputs'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_activity <- function(cell) {
  file <- file.path(activity_dir, paste0(cell, '_enhancer_activity.tsv'))
  if (!file.exists(file)) stop('Missing STARR-seq output: ', file)
  d <- read_tsv(file, show_col_types = FALSE)
  if (!all(c('element_id','active_call') %in% names(d))) stop('Missing columns in ', file)
  if (anyDuplicated(d$element_id)) stop('Duplicate element_id in ', file)
  values <- toupper(as.character(d$active_call))
  if (any(!is.na(values) & !values %in% c('TRUE','FALSE','1','0'))) stop('Invalid active_call in ', file)
  tibble(element_id = as.character(d$element_id), active = ifelse(is.na(values), NA, values %in% c('TRUE','1')))
}

if (!file.exists(seq_file)) stop('Missing oligonucleotide sequences: ', seq_file, '\nExpected TSV: element_id, sequence')
sequences <- read_tsv(seq_file, show_col_types = FALSE) %>%
  transmute(element_id = as.character(element_id), sequence = toupper(gsub('\\s+', '', as.character(sequence))))
if (anyDuplicated(sequences$element_id)) stop('Duplicated IDs in oligo_sequences.tsv')
if (anyNA(sequences$element_id) || any(!grepl('^[ACGTN]+$', sequences$sequence))) stop('Missing IDs or invalid sequences in oligo_sequences.tsv')

mec <- read_activity('MEC1') %>% rename(active_MEC1 = active)
hg3 <- read_activity('HG3') %>% rename(active_HG3 = active)
classified <- inner_join(mec, hg3, by = 'element_id') %>%
  left_join(sequences, by = 'element_id') %>%
  mutate(group = case_when(
    is.na(active_MEC1) | is.na(active_HG3) ~ 'not_assessed',
    active_MEC1 & active_HG3 ~ 'active_both',
    !active_MEC1 & !active_HG3 ~ 'inactive_both',
    TRUE ~ 'mixed_activity'
  ))
if (anyNA(classified$sequence)) stop('Missing oligonucleotide sequences for STARR-seq elements: ', paste(head(classified$element_id[is.na(classified$sequence)], 10), collapse=', '))

write_fasta <- function(tbl, filename) {
  if (nrow(tbl) < 1) stop('No sequences for ', filename)
  xs <- DNAStringSet(tbl$sequence); names(xs) <- tbl$element_id
  writeXStringSet(xs, file.path(out_dir, filename), format = 'fasta')
}
write_fasta(filter(classified, group == 'active_both'), 'Both_active.fa')
write_fasta(filter(classified, group == 'inactive_both'), 'Both_inactive.fa')
write_tsv(classified, file.path(out_dir, 'SEA_sequence_classification.tsv'))
write_tsv(count(classified, group, name = 'n_sequences'), file.path(out_dir, 'SEA_sequence_summary.tsv'))
message('SEA FASTA preparation complete: ', out_dir)
