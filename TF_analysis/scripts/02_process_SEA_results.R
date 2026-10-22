#!/usr/bin/env Rscript
# Process the single SEA analysis: shared active sequences against shuffled dinucleotide-preserving controls.
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
sea_default <- if (file.exists('results/SEA/SEA_Both_active_dinuc/sea.tsv')) 'results/SEA/SEA_Both_active_dinuc/sea.tsv' else 'inputs/SEA/sea.tsv'
sea_file <- Sys.getenv('TF_SEA_RESULTS', unset = sea_default)
motifs_file <- Sys.getenv('TF_SEA_MOTIFS', unset = 'inputs/motif_annotations/H13CORE_motifs.tsv')
expression_file <- Sys.getenv('TF_EXPRESSION_FILE', unset = 'results/RNAseq/MEC1_HG3_gene_expression_TPM.tsv')
out <- 'results/SEA/processed'; dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (f in c(sea_file, motifs_file, expression_file)) if (!file.exists(f)) stop('Missing SEA processing input: ', f)

sea <- read_tsv(sea_file, show_col_types = FALSE, comment = '#')
assert <- function(cols, df, f) {missing <- setdiff(cols, names(df)); if(length(missing))stop('Missing ', paste(missing,collapse=', '), ' in ', f)}
assert(c('ID','QVALUE'), sea, sea_file)

motifs <- read_tsv(motifs_file, show_col_types = FALSE, name_repair = 'unique')
get_column <- function(df, choices) {found <- intersect(choices, names(df)); if (!length(found)) NA_character_ else found[1]}
mc <- get_column(motifs, c('Motif','motif_id','motif'))
gc <- get_column(motifs, c('Gene (human)','Gene..human.','gene_symbol','Gene'))
fc <- get_column(motifs, c('TF family','TF.family','tf_family'))
if (is.na(mc) || is.na(gc)) stop('HOCOMOCO annotations must contain motif IDs and human TF symbols.')
annot <- tibble(motif_key = sub('\\..*$', '', as.character(motifs[[mc]])),
                annotated_tf = as.character(motifs[[gc]]),
                tf_family = if (is.na(fc)) NA_character_ else as.character(motifs[[fc]])) %>%
  distinct(motif_key, .keep_all = TRUE)
expr <- read_tsv(expression_file, show_col_types = FALSE)
assert(c('gene_symbol', 'MEC1_TPM', 'HG3_TPM'), expr, expression_file)
expr <- expr %>% filter(!is.na(gene_symbol), gene_symbol != '') %>%
  mutate(gene_key = toupper(trimws(gene_symbol))) %>%
  arrange(desc(pmax(MEC1_TPM, HG3_TPM, na.rm = TRUE))) %>%
  distinct(gene_key, .keep_all = TRUE)

prepared <- sea %>%
  mutate(motif_key = sub('\\..*$', '', as.character(ID)),
         q_value = as.numeric(QVALUE),
         primary_hits = if ('TP' %in% names(.)) as.numeric(TP) else NA_real_) %>%
  left_join(annot, by = 'motif_key') %>%
  mutate(tf_symbol = if_else(!is.na(annotated_tf) & annotated_tf != '', annotated_tf, motif_key),
         gene_key = toupper(trimws(tf_symbol))) %>%
  left_join(expr %>% select(gene_key, gene_symbol, MEC1_TPM, HG3_TPM), by='gene_key') %>%
  mutate(significant = !is.na(q_value) & q_value < 0.05,
         expressed_both = !is.na(MEC1_TPM) & !is.na(HG3_TPM) & MEC1_TPM >= 20 & HG3_TPM >= 20)
best <- prepared %>% filter(!is.na(tf_symbol), tf_symbol != '', !is.na(q_value)) %>%
  arrange(q_value, desc(primary_hits)) %>% group_by(tf_symbol) %>% slice_head(n=1) %>% ungroup()
write_tsv(prepared, file.path(out, 'SEA_motifs_annotated.tsv'))
write_tsv(best, file.path(out, 'SEA_best_motif_per_TF.tsv'))
write_tsv(filter(best, significant, expressed_both), file.path(out, 'SEA_significant_TFs_expressed_in_both.tsv'))
write_tsv(tibble(total_motifs = nrow(prepared), significant_motifs = sum(prepared$significant),
                 unique_TFs = nrow(best), significant_expressed_TFs = sum(best$significant & best$expressed_both)),
          file.path(out, 'SEA_summary.tsv'))
message('Processed SEA results: ', out)
