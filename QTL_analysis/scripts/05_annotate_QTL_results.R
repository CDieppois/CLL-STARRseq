#!/usr/bin/env Rscript
# Select predefined TAD-linked variant-feature pairs, then apply BH correction
# over that selected set (NOT the full Matrix eQTL Cartesian product).
source('scripts/utils.R')

annotate_qtl <- function(type) {
  raw <- read_dt(work_path(type, 'MatrixEQTL_all_results.tsv'))
  selected <- read_dt(work_path(type, 'MatrixEQTL_selected_pairs.tsv'))
  require_fields(raw, c('snps','gene','beta','statistic','pvalue'), 'Matrix eQTL output')
  require_fields(selected, c('snps','gene'), 'Selected pairs')
  selected <- unique(selected[, .(snps = as.character(snps), gene = as.character(gene))])
  raw[, `:=`(snps = as.character(snps), gene = as.character(gene))]
  if (anyDuplicated(raw, by = c('snps','gene'))) stop('Matrix eQTL result has duplicated variant-feature pairs.')
  setkey(raw, snps, gene)
  setkey(selected, snps, gene)
  result <- raw[selected, nomatch = 0]
  if (!nrow(result)) stop(type, ': no selected tests could be recovered.')
  if (nrow(result) != nrow(selected)) warning(type, ': ', nrow(selected) - nrow(result), ' selected pairs were missing from Matrix eQTL output.')
  result[, fdr_selected_pairs := p.adjust(pvalue, method = 'BH')]
  features <- read_dt(work_path(type, 'feature_annotations.tsv'))
  variants <- read_variants()
  # Never export per-patient genotype columns as part of summary association results.
  geno_cols <- sample_fields(variants)
  variant_info <- unique(variants[, setdiff(names(variants), geno_cols), with = FALSE], by = 'variant_id')
  setnames(variant_info, 'variant_id', 'snps')
  setnames(features, 'feature_id', 'gene')
  for (col in intersect(c('chromosome', 'start', 'end'), names(features))) setnames(features, col, paste0('feature_', col))
  for (col in intersect(c('chromosome', 'start_hg38', 'end_hg38'), names(variant_info))) setnames(variant_info, col, paste0('variant_', col))
  if (anyDuplicated(features$gene)) stop('Duplicated feature annotations')
  result <- merge(result, features, by = 'gene', all.x = TRUE, sort = FALSE)
  result <- merge(result, variant_info, by = 'snps', all.x = TRUE, sort = FALSE)
  setnames(result, c('snps','gene','fdr_selected_pairs'), c('variant_id','feature_id','fdr'))
  setcolorder(result, c('variant_id','feature_id','beta','statistic','pvalue','fdr', setdiff(names(result), c('variant_id','feature_id','beta','statistic','pvalue','fdr'))))
  setorder(result, pvalue)
  write_dt(result, work_path(type, paste0(type, '_annotated_results.tsv')))
  write_dt(result[fdr < 0.05], work_path(type, paste0(type, '_significant_results.tsv')))
  write_summary(type, 'QTL_results', list(selected_pairs = nrow(selected), tested_selected_pairs = nrow(result), nominal_p_005 = sum(result$pvalue < 0.05, na.rm = TRUE), FDR_005 = sum(result$fdr < 0.05, na.rm = TRUE), significant_variants = uniqueN(result[fdr < 0.05, variant_id])))
}
for (type in chosen_analyses()) {
  message('Annotating ', type)
  annotate_qtl(type)
}
