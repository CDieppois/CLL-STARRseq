#!/usr/bin/env Rscript
# Build hg38 SNP VCF for FABIAN from significant allelic effects in enhancer-active constructs.
# Run from TF_analysis/; requires reference/alternative oligo sequences and hg38 position metadata.
suppressPackageStartupMessages({library(dplyr); library(readr)})
results_dir <- '../STARR-seq_analysis/results/integrated'
seq_path <- Sys.getenv('TF_SEQUENCE_FILE', unset = 'inputs/oligo_sequences.tsv')
pos_path <- Sys.getenv('TF_VARIANT_POSITIONS', unset = 'inputs/variant_positions_hg38.tsv')
out <- 'results/FABIAN/inputs'; dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_ae <- function(cell) {
  path <- file.path(results_dir, paste0(cell, '_allelic_effects_annotated.tsv'))
  if (!file.exists(path)) stop('Missing STARR-seq analysis output: ', path)
  x <- read_tsv(path, show_col_types = FALSE)
  req <- c('contrast_id', 'variant_core', 'variant_id', 'allelic_effect_call', 'sig_enhancer')
  if (!all(req %in% names(x))) stop('Missing columns in ', path, ': ', paste(setdiff(req, names(x)),collapse=', '))
  x %>% transmute(contrast_id = as.character(contrast_id),
                   variant_core = as.character(variant_core), variant_id = as.character(variant_id),
                   ae = toupper(as.character(allelic_effect_call)) %in% c('TRUE','1'),
                   enhancer = toupper(as.character(sig_enhancer)) %in% c('TRUE','1'), cell_line = cell)
}
for (p in c(seq_path,pos_path)) if(!file.exists(p)) stop('Missing input: ', p)
sequences <- read_tsv(seq_path, show_col_types=FALSE) %>%
  transmute(element_id = as.character(element_id), sequence = toupper(gsub('\\s+', '', as.character(sequence))))
if (anyDuplicated(sequences$element_id)) stop('Duplicated element IDs in ', seq_path)
positions <- read_tsv(pos_path, show_col_types=FALSE)
needed <- c('variant_core', 'variant_id', 'chromosome', 'position_hg38', 'genomic_ref_hg38')
if (!all(needed %in% names(positions))) stop('Missing position fields: ', paste(setdiff(needed,names(positions)),collapse=', '))
if (anyDuplicated(positions$variant_core)) stop('Variant positions need one record per variant_core.')

eligible <- bind_rows(read_ae('MEC1'),read_ae('HG3')) %>%
  group_by(contrast_id,variant_core,variant_id) %>%
  summarise(significant_allelic_effect = any(ae), enhancer_either_cell = any(enhancer), .groups='drop') %>%
  filter(significant_allelic_effect, enhancer_either_cell) %>%
  left_join(positions, by = c('variant_core','variant_id'))

make_variant <- function(row) {
  core <- as.character(row$variant_core)
  contrast <- as.character(row$contrast_id)
  match <- regexec('^(.+)_ref_vs_(alt[0-9]*)$', contrast)
  parsed <- regmatches(contrast, match)[[1]]
  problem <- function(reason) tibble(contrast_id=contrast,variant_core=core,variant_id=row$variant_id,
                                     chromosome=row$chromosome,position_hg38=row$position_hg38,
                                     REF=NA_character_,ALT=NA_character_,reason=reason)
  if(length(parsed)!=3 || parsed[2]!=core) return(problem('unrecognized_contrast_id'))
  ref_ids <- sequences$element_id[startsWith(sequences$element_id,paste0(core,'_ref_'))]
  alt_ids <- sequences$element_id[startsWith(sequences$element_id,paste0(core,'_',parsed[3],'_'))]
  if(length(ref_ids)!=1 || length(alt_ids)!=1) return(problem('missing_or_ambiguous_oligonucleotide_sequence'))
  ref_seq <- sequences$sequence[match(ref_ids,sequences$element_id)]
  alt_seq <- sequences$sequence[match(alt_ids,sequences$element_id)]
  if(nchar(ref_seq)!=nchar(alt_seq)) return(problem('non_SNP_different_sequence_lengths'))
  a <- strsplit(ref_seq,'',fixed=TRUE)[[1]]; b <- strsplit(alt_seq,'',fixed=TRUE)[[1]]
  differences <- which(a!=b)
  if(length(differences)!=1) return(problem('non_SNP_expected_one_substitution'))
  p <- differences[[1]]
  if(p <= 30 || p > length(a)-30) return(problem('fewer_than_30_flanking_bases'))
  ref <- a[p]; alt <- b[p]
  if(!grepl('^[ACGT]$',ref) || !grepl('^[ACGT]$',alt)) return(problem('non_ACGT_allele'))
  if(is.na(row$genomic_ref_hg38) || toupper(as.character(row$genomic_ref_hg38))!=ref) return(problem('reference_allele_does_not_match_hg38'))
  if(is.na(row$chromosome) || is.na(row$position_hg38) || is.na(row$variant_id)) return(problem('missing_genomic_coordinates'))
  tibble(contrast_id=contrast, variant_core=core, variant_id=row$variant_id,
         chromosome=as.character(row$chromosome), position_hg38=as.integer(row$position_hg38),
         REF=ref, ALT=alt, reason='eligible')
}
if (!nrow(eligible)) stop('No enhancer-active, significant allele-specific contrasts available.')
rows <- bind_rows(lapply(seq_len(nrow(eligible)), function(i) make_variant(eligible[i,])))
write_tsv(rows, file.path(out, 'FABIAN_variant_eligibility.tsv'))
accepted <- rows %>% filter(reason=='eligible') %>%
  mutate(CHROM=sub('^chr','',chromosome,ignore.case=TRUE), POS=position_hg38,
         ID=variant_id, QUAL='.', FILTER='PASS', INFO='.', FORMAT='GT', SAMPLE='0/1',
         variant_key=paste0('chr',CHROM,':',POS,REF,'>',ALT)) %>%
  distinct(CHROM, POS, REF, ALT, .keep_all=TRUE)
if (!nrow(accepted)) stop('No eligible SNPs remain; inspect FABIAN_variant_eligibility.tsv')
write_tsv(accepted %>% select(variant_key, variant_id, variant_core, contrast_id, chromosome, position_hg38, REF, ALT),
          file.path(out, 'FABIAN_variant_metadata.tsv'))
vcf <- accepted %>% arrange(CHROM,POS) %>% select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE)
vcf_file <- file.path(out,'FABIAN_input_hg38.vcf')
writeLines(c('##fileformat=VCFv4.3', '##reference=GRCh38', paste0('#',paste(names(vcf),collapse='\t'))), vcf_file)
write.table(vcf,file=vcf_file,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,na='.')
message('FABIAN VCF variants: ', nrow(vcf), '; exclusions: ', sum(rows$reason!='eligible'))
