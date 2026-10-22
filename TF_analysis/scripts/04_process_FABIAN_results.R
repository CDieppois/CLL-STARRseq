#!/usr/bin/env Rscript
# Summarize FABIAN binding changes, favoring TFFM models and annotating TF expression.
# Input: original FABIAN 'fabian_KnownTFBS.data' (14 tab-separated fields, no header).
suppressPackageStartupMessages({library(dplyr);library(readr);library(tidyr);library(stringr)})
out <- 'results/FABIAN/processed'; dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw_file <- Sys.getenv('TF_FABIAN_RESULTS', unset = 'inputs/FABIAN/fabian_KnownTFBS.data')
motif_file <- Sys.getenv('TF_FABIAN_MOTIFS', unset = 'inputs/motif_annotations/H14CORE_motifs.tsv')
jaspar_file <- Sys.getenv('TF_JASPAR_META', unset = 'inputs/motif_annotations/ultimate_metadata_table_CORE.tsv')
expr_file <- Sys.getenv('TF_EXPRESSION_FILE', unset = 'results/RNAseq/MEC1_HG3_gene_expression_TPM.tsv')
meta_file <- 'results/FABIAN/inputs/FABIAN_variant_metadata.tsv'
for (file in c(raw_file,motif_file,jaspar_file,expr_file,meta_file)) if(!file.exists(file))stop('Missing FABIAN input: ',file)

fabian_columns <- c('variant','tf','model_id','database','wt_score','mt_score',
 'start_wt','end_wt','start_mt','end_mt','strand_wt','strand_mt','prediction','score')
raw <- read_tsv(raw_file,col_names=fabian_columns,show_col_types=FALSE,na=c('','NA')) %>%
  mutate(variant_key=sub('\\.[0-9]+$','',as.character(variant)),
         score=as.numeric(score),is_tffm=str_detect(str_to_upper(database),'TFFM'),
         motif_key=sub('\\..*$', '', as.character(model_id)),
         tf_from_model=if_else(str_detect(model_id,'@'),sub('.*@','',model_id),NA_character_))
combined <- raw %>% filter(!is.na(score)) %>% group_by(variant_key,tf) %>%
  filter(if(any(is_tffm)) is_tffm else TRUE) %>%
  summarise(combined_score=mean(score),n_models=n(),
            n_gain_models=sum(score>0),n_loss_models=sum(score<0),n_zero_models=sum(score==0),
            model_set=if_else(any(is_tffm),'TFFM','PWM'),
            contributing_models=paste(sort(unique(model_id)),collapse=';'),
            contributing_databases=paste(sort(unique(database)),collapse=';'),.groups='drop') %>%
  mutate(binding_prediction=case_when(combined_score>0~'gain',combined_score<0~'loss',TRUE~'none'))
write_tsv(combined,file.path(out,'FABIAN_combined_scores.tsv'))
retained <- filter(combined,abs(combined_score)>=0.6)

# HOCOMOCO v14 motif annotations and TF name-to-gene fallback.
motifs <- read_tsv(motif_file,show_col_types=FALSE,name_repair='unique')
get_col <- function(df,options){x<-intersect(options,names(df));if(length(x))x[1] else NA_character_}
mc <- get_col(motifs,c('Motif','motif_id'))
gc <- get_col(motifs,c('Gene (human)','Gene..human.','Gene','gene_symbol'))
fc <- get_col(motifs,c('TF family','TF.family','tf_family'))
if(is.na(mc)||is.na(gc))stop('HOCOMOCO motif file needs Motif and human gene fields')
map <- tibble(motif_key=sub('\\..*$', '',as.character(motifs[[mc]])),
              gene=as.character(motifs[[gc]]),
              tf_family=if(is.na(fc))NA_character_ else as.character(motifs[[fc]])) %>%
  separate_rows(gene,sep='::|;|,') %>% mutate(gene=str_trim(gene)) %>%
  filter(!is.na(gene),gene!='') %>% distinct()
expr <- read_tsv(expr_file,show_col_types=FALSE)
if(!all(c('gene_symbol','MEC1_TPM','HG3_TPM','expression_category_MEC1','expression_category_HG3') %in% names(expr)))
  stop('Expression table lacks required gene, TPM, or category fields')
expr <- expr %>% filter(!is.na(gene_symbol),gene_symbol!='') %>%
  mutate(gene_key=str_to_upper(str_trim(gene_symbol))) %>% distinct(gene_key,.keep_all=TRUE)

candidates <- raw %>% semi_join(retained,by=c('variant_key','tf')) %>%
  left_join(map,by='motif_key') %>%
  transmute(variant_key,tf,model_id,motif_key,tf_family,
            direct_gene=tf,model_gene=tf_from_model,motif_gene=gene) %>%
  pivot_longer(c(direct_gene,model_gene,motif_gene),names_to='gene_source',values_to='gene_symbol') %>%
  filter(!is.na(gene_symbol),gene_symbol!='') %>%
  separate_rows(gene_symbol,sep='::|;|,') %>%
  mutate(gene_key=str_to_upper(str_trim(gene_symbol))) %>%
  left_join(expr %>% select(gene_key,MEC1_TPM,HG3_TPM,expression_category_MEC1,expression_category_HG3),by='gene_key') %>%
  distinct()
expr_summary <- candidates %>% group_by(variant_key,tf) %>%
  summarise(mapped_genes=paste(sort(unique(gene_symbol)),collapse=';'),
    tf_families=paste(sort(unique(na.omit(tf_family))),collapse=';'),
    MEC1_TPM_max=if(all(is.na(MEC1_TPM)))NA_real_ else max(MEC1_TPM,na.rm=TRUE),
    HG3_TPM_max=if(all(is.na(HG3_TPM)))NA_real_ else max(HG3_TPM,na.rm=TRUE),
    expressed_MEC1=any(expression_category_MEC1 %in% c('high','very_high'),na.rm=TRUE),
    expressed_HG3=any(expression_category_HG3 %in% c('high','very_high'),na.rm=TRUE),.groups='drop') %>%
  mutate(expressed_either=expressed_MEC1|expressed_HG3)

annotated <- retained %>% left_join(expr_summary,by=c('variant_key','tf')) %>%
  mutate(expressed_either=coalesce(expressed_either,FALSE))
metadata <- read_tsv(meta_file,show_col_types=FALSE) %>% distinct(variant_key,.keep_all=TRUE)
annotated <- annotated %>% left_join(metadata,by='variant_key')

# JASPAR CORE human TF families (optional annotations within required metadata file).
jaspar <- read_tsv(jaspar_file,show_col_types=FALSE)
if(all(c('tax_id','name','class','family') %in% names(jaspar))) {
  jaspar <- jaspar %>% filter(as.character(tax_id)=='9606') %>%
    transmute(gene_key=str_to_upper(str_trim(name)),jaspar_class=as.character(class),jaspar_family=as.character(family)) %>% distinct()
  jaspar_map <- candidates %>% inner_join(jaspar,by='gene_key') %>% group_by(variant_key,tf) %>%
    summarise(jaspar_class=paste(sort(unique(na.omit(jaspar_class))),collapse=';'),
              jaspar_family=paste(sort(unique(na.omit(jaspar_family))),collapse=';'),.groups='drop')
  annotated <- left_join(annotated,jaspar_map,by=c('variant_key','tf'))
} else warning('JASPAR metadata columns not recognized, leaving family annotation to HOCOMOCO.')

write_tsv(annotated,file.path(out,'FABIAN_thresholded_annotated.tsv'))
write_tsv(filter(annotated,expressed_either),file.path(out,'FABIAN_significant_predicted_changes_expressed_either.tsv'))
write_tsv(filter(annotated,!expressed_either),file.path(out,'FABIAN_excluded_by_expression.tsv'))
write_tsv(tibble(total_variant_tf_pairs=nrow(combined),pairs_abs_score_ge_0_6=nrow(retained),
                 pairs_expressed_either=sum(annotated$expressed_either),
                 unique_variants_expressed_either=n_distinct(annotated$variant_key[annotated$expressed_either])),
          file.path(out,'FABIAN_summary.tsv'))
message('Processed FABIAN results: ',out)
