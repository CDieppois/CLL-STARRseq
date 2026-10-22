#!/usr/bin/env Rscript
# Reproduce gene-level TPM processing from LL100_RNA-seq.Rmd for MEC-1 and HG-3.
# Run from TF_analysis/. The Ensembl/BioMart gene mapping is queried live, as in the original notebook.
suppressPackageStartupMessages({library(tximport); library(GenomicFeatures); library(AnnotationDbi); library(biomaRt); library(readr)})

root <- Sys.getenv('TF_RNASEQ_WORK_DIR', unset = file.path(getwd(), 'local_RNAseq_data'))
gtf <- file.path(root, 'reference', 'gencode.v26.annotation.gtf')
files <- c(MEC1 = file.path(root, 'MEC1', 'MEC1_salmon_quant', 'quant.sf'),
           HG3 = file.path(root, 'HG3', 'HG3_salmon_quant', 'quant.sf'))
if (!file.exists(gtf) || !all(file.exists(files))) stop('Missing GENCODE v26 GTF or Salmon quant.sf. Run 01_process_LL100_RNAseq.sh first.')

# Same GENCODE v26 transcript-to-gene mapping used in the original notebook.
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
tx_names <- AnnotationDbi::keys(txdb, keytype = 'TXNAME')
tx2gene <- AnnotationDbi::select(txdb, keys = tx_names, keytype = 'TXNAME', columns = 'GENEID')
tx2gene$TXNAME <- sub('\\..*', '', tx2gene$TXNAME)
quant <- tximport::tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
tpm <- quant$abundance
if (!all(c('MEC1', 'HG3') %in% colnames(tpm))) stop('The Salmon data must include both MEC1 and HG3.')

# Preserve the original biomaRt HGNC annotation method and INNER merge behavior.
# Annotation can change with live Ensembl releases; compare symbols with the original exported table.
gene_ids <- sub('\\..*', '', rownames(tpm))
mart <- biomaRt::useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl')
annot <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                        filters = 'ensembl_gene_id', values = gene_ids, mart = mart)
rownames(tpm) <- gene_ids
expression <- merge(annot, tpm, by.x = 'ensembl_gene_id', by.y = 'row.names')
expression <- expression[, c('ensembl_gene_id', 'hgnc_symbol', 'MEC1', 'HG3')]
names(expression) <- c('ensembl_gene_id', 'gene_symbol', 'MEC1_TPM', 'HG3_TPM')

expression_category <- function(values) as.character(cut(values,
  breaks = c(-Inf, 1, 10, 20, 30, 100, Inf),
  labels = c('no', 'very_low', 'low', 'moderate', 'high', 'very_high'),
  right = FALSE, ordered_result = TRUE))
expression$expression_category_MEC1 <- expression_category(expression$MEC1_TPM)
expression$expression_category_HG3 <- expression_category(expression$HG3_TPM)

out <- 'results/RNAseq'; dir.create(out, recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(expression, file.path(out, 'MEC1_HG3_gene_expression_TPM.tsv'))
message('Gene-level TPM table saved. Rows: ', nrow(expression))
