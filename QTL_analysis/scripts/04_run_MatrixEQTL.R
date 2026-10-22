#!/usr/bin/env Rscript
# Original genome-wide variant-by-feature Matrix eQTL linear model, before pair filtering.
source('scripts/utils.R')
if (!requireNamespace('MatrixEQTL', quietly = TRUE)) stop('Install MatrixEQTL.')

load_sliced <- function(path) {
  x <- MatrixEQTL::SlicedData$new()
  x$fileDelimiter <- '\t'
  x$fileOmitCharacters <- 'NA'
  x$fileSkipRows <- 1L
  x$fileSkipColumns <- 1L
  x$fileSliceSize <- 2000L
  x$LoadFile(path)
  x
}

run_matrixeqtl <- function(type) {
  g <- read_dt(work_path(type, 'MatrixEQTL_genotype_dosage.tsv'))
  p <- read_dt(work_path(type, 'MatrixEQTL_phenotypes.tsv'))
  require_fields(g, 'snps', 'Genotype matrix')
  require_fields(p, 'gene', 'Phenotype matrix')
  if (!identical(names(g)[-1], names(p)[-1])) stop('Sample order differs between genotype and phenotype matrices.')
  # Computational cost is n_variants * n_features: this matches the original workflow.
  message(type, ': testing up to ', nrow(g) * nrow(p), ' unfiltered variant-feature combinations')
  snps <- load_sliced(work_path(type, 'MatrixEQTL_genotype_dosage.tsv'))
  gene <- load_sliced(work_path(type, 'MatrixEQTL_phenotypes.tsv'))
  cvrt <- MatrixEQTL::SlicedData$new()
  fit <- MatrixEQTL::Matrix_eQTL_main(
    snps = snps, gene = gene, cvrt = cvrt,
    output_file_name = NULL, pvOutputThreshold = 1,
    useModel = MatrixEQTL::modelLINEAR, errorCovariance = numeric(),
    verbose = TRUE, pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE
  )
  all <- as.data.table(fit$all$eqtls)
  if (!nrow(all)) stop('Matrix eQTL returned no results.')
  require_fields(all, c('snps','gene','beta','statistic','pvalue'), 'Matrix eQTL results')
  write_dt(all, work_path(type, 'MatrixEQTL_all_results.tsv'))
  write_summary(type, 'MatrixEQTL_test', list(variants = nrow(g), features = nrow(p), associations_returned = nrow(all)))
}
for (type in chosen_analyses()) {
  message('Running Matrix eQTL ', type)
  run_matrixeqtl(type)
}
