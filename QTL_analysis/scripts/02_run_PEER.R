#!/usr/bin/env Rscript
# PEER latent-factor correction, then inverse-normal transformation per feature.
# Uses original factor counts: expression 65, methylation 45, H3K27ac 10, ATAC 10.
source('scripts/utils.R')
if (!requireNamespace('peer', quietly = TRUE)) stop('Install the PEER R package (peer) in your analysis environment.')

run_peer <- function(type) {
  cfg <- qtl_config[[type]]
  mat <- read_sample_feature(work_path(type, paste0('PEER_matrix_', cfg$phenotype, '.txt')))
  if (anyNA(mat)) stop('PEER input contains missing values.')
  if (nrow(mat) <= cfg$factors) warning(type, ': number of samples <= PEER factor count; validate suitability for original model.')
  set.seed(123)
  model <- peer::PEER()
  peer::PEER_setNk(model, cfg$factors)
  peer::PEER_setPhenoMean(model, mat)
  peer::PEER_update(model)
  residuals <- as.matrix(peer::PEER_getResiduals(model))
  if (!identical(dim(residuals), dim(mat))) stop('Unexpected dimensions returned by PEER_getResiduals().')
  dimnames(residuals) <- dimnames(mat)
  factors <- as.matrix(peer::PEER_getX(model))
  if (nrow(factors) != nrow(mat)) stop('Unexpected PEER factor matrix orientation.')
  colnames(factors) <- paste0('PEER_factor_', seq_len(ncol(factors)))
  rownames(factors) <- rownames(mat)
  transformed <- vapply(seq_len(ncol(residuals)), function(i) invnorm(residuals[, i]), numeric(nrow(residuals)))
  transformed <- matrix(transformed, nrow = nrow(residuals), ncol = ncol(residuals), dimnames = dimnames(mat))
  if (anyNA(transformed)) stop('Inverse-normal transformation produced missing values; inspect inputs.')
  prefix <- paste0('PEER_', cfg$phenotype, '_K', cfg$factors)
  write_sample_feature(residuals, work_path(type, paste0(prefix, '_residuals.tsv')), rownames(mat), colnames(mat))
  write_sample_feature(transformed, work_path(type, paste0(prefix, '_residuals_INT.tsv')), rownames(mat), colnames(mat))
  write_sample_feature(factors, work_path(type, paste0(prefix, '_factors.tsv')), rownames(mat), colnames(factors))
  if ('PEER_getAlpha' %in% getNamespaceExports('peer')) {
    alpha <- as.numeric(peer::PEER_getAlpha(model))
    if (length(alpha) == ncol(factors)) {
      relevance <- data.table(factor_id = colnames(factors), alpha = alpha, relevance = 1 / alpha)
      write_dt(relevance, work_path(type, paste0(prefix, '_factor_relevance.tsv')))
    }
  }
  # PEER model contains patient-level signal: save only under git-ignored results/.
  tryCatch(saveRDS(model, work_path(type, paste0(prefix, '_model.rds'))), error = function(e) warning('PEER model object could not be serialized: ', conditionMessage(e)))
  write_summary(type, 'PEER', list(samples = nrow(mat), features = ncol(mat), factors = cfg$factors, seed = 123))
}
for (type in chosen_analyses()) {
  message('Running PEER ', type)
  run_peer(type)
}
