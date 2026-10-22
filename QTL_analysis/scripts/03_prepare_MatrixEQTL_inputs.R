#!/usr/bin/env Rscript
# Allele-oriented dosage matrices, aligned phenotype matrix, and selected pairs.
source('scripts/utils.R')

prepare_matrixeqtl <- function(type) {
  cfg <- qtl_config[[type]]
  variants <- read_variants()
  peer_input <- work_path(type, paste0('PEER_', cfg$phenotype, '_K', cfg$factors, '_residuals_INT.tsv'))
  residuals <- read_sample_feature(peer_input)
  geno_cols <- sample_fields(variants)
  samples <- rownames(residuals)[rownames(residuals) %in% geno_cols]
  if (length(samples) < 3L) stop(type, ': insufficient shared genotype/phenotype patients.')
  variants <- copy(variants)
  qc <- vector('list', nrow(variants))
  dosage <- matrix(NA_real_, nrow(variants), length(samples), dimnames = list(variants$variant_id, samples))
  for (i in seq_len(nrow(variants))) {
    var <- variants[i]
    codes <- allele_codes(var$genotype_ref, var$genotype_alt, var$tested_ref, var$tested_alt)
    values <- if (anyNA(codes)) rep(NA_real_, length(samples)) else genotype_dosage(unlist(var[, ..samples], use.names = FALSE), codes[1], codes[2])
    dosage[i, ] <- values
    observed <- unique(values[!is.na(values)])
    status <- if (anyNA(codes)) 'tested_alleles_not_resolved' else if (!length(observed)) 'no_valid_genotypes' else if (length(observed) < 2L) 'invariant_dosage' else 'success'
    qc[[i]] <- data.table(variant_id = var$variant_id, tested_ref_code = codes[1], tested_alt_code = codes[2], n_dosage_0 = sum(values == 0, na.rm = TRUE), n_dosage_1 = sum(values == 1, na.rm = TRUE), n_dosage_2 = sum(values == 2, na.rm = TRUE), n_missing = sum(is.na(values)), n_observed_classes = length(observed), status = status)
  }
  qc <- rbindlist(qc)
  write_dt(qc, work_path(type, 'genotype_dosage_QC.tsv'))
  valid <- qc[status == 'success', variant_id]
  if (!length(valid)) stop(type, ': no valid polymorphic tested-allele dosages.')
  pairs <- read_dt(work_path(type, 'selected_variant_feature_pairs.tsv'))
  require_fields(pairs, c('variant_id', 'feature_id'), 'Selected pairs')
  pairs <- unique(pairs[variant_id %in% valid & feature_id %in% colnames(residuals), .(variant_id, feature_id)])
  if (!nrow(pairs)) stop(type, ': no selected pairs remain after genotype QC.')
  valid <- valid[valid %in% pairs$variant_id]
  features <- colnames(residuals)[colnames(residuals) %in% pairs$feature_id]
  # Matrix eQTL requires SNPs/features in rows and samples in columns.
  geno_matrix <- dosage[valid, samples, drop = FALSE]
  pheno_matrix <- t(residuals[samples, features, drop = FALSE])
  geno_dt <- cbind(data.table(snps = rownames(geno_matrix)), as.data.table(geno_matrix))
  setnames(geno_dt, c('snps', samples))
  pheno_dt <- cbind(data.table(gene = rownames(pheno_matrix)), as.data.table(pheno_matrix))
  setnames(pheno_dt, c('gene', samples))
  matrix_pairs <- pairs[, .(snps = variant_id, gene = feature_id)]
  if (!identical(names(geno_dt)[-1], names(pheno_dt)[-1])) stop('Matrix sample order mismatch')
  if (any(!is.finite(pheno_matrix))) stop('Missing or nonfinite phenotype values')
  write_dt(geno_dt, work_path(type, 'MatrixEQTL_genotype_dosage.tsv'))
  write_dt(pheno_dt, work_path(type, 'MatrixEQTL_phenotypes.tsv'))
  write_dt(matrix_pairs, work_path(type, 'MatrixEQTL_selected_pairs.tsv'))
  write_summary(type, 'MatrixEQTL_inputs', list(samples = length(samples), valid_variants = nrow(geno_dt), selected_features = nrow(pheno_dt), selected_pairs = nrow(matrix_pairs)))
}
for (type in chosen_analyses()) {
  message('Preparing Matrix eQTL ', type)
  prepare_matrixeqtl(type)
}
