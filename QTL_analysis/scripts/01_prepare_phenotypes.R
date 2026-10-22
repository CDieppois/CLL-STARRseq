#!/usr/bin/env Rscript
# Select TAD-linked features, apply original phenotype/genotype-group presence filters,
# generate sample-by-feature PEER matrices. Expression/GEO coordinates: hg19; ATAC/H3K27ac/methylation: hg38.
source('scripts/utils.R')

prepare_phenotype <- function(type) {
  cfg <- qtl_config[[type]]
  variants <- read_variants()
  feature_data <- read_phenotype(type)
  samples <- intersect(sample_fields(variants), sample_fields(feature_data))
  if (length(samples) < 2L) stop(type, ': fewer than two matched patients.')
  pairs <- feature_overlap_pairs(variants, feature_data, cfg$overlap, type)
  pairs <- pairs[variant_id %in% variants$variant_id]
  if (!nrow(pairs)) stop(type, ': no TAD-linked features.')
  setkey(variants, variant_id)
  setkey(feature_data, feature_id)
  # Original eQTL, caQTL and hQTL filtering: in >=2 or >=10% of patients
  # in at least one genotype class (hom REF, heterozygote or hom ALT).
  # Methylation uses TAD CpG selection without this phenotype-presence filter.
  if (type != 'mQTL') {
    passed <- logical(nrow(pairs))
    for (i in seq_len(nrow(pairs))) {
      var <- variants[pairs$variant_id[i]]
      codes <- allele_codes(var$genotype_ref, var$genotype_alt, var$tested_ref, var$tested_alt)
      if (anyNA(codes)) next
      gt <- genotype_dosage(unlist(var[, ..samples], use.names = FALSE), codes[1], codes[2])
      values <- as.numeric(unlist(feature_data[pairs$feature_id[i], ..samples], use.names = FALSE))
      for (g in 0:2) {
        idx <- which(!is.na(gt) & gt == g)
        n <- length(idx)
        if (n < 2L) next
        required <- max(2L, ceiling(0.10 * n))
        if (sum(is.finite(values[idx]) & values[idx] > cfg$threshold) >= required) {
          passed[i] <- TRUE
          break
        }
      }
    }
    pairs[, passed_genotype_group_filter := passed]
    selected <- pairs[passed_genotype_group_filter == TRUE, .(variant_id, feature_id)]
  } else {
    pairs[, passed_genotype_group_filter := TRUE]
    selected <- pairs[, .(variant_id, feature_id)]
  }
  write_dt(pairs, work_path(type, 'variant_feature_selection_QC.tsv'))
  if (!nrow(selected)) stop(type, ': no pairs survived feature filtering.')
  ids <- unique(selected$feature_id)
  feature_sub <- feature_data[ids]
  if (anyNA(feature_sub$feature_id)) stop('Feature IDs lost during subsetting')
  measure <- as.matrix(feature_sub[, ..samples]); storage.mode(measure) <- 'numeric'
  if (type == 'mQTL') {
    eps <- 1e-6
    measure <- log2(pmin(pmax(measure, eps), 1 - eps) / (1 - pmin(pmax(measure, eps), 1 - eps)))
  }
  peer_matrix <- t(measure)
  dimnames(peer_matrix) <- list(samples, ids)
  if (any(!is.finite(peer_matrix))) {
    # Do not invent missing-value imputation; it was not part of the supplied pipeline.
    stop(type, ': missing/nonfinite values in PEER input. Resolve using authorized source data before rerunning.')
  }
  write_sample_feature(peer_matrix, work_path(type, paste0('PEER_matrix_', cfg$phenotype, '.txt')), samples, ids)
  annotation <- feature_sub[, intersect(c('feature_id','gene_symbol','chromosome','start','end'), names(feature_sub)), with = FALSE]
  write_dt(annotation, work_path(type, 'feature_annotations.tsv'))
  write_dt(unique(selected), work_path(type, 'selected_variant_feature_pairs.tsv'))
  write_summary(type, 'phenotype', list(samples = length(samples), input_features = nrow(feature_data), TAD_linked_pairs = nrow(pairs), selected_pairs = nrow(selected), selected_features = length(ids), peer_factors = cfg$factors))
}

for (type in chosen_analyses()) {
  message('Preparing ', type)
  prepare_phenotype(type)
}
