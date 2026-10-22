# Shared functions for the controlled-access QTL analyses.
# Run all scripts from the QTL_analysis module root. No individual-level data are distributed.
suppressPackageStartupMessages(library(data.table))

qtl_config <- list(
  eQTL  = list(phenotype = 'gene_expression',  factors = 65L, value_type = 'GCRMA', threshold = 4.5,  overlap = 'within', genome_build = 'hg19'),
  mQTL  = list(phenotype = 'DNA_methylation', factors = 45L, value_type = 'beta',  threshold = NA_real_, overlap = 'within', genome_build = 'hg38'),
  hQTL  = list(phenotype = 'H3K27ac',         factors = 10L, value_type = 'signal', threshold = 0, overlap = 'any', genome_build = 'hg38'),
  caQTL = list(phenotype = 'ATACseq',         factors = 10L, value_type = 'signal', threshold = 0, overlap = 'any', genome_build = 'hg38')
)

chosen_analyses <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || identical(args[1], 'all')) return(names(qtl_config))
  if (length(args) != 1L || !args[1] %in% names(qtl_config)) {
    stop('Usage: Rscript <script> [all|eQTL|mQTL|hQTL|caQTL]')
  }
  args[1]
}

# Read private source data from a local directory outside the GitHub repository.
# Configure QTL_INPUT_DIR in R session/shell/.Renviron to an external local folder.
# Keep real genotype/phenotype matrices out of version control; include only synthetic templates.
input_path <- function(name) {
  data_dir <- Sys.getenv('QTL_INPUT_DIR', unset = '')
  if (!nzchar(data_dir)) {
    stop('QTL_INPUT_DIR is not set. Point it to the local directory containing your authorized QTL input TSV files. See example_data/ for synthetic format examples.')
  }
  file.path(data_dir, name)
}
work_dir <- function(type) file.path('results', type)
work_path <- function(type, name) file.path(work_dir(type), name)
read_dt <- function(path) {
  if (!file.exists(path)) stop('Missing input: ', path)
  fread(path, sep = '\t', header = TRUE, check.names = FALSE, na.strings = c('NA', ''))
}
write_dt <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, sep = '\t', quote = FALSE, na = 'NA')
  message('Saved ', path, ' [', nrow(dt), ' rows]')
}
require_fields <- function(x, fields, label) {
  missing <- setdiff(fields, names(x))
  if (length(missing)) stop(label, ' is missing columns: ', paste(missing, collapse = ', '))
}
require_ids <- function(ids, label) {
  if (anyNA(ids) || any(!nzchar(trimws(as.character(ids))))) stop('Empty IDs in ', label)
  if (anyDuplicated(ids)) stop('Duplicate IDs in ', label, ': ', paste(head(unique(ids[duplicated(ids)]), 5), collapse = ', '))
}
sample_fields <- function(x) grep('^CLL[0-9]+$', names(x), value = TRUE)

read_variants <- function() {
  v <- read_dt(input_path('QTL_variant_genotypes.tsv'))
  required <- c('variant_id', 'genotype_ref', 'genotype_alt', 'tested_ref', 'tested_alt', 'tad_chromosome', 'tad_start_hg38', 'tad_end_hg38')
  require_fields(v, required, 'Genotype/variant table')
  v[, variant_id := trimws(as.character(variant_id))]
  require_ids(v$variant_id, 'variant table')
  s <- sample_fields(v)
  if (!length(s)) stop('No CLL patient genotype columns in variant table.')
  v[, `:=`(tad_chromosome = as.character(tad_chromosome), tad_start_hg38 = as.numeric(tad_start_hg38), tad_end_hg38 = as.numeric(tad_end_hg38))]
  # Original expression probe coordinates are hg19, other phenotype coordinates hg38.
  hg19 <- c('tad_chromosome_hg19', 'tad_start_hg19', 'tad_end_hg19')
  if (all(hg19 %in% names(v))) v[, `:=`(tad_chromosome_hg19 = as.character(tad_chromosome_hg19), tad_start_hg19 = as.numeric(tad_start_hg19), tad_end_hg19 = as.numeric(tad_end_hg19))]
  if (any(!is.na(v$tad_start_hg38) & !is.na(v$tad_end_hg38) & v$tad_start_hg38 > v$tad_end_hg38)) stop('Invalid TAD coordinates')
  v
}

# Read the exact previously published molecular matrices when available.
# A pre-standardized <phenotype>.tsv is also supported. All files are local only.
# Expression GEO alignments are hg19 in the original notebooks; the other assays hg38.
published_source_files <- c(
  eQTL = 'single_probes_all_cases_RMA_expression.txt',
  mQTL = 'cll_methylation_full_matrix_idibaps_norm_Beta_grch38.tsv.gz',
  hQTL = 'H3K27ac_consensus_sample_count_corrected.txt',
  caQTL = 'atac_consensus_sample_count_corrected.txt'
)

read_source_table <- function(path, ...) {
  if (!file.exists(path)) stop('Missing published molecular matrix: ', path)
  if (grepl('\\.gz$', path, ignore.case = TRUE) && !requireNamespace('R.utils', quietly = TRUE)) {
    stop('Reading gzipped methylation TSV requires R.utils. Install with install.packages("R.utils").')
  }
  fread(path, sep = '\t', header = TRUE, check.names = FALSE, na.strings = c('NA',''), ...)
}

read_geo_annotation <- function(path) {
  if (!file.exists(path)) stop('Missing GEO platform annotation file: ', path)
  header <- grep('^ID\t', readLines(path, warn = FALSE))
  if (!length(header)) stop('Cannot locate GEO annotation header beginning with ID in ', path)
  annot <- fread(path, sep = '\t', header = TRUE, skip = header[1] - 1L,
                 fill = TRUE, quote = '', check.names = FALSE)
  require_fields(annot, c('ID', 'Gene Symbol', 'Alignments'), 'GEO GPL13667 annotation')
  annot[, .(feature_id = trimws(as.character(ID)),
            gene_symbol = trimws(as.character(`Gene Symbol`)),
            alignments = as.character(Alignments))][!duplicated(feature_id)]
}

# Extract the first reported probe alignment, exactly as the original eQTL notebook.
extract_probe_coordinates <- function(alignment) {
  x <- as.character(alignment)
  matches <- regmatches(x, regexec('chr([^ :;]+):([0-9]+)-([0-9]+)', x, perl = TRUE))
  data.table(chromosome = vapply(matches, function(m) if (length(m) == 4L) paste0('chr', m[2]) else NA_character_, character(1)),
             start = as.numeric(vapply(matches, function(m) if (length(m) == 4L) m[3] else NA_character_, character(1))),
             end = as.numeric(vapply(matches, function(m) if (length(m) == 4L) m[4] else NA_character_, character(1))))
}

read_published_phenotype <- function(type) {
  path <- input_path(published_source_files[[type]])
  # Preserve the source filename as stated in the original notebooks.
  if (type == 'eQTL') {
    expression <- read_source_table(path)
    require_fields(expression, 'Probe_Set_ID', 'Published expression matrix')
    setnames(expression, 'Probe_Set_ID', 'feature_id')
    expression[, feature_id := trimws(as.character(feature_id))]
    annotation <- read_geo_annotation(input_path('GPL13667-15572.txt'))
    expression <- merge(expression, annotation, by = 'feature_id', all.x = TRUE, sort = FALSE)
    expression <- expression[!is.na(gene_symbol) & nzchar(gene_symbol) & gene_symbol != '---']
    positions <- extract_probe_coordinates(expression$alignments)
    expression[, `:=`(chromosome = positions$chromosome, start = positions$start, end = positions$end)]
    expression[, alignments := NULL]
    expression <- expression[!is.na(chromosome) & !is.na(start) & !is.na(end)]
    return(expression)
  }
  source <- read_source_table(path)
  require_fields(source, c('chr','start','end'), paste0(type, ' published matrix'))
  setnames(source, 'chr', 'chromosome')
  if (type == 'mQTL') {
    require_fields(source, 'name', 'Published DNA methylation matrix')
    setnames(source, 'name', 'feature_id')
    numeric_samples <- grep('^[0-9]+$', names(source), value = TRUE)
    if (!length(numeric_samples)) stop('No numerical patient columns in published methylation matrix')
    new_names <- paste0('CLL', as.integer(numeric_samples))
    if (anyDuplicated(new_names)) stop('Duplicate methylation CLL IDs after renaming')
    setnames(source, numeric_samples, new_names)
  } else {
    source[, feature_id := paste0(chromosome, ':', start, '-', end)]
  }
  source
}

read_phenotype <- function(type) {
  cfg <- qtl_config[[type]]
  standardized <- input_path(paste0(cfg$phenotype, '.tsv'))
  # Prefer standardized files when users already converted the published matrices.
  x <- if (file.exists(standardized)) {
    message(type, ': reading locally standardized phenotype file ', standardized)
    read_dt(standardized)
  } else {
    message(type, ': reading original published matrix ', published_source_files[[type]])
    read_published_phenotype(type)
  }
  require_fields(x, c('feature_id', 'chromosome', 'start', 'end'), cfg$phenotype)
  x[, `:=`(feature_id = trimws(as.character(feature_id)), chromosome = as.character(chromosome),
           start = as.numeric(start), end = as.numeric(end))]
  x <- x[!is.na(feature_id) & nzchar(feature_id) & !is.na(chromosome) &
           !is.na(start) & !is.na(end) & start <= end]
  require_ids(x$feature_id, paste0(type, ' features'))
  if (type == 'eQTL' && !'gene_symbol' %in% names(x)) x[, gene_symbol := NA_character_]
  samples <- sample_fields(x)
  if (!length(samples)) stop('No CLL patient columns in ', cfg$phenotype)
  x[, (samples) := lapply(.SD, function(y) suppressWarnings(as.numeric(y))), .SDcols = samples]
  values <- as.matrix(x[, ..samples])
  if (any(is.infinite(values))) stop('Infinite phenotype values in ', cfg$phenotype)
  if (type == 'mQTL' && any(values < 0 | values > 1, na.rm = TRUE)) stop('Methylation beta values must be in [0,1]')
  x
}

clean_allele <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c('', '.', 'NA')] <- NA_character_
  x
}
normalize_alleles <- function(a, b) {
  alleles <- clean_allele(c(a, b))
  if (anyNA(alleles)) return(rep(NA_character_, 2))
  a <- alleles[1]; b <- alleles[2]
  while (nchar(a) > 1 && nchar(b) > 1 && substr(a, 1, 1) == substr(b, 1, 1)) {
    a <- substring(a, 2); b <- substring(b, 2)
  }
  while (nchar(a) > 1 && nchar(b) > 1 && substr(a, nchar(a), nchar(a)) == substr(b, nchar(b), nchar(b))) {
    a <- substr(a, 1, nchar(a) - 1); b <- substr(b, 1, nchar(b) - 1)
  }
  c(a, b)
}

allele_codes <- function(geno_ref, geno_alt, tested_ref, tested_alt) {
  alleles <- clean_allele(c(geno_ref, strsplit(as.character(geno_alt), ',', fixed = TRUE)[[1]]))
  if (anyNA(alleles) || length(alleles) < 2L) return(c(NA_integer_, NA_integer_))
  tested <- clean_allele(c(tested_ref, tested_alt))
  if (anyNA(tested)) return(c(NA_integer_, NA_integer_))
  exact <- match(tested, alleles) - 1L
  if (!anyNA(exact) && length(unique(exact)) == 2L) return(as.integer(exact))
  norm_test <- normalize_alleles(tested[1], tested[2])
  for (r in seq_along(alleles)) for (a in seq_along(alleles)) {
    if (r == a) next
    if (identical(normalize_alleles(alleles[r], alleles[a]), norm_test)) return(c(r - 1L, a - 1L))
  }
  c(NA_integer_, NA_integer_)
}

genotype_dosage <- function(genotypes, ref_code, alt_code) {
  vapply(as.character(genotypes), function(g) {
    if (is.na(g) || g %in% c('', '.', './.', '.|.')) return(NA_real_)
    codes <- suppressWarnings(as.integer(strsplit(gsub('\\|', '/', g), '/', fixed = TRUE)[[1]]))
    if (length(codes) != 2L || anyNA(codes)) return(NA_real_)
    if (all(codes == ref_code)) return(0)
    if (all(codes == alt_code)) return(2)
    if (setequal(codes, c(ref_code, alt_code))) return(1)
    NA_real_
  }, numeric(1))
}

invnorm <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) < 2L) return(out)
  out[ok] <- qnorm((rank(x[ok], ties.method = 'average') - 0.5) / sum(ok))
  out
}

# Standardized matrices always use the first column as the identifier.
write_sample_feature <- function(mat, path, sample_ids, feature_ids) {
  if (!identical(dim(mat), c(length(sample_ids), length(feature_ids)))) stop('Sample/feature matrix dimensions do not match IDs.')
  tab <- as.data.table(mat)
  setnames(tab, feature_ids)
  tab[, sample_id := sample_ids]
  setcolorder(tab, c('sample_id', feature_ids))
  write_dt(tab, path)
}
read_sample_feature <- function(path) {
  dt <- read_dt(path)
  require_fields(dt, 'sample_id', basename(path))
  ids <- as.character(dt$sample_id)
  require_ids(ids, 'sample matrix')
  features <- setdiff(names(dt), 'sample_id')
  if (!length(features)) stop('No phenotype features in ', path)
  mat <- as.matrix(dt[, ..features])
  storage.mode(mat) <- 'numeric'
  if (any(is.infinite(mat))) stop('Infinite values in ', path)
  dimnames(mat) <- list(ids, features)
  mat
}

feature_overlap_pairs <- function(v, features, mode, type) {
  if (type == 'eQTL') {
    required <- c('tad_chromosome_hg19','tad_start_hg19','tad_end_hg19')
    require_fields(v, required, 'Expression hg19 TAD coordinates')
    keep <- which(!is.na(v$tad_chromosome_hg19) & !is.na(v$tad_start_hg19) & !is.na(v$tad_end_hg19))
    if (!length(keep)) stop('No valid hg19 TAD intervals for expression analysis')
    tads <- v[keep, .(variant_id, chromosome = tad_chromosome_hg19, start = tad_start_hg19, end = tad_end_hg19)]
  } else {
    keep <- which(!is.na(v$tad_chromosome) & !is.na(v$tad_start_hg38) & !is.na(v$tad_end_hg38))
    if (!length(keep)) stop('No variants have valid hg38 TAD coordinates')
    tads <- v[keep, .(variant_id, chromosome = tad_chromosome, start = tad_start_hg38, end = tad_end_hg38)]
  }
  f <- features[, .(feature_id, chromosome, start, end)]
  setkey(tads, chromosome, start, end)
  match <- foverlaps(f, tads, type = mode, nomatch = NULL)
  if (!nrow(match)) stop('No variant-feature TAD overlaps found; check genome builds and coordinates.')
  unique(match[, .(variant_id, feature_id)])
}

write_summary <- function(type, stage, metrics) {
  x <- data.table(metric = names(metrics), value = as.character(unlist(metrics)))
  write_dt(x, work_path(type, paste0(stage, '_summary.tsv')))
}
