#!/usr/bin/env Rscript
# Record the R software environment used for the STARR-seq analysis module.
# Run from the STARR-seq_analysis/ directory: Rscript scripts/05_record_environment.R

packages <- c("MPRAnalyze", "dplyr", "readr", "tidyr", "tibble")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

output_dir <- "docs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

package_versions <- vapply(packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1))
bioc_version <- if (requireNamespace("BiocManager", quietly = TRUE)) {
  as.character(BiocManager::version())
} else {
  "BiocManager not installed; Bioconductor version unavailable"
}

report <- c(
  paste("Recorded on:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Bioconductor:", bioc_version),
  "",
  "Direct package dependencies:",
  paste(names(package_versions), package_versions, sep = ": "),
  "",
  "Complete R session information:",
  capture.output(sessionInfo())
)
writeLines(report, con = file.path(output_dir, "software_environment.txt"))
message("Recorded R environment: ", file.path(output_dir, "software_environment.txt"))
