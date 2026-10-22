# QTL pipeline validation checklist

**Status: not validated against the original controlled-access datasets.** The script structure and data-reader assumptions were reviewed against the supplied notebooks, but no execution with the real inputs was possible.

Before a research release, for each QTL class:

1. Check that `01_prepare_phenotypes.R` reads the downloaded original publication files without changing the source measurements. Verify sample IDs and retained features.
2. Confirm coordinates: eQTL/GPL13667 probe alignments and their selected TADs **hg19**; H3K27ac/ATAC/CpGs and their TADs **hg38**. Confirm lifted TAD endpoints used in the original notebook and interval-coordinate conventions.
3. Compare TAD-linked pair counts and the genotype-group filter with original notebooks (eQTL >4.5; ATAC/H3K27ac >0, at least two or 10% in a genotype class; CpG no such presence filter).
4. Compare PEER input dimensions, 65/45/10/10 factors, residuals and transformed residuals. Record original PEER version and RNG seed; the draft currently uses seed 123, not verified against the original Linux run.
5. Confirm genotype allele coding and dosage orientation, including indels and multiallelic genotypes; check genotype QC and patient order.
6. Confirm the Matrix eQTL model, sample matching, p-values, tested pairs and BH FDR universe. The draft tests all variant–feature combinations and then retains the selected TAD pairs; this can be memory intensive.
7. Compare association coefficients, p-values, adjusted p-values, significant pairs and unique variant totals against the original QTL outputs.
8. Check the output tables, saved PEER model/residuals and logs for patient-level or derived restricted information; do not commit them to GitHub.

The provided synthetic TSVs demonstrate input layouts; they are **not** a powered or end-to-end test dataset, and cannot demonstrate scientific equivalence. Do not claim verified reproducibility until these checks have been completed on the authorized environment.
