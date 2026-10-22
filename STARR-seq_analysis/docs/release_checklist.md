# STARR-seq module validation checklist

Use this checklist for this module; licensing, citation and Zenodo release decisions belong to the main repository.

- [ ] Obtain approval to release the four processed input TSV files, or document their access restrictions.
- [ ] Supply the intended pre-QC `MEC1_counts.tsv`, `HG3_counts.tsv`, `input_counts.tsv` and `sequence_annotations.tsv` locally for testing.
- [ ] Verify genome build, genomic coordinate convention, oligonucleotide ID provenance and upstream count processing in [the data dictionary](../data/README.md).
- [ ] Run all scripts from the `STARR-seq_analysis/` directory in a clean checkout.
- [ ] Inspect filtering, MPRAnalyze exclusion, allele-eligibility and unparsed-ID reports.
- [ ] Compare original and refactored counts and exact enhancer statistics (`alpha`, `pval.mad`, `padj`) and allelic statistics (`logFC`, `fdr`) by ID, not only summary totals.
- [ ] Confirm unique significant rsIDs across MEC-1 and HG-3 using `variant_id`.
- [ ] Run `05_record_environment.R`; review and commit the generated `docs/software_environment.txt`.
- [ ] Confirm the versions, scientific methods and data availability statements match the final manuscript.
- [ ] Confirm this module does not contain restricted patient data, credentials or unapproved third-party inputs.
