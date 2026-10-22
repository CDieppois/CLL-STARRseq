# FABIAN-variant analysis

## Variant selection and genomic input

From both STARR-seq cell lines, select allele-comparison contrasts with a significant allelic effect (`allelic_effect_call`) and enhancer activity in either line (`sig_enhancer`). Deduplicate by `contrast_id` across cell lines, then derive the reference and alternative alleles from the matched oligonucleotide sequences. Support both `_alt_` and `_alt1_` / `_alt2_` labels, and require one and only one base substitution with at least 30 flanking nucleotides either side, matching the original preparation. Output VCF v4.3 on **GRCh38/hg38** only for eligible SNPs.

VCF REF must be the reference genome allele, not merely the allele named `ref` in an assay. Therefore the program requires an explicitly validated `genomic_ref_hg38` base in `variant_positions_hg38.tsv` and excludes discrepancies. Document how that base was obtained from GRCh38. No genotype data are included in the VCF; the sample `0/1` is an example/dummy field required by the original FABIAN submission format, **not** a patient's genotype.

## Online submission and postprocessing

Submit the VCF to FABIAN-variant and save the raw `fabian_KnownTFBS.data` result, with the version and online settings used. The provided postprocessing computes the average model score per variant–TF pair, preferring TFFM models when available and falling back to PWMs otherwise. Retain abs(combined_score) >= 0.6. Annotate HOCOMOCO v14 motif/gene mappings, then gene-level expression. Retain TFs classified as `high` or `very_high` in MEC-1 or HG-3, corresponding to TPM >= 30 from the original gene-expression categories. JASPAR CORE human metadata are used to annotate TF families/classes.

This is a motif-binding **prediction**, not experimental evidence of TF occupancy. Combined scores from different model versions may not be comparable.
