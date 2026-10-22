# Validation against original analyses

**Status: not numerically validated.** R execution and real reference/data files were not available during refactoring. Before publication, compare the cleaned outputs against the saved original notebooks and analysis results.

1. RNA-seq: compare gene identifiers and MEC-1/HG-3 TPM values with `LL100_TPM_from_Salmon_MEC1_HG3.xlsx`. Verify HG-3 was included and that live biomaRt versus GENCODE v26 symbols have been reconciled.
2. SEA: compare numbers of common active/inactive oligonucleotides with `Both_active.fa` and `Both_inactive.fa`; check sequence lengths, IDs and duplicates. Reproduce the exact `--order 1` shuffled-control SEA command; do not silently switch to using an inactive control FASTA.
3. SEA annotations: verify motif IDs from HOCOMOCO motif library and the annotation table, uniqueness of TF labels, SEA q-values and number of expressed TFs (>=20 TPM in both lines).
4. FABIAN input: compare the VCF SNP list to original `FABIAN_input_hg38.vcf`; record rows excluded due to missing sequence, indels, allele/reference mismatches or incomplete annotations. Validate every VCF REF base against GRCh38.
5. FABIAN postprocessing: compare combined scores, model preference, >=0.6 cutoff, and gene expression filters (>=30 TPM equivalent category in either cell line). Compare JASPAR family assignments and variant IDs.
6. Archive versions of Salmon, fastp, SRA Toolkit, R/Bioconductor packages, MEME Suite/SEA, FABIAN-variant, HOCOMOCO, JASPAR, and the transcript/GTF reference used in the final study.
