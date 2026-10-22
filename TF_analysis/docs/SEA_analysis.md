# SEA motif enrichment analysis

## Sequence classification

Read `MEC1_enhancer_activity.tsv` and `HG3_enhancer_activity.tsv` from the adjacent STARR-seq module and join on `element_id` with the actual oligonucleotide sequences from the Twist library design. A sequence is **shared active** if tested and called active in both cell lines, and **shared inactive** if tested and not called active in both cell lines. Mixed-status, missing, and untested elements are not assigned to either shared group. The two FASTA outputs are `Both_active.fa` and `Both_inactive.fa`.

## Actual enrichment comparison

The original analysis used the command:

```bash
sea --oc SEA_Both_active_dinuc --p Both_active.fa --order 1 --m motif_databases/HOCOMOCO/H12CORE_meme_format.meme
```

This is **not an active-versus-inactive sequence comparison** because `--n Both_inactive.fa` is not provided. SEA constructs a shuffled control for every primary sequence, preserving dinucleotide frequencies (order 1 = preserved 2-mers). Changing the command to use `--n` would change the statistical analysis and should not be represented as reproducing the original run.

## Enrichment annotations and filters

The SEA output `sea.tsv` is annotated with the HOCOMOCO motif-to-TF mapping, and LL-100-derived TPM values in MEC-1/HG-3. The postprocessing retains one best motif per TF using smallest SEA q-value (with positive-sequence hits as secondary ranking), then selects TFs with SEA q-value < 0.05 and TPM >= 20 in **both** cell lines. Use the exact motif library in the original `--m` call (`H12CORE_meme_format.meme`); your Rmd also used `H13CORE_motifs.tsv` for annotations, and motif identifiers need to be checked for consistency if those downloads originate from different HOCOMOCO releases.

## Sources

- MEME Suite SEA reference: https://meme-suite.org/meme/doc/sea.html
- HOCOMOCO: https://hocomoco11.autosome.org/ (consult the original versioned HOCOMOCO download used in the analysis rather than silently switching releases)
