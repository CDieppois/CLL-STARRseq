# Input file formats and locations

Run commands from `TF_analysis/`. The program intentionally does not include original RNA-seq FASTQ, oligonucleotide design, motif database downloads, online-tool results, or final patient-level/clinical information. Local input filenames are shown below; paths can be overridden where specified.

| Local file | Required fields or format | Original source |
|---|---|---|
| `inputs/oligo_sequences.tsv` | `element_id`, `sequence` (one uppercase or lowercase DNA string per oligonucleotide) | Twist/STARR-seq library design |
| `inputs/variant_positions_hg38.tsv` | `variant_core`, `variant_id`, `chromosome`, `position_hg38`, `genomic_ref_hg38` | Validated GRCh38 variant annotation and reference FASTA |
| `inputs/motif_annotations/H12CORE_meme_format.meme` | MEME motif format | Versioned HOCOMOCO download used for SEA |
| `inputs/motif_annotations/H13CORE_motifs.tsv` | HOCOMOCO motif ID, human gene, TF family | HOCOMOCO SEA annotation table from original notebook |
| `inputs/motif_annotations/H14CORE_motifs.tsv` | HOCOMOCO motif ID, human gene, TF family | HOCOMOCO FABIAN annotation table from original notebook |
| `inputs/motif_annotations/ultimate_metadata_table_CORE.tsv` | JASPAR `tax_id`, `name`, `class`, `family` | JASPAR CORE motif metadata |
| `inputs/SEA/sea.tsv` | Raw SEA output file | SEA online or command-line analysis |
| `inputs/FABIAN/fabian_KnownTFBS.data` | Original 14-column TSV-like FABIAN output without header | FABIAN-variant online analysis |

The original STARR-seq integrated allelic effects and enhancer activity results are read from `../STARR-seq_analysis/results/` (the neighboring analysis module); original private/data source locations are not hard-coded inside TF scripts.

Override local input paths without changing code using `TF_SEQUENCE_FILE`, `TF_VARIANT_POSITIONS`, `TF_SEA_RESULTS`, `TF_SEA_MOTIFS`, `TF_FABIAN_RESULTS`, `TF_FABIAN_MOTIFS`, `TF_JASPAR_META`, `TF_EXPRESSION_FILE`, and `TF_RNASEQ_WORK_DIR`.

`example_data/` contains synthetic illustrative tables, not scientific inputs or a full end-to-end test dataset.
