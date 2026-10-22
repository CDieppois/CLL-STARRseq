# Public LL-100 RNA-seq reprocessing

The original datasets are from the DSMZ LL-100 RNA-seq panel, viewable via [CellDive](https://celldive.dsmz.de/rna/ll-100) and [ENA project PRJEB30312](https://www.ebi.ac.uk/ena/browser/view/PRJEB30312). The run accessions used were ERR3003563 (MEC-1) and ERR3003538 (HG-3).

Original reported workflow: SRA Toolkit (`fasterq-dump` paired reads); fastp (paired-end adapter detection and trimming); Salmon **1.4.0**, GENCODE **v26** transcript index with `-k 31 --gencode`, paired-end library setting `-l ISR`, `--validateMappings`, `--gcBias`; followed by `tximport(type="salmon", tx2gene=...)` to obtain gene-level abundance TPM and HGNC gene annotations. The script reproduces these stated parameters. `fasterq-dump` creates plain FASTQ; fastp creates the compressed trimmed FASTQ. Original notes contained a commented HG-3 import and inconsistent compression/option syntax, repaired here.

TPM categories are `[0,1): no`, `[1,10): very_low`, `[10,20): low`, `[20,30): moderate`, `[30,100): high`, and `>=100: very_high` (zero is included in `no`).

**Important reproducibility note:** The original notebook used a live `biomaRt::getBM()` query for `hgnc_symbol` and an inner merge with the TPM table. The reorganized script preserves this approach. Because the original Ensembl archive snapshot was not recorded, gene-symbol mappings may change across reruns; compare against the original `LL100_TPM_from_Salmon_MEC1_HG3.xlsx` when available. Neither the original fastp/SRA Toolkit versions nor the original Ensembl archive snapshot were recorded in the notes provided. The script does not claim to replicate all DSMZ downstream normalization steps; it replicates the documented Salmon-to-TPM processing.
