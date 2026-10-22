# CLL-STARRseq

Code and documentation for analyses of regulatory genetic variation in chronic lymphocytic leukemia (CLL) assessed by allele-specific STARR-seq in the CLL-derived MEC-1 and HG-3 cell lines. The repository includes workflows for identifying active regulatory sequences, measuring allele-specific regulatory effects, integrating molecular quantitative trait loci (QTLs), and assessing transcription factor (TF) motifs and predicted allele-dependent TF binding.

The analyses are organized into independent modules. Each module has its own README describing the required inputs, methods, execution steps, and outputs.

## Repository structure

```text
CLL-STARRseq/
├── STARR-seq_analysis/
│   ├── README.md
│   ├── scripts/
│   ├── data/                  # input formats and locally supplied data
│   └── docs/
├── QTL_analysis/
│   ├── README.md
│   ├── scripts/
│   ├── docs/
│   └── example_data/          # synthetic examples only
├── TF_analysis/
│   ├── README.md
│   ├── RNAseq_processing/
│   ├── scripts/
│   ├── docs/
│   ├── example_data/          # synthetic examples only
│   └── inputs/                # local external inputs, not distributed
├── README.md
└── .gitignore
```

## Analysis modules

| Module | Scope | Documentation |
| --- | --- | --- |
| **STARR-seq analysis** | Oligonucleotide count QC; MPRAnalyze enhancer activity and allelic-effect testing in MEC-1 and HG-3; integration of enhancer activity with allelic effects | [STARR-seq_analysis/README.md](STARR-seq_analysis/README.md) |
| **QTL analysis** | eQTL, mQTL, hQTL and caQTL analyses; phenotype preparation, PEER correction, Matrix eQTL testing, and annotation of variant–feature associations | [QTL_analysis/README.md](QTL_analysis/README.md) |
| **TF analysis** | MEME Suite SEA motif enrichment; FABIAN-variant binding predictions; processing of LL-100 RNA-seq data for MEC-1 and HG-3 TF-expression annotation | [TF_analysis/README.md](TF_analysis/README.md) |

### STARR-seq analysis

The STARR-seq workflow starts from **processed oligonucleotide DNA/RNA count matrices** and construct annotations. It applies input/DNA count filtering, tests enhancer activity using MPRAnalyze, and evaluates allelic effects between tested reference and alternative sequences. Detailed models and significance criteria are documented in [the STARR-seq methods](STARR-seq_analysis/docs/methods.md).

Run from **inside** `STARR-seq_analysis/`:

```bash
cd STARR-seq_analysis
Rscript scripts/01_quality_control.R
Rscript scripts/02_enhancer_activity.R
Rscript scripts/03_allelic_effects.R
Rscript scripts/04_annotate_allelic_effects.R
Rscript scripts/05_record_environment.R
```

The required count matrices and annotations must be supplied locally. See the [STARR-seq data dictionary](STARR-seq_analysis/data/README.md).

### QTL analysis

The QTL workflows evaluate associations between tested variants and gene expression, DNA methylation, H3K27ac and ATAC-seq phenotypes, including region-specific selection of variant–feature pairs. The scripts use PEER correction and Matrix eQTL. Input molecular phenotype data originate from a previously published integrative CLL epigenomics study, while relevant genotype data are subject to the access requirements of their originating datasets. See [QTL data access](QTL_analysis/docs/data_access.md) and [input formats](QTL_analysis/docs/input_formats.md).

Run from **inside** `QTL_analysis/` after configuring the private `QTL_INPUT_DIR`:

```bash
cd QTL_analysis
Rscript scripts/01_prepare_phenotypes.R all
Rscript scripts/02_run_PEER.R all
Rscript scripts/03_prepare_MatrixEQTL_inputs.R all
Rscript scripts/04_run_MatrixEQTL.R all
Rscript scripts/05_annotate_QTL_results.R all
```

The supplied `example_data/` tables are synthetic illustrations of required formats; they are not study measurements.

### Transcription factor analysis

The TF module prepares sequence inputs for **SEA** (MEME Suite), processes motif enrichment, prepares genomic variant inputs for **FABIAN-variant**, and annotates predicted allele-dependent binding changes. It also documents reprocessing of MEC-1 and HG-3 transcriptome data from the **LL-100** panel for TF-expression annotation. The publicly available RNA-seq source is the [DSMZ CellDive LL-100 collection](https://celldive.dsmz.de/rna/ll-100), also available via [ENA PRJEB30312](https://www.ebi.ac.uk/ena/browser/view/PRJEB30312). The analyses use original motif resources from HOCOMOCO and JASPAR collections.

The SEA and FABIAN computations include steps carried out with their respective external tools. Read [the TF module README](TF_analysis/README.md) and [its input documentation](TF_analysis/docs/input_formats.md) for the precise commands, input requirements, settings and expected outputs.

## Data availability and reproducibility

- **STARR-seq:** Processed experimental count matrices, construct annotations and oligonucleotide sequences are not included in this code distribution. Their release is subject to the associated study's data-sharing arrangements.
- **QTL:** Individual-level genotypes and molecular phenotype matrices are not redistributed here. Relevant published data sources, accession numbers and access restrictions are documented in [QTL data access](QTL_analysis/docs/data_access.md).
- **TF:** Public LL-100 sequencing reads can be obtained from DSMZ/ENA. The module does not redistribute raw FASTQ files, motif reference databases or unpublished variant-level results.
- **Example inputs:** Files under `example_data/` are illustrative/synthetic, not experimental observations.

Locally generated result tables, intermediate files and manuscript submission materials are excluded from the repository.

OpenAI's ChatGPT was used to improve code structure, readability, and documentation. The methodology was defined by the authors, who reviewed and tested the resulting scripts.

## Citation

A citation for the associated study and a versioned software archive will be supplied when the manuscript and code release are finalized.

## Contact

For questions regarding the code, analyses, or reproducibility of the results, please contact:

Christelle Dieppois

Aix-Marseille University, INSERM, TAGC, UMR 1090, Marseille, France

Email: christelle.dieppois@inserm.fr
