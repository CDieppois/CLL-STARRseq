# QTL data sources and access

The analysis code can be shared without redistributing the source genotype or molecular phenotype matrices. **No real patient-level data are included in this repository.**

## Original publication and processed molecular datasets

Speedy HE, Beekman R, Chapaprieta V et al. *Insight into genetic predisposition to chronic lymphocytic leukemia from integrative epigenomics*. Nature Communications 10, 3615 (2019). DOI: https://doi.org/10.1038/s41467-019-11582-2

The publication's Data availability statement identifies a separate page for **processed matrices**:

https://resources.idibaps.org/paper/insight-into-genetic-predisposition-to-chronic-lymphocytic-leukemia-from-integrative-epigenomics

**Availability note:** This external download page could not be accessed during preparation of this repository documentation (October 2026); file-level download availability was not independently verified. Consult the publication or its authors if it is unavailable.

The R Markdown notebooks that produced the results refer to the following exact local source filenames:

| QTL | Published molecular matrix filename | Additional annotation |
|---|---|---|
| eQTL | `single_probes_all_cases_RMA_expression.txt` | `GPL13667-15572.txt` (GEO GPL13667) |
| caQTL | `atac_consensus_sample_count_corrected.txt` | None beyond peak intervals in source matrix |
| hQTL | `H3K27ac_consensus_sample_count_corrected.txt` | None beyond peak intervals in source matrix |
| mQTL | `cll_methylation_full_matrix_idibaps_norm_Beta_grch38.tsv.gz` | None beyond CpG intervals in source matrix |

GEO GPL13667: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13667

The code does **not** download or republish the source matrices. Researchers may need to verify their access rights and relevant data reuse conditions.

## Underlying EGA datasets

The publication lists the following EGA accessions: `EGAS00000000092`, `EGAD00001004046`, `EGAS00001000272`, and `EGAS00001001911`.

EGA portal: https://ega-archive.org/ . Access to controlled datasets is granted by the responsible Data Access Committee (DAC), subject to the study's conditions. See https://www.ebi.ac.uk/training/online/courses/ega-quick-tour/accessing-the-data-in-the-ega/ . **EGA accessions do not imply that this repository has permission to redistribute sample-level genotypes.**

## Analysis-ready genotype input

`QTL_variant_genotypes.tsv` is a **user-prepared** analysis-ready file derived from authorized genotype data and a list of variants selected for the STARR-seq study. It is not the name of a downloadable EGA file. The public `example_data/QTL_variant_genotypes_example.tsv` only demonstrates the tabular format with invented values.

To construct the actual genotype TSV on an authorized machine, include unique variant IDs, genome REF/ALT, STARR-seq tested REF/ALT, the selected TAD intervals in hg38 **and** the lifted hg19 TAD coordinates needed for expression, and coded diploid sample genotype columns named `CLL1`, `CLL2`, etc. See `input_formats.md`.

## Where to store source files

Keep the original source matrices and your genotype TSV in a local directory **outside the GitHub repository**, and point `QTL_INPUT_DIR` to that directory. The scripts read them; no `restricted/` folder inside GitHub is required.
