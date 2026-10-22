#!/usr/bin/env bash
# Run from TF_analysis/; downloads public ENA/SRA reads and runs Salmon quantification.
# Requires fasterq-dump, fastp, salmon, gzip, wget, and sufficient local disk space.
set -euo pipefail
THREADS="${THREADS:-8}"
WORK_DIR="${TF_RNASEQ_WORK_DIR:-$PWD/local_RNAseq_data}"
REF_DIR="$WORK_DIR/reference"
mkdir -p "$WORK_DIR" "$REF_DIR"
for tool in fasterq-dump fastp salmon wget gzip; do command -v "$tool" >/dev/null || { echo "Missing: $tool" >&2; exit 1; }; done

TX_FA_GZ="$REF_DIR/gencode.v26.transcripts.fa.gz"
GTF_GZ="$REF_DIR/gencode.v26.annotation.gtf.gz"
[[ -s "$TX_FA_GZ" ]] || wget -O "$TX_FA_GZ" 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz'
[[ -s "$GTF_GZ" ]] || wget -O "$GTF_GZ" 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz'
[[ -s "$REF_DIR/gencode.v26.transcripts.fa" ]] || gzip -dc "$TX_FA_GZ" > "$REF_DIR/gencode.v26.transcripts.fa"
[[ -s "$REF_DIR/gencode.v26.annotation.gtf" ]] || gzip -dc "$GTF_GZ" > "$REF_DIR/gencode.v26.annotation.gtf"
INDEX="$REF_DIR/salmon_index_v26"
if [[ ! -d "$INDEX" ]]; then salmon index -t "$REF_DIR/gencode.v26.transcripts.fa" -i "$INDEX" -k 31 --gencode; fi

process_sample() {
  local sample="$1" accession="$2" out="$WORK_DIR/$1"
  mkdir -p "$out"
  (
    cd "$out"
    [[ -s "${accession}_1.fastq" && -s "${accession}_2.fastq" ]] || fasterq-dump "$accession" --split-files --threads "$THREADS"
    fastp -i "${accession}_1.fastq" -I "${accession}_2.fastq" \
      -o "${sample}_trimmed_1.fastq.gz" -O "${sample}_trimmed_2.fastq.gz" \
      -h "${sample}_fastp.html" -j "${sample}_fastp.json" \
      --detect_adapter_for_pe --thread "$THREADS"
  )
  salmon quant -i "$INDEX" -l ISR \
    -1 "$out/${sample}_trimmed_1.fastq.gz" -2 "$out/${sample}_trimmed_2.fastq.gz" \
    --validateMappings --gcBias -p "$THREADS" -o "$out/${sample}_salmon_quant"
}
process_sample MEC1 ERR3003563
process_sample HG3 ERR3003538
printf 'Done. Quantifications: %s/{MEC1,HG3}/*_salmon_quant/quant.sf\n' "$WORK_DIR"
