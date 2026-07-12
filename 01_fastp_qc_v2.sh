#!/bin/bash
set -e
source activate biomage

# Input directory as argument
if [ -z "$1" ]; then
    echo "Usage: $0 <input_fastq_directory>"
    exit 1
fi

INPUT_DIR="$1"
OUT_DIR="/data/projects/biomage/data/processed/fastp"
THREADS=102

mkdir -p $OUT_DIR

for fq1 in ${INPUT_DIR}/*_R1*.fastq.gz; do
    fq2=${fq1/_R1/_R2}
    base=$(basename $fq1 | cut -d'_' -f1)
    fastp -i $fq1 -I $fq2 \
          -o ${OUT_DIR}/${base}_R1.trimmed.fastq.gz \
          -O ${OUT_DIR}/${base}_R2.trimmed.fastq.gz \
          --detect_adapter_for_pe --qualified_quality_phred 20 \
          --length_required 50 --thread $THREADS \
          --html ${OUT_DIR}/${base}.html \
          --json ${OUT_DIR}/${base}.json
done

