#!/bin/bash
source activate biomage
mkdir -p /data/projects/biomage/data/processed/filtered
THREADS=102

for fq1 in /data/projects/biomage/data/processed/fastp/*_R1.trimmed.fastq.gz; do
    fq2=${fq1/_R1/_R2}
    base=$(basename $fq1 | cut -d'_' -f1)
    bbduk.sh in1=$fq1 in2=$fq2 \
            out1=/data/projects/biomage/data/processed/filtered/${base}_R1.filtered.fastq.gz \
            out2=/data/projects/biomage/data/processed/filtered/${base}_R2.filtered.fastq.gz \
            entropy=0.7 entropywindow=50 entropyk=5 threads=$THREADS
done

