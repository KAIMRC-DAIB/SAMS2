#!/bin/bash
set -euo pipefail
source activate biomage

IN_DIR="/data/projects/biomage/data/processed/filtered"
OUT_DIR="/data/projects/biomage/data/processed/host_removed"
mkdir -p "$OUT_DIR"

# ~80% of 128 cores
THREADS=102
PIGZ_THREADS=$THREADS   # use same for compression; adjust if you prefer

bt2_idx="/data/projects/biomage/host/hg38"  # Bowtie2 index prefix (no .fa)

compress_file () {
  # compress a file in place with pigz if available, else gzip
  local f="$1"
  if command -v pigz >/dev/null 2>&1; then
    pigz -f -p "$PIGZ_THREADS" "$f"
  else
    gzip -f "$f"
  fi
}

for fq1 in "${IN_DIR}"/*_R1.filtered.fastq.gz; do
  fq2="${fq1/_R1/_R2}"
  base="$(basename "$fq1" | cut -d'_' -f1)"

  sam="${OUT_DIR}/${base}.sam"
  bam="${OUT_DIR}/${base}.bam"
  sbam="${OUT_DIR}/${base}.sorted.bam"
  r1="${OUT_DIR}/${base}_R1.clean.fastq"
  r2="${OUT_DIR}/${base}_R2.clean.fastq"

    # Skip if already completed
  if [[ -s "${r1}.gz" && -s "${r2}.gz" ]]; then
    echo "[SKIP] ${base} already processed"
    continue
  fi

  echo "[HOST-FILTER] ${base} : $(date)"

  # 1) Align to host
  bowtie2 -x "$bt2_idx" -1 "$fq1" -2 "$fq2" --very-sensitive -X 2000 -p "$THREADS" -S "$sam"

  # 2) Keep only pairs where BOTH reads are unmapped -> BAM
  samtools view -@ "$THREADS" -b -f 12 -F 256 "$sam" > "$bam"
  rm -f "$sam"   # free space ASAP

  # 3) Name-sort BAM for bamtofastq
  samtools sort -@ "$THREADS" -n "$bam" -o "$sbam"
  rm -f "$bam"   # free space ASAP

  # 4) Convert to paired FASTQ (uncompressed temp files)
  bedtools bamtofastq -i "$sbam" -fq "$r1" -fq2 "$r2"
  rm -f "$sbam"  # no longer needed

  # 5) Compress this sample’s FASTQs immediately
  compress_file "$r1"
  compress_file "$r2"

  # If you want to SKIP compressing zero-byte files, replace the two lines above with:
  # [[ -s "$r1" ]] && compress_file "$r1" || rm -f "$r1"
  # [[ -s "$r2" ]] && compress_file "$r2" || rm -f "$r2"

  echo "[DONE] ${base} : $(date)"
done

