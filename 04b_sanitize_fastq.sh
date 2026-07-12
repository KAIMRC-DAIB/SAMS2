#!/bin/bash
set -euo pipefail

# BBMap lives in 'biomage' for you
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
fi
conda activate biomage

IN_DIR="/data/projects/biomage/data/processed/deduplicated"
TMP_DIR="/data/projects/biomage/data/processed/deduplicated_sanitized/tmp_clean"
OUT_DIR="/data/projects/biomage/data/processed/deduplicated_sanitized"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"

THREADS=102
shopt -s nullglob
inputs=( "${IN_DIR}"/*_R1.dedup.fastq.gz )
if ((${#inputs[@]}==0)); then
  echo "[ERROR] No inputs in ${IN_DIR} (*_R1.dedup.fastq.gz)"; exit 1
fi

for fq1 in "${inputs[@]}"; do
  fq2="${fq1/_R1/_R2}"
  base="$(basename "$fq1" | cut -d'_' -f1)"

  tmp1="${TMP_DIR}/${base}_R1.clean.fastq.gz"
  tmp2="${TMP_DIR}/${base}_R2.clean.fastq.gz"
  out1="${OUT_DIR}/${base}_R1.sanitized.fastq.gz"
  out2="${OUT_DIR}/${base}_R2.sanitized.fastq.gz"

  log1="${LOG_DIR}/${base}.reformat.R1.log"
  log2="${LOG_DIR}/${base}.reformat.R2.log"
  logr="${LOG_DIR}/${base}.repair.log"

  if [[ -s "$out1" && -s "$out2" ]]; then
    echo "[SKIP] ${base} already sanitized."; continue
  fi

  echo "[SANITIZE:R1] ${base} : $(date)"
  reformat.sh in="$fq1" out="$tmp1" tossbrokenreads=t ow=t threads="$THREADS" &> "$log1"

  echo "[SANITIZE:R2] ${base} : $(date)"
  reformat.sh in="$fq2" out="$tmp2" tossbrokenreads=t ow=t threads="$THREADS" &> "$log2"

  echo "[REPAIR/PAIR-SYNC] ${base} : $(date)"
  repair.sh in="$tmp1" in2="$tmp2" out="$out1" out2="$out2" outs=/dev/null ow=t &> "$logr"

  echo "[DONE] ${base} : $(date)"
done

echo "[CLEANUP] Remove $TMP_DIR if you need space."

