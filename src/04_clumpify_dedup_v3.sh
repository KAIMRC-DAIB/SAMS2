#!/bin/bash
set -euo pipefail
source activate biomage

IN_DIR="/data/projects/biomage/data/processed/host_removed"
OUT_DIR="/data/projects/biomage/data/processed/deduplicated"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

# ~80% of 128 cores
THREADS=102

# --- Auto-size Java heap (safe defaults) ---
MEM_AVAIL_KB=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
AVAIL_GB=$(( MEM_AVAIL_KB / 1024 / 1024 ))
RESERVE_GB=48
CAPHARD_GB=180
JAVA_GB=$(( AVAIL_GB - RESERVE_GB ))
(( JAVA_GB < 16 )) && JAVA_GB=16
(( JAVA_GB > CAPHARD_GB )) && JAVA_GB=$CAPHARD_GB
JAVA_MEM="-Xmx${JAVA_GB}g"
echo "[INFO] MemAvailable=${AVAIL_GB}G  Using JAVA_MEM=${JAVA_MEM}  (reserve=${RESERVE_GB}G cap=${CAPHARD_GB}G)"

shopt -s nullglob
inputs=( "${IN_DIR}"/*_R1.clean.fastq.gz )
if ((${#inputs[@]} == 0)); then
  echo "[ERROR] No R1 inputs found in ${IN_DIR} (*_R1.clean.fastq.gz)."
  exit 1
fi

for fq1 in "${inputs[@]}"; do
  fq2="${fq1/_R1/_R2}"
  base="$(basename "$fq1" | cut -d'_' -f1)"

  out1="${OUT_DIR}/${base}_R1.dedup.fastq.gz"
  out2="${OUT_DIR}/${base}_R2.dedup.fastq.gz"
  log="${LOG_DIR}/${base}.clumpify.log"

  if [[ -s "$out1" && -s "$out2" ]]; then
    echo "[SKIP] ${base} already deduplicated."
    continue
  fi

  echo "[CLUMPIFY] ${base} : $(date)"
  clumpify.sh ${JAVA_MEM} \
    in="$fq1" in2="$fq2" \
    out="$out1" out2="$out2" \
    dedupe subs=0 threads="${THREADS}" zl=6 ow=t \
    fixjunk=t 2> "$log"

  if [[ ! -s "$out1" || ! -s "$out2" ]]; then
    echo "[WARN] ${base} produced empty dedup outputs. Check ${log}"
  fi
  echo "[DONE] ${base} : $(date)"
done

