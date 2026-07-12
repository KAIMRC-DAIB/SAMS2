#!/bin/bash
set -euo pipefail

# --- robust conda activation ---
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
else
  for c in "$HOME/miniconda3/etc/profile.d/conda.sh" \
           "$HOME/anaconda3/etc/profile.d/conda.sh" \
           "/opt/conda/etc/profile.d/conda.sh"; do
    [ -f "$c" ] && source "$c" && break
  done
fi
conda activate mpa  # MetaPhlAn 4.2.2 env

# --- paths ---
BASE="/data/projects/biomage/data/processed"
IN_DIR="${BASE}/deduplicated_sanitized"     # << use sanitized pairs
OUT_DIR="${BASE}/metaphlan"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

# --- params ---
THREADS=102
READ_MIN_LEN=100
TAX_LEVEL="a"
DB_DIR="/data/projects/biomage/data/db/metaphlan4/vJan25"
OFFLINE="--offline"   # set to "" if you want online checks

# --- gather inputs (sanitized preferred; fallback to dedup if none) ---
shopt -s nullglob
inputs=( "${IN_DIR}"/*_R1.sanitized.fastq.gz )
if ((${#inputs[@]} == 0)); then
  IN_DIR="${BASE}/deduplicated"
  inputs=( "${IN_DIR}"/*_R1.dedup.fastq.gz )
fi
if ((${#inputs[@]} == 0)); then
  echo "[ERROR] No R1 inputs found in ${IN_DIR} (expected *_R1.sanitized.fastq.gz or *_R1.dedup.fastq.gz)."
  exit 1
fi

# --- helper: safer completeness check ---
is_complete() {
  local profile="$1" mapout="$2"
  [[ -s "$profile" ]] || return 1
  grep -qvE '^\s*#' "$profile" || return 1
  [[ -s "$mapout" ]] || return 1
  bzip2 -t "$mapout" 2>/dev/null || return 1
  return 0
}

# --- main loop ---
for fq1 in "${inputs[@]}"; do
  file="$(basename "$fq1")"

  if [[ "$file" == *_R1.sanitized.fastq.gz ]]; then
    fq2="${IN_DIR}/${file/_R1.sanitized.fastq.gz/_R2.sanitized.fastq.gz}"
    base="${file%_R1.sanitized.fastq.gz}"
  else
    fq2="${IN_DIR}/${file/_R1.dedup.fastq.gz/_R2.dedup.fastq.gz}"
    base="${file%_R1.dedup.fastq.gz}"
  fi

  [[ -s "$fq2" ]] || { echo "[WARN] Missing mate for $file; skipping."; continue; }

  out_profile="${OUT_DIR}/${base}_profile.txt"
  out_map="${OUT_DIR}/${base}.mapout.bz2"
  log="${LOG_DIR}/${base}.metaphlan.log"

  # skip if already complete
  if is_complete "$out_profile" "$out_map"; then
    echo "[SKIP] ${base} already complete in ${OUT_DIR}."
    continue
  fi

  echo "[METAPHLAN] ${base} : $(date)"
  metaphlan \
    "${IN_DIR}/${file},${fq2}" \
    --input_type fastq \
    --nproc "$THREADS" \
    --tax_lev "$TAX_LEVEL" \
    --read_min_len "$READ_MIN_LEN" \
    --db_dir "$DB_DIR" \
    $OFFLINE \
    --mapout "$out_map" \
    -o "$out_profile" \
    --sample_id "$base" \
    &> "$log" || { echo "[ERROR] Failed: ${base}. See $log"; continue; }

  if is_complete "$out_profile" "$out_map"; then
    echo "[DONE] ${base} : $(date)"
  else
    echo "[WARN] Incomplete outputs for ${base}. See $log"
  fi
done

