#!/bin/bash
set -euo pipefail

# 06c_metaphlan_for_humann_vJun23.sh
# - Count read pairs per sample (R1 sanitized)
# - Build runlist using CUTOFF_PAIRS
# - Run MetaPhlAn with vJun23 DB only for samples in the runlist
#   (HUMAnN 3.9-compatible taxonomic profiles)

# ======================
# USER SETTINGS (edit or override at runtime)
# ======================
TOTAL_CORES="${TOTAL_CORES:-120}"
CUTOFF_PAIRS="${CUTOFF_PAIRS:-5000000}"        # <-- change anytime

# Counting stage (I/O heavy): total ~ COUNT_JOBS * COUNT_PIGZ_THREADS
COUNT_JOBS="${COUNT_JOBS:-24}"
COUNT_PIGZ_THREADS="${COUNT_PIGZ_THREADS:-4}"

# MetaPhlAn stage (CPU heavy): total ~ MP_JOBS * MP_THREADS
MP_JOBS="${MP_JOBS:-6}"
MP_THREADS="${MP_THREADS:-20}"

READ_MIN_LEN="${READ_MIN_LEN:-100}"
TAX_LEVEL="${TAX_LEVEL:-a}"
OFFLINE="${OFFLINE:---offline}"                 # set to "" if you want online

# ======================
# PATHS
# ======================
BASE="/data/projects/biomage/data/processed"
IN_DIR="${BASE}/deduplicated_sanitized"

QC_TSV="${BASE}/qc_humann_readpairs.tsv"
RUNLIST="${BASE}/humann_runlist_${CUTOFF_PAIRS}.txt"

OUT_DIR="${BASE}/metaphlan_vJun23"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

DB_DIR="/data/projects/biomage/data/db/metaphlan4/vJun23"
INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"

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
conda activate mpa

echo "[CONFIG] TOTAL_CORES=${TOTAL_CORES}"
echo "[CONFIG] CUTOFF_PAIRS=${CUTOFF_PAIRS}"
echo "[CONFIG] COUNT_JOBS=${COUNT_JOBS} COUNT_PIGZ_THREADS=${COUNT_PIGZ_THREADS} (total~$((COUNT_JOBS*COUNT_PIGZ_THREADS)))"
echo "[CONFIG] MP_JOBS=${MP_JOBS} MP_THREADS=${MP_THREADS} (total~$((MP_JOBS*MP_THREADS)))"
echo "[CONFIG] DB_DIR=${DB_DIR} INDEX=${INDEX} OFFLINE='${OFFLINE}'"

# --- sanity: input files ---
shopt -s nullglob
r1s=( "${IN_DIR}"/*_R1.sanitized.fastq.gz )
if ((${#r1s[@]} == 0)); then
  echo "[ERROR] No inputs in ${IN_DIR} (*_R1.sanitized.fastq.gz)"
  exit 1
fi

# --- sanity: DB looks installed ---
ls "${DB_DIR}/${INDEX}"*.pkl >/dev/null 2>&1 || {
  echo "[ERROR] Missing ${INDEX}*.pkl in ${DB_DIR} (DB not installed/extracted)."
  exit 1
}
ls "${DB_DIR}/${INDEX}".*.bt2l >/dev/null 2>&1 || {
  echo "[ERROR] Missing ${INDEX}.*.bt2l in ${DB_DIR} (Bowtie2 index missing)."
  exit 1
}

# ---------------------------
# Helper: completeness check like 06
# ---------------------------
is_complete () {
  local profile="$1" mapout="$2"
  [[ -s "$profile" ]] || return 1
  grep -qvE '^\s*#' "$profile" || return 1
  [[ -s "$mapout" ]] || return 1
  bzip2 -t "$mapout" 2>/dev/null || return 1
  return 0
}

# ---------------------------
# 1) COUNT READ PAIRS (R1)
# ---------------------------
echo "[QC] Counting read pairs from R1 files..."

export COUNT_PIGZ_THREADS
count_one () {
  local r1="$1"
  local base n
  base="$(basename "$r1" | sed 's/_R1\.sanitized\.fastq\.gz$//')"

  # Robust counting:
  # - try pigz (fast)
  # - if pigz complains, fallback to gzip
  if command -v pigz >/dev/null 2>&1; then
    if n=$(( $(pigz -dc -p "$COUNT_PIGZ_THREADS" "$r1" 2>/dev/null | wc -l) / 4 )); then
      :
    else
      n=$(( $(gzip -dc "$r1" | wc -l) / 4 ))
    fi
  else
    n=$(( $(gzip -dc "$r1" | wc -l) / 4 ))
  fi

  printf "%s\t%s\n" "$base" "$n"
}
export -f count_one

printf '%s\0' "${r1s[@]}" \
  | xargs -0 -P "$COUNT_JOBS" -I{} bash -c 'count_one "$1"' _ {} \
  | sort -k2,2nr > "$QC_TSV"

echo "[QC] Wrote: $QC_TSV"
echo "[QC] Top 5:"
head -n 5 "$QC_TSV"

# ---------------------------
# 2) BUILD RUNLIST (>= cutoff)
# ---------------------------
awk -F'\t' -v C="$CUTOFF_PAIRS" '$2+0 >= C {print $1}' "$QC_TSV" > "$RUNLIST"
echo "[RUNLIST] Cutoff=$CUTOFF_PAIRS pairs"
echo "[RUNLIST] Wrote: $RUNLIST"
echo "[RUNLIST] N passing: $(wc -l < "$RUNLIST")"
echo "[RUNLIST] First 10:"
head "$RUNLIST"

# ---------------------------
# 3) RUN MetaPhlAn vJun23 (only runlist samples)
# ---------------------------
run_one () {
  local base="$1"
  local fq1 fq2 out_profile out_map log

  fq1="${IN_DIR}/${base}_R1.sanitized.fastq.gz"
  fq2="${IN_DIR}/${base}_R2.sanitized.fastq.gz"
  [[ -s "$fq1" && -s "$fq2" ]] || { echo "[WARN] Missing FASTQs for $base; skipping."; return 0; }

  out_profile="${OUT_DIR}/${base}_profile.txt"
  out_map="${OUT_DIR}/${base}.mapout.bz2"
  log="${LOG_DIR}/${base}.metaphlan.log"

  if is_complete "$out_profile" "$out_map"; then
    echo "[SKIP] ${base} already complete in ${OUT_DIR}."
    return 0
  fi

  echo "[MP-vJun23] ${base} : $(date)"
  metaphlan \
    "${fq1},${fq2}" \
    --input_type fastq \
    --nproc "$MP_THREADS" \
    --tax_lev "$TAX_LEVEL" \
    --read_min_len "$READ_MIN_LEN" \
    --db_dir "$DB_DIR" \
    -x "$INDEX" \
    $OFFLINE \
    --mapout "$out_map" \
    -o "$out_profile" \
    --sample_id "$base" \
    &> "$log" || { echo "[ERROR] Failed: ${base}. See $log"; return 0; }

  if is_complete "$out_profile" "$out_map"; then
    echo "[DONE] ${base} : $(date)"
  else
    echo "[WARN] Incomplete outputs for ${base}. See $log"
  fi
}

export -f run_one is_complete
export IN_DIR OUT_DIR LOG_DIR DB_DIR INDEX OFFLINE MP_THREADS READ_MIN_LEN TAX_LEVEL

echo "[MP] Running MetaPhlAn vJun23 for samples in runlist..."
tr -d '\r' < "$RUNLIST" | xargs -P "$MP_JOBS" -I{} bash -c 'run_one "$1"' _ {}

echo "[DONE] vJun23 profiles in: $OUT_DIR"
echo "[DONE] Runlist saved to: $RUNLIST"
