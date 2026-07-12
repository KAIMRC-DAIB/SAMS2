#!/bin/bash
set -euo pipefail

# 07_humann3_run_5M_vJun23.sh
# Run HUMAnN 3.9 ONLY on samples in RUNLIST, using MetaPhlAn vJun23 profiles.
# Interleaves paired FASTQs on-the-fly.

# ======================
# USER SETTINGS (override at runtime if you want)
# ======================
TOTAL_CORES="${TOTAL_CORES:-120}"
THREADS="${THREADS:-24}"              # threads per HUMAnN sample
JOBS="${JOBS:-5}"                     # parallel samples (JOBS*THREADS ~= total cores)

RUNLIST="${RUNLIST:-/data/projects/biomage/data/processed/humann_runlist_5000000.txt}"

# ======================
# PATHS
# ======================
BASE="/data/projects/biomage/data/processed"
IN_DIR="${BASE}/deduplicated_sanitized"
TPROF_DIR="${BASE}/metaphlan_vJun23"          # vJun23 profiles from 06c
OUT_DIR="${BASE}/humann3"
LOG_DIR="${OUT_DIR}/logs"
TMP_BASE="${TMPDIR:-/tmp}/humann3_interleaved"
mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_BASE"

# --- conda ---
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
else
  for c in "$HOME/miniconda3/etc/profile.d/conda.sh" \
           "$HOME/anaconda3/etc/profile.d/conda.sh" \
           "/opt/conda/etc/profile.d/conda.sh"; do
    [ -f "$c" ] && source "$c" && break
  done
fi
conda activate humann3

HUMANN_BIN="$(which humann)"
REFORMAT_BIN="$(which reformat.sh || true)"
[[ -x "$HUMANN_BIN" ]] || { echo "[ERROR] humann not found in humann3 env"; exit 1; }
[[ -x "$REFORMAT_BIN" ]] || { echo "[ERROR] reformat.sh not found (install bbmap in humann3)"; exit 1; }

# Ensure runlist exists
[[ -s "$RUNLIST" ]] || { echo "[ERROR] RUNLIST missing/empty: $RUNLIST"; exit 1; }

echo "[CONFIG] TOTAL_CORES=$TOTAL_CORES  JOBS=$JOBS  THREADS=$THREADS  (total~$((JOBS*THREADS)))"
echo "[CONFIG] RUNLIST=$RUNLIST"
echo "[CONFIG] IN_DIR=$IN_DIR"
echo "[CONFIG] TPROF_DIR=$TPROF_DIR"
echo "[CONFIG] OUT_DIR=$OUT_DIR"
echo "[CONFIG] TMP_BASE=$TMP_BASE"

is_complete () {
  local outdir="$1" base="$2"
  [[ -s "${outdir}/${base}_genefamilies.tsv" ]] && \
  [[ -s "${outdir}/${base}_pathabundance.tsv" ]] && \
  [[ -s "${outdir}/${base}_pathcoverage.tsv"  ]]
}

run_one () {
  local base="$1"

  local r1="${IN_DIR}/${base}_R1.sanitized.fastq.gz"
  local r2="${IN_DIR}/${base}_R2.sanitized.fastq.gz"
  local prof="${TPROF_DIR}/${base}_profile.txt"

  [[ -s "$r1" && -s "$r2" ]] || { echo "[WARN] Missing FASTQs for $base; skip"; return 0; }
  [[ -s "$prof" ]] || { echo "[WARN] Missing vJun23 profile for $base ($prof); skip"; return 0; }

  local sample_out="${OUT_DIR}/${base}"
  mkdir -p "$sample_out"

  local log="${LOG_DIR}/${base}.humann.log"
  local ilv="${TMP_BASE}/${base}.interleaved.fastq.gz"

  if is_complete "$sample_out" "$base"; then
    echo "[SKIP] $base already complete."
    return 0
  fi

  echo "[INTERLEAVE] $base : $(date)"
  "$REFORMAT_BIN" in1="$r1" in2="$r2" out="$ilv" interleaved=t ow=t threads="$THREADS" &>> "$log" \
    || { echo "[ERROR] Interleave failed: $base (see $log)"; rm -f "$ilv"; return 0; }

  echo "[HUMANN] $base : $(date)"
  "$HUMANN_BIN" \
    --input "$ilv" \
    --input-format fastq \
    --threads "$THREADS" \
    --taxonomic-profile "$prof" \
    --output "$sample_out" \
    --output-basename "$base" \
    &>> "$log" || { echo "[ERROR] HUMAnN failed: $base (see $log)"; rm -f "$ilv"; return 0; }

  rm -f "$ilv"

  if is_complete "$sample_out" "$base"; then
    echo "[DONE] $base"
  else
    echo "[WARN] Incomplete outputs for $base (see $log)"
  fi
}

export -f run_one is_complete
export IN_DIR TPROF_DIR OUT_DIR LOG_DIR TMP_BASE THREADS HUMANN_BIN REFORMAT_BIN

# Run in parallel
tr -d '\r' < "$RUNLIST" | xargs -P "$JOBS" -I{} bash -c 'run_one "$1"' _ {}

echo "[DONE] HUMAnN finished. Outputs in: $OUT_DIR"
