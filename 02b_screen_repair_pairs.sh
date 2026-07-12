#!/usr/bin/env bash
set -euo pipefail
source activate biomage

IN_DIR="/data/projects/biomage/data/processed/filtered"
LOG_DIR="/data/projects/biomage/data/qc/02b_pairscreen"
mkdir -p "$LOG_DIR"

CPU_TOTAL=120

# Parallelization knobs (you can override when running)
PIGZ_THREADS="${PIGZ_THREADS:-4}"   # threads per file for counting
JOBS="${JOBS:-$(( CPU_TOTAL / PIGZ_THREADS ))}"

# cap jobs to avoid destroying disk IO
if (( JOBS > 24 )); then JOBS=24; fi
if (( JOBS < 1 )); then JOBS=1; fi

MODE="${MODE:-repair}"  # repair | check

SUMMARY_TMP="$LOG_DIR/summary.tmp.tsv"
SUMMARY="$LOG_DIR/summary.tsv"
FAIL_LIST="$LOG_DIR/fail.list"
: > "$SUMMARY_TMP"
: > "$FAIL_LIST"

echo "[02b] IN_DIR=$IN_DIR"
echo "[02b] MODE=$MODE  JOBS=$JOBS  PIGZ_THREADS=$PIGZ_THREADS"
echo -e "sample\tstatus\tr1_reads\tr2_reads\taction\tnote" > "$SUMMARY"

count_reads () {
  local f="$1"
  local lines

  if command -v pigz >/dev/null 2>&1; then
    lines=$(pigz -dc -p "$PIGZ_THREADS" "$f" | wc -l)
  else
    lines=$(gzip -dc "$f" | wc -l)
  fi

  echo $(( lines / 4 ))
}

process_one () {
  local base="$1"
  local r1="$IN_DIR/${base}_R1.filtered.fastq.gz"
  local r2="$IN_DIR/${base}_R2.filtered.fastq.gz"

  local status="OK"
  local action="none"
  local note=""
  local n1="NA"
  local n2="NA"

  # existence
  if [[ ! -f "$r1" ]]; then
    status="MISSING"; note="missing R1"
    echo -e "${base}\t${status}\t${n1}\t${n2}\t${action}\t${note}"
    echo "$base" >> "$FAIL_LIST"
    return 0
  fi
  if [[ ! -f "$r2" ]]; then
    status="MISSING"; note="missing R2"
    echo -e "${base}\t${status}\t${n1}\t${n2}\t${action}\t${note}"
    echo "$base" >> "$FAIL_LIST"
    return 0
  fi

  # gzip integrity
  if ! gzip -t "$r1" 2>/dev/null; then
    status="CORRUPT"; note="R1 not valid gzip"
    echo -e "${base}\t${status}\t${n1}\t${n2}\t${action}\t${note}"
    echo "$base" >> "$FAIL_LIST"
    return 0
  fi
  if ! gzip -t "$r2" 2>/dev/null; then
    status="CORRUPT"; note="R2 not valid gzip"
    echo -e "${base}\t${status}\t${n1}\t${n2}\t${action}\t${note}"
    echo "$base" >> "$FAIL_LIST"
    return 0
  fi

  # read counts
  n1="$(count_reads "$r1")"
  n2="$(count_reads "$r2")"

  if [[ "$n1" -ne "$n2" ]]; then
    status="MISMATCH"
    note="pairs desynced"

    if [[ "$MODE" == "repair" ]]; then
      action="repair"

      local out1="$IN_DIR/${base}_R1.filtered.repaired.fastq.gz"
      local out2="$IN_DIR/${base}_R2.filtered.repaired.fastq.gz"
      local outs="$IN_DIR/${base}.singletons.fastq.gz"

      if ! repair.sh in1="$r1" in2="$r2" out1="$out1" out2="$out2" outs="$outs" >/dev/null 2>&1; then
        status="REPAIR_FAIL"
        action="repair_failed"
        note="repair.sh failed"
        echo -e "${base}\t${status}\t${n1}\t${n2}\t${action}\t${note}"
        echo "$base" >> "$FAIL_LIST"
        return 0
      fi

      # backup originals (don’t overwrite if backups already exist)
      [[ -f "${r1}.bak" ]] || mv -f "$r1" "${r1}.bak"
      [[ -f "${r2}.bak" ]] || mv -f "$r2" "${r2}.bak"

      # put repaired in place
      mv -f "$out1" "$r1"
      mv -f "$out2" "$r2"

      # verify counts after repair
      local nn1 nn2
      nn1="$(count_reads "$r1")"
      nn2="$(count_reads "$r2")"

      if [[ "$nn1" -ne "$nn2" ]]; then
        status="REPAIR_BAD"
        note="still mismatched after repair"
        echo -e "${base}\t${status}\t${nn1}\t${nn2}\t${action}\t${note}"
        echo "$base" >> "$FAIL_LIST"
        return 0
      fi

      status="REPAIRED"
      note="ok after repair"
      echo -e "${base}\t${status}\t${nn1}\t${nn2}\t${action}\t${note}"
      return 0
    fi
  fi

  echo -e "${base}\t${status}\t${n1}\t${n2}\t${action}\t${note}"
}

export -f process_one count_reads
export IN_DIR MODE PIGZ_THREADS FAIL_LIST

# List samples
mapfile -t BASES < <(cd "$IN_DIR" && ls *_R1.filtered.fastq.gz 2>/dev/null | sed 's/_R1.filtered.fastq.gz$//' | sort)

if (( ${#BASES[@]} == 0 )); then
  echo "ERROR: no *_R1.filtered.fastq.gz found in $IN_DIR" >&2
  exit 1
fi

# Run in parallel (no -n1 with -I)
printf "%s\n" "${BASES[@]}" \
  | xargs -P "$JOBS" -I{} bash -c 'process_one "$1"' _ {} \
  >> "$SUMMARY_TMP"

cat "$SUMMARY_TMP" >> "$SUMMARY"

echo
echo "[02b] Done."
echo "[02b] Summary: $SUMMARY"
echo "[02b] Fail list (if any): $FAIL_LIST"

