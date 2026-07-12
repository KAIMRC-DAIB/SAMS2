#!/bin/bash
set -euo pipefail

# 07b_humann3_postprocess_v6.sh
# End-to-end HUMAnN postprocess for per-sample subfolders, plus MetaCyc links (ID-only).

BASE="/data/projects/biomage/data/processed/humann3"
OUT="${BASE}/postprocess"
FLAT="${OUT}/tables_flat"

FORCE="${FORCE:-0}"          # FORCE=1 rebuild even if outputs exist
RENAME_GENE="${RENAME_GENE:-0}"  # RENAME_GENE=1 also rename genefamilies to uniref50 names
LINK_BASE="${LINK_BASE:-https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=}"  # or https://biocyc.org/META/new-image?type=PATHWAY&object=

if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
fi
conda activate humann3

mkdir -p "$OUT"

echo "[INFO] BASE=$BASE"
echo "[INFO] OUT=$OUT"
echo "[INFO] FORCE=$FORCE  RENAME_GENE=$RENAME_GENE"
echo "[INFO] LINK_BASE=$LINK_BASE"

need_make () {
  local f="$1"
  if [[ "$FORCE" == "1" ]]; then return 0; fi
  [[ -s "$f" ]] && return 1
  return 0
}

qc_sum () {
  local in="$1"; local out="$2"
  awk -F'\t' '
    NR==1 {for(i=2;i<=NF;i++) h[i]=$i; n=NF; next}
    {for(i=2;i<=n;i++){v=($i==""?0:$i); s[i]+=v}}
    END{
      print "Sample\tSum";
      for(i=2;i<=n;i++) printf "%s\t%.6f\n", h[i], s[i]+0
    }' "$in" > "$out"
}

# 0) completeness
echo "[CHECK] completeness"
ok=0; bad=0
shopt -s nullglob
for d in "$BASE"/SAMS*; do
  s=$(basename "$d")
  if [[ -s "$d/${s}_genefamilies.tsv" && -s "$d/${s}_pathabundance.tsv" && -s "$d/${s}_pathcoverage.tsv" ]]; then
    ok=$((ok+1))
  else
    echo "INCOMPLETE: $s"
    bad=$((bad+1))
  fi
done
echo "[CHECK] COMPLETE=$ok INCOMPLETE=$bad"
if (( bad > 0 )); then
  echo "[ERROR] Some samples incomplete. Fix them before postprocess."
  exit 1
fi

# 1) flat symlinks
echo "[LINK] building flat table folder (symlinks)"
rm -rf "$FLAT"
mkdir -p "$FLAT"

for d in "$BASE"/SAMS*; do
  s=$(basename "$d")
  ln -sf "$d/${s}_genefamilies.tsv" "$FLAT/" || true
  ln -sf "$d/${s}_pathabundance.tsv" "$FLAT/" || true
  ln -sf "$d/${s}_pathcoverage.tsv" "$FLAT/" || true
done

n_gf=$(ls -1 "$FLAT"/*_genefamilies.tsv 2>/dev/null | wc -l)
n_pa=$(ls -1 "$FLAT"/*_pathabundance.tsv 2>/dev/null | wc -l)
n_pc=$(ls -1 "$FLAT"/*_pathcoverage.tsv 2>/dev/null | wc -l)
echo "[LINK] Found genefamilies=$n_gf pathabundance=$n_pa pathcoverage=$n_pc (expect $ok each)"
if (( n_gf != ok || n_pa != ok || n_pc != ok )); then
  echo "[ERROR] Missing some tables in FLAT folder. Check naming."
  exit 1
fi

# 2) join
GF_M="$OUT/genefamilies_merged.tsv"
PA_M="$OUT/pathabundance_merged.tsv"
PC_M="$OUT/pathcoverage_merged.tsv"

if need_make "$GF_M"; then
  echo "[JOIN] genefamilies"
  humann_join_tables --input "$FLAT" --output "$GF_M" --file_name genefamilies
else
  echo "[SKIP] Found $GF_M"
fi

if need_make "$PA_M"; then
  echo "[JOIN] pathabundance"
  humann_join_tables --input "$FLAT" --output "$PA_M" --file_name pathabundance
else
  echo "[SKIP] Found $PA_M"
fi

if need_make "$PC_M"; then
  echo "[JOIN] pathcoverage"
  humann_join_tables --input "$FLAT" --output "$PC_M" --file_name pathcoverage
else
  echo "[SKIP] Found $PC_M"
fi

# 3) renorm
GF_RELAB="$OUT/genefamilies_relab.tsv"
GF_CPM="$OUT/genefamilies_cpm.tsv"
PA_RELAB="$OUT/pathabundance_relab.tsv"
PC_RELAB="$OUT/pathcoverage_relab.tsv"

if need_make "$GF_RELAB"; then
  echo "[RENORM] genefamilies relab"
  humann_renorm_table --input "$GF_M" --output "$GF_RELAB" --units relab
else
  echo "[SKIP] Found $GF_RELAB"
fi

if need_make "$GF_CPM"; then
  echo "[RENORM] genefamilies cpm"
  humann_renorm_table --input "$GF_M" --output "$GF_CPM" --units cpm
else
  echo "[SKIP] Found $GF_CPM"
fi

if need_make "$PA_RELAB"; then
  echo "[RENORM] pathabundance relab"
  humann_renorm_table --input "$PA_M" --output "$PA_RELAB" --units relab
else
  echo "[SKIP] Found $PA_RELAB"
fi

if need_make "$PC_RELAB"; then
  echo "[RENORM] pathcoverage relab"
  humann_renorm_table --input "$PC_M" --output "$PC_RELAB" --units relab
else
  echo "[SKIP] Found $PC_RELAB"
fi

# 4) rename pathways (correct key is metacyc-pwy)
PA_NAMED="$OUT/pathabundance_relab_named.tsv"
if need_make "$PA_NAMED"; then
  echo "[RENAME] pathways (metacyc-pwy)"
  humann_rename_table --input "$PA_RELAB" --output "$PA_NAMED" --names metacyc-pwy
else
  echo "[SKIP] Found $PA_NAMED"
fi

if [[ "$RENAME_GENE" == "1" ]]; then
  GF_NAMED="$OUT/genefamilies_relab_named.tsv"
  if need_make "$GF_NAMED"; then
    echo "[RENAME] genefamilies (uniref50)"
    humann_rename_table --input "$GF_RELAB" --output "$GF_NAMED" --names uniref50
  else
    echo "[SKIP] Found $GF_NAMED"
  fi
fi

# 5) QC
QC_PA="$OUT/qc_pathabundance_relab_colsum.tsv"
QC_GF="$OUT/qc_genefamilies_relab_colsum.tsv"

if need_make "$QC_PA"; then
  echo "[QC] pathabundance relab column sums"
  qc_sum "$PA_RELAB" "$QC_PA"
else
  echo "[SKIP] Found $QC_PA"
fi

if need_make "$QC_GF"; then
  echo "[QC] genefamilies relab column sums"
  qc_sum "$GF_RELAB" "$QC_GF"
else
  echo "[SKIP] Found $QC_GF"
fi

# 6) Build final table with ID-only + link (FIXED2 logic)
PA_LINKED="$OUT/pathabundance_relab_named_with_links.tsv"
if need_make "$PA_LINKED"; then
  echo "[LINKS] building named+ID+link table (ID-only)"
  awk -F'\t' -v BASEURL="$LINK_BASE" '
    NR==FNR{
      if(FNR==1) next
      feat[FNR]=$1
      c=index($1,":")
      if(c>0){
        nm=substr($1,c+1)
        sub(/^[[:space:]]+/,"",nm)
      } else nm=""
      pname[FNR]=nm
      next
    }
    FNR==1{
      printf "feature\tpathway_id\tpathway_name\tmetacyc_link"
      for(i=2;i<=NF;i++) printf "\t%s",$i
      printf "\n"
      next
    }
    {
      raw=$1
      pid=raw
      c=index(raw,":")
      if(c>0) pid=substr(raw,1,c-1)

      enc=pid
      gsub(/\+/,"%2B",enc)

      link=""
      if(pid!="UNMAPPED" && pid!="UNINTEGRATED" && pid!=""){
        link=BASEURL enc
      }

      printf "%s\t%s\t%s\t%s", feat[FNR], pid, pname[FNR], link
      for(i=2;i<=NF;i++) printf "\t%s",$i
      printf "\n"
    }
  ' "$PA_NAMED" "$PA_RELAB" > "$PA_LINKED"
else
  echo "[SKIP] Found $PA_LINKED"
fi

echo "[DONE] Postprocess complete."
echo "Key pathway file (named): $PA_NAMED"
echo "Key pathway file (named+ID+links): $PA_LINKED"
