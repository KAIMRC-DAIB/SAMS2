#!/bin/bash
set -euo pipefail

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
conda activate mpa

MP_OUT="/data/projects/biomage/data/processed/metaphlan"
LOG_DIR="${MP_OUT}/logs"
mkdir -p "$LOG_DIR"

MERGED="${MP_OUT}/metaphlan_merged_abundance.tsv"
CLEAN="${MP_OUT}/metaphlan_merged_abundance.clean.tsv"

# 1) Build merged if missing
shopt -s nullglob
profiles=( "${MP_OUT}"/*_profile.txt )
if ((${#profiles[@]} == 0)); then
  echo "[ERROR] No *_profile.txt in ${MP_OUT}"; exit 1
fi
if [[ ! -s "$MERGED" ]]; then
  echo "[MERGE] ${MERGED}"
  merge_metaphlan_tables.py "${profiles[@]}" > "$MERGED"
else
  echo "[SKIP] Found ${MERGED}"
fi

# 2) Make a clean merged (drop comment lines, keep first non-comment as header)
echo "[CLEAN] ${CLEAN}"
awk 'BEGIN{FS=OFS="\t"}
     /^#/ {next}
     !hdr {print; hdr=1; next}
     {print}' "$MERGED" > "$CLEAN"

# 3) Split by rank from the CLEAN file (keep header + rows whose last rank matches)
echo "[SPLIT] per-rank tables"
in="$CLEAN"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)k__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_kingdom.tsv"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)p__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_phylum.tsv"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)c__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_class.tsv"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)o__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_order.tsv"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)f__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_family.tsv"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)g__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_genus.tsv"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)s__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_species.tsv"
awk -F'\t' 'NR==1||$1 ~ /(^|\|)t__[^|]+$/ {print}' "$in" > "${MP_OUT}/metaphlan_sgb.tsv"

# Convenience copies
cp "${MP_OUT}/metaphlan_species.tsv" "${MP_OUT}/metaphlan_species_only.tsv"
cp "${MP_OUT}/metaphlan_sgb.tsv"     "${MP_OUT}/metaphlan_sgb_only.tsv"

# 4) Proportions versions (percent / 100)
mk_prop () {
  local src="$1"; local dst="$2"
  awk -F'\t' 'NR==1{print; next}{
    printf "%s",$1;
    for(i=2;i<=NF;i++){v=($i==""?0:$i)/100.0; printf "\t%.10f", v}
    printf "\n"
  }' OFS='\t' "$src" > "$dst"
}
echo "[PROP] build *_prop.tsv"
mk_prop "${MP_OUT}/metaphlan_species_only.tsv" "${MP_OUT}/metaphlan_species_only_prop.tsv"
mk_prop "${MP_OUT}/metaphlan_sgb_only.tsv"     "${MP_OUT}/metaphlan_sgb_only_prop.tsv"
mk_prop "${MP_OUT}/metaphlan_genus.tsv"        "${MP_OUT}/metaphlan_genus_prop.tsv"
mk_prop "${MP_OUT}/metaphlan_phylum.tsv"       "${MP_OUT}/metaphlan_phylum_prop.tsv"

# 5) QC: column sums (header-driven, robust)
sum_cols () {
  local src="$1"; local out="$2"
  awk -F'\t' '
    NR==1{for(i=2;i<=NF;i++){hdr[i]=$i}; ncol=NF; next}
    {for(i=2;i<=ncol;i++){v=($i==""?0:$i); sum[i]+=v}}
    END{
      print "Sample\tSum";
      for(i=2;i<=ncol;i++){printf "%s\t%.6f\n", hdr[i], (sum[i]+0)}
    }' "$src" > "$out"
}
echo "[QC] column sums"
sum_cols "${MP_OUT}/metaphlan_species_only.tsv"        "${MP_OUT}/qc_species_colsum.tsv"
sum_cols "${MP_OUT}/metaphlan_species_only_prop.tsv"   "${MP_OUT}/qc_species_colsum_prop.tsv"

# 6) Sample IDs (from CLEAN header)
awk -F'\t' 'NR==1{for(i=2;i<=NF;i++) print $i}' "$CLEAN" > "${MP_OUT}/sample_ids.txt"

echo "[DONE] Wrote per-rank, *_prop, QC, and sample_ids to ${MP_OUT}"


# === Outputs (written to ${MP_OUT}) ===
# metaphlan_merged_abundance.tsv            # merged raw MetaPhlAn table (may start with '#')
# metaphlan_merged_abundance.clean.tsv      # same, comments removed; real header retained
# metaphlan_{kingdom,phylum,class,order,family,genus,species,sgb}.tsv   # per-rank % abundance
# metaphlan_species_only.tsv                # alias of species table (%)
# metaphlan_sgb_only.tsv                    # alias of SGB table (%; t__SGB IDs)
# metaphlan_{phylum,genus}_prop.tsv         # per-rank abundances as proportions (0–1)
# metaphlan_species_only_prop.tsv           # species abundances as proportions (0–1)
# metaphlan_sgb_only_prop.tsv               # SGB abundances as proportions (0–1)
# qc_species_colsum.tsv                     # column sums for species % (~100 per sample)
# qc_species_colsum_prop.tsv                # column sums for species proportions (~1.0 per sample)
# sample_ids.txt                            # sample names in column order

