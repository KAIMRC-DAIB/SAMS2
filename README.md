# SAMS2_Metagenomics

## Healthy Ageing Shotgun Metagenomics Project to compare between Healthy vs Resilient Vs Unhealthy groups of individuals aged 55 and above

This repository contains the shotgun metagenomics workflow used for the **SAMS2 healthy ageing gut microbiome project** 

The workflow processes paired-end shotgun metagenomic FASTQ files from raw reads through preprocessing, host-read removal, taxonomic profiling, functional profiling, quality-control summaries, and downstream visualization-ready abundance tables.

## Important Notes

- The final taxonomic profiling workflow uses **MetaPhlAn4**.
- The script `05_kraken2_bracken_v2.sh` is present in the repository, but Kraken2/Bracken was **not used in the final analysis for this study**.
- Most paths in this README reflect the local project layout used for this analysis:

```bash
/data/projects/biomage
```

If the repository is cloned or moved elsewhere, update the paths inside the scripts before running the workflow.

---

## Repository Contents

Scripts are located in:

```bash
/data/projects/biomage/src
```

Current contents:

```text
01_fastp_qc_v2.sh
02_bbduk_low_complexity.sh
02b_screen_repair_pairs.sh
03_bowtie2_host_removal_v3.sh
04_clumpify_dedup_v3.sh
04b_sanitize_fastq.sh
05_kraken2_bracken_v2.sh
06_metaphlan_v2.sh
06b_metaphlan_postprocess.sh
06c_metaphlan_for_humann_vJun23.sh
07_humann3_run.sh
07_humann3_run_v3_serial.sh
07_humann3_run_5M_vJun23.sh
07b_humann3_postprocess_v6.sh
07c_humann3_plots.py
rename_fastq.sh
R/
old_script/
README.md
```

---

## Project Directory Structure

```text
/data/projects/biomage/
├── src/
│   ├── 01_fastp_qc_v2.sh
│   ├── 02_bbduk_low_complexity.sh
│   ├── 02b_screen_repair_pairs.sh
│   ├── 03_bowtie2_host_removal_v3.sh
│   ├── 04_clumpify_dedup_v3.sh
│   ├── 04b_sanitize_fastq.sh
│   ├── 05_kraken2_bracken_v2.sh
│   ├── 06_metaphlan_v2.sh
│   ├── 06b_metaphlan_postprocess.sh
│   ├── 06c_metaphlan_for_humann_vJun23.sh
│   ├── 07_humann3_run.sh
│   ├── 07_humann3_run_v3_serial.sh
│   ├── 07_humann3_run_5M_vJun23.sh
│   ├── 07b_humann3_postprocess_v6.sh
│   ├── 07c_humann3_plots.py
│   ├── rename_fastq.sh
│   ├── R/
│   └── old_script/
│
├── data/
│   ├── raw/                         # Raw FASTQ files
│   ├── processed/
│   │   ├── fastp/                   # fastp-trimmed reads and reports
│   │   ├── filtered/                # Low-complexity filtered reads
│   │   ├── host_removed/            # Non-human reads after host removal
│   │   ├── deduplicated/            # Deduplicated FASTQ files
│   │   ├── deduplicated_sanitized/  # Sanitized and pair-repaired FASTQ files
│   │   ├── metaphlan/               # Main MetaPhlAn taxonomic profiles
│   │   ├── metaphlan_vJun23/        # HUMAnN-compatible MetaPhlAn profiles
│   │   └── humann3/                 # HUMAnN output folders
│   │       └── postprocess/         # HUMAnN merged and normalized tables
│   │
│   └── qc/
│       └── 02b_pairscreen/          # Pair-screening summaries
│
└── host/
    └── hg38                         # Bowtie2 human genome index prefix
```

---

## Software Requirements

The workflow uses several Conda environments.

### `biomage`

Used for preprocessing, read filtering, host-read removal, duplicate removal, and FASTQ repair.

Required tools:

```text
fastp
bbmap
bbduk.sh
repair.sh
reformat.sh
clumpify.sh
bowtie2
samtools
bedtools
pigz
gzip
```

### `mpa`

Used for MetaPhlAn taxonomic profiling and MetaPhlAn table postprocessing.

Required tools:

```text
MetaPhlAn 4.2.2
merge_metaphlan_tables.py
Bowtie2
bzip2
```

Main MetaPhlAn database:

```text
/data/projects/biomage/data/db/metaphlan4/vJan25
```

HUMAnN-compatible MetaPhlAn database:

```text
/data/projects/biomage/data/db/metaphlan4/vJun23
```

HUMAnN-compatible MetaPhlAn index:

```text
mpa_vJun23_CHOCOPhlAnSGB_202403
```

### `humann3`

Used for HUMAnN3 functional profiling and HUMAnN postprocessing.

Required tools:

```text
humann
humann_join_tables
humann_renorm_table
humann_rename_table
reformat.sh
```

### Python

Used for HUMAnN plotting.

Required Python packages:

```text
python >= 3.8
pandas
numpy
matplotlib
```

---

## Workflow Summary

| Step | Script | Purpose | Main tools |
|---:|---|---|---|
| 0 | `rename_fastq.sh` | Optional FASTQ renaming before analysis | Bash |
| 1 | `01_fastp_qc_v2.sh` | Adapter trimming and read quality filtering | fastp |
| 2 | `02_bbduk_low_complexity.sh` | Remove low-complexity reads | BBDuk |
| 3 | `02b_screen_repair_pairs.sh` | Check, repair, and synchronize paired reads | BBMap `repair.sh` |
| 4 | `03_bowtie2_host_removal_v3.sh` | Remove human host reads | Bowtie2, SAMtools, BEDTools |
| 5 | `04_clumpify_dedup_v3.sh` | Remove duplicate reads | Clumpify |
| 6 | `04b_sanitize_fastq.sh` | Remove broken FASTQ records and repair pairs | BBMap `reformat.sh`, `repair.sh` |
| 7 | `06_metaphlan_v2.sh` | Main taxonomic profiling | MetaPhlAn4 |
| 8 | `06b_metaphlan_postprocess.sh` | Merge and split MetaPhlAn tables | MetaPhlAn utilities, awk |
| 9 | `06c_metaphlan_for_humann_vJun23.sh` | Create HUMAnN-compatible MetaPhlAn profiles | MetaPhlAn vJun23 |
| 10 | `07_humann3_run_5M_vJun23.sh` | Run HUMAnN3 on selected samples | HUMAnN3 |
| 11 | `07b_humann3_postprocess_v6.sh` | Merge, normalize, rename, and QC HUMAnN tables | HUMAnN utilities |
| 12 | `07c_humann3_plots.py` | Generate HUMAnN QC and pathway plots | Python |

---

## Recommended Run Order

Start from the source directory:

```bash
cd /data/projects/biomage/src
```

Run the scripts in order:

```bash
# Optional: rename raw FASTQ files if needed
bash rename_fastq.sh

# 1. Adapter trimming and read quality filtering
bash 01_fastp_qc_v2.sh /path/to/raw_fastq_directory

# 2. Low-complexity filtering
bash 02_bbduk_low_complexity.sh

# 3. Pair checking and repair
bash 02b_screen_repair_pairs.sh

# 4. Human host-read removal
bash 03_bowtie2_host_removal_v3.sh

# 5. Duplicate removal
bash 04_clumpify_dedup_v3.sh

# 6. FASTQ sanitization and pair repair
bash 04b_sanitize_fastq.sh

# 7. Main MetaPhlAn taxonomic profiling
bash 06_metaphlan_v2.sh

# 8. MetaPhlAn postprocessing
bash 06b_metaphlan_postprocess.sh

# 9. HUMAnN-compatible MetaPhlAn profiling using the vJun23 database
bash 06c_metaphlan_for_humann_vJun23.sh

# 10. HUMAnN3 functional profiling
bash 07_humann3_run_5M_vJun23.sh

# 11. HUMAnN3 postprocessing
bash 07b_humann3_postprocess_v6.sh

# 12. HUMAnN3 plotting
python 07c_humann3_plots.py
```

---

## Detailed Workflow

### 0. Optional FASTQ Renaming

**Script:** `rename_fastq.sh`

Use this script only if raw FASTQ files need to be renamed into the format expected by the pipeline.

Most downstream scripts expect paired-end files to contain `_R1` and `_R2` in their names.

Expected format:

```text
SAMPLE_R1.fastq.gz
SAMPLE_R2.fastq.gz
```

Run only if needed:

```bash
bash rename_fastq.sh
```

---

### 1. Adapter Trimming and Quality Filtering

**Script:** `01_fastp_qc_v2.sh`

This step uses `fastp` for adapter detection, adapter trimming, quality filtering, and length filtering.

Command:

```bash
bash 01_fastp_qc_v2.sh /path/to/raw_fastq_directory
```

Main parameters:

```bash
--detect_adapter_for_pe
--qualified_quality_phred 20
--length_required 50
--thread 102
```

Input:

```text
/path/to/raw_fastq_directory/*_R1*.fastq.gz
/path/to/raw_fastq_directory/*_R2*.fastq.gz
```

Output:

```text
/data/projects/biomage/data/processed/fastp/
├── SAMPLE_R1.trimmed.fastq.gz
├── SAMPLE_R2.trimmed.fastq.gz
├── SAMPLE.html
└── SAMPLE.json
```

The `.html` and `.json` files are per-sample `fastp` quality-control reports.

---

### 2. Low-Complexity Filtering

**Script:** `02_bbduk_low_complexity.sh`

This step removes low-complexity reads using `bbduk.sh`.

Command:

```bash
bash 02_bbduk_low_complexity.sh
```

Main parameters:

```bash
entropy=0.7
entropywindow=50
entropyk=5
threads=102
```

Input:

```text
/data/projects/biomage/data/processed/fastp/*_R1.trimmed.fastq.gz
/data/projects/biomage/data/processed/fastp/*_R2.trimmed.fastq.gz
```

Output:

```text
/data/projects/biomage/data/processed/filtered/
├── SAMPLE_R1.filtered.fastq.gz
└── SAMPLE_R2.filtered.fastq.gz
```

---

### 3. Paired-End Screening and Repair

**Script:** `02b_screen_repair_pairs.sh`

This step checks whether paired-end FASTQ files are synchronized.

The script checks:

- whether R1 exists
- whether R2 exists
- whether each file is a valid gzip file
- whether R1 and R2 have the same number of reads
- whether mismatched pairs can be repaired

By default, the script runs in repair mode:

```bash
MODE=repair
```

Command:

```bash
bash 02b_screen_repair_pairs.sh
```

Check-only mode:

```bash
MODE=check bash 02b_screen_repair_pairs.sh
```

Control parallelization:

```bash
PIGZ_THREADS=4 JOBS=24 bash 02b_screen_repair_pairs.sh
```

Input:

```text
/data/projects/biomage/data/processed/filtered/*_R1.filtered.fastq.gz
/data/projects/biomage/data/processed/filtered/*_R2.filtered.fastq.gz
```

Output:

```text
/data/projects/biomage/data/qc/02b_pairscreen/
├── summary.tsv
├── summary.tmp.tsv
└── fail.list
```

Summary columns:

```text
sample
status
r1_reads
r2_reads
action
note
```

Possible statuses:

```text
OK
MISSING
CORRUPT
MISMATCH
REPAIRED
REPAIR_FAIL
REPAIR_BAD
```

If repair is successful, original files are backed up as:

```text
SAMPLE_R1.filtered.fastq.gz.bak
SAMPLE_R2.filtered.fastq.gz.bak
```

---

### 4. Human Host-Read Removal

**Script:** `03_bowtie2_host_removal_v3.sh`

This step aligns paired reads against the human reference genome and keeps only pairs where both mates are unmapped.

Human Bowtie2 index prefix:

```text
/data/projects/biomage/host/hg38
```

Command:

```bash
bash 03_bowtie2_host_removal_v3.sh
```

Main Bowtie2 parameters:

```bash
--very-sensitive
-X 2000
-p 102
```

The script keeps pairs where both reads are unmapped:

```bash
samtools view -@ "$THREADS" -b -f 12 -F 256
```

Input:

```text
/data/projects/biomage/data/processed/filtered/*_R1.filtered.fastq.gz
/data/projects/biomage/data/processed/filtered/*_R2.filtered.fastq.gz
```

Output:

```text
/data/projects/biomage/data/processed/host_removed/
├── SAMPLE_R1.clean.fastq.gz
└── SAMPLE_R2.clean.fastq.gz
```

Temporary SAM and BAM files are removed after each sample to save disk space.

---

### 5. Duplicate Removal

**Script:** `04_clumpify_dedup_v3.sh`

This step removes duplicate reads using `clumpify.sh`.

The script automatically estimates available system memory and selects a Java heap size while reserving memory for the system.

Command:

```bash
bash 04_clumpify_dedup_v3.sh
```

Main parameters:

```bash
dedupe
subs=0
threads=102
zl=6
fixjunk=t
```

Input:

```text
/data/projects/biomage/data/processed/host_removed/*_R1.clean.fastq.gz
/data/projects/biomage/data/processed/host_removed/*_R2.clean.fastq.gz
```

Output:

```text
/data/projects/biomage/data/processed/deduplicated/
├── SAMPLE_R1.dedup.fastq.gz
├── SAMPLE_R2.dedup.fastq.gz
└── logs/
    └── SAMPLE.clumpify.log
```

---

### 6. FASTQ Sanitization and Pair Repair

**Script:** `04b_sanitize_fastq.sh`

This step removes broken FASTQ records and repairs paired-end synchronization after deduplication.

Although this script is named `04b`, it should be run **after** `04_clumpify_dedup_v3.sh`.

Command:

```bash
bash 04b_sanitize_fastq.sh
```

Main parameters:

```bash
tossbrokenreads=t
ow=t
threads=102
```

Input:

```text
/data/projects/biomage/data/processed/deduplicated/*_R1.dedup.fastq.gz
/data/projects/biomage/data/processed/deduplicated/*_R2.dedup.fastq.gz
```

Output:

```text
/data/projects/biomage/data/processed/deduplicated_sanitized/
├── SAMPLE_R1.sanitized.fastq.gz
├── SAMPLE_R2.sanitized.fastq.gz
├── tmp_clean/
└── logs/
    ├── SAMPLE.reformat.R1.log
    ├── SAMPLE.reformat.R2.log
    └── SAMPLE.repair.log
```

Temporary files are written to:

```text
/data/projects/biomage/data/processed/deduplicated_sanitized/tmp_clean
```

After checking the outputs, the temporary directory can be removed if disk space is needed:

```bash
rm -rf /data/projects/biomage/data/processed/deduplicated_sanitized/tmp_clean
```

---

### 7. Taxonomic Profiling with MetaPhlAn4

**Script:** `06_metaphlan_v2.sh`

This step performs taxonomic profiling using MetaPhlAn4.

The script uses sanitized FASTQ files when available. If sanitized files are missing, it falls back to deduplicated FASTQ files.

Main MetaPhlAn database:

```text
/data/projects/biomage/data/db/metaphlan4/vJan25
```

Command:

```bash
bash 06_metaphlan_v2.sh
```

Main parameters:

```bash
--input_type fastq
--nproc 102
--tax_lev a
--read_min_len 100
--offline
```

Preferred input:

```text
/data/projects/biomage/data/processed/deduplicated_sanitized/*_R1.sanitized.fastq.gz
/data/projects/biomage/data/processed/deduplicated_sanitized/*_R2.sanitized.fastq.gz
```

Fallback input:

```text
/data/projects/biomage/data/processed/deduplicated/*_R1.dedup.fastq.gz
/data/projects/biomage/data/processed/deduplicated/*_R2.dedup.fastq.gz
```

Output:

```text
/data/projects/biomage/data/processed/metaphlan/
├── SAMPLE_profile.txt
├── SAMPLE.mapout.bz2
└── logs/
    └── SAMPLE.metaphlan.log
```

The script checks output completeness before skipping previously processed samples.

---

### 8. MetaPhlAn Postprocessing

**Script:** `06b_metaphlan_postprocess.sh`

This step merges all MetaPhlAn sample profiles into one abundance table, removes comment lines, splits the merged table by taxonomic rank, converts percentage abundances to proportions, and generates QC files.

Command:

```bash
bash 06b_metaphlan_postprocess.sh
```

Input:

```text
/data/projects/biomage/data/processed/metaphlan/*_profile.txt
```

Main outputs:

```text
/data/projects/biomage/data/processed/metaphlan/
├── metaphlan_merged_abundance.tsv
├── metaphlan_merged_abundance.clean.tsv
├── metaphlan_kingdom.tsv
├── metaphlan_phylum.tsv
├── metaphlan_class.tsv
├── metaphlan_order.tsv
├── metaphlan_family.tsv
├── metaphlan_genus.tsv
├── metaphlan_species.tsv
├── metaphlan_sgb.tsv
├── metaphlan_species_only.tsv
├── metaphlan_sgb_only.tsv
├── metaphlan_species_only_prop.tsv
├── metaphlan_sgb_only_prop.tsv
├── metaphlan_genus_prop.tsv
├── metaphlan_phylum_prop.tsv
├── qc_species_colsum.tsv
├── qc_species_colsum_prop.tsv
└── sample_ids.txt
```

Abundance units:

```text
metaphlan_species_only.tsv        # percentage abundance
metaphlan_species_only_prop.tsv   # proportional abundance, 0–1
```

---

### 9. HUMAnN-Compatible MetaPhlAn Profiling

**Script:** `06c_metaphlan_for_humann_vJun23.sh`

This script generates MetaPhlAn profiles using the vJun23 database for HUMAnN compatibility.

Before running MetaPhlAn, the script counts read pairs per sample and creates a runlist of samples passing a minimum read-depth threshold.

Default cutoff:

```bash
CUTOFF_PAIRS=5000000
```

Only samples with at least 5 million read pairs are included in the HUMAnN workflow by default.

Command:

```bash
bash 06c_metaphlan_for_humann_vJun23.sh
```

Custom cutoff example:

```bash
CUTOFF_PAIRS=3000000 bash 06c_metaphlan_for_humann_vJun23.sh
```

Parallel settings example:

```bash
COUNT_JOBS=24 COUNT_PIGZ_THREADS=4 MP_JOBS=6 MP_THREADS=20 bash 06c_metaphlan_for_humann_vJun23.sh
```

Input:

```text
/data/projects/biomage/data/processed/deduplicated_sanitized/*_R1.sanitized.fastq.gz
/data/projects/biomage/data/processed/deduplicated_sanitized/*_R2.sanitized.fastq.gz
```

Read-pair QC and runlist outputs:

```text
/data/projects/biomage/data/processed/
├── qc_humann_readpairs.tsv
└── humann_runlist_5000000.txt
```

MetaPhlAn vJun23 outputs:

```text
/data/projects/biomage/data/processed/metaphlan_vJun23/
├── SAMPLE_profile.txt
├── SAMPLE.mapout.bz2
└── logs/
    └── SAMPLE.metaphlan.log
```

---

### 10. HUMAnN3 Functional Profiling

**Main script used:** `07_humann3_run_5M_vJun23.sh`

This script runs HUMAnN3 functional profiling on the samples listed in the HUMAnN runlist.

It uses:

- sanitized paired FASTQ files
- MetaPhlAn vJun23 taxonomic profiles
- an on-the-fly interleaved FASTQ generated using `reformat.sh`

Command:

```bash
bash 07_humann3_run_5M_vJun23.sh
```

Optional resource settings:

```bash
JOBS=5 THREADS=24 bash 07_humann3_run_5M_vJun23.sh
```

Optional custom runlist:

```bash
RUNLIST=/data/projects/biomage/data/processed/humann_runlist_3000000.txt bash 07_humann3_run_5M_vJun23.sh
```

Input FASTQ files:

```text
/data/projects/biomage/data/processed/deduplicated_sanitized/
├── SAMPLE_R1.sanitized.fastq.gz
└── SAMPLE_R2.sanitized.fastq.gz
```

Input MetaPhlAn vJun23 profile:

```text
/data/projects/biomage/data/processed/metaphlan_vJun23/SAMPLE_profile.txt
```

Input runlist:

```text
/data/projects/biomage/data/processed/humann_runlist_5000000.txt
```

Per-sample output:

```text
/data/projects/biomage/data/processed/humann3/SAMPLE/
├── SAMPLE_genefamilies.tsv
├── SAMPLE_pathabundance.tsv
└── SAMPLE_pathcoverage.tsv
```

Log files:

```text
/data/projects/biomage/data/processed/humann3/logs/
└── SAMPLE.humann.log
```

The repository also contains older or alternative HUMAnN scripts:

```text
07_humann3_run.sh
07_humann3_run_v3_serial.sh
```

The recommended script for the final workflow is:

```text
07_humann3_run_5M_vJun23.sh
```

---

### 11. HUMAnN3 Postprocessing

**Script:** `07b_humann3_postprocess_v6.sh`

This script performs HUMAnN postprocessing across completed sample folders.

It performs:

1. Completeness checking for each sample
2. Creation of a flat symlink folder
3. Merging per-sample HUMAnN tables
4. Renormalization to relative abundance and CPM
5. Renaming MetaCyc pathways
6. Column cleanup
7. QC file generation

Command:

```bash
bash 07b_humann3_postprocess_v6.sh
```

Input:

```text
/data/projects/biomage/data/processed/humann3/SAMPLE/
├── SAMPLE_genefamilies.tsv
├── SAMPLE_pathabundance.tsv
└── SAMPLE_pathcoverage.tsv
```

Output directory:

```text
/data/projects/biomage/data/processed/humann3/postprocess/
```

Main outputs:

```text
humann_completed_samples.txt
humann_incomplete_samples.txt
humann_genefamilies_merged.tsv
humann_pathabundance_merged.tsv
humann_pathcoverage_merged.tsv
humann_genefamilies_relab.tsv
humann_pathabundance_relab.tsv
humann_genefamilies_cpm.tsv
humann_pathabundance_cpm.tsv
humann_pathabundance_relab_named.tsv
humann_pathabundance_cpm_named.tsv
humann_pathcoverage_clean.tsv
humann_pathabundance_relab_named_stratified.tsv
humann_pathabundance_relab_named_unstratified.tsv
humann_pathabundance_cpm_named_stratified.tsv
humann_pathabundance_cpm_named_unstratified.tsv
qc_pathabundance_relab_colsum.tsv
qc_genefamilies_relab_colsum.tsv
qc_pathcoverage_stats.tsv
```

Key pathway files:

```text
humann_pathabundance_relab_named_unstratified.tsv
humann_pathabundance_relab_named_stratified.tsv
humann_pathabundance_cpm_named_unstratified.tsv
humann_pathabundance_cpm_named_stratified.tsv
```

---

### 12. HUMAnN3 Visualization

**Script:** `07c_humann3_plots.py`

This script generates HUMAnN QC and pathway visualizations from postprocessed HUMAnN outputs.

Command:

```bash
python 07c_humann3_plots.py
```

Input:

```text
/data/projects/biomage/data/processed/humann3/postprocess/
```

Output:

```text
/data/projects/biomage/data/processed/humann3/postprocess/figures/
```

Expected figure files:

```text
humann_pathway_top20_barplot.png
humann_pathway_top20_heatmap.png
humann_qc_pathabundance_colsum.png
humann_qc_genefamilies_colsum.png
humann_qc_pathcoverage_mean.png
```

---

## Optional / Unused Script

### Kraken2 and Bracken

**Script:** `05_kraken2_bracken_v2.sh`

This script is present in the repository but was **not used in the final analysis for this study**.

The final taxonomic profiling workflow used:

```text
MetaPhlAn4
```

not Kraken2/Bracken.

---

## Main Output Summary

### Cleaned FASTQ files

```text
/data/projects/biomage/data/processed/deduplicated_sanitized/
├── SAMPLE_R1.sanitized.fastq.gz
└── SAMPLE_R2.sanitized.fastq.gz
```

### Taxonomic profiles

```text
/data/projects/biomage/data/processed/metaphlan/
├── metaphlan_merged_abundance.tsv
├── metaphlan_species_only.tsv
├── metaphlan_species_only_prop.tsv
├── metaphlan_sgb_only.tsv
├── metaphlan_sgb_only_prop.tsv
├── metaphlan_genus_prop.tsv
└── metaphlan_phylum_prop.tsv
```

### Functional profiles

```text
/data/projects/biomage/data/processed/humann3/postprocess/
├── humann_genefamilies_relab.tsv
├── humann_genefamilies_cpm.tsv
├── humann_pathabundance_relab_named_unstratified.tsv
├── humann_pathabundance_relab_named_stratified.tsv
├── humann_pathabundance_cpm_named_unstratified.tsv
├── humann_pathabundance_cpm_named_stratified.tsv
└── humann_pathcoverage_clean.tsv
```

### Figures

```text
/data/projects/biomage/data/processed/humann3/postprocess/figures/
```

---

## Downstream Analysis

The resulting taxonomic and functional abundance tables can be used for downstream microbiome analyses, including:

- alpha diversity
- beta diversity
- taxonomic composition
- functional composition
- association analysis with healthy ageing metadata

R scripts are stored in:

```text
/data/projects/biomage/src/R/
```

---

## Resource Notes

Several scripts are configured for a high-resource workstation and use up to 102 threads by default. If running on a smaller machine or an HPC job with lower resource allocation, edit the thread and job settings before running.

Examples of configurable settings used by the workflow include:

```bash
THREADS=24
JOBS=5
PIGZ_THREADS=4
COUNT_JOBS=24
COUNT_PIGZ_THREADS=4
MP_JOBS=6
MP_THREADS=20
```

---

## Path Notes

The workflow uses absolute paths under:

```text
/data/projects/biomage/
```

Before running the workflow on a new system, check and update paths for:

```text
/data/projects/biomage/data/raw
/data/projects/biomage/data/processed
/data/projects/biomage/data/db/metaphlan4/vJan25
/data/projects/biomage/data/db/metaphlan4/vJun23
/data/projects/biomage/host/hg38
```

---

## Troubleshooting

### No input files found

Check that files follow the expected paired-end naming pattern:

```text
*_R1*.fastq.gz
*_R2*.fastq.gz
```

### Missing mate pair

Confirm that each sample has both R1 and R2 files.

### Corrupt gzip file

Test FASTQ files with:

```bash
gzip -t file.fastq.gz
```

### HUMAnN sample skipped

Check whether the sample passed the read-pair cutoff used to create the HUMAnN runlist.

Default cutoff:

```bash
CUTOFF_PAIRS=5000000
```

Check the runlist:

```bash
cat /data/projects/biomage/data/processed/humann_runlist_5000000.txt
```

### Incomplete HUMAnN outputs

Each completed HUMAnN sample folder should contain:

```text
SAMPLE_genefamilies.tsv
SAMPLE_pathabundance.tsv
SAMPLE_pathcoverage.tsv
```

Check logs:

```bash
less /data/projects/biomage/data/processed/humann3/logs/SAMPLE.humann.log
```

---

## Final Project Summary

This workflow processes paired-end shotgun metagenomic sequencing reads from raw FASTQ files to taxonomic and functional profiles.

The final pipeline includes:

1. Adapter trimming and quality filtering
2. Low-complexity read filtering
3. Paired-end read screening and repair
4. Human host-read removal
5. Duplicate read removal
6. FASTQ sanitization
7. MetaPhlAn4 taxonomic profiling
8. MetaPhlAn postprocessing
9. HUMAnN-compatible MetaPhlAn profiling
10. HUMAnN3 functional profiling
11. HUMAnN postprocessing
12. HUMAnN visualization

The resulting taxonomic and functional abundance tables are suitable for downstream analyses of gut microbiome composition, microbial diversity, functional potential, and associations with healthy ageing.
