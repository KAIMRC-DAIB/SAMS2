# SAMS2 Metagenomics

Shotgun metagenomics workflow developed for the **SAMS2 Healthy Ageing Gut Microbiome Project**, investigating gut microbiome composition and function in individuals aged 55 years and above.

The pipeline processes paired-end shotgun metagenomic sequencing data from raw FASTQ files through quality control, host-read removal, taxonomic profiling, and functional profiling.

## Workflow Overview

The main workflow consists of:

1. **Quality control and adapter trimming** — `fastp`
2. **Low-complexity filtering** — `BBDuk`
3. **Paired-end read validation and repair** — `BBMap`
4. **Human host-read removal** — `Bowtie2`
5. **Read deduplication** — `Clumpify`
6. **FASTQ sanitization and pair repair**
7. **Taxonomic profiling** — `MetaPhlAn 4`
8. **Taxonomic table postprocessing**
9. **Functional profiling** — `HUMAnN 3`
10. **postprocessing and visualization**

## Repository Structure

```text
.
├── 01_fastp_qc_v2.sh
├── 02_bbduk_low_complexity.sh
├── 02b_screen_repair_pairs.sh
├── 03_bowtie2_host_removal_v3.sh
├── 04_clumpify_dedup_v3.sh
├── 04b_sanitize_fastq.sh
├── 06_metaphlan_v2.sh
├── 06b_metaphlan_postprocess.sh
├── 06c_metaphlan_for_humann_vJun23.sh
├── 07_humann3_run_5M_vJun23.sh
├── 07b_humann3_postprocess_v6.sh
├── 07c_humann3_plots.py
├── rename_fastq.sh
└── R/ Combined_Analysis_Script_FINAL.R
```

## Requirements

The workflow uses the following main tools:

- `fastp`
- `BBMap`
- `Bowtie2`
- `SAMtools`
- `BEDTools`
- `MetaPhlAn 4`
- `HUMAnN 3`
- Python 3
- R

Python dependencies used for visualization include:

```text
pandas
numpy
matplotlib
```

Conda environments are recommended for managing dependencies.

## Input

The workflow expects paired-end shotgun metagenomic FASTQ files using an R1/R2 naming convention, for example:

```text
SAMPLE_R1.fastq.gz
SAMPLE_R2.fastq.gz
```

## Running the Workflow

Scripts are designed to be executed sequentially:

```bash
# Optional FASTQ renaming
bash rename_fastq.sh

# Read preprocessing
bash 01_fastp_qc_v2.sh /path/to/raw_fastq
bash 02_bbduk_low_complexity.sh
bash 02b_screen_repair_pairs.sh

# Host removal and read cleanup
bash 03_bowtie2_host_removal_v3.sh
bash 04_clumpify_dedup_v3.sh
bash 04b_sanitize_fastq.sh

# Taxonomic profiling
bash 06_metaphlan_v2.sh
bash 06b_metaphlan_postprocess.sh

# Functional profiling
bash 06c_metaphlan_for_humann_vJun23.sh
bash 07_humann3_run_5M_vJun23.sh
bash 07b_humann3_postprocess_v6.sh

# Visualization
R Combined_Analysis_Script_FINAL.R
```

Paths, database locations, and computational resource settings should be configured for the target environment before running the workflow.

## Outputs

The pipeline produces:

- Quality-controlled and host-filtered FASTQ files
- MetaPhlAn taxonomic abundance profiles
- Species-, genus-, and higher-level abundance tables
- HUMAnN gene-family profiles
- HUMAnN pathway abundance and coverage tables
- Quality-control summaries and visualizations

These outputs can be used for downstream microbiome analyses such as:

- Alpha and beta diversity
- Taxonomic composition
- Functional composition
- Differential abundance analysis
- Association analyses with phenotype or clinical metadata

## Notes

MetaPhlAn 4 is used for the final taxonomic profiling workflow.

HUMAnN-compatible MetaPhlAn profiles are generated separately where required to maintain database compatibility between MetaPhlAn and HUMAnN.

The scripts were originally developed for a high-performance computational environment. Users should adjust thread counts, paths, and resource settings according to their system.

## Project

**SAMS2 Healthy Ageing Shotgun Metagenomics Project**

This repository provides the computational workflow used to process shotgun metagenomic sequencing data and generate taxonomic and functional microbiome profiles for downstream healthy-ageing analyses.
