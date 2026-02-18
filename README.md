# 🧬 RNA-seq Universal Pipeline

An automated RNA-seq pipeline built with Nextflow DSL2 for gene expression analysis.
Supports both conventional paired-end reads and interleaved FASTQ files, with automatic
input detection and dynamic parameter adjustment.

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Parameters](#parameters)
- [Output](#output)
- [Author](#author)

## 🔍 Overview

This pipeline performs a complete RNA-seq analysis — from read quality control to gene
quantification — generating consolidated reports with MultiQC. The entire process is
automated: simply organize your files inside the `data/` folder and run the pipeline.
It will automatically locate the genome, annotation, and reads.
```
Reads (paired / interleaved)
        │
        ▼
  [ FASTQC ] ──────────────────────────────────┐
        │                                       │
        ▼                                       │
  [ STAR Index ] → [ STAR Align ] → [ Index BAM ] → [ featureCounts ] → [ MultiQC ]
```

## ✨ Features

- Automatic detection of genome (`.fa`, `.fasta`, `.fna`) and annotation
  (`.gtf`, `.gff`, `.gff3`) files inside `data/`
- Dual support for **paired-end** reads (`_1`/`_2` files) and **interleaved**
  reads (single interleaved FASTQ)
- Automatic splitting of interleaved reads via dedicated `SPLIT_INTERLEAVED` process
- Dynamic calculation of `genomeSAindexNbases` based on genome size
- Intelligent RAM allocation for STAR Aligner (80% of available system memory)
- Consolidated MultiQC report aggregating QC, mapping, and count metrics
- Conda-based environment management

## 📦 Requirements

| Tool | Recommended Version | Description |
|------|-------------------|-------------|
| Nextflow | ≥ 22.10 | Pipeline orchestrator |
| Conda | ≥ 4.12 | Environment manager |
| STAR | ≥ 2.7 | RNA-seq aligner |
| Samtools | ≥ 1.15 | BAM file manipulation |
| FastQC | ≥ 0.11 | Quality control |
| featureCounts | ≥ 2.0 (Subread) | Gene quantification |
| MultiQC | ≥ 1.14 | Aggregated report |

> All dependencies should be listed in `envs/bioinfo.yml`.

## 🗂 Directory Structure
```
project/
├── data/
│   ├── genome/
│   │   └── genome.fasta             # Reference genome
│   ├── annotation/
│   │   └── annotation.gtf           # Genome annotation
│   └── reads/
│       ├── sample1_1.fastq.gz       # Paired-end R1
│       ├── sample1_2.fastq.gz       # Paired-end R2
│       └── sample2.fastq.gz         # Interleaved
├── envs/
│   └── bioinfo.yml                  # Conda environment
├── setup.sh                         # Directory setup script
├── main.nf                          # Main pipeline script
└── results_rnaseq/                  # Auto-generated output
```

## 🚀 Usage

### 1. Clone the repository
```bash
git clone https://github.com/Lucas-Guerra1/rnaseq-star-nextflow.git
cd rnaseq-star-nextflow
```

### 2. Set up the directory structure

A setup script is provided to create the required folder structure automatically:
```bash
bash setup.sh
```

This will create the following directories:
```
data/
├── genome/       # Place your reference genome here (.fa, .fasta, .fna)
├── annotation/   # Place your annotation file here (.gtf, .gff, .gff3)
└── reads/        # Place your FASTQ read files here
```

> After running the script, add your files to the appropriate folders before
> executing the pipeline.

### 3. Organize your data
```bash
# Place genome, annotation, and reads inside their respective data/ subdirectories
# Subdirectory structures are supported (recursive search)
```

### 4. Run the pipeline
```bash
nextflow run main.nf
```

### 5. (Optional) Customize the output directory
```bash
nextflow run main.nf --outdir my_results
```

## ⚙️ Pipeline Steps

### 1. SPLIT_INTERLEAVED
Detects and splits interleaved FASTQ files into two independent files
(`_1.fastq.gz` and `_2.fastq.gz`), validating read count parity before processing.

### 2. FASTQC
Evaluates read quality for all samples, generating `.html` and `.zip` reports.

### 3. STAR_INDEX
Builds the reference genome index. The `--genomeSAindexNbases` parameter is
calculated automatically based on genome size:

| Genome Size | genomeSAindexNbases |
|-------------|-------------------|
| < 10 Mb | 10 |
| 10 – 100 Mb | 12 |
| > 100 Mb | 14 |

### 4. STAR_ALIGN
Aligns reads to the reference genome using splice-aware parameters optimized for
RNA-seq. RAM limit is calculated dynamically (80% of total system memory).

### 5. INDEX_BAM
Indexes sorted BAM files with `samtools index` for efficient random access.

### 6. FEATURE_COUNTS
Quantifies reads per gene using featureCounts. Attempts paired-end mode first
(`-p -B -C`); falls back to single-end mode if needed.

### 7. MULTIQC
Aggregates all QC reports, mapping logs, flagstats, and count summaries into a
single interactive HTML report.

## 🛠️ Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results_rnaseq` | Output directory |
| `--reads_paired` | `data/**/*_{1,2}.f*q{,.gz}` | Glob pattern for paired-end reads |
| `--reads_interleaved` | `data/**/*.f*q{,.gz}` | Glob pattern for interleaved reads |

## 📊 Output
```
results_rnaseq/
├── 00_preprocessed/   # Split interleaved reads
├── 01_fastqc/         # FastQC quality reports
├── 02_star_index/     # STAR genome index
├── 03_mapping/        # Aligned BAMs, indexes, flagstat and STAR logs
├── 04_counts/         # counts.txt and counts.txt.summary (featureCounts)
└── 05_multiqc/        # Consolidated MultiQC report
```

## 👤 Author

**Lucas Guerra**
Federal University of Lavras (UFLA) — Brazil
Ph.D. candidate | M.Sc. Plant Biotechnology
[GitHub](https://github.com/Lucas-Guerra1) ·
[LinkedIn](https://www.linkedin.com/in/lucas-ribeiro-de-souza-guerra-082621186/)
