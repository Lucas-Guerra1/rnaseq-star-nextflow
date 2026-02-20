# RNA-Seq Pipeline — Nextflow DSL2 (STAR)

An automated pipeline for RNA-seq data analysis, implemented in Nextflow DSL2. The pipeline automatically discovers paired-end and interleaved FASTQ files, performs quality control, aligns reads to a reference genome with STAR, quantifies gene expression with featureCounts, and aggregates all reports with MultiQC.

---

## Table of Contents

- [Overview](#overview)
- [Workflow](#workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Format](#input-format)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
- [Parameters](#parameters)
- [Output Structure](#output-structure)
- [Output Files](#output-files)
- [Technical Description](#technical-description)
- [Limitations](#limitations)
- [License](#license)

---

## Overview

This pipeline performs end-to-end RNA-seq analysis:

- Automatic detection of paired-end and interleaved FASTQ files
- Quality control of raw reads
- Genome index generation with STAR (with automatic `genomeSAindexNbases` adjustment based on genome size)
- Paired-end alignment with STAR
- BAM indexing and flagstat reporting
- Gene expression quantification with featureCounts
- Aggregated QC and mapping report with MultiQC

---

## Workflow

```
FASTQ (paired-end or interleaved)
        │
        ▼
Interleaved split (if needed)
        │
        ▼
FastQC → STAR Index → STAR Align → BAM Index → featureCounts → MultiQC
```

---

## Requirements

### Software

| Tool | Version | Purpose |
|---|---|---|
| Nextflow | ≥ 22 | Workflow manager |
| Java | ≥ 11 | Required by Nextflow |
| Conda or Mamba | any | Dependency management |

### Bioinformatics tools (managed via Conda)

- `STAR` ≥ 2.7
- `FastQC`
- `MultiQC`
- `Samtools`
- `featureCounts` (Subread package)

All tools are installed automatically via the `envs/bioinfo.yml` Conda environment.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Lucas-Guerra1/rnaseq-star-nextflow.git
cd rnaseq-star-nextflow
```

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

> Requires Java ≥ 11. Check with `java -version`.

### 3. Create the Conda environment

```bash
conda env create -f envs/bioinfo.yml
```

> If using Mamba: `mamba env create -f envs/bioinfo.yml`

### 4. Set up the directory structure

```bash
bash setup.sh
```

This creates the expected `data/` structure for your input files.

---

## Input Format

### Reference genome

Place your genome FASTA file anywhere inside `data/`. The pipeline discovers it automatically:

```
data/genome/genome.fa        # or .fasta / .fna
```

### Annotation

Place your GTF or GFF annotation file inside `data/`. The pipeline discovers it automatically:

```
data/annotation/annotation.gtf   # or .gff / .gff3
```

### FASTQ reads

The pipeline supports two input formats, detected automatically:

**Paired-end** (separate R1/R2 files):
```
data/reads/sample_1.fastq.gz
data/reads/sample_2.fastq.gz
```

**Interleaved** (R1 and R2 in the same file):
```
data/reads/sample.fastq.gz
```

> Interleaved files are split automatically before alignment. Files without `_1` / `_2` in their names are treated as interleaved.

---

## Directory Structure

```
rnaseq-star-nextflow/
├── data/
│   ├── genome/          # Reference genome FASTA (.fa, .fasta, .fna)
│   ├── annotation/      # Gene annotation (.gtf, .gff, .gff3)
│   └── reads/           # FASTQ files (paired-end or interleaved)
├── envs/
│   └── bioinfo.yml      # Conda environment definition
├── main.nf              # Main Nextflow pipeline
├── nextflow.config      # Resource and executor configuration
├── setup.sh             # Directory structure setup script
└── README.md
```

---

## Usage

### Basic run

```bash
nextflow run main.nf
```

The pipeline will automatically find the genome, annotation, and FASTQ files inside `data/`.

### Run with explicit parameters

```bash
nextflow run main.nf \
  --outdir my_results
```

### Run with multiple CPUs

Configure resources in `nextflow.config`:

```groovy
process.cpus   = 8
process.memory = '16 GB'
```

### Resume a failed run

```bash
nextflow run main.nf -resume
```

---

## Parameters

| Parameter | Description | Default |
|---|---|---|
| `--reads_paired` | Glob pattern for paired-end FASTQ files | `data/**/*_{1,2}.f*q{,.gz}` |
| `--reads_interleaved` | Glob pattern for interleaved FASTQ files | `data/**/*.f*q{,.gz}` |
| `--outdir` | Output directory | `results_rnaseq` |

> The reference genome and annotation are discovered automatically from `data/`. No manual path configuration is required.

---

## Output Structure

```
results_rnaseq/
├── 00_preprocessed/     # Split reads from interleaved files
├── 01_fastqc/           # FastQC quality reports
├── 02_star_index/       # STAR genome index
├── 03_mapping/          # Sorted BAM files, indexes, flagstat, and STAR logs
├── 04_counts/           # featureCounts output (counts.txt)
└── 05_multiqc/          # Aggregated MultiQC report
```

---

## Output Files

| File | Description |
|---|---|
| `*_fastqc.html` | Per-sample FastQC quality report |
| `star_index/` | STAR genome index directory |
| `*.sorted.bam` | Coordinate-sorted BAM file |
| `*.sorted.bam.bai` | BAM index |
| `*.Log.final.out` | STAR mapping statistics summary |
| `*.flagstat.txt` | Samtools flagstat alignment summary |
| `counts.txt` | Gene-level read counts (featureCounts) |
| `counts.txt.summary` | featureCounts assignment summary |
| `multiqc_report.html` | Aggregated QC and mapping report |

---

## Technical Description

### Automatic genome size detection

The STAR indexing step automatically adjusts `--genomeSAindexNbases` based on genome size:

| Genome size | `genomeSAindexNbases` |
|---|---|
| < 10 Mb | 10 |
| 10–100 Mb | 12 |
| > 100 Mb | 14 (STAR default) |

### Automatic RAM detection

The alignment step calculates `--limitBAMsortRAM` dynamically as 80% of available system memory (minimum 2 GB), avoiding out-of-memory failures on different machines.

### Interleaved detection

Files are classified as interleaved if they do not contain `_1` or `_2` before the extension. These are split into R1/R2 using an `awk`-based approach before alignment.

### featureCounts fallback

featureCounts is first run in paired-end mode (`-p --countReadPairs`). If this fails (e.g., for single-end libraries), the pipeline automatically retries in single-end mode.

---

## Limitations

- Designed for **Linux** only; not tested on macOS or Windows
- Reference genome and annotation must be placed inside `data/` for automatic discovery
- Only supports standard Illumina short reads; long-read data is not supported
- Stranded library protocols are not configured by default — set the `-s` flag in featureCounts manually if needed
- Does not perform differential expression analysis; `counts.txt` is the final quantification output

---

## License

MIT License — Copyright (c) 2026. Feel free to use, modify, and distribute with attribution.
