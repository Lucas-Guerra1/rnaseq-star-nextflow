# SNP Phylogenomics Pipeline (Nextflow DSL2)

An automated pipeline for SNP detection, genome alignment construction, and phylogenetic inference from paired-end Illumina sequencing data, implemented in Nextflow DSL2.

---

## Table of Contents

- [Overview](#overview)
- [Workflow](#workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Format](#input-format)
- [Usage](#usage)
- [Parameters](#parameters)
- [Output Structure](#output-structure)
- [Output Files](#output-files)
- [Technical Description](#technical-description)
- [Use Cases](#use-cases)
- [Limitations](#limitations)
- [Suggested Improvements](#suggested-improvements)
- [License](#license)

---

## Overview

This pipeline performs end-to-end SNP-based phylogenomics:

- Read quality control and trimming
- Reference genome alignment
- Per-sample variant calling
- SNP filtering
- Multi-sample VCF merging
- FASTA alignment generation from SNPs
- Phylogenetic tree inference

The pipeline is fully reproducible, modular, and scalable across local machines and HPC environments.

---

## Workflow

```
Reads → QC (fastp) → Alignment (BWA) → Variant Calling (bcftools) → SNP Filter → Merge VCFs → SNP Alignment → Phylogeny (FastTree)
```

**Detailed steps:**

1. **Reference indexing** — `bwa index`
2. **Quality control and trimming** — `fastp` with automatic adapter detection
3. **Genome alignment** — `bwa mem`, sorting and indexing with `samtools`
4. **Per-sample variant calling** — `bcftools mpileup` + `bcftools call`
5. **SNP filtering** — retains only biallelic SNPs
6. **Multi-sample VCF merging** — `bcftools merge`
7. **SNP-to-FASTA alignment** — concatenated SNP matrix with corrected sample names
8. **Phylogenetic inference** — `FastTree` with GTR model

---

## Requirements

### Software

| Tool | Version |
|---|---|
| Nextflow | ≥ 22 |
| Conda or Mamba | any |

### Bioinformatics tools (managed via Conda)

- `bwa`
- `fastp`
- `samtools`
- `bcftools`
- `FastTree`

All tools are installed automatically via the `envs/bioinfo.yml` Conda environment.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/snp-phylogenomics-nextflow.git
cd snp-phylogenomics-nextflow
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

---

## Input Format

### Paired-end FASTQ reads

The pipeline expects paired-end FASTQ files following these naming conventions:

```
data/sample_1.fastq.gz    data/sample_2.fastq.gz
data/SAMPLEID_R1.fastq.gz data/SAMPLEID_R2.fastq.gz
```

Configurable via `--reads` (default: `data/*_{1,2}.fastq.gz`).

### Reference genome

A single FASTA file:

```
data/ref.fa
```

Configurable via `--ref_fasta`.

---

## Usage

### Basic run

```bash
nextflow run main.nf
```

### Run with explicit parameters

```bash
nextflow run main.nf \
  --reads "data/*_{1,2}.fastq.gz" \
  --ref_fasta data/ref.fa \
  --outdir results
```

### Run with multiple CPUs

```bash
nextflow run main.nf -profile local --max_cpus 8
```

> Configure CPU and memory limits in `nextflow.config`.

### Resume a failed run

```bash
nextflow run main.nf -resume
```

---

## Parameters

| Parameter | Description | Default |
|---|---|---|
| `--reads` | Glob pattern for paired-end FASTQ files | `data/*_{1,2}.fastq.gz` |
| `--ref_fasta` | Reference genome FASTA file | `data/ref.fa` |
| `--outdir` | Output directory | `results` |

The pipeline automatically validates that `--reads` and `--ref_fasta` are provided.

---

## Output Structure

```
results/
├── 00_ref_index/       # BWA reference indices
├── 01_fastp/           # Trimmed reads
├── 02_bam/             # Sorted and indexed BAM files
├── 03_vcf/             # Per-sample VCF files
├── 04_snps/            # Biallelic SNP-only VCF files
├── 05_merged_vcf/      # Multi-sample merged VCF
├── 06_alignment/       # SNP FASTA alignment
└── 07_tree/            # Phylogenetic tree (Newick)
```

---

## Output Files

| File | Description |
|---|---|
| `*.sorted.bam` | Sorted alignment file |
| `*.sorted.bam.bai` | BAM index |
| `*.vcf.gz` | Per-sample variant calls |
| `*.snps.vcf.gz` | Biallelic SNPs only |
| `all_samples.vcf.gz` | Multi-sample merged VCF |
| `snp_alignment.fasta` | Concatenated SNP matrix (FASTA) |
| `tree.nwk` | Phylogenetic tree in Newick format |

---

## Technical Description

### SNP alignment strategy

The pipeline:

1. Extracts only SNP positions from the merged VCF
2. Removes phasing symbols (`/` and `|`)
3. Uses the first allele per position
4. Replaces missing data with `N`
5. Concatenates all SNPs per sample
6. Generates a multi-sequence FASTA alignment with correct sample identifiers

### Pipeline architecture

Implemented in Nextflow DSL2 with:

- Independent modular processes
- Typed channels
- Automatic parallelization
- Reproducibility via versioned Conda environments
- Compatibility with local and HPC execution

---

## Use Cases

- SNP-based bacterial phylogeny
- Comparative genomics
- Epidemiological surveillance
- Population structure studies
- Evolutionary analysis

---

## Limitations

- Assumes Illumina paired-end sequencing data only
- Retains only biallelic SNPs; multiallelic variants are excluded
- No advanced variant quality filtering (DP, MQ, QUAL thresholds)
- Does not mask repetitive regions
- Does not include recombination filtering

---

## Suggested Improvements

- SNP quality filtering by depth (DP), mapping quality (MQ), and QUAL score
- Repetitive region masking (e.g., with RepeatMasker)
- Recombination detection (e.g., ClonalFrameML, Gubbins)
- Long-read data support (PacBio, Oxford Nanopore)
- Additional phylogenetic models and bootstrap support

---

## License

MIT License — Copyright (c) 2026. Feel free to use, modify, and distribute with attribution.
