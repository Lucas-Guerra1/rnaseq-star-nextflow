# rnaseq-star-nextflow# 🧬 RNA-seq Universal Pipeline

<p align="center">
  <img src="https://img.shields.io/badge/Nextflow-DSL2-brightgreen?style=for-the-badge&logo=nextflow" />
  <img src="https://img.shields.io/badge/Conda-Environment-blue?style=for-the-badge&logo=anaconda" />
  <img src="https://img.shields.io/badge/STAR-Aligner-orange?style=for-the-badge" />
  <img src="https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge" />
</p>

<p align="center">
  Pipeline automatizado de RNA-seq desenvolvido em <strong>Nextflow DSL2</strong> para análise de expressão gênica. Suporta leituras <em>paired-end</em> convencionais e arquivos <em>interleaved</em>, com detecção automática de arquivos de entrada e ajuste dinâmico de parâmetros.
</p>

---

## 📋 Índice

- [Visão Geral](#-visão-geral)
- [Funcionalidades](#-funcionalidades)
- [Requisitos](#-requisitos)
- [Estrutura de Diretórios](#-estrutura-de-diretórios)
- [Como Usar](#-como-usar)
- [Etapas do Pipeline](#-etapas-do-pipeline)
- [Parâmetros](#️-parâmetros)
- [Resultados](#-resultados)
- [Autor](#-autor)

---

## 🔍 Visão Geral

Este pipeline realiza análise completa de RNA-seq, desde o controle de qualidade das leituras até a quantificação de genes, gerando relatórios consolidados com MultiQC. Todo o processo é **automatizado**: basta organizar os arquivos na pasta `data/` e executar o pipeline — ele encontra o genoma, a anotação e as leituras automaticamente.

```
Reads (paired / interleaved)
        │
        ▼
  [ FASTQC ] ──────────────────────────────────┐
        │                                       │
        ▼                                       │
  [ STAR Index ] → [ STAR Align ] → [ Index BAM ] → [ featureCounts ] → [ MultiQC ]
```

---

## ✨ Funcionalidades

- **Detecção automática** de arquivos de genoma (`.fa`, `.fasta`, `.fna`) e anotação (`.gtf`, `.gff`, `.gff3`) dentro de `data/`
- **Suporte duplo** a reads *paired-end* (arquivos `_1`/`_2`) e *interleaved* (arquivo único intercalado)
- **Separação automática** de leituras interleaved via processo dedicado (`SPLIT_INTERLEAVED`)
- **Ajuste dinâmico** do parâmetro `genomeSAindexNbases` conforme o tamanho do genoma
- **Alocação inteligente** de memória RAM para o STAR Aligner (80% da RAM disponível)
- **Relatório consolidado** com MultiQC agregando QC, mapeamento e contagem
- **Controle de ambiente** via Conda

---

## 📦 Requisitos

| Ferramenta      | Versão Recomendada | Descrição                        |
|-----------------|--------------------|----------------------------------|
| [Nextflow](https://www.nextflow.io/) | ≥ 22.10 | Orquestrador do pipeline |
| [Conda](https://docs.conda.io/)     | ≥ 4.12  | Gerenciador de ambientes         |
| STAR            | ≥ 2.7              | Alinhador de RNA-seq             |
| Samtools        | ≥ 1.15             | Manipulação de arquivos BAM      |
| FastQC          | ≥ 0.11             | Controle de qualidade            |
| featureCounts   | ≥ 2.0 (Subread)    | Quantificação de genes           |
| MultiQC         | ≥ 1.14             | Relatório agregado               |

> O arquivo `envs/bioinfo.yml` deve conter todas as dependências acima.

---

## 🗂 Estrutura de Diretórios

```
project/
├── data/
│   ├── genome.fasta          # Genoma de referência
│   ├── annotation.gtf        # Anotação genômica
│   └── reads/
│       ├── sample1_1.fastq.gz   # Paired-end R1
│       ├── sample1_2.fastq.gz   # Paired-end R2
│       └── sample2.fastq.gz     # Interleaved
├── envs/
│   └── bioinfo.yml           # Ambiente Conda
├── main.nf                   # Pipeline principal
└── results_rnaseq/           # Saída gerada automaticamente
```

---

## 🚀 Como Usar

**1. Clone o repositório**
```bash
git clone https://github.com/Lucas-Guerra1/rnaseq-pipeline.git
cd rnaseq-pipeline
```

**2. Organize seus dados**
```bash
# Coloque o genoma, anotação e reads dentro de data/
# A estrutura de subpastas é suportada (busca recursiva)
```

**3. Execute o pipeline**
```bash
nextflow run main.nf
```

**4. (Opcional) Personalize o diretório de saída**
```bash
nextflow run main.nf --outdir meus_resultados
```

---

## ⚙️ Etapas do Pipeline

### 1. `SPLIT_INTERLEAVED`
Detecta e separa arquivos FASTQ no formato interleaved em dois arquivos independentes (`_1.fastq.gz` e `_2.fastq.gz`), validando a paridade do número de leituras antes do processamento.

### 2. `FASTQC`
Avalia a qualidade das leituras de todas as amostras (paired e interleaved separados), gerando relatórios `.html` e `.zip`.

### 3. `STAR_INDEX`
Constrói o índice do genoma de referência. O parâmetro `--genomeSAindexNbases` é calculado automaticamente com base no tamanho do genoma:

| Tamanho do Genoma | `genomeSAindexNbases` |
|-------------------|-----------------------|
| < 10 Mb           | 10                    |
| 10 – 100 Mb       | 12                    |
| > 100 Mb          | 14                    |

### 4. `STAR_ALIGN`
Alinha as leituras ao genoma de referência com parâmetros otimizados para RNA-seq (splice-aware). O limite de RAM é calculado dinamicamente (80% da memória total do sistema).

### 5. `INDEX_BAM`
Indexa os arquivos BAM ordenados com `samtools index`, possibilitando acesso aleatório eficiente.

### 6. `FEATURE_COUNTS`
Quantifica leituras por gene usando `featureCounts`. Tenta primeiro o modo *paired-end* com validação de pares (`-p -B -C`); em caso de falha, executa em modo simples como fallback.

### 7. `MULTIQC`
Agrega todos os relatórios de QC, logs de mapeamento, flagstats e sumários de contagem em um único relatório HTML interativo.

---

## 🛠️ Parâmetros

| Parâmetro            | Padrão                          | Descrição                                      |
|----------------------|---------------------------------|------------------------------------------------|
| `--outdir`           | `results_rnaseq`                | Diretório de saída dos resultados              |
| `--reads_paired`     | `data/**/*_{1,2}.f*q{,.gz}`    | Padrão glob para reads paired-end              |
| `--reads_interleaved`| `data/**/*.f*q{,.gz}`          | Padrão glob para reads interleaved             |

---

## 📊 Resultados

```
results_rnaseq/
├── 00_preprocessed/     # Reads interleaved separados
├── 01_fastqc/           # Relatórios de qualidade (FastQC)
├── 02_star_index/       # Índice do genoma (STAR)
├── 03_mapping/          # BAMs alinhados, indexados, flagstat e logs do STAR
├── 04_counts/           # counts.txt e counts.txt.summary (featureCounts)
└── 05_multiqc/          # Relatório consolidado (MultiQC)
```

---

## 👤 Autor

**Lucas Guerra**
Universidade Federal de Lavras (UFLA)

[![GitHub](https://img.shields.io/badge/GitHub-Lucas--Guerra1-181717?style=flat&logo=github)](https://github.com/Lucas-Guerra1)

---

<p align="center">
  Desenvolvido com 🧬 para análises de transcriptômica reproduzíveis e automatizadas.
</p>
