nextflow.enable.dsl=2

/*
 * CONFIGURAÇÃO AUTOMÁTICA (Lucas-Guerra1)
 * Pipeline RNA-seq universal - busca arquivos recursivamente
 */

// 1. Busca automática de Referência (recursiva em data/)
def genome_file = file("data/**/*.{fa,fasta,fna}").find { it.exists() }

// 2. Busca automática de Anotação (recursiva em data/)
def gtf_file = file("data/**/*.{gtf,gff,gff3}").find { it.exists() }

// 3. Padrões de busca para reads (recursivo em data/)
params.reads_paired = "data/**/*_{1,2}.f*q{,.gz}"
params.reads_interleaved = "data/**/*.f*q{,.gz}"
params.outdir = "results_rnaseq"

process SPLIT_INTERLEAVED {
    tag "$sample_id"
    publishDir "${params.outdir}/00_preprocessed", mode: 'copy'
    conda "envs/bioinfo.yml"
    
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz"), emit: reads
    
    script:
    """
    echo "Separando reads paired-end de: ${sample_id}"
    
    # Conta total de reads
    if [[ ${fastq} == *.gz ]]; then
        TOTAL_LINES=\$(zcat ${fastq} | wc -l)
    else
        TOTAL_LINES=\$(wc -l < ${fastq})
    fi
    
    TOTAL_READS=\$((TOTAL_LINES / 4))
    echo "Total de reads: \$TOTAL_READS"
    
    # Se número ímpar, não é paired
    if [ \$((TOTAL_READS % 2)) -ne 0 ]; then
        echo "✗ ERRO: Número ímpar de reads (\$TOTAL_READS) - não é paired-end"
        exit 1
    fi
    
    # Separar reads: linhas ímpares = R1, linhas pares = R2
    if [[ ${fastq} == *.gz ]]; then
        zcat ${fastq} | paste - - - - | awk 'NR%2==1' | tr '\\t' '\\n' | gzip > ${sample_id}_1.fastq.gz
        zcat ${fastq} | paste - - - - | awk 'NR%2==0' | tr '\\t' '\\n' | gzip > ${sample_id}_2.fastq.gz
    else
        cat ${fastq} | paste - - - - | awk 'NR%2==1' | tr '\\t' '\\n' | gzip > ${sample_id}_1.fastq.gz
        cat ${fastq} | paste - - - - | awk 'NR%2==0' | tr '\\t' '\\n' | gzip > ${sample_id}_2.fastq.gz
    fi
    
    echo "✓ Separação concluída: ${sample_id}_1.fastq.gz e ${sample_id}_2.fastq.gz"
    """
}

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    conda "envs/bioinfo.yml"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.{html,zip}", emit: fastqc_files
    
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

process STAR_INDEX {
    tag "index"
    publishDir "${params.outdir}/02_star_index", mode: 'copy'
    conda "envs/bioinfo.yml"
    
    input:
    path genome
    path gtf
    
    output:
    path "star_index", emit: index
    
    script:
    """
    mkdir star_index
    
    # Calcular genomeSAindexNbases baseado no tamanho do genoma
    GENOME_SIZE=\$(grep -v ">" ${genome} | tr -d '\\n' | wc -c)
    
    # Fórmula simplificada para diferentes tamanhos de genoma
    if [ \$GENOME_SIZE -lt 10000000 ]; then
        SA_INDEX=10
    elif [ \$GENOME_SIZE -lt 100000000 ]; then
        SA_INDEX=12
    else
        SA_INDEX=14
    fi
    
    echo "Genoma: \$GENOME_SIZE bp, usando --genomeSAindexNbases \$SA_INDEX"
    
    STAR --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles ${genome} \
         --sjdbGTFfile ${gtf} \
         --runThreadN ${task.cpus} \
         --genomeSAindexNbases \$SA_INDEX \
         --genomeChrBinNbits 10
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping", mode: 'copy'
    conda "envs/bioinfo.yml"
    
    input:
    tuple val(sample_id), path(reads)
    path index
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    path "${sample_id}.Log.final.out", emit: logs
    path "${sample_id}.flagstat.txt", emit: flagstat
    
    script:
    """
    # Calcular RAM disponível (80% da memória total, mínimo 2GB)
    TOTAL_RAM=\$(free -b | awk 'NR==2{print \$2}')
    RAM_80=\$(echo "\$TOTAL_RAM * 0.8 / 1" | bc)
    RAM_2GB=2000000000
    
    if [ \$RAM_80 -gt \$RAM_2GB ]; then
        RAM_LIMIT=\$RAM_80
    else
        RAM_LIMIT=\$RAM_2GB
    fi
    
    echo "RAM total: \$TOTAL_RAM bytes"
    echo "Usando --limitBAMsortRAM \$RAM_LIMIT bytes"
    
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${index} \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS NM MD \
         --outFilterType BySJout \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --limitBAMsortRAM \$RAM_LIMIT \
         --outFileNamePrefix ${sample_id}.
    
    mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}

process INDEX_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping", mode: 'copy'
    conda "envs/bioinfo.yml"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_bam
    
    script:
    """
    samtools index ${bam}
    """
}

process FEATURE_COUNTS {
    tag "quant"
    publishDir "${params.outdir}/04_counts", mode: 'copy'
    conda "envs/bioinfo.yml"
    
    input:
    tuple val(sample_ids), path(bams), path(bais)
    path gtf
    
    output:
    path "counts.txt", emit: counts
    path "counts.txt.summary", emit: summary
    
    script:
    def bam_list = bams.join(' ')
    """
    featureCounts -T ${task.cpus} \
                  -p --countReadPairs \
                  -B -C \
                  -a ${gtf} \
                  -o counts.txt \
                  ${bam_list} 2>&1 | tee featureCounts.log || \
    featureCounts -T ${task.cpus} \
                  -a ${gtf} \
                  -o counts.txt \
                  ${bam_list}
    """
}

process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/05_multiqc", mode: 'copy'
    conda "envs/bioinfo.yml"
    
    input:
    path all_files
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc .
    """
}

workflow {
    if (!genome_file) { 
        error "❌ Genoma (.fa/.fasta/.fna) não encontrado em 'data/'" 
    }
    if (!gtf_file) { 
        error "❌ Anotação (.gtf/.gff) não encontrada em 'data/'" 
    }
    
    log.info """
    ======================================================
    🧬 RNA-SEQ AUTOMÁTICO: Lucas-Guerra1 (UFLA)
    ======================================================
    Genome: ${genome_file}
    GTF   : ${gtf_file}
    ======================================================
    """
    
    // Detectar arquivos paired-end (com _1/_2 no nome)
    paired_ch = Channel.fromFilePairs(params.reads_paired, size: 2, checkIfExists: false)
        .filter { id, files -> files.size() == 2 }
    
    // Detectar arquivos interleaved (sem _1/_2)
    interleaved_ch = Channel.fromPath(params.reads_interleaved, checkIfExists: false)
        .filter { file -> 
            !file.name.matches('.*_[12]\\.f.*q.*') && 
            (file.name.endsWith('.fastq') || 
             file.name.endsWith('.fq') || 
             file.name.endsWith('.fastq.gz') || 
             file.name.endsWith('.fq.gz'))
        }
        .map { file -> 
            def sample_id = file.baseName
                .replaceAll(/\.fastq$/, '')
                .replaceAll(/\.fq$/, '')
                .replaceAll(/\.gz$/, '')
            [sample_id, file]
        }
    
    // Contar e reportar
    paired_ch.count().subscribe { n -> 
        if (n > 0) log.info "✓ ${n} amostra(s) PAIRED encontradas"
    }
    
    interleaved_ch.count().subscribe { n -> 
        if (n > 0) log.info "✓ ${n} amostra(s) INTERLEAVED encontradas"
    }
    
    // Recriar canais
    paired_ch = Channel.fromFilePairs(params.reads_paired, size: 2, checkIfExists: false)
        .filter { id, files -> files.size() == 2 }
    
    interleaved_ch = Channel.fromPath(params.reads_interleaved, checkIfExists: false)
        .filter { file -> 
            !file.name.matches('.*_[12]\\.f.*q.*') && 
            (file.name.endsWith('.fastq') || 
             file.name.endsWith('.fq') || 
             file.name.endsWith('.fastq.gz') || 
             file.name.endsWith('.fq.gz'))
        }
        .map { file -> 
            def sample_id = file.baseName
                .replaceAll(/\.fastq$/, '')
                .replaceAll(/\.fq$/, '')
                .replaceAll(/\.gz$/, '')
            [sample_id, file]
        }
    
    // Separar interleaved
    split_reads = SPLIT_INTERLEAVED(interleaved_ch)
    
    // Unir todos os reads
    all_reads = paired_ch
        .map { id, files -> [id, [files[0], files[1]]] }
        .mix(split_reads.reads.map { id, f1, f2 -> [id, [f1, f2]] })
    
    // Pipeline principal
    fastqc_out = FASTQC(all_reads)
    index_ch = STAR_INDEX(genome_file, gtf_file)
    mapping_out = STAR_ALIGN(all_reads, index_ch)
    indexed_bams = INDEX_BAM(mapping_out.bam)
    
    bams_collected = indexed_bams.indexed_bam
        .map { sample_id, bam, bai -> [bam, bai] }
        .collect()
        .map { files -> 
            def bams = files.findAll { it.name.endsWith('.bam') }
            def bais = files.findAll { it.name.endsWith('.bai') }
            def ids = bams.collect { it.baseName.replaceAll(/\.sorted$/, '') }
            [ids, bams, bais]
        }
    
    counts_out = FEATURE_COUNTS(bams_collected, gtf_file)
    
    report_files = fastqc_out.fastqc_files.collect()
        .mix(mapping_out.logs.collect())
        .mix(mapping_out.flagstat.collect())
        .mix(counts_out.counts)
        .mix(counts_out.summary)
        .collect()
    
    MULTIQC(report_files)
}
