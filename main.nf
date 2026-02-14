nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: fastqc_files

    script:
    "fastqc ${reads}"
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
    STAR --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles ${genome} \
         --sjdbGTFfile ${gtf} \
         --runThreadN ${task.cpus} \
         --genomeSAindexNbases 6 \
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
    path "${sample_id}.sorted.bam"
    path "${sample_id}.Log.final.out", emit: logs

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${index} \
         --readFilesIn ${reads} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${sample_id}.
    
    mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.sorted.bam
    """
}

process MULTIQC {
    publishDir "${params.outdir}/04_multiqc", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    path files

    output:
    path "multiqc_report.html"

    script:
    "multiqc ."
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    genome_file = file(params.genome, checkIfExists: true)
    gtf_file = file(params.gtf, checkIfExists: true)
    
    fastqc_out = FASTQC(reads_ch)
    index_ch = STAR_INDEX(genome_file, gtf_file)
    mapping_out = STAR_ALIGN(reads_ch, index_ch)

    // Coletar todos os arquivos de log e QC para o MultiQC
    all_reports = fastqc_out.fastqc_files.collect().combine(mapping_out.logs.collect())
    MULTIQC(all_reports)
}
