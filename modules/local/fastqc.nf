// ========================================================================================
// FASTQC - read-level QC for HiFi reads
// ========================================================================================

process FASTQC {
    tag "${meta.sample}"
    label 'process_low'
    publishDir "${params.outdir}/qc/fastqc", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample}/*.html"), emit: html
        tuple val(meta), path("${meta.sample}/*.zip"),  emit: zip
        path "versions.yml",                            emit: versions

    script:
        // FastQC --memory expects MB-per-thread, clamped to 250..10000
        def per_thread_mb = ((task.memory.toMega() as long) / Math.max(task.cpus as int, 1)) as long
        def fq_mem = Math.max(250L, Math.min(10000L, per_thread_mb))
        """
        mkdir -p ${meta.sample}
        fastqc \\
            --threads ${task.cpus} \\
            --memory ${fq_mem} \\
            --outdir ${meta.sample} \\
            ${reads}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$(fastqc --version | sed 's/^FastQC v//')
        END_VERSIONS
        """
}
