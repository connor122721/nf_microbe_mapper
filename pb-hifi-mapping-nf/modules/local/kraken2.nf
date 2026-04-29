// ========================================================================================
// KRAKEN2 - taxonomic classification of HiFi reads
// ========================================================================================

process KRAKEN2 {
    tag "${meta.sample}"
    label 'process_high_memory'
    publishDir "${params.outdir}/taxonomy/kraken2", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_4'

    input:
        tuple val(meta), path(reads)
        path  kraken2_db

    output:
        tuple val(meta), path("*.kraken2.report.txt"), emit: report
        tuple val(meta), path("*.kraken2.out.txt"),    emit: classified
        path "versions.yml",                           emit: versions

    script:
        def prefix = "${meta.sample}"
        // Decompress on read if needed (HiFi fastq is often .gz)
        def gz_arg = reads.toString().endsWith('.gz') ? '--gzip-compressed' : ''
        """
        kraken2 \\
            --db ${kraken2_db} \\
            --threads ${task.cpus} \\
            --report ${prefix}.kraken2.report.txt \\
            --output ${prefix}.kraken2.out.txt \\
            --use-names \\
            --report-minimizer-data \\
            ${gz_arg} \\
            ${reads}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kraken2: \$(kraken2 --version | head -n1 | sed 's/Kraken version //')
        END_VERSIONS
        """
}
