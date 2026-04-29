// ========================================================================================
// MULTIQC - aggregate FastQC, samtools stats/flagstat, kraken2 reports
// ========================================================================================

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/qc/multiqc", mode: params.publish_dir_mode

    // container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'
    container 'docker://staphb/multiqc:latest'

    input:
        path '*'                            // anything multiqc can parse
        path multiqc_config, stageAs: 'multiqc_config.yaml'

    output:
        path "multiqc_report.html",         emit: report
        // path "multiqc_data",                emit: data
        // path "versions.yml",                emit: versions

    script:
        def cfg_arg = multiqc_config ? "--config ${multiqc_config}" : ''
        """
        multiqc \\
            --force \\
            --filename multiqc_report.html \\
            ${cfg_arg} \\
            .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$(multiqc --version | sed 's/^multiqc, version //')
        END_VERSIONS
        """
}
