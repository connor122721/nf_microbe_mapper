// ========================================================================================
// SAMPLE_REPORT - render a per-sample HTML summary combining mapping stats and top taxa
// The python script auto-discovers files by name suffix, so we just stage all per-sample
// files into one directory and pass the directory.
// ========================================================================================

process SAMPLE_REPORT {
    tag "${meta.sample}"
    label 'process_low'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    // Plain python from Docker Hub - stdlib only, no extra deps needed
    container 'docker.io/library/python:3.11-slim'

    input:
        tuple val(meta), path('inputs/*')
        path  report_script

    output:
        tuple val(meta), path("${meta.sample}.report.html"), emit: html
        tuple val(meta), path("${meta.sample}.summary.tsv"), emit: summary
        path "versions.yml",                                  emit: versions

    script:
        """
        python ${report_script} \\
            --sample ${meta.sample} \\
            --inputs-dir inputs \\
            --html-out ${meta.sample}.report.html \\
            --tsv-out  ${meta.sample}.summary.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //')
        END_VERSIONS
        """
}
