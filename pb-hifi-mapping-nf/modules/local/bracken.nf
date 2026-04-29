// ========================================================================================
// BRACKEN - re-estimate species/genus abundance from a Kraken2 report
//           Requires that the kraken2 DB was built with bracken-database (most prebuilt
//           DBs from the Langmead lab include the .kmer_distrib files).
// ========================================================================================

process BRACKEN {
    tag "${meta.sample}"
    label 'process_low'
    publishDir "${params.outdir}/taxonomy/bracken", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/bracken:3.0--h9948957_0'

    input:
        tuple val(meta), path(kraken_report)
        path  kraken2_db

    output:
        tuple val(meta), path("*.bracken.tsv"),     emit: bracken
        tuple val(meta), path("*.bracken_report.txt"), emit: report, optional: true
        path "versions.yml",                          emit: versions

    when:
        params.run_bracken

    script:
        def prefix    = "${meta.sample}"
        def read_len  = params.bracken_read_length ?: 10000  // HiFi typical
        def level     = params.bracken_level       ?: 'S'
        def threshold = params.bracken_threshold   ?: 10
        """
        # bracken refuses to overwrite, and the binary echoes plenty of info
        bracken \\
            -d ${kraken2_db} \\
            -i ${kraken_report} \\
            -o ${prefix}.bracken.tsv \\
            -w ${prefix}.bracken_report.txt \\
            -r ${read_len} \\
            -l ${level} \\
            -t ${threshold}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bracken: \$(bracken -v 2>&1 | head -n1 | sed 's/Bracken v//')
        END_VERSIONS
        """
}
