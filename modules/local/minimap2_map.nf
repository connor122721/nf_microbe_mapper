// ========================================================================================
// MINIMAP2_MAP - map HiFi reads to a single reference, output sorted+indexed BAM
//                runs once per (sample, reference) pair (see workflow combine())
// ========================================================================================

process MINIMAP2_MAP {
    tag "${meta.sample}__${meta.ref_id}"
    label 'process_high'
    publishDir "${params.outdir}/mapping/${meta.sample}", mode: params.publish_dir_mode

    // Mulled biocontainer: minimap2 2.28 + samtools 1.21
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'

    input:
        tuple val(meta), path(reads), path(reference)

    output:
        tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
        tuple val(meta), path("*.flagstat.txt"),                          emit: flagstat
        tuple val(meta), path("*.stats.txt"),                             emit: stats
        tuple val(meta), path("*.idxstats.txt"),                          emit: idxstats
        path "versions.yml",                                              emit: versions

    script:
        def prefix     = "${meta.sample}.vs.${meta.ref_id}"
        def preset     = params.minimap2_preset ?: 'map-hifi'
        def extra_args = params.minimap2_args   ?: ''
        def sort_threads = Math.max(1, (task.cpus as int).intdiv(4))
        def map_threads  = Math.max(1, (task.cpus as int) - sort_threads)
        """
        # Stream minimap2 -> samtools sort to keep memory bounded
        minimap2 \\
            -t ${map_threads} \\
            -ax ${preset} \\
            ${extra_args} \\
            ${reference} \\
            ${reads} \\
        | samtools sort \\
            -@ ${sort_threads} \\
            -O bam \\
            -o ${prefix}.sorted.bam -

        samtools index -@ ${task.cpus} ${prefix}.sorted.bam

        samtools flagstat -@ ${task.cpus} ${prefix}.sorted.bam > ${prefix}.flagstat.txt
        samtools stats    -@ ${task.cpus} ${prefix}.sorted.bam > ${prefix}.stats.txt
        samtools idxstats ${prefix}.sorted.bam              > ${prefix}.idxstats.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
            samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
        END_VERSIONS
        """
}
