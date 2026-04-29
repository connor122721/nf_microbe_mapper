#!/usr/bin/env nextflow

/*
========================================================================================
    pb-hifi-mapping-nf
========================================================================================
    A simple, robust Nextflow pipeline for PacBio HiFi metagenomic-style work:

      samplesheet.csv  +  one or more reference FASTAs
              │
              ├─► FastQC               (read QC)
              ├─► minimap2  (sample × reference)  ─► sorted/indexed BAM + samtools stats
              ├─► Kraken2  + Bracken   (taxonomic profile)
              ├─► Per-sample HTML report
              └─► MultiQC              (aggregate report)

    Author : Connor S. Murray, PhD
    License: MIT
========================================================================================
*/

nextflow.enable.dsl = 2

// ========================================================================================
// HELP
// ========================================================================================

def helpMessage() {
    log.info"""
    ========================================================================================
     pb-hifi-mapping-nf  v${workflow.manifest.version ?: 'dev'}
    ========================================================================================

    Usage:
        nextflow run main.nf \\
            --samplesheet samplesheet.csv \\
            --references  'refs/*.{fa,fasta,fna}{,.gz}' \\
            --kraken2_db  /path/to/kraken2_db \\
            --outdir      results \\
            -profile      slurm,singularity

    Required:
        --samplesheet     CSV with columns: sample,fastq
        --references      Glob OR comma-separated list of reference FASTAs.
                          Accepts .fa .fasta .fna and their .gz variants.

    Optional:
        --kraken2_db      Pre-built Kraken2 database directory. If omitted,
                          taxonomic profiling is skipped.
        --run_bracken     Run Bracken on Kraken2 reports [default: true]
        --bracken_read_length  Read length for Bracken [default: 10000]
        --bracken_level   Tax level (D,P,C,O,F,G,S) [default: S]
        --bracken_threshold    Min reads per taxon [default: 10]
        --minimap2_preset Minimap2 preset [default: map-hifi]
        --minimap2_args   Extra minimap2 args (quoted string) [default: '']
        --outdir          Output directory [default: ./results]
        --publish_dir_mode How to publish outputs (copy|link|symlink) [default: copy]

    Profiles:
        -profile slurm,singularity   Submit jobs via SLURM with Singularity (zurada HPC)
        -profile slurm,apptainer     Same but Apptainer
        -profile local,docker        Run on the current machine with Docker
        -profile test                Tiny smoke-test (uses bundled config)

    Examples:
        # Standard run on the zurada HPC
        nextflow run main.nf \\
            --samplesheet samples.csv \\
            --references  'genomes/*.fasta.gz' \\
            --kraken2_db  /scratch/c0murr09/kraken2/k2_pluspf_20240605 \\
            -profile slurm,singularity

        # No taxonomy (mapping + QC only)
        nextflow run main.nf \\
            --samplesheet samples.csv \\
            --references  ref1.fa,ref2.fa \\
            -profile local,docker
    ========================================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ========================================================================================
// VALIDATE INPUTS
// ========================================================================================

if (!params.samplesheet) { exit 1, "ERROR: --samplesheet is required.  Run with --help for usage." }
if (!params.references)  { exit 1, "ERROR: --references is required.  Run with --help for usage." }

log.info ""
log.info "========================================================================================"
log.info " pb-hifi-mapping-nf"
log.info "========================================================================================"
log.info " Samplesheet     : ${params.samplesheet}"
log.info " References      : ${params.references}"
log.info " Kraken2 DB      : ${params.kraken2_db ?: '(skipped)'}"
log.info " Run Bracken     : ${params.run_bracken && params.kraken2_db}"
log.info " Minimap2 preset : ${params.minimap2_preset}"
log.info " Output dir      : ${params.outdir}"
log.info " Profile         : ${workflow.profile}"
log.info "========================================================================================"
log.info ""

// ========================================================================================
// IMPORT MODULES
// ========================================================================================

include { FASTQC        } from './modules/local/fastqc.nf'
include { MINIMAP2_MAP  } from './modules/local/minimap2_map.nf'
include { KRAKEN2       } from './modules/local/kraken2.nf'
include { BRACKEN       } from './modules/local/bracken.nf'
include { MULTIQC       } from './modules/local/multiqc.nf'
include { SAMPLE_REPORT } from './modules/local/sample_report.nf'

// ========================================================================================
// HELPERS
// ========================================================================================

// Build the reference channel. The user may pass a glob ("refs/*.fa") OR a comma-
// separated list of paths, OR a single path. Each is normalized to (ref_id, file)
// where ref_id is the basename minus its FASTA extension(s).
def buildReferenceChannel(spec) {
    def fasta_exts = ~/(?i)\.(fa|fasta|fna)(\.gz)?$/

    // Comma-separated list?
    def items = spec.toString().contains(',') \
        ? spec.toString().split(',').collect { it.trim() }.findAll { it } \
        : [ spec.toString() ]

    // Each item may itself be a glob - resolve via Channel.fromPath(...)
    return Channel
        .fromList(items)
        .flatMap { pattern -> file(pattern, checkIfExists: false).with { f ->
            f instanceof List ? f : [f]
        } }
        .filter { it.exists() }
        .filter { it.name ==~ /(?i).*\.(fa|fasta|fna)(\.gz)?$/ }
        .map    { f ->
            def ref_id = f.name.replaceAll(fasta_exts, '')
            tuple(ref_id, f)
        }
        .ifEmpty { exit 1, "ERROR: no FASTA files matched --references '${spec}'.  Accepted extensions: .fa .fasta .fna (optionally .gz)" }
}

// ========================================================================================
// WORKFLOW
// ========================================================================================

workflow {

    // ----------------------------------------------------------------------------------
    // 1.  Read samplesheet  →  channel of (meta, reads_file)
    // ----------------------------------------------------------------------------------
    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq) {
                exit 1, "ERROR: samplesheet row missing 'sample' or 'fastq': ${row}"
            }
            def meta = [ sample: row.sample.toString().trim() ]
            def reads = file(row.fastq.toString().trim(), checkIfExists: true)
            tuple(meta, reads)
        }

    // ----------------------------------------------------------------------------------
    // 2.  Reference channel  →  (ref_id, ref_file)
    // ----------------------------------------------------------------------------------
    refs_ch = buildReferenceChannel(params.references)

    // ----------------------------------------------------------------------------------
    // 3.  FastQC (per sample)
    // ----------------------------------------------------------------------------------
    FASTQC(samples_ch)

    // ----------------------------------------------------------------------------------
    // 4.  Cartesian product (sample × reference) for parallel mapping
    // ----------------------------------------------------------------------------------
    map_input_ch = samples_ch
        .combine(refs_ch)
        .map { meta, reads, ref_id, ref ->
            // augment meta with ref_id so the BAM filename is unique per pair
            def meta2 = meta + [ ref_id: ref_id ]
            tuple(meta2, reads, ref)
        }

    MINIMAP2_MAP(map_input_ch)

    // ----------------------------------------------------------------------------------
    // 5.  Taxonomic profile (Kraken2 + Bracken) - skipped if --kraken2_db not set
    // ----------------------------------------------------------------------------------
    kraken_db_ch = params.kraken2_db \
        ? Channel.value( file(params.kraken2_db, checkIfExists: true) ) \
        : Channel.empty()

    if (params.kraken2_db) {

        KRAKEN2(samples_ch, kraken_db_ch)

        if (params.run_bracken) {
            BRACKEN(KRAKEN2.out.report, kraken_db_ch)
            bracken_out_ch = BRACKEN.out.bracken
        } else {
            bracken_out_ch = Channel.empty()
        }

        kraken_out_ch = KRAKEN2.out.report

    } else {
        kraken_out_ch  = Channel.empty()
        bracken_out_ch = Channel.empty()
    }

    // ----------------------------------------------------------------------------------
    // 6.  Per-sample HTML report
    //     Collect all per-sample stat files (across all references) + tax outputs
    //     into one bundle keyed by meta.sample.
    // ----------------------------------------------------------------------------------

    // strip ref_id from meta so we can group by sample
    flagstat_per_sample = MINIMAP2_MAP.out.flagstat
        .map { meta, f -> tuple([sample: meta.sample], f) }
        .groupTuple(by: 0)

    idxstats_per_sample = MINIMAP2_MAP.out.idxstats
        .map { meta, f -> tuple([sample: meta.sample], f) }
        .groupTuple(by: 0)

    // Combine all per-sample input files into a single list for SAMPLE_REPORT.
    // Use .join() so we only emit when both flagstat AND idxstats are ready,
    // and merge in kraken/bracken (or empty markers) by sample.
    report_bundle_ch = flagstat_per_sample
        .join(idxstats_per_sample, by: 0)
        .map { meta, flagstats, idxstats -> tuple(meta.sample, meta, flagstats, idxstats) }

    // Pull kraken / bracken into per-sample lists keyed by sample string.
    // If upstream is empty (no DB / bracken disabled), the channel stays empty;
    // .join(..., remainder: true) below will pass mapping rows through with nulls
    // in those slots.
    kraken_keyed  = kraken_out_ch .map { meta, f -> tuple(meta.sample, [f]) }
    bracken_keyed = bracken_out_ch.map { meta, f -> tuple(meta.sample, [f]) }

    // Left-join taxonomy onto mapping stats - keep sample even when tax missing
    report_input_ch = report_bundle_ch
        .map { sample, meta, flagstats, idxstats -> tuple(sample, meta, flagstats, idxstats) }
        .join(kraken_keyed,  by: 0, remainder: true)
        .join(bracken_keyed, by: 0, remainder: true)
        .map { row ->
            def sample    = row[0]
            def meta      = row[1]
            def flagstats = row[2]
            def idxstats  = row[3]
            def kraken    = row[4] ?: []
            def bracken   = row[5] ?: []
            def all_files = (flagstats + idxstats + kraken + bracken).flatten().findAll { it != null }
            tuple(meta, all_files)
        }

    SAMPLE_REPORT(
        report_input_ch,
        Channel.value( file("${projectDir}/bin/make_sample_report.py", checkIfExists: true) )
    )

    // ----------------------------------------------------------------------------------
    // 7.  MultiQC across the whole run
    // ----------------------------------------------------------------------------------
    multiqc_input_ch = Channel.empty()
        .mix( FASTQC.out.zip.map        { meta, f -> f } )
        .mix( MINIMAP2_MAP.out.flagstat .map { meta, f -> f } )
        .mix( MINIMAP2_MAP.out.stats    .map { meta, f -> f } )
        .mix( MINIMAP2_MAP.out.idxstats .map { meta, f -> f } )
        .mix( kraken_out_ch             .map { meta, f -> f } )
        .collect()

    multiqc_config_ch = file("${projectDir}/assets/multiqc_config.yaml", checkIfExists: true)

    MULTIQC(multiqc_input_ch, multiqc_config_ch)
}

// ========================================================================================
// COMPLETION
// ========================================================================================

workflow.onComplete {
    log.info ""
    log.info "========================================================================================"
    log.info " pb-hifi-mapping-nf  -  COMPLETE"
    log.info "========================================================================================"
    log.info " Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info " Duration : ${workflow.duration}"
    log.info " Outdir   : ${params.outdir}"
    log.info ""
    log.info " Key outputs:"
    log.info "   ${params.outdir}/mapping/<sample>/                 sorted+indexed BAMs (one per reference)"
    log.info "   ${params.outdir}/qc/fastqc/                        FastQC reports"
    log.info "   ${params.outdir}/qc/multiqc/multiqc_report.html    Aggregated MultiQC"
    log.info "   ${params.outdir}/taxonomy/kraken2/                 Kraken2 reports + classifications"
    log.info "   ${params.outdir}/taxonomy/bracken/                 Bracken abundance estimates"
    log.info "   ${params.outdir}/reports/<sample>.report.html      Per-sample summary"
    log.info "   ${params.outdir}/pipeline_info/                    Trace, timeline, DAG, report"
    log.info "========================================================================================"
}

workflow.onError {
    log.info ""
    log.info "Pipeline failed: ${workflow.errorMessage}"
}
