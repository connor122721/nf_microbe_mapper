# pb-hifi-mapping-nf

A small, robust [Nextflow](https://www.nextflow.io/) DSL2 pipeline for mapping
PacBio HiFi reads against one or more reference genomes in parallel, with
read-level QC and a taxonomic profile per sample.

```
   samplesheet.csv   ──┐
                       │
   refs/*.{fa,fasta,fna}{,.gz} ──► minimap2 (sample × reference)
                       │             ├─► sorted/indexed BAM (per pair)
                       │             └─► samtools flagstat / stats / idxstats
                       │
                       ├─► FastQC                              (per sample)
                       ├─► Kraken2  → Bracken                  (per sample)
                       ├─► Per-sample HTML report              (per sample)
                       └─► MultiQC                             (aggregate)
```

The minimap2 + iterative-mapping pattern was inspired by the
[PacificBiosciences/pb-metagenomics-tools](https://github.com/PacificBiosciences/pb-metagenomics-tools)
collection; this pipeline is intentionally simpler — one BAM per
`(sample × reference)` pair, plus Kraken2/Bracken for taxonomy.

---

## Quick start

```bash
# 1.  Get the code
git clone https://github.com/connor122721/nf_microbe_mapper.git
cd nf_microbe_mapper

# 2.  Edit the samplesheet
cp assets/samplesheet.csv my_samples.csv
nano my_samples.csv

# 3.  Run on the zurada HPC
nextflow run main.nf \
    --samplesheet my_samples.csv \
    --references 'refs/*.fasta.gz' \
    --kraken2_db  /work/c0murr09/kraken2/k2_pluspf_20240605 \
    --outdir ./results \
    -profile slurm,apptainer
```
---

## Inputs

### Samplesheet (`--samplesheet`)

A CSV with two required columns:

```csv
sample,fastq
sample01,/data/runs/sample01.hifi_reads.fastq.gz
sample02,/data/runs/sample02.hifi_reads.fastq.gz
```

* `sample` — short, unique sample ID (used in filenames).
* `fastq`  — absolute path to the HiFi reads.  `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz` accepted.

Validate before running:

```bash
python bin/check_samplesheet.py my_samples.csv
```

### References (`--references`)

Pass either a glob or a comma-separated list.  Recognized extensions:
`.fa  .fasta  .fna  .fa.gz  .fasta.gz  .fna.gz`.

```bash
--references 'genomes/*.fasta.gz'
--references ref1.fa,ref2.fasta,ref3.fna.gz
```

The basename (minus FASTA extension) becomes the **`ref_id`** that appears in
output filenames, e.g. `sample01.vs.GRCh38.sorted.bam`.

### Taxonomic database (`--kraken2_db`)

A pre-built Kraken2 database **directory**.  We recommend the standard
[`k2_pluspf`](https://benlangmead.github.io/aws-indexes/k2) build for HiFi
metagenomic samples.  Bracken's `.kmer_distrib` files ship inside the same
directory, so no separate Bracken DB is needed.

If `--kraken2_db` is omitted, taxonomic profiling is skipped and the per-sample
report contains only mapping stats.

---

## Outputs

```
results/
├── mapping/
│   └── <sample>/
│       ├── <sample>.vs.<ref_id>.sorted.bam        ← primary deliverable
│       ├── <sample>.vs.<ref_id>.sorted.bam.bai
│       ├── <sample>.vs.<ref_id>.flagstat.txt
│       ├── <sample>.vs.<ref_id>.stats.txt
│       └── <sample>.vs.<ref_id>.idxstats.txt
├── qc/
│   ├── fastqc/<sample>/                           ← FastQC HTML + zip
│   └── multiqc/multiqc_report.html                ← aggregated report
├── taxonomy/
│   ├── kraken2/
│   │   ├── <sample>.kraken2.report.txt
│   │   └── <sample>.kraken2.out.txt
│   └── bracken/
│       ├── <sample>.bracken.tsv
│       └── <sample>.bracken_report.txt
├── reports/
│   ├── <sample>.report.html                       ← per-sample summary
│   └── <sample>.summary.tsv
└── pipeline_info/
    ├── execution_report.html
    ├── timeline.html
    ├── trace.tsv
    └── dag.html
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet`         | _required_ | CSV with `sample,fastq` columns |
| `--references`          | _required_ | Glob or comma list of FASTAs |
| `--kraken2_db`          | `null` | Kraken2 DB directory (skips taxonomy if unset) |
| `--run_bracken`         | `true` | Run Bracken on Kraken2 reports |
| `--bracken_read_length` | `300` | HiFi-typical read length for Bracken |
| `--bracken_level`       | `S` | Tax level (D/P/C/O/F/G/S) |
| `--bracken_threshold`   | `10` | Min reads per taxon |
| `--minimap2_preset`     | `map-hifi` | Minimap2 `-ax` preset |
| `--minimap2_args`       | `''` | Extra minimap2 args (quoted) |
| `--outdir`              | `./results` | Output root |
| `--publish_dir_mode`    | `copy` | `copy`, `link`, or `symlink` |
| `--partition`           | `cpu384g` | SLURM partition (zurada default) |

Per-process container images are pinned in each `modules/local/*.nf` file.
Override any of them with the matching `--<tool>_container <image>` parameter
(e.g. `--minimap2_container /home/c0murr09/sif/minimap2.sif`).

---

## Profiles

| Profile | Effect |
|---|---|
| `slurm`         | Submit jobs through SLURM (zurada `cpu384g` partition) |
| `local`         | Run on the current machine |
| `singularity`   | Use Singularity to run containers |
| `apptainer`     | Use Apptainer to run containers |
| `docker`        | Use Docker to run containers |
| `podman`        | Use Podman to run containers |
| `test`          | Tiny smoke-test config (see `conf/test.config`) |
| `debug`         | Disable cleanup, dump task hashes |

Combine profiles with commas: `-profile slurm,apptainer`.

---

## Containers

All containers are pulled from public registries.  Defaults (set in each module):

| Tool | Image |
|---|---|
| FastQC          | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| minimap2 + samtools | `quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0` |
| Kraken2         | `quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_4` |
| Bracken         | `quay.io/biocontainers/bracken:3.0--h9948957_0` |
| MultiQC         | `quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0` |
| Python (report) | `docker.io/library/python:3.11-slim` |

If your HPC is offline-only, set `NXF_SINGULARITY_CACHEDIR` (or
`NXF_APPTAINER_CACHEDIR`) so Nextflow caches one copy of each `.sif`, and run
once on a node with internet access to populate it.

---

## Resource configuration

Per-process resources are controlled by **labels** in `conf/base.config`:

* `process_low`         – 2 CPU / 8 GB / 2 h
* `process_medium`      – 8 CPU / 16 GB / 4 h
* `process_high`        – 16 CPU / 32 GB / 8 h  (minimap2)
* `process_high_memory` – 8 CPU / 80 GB / 6 h  (Kraken2 — DB lives in RAM)

Edit those numbers to match your queue.  The `MINIMAP2_MAP` process is capped
at `maxForks = 16` so a 100-sample × 5-reference run doesn't flood the
scheduler — bump it if your cluster can absorb more.

---

## Recipe: building / sourcing a Kraken2 DB

The most popular pre-built DB for HiFi work is the
[Langmead "PlusPF" index](https://benlangmead.github.io/aws-indexes/k2)
(bacteria + archaea + viral + protozoa + fungi + human).  It is ~80 GB on disk
and ships with the Bracken `.kmer_distrib` files.

```bash
mkdir -p k2_pluspf && cd k2_pluspf
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz
tar -xzvf k2_standard_08_GB_20260226.tar.gz
# point the pipeline at this directory:
nextflow run ... --kraken2_db $PWD
```

---

## Citation

If you use this pipeline, please also cite the underlying tools:

* Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34: 3094–3100.
* Wood, D. E., Lu, J., Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology* 20: 257.
* Lu, J. _et al._ (2017). Bracken: estimating species abundance in metagenomics samples. *PeerJ Computer Science* 3: e104.
* Andrews S. (2010). FastQC: A quality control tool for high throughput sequence data.
* Ewels, P. _et al._ (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 32: 3047–3048.

---

## Author

**Connor S. Murray, PhD** — Bioinformatician, University of Louisville School of Medicine.

License: MIT (see `LICENSE`).
