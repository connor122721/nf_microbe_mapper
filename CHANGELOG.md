# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-04-28

### Added
- Initial release.
- Samplesheet-driven HiFi mapping with `minimap2 -ax map-hifi`.
- Parallel `(sample × reference)` mapping; one sorted/indexed BAM per pair.
- FastQC + samtools flagstat/stats/idxstats per BAM.
- Optional Kraken2 + Bracken taxonomic profiling.
- Per-sample HTML report (mapping rates + top taxa).
- MultiQC aggregation.
- SLURM, local, Singularity, Apptainer, Docker, Podman profiles.
- `bin/check_samplesheet.py` validator.
