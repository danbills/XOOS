<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to STR Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0]

### Added

- Genotyping and detection of repeat expansions in short tandem repeats (STR) from SBX sequencing data.
- VCF and JSON output describing alternate repeat structures.
- Metrics TSV output with per-locus statistics.
- STR catalog input via `--str-catalog` for defining target loci.
- Support for SBX duplex reads with lossless encoding to recover R1 and R2 reads for additional support.
- Configurable region padding, read filtering, and confidence thresholds.
- FXN locus support with configurable spanning-read fraction cut-offs for HOM and low-confidence calls.
- Support for compressed reference files.

<!-- Version comparison links -->
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.1.0
