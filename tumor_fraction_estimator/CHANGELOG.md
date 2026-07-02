<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Tumor Fraction Estimator will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0]

### Added

- Tumor fraction estimation from SBX duplex sequencing data via the `tumor-fraction` subcommand.
- Contamination estimation via the `contamination` subcommand using population SNP VCFs.
- TSV metrics output with per-subcommand prefixed filenames (`tumor_fraction.results.tsv`, `contamination.results.tsv`).
- Input validation with file extension and readability checks, and fail-fast errors when output files already exist.
- Read filtering controls grouped into read-level, base-level, pileup, and noise-estimation options.
- `--save-read-features` developer option for emitting per-read feature data.

<!-- Version comparison links -->
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.1.0
