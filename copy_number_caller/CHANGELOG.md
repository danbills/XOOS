<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Copy Number Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0]

### Added

- Genomic copy number variation detection by analyzing coverage patterns in DNA sequencing alignments.
- Germline WGS workflow (`germline-wgs`) for detecting copy number variants in normal samples.
- Somatic tumor-normal WGS workflow (`tumor-normal-wgs`) for characterizing allele-specific somatic CNVs in matched tumor-normal pairs.
- Segmentation workflow (`segmentation`) for standalone genome segmentation.
- Purity/ploidy estimation with grid search and whole genome doubling (WGD) QC.
- Somatic CNV prediction via `predict-somatic-cnv` subcommand.
- Interval coverage calculation with GC-content normalization.
- Intermediate output files in standard BED format with metadata headers.
- VCF and SEG output with version and command line metadata for traceability.
- Mirrored BAF filter and low interval density filter for purity/ploidy estimation.
- Decoy interval filtering via `hg38_decoy_intervals.bed` to reduce false positive duplications.
- Low-coverage early exit: when median coverage is below 10x, the module writes valid header-only output files.
- Graceful exit on VCF data-quality issues in `tumor-normal-wgs` (low het-SNP fraction, zero variants after filtering).
- Warning and skip for calls in pseudo-autosomal regions (PAR) if chrY PAR alignments are detected.
- Configurable blocklist BED for excluding problematic regions.

<!-- Version comparison links -->
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.1.0
