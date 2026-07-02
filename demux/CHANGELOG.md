<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Demux will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0]

### Added

- Added `--read-length-mode all-split` option that writes full and partial simplex reads into separate subdirectories (`<sample>/full/` and `<sample>/partial/`).
- Added `--suppress-simplex-qual-override` flag to preserve original quality scores for simplex adapter reads instead of overriding with a fixed value.
- Added new simplex adapter type (`YS-NEW`) for improved SID assignment accuracy.
- Added new `SIMPLEX-10X` adapter type for 10x-linked simplex libraries.
- Added `--min-score` CLI parameter to set a minimum log-odds score threshold for adapter matching (default `30`, experimental).
- Added metrics `raw_bases`, `unassigned_bases`, `failed_assigned_bases`, and `trimmed_bases` to duplex metrics.
- Added sample sheet validation to reject the reserved sample name `raw_failed`.
- Added `file_provenance.tsv` in each sample output directory, mapping each output FASTQ file to its source input file.
- Added optional `prefix` and `suffix` fields to LUT definitions in the adapter design manifest.
- Added TSV support for `--sample-sheet` with auto-detected delimiter (comma or tab).
- Added `--worker-threads-per-input` CLI option to control concurrent input file processing (default `4`).
- Added optional `bait_5p` and `bait_3p` fields in adapter design manifest.
- Added `--discordant-sid-mode` CLI option to control handling of reads with discordant 5'/3' SIDs (`discard-tied`, `discard-all`, `keep`).
- Added `sid_discordant_discarded_reads` and `failed_assigned_reads` metrics to simplex run metrics.

### Changed

- Changed output FASTQ filename format to `<sample_name>-<part_index>.fastq[.gz|.zst]`. Part indices are 1-based, zero-padded, contiguous, and assigned per sample in input file order.
- Input files are now sorted lexicographically before processing to ensure deterministic part index assignment.
- Migrated TSV/histogram metadata to structured `##RocheCommandLine=<...>` header format.
- Using `NA` instead of `0` for metrics not applicable to duplex UMI adapter types and for strand metrics when strand detection is disabled.
- Sample names are determined from prefixes instead of directory names, supporting nested output directories.
- CLI `--help` output now groups adapter-specific options under "Simplex Options" and "Duplex Options" headings.
- Changed discordant SID behavior of simplex reads to no longer partial trim or trim both SIDs when discordant.

### Fixed

- Fixed `IsPartialRead` for YSU adapter: reads with one SID and zero UMIs are now correctly classified as partial.
- Improved duplex hairpin detection: fixed out-of-bounds search window near read ends, corrected SID search anchoring for indel-containing loop matches, and scaled search windows to loop length.
- Fixed non-deterministic `file_sid` assignment caused by unordered map iteration in sample sheet parsing.
- Fixed `raw_failed` FASTQ files being written to the output root instead of the `raw_failed/` subdirectory.
- Fixed simplex metrics counting: `preassignment_passing_reads` now includes reads filtered after trimming, and `too_short_reads` vs `too_short_trimmed_reads` are tracked separately.
- Fixed YSU false positive UMI matches by extending the LUT match region and gating SID trimming on UMI confirmation.
- Fixed `SeqMatcher::FindNextBarcode` failing to find exact barcode matches when shorter approximate matches existed at adjacent start positions.

### Removed

- Removed unused `stem-5p.fa` and `stem-3p.fa` from `YSU` adapter design bundle.

## [1.0.0]

### Added

- Adapter trimming and sample demultiplexing for SBX reads.
- Adapter designs: SBX-D (duplex), SBX-DM (TAPS+ methylation with XM tag output), SBX-FAST (simplex), YS.
- SID identification and UMI extraction from adapter sequences.
- YC tag annotation (v1.5) with base-type quality scores.
- Duplex and simplex read processing with configurable read-length filtering.
- FASTQ input validation for sequence/quality length consistency and valid characters.
- Compressed output via gzip (default) or zstd.
- Per-sample output directories with configurable writing threads.
- Run-level and sample-level metrics output (`run_stats.tsv`, `sample_stats.tsv`).
- Strand index generation via `demux_strand_index`.

<!-- Version comparison links -->
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.0.0
