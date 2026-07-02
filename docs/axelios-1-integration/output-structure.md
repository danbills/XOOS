<!-- markdownlint-disable MD024 -->
# Sequencer output directory structure

This page describes the standard storage folder structure for sequencing
data produced on the AXELIOS platform.

## Overall structure

Each AXELIOS platform writes all output under a single top-level
`customer-storage/` folder on the network storage volume root:

```text
<volume-root>/
└── customer-storage/
    ├── runs/             ← Sequencing run output
    ├── analyses/         ← Analysis task output
    ├── problem-reports/  ← Diagnostic packages
    ├── logs/             ← Platform logs
    └── resources/        ← User-provided reference files
```

`<volume-root>` is the NFS mount root for one deployed platform
instance, for example `roche-seq-platform-1/`.

| Folder             | Purpose                                                                          |
| ------------------ | -------------------------------------------------------------------------------- |
| `runs/`            | Output data from sequencing runs and primary analysis (basecall and demux data). |
| `analyses/`        | Results of analysis tasks (for example, RDB-to-FASTQ conversion and alignment).  |
| `problem-reports/` | Encrypted diagnostic packages used for troubleshooting.                          |
| `logs/`            | Audit, event, and service logs.                                                  |
| `resources/`       | User-provided reference files reused across runs and analyses.                   |

## Accessing the storage

The storage volume is served over NFS.
Mount the network share at `<volume-root>/` on any system that needs to read
the output; all platform output lives under `<volume-root>/customer-storage/`.

Grant read access to downstream consumers (for example, analysis servers and
data pipelines) and restrict write access to the platform itself.

## `runs/` folder

`runs/` is the repository for all output data generated from sequencing runs
and primary analysis.

### Naming convention

Each run resides in a dedicated subfolder:

```text
<YYYYMMDD>_<device-name>_Q<queue-number-for-the-day>_R<run-number-in-the-queue>_<workflow-order-name>_<synthesis-output-tubeId>/
```

| Token                          | Description                                                    | Example        |
| ------------------------------ | -------------------------------------------------------------- | -------------- |
| `<YYYYMMDD>`                   | Date the run started.                                          | `20250609`     |
| `<device-name>`                | User-assigned identifier of the sequencing instrument.         | `HTP-Seq1`     |
| `Q<queue-number-for-the-day>`  | Queue number for the run on that day.                          | `Q1`           |
| `R<run-number-in-the-queue>`   | Run number within that queue.                                  | `R2`           |
| `<workflow-order-name>`        | Workflow order name assigned in the platform management software (SPM). | `KinaseStudy3` |
| `<synthesis-output-tubeId>`    | Human-readable portion of the synthesis output tube.           | `987654B`      |

Example run folder path:
`roche-seq-platform-1/customer-storage/runs/20250609_HTP-Seq1_Q1_R2_KinaseStudy3_987654B/`

### Contents

A run folder holds the sequencing data produced during the run, plus a workflow
configuration file.
Only a defined subset of the sequencer's output files is included.

The sequencing data is delivered in RDB (Roche Data Block) format.
A run operates in one of two modes — raw basecall reads or demultiplexed
consensus reads — and the output layout differs accordingly.
FASTQ files are not present in the run folder; FASTQ is produced separately by
the RDB-to-FASTQ conversion analysis task and lands under `analyses/`.

#### Output directory layout

For both modes, sequencing data files are nested under a `capture/` directory
at the root of the run folder.
Output is partitioned by GPU into four subdirectories — `device_0` through
`device_3` (zero-based index) — each containing a `main/` subdirectory with
the RDB data files.
A run-level metrics file, *run_key_demux_metrics.csv*, is written directly
under `capture/`.

#### Basecall (raw read) output

| Item         | Detail                                                                                                             |
| ------------ | ------------------------------------------------------------------------------------------------------------------ |
| File pattern | `PrimaryAnalysis_<M>.rdb`                                                                                          |
| `<M>`        | Six-digit, zero-padded, one-based time-chunk index (for example, `000001`).                                        |
| Content      | Base-called molecular traces and quality scores for one fixed-duration time chunk; each chunk is independent.      |

```text
<run-folder>/
└── capture/
    ├── device_0/
    │   └── main/
    │       ├── PrimaryAnalysis_000001.rdb
    │       └── PrimaryAnalysis_000002.rdb
    ├── device_1/ … device_3/   (same structure)
    └── run_key_demux_metrics.csv
```

#### Demultiplexed consensus output

| Item         | Detail                                                                                                             |
| ------------ | ------------------------------------------------------------------------------------------------------------------ |
| File pattern | `DemuxAnalysis_<M>.rdb`                                                                                            |
| `<M>`        | Six-digit, zero-padded, one-based time-chunk index.                                                                |
| Content      | Reads assigned to samples (by sample identification) for one time chunk; one RDB per time chunk per GPU device.     |
| Extra file   | `DemuxParams.rdb` — written under `capture/` at run start, holding global demux metadata.                          |

A run produces either `PrimaryAnalysis_*` (basecall) or `DemuxAnalysis_*`
(consensus) files depending on the configured workflow, not both.

```text
<run-folder>/
└── capture/
    ├── device_0/
    │   └── main/
    │       ├── DemuxAnalysis_000001.rdb
    │       └── DemuxAnalysis_000002.rdb
    ├── device_1/ … device_3/   (same structure)
    ├── DemuxParams.rdb
    └── run_key_demux_metrics.csv
```

#### Ancillary run files

The run folder also includes declaration, metadata, and metrics files:

| File                                                    | Description                                              |
| ------------------------------------------------------- | -------------------------------------------------------- |
| `AXELIOS-1-Platform_metrics-glossary_customer-facing.txt` | Glossary of platform metrics.                          |
| `{RunID}.summary.json`                                  | Run summary, where `{RunID}` is the run identifier.      |
| `run_key_quality_metrics.csv`                           | Per-run-key quality metrics.                             |
| `failed_psp_run_key_quality_metrics.csv`                | Quality metrics for failed probe-set pool (PSP) run keys. |

#### Workflow configuration file

*workflow-configuration.yaml* captures the parameters of the synthesis,
sequencing, and analysis steps for use by downstream tools.

```text
roche-seq-platform-1/customer-storage/runs/20250609_HTP-Seq1_Q1_R2_KinaseStudy3_987654B/
├── capture/
│   ├── device_0/
│   │   └── main/
│   │       └── <PrimaryAnalysis_* or DemuxAnalysis_*>.rdb
│   └── … device_1 … device_3 (same structure)
└── workflow-configuration.yaml
```

## `analyses/` folder

Each analysis folder stores the results of an analysis task.
The types, content, and filenames of results depend on the analysis task that
produced them.
A *workflow-configuration.yaml* file is always included.

### Naming convention

```text
<analysis-task-name>_<YYYYMMDD>_<device-name>_A<analysis-number-for-the-day>/
```

| Token                            | Description                                                                                              |
| -------------------------------- | -------------------------------------------------------------------------------------------------------- |
| `<analysis-task-name>`           | Workflow order name for tasks initiated automatically by the platform; otherwise specified by the user via the Analysis Controller CLI (AC-CLI). |
| `<YYYYMMDD>`                     | Date the analysis task started.                                                                          |
| `<device-name>`                  | Analysis device/system where the analysis ran (for example, `sas-1`).                                   |
| `A<analysis-number-for-the-day>` | Sequential analysis number for that day, starting at 1 (for example, `A1`).                             |

Example analyses folder path:
`roche-seq-platform-1/customer-storage/analyses/KinaseStudy3_20250609_sas-1_A1/`

### Contents

The result files depend on the analysis task type.
The two task types are RDB-to-FASTQ conversion and alignment.
A *workflow-configuration.yaml* file is always present.

Common conventions across RDB-to-FASTQ output filenames:

- **`<analysis-task-name>`** — the name of the analysis task that produced the
  output.
  Allowed characters: alphanumeric, dashes, underscores; non-empty; maximum
  50 characters.
- **`<part-index>`** — FASTQ output is split across many files to ease
  transfer and parallel processing.
  The part index is a zero-padded six-digit integer starting at `000001`.
  Each part is approximately 2 GB, except the final part of a task, which may
  be smaller.

#### RDB-to-FASTQ conversion output

The conversion task takes RDB input and produces FASTQ.
You select compressed or uncompressed output, which changes the file
extension (*.fastq.gz* vs *.fastq*).

**Basecall FASTQ** (from basecall RDB input) — all reads from the input,
re-chunked independently of the RDB chunking:

| Compression       | Filename pattern                                   |
| ----------------- | -------------------------------------------------- |
| 1 or above (gzip) | `<analysis-task-name>-<part-index>.fastq.gz`       |
| 0 (uncompressed)  | `<analysis-task-name>-<part-index>.fastq`          |

**Consensus FASTQ** (from consensus RDB input) — chunked FASTQ files produced
per detected sample identity (SID).
The SID nucleotide sequence is embedded in the filename via the
`<SID-sequence>` token (uppercase IUPAC codes A/C/G/T).
Reads with no recognized SID use the literal token `Undetermined`:

| Compression       | Filename pattern                                                         |
| ----------------- | ------------------------------------------------------------------------ |
| 1 or above (gzip) | `<analysis-task-name>-<SID-sequence>-<part-index>.fastq.gz`              |
| 0 (uncompressed)  | `<analysis-task-name>-<SID-sequence>-<part-index>.fastq`                 |

#### Alignment output

Alignment tasks produce aligned data and index files (for example, BAM and
associated indices).
See the task-specific output specification for per-file format details.

#### Workflow configuration file

*workflow-configuration.yaml* contains synthesis, sequencing, and analysis
parameters for downstream tools.
It is included in every analysis output folder.

#### Examples

**Compressed basecall FASTQ:**

```text
KinaseStudy5_20250609_SAS-399_A2/
├── KinaseStudy5-000001.fastq.gz
├── KinaseStudy5-000002.fastq.gz
├── KinaseStudy5-000003.fastq.gz
└── workflow-configuration.yaml
```

**Compressed consensus FASTQ** (per SID, plus `Undetermined`):

```text
Trio187_20250616_SAS-399_A1/
├── Trio187-ATGCAGATA-000001.fastq.gz
├── Trio187-ATGCAGATA-000002.fastq.gz
├── Trio187-GGAATCGTT-000001.fastq.gz
├── Trio187-CTGTCCAAG-000001.fastq.gz
├── Trio187-Undetermined-000001.fastq.gz
└── workflow-configuration.yaml
```

## `problem-reports/` folder

`problem-reports/` stores encrypted packages of logs and diagnostic
information used by Roche Service for troubleshooting.

### Layout

Sequencer problem reports are written directly under `problem-reports/` as
`PR_[A/M/Q]_*.package.zip` files.
Analysis server (SAS) problem reports are organized into per-device subfolders
named `<sas-device-name>/`, where `<sas-device-name>` is the user-assigned
identifier of the originating SAS device (for example, `sas-1`).
The naming convention for individual package files typically includes the date
and device identifiers.

These packages are encrypted and intended for Roche Service Representatives
for troubleshooting only.

The sequencer uses the following package-name conventions:

| Pattern                | Trigger                                                          |
| ---------------------- | --------------------------------------------------------------- |
| `PR_A_*.package.zip`   | Automatic (abnormal shutdown / system-issue reports).           |
| `PR_M_*.package.zip`   | Manual / on-demand (triggered via Service software).            |
| `PR_Q_*.package.zip`   | Run-queue report (auto-generated at the end of each run queue). |

```text
roche-seq-platform-1/customer-storage/problem-reports/
├── PR_A_xxx.package.zip
├── PR_M_xxx.package.zip
├── PR_Q_xxx.package.zip
└── sas-1/
    └── 2025052955_sas-1_problemreport.zip
```

## `logs/` folder

`logs/` stores customer-facing platform logs — audit logs, event logs, and
service logs.
These include summaries of workflow execution, user activity, and other
platform events.

Logs are organized by log type into separate subfolders (`AuditLog/`,
`EventLog/`, `ServiceLog/`).
Log files are gzip-compressed (`.log.gz`):

| Subfolder     | Content                                                                                |
| ------------- | ------------------------------------------------------------------------------------- |
| `AuditLog/`   | Audit logs capturing user activities.                                                 |
| `EventLog/`   | Event logs capturing system events with severity.                                     |
| `ServiceLog/` | Service logs capturing maintenance status, actions, and other service-related activity. |

```text
roche-seq-platform-1/customer-storage/logs/
├── AuditLog/
│   └── AuditLog_GTIN-07613336223307-361_xxx.log.gz
├── EventLog/
│   └── EventLog_GTIN-07613336223307-361_xxx.log.gz
└── ServiceLog/
    └── ServiceLog_GTIN-07613336223307-361_xxx.log.gz
```

## `resources/` folder

`resources/` holds user-provided reference files used by sequencing and
analysis workflows.
These are typically static files reused across multiple runs and analyses, such
as reference genomes (FASTA) and target-region definitions (BED).

```text
roche-seq-platform-1/customer-storage/resources/
├── homo_sapiens.GRCh38.dna.alt.fa
└── Target_regions.bed
```

## General notes and best practices

- **Immutability of run folders:** Once a run completes and its data is
  published to `runs/`, the folder contents are considered immutable.
  Re-processing produces a new run or analysis folder rather than modifying
  existing data.
- **Permissions:** Apply file-system permissions to protect data integrity.
  Grant read-only access to downstream consumers (for example, analysis
  pipelines and data integration tools) and restrict write access to the
  platform itself.
- **Local network storage is required:** Network storage must be provisioned
  before the platform can operate.
  Cloud storage (for example, Google Cloud Storage), if configured, can be an
  additional destination for analysis output but does not replace local network
  storage.
- **Scalability:** The hierarchy is designed to scale with data volume and run
  count; plan capacity primarily around `runs/` and `analyses/`, which grow per
  run.

## Consolidated layout

```text
<volume-root>/
└── customer-storage/
    ├── runs/
    │   └── <YYYYMMDD>_<device-name>_Q<queue-number-for-the-day>_R<run-number-in-the-queue>_<workflow-order-name>_<synthesis-output-tubeId>/
    │       ├── capture/
    │       │   ├── device_0/ … device_3/
    │       │   │   └── main/
    │       │   │       └── <PrimaryAnalysis_* or DemuxAnalysis_*>.rdb
    │       │   └── run_key_demux_metrics.csv
    │       └── workflow-configuration.yaml
    ├── analyses/
    │   └── <analysis-task-name>_<YYYYMMDD>_<device-name>_A<analysis-number-for-the-day>/
    │       ├── <Files produced by analysis tasks>
    │       └── workflow-configuration.yaml
    ├── problem-reports/
    │   ├── PR_[A/M/Q]_*.package.zip
    │   └── <sas-device-name>/
    │       └── <filename>.zip
    ├── logs/
    │   ├── AuditLog/
    │   │   └── <filename>.log.gz
    │   ├── EventLog/
    │   │   └── <filename>.log.gz
    │   └── ServiceLog/
    │       └── <filename>.log.gz
    └── resources/
        └── <resource filename>
```

Concrete example:

```text
roche-seq-platform-1/
└── customer-storage/
    ├── runs/
    │   └── 20250609_HTP-Seq1_Q1_R2_KinaseStudy3_987654B/
    │       ├── capture/
    │       │   ├── device_0/
    │       │   │   └── main/
    │       │   │       ├── PrimaryAnalysis_000001.rdb
    │       │   │       └── PrimaryAnalysis_000002.rdb
    │       │   └── run_key_demux_metrics.csv
    │       └── workflow-configuration.yaml
    ├── analyses/
    │   └── KinaseStudy3_20250609_sas-1_A1/
    │       ├── KinaseStudy3-000001.fastq.gz
    │       ├── KinaseStudy3-000002.fastq.gz
    │       └── workflow-configuration.yaml
    ├── problem-reports/
    │   ├── PR_A_xxx.package.zip
    │   ├── PR_M_xxx.package.zip
    │   ├── PR_Q_xxx.package.zip
    │   └── sas-1/
    │       └── 2025052955_sas-1_problemreport.zip
    ├── logs/
    │   ├── AuditLog/
    │   │   └── AuditLog_GTIN-07613336223307-361_xxx.log.gz
    │   ├── EventLog/
    │   │   └── EventLog_GTIN-07613336223307-361_xxx.log.gz
    │   └── ServiceLog/
    │       └── ServiceLog_GTIN-07613336223307-361_xxx.log.gz
    └── resources/
        ├── homo_sapiens.GRCh38.dna.alt.fa
        └── Target_regions.bed
```
