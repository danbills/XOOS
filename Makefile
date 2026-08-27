# ============================================================================
# XOOS CUDA Acceleration Developer Makefile
# File: Makefile
#
# Provides high-level build, benchmark, and Nix layered container commands.
# ============================================================================

SHELL := /usr/bin/env bash
NIX_BIN := $(shell which nix 2>/dev/null || echo /nix/var/nix/profiles/default/bin/nix)
PODMAN_BIN := $(shell which podman 2>/dev/null || echo podman)

.PHONY: all help build aligner demux docker-aligner docker-demux docker-all \
        stream-aligner stream-demux run-aligner run-demux bench-aligner bench-demux \
        dev check clean

all: help

help:
	@echo "================================================================================"
	@echo "  XOOS CUDA Acceleration Build & Execution System (sm_120 / CUDA 13.x)"
	@echo "================================================================================"
	@echo "  make aligner           - Build CUDA Aligner package via Nix Flake"
	@echo "  make demux             - Build CUDA Demux package via Nix Flake"
	@echo "  make bench-aligner     - Run native CUDA Aligner benchmark (500k reads)"
	@echo "  make bench-demux       - Run native CUDA Demux benchmark (500k reads)"
	@echo ""
	@echo "  [Nix Layered Container Targets]"
	@echo "  make docker-aligner    - Build 16-layer OCI container image for Aligner"
	@echo "  make docker-demux      - Build 16-layer OCI container image for Demux"
	@echo "  make stream-aligner    - Stream & load Aligner image into Podman"
	@echo "  make stream-demux      - Stream & load Demux image into Podman"
	@echo "  make run-aligner       - Execute containerized Aligner in Podman with GPU"
	@echo "  make run-demux         - Execute containerized Demux in Podman with GPU"
	@echo ""
	@echo "  [Development]"
	@echo "  make dev               - Enter Nix development shell (CUDA, GCC, NVCC)"
	@echo "  make check             - Validate Nix Flake outputs"
	@echo "  make clean             - Remove local build symlinks and caches"
	@echo "================================================================================"

# ----------------------------------------------------------------------------
# Compilation Targets
# ----------------------------------------------------------------------------

aligner:
	@echo ">>> Building CUDA Aligner with Nix Flake..."
	$(NIX_BIN) build .#aligner --print-build-logs

demux:
	@echo ">>> Building CUDA Demux with Nix Flake..."
	$(NIX_BIN) build .#demux --print-build-logs

read-collapser:
	@echo ">>> Building CUDA Read Collapser with Nix Flake..."
	$(NIX_BIN) build .#read-collapser --print-build-logs

metrics:
	@echo ">>> Building CUDA Alignment Metrics with Nix Flake..."
	$(NIX_BIN) build .#metrics --print-build-logs

svc:
	@echo ">>> Building CUDA Small Variant Caller with Nix Flake..."
	$(NIX_BIN) build .#svc --print-build-logs

cnv:
	@echo ">>> Building CUDA Copy Number Caller with Nix Flake..."
	$(NIX_BIN) build .#cnv --print-build-logs

str:
	@echo ">>> Building CUDA Short Tandem Repeat (STR) Caller with Nix Flake..."
	$(NIX_BIN) build .#str --print-build-logs

tfe:
	@echo ">>> Building CUDA Tumor Fraction Estimator with Nix Flake..."
	$(NIX_BIN) build .#tfe --print-build-logs

pangenome:
	@echo ">>> Building CUDA Pangenome Consensus Caller with Nix Flake..."
	$(NIX_BIN) build .#pangenome --print-build-logs

fused:
	@echo ">>> Building CUDA Fused Super-Kernel Engine (AOT + JIT) with Nix Flake..."
	$(NIX_BIN) build .#fused --print-build-logs

# ----------------------------------------------------------------------------
# Benchmarking & Profiling Targets
# ----------------------------------------------------------------------------

bench-aligner: aligner
	@echo ">>> Executing Native CUDA Aligner Benchmark..."
	./result/bin/aligner_cuda_bench 500000

bench-demux: demux
	@echo ">>> Executing Native CUDA Demux Benchmark..."
	./result/bin/demux_cuda_bench 500000

bench-read-collapser: read-collapser
	@echo ">>> Executing Native CUDA Read Collapser Benchmark..."
	./result/bin/read_collapser_cuda_bench 500000

bench-metrics: metrics
	@echo ">>> Executing Native CUDA Alignment Metrics Benchmark..."
	./result/bin/alignment_metrics_cuda_bench 1000000

bench-svc: svc
	@echo ">>> Executing Native CUDA Small Variant Caller Benchmark..."
	./result/bin/small_variant_caller_cuda_bench 100000

bench-cnv: cnv
	@echo ">>> Executing Native CUDA Copy Number Caller Benchmark..."
	./result/bin/copy_number_caller_cuda_bench 500000

bench-str: str
	@echo ">>> Executing Native CUDA STR Caller Benchmark..."
	./result/bin/str_caller_cuda_bench 2000

bench-tfe: tfe
	@echo ">>> Executing Native CUDA Tumor Fraction Estimator Benchmark..."
	./result/bin/tumor_fraction_estimator_cuda_bench 5000

bench-pangenome: pangenome
	@echo ">>> Executing Native CUDA Pangenome Consensus Caller Benchmark..."
	./result/bin/pangenome_consensus_caller_cuda_bench 1000000

bench-fused: fused
	@echo ">>> Executing Native CUDA Fused Super-Kernel (AOT vs JIT) Benchmark..."
	./result/bin/fused_stage2_cuda_bench 1000000

bench-realdata: fused
	@echo ">>> Executing Real Roche SBX Chr20 BAM Pipeline Benchmark..."
	./result/bin/roche_sbx_cuda_eval 1000000 input/HG002.roche_sbx.chr20.bam

profile-realdata: fused
	@echo ">>> Profiling Real Roche SBX Execution with Nsight Systems (nsys)..."
	nsys profile --stats=true --force-overwrite=true -o /tmp/roche_sbx_nsys_profile ./result/bin/roche_sbx_cuda_eval 1000000 input/HG002.roche_sbx.chr20.bam

ncu-realdata: fused
	@echo ">>> Profiling Real Roche SBX Kernels with Nsight Compute (ncu)..."
	ncu --metrics sm__warps_active.avg.pct_of_peak_sustained_active,gpu__time_duration.sum,dram__throughput.avg.pct_of_peak_sustained_elapsed ./result/bin/roche_sbx_cuda_eval 100000 input/HG002.roche_sbx.chr20.bam

profile-aligner: aligner
	@echo ">>> Profiling CUDA Aligner with Nsight Systems (nsys)..."
	nsys profile --stats=true ./result/bin/aligner_cuda_bench 500000

ncu-aligner: aligner
	@echo ">>> Profiling CUDA Aligner Kernels with Nsight Compute (ncu)..."
	ncu --set full ./result/bin/aligner_cuda_bench 100000

# ----------------------------------------------------------------------------
# Layered Container Targets
# ----------------------------------------------------------------------------

docker-aligner:
	@echo ">>> Building Layered OCI Container for Aligner..."
	$(NIX_BIN) build .#docker-aligner --print-build-logs

docker-demux:
	@echo ">>> Building Layered OCI Container for Demux..."
	$(NIX_BIN) build .#docker-demux --print-build-logs

docker-read-collapser:
	@echo ">>> Building Layered OCI Container for Read Collapser..."
	$(NIX_BIN) build .#docker-read-collapser --print-build-logs

docker-metrics:
	@echo ">>> Building Layered OCI Container for Alignment Metrics..."
	$(NIX_BIN) build .#docker-metrics --print-build-logs

docker-svc:
	@echo ">>> Building Layered OCI Container for Small Variant Caller..."
	$(NIX_BIN) build .#docker-svc --print-build-logs

docker-cnv:
	@echo ">>> Building Layered OCI Container for Copy Number Caller..."
	$(NIX_BIN) build .#docker-cnv --print-build-logs

docker-str:
	@echo ">>> Building Layered OCI Container for STR Caller..."
	$(NIX_BIN) build .#docker-str --print-build-logs

docker-tfe:
	@echo ">>> Building Layered OCI Container for Tumor Fraction Estimator..."
	$(NIX_BIN) build .#docker-tfe --print-build-logs

docker-pangenome:
	@echo ">>> Building Layered OCI Container for Pangenome Consensus Caller..."
	$(NIX_BIN) build .#docker-pangenome --print-build-logs

docker-fused:
	@echo ">>> Building Layered OCI Container for Fused Super-Kernel Engine..."
	$(NIX_BIN) build .#docker-fused --print-build-logs

docker-all:
	@echo ">>> Building Unified Layered OCI Container..."
	$(NIX_BIN) build .#docker-all --print-build-logs

stream-aligner:
	@echo ">>> Streaming Aligner Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-aligner | $(PODMAN_BIN) load

stream-demux:
	@echo ">>> Streaming Demux Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-demux | $(PODMAN_BIN) load

stream-read-collapser:
	@echo ">>> Streaming Read Collapser Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-read-collapser | $(PODMAN_BIN) load

stream-metrics:
	@echo ">>> Streaming Alignment Metrics Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-metrics | $(PODMAN_BIN) load

stream-svc:
	@echo ">>> Streaming Small Variant Caller Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-svc | $(PODMAN_BIN) load

stream-cnv:
	@echo ">>> Streaming Copy Number Caller Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-cnv | $(PODMAN_BIN) load

stream-str:
	@echo ">>> Streaming STR Caller Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-str | $(PODMAN_BIN) load

stream-tfe:
	@echo ">>> Streaming Tumor Fraction Estimator Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-tfe | $(PODMAN_BIN) load

stream-pangenome:
	@echo ">>> Streaming Pangenome Consensus Caller Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-pangenome | $(PODMAN_BIN) load

stream-fused:
	@echo ">>> Streaming Fused Super-Kernel Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-fused | $(PODMAN_BIN) load

run-aligner: stream-aligner
	@echo ">>> Running Containerized CUDA Aligner with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-aligner-cuda:latest 500000

run-demux: stream-demux
	@echo ">>> Running Containerized CUDA Demux with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-demux-cuda:latest 500000

run-read-collapser: stream-read-collapser
	@echo ">>> Running Containerized CUDA Read Collapser with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-read-collapser-cuda:latest 500000

run-metrics: stream-metrics
	@echo ">>> Running Containerized CUDA Alignment Metrics with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-alignment-metrics-cuda:latest 1000000

run-svc: stream-svc
	@echo ">>> Running Containerized CUDA Small Variant Caller with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-small-variant-caller-cuda:latest 100000

run-cnv: stream-cnv
	@echo ">>> Running Containerized CUDA Copy Number Caller with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-copy-number-caller-cuda:latest 500000

run-str: stream-str
	@echo ">>> Running Containerized CUDA STR Caller with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-str-caller-cuda:latest 2000

run-tfe: stream-tfe
	@echo ">>> Running Containerized CUDA Tumor Fraction Estimator with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-tumor-fraction-estimator-cuda:latest 5000

run-pangenome: stream-pangenome
	@echo ">>> Running Containerized CUDA Pangenome Consensus Caller with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-pangenome-consensus-caller-cuda:latest 1000000

run-fused: stream-fused
	@echo ">>> Running Containerized CUDA Fused Super-Kernel with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-fused-engine-cuda:latest 1000000

# ----------------------------------------------------------------------------
# Development Utilities
# ----------------------------------------------------------------------------

dev:
	$(NIX_BIN) develop

check:
	$(NIX_BIN) flake show

clean:
	rm -rf result result-* build
