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

# ----------------------------------------------------------------------------
# Development Utilities
# ----------------------------------------------------------------------------

dev:
	$(NIX_BIN) develop

check:
	$(NIX_BIN) flake show

clean:
	rm -rf result result-* build
