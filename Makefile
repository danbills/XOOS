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

# ----------------------------------------------------------------------------
# Benchmarking Targets
# ----------------------------------------------------------------------------

bench-aligner: aligner
	@echo ">>> Executing Native CUDA Aligner Benchmark..."
	./result/bin/aligner_cuda_bench 500000

bench-demux: demux
	@echo ">>> Executing Native CUDA Demux Benchmark..."
	./result/bin/demux_cuda_bench 500000

# ----------------------------------------------------------------------------
# Layered Container Targets
# ----------------------------------------------------------------------------

docker-aligner:
	@echo ">>> Building Layered OCI Container for Aligner..."
	$(NIX_BIN) build .#docker-aligner --print-build-logs

docker-demux:
	@echo ">>> Building Layered OCI Container for Demux..."
	$(NIX_BIN) build .#docker-demux --print-build-logs

docker-all:
	@echo ">>> Building Unified Layered OCI Container..."
	$(NIX_BIN) build .#docker-all --print-build-logs

stream-aligner:
	@echo ">>> Streaming Aligner Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-aligner | $(PODMAN_BIN) load

stream-demux:
	@echo ">>> Streaming Demux Container directly into Podman..."
	$(NIX_BIN) run .#stream-docker-demux | $(PODMAN_BIN) load

run-aligner: stream-aligner
	@echo ">>> Running Containerized CUDA Aligner with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-aligner-cuda:latest 500000

run-demux: stream-demux
	@echo ">>> Running Containerized CUDA Demux with GPU Passthrough..."
	$(PODMAN_BIN) run --rm --device nvidia.com/gpu=all localhost/roche-axelios/xoos-demux-cuda:latest 500000

# ----------------------------------------------------------------------------
# Development Utilities
# ----------------------------------------------------------------------------

dev:
	$(NIX_BIN) develop

check:
	$(NIX_BIN) flake show

clean:
	rm -rf result result-* build
