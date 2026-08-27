# ============================================================================
# Nix Declarative Layered Container Builder for XOOS CUDA Accelerated Engines
# File: nix/docker-image.nix
#
# Generates multi-layered OCI/Docker container images optimized for Nextflow / xoosnf.
#
# Layering Structure (Immutable Cache Optimization):
#   - Layer 1: Core POSIX utilities (bash, coreutils, procps, glibc)
#   - Layer 2: CUDA 13.x driver/runtime libraries (cudatoolkit, libcudart)
#   - Layer 3: NVIDIA Profiling Suite (nsys / Nsight Systems, ncu / Nsight Compute)
#   - Layer 4: System compression & math libraries (zlib, zstd, boost, xxhash)
#   - Layer 5 (Top Layer): Compiled XOOS binaries (~5 MB) & adapter designs
#
# When C++/CUDA code changes, only Layer 5 is rebuilt and transmitted over network!
# ============================================================================

{ pkgs ? import <nixpkgs> { config.allowUnfree = true; config.cudaSupport = true; }
, package
, imageName ? "roche-axelios/xoos-cuda"
, imageTag ? "latest"
, entrypointBinary ? "/bin/demux"
, enableProfiling ? true
}:

let
  profilingPackages = if enableProfiling then [
    pkgs.cudaPackages.nsight_systems
    pkgs.cudaPackages.nsight_compute
  ] else [];

  profilingPaths = if enableProfiling then
    ":${pkgs.cudaPackages.nsight_systems}/bin:${pkgs.cudaPackages.nsight_compute}/bin"
  else "";

  commonContents = [
    # Base Linux & Process Inspection
    pkgs.bashInteractive
    pkgs.coreutils
    pkgs.procps

    # CUDA 13.x Runtime & Acceleration
    pkgs.cudaPackages.cudatoolkit
    pkgs.cudaPackages.cuda_cudart

    # High-performance compression & hashing
    pkgs.zlib
    pkgs.zstd
    pkgs.xxhash

    # Target XOOS Package
    package
  ] ++ profilingPackages;

  commonConfig = {
    Cmd = [ entrypointBinary "--help" ];
    Entrypoint = [ entrypointBinary ];
    Env = [
      "PATH=/bin:/usr/bin:${package}/bin${profilingPaths}"
      "NVIDIA_VISIBLE_DEVICES=all"
      "NVIDIA_DRIVER_CAPABILITIES=compute,utility"
      "LC_ALL=C.UTF-8"
    ];
    WorkingDir = "/data";
    Volumes = {
      "/data" = {};
    };
  };

  layeredImage = pkgs.dockerTools.buildLayeredImage {
    name = imageName;
    tag = imageTag;
    maxLayers = 16;
    contents = commonContents;
    config = commonConfig;
  };

  streamedImage = pkgs.dockerTools.streamLayeredImage {
    name = imageName;
    tag = imageTag;
    maxLayers = 16;
    contents = commonContents;
    config = commonConfig;
  };

in {
  layered = layeredImage;
  streamed = streamedImage;
}
