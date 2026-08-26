# ============================================================================
# Nix Declarative Layered Container Builder for XOOS CUDA Accelerated Engines
# File: nix/docker-image.nix
#
# Generates multi-layered OCI/Docker container images optimized for Nextflow / xoosnf.
#
# Layering Structure (Immutable Cache Optimization):
#   - Layer 1: Core POSIX utilities (bash, coreutils, procps, glibc)
#   - Layer 2: CUDA 13.x driver/runtime libraries (cudatoolkit, libcudart)
#   - Layer 3: System compression & math libraries (zlib, zstd, boost, xxhash)
#   - Layer 4 (Top Layer): Compiled XOOS binaries (~5 MB) & adapter designs
#
# When C++/CUDA code changes, only Layer 4 is rebuilt and transmitted over network!
# ============================================================================

{ pkgs ? import <nixpkgs> { config.allowUnfree = true; config.cudaSupport = true; }
, package
, imageName ? "roche-axelios/xoos-cuda"
, imageTag ? "latest"
, entrypointBinary ? "/bin/demux"
}:

let
  layeredImage = pkgs.dockerTools.buildLayeredImage {
    name = imageName;
    tag = imageTag;
    maxLayers = 16;

    contents = [
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
    ];

    config = {
      Cmd = [ entrypointBinary "--help" ];
      Entrypoint = [ entrypointBinary ];
      Env = [
        "PATH=/bin:/usr/bin:${package}/bin"
        "NVIDIA_VISIBLE_DEVICES=all"
        "NVIDIA_DRIVER_CAPABILITIES=compute,utility"
        "LC_ALL=C.UTF-8"
      ];
      WorkingDir = "/data";
      Volumes = {
        "/data" = {};
      };
    };
  };

  streamedImage = pkgs.dockerTools.streamLayeredImage {
    name = imageName;
    tag = imageTag;
    maxLayers = 16;

    contents = [
      pkgs.bashInteractive
      pkgs.coreutils
      pkgs.procps
      pkgs.cudaPackages.cudatoolkit
      pkgs.cudaPackages.cuda_cudart
      pkgs.zlib
      pkgs.zstd
      pkgs.xxHash
      package
    ];

    config = {
      Cmd = [ entrypointBinary "--help" ];
      Entrypoint = [ entrypointBinary ];
      Env = [
        "PATH=/bin:/usr/bin:${package}/bin"
        "NVIDIA_VISIBLE_DEVICES=all"
        "NVIDIA_DRIVER_CAPABILITIES=compute,utility"
        "LC_ALL=C.UTF-8"
      ];
      WorkingDir = "/data";
      Volumes = {
        "/data" = {};
      };
    };
  };

in {
  layered = layeredImage;
  streamed = streamedImage;
}
