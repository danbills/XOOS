# ============================================================================
# Nix Declarative Layered Docker Image Builder for XOOS CUDA Demux
# File: nix/docker-image.nix
#
# Generates a multi-layered OCI/Docker container image optimized for Nextflow / xoosnf.
# Dependencies (glibc, CUDA runtime, zlib) are cached in lower immutable layers,
# while the top application layer (~5 MB) rebuilds instantly on code updates.
# ============================================================================

{ pkgs ? import <nixpkgs> { config.allowUnfree = true; config.cudaSupport = true; }
, xoosDemuxPackage
}:

pkgs.dockerTools.buildLayeredImage {
  name = "roche-axelios/xoos-demux-cuda";
  tag = "latest";

  # Maximize layer sharing across container pulls
  maxLayers = 8;

  contents = [
    # System base utilities
    pkgs.bashInteractive
    pkgs.coreutils
    pkgs.procps

    # CUDA and compression runtimes
    pkgs.cudaPackages.cudatoolkit
    pkgs.zlib
    pkgs.zstd

    # The compiled XOOS binary package
    xoosDemuxPackage
  ];

  config = {
    Cmd = [ "/bin/demux" "--help" ];
    Entrypoint = [ "/bin/demux" ];
    Env = [
      "PATH=/bin:/usr/bin:${xoosDemuxPackage}/bin"
      "NVIDIA_VISIBLE_DEVICES=all"
      "NVIDIA_DRIVER_CAPABILITIES=compute,utility"
      "LC_ALL=C.UTF-8"
    ];
    WorkingDir = "/data";
    Volumes = {
      "/data" = {};
    };
  };
}
