{
  description = "XOOS: High-Performance CUDA-Accelerated Secondary Genomics Pipeline for NVIDIA Blackwell (SM120) & CUDA 13.x";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachSystem [ "x86_64-linux" "aarch64-linux" ] (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          config.cudaSupport = true;
        };

        # Base derivation for XOOS CUDA components
        mkXoosCudaPackage = { pname, subDir, cmakeExtraFlags ? [] }: pkgs.stdenv.mkDerivation {
          inherit pname;
          version = "1.0.0";
          src = ./.;

          nativeBuildInputs = [
            pkgs.cmake
            pkgs.ninja
            pkgs.cudaPackages.cuda_nvcc
            pkgs.cudaPackages.cudatoolkit
          ];

          buildInputs = [
            pkgs.cudaPackages.cuda_cudart
            pkgs.zlib
            pkgs.zstd
            pkgs.boost
            pkgs.nlohmann_json
            pkgs.fmt
            pkgs.spdlog
            pkgs.xxhash
            pkgs.cli11
            pkgs.magic-enum
            pkgs.microsoft-gsl
            pkgs.taskflow
          ];

          cmakeFlags = [
            "-S" "../${subDir}"
            "-DCMAKE_BUILD_TYPE=Release"
            "-DENABLE_CUDA=ON"
            "-DLINT_ENABLE=OFF"
            "-DCMAKE_CUDA_ARCHITECTURES=120;90;89;86"
            "-DSTATIC_LINK_DISABLE=ON"
            "-DCODE_COVERAGE_ENABLE=OFF"
          ] ++ cmakeExtraFlags;

          installPhase = ''
            mkdir -p $out/bin
            find . -type f -executable -exec cp {} $out/bin/ \; 2>/dev/null || true
          '';
        };

        xoosDemuxPkg = mkXoosCudaPackage {
          pname = "xoos-demux-cuda";
          subDir = "demux";
        };

        xoosAlignerPkg = mkXoosCudaPackage {
          pname = "xoos-aligner-cuda";
          subDir = "aligner";
        };

        xoosReadCollapserPkg = mkXoosCudaPackage {
          pname = "xoos-read-collapser-cuda";
          subDir = "read_collapser";
        };

        # Layered Docker Images
        demuxDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosDemuxPkg;
          imageName = "roche-axelios/xoos-demux-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/demux";
        };

        alignerDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosAlignerPkg;
          imageName = "roche-axelios/xoos-aligner-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/aligner_cuda_bench";
        };

        readCollapserDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosReadCollapserPkg;
          imageName = "roche-axelios/xoos-read-collapser-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/read_collapser_cuda_bench";
        };

        unifiedDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosDemuxPkg;
          imageName = "roche-axelios/xoos-cuda-all";
          imageTag = "latest";
          entrypointBinary = "/bin/demux";
        };

      in {
        packages = {
          default = xoosDemuxPkg;
          demux = xoosDemuxPkg;
          aligner = xoosAlignerPkg;
          read-collapser = xoosReadCollapserPkg;
          
          # Container outputs
          docker-demux = demuxDocker.layered;
          docker-aligner = alignerDocker.layered;
          docker-read-collapser = readCollapserDocker.layered;
          docker-all = unifiedDocker.layered;
          
          stream-docker-demux = demuxDocker.streamed;
          stream-docker-aligner = alignerDocker.streamed;
          stream-docker-read-collapser = readCollapserDocker.streamed;
          stream-docker-all = unifiedDocker.streamed;
        };

        devShells.default = pkgs.mkShell {
          name = "xoos-cuda-dev-shell";
          packages = with pkgs; [
            cmake
            ninja
            gnumake
            gcc
            gdb
            cudaPackages.cudatoolkit
            cudaPackages.cuda_nvcc
            cudaPackages.cuda_cudart
            cudaPackages.nsight_systems
            cudaPackages.nsight_compute
            zlib
            zstd
            boost
            xxHash
          ];

          shellHook = ''
            export CUDA_PATH=${pkgs.cudaPackages.cudatoolkit}
            export PATH=${pkgs.cudaPackages.cuda_nvcc}/bin:$PATH
            echo "🚀 XOOS CUDA 13.x DevShell Activated (Blackwell SM120 Target)"
          '';
        };
      }
    );
}
