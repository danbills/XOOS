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

        xoosMetricsPkg = mkXoosCudaPackage {
          pname = "xoos-alignment-metrics-cuda";
          subDir = "alignment_metrics";
        };

        xoosSvcPkg = mkXoosCudaPackage {
          pname = "xoos-small-variant-caller-cuda";
          subDir = "small_variant_caller";
        };

        xoosCnvPkg = mkXoosCudaPackage {
          pname = "xoos-copy-number-caller-cuda";
          subDir = "copy_number_caller";
        };

        xoosStrPkg = mkXoosCudaPackage {
          pname = "xoos-str-caller-cuda";
          subDir = "str_caller";
        };

        xoosTfePkg = mkXoosCudaPackage {
          pname = "xoos-tumor-fraction-estimator-cuda";
          subDir = "tumor_fraction_estimator";
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

        metricsDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosMetricsPkg;
          imageName = "roche-axelios/xoos-alignment-metrics-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/alignment_metrics_cuda_bench";
        };

        svcDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosSvcPkg;
          imageName = "roche-axelios/xoos-small-variant-caller-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/small_variant_caller_cuda_bench";
        };

        cnvDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosCnvPkg;
          imageName = "roche-axelios/xoos-copy-number-caller-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/copy_number_caller_cuda_bench";
        };

        strDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosStrPkg;
          imageName = "roche-axelios/xoos-str-caller-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/str_caller_cuda_bench";
        };

        tfeDocker = import ./nix/docker-image.nix {
          inherit pkgs;
          package = xoosTfePkg;
          imageName = "roche-axelios/xoos-tumor-fraction-estimator-cuda";
          imageTag = "latest";
          entrypointBinary = "/bin/tumor_fraction_estimator_cuda_bench";
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
          metrics = xoosMetricsPkg;
          alignment-metrics = xoosMetricsPkg;
          svc = xoosSvcPkg;
          small-variant-caller = xoosSvcPkg;
          cnv = xoosCnvPkg;
          copy-number-caller = xoosCnvPkg;
          str = xoosStrPkg;
          str-caller = xoosStrPkg;
          tfe = xoosTfePkg;
          tumor-fraction = xoosTfePkg;
          tumor-fraction-estimator = xoosTfePkg;
          
          # Container outputs
          docker-demux = demuxDocker.layered;
          docker-aligner = alignerDocker.layered;
          docker-read-collapser = readCollapserDocker.layered;
          docker-metrics = metricsDocker.layered;
          docker-svc = svcDocker.layered;
          docker-cnv = cnvDocker.layered;
          docker-str = strDocker.layered;
          docker-tfe = tfeDocker.layered;
          docker-all = unifiedDocker.layered;
          
          stream-docker-demux = demuxDocker.streamed;
          stream-docker-aligner = alignerDocker.streamed;
          stream-docker-read-collapser = readCollapserDocker.streamed;
          stream-docker-metrics = metricsDocker.streamed;
          stream-docker-svc = svcDocker.streamed;
          stream-docker-cnv = cnvDocker.streamed;
          stream-docker-str = strDocker.streamed;
          stream-docker-tfe = tfeDocker.streamed;
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
