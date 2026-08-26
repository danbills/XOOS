{
  description = "XOOS: SBX Optimized Open Source secondary genomics analysis with CUDA acceleration";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          config.cudaSupport = true;
        };

        xoosDemuxPkg = pkgs.stdenv.mkDerivation {
          pname = "xoos-demux-cuda";
          version = "1.0.0";
          src = ./.;

          nativeBuildInputs = [
            pkgs.cmake
            pkgs.ninja
            pkgs.cudaPackages.cuda_nvcc
            pkgs.cudaPackages.cudatoolkit
          ];

          buildInputs = [
            pkgs.zlib
            pkgs.zstd
            pkgs.boost
            pkgs.nlohmann_json
            pkgs.fmt
            pkgs.spdlog
            pkgs.xxHash
          ];

          cmakeFlags = [
            "-DCMAKE_BUILD_TYPE=Release"
            "-DENABLE_CUDA=ON"
            "-DSTATIC_LINK_DISABLE=ON"
            "-DCODE_COVERAGE_ENABLE=OFF"
          ];
        };

        dockerImg = import ./nix/docker-image.nix {
          inherit pkgs;
          xoosDemuxPackage = xoosDemuxPkg;
        };

      in {
        packages = {
          default = xoosDemuxPkg;
          demux = xoosDemuxPkg;
          docker = dockerImg;
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
            cudaPackages.nsight_systems
            cudaPackages.nsight_compute
            zlib
            zstd
            boost
          ];

          shellHook = ''
            export CUDA_PATH=${pkgs.cudaPackages.cudatoolkit}
            export PATH=${pkgs.cudaPackages.cuda_nvcc}/bin:$PATH
            echo "🚀 XOOS CUDA DevShell Activated (Blackwell SM120 / CUDA 13.x)"
          '';
        };
      }
    );
}
