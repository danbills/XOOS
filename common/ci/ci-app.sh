#!/bin/bash

set -eu -o pipefail

code_coverage_enable="${code_coverage_enable:-OFF}"
build_type="${build_type:-Release}"
version="${version:-unknown}"
build_testing="${build_testing:-OFF}"

cmake_opts=()
use_sccache=false
if command -v sccache >/dev/null 2>&1 ; then
  # Verify sccache can start before using it as a compiler launcher.
  # If the cache backend (e.g. S3) is unreachable, sccache fails fatally
  # and would break the build.
  if sccache --start-server 2>/dev/null ; then
    cmake_opts+=("-DCMAKE_C_COMPILER_LAUNCHER=sccache")
    cmake_opts+=("-DCMAKE_CXX_COMPILER_LAUNCHER=sccache")
    use_sccache=true
  else
    echo "WARNING: sccache server failed to start, building without cache" >&2
    unset SCCACHE_BUCKET SCCACHE_REGION AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN
  fi
fi

cmake --preset conan-default \
  -DLINT_ENABLE=OFF \
  -DBUILD_TESTING="${build_testing}" \
  -DCODE_COVERAGE_ENABLE="${code_coverage_enable}" \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DVERSION="${version}" \
  "${cmake_opts[@]}"

cmake --build --preset "conan-${build_type,,}"

if [ "${use_sccache}" = "true" ] ; then
  sccache --show-stats
fi
