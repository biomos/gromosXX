#!/usr/bin/env bash
#
# Build md++ for one asv environment and stage the result in asv's build cache.
#
# asv runs build commands without a shell -- it shlex-splits the string and
# only substitutes {placeholder} variables -- so "$VAR" in asv.conf.json would
# never expand. Configuration therefore arrives two ways: positional arguments
# for asv's own paths, and environment variables for the matrix settings, which
# asv passes through to build commands.
#
# This script is invoked as {conf_dir}/build/build_gromos.sh, i.e. from the
# working checkout, NOT from the commit being built. That is deliberate: the
# commits in the timeline predate this benchmark suite and do not contain it.
#
# Usage: asv_build.sh <build_dir> <build_cache_dir>
#
# Environment:
#   GROMOS_VARIANT      serial | omp | cuda | mpi (selects the -D flags)
#   GROMOS_CMAKE_ARGS   optional extra -D flags, appended last
#   GROMOS_BUILD_JOBS   parallel build jobs (default: nproc)
#   GROMOS_NVCC         path to nvcc, for the CUDA variant
#   GROMOS_CUDA_HOST_COMPILER  host compiler nvcc drives (default: g++ on PATH)

set -euo pipefail

BUILD_DIR="${1:?build_dir required}"
CACHE_DIR="${2:?build_cache_dir required}"

JOBS="${GROMOS_BUILD_JOBS:-$(nproc)}"
PREFIX="${CACHE_DIR}/prefix"

# md++ lives in a subdirectory of the repository. Commits older than the
# 2023-10-21 rename (284700de) have it at gromosXX/ instead and are outside the
# supported range; say so clearly rather than failing inside cmake.
if [ -f "${BUILD_DIR}/md++/CMakeLists.txt" ]; then
    SRC="${BUILD_DIR}/md++"
elif [ -f "${BUILD_DIR}/gromosXX/CMakeLists.txt" ]; then
    echo "This commit predates the gromosXX -> md++ rename (284700de," \
         "2023-10-21). The benchmark suite does not support it." >&2
    exit 1
else
    echo "No CMakeLists.txt under ${BUILD_DIR}/md++ -- commit too old to build" >&2
    exit 1
fi

CMAKE_ARGS=(
    -S "${SRC}"
    -B "${BUILD_DIR}/_asvbuild"
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX="${PREFIX}"
)

# ccache makes neighbouring commits cheap to build, which matters because the
# build dominates: ~230 translation units versus about a minute of benchmarking.
if command -v ccache >/dev/null 2>&1; then
    CMAKE_ARGS+=(-DCMAKE_CXX_COMPILER_LAUNCHER=ccache)
fi
if command -v ninja >/dev/null 2>&1; then
    CMAKE_ARGS+=(-GNinja)
fi
if [ -n "${GROMOS_NVCC:-}" ]; then
    if [ ! -x "${GROMOS_NVCC}" ]; then
        echo "GROMOS_NVCC=${GROMOS_NVCC} is not executable" >&2
        exit 1
    fi
    CMAKE_ARGS+=(-DCMAKE_CUDA_COMPILER="${GROMOS_NVCC}")

    # Pin the host compiler nvcc drives. When nvcc comes from a conda
    # environment (the usual case here, since there is no system CUDA
    # toolkit), it otherwise defaults to that environment's toolchain and
    # picks up its sysroot headers. Those disagree with the host glibc and
    # configuration dies with hundreds of "identifier _Float32 is undefined"
    # errors out of bits/mathcalls.h -- an error that points at the standard
    # library and says nothing about the actual cause.
    CUDA_HOST_CXX="${GROMOS_CUDA_HOST_COMPILER:-$(command -v g++ || true)}"
    if [ -n "${CUDA_HOST_CXX}" ]; then
        CMAKE_ARGS+=(-DCMAKE_CUDA_HOST_COMPILER="${CUDA_HOST_CXX}")
    fi
fi

# The variant is the single source of truth for how a build is configured;
# deriving the flags here keeps the asv matrix down to one dimension and keeps
# the mutually exclusive combinations from being expressible at all.
# Constraints are enforced by md++ itself in cmake/options.cmake: OMP and MPI
# may not both be on, and CUKERNEL requires OMP.
VARIANT="${GROMOS_VARIANT:-omp}"
case "${VARIANT}" in
    serial) VARIANT_ARGS=() ;;
    omp)    VARIANT_ARGS=(-DOMP=ON) ;;
    cuda)   VARIANT_ARGS=(-DOMP=ON -DCUKERNEL=ON) ;;
    mpi)    VARIANT_ARGS=(-DMPI=ON) ;;
    *)      echo "unknown GROMOS_VARIANT '${VARIANT}'" >&2; exit 1 ;;
esac

if [ "${VARIANT}" = "cuda" ] && [ -z "${GROMOS_NVCC:-}" ] \
   && ! command -v nvcc >/dev/null 2>&1; then
    echo "variant 'cuda' needs nvcc: set GROMOS_NVCC to its path" >&2
    exit 1
fi

# Optional extra flags, appended last so they win. Unquoted on purpose:
# this carries several independent -D flags.
# shellcheck disable=SC2206
EXTRA=(${GROMOS_CMAKE_ARGS:-})

rm -rf "${PREFIX}"
mkdir -p "${PREFIX}"

echo "=== configure variant=${VARIANT} ${VARIANT_ARGS[*]:-} ${GROMOS_CMAKE_ARGS:-}"
# No ":-" defaults here: on an empty array (the serial variant has no extra
# flags) "${arr[@]:-}" expands to one empty word rather than to nothing, and
# cmake reads that empty argument as a source directory.
cmake "${CMAKE_ARGS[@]}" ${VARIANT_ARGS[@]+"${VARIANT_ARGS[@]}"} ${EXTRA[@]+"${EXTRA[@]}"}

echo "=== build (-j ${JOBS})"
cmake --build "${BUILD_DIR}/_asvbuild" -j "${JOBS}"

echo "=== install -> ${PREFIX}"
cmake --install "${BUILD_DIR}/_asvbuild"

ls "${PREFIX}/bin"
