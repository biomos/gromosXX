#!/usr/bin/env bash
#
# Install a cached md++ build into an asv environment.
#
# asv's default install step is "pip install <wheel>", which has nothing to do
# with a C++ project. Copying the install prefix into the environment root
# instead puts the binary at <env>/bin/md, so benchmarks can find it via
# sys.prefix without any path configuration.
#
# Usage: asv_install.sh <build_cache_dir> <env_dir>

set -euo pipefail

CACHE_DIR="${1:?build_cache_dir required}"
ENV_DIR="${2:?env_dir required}"
PREFIX="${CACHE_DIR}/prefix"

if [ ! -d "${PREFIX}" ]; then
    echo "no build prefix at ${PREFIX} -- the build step did not complete" >&2
    exit 1
fi

mkdir -p "${ENV_DIR}"
cp -a "${PREFIX}/." "${ENV_DIR}/"

echo "installed: $(ls "${ENV_DIR}/bin" | tr '\n' ' ')"
