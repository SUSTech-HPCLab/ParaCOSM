#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"

# ── 1. Setup Intel oneAPI environment ──
ONEAPI_ROOT="${ONEAPI_ROOT:-$HOME/intel/oneapi}"
SETVARS="${ONEAPI_ROOT}/setvars.sh"

if [[ -f "${SETVARS}" ]]; then
    echo "[build] Sourcing oneAPI environment: ${SETVARS}"
    set +u  # oneAPI setvars.sh may reference unbound variables
    source "${SETVARS}" --force 2>&1 | tail -1
    set -u
else
    echo "[build] WARNING: ${SETVARS} not found, assuming environment is already set up"
fi

# ── 2. Ensure external/taskflow exists ──
TASKFLOW_DIR="${SCRIPT_DIR}/external/taskflow"
if [[ ! -d "${TASKFLOW_DIR}" ]]; then
    echo "[build] taskflow not found at ${TASKFLOW_DIR}"
    echo "[build] Attempting to clone taskflow from GitHub..."
    mkdir -p "${SCRIPT_DIR}/external"
    git clone --depth 1 https://github.com/taskflow/taskflow.git /tmp/taskflow_clone \
        && cp -r /tmp/taskflow_clone/taskflow "${TASKFLOW_DIR}" \
        && rm -rf /tmp/taskflow_clone \
        && echo "[build] taskflow downloaded successfully" \
        || { echo "[build] ERROR: Failed to download taskflow. Please place it manually at ${TASKFLOW_DIR}"; exit 1; }
fi

# ── 3. Configure & Build ──
JOBS="${JOBS:-$(nproc)}"
BUILD_TYPE="${BUILD_TYPE:-Release}"

echo "[build] Build directory : ${BUILD_DIR}"
echo "[build] Build type      : ${BUILD_TYPE}"
echo "[build] Parallel jobs   : ${JOBS}"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake "${SCRIPT_DIR}" \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    2>&1

cmake --build . -j"${JOBS}" 2>&1

echo ""
echo "[build] ✅ Build complete. Binaries:"
ls -lh "${BUILD_DIR}/bin/" 2>/dev/null || echo "  (check ${BUILD_DIR} for outputs)"
