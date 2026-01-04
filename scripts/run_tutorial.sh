#!/bin/bash
set -e

# Get the project root directory
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
TUTORIAL_RUN_DIR="${PROJECT_ROOT}/_tutorial_run"

echo "Project Root: ${PROJECT_ROOT}"
echo "Tutorial Run Dir: ${TUTORIAL_RUN_DIR}"

# Clean up previous run
if [ -d "${TUTORIAL_RUN_DIR}" ]; then
    rm -rf "${TUTORIAL_RUN_DIR}"
fi
mkdir -p "${TUTORIAL_RUN_DIR}"

# Unpack Tutorial data
echo "Unpacking Tutorial data..."
tar -xzf "${PROJECT_ROOT}/Tutorial.tar.gz" -C "${TUTORIAL_RUN_DIR}"

# The tarball contains a "Tutorial" directory.
WORK_DIR="${TUTORIAL_RUN_DIR}/Tutorial"
cd "${WORK_DIR}"

echo "Running in: $(pwd)"

# Run Clustering via CLI
echo "----------------------------------------------------------------"
echo "Running K-Centers Clustering (via hkdm)..."
# Using absolute path to hkdm if it's not in PATH, or assuming conda env is active
# Since we are running in a shell that might not have the env activated if called directly,
# but usually we run this *inside* the env.
hkdm cluster kcenters \
    --trajlist trajlist \
    --atomlist atom_indices \
    --topology native.pdb \
    --iext xtc \
    --n-clusters 100 \
    --output-dir .

# Run APLoD Clustering via CLI (Demonstration)
echo "----------------------------------------------------------------"
echo "Running APLoD Clustering (via hkdm)..."
hkdm cluster aplod \
    --trajlist trajlist \
    --atomlist atom_indices \
    --topology native.pdb \
    --iext xtc \
    --rho-cutoff 0.1 \
    --delta-cutoff 0.1 \
    --n-neighbors 10 \
    --output-dir .

# Run Lumping via CLI
echo "----------------------------------------------------------------"
echo "Running PCCA Lumping (via hkdm)..."
# Find assignment file from k-centers
ASSIGNMENT_FILE=$(ls assignments_kcenters_n_*.txt | head -n 1)
echo "Using assignment file: ${ASSIGNMENT_FILE}"

hkdm lump pcca \
    --assignments "${ASSIGNMENT_FILE}" \
    --n-macro 4 \
    --homedir .

echo "----------------------------------------------------------------"
echo "Tutorial run completed successfully."
echo "Outputs are in: ${WORK_DIR}"
ls -l "${WORK_DIR}"