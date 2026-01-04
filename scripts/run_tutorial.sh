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

# Run Clustering
echo "----------------------------------------------------------------"
echo "Running K-Centers Clustering..."
python "${PROJECT_ROOT}/scripts/test_kcenters.py" \
    -t trajlist \
    -a atom_indices \
    -g native.pdb \
    -e xtc \
    -n 100 \
    -o .

# The output of clustering is assignments_kcenters_n_?.txt
# We need to find the exact name. It depends on n_microstates found.
# test_kcenters.py prints "Estimated number of clusters: ..."
# And saves "assignments_kcenters_n_X.txt"

# Let's find the assignment file
ASSIGNMENT_FILE=$(ls assignments_kcenters_n_*.txt | head -n 1)
echo "Found assignment file: ${ASSIGNMENT_FILE}"

# Run Lumping
echo "----------------------------------------------------------------"
echo "Running PCCA Lumping..."
python "${PROJECT_ROOT}/scripts/doLumping.py" \
    -c "${ASSIGNMENT_FILE}" \
    -m 4 \
    -l traj_len.txt

echo "----------------------------------------------------------------"
echo "Tutorial run completed successfully."
echo "Outputs are in: ${WORK_DIR}"
ls -l "${WORK_DIR}"
