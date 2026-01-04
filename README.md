# HKDataMiner

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)

**HKDataMiner** is a comprehensive Python library designed for the statistical analysis of biomolecular dynamics simulations. It provides a suite of algorithms for clustering conformational states and constructing Markov State Models (MSMs) to elucidate the kinetic properties of molecular systems.

Developed by **Prof. Xuhui Huang's Group** at HKUST.

## Key Features

HKDataMiner integrates advanced data mining techniques tailored for Molecular Dynamics (MD) trajectories:

### 1. Clustering Algorithms (Microstates)
*   **K-Centers**: A geometric clustering algorithm that minimizes the maximum radius of clusters. Effective for ensuring full coverage of the conformational space.
*   **APLoD (Adaptive Partitioning by Local Density-peaks)**: An efficient, density-based clustering algorithm. APLoD estimates the local density of MD conformations using k-nearest neighbors and groups high-density regions into clusters. It features **adaptive resolution**: generating larger clusters in low-density regions and finer clusters in high-density regions, significantly reducing computational cost and memory usage compared to traditional methods.

### 2. Lumping Algorithms (Macrostates)
*   **PCCA (Perron Cluster Cluster Analysis)**: Uses the eigenspectrum of the transition probability matrix to aggregate microstates into kinetically metastable macrostates.
*   **PCCA+**: An improved version of PCCA with better fuzzy membership handling.
*   **Spectral Clustering**: Graph-based clustering on the transition matrix.
*   **Ward**: Hierarchical clustering minimizing variance.

### 3. Template Matching
A specialized module for analyzing cryo-EM data distributions using template matching techniques (requires Xmipp).

## Installation

We recommend using `mamba` (or `conda`) to manage the scientific Python stack.

### Prerequisites
*   Linux or macOS
*   Python 3.10+
*   Mamba/Conda

### Steps

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/liusong299/HK_DataMiner.git
    cd HK_DataMiner
    ```

2.  **Create and activate the environment:**
    ```bash
    mamba env create -f environment.yml
    mamba activate hkdataminer-py310
    ```

3.  **Install the package:**
    ```bash
    pip install -e .
    ```

## Quick Start: Alanine Dipeptide Tutorial

We provide a one-line command to run a full end-to-end analysis on the Alanine Dipeptide dataset. This workflow includes:
1.  Unpacking tutorial data.
2.  **Clustering** trajectory frames into microstates (using K-Centers or APLoD).
3.  **Lumping** microstates into macrostates using PCCA.
4.  Visualizing the results.

**Run the tutorial:**
```bash
hkdm tutorial
```

**Output:**
Results are saved in `_tutorial_run/Tutorial/`. Key files include:
*   `kcenters_n_100.png`: Microstate distribution on the Ramachandran plot.
*   `PCCA_100_to_4_states_Matrix.png`: Transition probability matrix.
*   `assignments_*.txt`: Frame-to-state assignments.

## Usage Guide

HKDataMiner provides a unified CLI tool `hkdm`.

### 1. Clustering with K-Centers
```bash
hkdm cluster kcenters \
    --trajlist trajlist.txt \
    --atomlist atom_indices.txt \
    --topology native.pdb \
    --n-clusters 100 \
    --output-dir ./results
```

### 2. Clustering with APLoD
APLoD is recommended for large datasets.
```bash
hkdm cluster aplod \
    --trajlist trajlist.txt \
    --atomlist atom_indices.txt \
    --topology native.pdb \
    --rho-cutoff 1.0 \
    --delta-cutoff 1.0 \
    --n-neighbors 100 \
    --output-dir ./results
```

### 3. Lumping (MSM Construction)
```bash
hkdm lump pcca \
    --assignments ./results/assignments_kcenters_n_100.txt \
    --traj-len ./results/traj_len.txt \
    --n-macro 4 \
    --homedir ./results
```

## References

If you use **APLoD** in your research, please cite:

> Song Liu, Lizhe Zhu and Xuhui Huang, "Adaptive Partitioning by Local Density-peaks (APLoD): An efficient density-based clustering algorithm for analyzing molecular dynamics trajectories". *Journal of Computational Chemistry*, 2016.

For **HKDataMiner** general usage:

> HKUST Huang Group (http://compbio.ust.hk)

## License

This project is licensed under the **Apache 2.0 License**. See the [LICENSE](LICENSE) file for details.

## Authors

*   **Prof. Xuhui Huang** - *Project Leader*
*   **Dr. Song Liu** - *Developer*
*   **Mr. Hanlin Gu** - *Developer*
