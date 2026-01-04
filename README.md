# HKDataMiner

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)

**HKDataMiner** is a python library for constructing statistical models for biomolecular dynamics data. It is developed by **Prof. Xuhui Huang's Group** at HKUST.

## Project History & Modernization

**HKDataMiner** was originally initiated in **2014** and built upon the Python 2 ecosystem. After serving the community for several years, the project entered a dormant phase without active maintenance.

**Recently, we have revitalized the project with a major modernization:**
*   **Python 3.10+ Support**: The codebase has been fully ported from Python 2 to Python 3.10, ensuring compatibility with modern scientific computing stacks.
*   **Modern Engineering**: We reorganized the package structure, added `mamba` environment support, and introduced a unified CLI (`hkdm`) for better usability.
*   **Active Maintenance**: We are committed to further maintaining and updating HKDataMiner to support modern computational workflows and new research needs.

The core highlight of this package is **APLoD** (Adaptive Partitioning by Local Density-peaks), a highly efficient clustering algorithm designed specifically for analyzing large-scale Molecular Dynamics (MD) trajectories and constructing Markov State Models (MSMs).

## Why APLoD?

**Adaptive Partitioning by Local Density-peaks (APLoD)** addresses the challenge of clustering ultra-large MD datasets containing millions of conformations.

### Key Advantages:
1.  **Extreme Efficiency**: APLoD reduces running time and memory usage by **2–3 orders of magnitude** compared to standard Density Peaks (DP) algorithms. It achieves a temporal complexity of $O(N \log N)$ and spatial complexity of $O(N)$, making it feasible to run on standard desktops even for massive datasets.
2.  **Adaptive Resolution**: Unlike geometric clustering (e.g., K-Centers) which tends to partition space uniformly, APLoD is density-based. It produces **clusters with adaptive sizes**:
    *   **High-density regions** (energy minima): Finer resolution with smaller clusters.
    *   **Low-density regions** (transition states): Coarser resolution with larger clusters.
    *   This automatically minimizes statistical error within clusters, preserving kinetic boundaries.
3.  **Local Density Estimation**: By utilizing **k-Nearest-Neighbors (kNN)** search (via VP-trees), APLoD estimates density locally. This avoids the quadratic $O(N^2)$ cost of global density estimation found in traditional methods.

## Features

### 1. Clustering Algorithms (Microstate Identification)
*   **APLoD (Recommended)**: Best for large-scale datasets. Groups conformations based on local density peaks.
*   **K-Centers**: Standard geometric clustering minimizing the maximum cluster radius. Useful for ensuring uniform coverage of conformational space.

### 2. Lumping Algorithms (Macrostate Construction)
Once microstates are identified, HKDataMiner provides algorithms to lump them into kinetically metastable macrostates:
*   **PCCA (Perron Cluster Cluster Analysis)**: Spectral method based on the transition probability matrix eigenvalues.
*   **PCCA+**: Robust extension of PCCA with fuzzy memberships.
*   **Spectral Clustering**: Graph-based clustering.
*   **Ward**: Hierarchical clustering.

### 3. Template Matching
A specialized module for cryo-EM data distribution analysis (requires Xmipp).

## Installation

We recommend using `mamba` (or `conda`) to manage the environment.

```bash
# 1. Clone the repository
git clone https://github.com/liusong299/HK_DataMiner.git
cd HK_DataMiner

# 2. Create the environment (Python 3.10)
mamba env create -f environment.yml
mamba activate hkdataminer-py310

# 3. Install in editable mode
pip install -e .
```

## Quick Start

We provide a one-line command to run a full tutorial on the **Alanine Dipeptide** dataset. This workflow unpacks data, performs clustering, builds an MSM via lumping, and plots the results.

```bash
hkdm tutorial
```

**Results** will be saved in `_tutorial_run/Tutorial/`, including:
*   `assignments_*.txt`: Microstate assignments for every frame.
*   `cluster_centers.pdb`: Structures of cluster centers.
*   `PCCA_*_Matrix.png`: Transition probability matrix of the macrostates.
*   `*_Metastability_Modularity.txt`: MSM quality metrics.

## CLI Usage

HKDataMiner provides a unified command-line interface `hkdm`.

### Clustering with APLoD (Recommended)
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

### Clustering with K-Centers
```bash
hkdm cluster kcenters \
    --trajlist trajlist.txt \
    --atomlist atom_indices.txt \
    --topology native.pdb \
    --n-clusters 100 \
    --output-dir ./results
```

### Lumping (MSM Construction)
```bash
hkdm lump pcca \
    --assignments ./results/assignments_kcenters_n_100.txt \
    --traj-len ./results/traj_len.txt \
    --n-macro 4 \
    --homedir ./results
```

## Citation

If you use **APLoD** or **HKDataMiner** in your research, please cite:

```bibtex
@article{liu2017adaptive,
  title={Adaptive partitioning by local density-peaks: An efficient density-based clustering algorithm for analyzing molecular dynamics trajectories},
  author={Liu, Song and Zhu, Lizhe and Sheong, Fu Kit and Wang, Wei and Huang, Xuhui},
  journal={Journal of Computational Chemistry},
  volume={38},
  number={3},
  pages={152--160},
  year={2017},
  publisher={Wiley Online Library}
}
```

## License

This project is licensed under the **Apache 2.0 License**. See the [LICENSE](LICENSE) file for details.

## Authors

*   **Prof. Xuhui Huang** - *Project Leader*
*   **Dr. Song Liu** - *Developer*