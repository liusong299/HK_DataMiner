# HKDataMiner

**HKDataMiner** is a Python library for constructing statistical models from biomolecular dynamics data. It provides algorithms for clustering conformational states and building Markov State Models (MSMs) via lumping techniques.

This repository has been modernized to support **Python 3.10+** and includes a CLI and reproducible tutorial.

## Quickstart

### 1. Installation

We recommend using `mamba` (or `conda`) to manage dependencies.

```bash
# Clone the repository
git clone https://github.com/liusong299/HK_DataMiner.git
cd HK_DataMiner

# Create the environment
mamba env create -f environment.yml

# Activate the environment
mamba activate hkdataminer-py310

# Install the package in editable mode
pip install -e .
```

### 2. Run the Tutorial (One Command)

To run the end-to-end tutorial (Alanine Dipeptide):

```bash
hkdm tutorial
```

This command will:
1.  Unpack the included tutorial data (`Tutorial.tar.gz`).
2.  Run K-Centers clustering on the trajectory data.
3.  Run PCCA lumping to identify macrostates.
4.  Generate plots and assignment files in `_tutorial_run/Tutorial`.

## Tutorial Walkthrough

The tutorial analyzes Alanine Dipeptide dynamics. The raw data includes trajectory files (`.xtc`) and topology (`native.pdb`).

### Step 1: Clustering (Microstates)

The first step clusters the conformations into microstates using the K-Centers algorithm based on RMSD.

**CLI Command:**
```bash
hkdm cluster kcenters --n-clusters 100
```

**Key Outputs:**
*   `assignments_kcenters_n_100.txt`: Cluster assignment for each frame.
*   `cluster_centers.pdb`: PDB file containing the center conformation of each cluster.
*   `kcenters_n_100.png`: Plot of microstates on the Phi-Psi landscape.

### Step 2: Lumping (Macrostates)

The second step lumps microstates into kinetically metastable macrostates using PCCA.

**CLI Command:**
```bash
hkdm lump pcca --assignments assignments_kcenters_n_100.txt --n-macro 4
```

**Key Outputs:**
*   `PCCA_100_to_4_states_MacroAssignments.txt`: Macrostate assignment for each frame.
*   `PCCA_100_to_4_states_Matrix.png`: Transition probability matrix plot.
*   `PCCA_100_to_4_states_Metastability_Modularity.txt`: Quality metrics.

### Step 3: Visualization

You can view the generated `.png` files to see the energy landscape and cluster distributions. A Jupyter Notebook is also provided in `notebooks/AlanineDipeptide_Tutorial.ipynb`.

## Troubleshooting

*   **Missing Dependencies:** Ensure you have activated the conda environment (`mamba activate hkdataminer-py310`).
*   **Compilation Errors:** If `pip install -e .` fails, ensure you have a C compiler installed (usually present in standard Linux/macOS environments).
*   **Import Errors:** If you see `ModuleNotFoundError`, try reinstalling the package: `pip install -e .`.

## Optional: Cryo-EM Template Matching

The repository contains a `template_matching` module for analyzing cryo-EM data. This requires **Xmipp** and is currently optional.

To run the template matching workflow (Linux + Xmipp required):
1.  Navigate to `template_matching/`.
2.  Follow instructions in the original documentation or scripts.
3.  Run `python TSTM.py ...` (Note: This part has not been fully modernized for Python 3.10 and may require adjustments).

## License

This project is licensed under the Apache-2.0 License. See the [LICENSE](LICENSE) file for details.

## Authors

*   **Prof. Xuhui Huang** - *Project Leader*
*   **Mr. Song Liu** - *Developer*
*   **Mr. Hanlin Gu** - *Developer*