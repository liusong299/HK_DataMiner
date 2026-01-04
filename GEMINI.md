# HKDataMiner

## Project Overview
HKDataMiner is a Python library and toolkit designed for constructing statistical models from biomolecular dynamics data. It provides algorithms for clustering conformational states and building Markov State Models (MSMs) via lumping techniques. Additionally, it includes a specialized "Template Matching" module for analyzing cryo-EM data distributions.

## Key Features

### 1. Clustering (`cluster/`)
Algorithms to group molecular conformations into microstates.
*   **K-Centers:** Standard geometric clustering.
*   **DBSCAN:** Density-based clustering.
*   **Faiss DBSCAN:** Accelerated DBSCAN using the Facebook AI Similarity Search (Faiss) library for large datasets.

### 2. Lumping (`lumping/`)
Techniques to aggregate microstates into kinetically meaningful macrostates.
*   **PCCA / PCCA+:** Perron Cluster Cluster Analysis.
*   **Spectral Clustering:** Graph-based spectral methods.
*   **Ward:** Hierarchical clustering.

### 3. Template Matching (`template_matching/`)
A workflow for determining population distributions in multiple conformations, particularly for cryo-EM datasets.
*   **Workflow:**
    1.  **Select Angle:** PCA-based view angle selection.
    2.  **Generate 2D Images:** Project 3D volumes.
    3.  **Two-Stage Matching:** Algorithm to match experimental data to templates.

## Directory Structure
*   `cluster/`: Implementations of clustering algorithms (`dbscan_.py`, `kcenters_.py`, `faiss_dbscan_.py`).
*   `lumping/`: Implementations of lumping algorithms (`pcca_.py`, `spectral_.py`, `ward_.py`).
*   `metrics/`: Distance metrics (e.g., RMSD) used for clustering.
*   `scripts/`: **Main Entry Points.** Contains executable scripts for running analyses (e.g., `doLumping.py`, `test_kcenters.py`).
*   `template_matching/`: Self-contained module for the template matching workflow.
    *   `TSTM.py`: Main driver script for this module.
*   `utils/`: Utility functions for file I/O (`reader_.py`) and plotting (`plot_.py`).

## Usage & Workflow

### Standard MSM Construction
1.  **Clustering:** Use scripts in `scripts/` to cluster trajectory data.
    ```bash
    python ../scripts/test_kcenters.py -n 500
    ```
    *   *Input:* Trajectory files (often handled via `mdtraj`).
    *   *Output:* `assignments.txt` (cluster labels), `generators.txt` (centers).

2.  **Lumping:** Lump microstates into macrostates.
    ```bash
    python ../scripts/doLumping.py -c assignments.txt -m 4
    ```
    *   *Input:* Cluster assignments from the previous step.
    *   *Output:* Macrostate definitions and plots.

### Template Matching
Run the main script from the `template_matching` directory.
```bash
python TSTM.py --datatype='sim' --vol_size=128
```

## Dependencies
*   **Python:** 3.7+
*   **Core:** `numpy`, `scipy`, `scikit-learn`, `matplotlib`
*   **Bio/MD:** `mdtraj`
*   **Image/Vol:** `xmipp`, `mrcfile`
*   **Other:** `sh`, `faiss` (optional but recommended for performance)

## Notes
*   **File Sync Artifacts:** You may see files like `... (stephen-desktop-linux's conflicted copy ...).py`. These are likely synchronization artifacts and should generally be ignored in favor of the clean filenames.
*   **Configuration:** Most scripts use `argparse` and can be configured via command-line flags. Check individual scripts for specific options.
