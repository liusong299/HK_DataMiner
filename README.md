# HKDataminer Tutorial - Data Mining Tool for Biomolecular Dynamics
**Version 0.9.1, (c) Song LIU**
## Overview
HK Dataminer is a python library for constructing statistical models for bimolecular dynamics data. The package includes complementary data mining algorithms such as clustering and Markov state models (MSMs).

# Installation guide
## I. Package requirement
* [Python v3.7](https://www.python.org)
* [NumPy v1.1.8](https://numpy.org)
* [SciPy v1.17.0](https://www.scipy.org)
* [MDTraj 1.9.3](mdtraj.org)
* [Scikit-Learn v0.21.3](https://scikit-learn.org)
* [Matplotlib](https://matplotlib.org)
## II.Installing
We highly recommend that you download the Python 3.7 version of Anaconda, which is a completely free enterprise-ready Python distribution for large-scale data processing, predictive analytics, and scientific computing.

# A step-by-step guide of running HKDataminer
## 1.Introduction and Alanine Dipeptide

## 2.Move to tutorial directory Assuming you’ve in the MD Data Miner folder, move to the Tutorial directory.
`cd Tutorial`

## 3.Cluster your data The following command clusters your data using the RMSD metric by k-centers clustering method.
`python ../scripts/test_kcenters.py  -n 500`

After Clustering, the assignments of each conformation are stored as assignments.txt, the cluster centers are sotred as generators.txt.

## 4.Lump microstates into macrostates
Once we have determined the number of macrostates in the system, we will use Perron Cluster CLuster Analysis (PCCA) algorithm to lump microstates into macrostates. We use the doLumping.py script to lump the microstates into macrostates. We could define the number of macrostates using -n option.The command below will build 4 macrostates.

`python ../scripts/doLumping.py -c assignments_kcenters_n_500.txt -m 6`

Examining the macrostate decomposition It is known that the relevant degrees of freedom for alanine dipeptide are the phi and psi backbone angles. Let’s compute the dihedrals and plot the conformations of each macrostate.

# Deployment
HKUST method is developed by [Prof. Xuhui Huang's group](http://compbio.ust.hk)

# Authors
* **Prof. Xuhui Huang** - *Project leader* - [xuhuihuang](http://compbio.ust.hk/public_html/pmwiki-2.2.8/pmwiki.php?n=People.XuhuiHuang)
* **Mr. Song Liu** - *Developer* -[liusong299](https://github.com/liusong299/)

See also the list of [contributors](https://github.com/liusong299/gromacs-2019-CWBSol/graphs/contributors) who participated in this project.

# License
This project is licensed under the GPL License - see the [LICENSE](LICENSE) file for details

# Acknowledgments

# References:

