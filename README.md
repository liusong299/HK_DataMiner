HK DataMiner
===============
New version of HK DataMiner, using a lot of features of scikit-learn

Requirements
--------------
- python 2.7.8
- numpy 1.9.2
- scipy 0.15.1
- MDTraj 1.3.0
- scikit-learn 0.16.1
- matplotlib 1.4.3

**OBS:**
- These requirements were tested and proved to work. Feel free to test
older versions of the requirements.

Documentation
--------------
The Doxygen documentation for this code can be found at
[doc](http://chz382.ust.hk/doc/stephen/HK_DataMiner_new/).

Install Python and Python Packages
----------------------------------
We highly recommend that you download the Python 2.7 version of Anaconda, which is a completely free enterprise-ready Python distribution for large-scale data processing, predictive analytics, and scientific computing.

Design and Usage
-----------------
MD Data Miner is a package of python scripts to build an Markov State Models (MSMs) from a simulation dataset.

MD DataMiner Tutorial
---------------------
1.Introduction and Alanine Dipeptide

2.Move to tutorial directory Assuming you’ve in the MD Data Miner folder, move to the Tutorial directory.

    cd Tutorial

3.Cluster your data The following command clusters your data using the RMSD metric by k-centers clustering method.
  ,lumping by PCCA method and plot the lumping results
    
    python ../scripts/test_kcenters.py  -n 500

After Clustering, the assignments of each conformation are stored as assignments.txt, the cluster centers are sotred as generators.txt.

5.Lump microstates into macrostates Once we have determined the number of macrostates in the system, we will use Perron Cluster CLuster Analysis (PCCA) algorithm to lump microstates into macrostates. We use the doLumping.py script to lump the microstates into macrostates. We could define the number of macrostates using -n option.The command below will build 4 macrostates.

    python ../scripts/doLumping.py -c assignments_kcenters_n_500.txt -m 6
Examining the macrostate decomposition It is known that the relevant degrees of freedom for alanine dipeptide are the phi and psi backbone angles. Let’s compute the dihedrals and plot the conformations of each macrostate.
    

TODO
--------------
- [ ] Update Tutorial
- [ ] Test more

labels:
- [DON] - Done
- [IMP] - To be improoved
- [BUG] - Buggy and experimental

Credits
--------------
**Author:** [Stephen Liu](http://chz382.ust.hk/u/stephen): *liusong299@gmail.com*

**Contributors:** Maybe you !
