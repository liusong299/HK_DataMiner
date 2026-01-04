import os
import numpy as np
from hkdataminer.cluster import KCenters, APLoD
from hkdataminer.lumping import PCCA
from hkdataminer.utils import XTCReader, plot_cluster, utils

def run_clustering(
    trajListFns='trajlist',
    atomListFns='atom_indices',
    topology='native.pdb',
    homedir='.',
    iext='xtc',
    n_clusters=100,
    stride=None,
    output_dir='.'
):
    print(f"Running Clustering with n_clusters={n_clusters}, stride={stride}")
    
    # Check cache for phi/psi
    phi_file = os.path.join(output_dir, "phi_angles.txt")
    psi_file = os.path.join(output_dir, "psi_angles.txt")
    traj_len_file = os.path.join(output_dir, "traj_len.txt")
    
    # Always read trajs if using RMSD
    print("Reading trajectories...")
    trajreader = XTCReader(trajListFns, atomListFns, homedir, iext, topology, nSubSample=stride)
    trajs = trajreader.trajs
    traj_len = trajreader.traj_len
    np.savetxt(traj_len_file, traj_len, fmt="%d")
    
    if not (os.path.isfile(phi_file) and os.path.isfile(psi_file)):
        phi_angles, psi_angles = trajreader.get_phipsi(trajs, psi=[6, 8, 14, 16], phi=[4, 6, 8, 14])
        np.savetxt(phi_file, phi_angles, fmt="%f")
        np.savetxt(psi_file, psi_angles, fmt="%f")
    else:
        phi_angles = np.loadtxt(phi_file, dtype=np.float32)
        psi_angles = np.loadtxt(psi_file, dtype=np.float32)

    # do Clustering using KCenters method
    print(f"Clustering with KCenters (n={n_clusters})...")
    cluster = KCenters(n_clusters=n_clusters, metric="rmsd", random_state=0)
    cluster.fit(trajs)

    labels = cluster.labels_
    n_microstates = len(set(labels)) - (1 if -1 in labels else 0)
    print(f'Estimated number of clusters: {n_microstates}')

    cluster_centers_ = cluster.cluster_centers_
    
    clustering_name = "kcenters_n_" + str(n_microstates)
    
    # Save outputs
    assign_file = os.path.join(output_dir, f"assignments_{clustering_name}.txt")
    centers_file = os.path.join(output_dir, f"cluster_centers_{clustering_name}.txt")
    pdb_file = os.path.join(output_dir, "cluster_centers.pdb")
    
    np.savetxt(assign_file, labels, fmt="%d")
    np.savetxt(centers_file, cluster_centers_, fmt="%d")
    
    original_cwd = os.getcwd()
    try:
        os.chdir(output_dir)
        plot_cluster(labels=labels, phi_angles=phi_angles, psi_angles=psi_angles, name=clustering_name)
        trajs[cluster_centers_].save(pdb_file)
    finally:
        os.chdir(original_cwd)
        
    return assign_file

def run_aplod(
    trajListFns='trajlist',
    atomListFns='atom_indices',
    topology='native.pdb',
    homedir='.',
    iext='xtc',
    rho_cutoff=1.0,
    delta_cutoff=1.0,
    n_neighbors=100,
    stride=None,
    output_dir='.'
):
    print(f"Running APLoD Clustering with rho={rho_cutoff}, delta={delta_cutoff}")
    
    # Check cache for phi/psi
    phi_file = os.path.join(output_dir, "phi_angles.txt")
    psi_file = os.path.join(output_dir, "psi_angles.txt")
    traj_len_file = os.path.join(output_dir, "traj_len.txt")
    
    print("Reading trajectories...")
    trajreader = XTCReader(trajListFns, atomListFns, homedir, iext, topology, nSubSample=stride)
    trajs = trajreader.trajs
    traj_len = trajreader.traj_len
    np.savetxt(traj_len_file, traj_len, fmt="%d")
    
    if not (os.path.isfile(phi_file) and os.path.isfile(psi_file)):
        phi_angles, psi_angles = trajreader.get_phipsi(trajs, psi=[6, 8, 14, 16], phi=[4, 6, 8, 14])
        np.savetxt(phi_file, phi_angles, fmt="%f")
        np.savetxt(psi_file, psi_angles, fmt="%f")
    else:
        phi_angles = np.loadtxt(phi_file, dtype=np.float32)
        psi_angles = np.loadtxt(psi_file, dtype=np.float32)

    # Use Phi/Psi for APLoD because RMSD is broken/disabled in APLoD without knnn
    print("Clustering with APLoD (using Phi/Psi angles)...")
    data = np.column_stack((phi_angles, psi_angles))
    
    cluster = APLoD(rho_cutoff=rho_cutoff, delta_cutoff=delta_cutoff, n_neighbors=n_neighbors, metric='euclidean')
    cluster.fit(data)

    labels = cluster.labels_
    n_microstates = len(set(labels)) - (1 if -1 in labels else 0)
    print(f'Estimated number of clusters: {n_microstates}')

    cluster_centers_indices = cluster.cluster_centers_
    # For APLoD, cluster_centers_ are indices into the data
    
    clustering_name = "aplod_n_" + str(n_microstates)
    
    # Save outputs
    assign_file = os.path.join(output_dir, f"assignments_{clustering_name}.txt")
    centers_file = os.path.join(output_dir, f"cluster_centers_{clustering_name}.txt")
    pdb_file = os.path.join(output_dir, "cluster_centers.pdb")
    
    np.savetxt(assign_file, labels, fmt="%d")
    # Save indices of centers
    np.savetxt(centers_file, cluster_centers_indices, fmt="%d")
    
    original_cwd = os.getcwd()
    try:
        os.chdir(output_dir)
        plot_cluster(labels=labels, phi_angles=phi_angles, psi_angles=psi_angles, name=clustering_name)
        # Check if indices are valid before saving pdb
        if len(cluster_centers_indices) > 0:
             trajs[cluster_centers_indices].save(pdb_file)
    finally:
        os.chdir(original_cwd)
        
    return assign_file

def run_lumping(
    assignments_file,
    traj_len_file='traj_len.txt',
    n_macro_states=6,
    homedir='.'
):
    print(f"Running Lumping with n_macro_states={n_macro_states}")
    
    labels = np.loadtxt(assignments_file, dtype=np.int32)
    traj_len = np.loadtxt(traj_len_file, dtype=np.int32)
    
    n_microstates = len(set(labels)) - (1 if -1 in labels else 0)
    print(f'Estimated number of clusters: {n_microstates}')
    
    # We assume phi/psi angles are available in homedir/current dir for plotting?
    # doLumping.py reads them from ./phi_angles.txt
    
    print("Running PCCA...")
    PCCA_lumper = PCCA(n_macro_states=n_macro_states, traj_len=traj_len, cut_by_mean=False, homedir=homedir)
    algorithm = PCCA_lumper
    algorithm.fit(labels)
    
    lumping_name = 'micro_' + str(n_microstates) + "_PCCA_" + str(n_macro_states)
    
    # Plotting? doLumping.py does NOT call plot_cluster by default (commented out).
    # But tutorial step 4 says "Let’s compute the dihedrals and plot the conformations of each macrostate."
    # The script has commented out plots.
    # I will enable plotting if phi/psi available.
    
    phi_file = os.path.join(homedir, "phi_angles.txt")
    psi_file = os.path.join(homedir, "psi_angles.txt")
    
    if os.path.isfile(phi_file) and os.path.isfile(psi_file):
        phi_angles = np.loadtxt(phi_file, dtype=np.float32)
        psi_angles = np.loadtxt(psi_file, dtype=np.float32)
        
        # PCCA object saves results via OutputResult which calls Evaluate_Result which calls plot_matrix.
        # So some plotting is done.
        # But macrostate plotting (phi/psi) is not called in doLumping.py active code.
        # I'll add it.
        
        original_cwd = os.getcwd()
        try:
            os.chdir(homedir)
            plot_cluster(labels=algorithm.MacroAssignments_, phi_angles=phi_angles, psi_angles=psi_angles, name=lumping_name)
        finally:
            os.chdir(original_cwd)