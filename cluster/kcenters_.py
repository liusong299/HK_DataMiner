__author__ = 'stephen'
#===============================================================================
# GLOBAL IMPORTS:
import os, sys
import numpy as np
import time
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.utils import check_array, check_random_state
#===============================================================================
# LOCAL IMPORTS:
HK_DataMiner_Path = os.path.relpath(os.pardir)
sys.path.append(HK_DataMiner_Path)
from metrics.pairwise import pairwise_distances
#===============================================================================


def k_centers(X, n_clusters=8, metric='rmsd', random_state=None):
    """K-Centers clustering
    Cluster a vector or Trajectory dataset using a simple heuristic to minimize
    the maximum distance from any data point to its assigned cluster center.
    The runtime of this algorithm is O(kN), where k is the number of
    clusters and N is the size of the dataset, making it one of the least
    expensive clustering algorithms available.
    Parameters
    ----------
    n_clusters : int, optional, default: 8
        The number of clusters to form as well as the number of
        centroids to generate.
    metric : {"euclidean", "sqeuclidean", "cityblock", "chebyshev", "canberra",
              "braycurtis", "hamming", "jaccard", "cityblock", "rmsd"}
        The distance metric to use. metric = "rmsd" requires that sequences
        passed to ``fit()`` be ```md.Trajectory```; other distance metrics
        require ``np.ndarray``s.
    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.
    References
    ----------
    .. [1] Gonzalez, Teofilo F. "Clustering to minimize the maximum
       intercluster distance." Theor. Comput. Sci. 38 (1985): 293-306.
    .. [2] Beauchamp, Kyle A., et al. "MSMBuilder2: modeling conformational
       dynamics on the picosecond to millisecond scale." J. Chem. Theory.
       Comput. 7.10 (2011): 3412-3419.
    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features] or md.Trajectory
        Coordinates of cluster centers
    labels_ : array, [n_samples,]
        The label of each point is an integer in [0, n_clusters).
    """
    n_samples = len(X)
    if random_state is -1:
        seed = check_random_state(None).randint(0, n_samples)
    else:
        seed = random_state
    print "seed=", seed
    cluster_centers_ = []
    cluster_centers_.append(seed)  #seed = random
    distances_ = pairwise_distances(X, index=seed, metric=metric)
    labels_ = np.zeros(len(X), dtype=np.int32)

    for i in xrange(1, n_clusters):
        MaxIndex = np.argmax(distances_)
        cluster_centers_.append(MaxIndex)
        #set the furthest point from existing center as a new center

        if distances_[ MaxIndex ] < 0:
            break

        new_distance_list = pairwise_distances(X, index=MaxIndex, metric=metric)
        updated_indices = np.where(new_distance_list < distances_)[0]
        distances_[ updated_indices ] = new_distance_list[ updated_indices ]
        labels_[ updated_indices ] = i

    return cluster_centers_, labels_

class KCenters(BaseEstimator, ClusterMixin):
    """K-Centers clustering
    Cluster a vector or Trajectory dataset using a simple heuristic to minimize
    the maximum distance from any data point to its assigned cluster center.
    The runtime of this algorithm is O(kN), where k is the number of
    clusters and N is the size of the dataset, making it one of the least
    expensive clustering algorithms available.
    Parameters
    ----------
    n_clusters : int, optional, default: 8
        The number of clusters to form as well as the number of
        centroids to generate.
    metric : {"euclidean", "sqeuclidean", "cityblock", "chebyshev", "canberra",
              "braycurtis", "hamming", "jaccard", "cityblock", "rmsd"}
        The distance metric to use. metric = "rmsd" requires that sequences
        passed to ``fit()`` be ```md.Trajectory```; other distance metrics
        require ``np.ndarray``s.
    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.
    References
    ----------
    .. [1] Gonzalez, Teofilo F. "Clustering to minimize the maximum
       intercluster distance." Theor. Comput. Sci. 38 (1985): 293-306.
    .. [2] Beauchamp, Kyle A., et al. "MSMBuilder2: modeling conformational
       dynamics on the picosecond to millisecond scale." J. Chem. Theory.
       Comput. 7.10 (2011): 3412-3419.
    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features] or md.Trajectory
        Coordinates of cluster centers
    labels_ : array, [n_samples,]
        The label of each point is an integer in [0, n_clusters).
    """
    def __init__(self, n_clusters=8,  metric='rmsd', random_state=None):
        self.n_clusters = n_clusters
        self.random_state = random_state
        self.metric = metric
    def fit(self, X, y=None):
        """Perform clustering.
        Parameters
        -----------
        X : array-like, shape=[n_samples, n_features]
            Samples to cluster.
        """
        #X = check_array(X)
        t0 = time.time()
        self.cluster_centers_, self.labels_ = \
            k_centers(X, n_clusters=self.n_clusters,  metric=self.metric, random_state=self.random_state)
        t1 = time.time()
        print "K-Centers clustering Time Cost:", t1 - t0
        return self






