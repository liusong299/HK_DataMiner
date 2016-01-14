__author__ = 'stephen'
import numpy as np
import scipy.io
import scipy.sparse
import scipy.sparse.linalg
from sklearn.base import BaseEstimator, ClusterMixin


class MarkovStateModel(BaseEstimator, ClusterMixin):
    '''
    Reversible Markov State Model
    This model fits a first-order Markov model to a dataset of integer-valued
    timeseries. The key estimated attribute, ``transmat_`` is a matrix
    containing the estimated probability of transitioning between pairs
    of states in the duration specified by ``lag_time``.
    Unless otherwise specified, the model is constrained to be reversible
    (satisfy detailed balance), which is appropriate for equilibrium chemical
    systems.
    Parameters
    ----------
    lag_time : int
        The lag time of the model
    '''
    def __init__(self, lag_time=1, n_macro_states=None, traj_len=None):
        self.homedir='.'
        self.lag_time = lag_time
        self.traj_len = traj_len
        # Cached eigensystem
        self.eigenvalues = None
        self.right_eigenvectors = None
        self.n_micro_states = None
        self.n_macro_states = n_macro_states
        self.tCount_ = None
        self.tProb_ = None

    '''
    def getTransitionMatrix(self, assignments, lag_time=1, outlier=-1):
        """Get the transition matrix.
        """
        # *** Building Transition Count Matrix
        assignments_len = len(assignments)
        n_micro_states = np.max(assignments)+1
        self.n_micro_states = n_micro_states
        tCount_ = scipy.sparse.lil_matrix(n_micro_states, n_micro_states, dtype=np.int32)
        tProb_ = scipy.sparse.lil_matrix(n_micro_states, n_micro_states, dtype=np.float32)
        tCount_ = np.zeros((n_micro_states, n_micro_states))
        tProb_ = np.zeros((n_micro_states, n_micro_states))

        total_count = 0
        assignments_begin = 0
        new_assignments = []
        outlier_flag = False
        if self.traj_len is not None and isinstance(self.traj_len, list):
            #Compute transition count matrix ONLY in one traj
            for per_traj_len in self.traj_len:
                per_traj_assigments = []
                outlier_flag = False
                for i in xrange(0, per_traj_len):
                    index = assignments_begin + i
                    cur_item = int(assignments[index])
                    if cur_item is not outlier:
                        if outlier_flag is True:
                            outlier_flag = False
                            new_assignments.append(per_traj_assigments)
                            #print "outlier", len(per_traj_assigments), len(new_assignments)
                            per_traj_assigments = []
                        if outlier_flag is False:
                            per_traj_assigments.append(cur_item)
                    else:
                        outlier_flag = True
                    if i == assignments_len-1:
                        new_assignments.append(per_traj_assigments)
                        #print len(per_traj_assigments), len(new_assignments)
                        per_traj_assigments = []
                assignments_begin += per_traj_len
        else:   # no traj_len or just have only one traj
            #print "**DEBUG** no traj_len"
            per_traj_assigments = []
            outlier_flag = False
            #print "assignments_len", assignments_len
            for i in xrange(0, assignments_len):
                index = assignments_begin + i
                cur_item = int(assignments[index])
                if cur_item is not outlier:
                    if outlier_flag is True:
                        outlier_flag = False
                        if per_traj_assigments: #if per_traj_assignments is not empty
                            new_assignments.append(per_traj_assigments)
                            #print "outlier", len(per_traj_assigments), len(new_assignments)
                        per_traj_assigments = []
                    if outlier_flag is False:
                        per_traj_assigments.append(cur_item)
                        #if i == assignments_len-1:
                        #    print "i=", i, "cur_item:", cur_item
                        #    new_assignments.append(per_traj_assigments)
                        #    #print len(per_traj_assigments), len(new_assignments)
                        #    per_traj_assigments = []

                else:   # cur_item is outlier
                    outlier_flag = True

                if i == assignments_len-1:
                    new_assignments.append(per_traj_assigments)
                    #print len(per_traj_assigments), len(new_assignments)
                    per_traj_assigments = []
        #print "new_assignmenet shape:", len(new_assignments)
        #print new_assignments
        for i in xrange(0, len(new_assignments)):
            cur_assignments = new_assignments[i]
            for j in xrange(0, len(cur_assignments) - lag_time):
                origState = cur_assignments[j]
                newState = cur_assignments[j + lag_time]
                #if origState is not outlier and newState is not outlier:
                tCount_[origState, newState] += 1
                #print origState, "->", newState, "count:", tCount_[origState, newState]
                total_count += 1

        #print "tCount:"
        #print tCount_
        scipy.io.mmwrite("./tCount_test.mtx", scipy.sparse.csr_matrix(tCount_), field='integer')
        # *** Fulfilling Detailed Balance
        # determine the deviation from reversibility (or symmetry of count matrix)
        d = abs(tCount_ - tCount_.transpose()).sum() / 2.0
        #print " Deviation from reversibility/symmetry: %f" % (float(d)/total_count)

        # Symmetrize count matrix (optional)
        tCount_symmerized = (tCount_+tCount_.transpose())/2

        # reweight to get transition probability matrix
        tProb_ = tCount_symmerized.copy()
        weights_ = np.array(tCount_.sum(axis=1)).flatten()
        # adjust weights for low counts
        weights_[weights_==0] = 1

        # normalize to get tProb
        indices = np.array(tProb_.nonzero()).transpose()
        for ind in indices:
            tProb_[ind[0], ind[1]] /= weights_[ind[0]]
        return tCount_symmerized, tProb_
    '''

    def getTransitionMatrix(self, assignments, lag_time=1, outlier=-1):
        """Get the transition matrix.
        """
        # *** Building Transition Count Matrix
        assignments_len = len(assignments)
        n_micro_states = np.max(assignments)+1
        self.n_micro_states = n_micro_states
        tCount_ = scipy.sparse.lil_matrix(n_micro_states, n_micro_states, dtype=np.int32)
        tProb_ = scipy.sparse.lil_matrix(n_micro_states, n_micro_states, dtype=np.float128)
        tCount_ = np.zeros((n_micro_states, n_micro_states))
        tProb_ = np.zeros((n_micro_states, n_micro_states))

        total_count = 0
        assignments_begin = 0
        #print "traj_len", self.traj_len
        if self.traj_len is not None:
            #Compute transition count matrix ONLY in one traj
            for per_traj_len in self.traj_len:
                for i in range(0, per_traj_len - lag_time):
                    index = assignments_begin + i
                    origState = assignments[index]
                    newState = assignments[index + lag_time]
                    if origState is not outlier and newState is not outlier:
                        tCount_[origState, newState] += 1
                        total_count += 1
                assignments_begin += per_traj_len
        else:
            for i in range(0, assignments_len - lag_time):
                origState = assignments[i]
                newState = assignments[i + lag_time]
                if origState is not outlier and newState is not outlier:
                    tCount_[origState, newState] += 1
                    total_count += 1

        #This code from MSMBuilder 3.0
        #=============================================================================
        tCount_symmerized = 0.5 * (tCount_ + tCount_.T)

        populations = tCount_symmerized.sum(axis=0)
        populations /= populations.sum(dtype=float)
        tProb_ = tCount_symmerized.astype(float) / tCount_symmerized.sum(axis=1)[:, None]
        #=============================================================================
        '''
        # *** Fulfilling Detailed Balance
        # determine the deviation from reversibility (or symmetry of count matrix)
        d = abs(tCount_ - tCount_.transpose()).sum() / 2.0
        #print " Deviation from reversibility/symmetry: %f" % (float(d)/total_count)

        # Symmetrize count matrix (optional)
        tCount_symmerized = (tCount_+tCount_.transpose())/2.0

        # reweight to get transition probability matrix
        tProb_ = tCount_symmerized.copy()
        weights_ = np.array(tCount_symmerized.sum(axis=1)).flatten()
        # adjust weights for low counts
        #weights_[weights_==0] = 1

        # normalize to get tProb
        indices = np.array(tProb_.nonzero()).transpose()
        for ind in indices:
            tProb_[ind[0], ind[1]] /= weights_[ind[0]]
        '''
        return tCount_symmerized, tProb_

    def get_righteigenvectors(self, t_matrix, n_macro_states, tol=1E-30):
        """Get the right eigenvectors of a transition matrix, sorted by
        eigenvalue magnitude
        """
        #if scipy.sparse.issparse(t_matrix):
        #print "Is sparse matrix"
        #values, vectors = scipy.sparse.linalg.eigs(t_matrix, n_macro_states, which="LR", maxiter=100000, tol=tol)
        #t_matrix = t_matrix.astype(np.float64)
        values, vectors = scipy.sparse.linalg.eigs(t_matrix, n_macro_states, maxiter=1000000, tol=tol)
        #values, vectors = scipy.linalg.eig(t_matrix)
        #values, vectors = np.linalg.eig(t_matrix)
        #else:
        #    print "Is NOT sparse matrix"
        #    values, vectors = scipy.linalg.eig(t_matrix)

        #values, vectors = scipy.linalg.eigs(t_matrix)

        order = np.argsort(-np.real(values))
        e_lambda = values[order]
        e_vectors = vectors[:, order]

        # normalize the first eigenvector (populations)
        #e_vectors[:, 0] /= sum(e_vectors[:, 0])

        e_lambda = np.real(e_lambda)
        e_vectors = np.real(e_vectors)
        return e_lambda, e_vectors

    def get_implied_timescales(self, n_implied_times=50, homedir='.'):
        """
            Calculate implied timescales
        """
        for lag_time in range(1, n_implied_times):
            self.lagTime = lag_time  # a bug will cause if you remove this command, Stephen fixed 20141208
            print "Calculating implied timescales at lagtime ", lag_time
            print "- getting Transition Matrix ... ",
            self.getTransitionMatrix()
            n_eigenvectors = n_implied_times + 1
            print "Done."
            print "- getting Eigenvectors ...",
            e_values = self.get_eigenvectors(self.tProb, n_eigenvectors, epsilon=1)[0]
            print "Done."
            print "- Calculating Implied Timescales ...",
                # Correct for possible change in n_eigenvectors from trimming
            n_eigenvectors = len(e_values)
            n_implied_times = n_eigenvectors - 1

            # make sure to leave off equilibrium distribution
            lag_times = lag_time * np.ones(n_implied_times)
            imp_times = -lag_times / np.log(e_values[1: n_eigenvectors])
            print "Done."
            # save implied timescales results to file.
            fn = homedir + "/" + "lag." + str(lag_time) + ".tau.dat"
            print "- Saving results to file:", fn, " ...",
            f = open(fn, 'w')
            for i in range(n_implied_times):
                line = "%d %f\n" %(lag_time, imp_times[i])
                #line = "%d %f\n" %(lag_time, e_values[i])
                f.write(line)
            f.close()
            print "Done."
            print "Calculation Done."

    def fit(self, assignments, y=None):
        """Estimate model parameters.
        Parameters
        ----------
        assignments : list of array-like
            List of assignments, or a single assignments. Each assignment should be a
            1D iterable of state labels. Labels can be integers, strings, or
            other orderable objects.
        Returns
        -------
        self
        Notes
        -----
        `None` and `NaN` are recognized immediately as invalid labels.
        Therefore, transition counts from or to a sequence item which is NaN or
        None will not be counted. The mapping_ attribute will not include the
        NaN or None.
        """
        self.tCount_, self.tProb_ = self.getTransitionMatrix(assignments=assignments, lag_time=self.lag_time, outlier=-1)
        return self
