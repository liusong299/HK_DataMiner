import numpy as np
import scipy.sparse
import scipy.linalg

class MarkovStateModel:
    def __init__(self, lag_time=1, n_macro_states=None, traj_len=None):
        self.lag_time = lag_time
        self.n_macro_states = n_macro_states
        self.traj_len = traj_len
        self.tCount_ = None
        self.tProb_ = None
        self.n_micro_states = 0
        self.assignments = None

    def fit(self, assignments):
        """
        Construct the transition count matrix and transition probability matrix.
        """
        self.assignments = assignments
        # Determine number of microstates
        unique_states = np.unique(assignments)
        if -1 in unique_states:
            unique_states = unique_states[unique_states != -1]
        self.n_micro_states = len(unique_states)
        
        # We assume states are 0..N-1. If not, we might need mapping, but let's assume they are compact or use max+1
        max_state = np.max(assignments)
        n_states = max_state + 1
        
        self.tCount_ = scipy.sparse.lil_matrix((n_states, n_states), dtype=np.int32)
        
        # Handle trajectories
        if self.traj_len is None:
            # Assume single trajectory
            traj_len = [len(assignments)]
        else:
            traj_len = self.traj_len
            
        start_idx = 0
        for length in traj_len:
            end_idx = start_idx + length
            traj = assignments[start_idx:end_idx]
            # Count transitions with lag_time
            for i in range(len(traj) - self.lag_time):
                s_curr = traj[i]
                s_next = traj[i + self.lag_time]
                if s_curr != -1 and s_next != -1:
                    self.tCount_[s_curr, s_next] += 1
            start_idx = end_idx
            
        self.tCount_ = self.tCount_.tocsr()
        
        # Compute Transition Probability Matrix
        row_sums = np.array(self.tCount_.sum(axis=1)).flatten()
        # Avoid division by zero
        row_sums[row_sums == 0] = 1.0
        
        # Normalize
        self.tProb_ = scipy.sparse.diags(1.0 / row_sums).dot(self.tCount_)
        
        return self

    def get_righteigenvectors(self, tProb, n_modes):
        """
        Compute the top n_modes right eigenvectors.
        """
        if scipy.sparse.issparse(tProb):
            vals, vecs = scipy.sparse.linalg.eigs(tProb, k=n_modes, which='LR')
        else:
            vals, vecs = scipy.linalg.eig(tProb)
            
        # Sort by eigenvalue magnitude (descending)
        idx = np.argsort(np.abs(vals))[::-1]
        vals = vals[idx]
        vecs = vecs[:, idx]
        
        return vals[:n_modes].real, vecs[:, :n_modes].real
