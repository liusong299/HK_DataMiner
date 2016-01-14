__author__ = 'stephen'
import numpy as np
import mdtraj as md
import os, sys


def get_subindices(assignments=None, state=None, samples=10):
    'Get Subsamples assignments from same state'
    assignments = np.array(assignments)
    if state is not None:
        indices = np.where(np.array(assignments) == 1)[0]
    else:
        indices = range(0, len(assignments))
    if samples is not None:
        indices = np.random.choice(indices, size=samples)
    return indices

def output_trajs(trajs=None, assignments=None, n_states=None, samples=10, output_name='output', output_type='.pdb' ):
    for i in xrange(0, n_states):
        indices =get_subindices(assignments, i, samples)
        name = output_name + '_' + str(i) + '_' + output_type
        trajs[indices].save(name)


