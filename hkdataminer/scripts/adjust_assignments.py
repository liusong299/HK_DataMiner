import os,sys
import numpy as np
import argparse

def merge_assignments(new_assignments, old_assignments, remove_outliers=False):
    outliers = -1
    # Number of clusters in assignments, ignoring noise if present.
    #clusters_size = len(set(old_assignments)) - (1 if -1 in old_assignments else 0)
    clusters_size = np.max(old_assignments) + 1
    max_clust_id = clusters_size
    print("max_clust_id:", max_clust_id)
    count_first = [0] * clusters_size
    count_second = [0] * clusters_size

    old_assignments_size = len(old_assignments)
    # new_assignments_size = len(new_assignments)
    for i in range(0, old_assignments_size):
        if old_assignments[i] != outliers:
            if new_assignments[i] != outliers:
                count_first[old_assignments[i]] += 1
            count_second[old_assignments[i]] += 1

    # Percentage
    percentage = [0.0] * clusters_size
    for i in range(0, clusters_size):
        if count_second[i] is 0:
            percentage[i] = 0.0
        else:
            percentage[i] = float(count_first[i])/float(count_second[i])

    # Starting assignment
    assignments=np.copy(old_assignments)
    for i in range(0, old_assignments_size):
        if old_assignments[i] != outliers and percentage[old_assignments[i]] > 0.7:
            if new_assignments[i] != outliers:
                assignments[i] = new_assignments[i] + max_clust_id
                # print old_assignments[i]
            elif remove_outliers is True:  #if want to remove outliers in the iterations
                assignments[i] = outliers

    n_microstates = len(set(assignments)) - (1 if -1 in assignments else 0)
    adjust = np.max(assignments) - n_microstates
    for i in range(len(assignments)):
        if assignments[i] > 0:
            assignments[i] -= adjust

    return assignments

def main():
    cli = argparse.ArgumentParser()
    cli.add_argument('-n', '--new_assignments',  type=str)
    cli.add_argument('-o', '--old_assignments', type=str)
    cli.add_argument('-m', '--merged_assignments', type=str)
    cli.add_argument('-a', '--adjust', type=int)
    args = cli.parse_args()


    #new_assignments = np.loadtxt(args.new_assignments, dtype=np.int32)
    old_assignments = np.loadtxt(args.old_assignments, dtype=np.int32)
    print(old_assignments)
    adjust = args.adjust
    print(adjust)

    for i in range(len(old_assignments)):
        if old_assignments[i] >= 0:
            old_assignments[i] += adjust
    print(old_assignments)
    #merged_assignments = merge_assignments(new_assignments,old_assignments, remove_outliers=True)
    #n_microstates = len(set(merged_assignments)) - (1 if -1 in merged_assignments else 0)
    #print('Estimated number of clusters: %d' % n_microstates)
    np.savetxt(args.old_assignments, old_assignments, fmt="%d")

if __name__ == "__main__":
    main()

