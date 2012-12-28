#!/usr/bin/env python
# encoding: utf-8


import numpy
import numpy.matlib
import scipy.stats
import prody
import os
import itertools
import multiprocessing

prody.confProDy(verbosity='none')

n_process = 4


def compute_dist_matrix(points):
    numPoints = len(points)
    distMat = numpy.sqrt(numpy.sum((numpy.matlib.repmat(points, numPoints, 1) - numpy.matlib.repeat(points, numPoints, axis=0))**2, axis=1))
    return distMat.reshape((numPoints,numPoints))


def extract_matrix(inputs):
    pdb_file, secondary_file, ss_type, ss_length = inputs

    path = os.path.join('PDB', pdb_file)
    protein = prody.parsePDB(path, subset='CA')
    coords = protein.getCoords()
    n_residues = coords.shape[0]

    secondary_path = os.path.join('SS', secondary_file)
    ss = open(secondary_path).read()

    results = []

    # loop over all segments
    for start in range(0, n_residues - ss_length + 1):
        end_plus_one = start + ss_length

        if ss[start:end_plus_one] == ss_type * ss_length:
            if ss_length == 20 and ss_type == 'E':
                print pdb_file
            dist_mat = compute_dist_matrix( coords[start:end_plus_one,:] )
            results.append(dist_mat)
    return results


def main():
    # create process pool
    pool = multiprocessing.Pool(n_process)

    # read in list of pdbs
    pdb_list = open('effective_list.txt').readlines()
    pdb_list = [ line.strip() for line in pdb_list ]

    # create list of secondary structures
    ss_list = [ line[:-4]+'.dat' for line in pdb_list ]

    # loop over H, E
    for ss_type in ['H', 'E']:
        for ss_length in range(4, 21):
        #for ss_length in range(5, 6):
            params = zip( pdb_list, ss_list, itertools.repeat(ss_type), itertools.repeat(ss_length) )
            results = pool.map(extract_matrix, params)
            matrices = []
            for r in results:
                matrices.extend(r)

            print ss_type, ss_length, len(matrices)
            # if there are no fragments, then continue
            if not matrices:
                continue

            matrices = numpy.dstack(matrices)

            # extract info
            lower_bounds = numpy.zeros( (ss_length, ss_length) )
            upper_bounds = numpy.zeros( (ss_length, ss_length) )

            for i in range(ss_length):
                for j in range(i, ss_length):
                    lb = scipy.stats.scoreatpercentile( matrices[i,j,:], 5 )
                    ub = scipy.stats.scoreatpercentile( matrices[i,j,:], 95 )
                    lower_bounds[i,j] = lower_bounds[j,i] = lb
                    upper_bounds[i,j] = upper_bounds[j,i] = ub
            bounds = numpy.dstack( [lower_bounds, upper_bounds] )
            # write to disk
            numpy.save( 'Results/bounds_{}_{}.npy'.format(ss_type, ss_length), bounds )


if __name__ == '__main__':
    main()
