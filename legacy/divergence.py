#!/usr/bin/env python

from __future__ import division
from scipy.special import xlogy
import numpy as np

def jensen_shannon(counts, nstates=None, nsamples=None):
    """ Calculate Jensen-Shannon Divergence from count data.
    """
    #test data
    #np.random.seed(1)
    #counts = np.random.randint(0, 50, (10, N_SAMPLES*N_STATES)).astype(float)

    # We reshape the matrix to group the data into samples and exclude all loci that
    # are not covered by each and every sample by setting it to np.nan. The replacement
    # withth np.nan of values that belong to these excluded loci
    # is crucial and ensures that downstream calculations with numpy arrays do not
    # need to be modified, for example to check for division by zero.

    counts3d = counts.reshape(-1, nsamples, nstates)
    counts3d[~np.all(counts3d.sum(axis=2), axis=1)] = np.nan

    # calculate the mixture entropy
    counts_per_sample = counts3d.sum(axis=1)
    counts_per_locus = counts_per_sample.sum(axis=1)[:, np.newaxis]
    mixture = counts_per_sample/counts_per_locus
    mixture_entropy = -np.sum(xlogy(mixture, mixture), axis=1)

    # calculate the weighted average entropy
    probabilities = counts3d/counts3d.sum(axis=2)[:, :, np.newaxis]
    entropies = -np.sum(xlogy(probabilities, probabilities), axis=2)
    average_entropy = np.sum(entropies*counts3d.sum(axis=2)/counts_per_locus, axis=1)
    
    # Jensen-Shannon Divergence
    return  np.around(mixture_entropy - average_entropy, decimals=3) + 0.

if __name__ == '__main__':
    import sys

    ### INPUT ###

    data_in = sys.stdin
    context = sys.argv[1]
    NSTATES = 2
    NSAMPLES = 7

    dtype = 'a4, uint32, {shape}float32'.format(shape=NSTATES*NSAMPLES)
    chrom, end, counts = np.loadtxt(data_in, dtype=dtype, unpack=True)

    ### PROCESS ###

    jsd = jensen_shannon(counts, nstates=NSTATES, nsamples=NSAMPLES)
    dtype = [('chrom', '|S4'), ('start', 'uint32'), ('end', 'uint32'), ('div', 'float32')]
    X = np.array(zip(chrom, end - 1, end, jsd), dtype=dtype)
    X = X.compress(~np.isnan(jsd), axis=0)

    ### OUTPUT ###

    #data_out = 'total_by_identity.jsd.{context}.bedg'.format(context=context)
    data_out = 'tmp.bedg'
    with open(data_out, 'wb') as bedg:
        trackinfo = 'track type=bedGraph name="Cme Diversity" description="JSD for population in {context} context" db=hg19\n'.format(context=context)
        bedg.write(trackinfo)
    
    with open(data_out, 'a') as bedg:
        fmt = '%s %i %i %2.3f'
        np.savetxt(bedg, X, delimiter='\t', fmt=fmt) 
