#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from scipy.special import xlogy
import sys

def countmeth(basecalls): # maybe use np.fromregex if you read directly from filtered generator
    """
    """
    
    ref = basecalls[2]
    callstrings = basecalls[3:]
    if ref.upper() == 'C':
        strand = '+'
        counts_M = [call.count('.') for call in callstrings]
        counts_m = [call.count('T') for call in callstrings]
    elif ref.upper() == 'G':
        strand = '-'
        counts_M = [call.count(',') for call in callstrings]
        counts_m = [call.count('a') for call in callstrings]
    total_count = sum(counts_M) + sum(counts_m)
    counts = zip(counts_M, counts_m)
    return basecalls[:2] + [strand] + [total_count] + counts

def diversity(counts):
    """
    """

    position = counts[:3]
    
    total_count = counts[3]
    
    count_matrix = np.array(counts[4:])
    
    num_ind, dim_sampleSpace = count_matrix.shape
    
    wavg_methylation = (count_matrix.sum(axis=0)/total_count)[0]
    
    ind_counts = count_matrix.sum(axis=1)
    
    freq_matrix = count_matrix.T.compress(ind_counts != 0, axis=1)/ind_counts[ind_counts != 0]

    ind_weights = ind_counts[ind_counts != 0]/total_count

    #avg_count = round(np.average(ind_counts), 1)    
    #entropy_weight_dist = -xlogy(ind_weights, ind_weights).sum()
    #confidence_metric = round(avg_count/float(num_ind)*(np.log(dim_sampleSpace) - entropy_weight_dist), 1)
    
    diversity = jsd(distributions=freq_matrix, weights=ind_weights) # jsd could be vectorized/mapped over loci
    
    return position + [str(wavg_methylation), str(diversity)]


from scipy.special import xlogy

def jsd(distributions=None, weights=None):
    """ Calculates the Jensen-Shannon Divergence.
    
    Parameters
    ----------
    d : array_like
        d.shape = (k,n)
        n distributions over k-dimensional sample space
    w : array_like
        w.shape = (n,)
        weights for n distributions

    Returns
    -------
    j : float
        Jensen-Shannon Divergence between the distributions

    Notes
    -----
    The Jensen-Shannon Divergence is given by J = H(<P>) - <H>, the
    entropy of the weighted average distribution minus the weighted
    average of the distribution entropies. It is a measure of the
    overall difference between distributions based on information
    loss.

    """

    d = distributions
    w = weights

    dAvg = d.dot(w)
    h_dAvg = -xlogy(dAvg, dAvg).sum()
    hAvg = -xlogy(w*d, d).sum()
    j = h_dAvg - hAvg
    return j

if __name__ == '__main__':
    
    with open(sys.argv[1], 'r') as mpileup:
        pileup = (line.strip().split('\t') for line in mpileup)
        base_calls = (line[:3] + line[4:][::3] for line in pileup if line[2].upper() in ['C', 'G'])
        base_counts = (countmeth(calls) for calls in base_calls)
        covered_bases = (counts for counts in base_counts if counts[3] > 0)
        div = (diversity(counts) for counts in covered_bases)
        c = 0
        for line in div:
            print '\t'.join(line)
#            c += 1
#            if c==30:
#                break
