#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import csv
import glob
import gzip
import sys
from itertools import izip

import numpy as np
import yaml
from scipy.special import xlogy

### FUNCTIONS ###

def make_pmf(data):
    meth_prob = data.reshape(data.size, 1)
    nometh_prob = np.ones(meth_prob.shape) - meth_prob
    return np.concatenate((meth_prob, nometh_prob), axis=1)

def js_divergence(pmf_sample):
    """Returns the Jensen-Shannon Divergence"""
    pmf_sample = np.asarray(pmf_sample)
    sample_size = pmf_sample.shape[0]
    assert (pmf_sample.sum(axis=1) == np.ones(sample_size)).all()
    weight = 1.0/sample_size # TODO: add other weight options
    mix_pmf = (weight*pmf_sample).sum(axis=0)
    mix_entropy = -(xlogy(mix_pmf, mix_pmf)).sum()
    entropies = -(xlogy(pmf_sample, pmf_sample)).sum(axis=1)
    return mix_entropy - (weight*entropies).sum()

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

    assert len(d.shape) == 2, "Distributions must be a 2d-array"
    assert len(w.shape) == 1, "Weights must be a 1d-array"
    assert d.shape[1] == w.size, "The number of distributions must equal the number of weights."
    assert (w >= 0.).all(), "Your weight vector contains negative weights."
    assert np.isclose(w.sum(), 1.), "Your weight vector must sum to 1."

    dAvg = d.dot(w)
    h_dAvg = -xlogy(dAvg, dAvg).sum()
    hAvg = -xlogy(w*d, d).sum()
    j = h_dAvg - hAvg
    return j

### INPUT DATA ###

with open(sys.argv[1], 'rb') as yamlfile:
    config = yaml.safe_load(yamlfile)

inpath = config['input']['path']
outpath = config['output']['path']

incols = config['input']['columns']
outcols = config['output']['columns']

depth_idx = incols.index('read depth')
score_idx = incols.index('score')

### PROCESS and OUTPUT ###

in_counter = 0
out_counter = 0

# open files

if len(inpath) == 1:
    infiles = [open(f, 'rb') for f in glob.glob(inpath[0])]
elif len(inpath) > 1:
    infiles = [open(f, 'rb') for f in inpath]

if outpath.endswith('.gz'):
    outfile = gzip.open(outpath, 'wb')
else:
    outfile = open(outpath, 'wb')

writer = csv.DictWriter(outfile, fieldnames=outcols, delimiter='\t')

for lines in izip(*infiles):
    """see also http://stackoverflow.com/a/24108858"""
    # "lines" is a list of strings. The map function makes it into a
    # list of lists that have strings (the values) as elements. This
    # list is converted into a 2-d numpy.array "data" which is used to
    # calculate jsd if all samples are covered (have depth) at this
    # site. The dictionary "out" is populated and written as a row to
    # the outfile.
    in_counter += 1

    data = np.array(map(lambda x: x.rstrip('\n').split('\t'), lines))

    depth = data[:, depth_idx].astype(int)

    if depth.all():
        weights = depth/depth.sum()
        scores = data[:, score_idx].astype(float)
        pmfs  = make_pmf(scores)
        #jsd_value = js_divergence(pmfs)
        #jsd_value = jsd(distributions=pmfs.T, weights=np.ones_like(weights)/weights.size)
        jsd_value = jsd(distributions=pmfs.T, weights=weights)

        # get inherited fields and test them for equality in all samples
        mask = np.in1d(incols, outcols)
        inherited_fields = data.compress(mask, axis=1)
        try:
            n_samples, _ = data.shape
            test = [np.array_equal(inherited_fields[0,:], inherited_fields[i,:]) for i in range(n_samples)]
            assert all(test)
        except AssertionError:
            print inherited_fields

        # create output dictionary, fill it and write the data to the file
        out = dict.fromkeys(outcols)
        out['diversity'] = round(np.exp(jsd_value), 2)
        for i, col in enumerate(np.array(incols)[mask]):
            out[col] = inherited_fields[0, i]
        if out['end'] is None:
            out['end'] = int(out['start']) + 1

        writer.writerow(out)
        out_counter += 1
        
        if not out_counter%100000:
            print '{0: 8d} {1: 8d}'.format(in_counter, out_counter)
        #if out_counter == 30: break

# close files
for f in infiles:
    f.close()
outfile.close()
