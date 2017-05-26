# -*- coding: utf-8 -*-

from __future__ import division
from scipy.special import xlogy
import numpy as np

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


def jensen_shannon(bedcount, nstates=2):
    """ Calculate Jensen-Shannon Divergence from count data.
    """

    #test data
    #np.random.seed(1)
    #counts = np.random.randint(0, 50, (10, N_SAMPLES*N_STATES)).astype(float)

    # We reshape the matrix to group the data into samples and exclude
    # all loci that are not covered by each and every sample by
    # setting it to np.nan. The replacement withth np.nan of values
    # that belong to these excluded loci is crucial and ensures that
    # downstream calculations with numpy arrays do not need to be
    # modified, for example to check for division by zero.

    counts = bedcount.values
    nloci, ncounts = counts.shape
    nsamples = ncounts/nstates

    counts3d = counts.reshape(-1, nsamples, nstates)
    counts3d[~np.all(counts3d.sum(axis=2), axis=1)] = np.nan

    # calculate the entropy of the total population mixture
    counts_per_sample = counts3d.sum(axis=1)
    counts_per_locus = counts_per_sample.sum(axis=1)[:, np.newaxis]
    mixture = counts_per_sample/counts_per_locus
    mixture_entropy = entropy(mixture, axis=1)#-np.sum(xlogy(mixture, mixture), axis=1)

    # calculate the weighted average entropy
    probabilities = counts3d/counts3d.sum(axis=2)[:, :, np.newaxis]
    entropies = entropy(probabilities, axis=2)#-np.sum(xlogy(probabilities, probabilities), axis=2)
    average_entropy = np.sum(entropies*counts3d.sum(axis=2)/counts_per_locus, axis=1)

    # Jensen-Shannon Divergence
    jsd = mixture_entropy - average_entropy
    return  np.around(np.exp(jsd), decimals=3) + 0.

def entropy(pmf, axis=1):
    """Shannon Entropy
    TODO: add option for units (nats, bits, dits)."""
    # try:
    #     msg = 'The probability mass does not sum to 1!'
    #     np.testing.assert_allclose(pmf.sum(axis=axis), 1.0, rtol=1e-04, err_msg=msg)
    # except AssertionError, err_msg:
    #     print err_msg
    return -np.sum(xlogy(pmf, pmf), axis=axis)

# ============================================

def count_data(data, samples, num_states=2):
    """Extract counts from DataFrame.

    We use pandas' DataFrame.filter method to get the subset defined
    by the sample IDs followed by extracting the count values as
    ndarray. Reshape the count array for convenient vectorization of
    downstream calculations. Exclude the loci with 'missing data'
    (i.e. with 0 counts in all states for the given set of samples) by
    setting the respective subarrays to 'np.nan'. This ensures sane
    behaviour in downstream calculations and avoids
    e.g. division-by-zero warnings when calculating the mixture
    distributions. Also we do not exclude the NaNs until we write the
    output data to ensure that the calculated data maps to the
    chromosome and position information in the input.

    Parameters
    ----------
    data : pandas DataFrame
        Column labels for the count data should follow the convention
        "<samples>_<state>"
    samples : list-like
        Column labels to filter for

    Returns
    -------
    counts3d : ndarray
        axis 0 = locus
        axis 1 = sample
        axis 2 = state

    Examples
    --------

    """

    num_samples = len(samples)

    samples_regex = '|'.join([sample + '_' for sample in samples])

    counts = data.filter(regex=samples_regex).values.astype(float)

    num_loci, num_values = counts.shape

    assert num_samples*num_states == num_values

    counts3d = counts.reshape(num_loci, num_samples, num_states)

    missing_loci = ~np.all(counts3d.sum(axis=2), axis=1)

    counts3d[missing_loci] = np.nan

    return counts3d

def mixture_distribution(counts3d):
    """
    """

    counts_per_sample = counts3d.sum(axis=1)

    counts_per_locus = counts_per_sample.sum(axis=1)[:, np.newaxis]

    mixture = counts_per_sample/counts_per_locus

    return mixture, counts_per_locus[:, 0] # to allow broadcasting outside this function

def weighted_mixture(data, subset=None, sample_space=None):
    """Calculate the finite mixture distribution.

    This function is used to calculate the empirical finite mixture
    distribution from count data over a defined sample space. It can
    be applied to the whole set or a subset of the population. The
    computation of the mixture is done for each unit of observation
    given by the size of the first dimension of `data`.

    Parameters
    ----------
    data : pandas.DataFrame object
        Input data frame with a set of counts of observations for each
        state and sample.
    samples : list of str, optional
        Sample IDs of subpopulation if applicable, defaults to
        None. Must match the sample IDs used in `data`.
    states : list of str
        IDs of elements of the hypothesis space to which `data`
        applies. Must match IDs used within column names of `data`.

    Returns
    -------
    mix : numpy.ndarray
        The mixture distribution as a 2-d array. The shape is
        according to `data` in the first axis (i.e. the number of
        calculated distributions) and `states` in the second axis.
    cov : numpy.ndarray
        The sum of counts for each mixture distributions as a 1-d
        array.

    Notes
    -----
    The weight of each individual in the mixture equals the sum of
    counts measured in the individual divided by the sum of counts in
    all individuals.

    """
