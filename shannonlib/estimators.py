# -*- coding:utf-8 -*-
# estimators.py

"""Estimators of information-theoretic quantities.
"""

from functools import partial
from itertools import compress
from multiprocessing import Pool, cpu_count

import numexpr as ne
import numpy as np
import pandas as pd

import shannonlib.constants as constant
from shannonlib.preprocessing import partition


def shannon_entropy(countmatrix, axis=1, method='plug-in'):
    """Shannon entropy (in nat) of the feature frequency profile.
    """

    if method == 'plug-in':
        expression = ("sum(where(prob > 0, -prob * log(prob), 0), axis={})"
                      .format(axis))
        count_distribution = countmatrix.sum(axis)[..., np.newaxis]
        prob = countmatrix / count_distribution
        return ne.evaluate(expression)


def jsd_is(data, grouped=False):
    """
    """
    # BUG in JSD calculation!!!

    # if imputation:
    #     impute_value = impute(data, method='pseudocount')
    #     data.fillna(impute_value, inplace=True)

    if grouped:
        return data.apply(jsd_is)

    count_per_unit = data.sum(axis=1, level='sampling_unit')
    samplesize = count_per_unit.notnull().sum(axis=1)

    # QC filter
    min_samplesize = 1
    min_count = 1
    count_filter = (count_per_unit >= min_count).any(axis=1)
    samplesize_filter = (samplesize >= min_samplesize)
    combined_filter = (count_filter & samplesize_filter)
    data = data[combined_filter]

    if data.empty:
        print('...{:>5} % (skipped low-quality region)'.format(progress))
        return data  # continue
    else:
        data_unit = count_per_unit[combined_filter]
        data_feature = (data.sum(axis=1, level='feature').astype(np.int32))
        counts = data.values.reshape(
            data_unit.shape[0],
            data_unit.shape[1],
            data_feature.shape[1])

        mix_entropy = shannon_entropy(data_feature.values)
        avg_entropy = np.average(shannon_entropy(counts, axis=2),
                                 weights=data_unit.fillna(0),
                                 axis=1)
        # FIXME: do fillna at the end? ZeroDivisionError
        # possible if all weights become zero but that
        # only happens if all weights ar NaN which
        # should not happen if data is preprocessed
        # properly

        jsd = np.where(np.isclose(mix_entropy, avg_entropy, atol=1.e-2),
                       .0, mix_entropy - avg_entropy)
        # FIXME
        # if PDF binary and ratio (here meth. level) of all counts summed is close to 0 (or 1) set JSD=0
        # elif entropy_of_mixture close to zero set JSD=0
        # else set JSD = entropy_of_mixture - avg_of_entropy

        div = data_feature
        div['JSD (bit)'] = constant.LOG2E * jsd
        div['sample size'] = samplesize[combined_filter]
        # div.insert(2, 'HMIX (bit)', constant.LOG2E * mix_entropy)
        # div['members'] = (data_unit.apply(lambda x: ','.join(x.dropna().index), axis=1))

        return div


def jsd_average(divset, mask):
    """Returns average JS divergence over a subset given by mask.
    """

    subset = list(compress(divset, mask))

    if len(subset) == 1:
        return subset[0]['JSD (bit)']
    elif len(subset) > 1:
        div = pd.concat(subset, keys=None, axis=1)
        jsdiv = div.xs('JSD (bit)', axis=1)
        ssize = div.xs('sample size', axis=1).fillna(0)
        avg = np.average(jsdiv.values, weights=ssize.values, axis=1)
        return pd.Series(avg, index=div.index)
    # else return an exception for empty subset and let caller handle it


def js_divergence(data, meta, subdivision=None, weight=None):
    """Hierarchical Jensen-Shannon Divergence.

    Computes hierarchical Jensen-Shannon Divergence (hJSD) over a statistical
    ensemble at each index.

    Parameters
    ----------
    data : pandas.DataFrame
        Data values for population sample. The DataFrame is assumed to have a
        2-level index along axis 1. The first level indexes the population
        sample, the second level indexes the sample space. The values are counts
        or probabilities for each event.
    meta : pandas.DataFrame
        Metadata for the population sample. The index (along axis 0) must
        correspond with the population sample index in the data.
    subdivision : list of string objects
        Names of factors along which the population is subdivided. The order of
        factors matters as it defines the sampling hierarchy. The names must
        appear as columns in the metadata.
    weight : array-like
        Weights for the individual samples. If None (default), input is assumed
        to be count data. A value implies that input data is given in the form
        of probabilities (i.e. a probability distribution for each sample). The
        weight matrix must be of the same shape as data along axis 0 and axis 1,
        level 0 (sample). 

    Returns
    -------
    divergence : pandas.DataFrame
        Without subdivision, the function returns a single value (the classical
        JSD) for each index (axis 0). For n subdivisions, the function returns
        n+1 values for each index: the n between-group divergences at each level
        of subdivision plus the averaged within-group divergence for the
        terminal (i.e. non-subdivided) groups in the hierarchy.

    Notes
    -----
    The Jensen-Shannon Divergence (JSD) is a measure of heterogeneity for a set
    of probability distributions. This function implements a measure that
    generalizes JSD to the case of a hierarchically subdivided population
    (hJSD); JSD is the special case of hJSD without subdivision.

    TODO
    ----
    - remove loci without divergence for downstream nodes
    - diff = lambda x: [x[i] - x[i+1] for i,n in enumerate(x) if i < len(x)-1]
    """

    # preprocess
    if subdivision is None:
        subdivision = []

    # initialize key list
    template = 'J[{reference_id}:{divisor_id}{divisor_name}] (bit)'
    keys = [template.format(reference_id='0', divisor_id='i', divisor_name='')]

    if subdivision:
        quset = range(len(subdivision))
        for n, level in enumerate(subdivision):
            depth = n + 1
            subdivided = depth < len(subdivision)
            # add key
            key = template.format(
                reference_id=str(n),
                divisor_id=str(depth) if subdivided else 'i',
                divisor_name='({})'.format(level) if subdivided else '')
            keys.append(key)
            # add quotient set
            eqrel = subdivision[:depth]
            groups = meta.groupby(eqrel).groups.values()
            quset[n] = set(tuple(sorted(g)) for g in groups)  # if len(s) > 1 ?

        # Eliminate duplicate equivalence classes over all levels and keep track
        # at which level each unique class appears by using a boolean mask. This
        # way qusets can be recovered.
        uniq_classes = list(set.union(*quset))
        mask = [[1 if eclass in quset[n] else 0 for eclass in uniq_classes]
                for n, _ in enumerate(subdivision)]

        entropy = range(len(subdivision) + 2)

    # TODO for all sets in the terminal quset (leaves) do
    # ensuring proper merge (concat) along axis 0 at construction of the entropy data frame

    i_mixture = data.sum(axis=1, level='event')
    i_weights = data.sum(axis=1, level='sample')
    i_samples = i_weights.notnull().sum(axis=1)
    i_entropy = shannon_entropy(data.values.reshape(None), axis=None)
    
    g_mixture = np.average() # over all i!

    entropy[first] = shannon_entropy(g_mixture)
    entropy[last] = entropy[first] if entropy[first] close to 0 else np.average(i_entropy)

    entropy_actual = np.where(np.isclose(entropy[-2], 0.0), entropy[-2], VALUE)
    entropy.insert(-1, )

    # process entropy data
    mixture_entropy = entropy.iloc[:, :n].values
    entropy_mixture = entropy.iloc[:, 1:].values
    jsd = jsd_shannon(mixture_entropy, entropy_mixture)

    if subdivision:
        root_mixture_entropy = entropy.iloc[:, 0].values
        leaf_entropy_mixture = entropy.iloc[:, -1].values
        jsd_total = jsd_shannon(root_mixture_entropy, leaf_entropy_mixture)
        jsd = np.hstack([jsd_total, jsd])

    return pd.DataFrame(jsd, index=data.index, columns=keys)


def jsd_shannon(mixture_entropy, entropy_mixture):
    """Jensen-Shannon divergence from Shannon entropies.

    Jensen-Shannon divergence is calculated using pre-estimated Shannon
    entropies, mixture_entropy - entropy_mixture. Takes into account that
    mixture_entropy >= entropy_mixture must hold in order to avoid subtracting
    values close to zero.

    Parameters
    ----------
    mixture_entropy : array-like
    entropy_mixture : array-like

    Returns
    -------
    js_divergence : array-like

    TODO
    ----
    - return exception if arrays don't broadcast
    """

    a = mixture_entropy
    b = entropy_mixture
    return np.where(np.isclose(a, b), a, a - b)

    ###############
    # previous code
    ###############

    if not qusetunion:
        return jsd_is(data)
    else:
        dataset = [data[list(s)] for s in qusetunion]

        mask = [[1 if s in quset[level] else 0 for s in qusetunion]
                for level in subdivision]

        nproc = min(len(qusetunion), cpu_count())
        with Pool(processes=nproc) as pool:
            jsd_dataset = pool.map(jsd_is, dataset)

        nproc = min(len(mask), cpu_count())
        with Pool(processes=nproc) as pool:
            args = ((jsd_dataset, m) for m in mask)
            div = pool.starmap(jsd_average, args)

        import pdb
        pdb.set_trace()
        df = pd.concat(div, keys=title, axis=1)

        return df

    #         elif depth == 1:
    #             div_average(div_quset)
    #         elif depth > 1:
    #             keys = quset.keys()
    #             dist_to_root = depth
    #             while dist_to_root > 0:
    #                 if dist_to_root > 1:
    #                     groups = [list(g) for k, g in groupby(
    #                         sorted(keys, key=lambda x: x[-2]), lambda x: x[-2])]
    #                 elif dist_to_root == 1:
    #                     groups = keys
    #                 # apply a function to each group that averages or takes the DIVs
    #                 dist_to_root -= 1
    #                 keys =

    #         # averaging depth times based on parent nodes
    #         for k in range(depth):
    #             for each parent node do
    #             siblings = len() > 1
    #             if siblings:
    #                 d = average(d)

    #         # - set parent node based on previous quset; if depth==1, parent node is root
    #         # - group each parent node by current level; groups are child nodes
    #         # - exclude single-element children
    #         # - apply jsd_is to each child node
    #         # - if only one child, take jsd_is else: average over children = jsd of parent node
    #         # - do this 'depth' times until we get a single value per position (could also use a while loop although less explicit)
    #         # alternatives: pd.merge; a list of pd.Series followed by pd.concat and .fillna(0)
    #         div[title] =

    # return div
