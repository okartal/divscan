# -*- coding:utf-8 -*-
# # estimators.py

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


def jsd_average(mask, divset=None):
    """Returns average JS divergence over a subset of superset.
    """

    subset = list(compress(divset, mask))

    if len(subset) == 1:
        return subset[0]['JSD (bit)']
    elif len(subset) > 1:
        div = pd.concat(subset, keys=None, axis=1)
        jsdiv = div.xs('JSD (bit)', level='feature', axis=1)
        ssize = div.xs('sample size', level='feature', axis=1).fillna(0)
        avg = np.average(jsdiv.values, weights=ssize.values, axis=1)
        return pd.DataFrame(avg, index=div.index)
    # else return an exception for empty subset and let caller handle it


def js_divergence(data, meta, hierarchy=None):
    """Returns within-group divergences along the specified sampling hierarchy.

    compute the quset at each level; note that the same set can be present
    at different levels, so do not compute div_is more than once. that
    means, if we use a pool of workers to apply jsd_is to the *union* of
    qusets, we need to keep track of where each set belongs to, in order
    to make the right averaging;
    1. make qusets (q_level1, q_level2, etc.)
    2. make list of union of qusets (q)
    3. make lookup table for q (essentially given by qusets)
    columns: index(int), level1(bool), level2(bool), etc
    if q[i] in q_level1 set level1==True, etc.
    4. compute divIS list by applying div_is to q
    5. compress/mask the divIS result list with the level column and average the set

    FIXME: index_col must be given for meta, add index_col parameter in CLI
    TODO: only keep loci with divergence and sample size > 1 for downstream nodes
    """

    div = jsd_is(data)

    if hierarchy is None:
        hierarchy = []

    if not hierarchy or div.empty:
        return div

    title = ['JSD (bit)']
    quset = dict()

    for i, level in enumerate(hierarchy):
        depth = i + 1
        eqrel = hierarchy[:depth]
        groups = meta.groupby(eqrel).groups.values()
        title.append('JSD[{d}~{l}] (bit)'.format(d=depth, l=level))
        quset[level] = set(tuple(sorted(s)) for s in groups if len(s) > 1)

    qusetunion = list(set.union(*quset.values()))

    if not qusetunion:
        return div
    else:
        quset_data = [data[list(s)] for s in qusetunion]

        mask = [[1 if s in quset[level] else 0 for s in qusetunion]
                for level in hierarchy]
        
        nproc = min(len(qusetunion), cpu_count())

        with Pool(processes=nproc) as pool:
            div_qusetunion = pool.map(jsd_is, quset_data)
            # div_level = pool.map(partial(jsd_average, divset=div_qusetunion), mask)

        import pdb; pdb.set_trace()
        return pd.concat([div] + div_level, keys=title, axis=1)

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
