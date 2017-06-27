# -*- coding:utf-8 -*-
# # estimators.py

"""Estimators of information-theoretic quantities.
"""

from functools import partial
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
                                 # FIXME do fillna at the end!
                                 # ZeroDivisionError possible if all weights
                                 # are zero!
                                 weights=data_unit.fillna(0),
                                 axis=1)

        jsd = np.where(np.isclose(mix_entropy, avg_entropy, atol=1.e-2),
                       .0, mix_entropy - avg_entropy)

        div = data_feature
        div['JSD (bit)'] = constant.LOG2E * jsd
        div['sample size'] = samplesize[combined_filter]
        # div.insert(2, 'HMIX (bit)', constant.LOG2E * mix_entropy)
        # div['members'] = (data_unit.apply(lambda x: ','.join(x.dropna().index), axis=1))

        return div


def jsd_st(data):
    """Hierarchical JS divergence.
    """

    divIT = jsd_is(data)

    if divIT.empty:
        return False
    elif data.columns.nlevels > 2:
        hierarchy = data.columns.names[:-2]
        data_groups = [data.groupby(axis=1, level=hierarchy[:1 + i])
                       for i, _ in enumerate(hierarchy)]
        processes = min(len(hierarchy), cpu_count())
        with Pool(processes) as pool:
            divIS = pool.map(partial(jsd_is, grouped=True), data_groups)

        return divIS
    else:
        return divIT


def jsd_st_bottomup(parameter_list):
    pass