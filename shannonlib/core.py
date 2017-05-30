# -*- coding:utf-8 -*-
# core.py

"""Core functions.
"""

import collections
import csv
import functools
import glob
import io
import math
import os
import subprocess
import sys

import numexpr as ne
import numpy as np
import pandas as pd

import shannonlib.constants as constant

from .estimators import shannon_entropy
from .preprocessing import impute, groupname
from .gpf_utils import get_regions, get_data


def divergence(pop=None, chrom=None, filename=None, imputation=False,
               method='jsd', filetype='bismark_coverage'):
    """Computes within-group divergence for population.
    """

    # parameters
    min_samplesize = 2
    min_count = 3
    load = 1e3

    input_files = pop['url']
    labels = pop['label']
    data_columns = [[(4, 'count_mC', np.uint8), (5, 'count_C', np.uint8)]]

    regions_pct, regions = get_regions(input_files, chrom=chrom, exp_numsites=load)

    regions_data = get_data(input_files, labels=labels, data_columns=data_columns, regions=regions)

    for progress, indata in zip(regions_pct, regions_data):

        if indata.empty:
            print('...{:>5} % (skipped empty region)'.format(progress))
            continue

        if imputation:
            impute_value = impute(indata, method='pseudocount')
            indata.fillna(impute_value, inplace=True)

        # preprocess data
        count_per_unit = indata.sum(axis=1, level='sample')
        samplesize = count_per_unit.notnull().sum(axis=1)
        count_filter = (count_per_unit >= min_count).any(axis=1)
        samplesize_filter = (samplesize >= min_samplesize)
        combined_filter = (count_filter & samplesize_filter)

        data = indata[combined_filter]

        if not data.empty:
            # extract data
            data_unit = count_per_unit[combined_filter]
            data_feature = (data.sum(axis=1, level='feature').astype(np.int32))
            counts = data.values.reshape(
                data_unit.shape[0],
                data_unit.shape[1],
                data_feature.shape[1])
            # estimate JSD
            mix_entropy = shannon_entropy(data_feature.values)
            avg_entropy = np.average(
                shannon_entropy(counts, axis=2),
                weights=data_unit.fillna(0),
                axis=1)
            # output
            div = data_feature
            div.insert(0, 'JSD_bit_', constant.LOG2E *
                       (mix_entropy - avg_entropy))
            div.insert(1, 'sample size', samplesize[combined_filter])
            div.insert(2, 'HMIX_bit_', constant.LOG2E * mix_entropy)
            # div['members'] = (data_unit.apply(lambda x: ','.join(x.dropna().index), axis=1))

            if not os.path.isfile(filename):
                header = True
            elif os.stat(filename).st_size == 0:
                header = True
            else:
                header = False

            (div
             .round({'JSD_bit_': 3, 'HMIX_bit_': 3})
             .to_csv(filename, header=header, sep='\t', index=True, mode='a'))
        else:
            print('...{:>5} % (skipped region not meeting QC)'.format(progress))
            continue

        print('...{:>5} %'.format(progress))

    return None
