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
import shannonlib.gpf_utils as gpf

from .estimators import shannon_entropy


def divergence(sample, chrom=None, outfile=None, chunksize=1e3):
    """Computes within-group divergence for population.
    """

    data_columns = [[(4, 'count_mC', np.uint8), (5, 'count_C', np.uint8)]]

    regions_pct, regions = gpf.get_regions(
        sample['url'], chrom=chrom, exp_numsites=chunksize)

    regions_data = gpf.get_data(sample['url'], labels=sample['label'],
                                data_columns=data_columns, regions=regions)

    for progress, indata in zip(regions_pct, regions_data):

        if indata.empty:
            print('...{:>5} % (skipped empty region)'.format(progress))
            continue

        # if imputation:
        #     impute_value = impute(indata, method='pseudocount')
        #     indata.fillna(impute_value, inplace=True)

        count_per_unit = indata.sum(axis=1, level='sample')
        samplesize = count_per_unit.notnull().sum(axis=1)

        # QC filter
        min_samplesize = 2
        min_count = 3
        count_filter = (count_per_unit >= min_count).any(axis=1)
        samplesize_filter = (samplesize >= min_samplesize)
        combined_filter = (count_filter & samplesize_filter)
        data = indata[combined_filter]

        if data.empty:
            print('...{:>5} % (skipped region that failed QC)'.format(progress))
            continue

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

        if not os.path.isfile(outfile):
            header = True
        elif os.stat(outfile).st_size == 0:
            header = True
        else:
            header = False

        (div
            .round({'JSD_bit_': 3, 'HMIX_bit_': 3})
            .to_csv(outfile, header=header, sep='\t', index=True, mode='a'))

        print('...{:>5} %'.format(progress))

    return None
