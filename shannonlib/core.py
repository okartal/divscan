# -*- coding:utf-8 -*-
# core.py

"""Core functions.
"""

import os

import numpy as np
import pandas as pd

import shannonlib.estimators as est
import shannonlib.gpf_utils as gpf


def divergence(sample, chrom=None, data_columns=None, outfile=None,
               chunksize=None, hierarchy=None):
    """Output genome-wide divergence for a population.
    """

    regions_pct, regions = gpf.get_regions(
        sample['url'], chrom=chrom, exp_numsites=chunksize)

    regions_data = gpf.get_data(sample, hierarchy, data_columns, regions)

    for progress, data in zip(regions_pct, regions_data):

        if data.empty:
            print('...{:>5} % (skipped empty region)'.format(progress))
            continue
        
        div = est.js_divergence(data)

        if not div:
            continue
        
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
