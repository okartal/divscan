# -*- coding:utf-8 -*-
# core.py

"""Core functions.
"""

import os

import numpy as np
import pandas as pd

import shannonlib.estimators as est
import shannonlib.gpf_utils as gpf


def divergence(metadata, chrom=None, data_columns=None, outfile=None,
               chunksize=None, hierarchy=None):
    """Output genome-wide divergence for a population.
    """

    regions_pct, regions = gpf.get_regions(
        metadata['url'],
        chrom=chrom,
        exp_numsites=chunksize)

    regions_data = gpf.get_data(
        metadata['url'],
        keys=metadata.index,
        data_columns=data_columns,
        regions=regions)

    for progress, data in zip(regions_pct, regions_data):

        if data.empty:
            print('...{:>5} % (skipped empty region)'.format(progress))
            continue

        div = est.js_divergence(data, metadata, hierarchy=hierarchy)

        if div.empty:
            continue

        if not os.path.isfile(outfile):
            header = True
        elif os.stat(outfile).st_size == 0:
            header = True
        else:
            header = False

        (div
        #  .round({'JSD (bit)': 3})
         .to_csv(outfile, header=header, sep='\t', index=True, mode='a'))

        print('...{:>5} %'.format(progress))

    return None
