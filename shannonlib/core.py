# -*- coding:utf-8 -*-
# core.py

"""Core functions.
"""

import os

import numpy as np
import pandas as pd

import shannonlib.estimators as est
import shannonlib.gpf_utils as gpf


def divergence(sample, chrom=None, data_columns=None, outfile=None, chunksize=None):
    """Computes within-group divergence for population.
    """

    regions_pct, regions = gpf.get_regions(
        sample['url'], chrom=chrom, exp_numsites=chunksize)

    regions_data = gpf.get_data(sample['url'], labels=sample['label'],
                                data_columns=data_columns, regions=regions)

    for progress, data in zip(regions_pct, regions_data):

        if data.empty:
            print('...{:>5} % (skipped empty region)'.format(progress))
            continue

        div = est.jsd_is(data)

        if div.empty:
            continue

        # div_subgroups = gpf.groupby('stage', metadata=sample, data=data)
        # note: use multiprocessing if you have a list of groupby objects
        # div_is = div_subgroups.apply(est.jsd_is)
        # js = div_is.xs('JSD_bit_', level='feature', axis=1)
        # ss = div_is.xs('sample size', level='feature', axis=1)
        # avg = np.average(js.values, weights=ss.values, axis=1)
        # div.insert(1, 'JSD_is', avg)

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
