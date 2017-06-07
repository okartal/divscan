# -*- coding:utf-8 -*-
# core.py

"""Core functions.
"""

import os

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

        div = est.js_divergence(data) # this is div_it!
        # div_subgroups = gpf.groupby('stage', metadata=sample, data=data)
        # div_is = div_subgroups.apply(est.js_divergence) # alternatively map js_divergence to groups using multiprocessing
        # div_st = div_it - div_is_wavg # use np.average to over subgroups using 'sample size' as weights, see https://stackoverflow.com/a/33054358/2136626

        if div.empty:
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
