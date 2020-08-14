# -*- coding: utf-8 -*-
"""diversity.io
    ~~~~~~~~~~~~

    This module implements functions to load and transform data.

"""

import subprocess

import numpy as np
import pandas as pd


def bedcount_reader(bedcount, compression=None, chunksize=10000):
    """bedcount_reader returns a dataframe reader of the data."""
    reader = pd.read_table(bedcount, compression=compression,
                           chunksize=chunksize, header=0,
                           dtype={'#chrom': str, 'start': np.int})
    return reader


def population_filter(metadata, subset=None, relation=None):
    """Read metadata and return the population and the quotient set
    """

    pop = {'reference': None, 'qset': None}
    meta = pd.read_table(metadata, header=0)

    if subset is not None:
        pop['reference'] = list(meta.query(subset)['sample'])
    else:
        pop['reference'] = list(meta['sample'])

    if relation is not None:
        reference_meta = meta[meta['sample'].isin(pop['reference'])]
        group = reference_meta.groupby([relation])
        qset = []
        for _, df in group:
            qset.append(list(df['sample']))
        pop['qset'] = qset

    return pop

    # header_selection = '|'.join([s + '_' for s in population['sample']])
    # return dataframe.filter(regex=header_selection)


