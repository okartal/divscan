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


def check_arguments(args):
    """
    """

    if os.path.isfile(args.output) and not os.stat(args.output).st_size == 0:
        msg = "\nExecution stops! Output file exists and is not empty!\n"
        sys.exit(msg)


class InputMismatchError(Exception):
    pass


class MissingInputError(Exception):
    pass


def tbx_merge(files, labels=[], data_columns=[], regions=[], join='outer',
              preset='bed'):
    """Combines tabix-indexed genome position files.
    """

    # check input arguments

    if not labels:
        keys = ['sample_{}'.format(pos + 1) for pos, _ in enumerate(files)]
    elif len(labels) == len(files):
        keys = labels
    else:
        raise InputMismatchError('Number of files and labels must match!')

    if not data_columns:
        # FIXME: assume that all columns are taken and are of type float
        raise MissingInputError(
            'The list of data_colums must have at least one entry!')
    elif len(data_columns) == len(files):
        pass
    elif len(data_columns) == 1:
        data_columns = data_columns * len(files)
    else:
        raise InputMismatchError(
            'Either supply a single entry in data_columns or'
            'the number of entries must match the number of files!')

    if preset == 'bed':
        index = [
            (0, '#chrom', str),
            (1, 'start', np.int64),
            (2, 'end', np.int64)]
        index_col = [i[0] for i in index]
    if preset == 'gff':
        pass
    if preset == 'vcf':
        pass
    if preset == 'sam':
        pass

    # output columns
    names = ['sample', 'value']
    columns = [index + cols for cols in data_columns]

    for region in regions:
        query = '{0}:{1}-{2}'.format(*region)

        tabix = enumerate(
            subprocess.Popen(['tabix', file_, query], stdout=subprocess.PIPE,
                             universal_newlines=True)
            for file_ in files)

        # import pdb; pdb.set_trace()

        dframes = (
            pd.read_table(
                tbx.stdout, header=None, index_col=index_col, comment='#',
                usecols=[f[0] for f in columns[i]],
                names=[f[1] for f in columns[i]],
                dtype={f[1]: f[2] for f in columns[i]})
            for i, tbx in tabix)

        merged_dframe = pd.concat(
            dframes, axis=1, keys=keys, names=names, join=join)

        if not merged_dframe.empty:
            yield merged_dframe
