"""Core functions.
"""

import collections
import csv
import functools
import glob
import io
import logging
import math
import os
import subprocess
import sys

import numexpr as ne
import numpy as np
import pandas as pd

logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s",
                    level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

LOG2E = np.log2(math.e)


def divergence(pop=None, chrom=None, filename=None, imputation=False,
               method='jsd', filetype='bismark_coverage'):
    """Computes within-group divergence for population.
    """

    # parameters
    # load is measured in expected number of sites in the queried region
    min_sample_size = 2
    min_count = 3
    load = 1e5

    input_files = pop['url']

    labels = pop['label']

    column = collections.OrderedDict()

    if filetype == 'bismark_coverage':
        column[0] = '#chrom'
        column[1] = 'start'
        column[2] = 'end'
        column[4] = 'E'  # my single letter code for methylated C
        column[5] = 'C'
        types = [str, np.int64, np.int64, np.int32, np.int32]
        dtype = dict(zip(column.values(), types))
        alphabet = ['E', 'C']
        coordinate = [column[i] for i in range(3)]

    sup_chrom = chrom_supremum(input_files, chrom)
    if sup_chrom == None:
        logging.info("Skipping because chromosome is missing.")
        return False

    sup_nsites = nsites_supremum(input_files, chrom)

    if sup_nsites == None or sup_nsites == 0:
        logging.info("Skipping because there are no entries.")
        return False

    stepsize = math.ceil(sup_chrom / sup_nsites * load)

    if stepsize < sup_chrom:
        step = stepsize
        print("step size:", step)
    else:
        step = sup_chrom
        print("step size:", step, "(max. for contig {0})".format(chrom))

    pos_start = list(range(0, sup_chrom, step + 1))
    pos_end = list(range(step, sup_chrom, step + 1)) + [sup_chrom]
    regions = zip(pos_start, pos_end)

    # Print step size and start the calculation for each interval

    for interval in regions:

        # INPUT

        region = chrom + ':{0}-{1}'.format(*interval)

        tabix_query = (subprocess.Popen(['tabix', f, region],
                                        stdout=subprocess.PIPE,
                                        universal_newlines=True)
                       for f in input_files)

        dataframes = (pd.read_table(query.stdout, comment='#', header=None,
                                    usecols=list(column.keys()),
                                    names=list(column.values()), dtype=dtype,
                                    index_col=[0, 1, 2])
                      for query in tabix_query)

        indata = pd.concat(dataframes, axis=1, keys=labels,
                           names=['sample', 'feature'])

        # PROCESS
        print('processing', region, '...')

        if imputation:
            impute_value = impute(indata, method='pseudocount')
            indata.fillna(impute_value, inplace=True)

        count_per_unit = indata.sum(axis=1, level='sample')

        sample_size = count_per_unit.notnull().sum(axis=1)

        # filter
        pass_count = (count_per_unit >= min_count).any(axis=1)
        pass_sample_size = (sample_size >= min_sample_size)
        pass_filter = (pass_count & pass_sample_size)

        data = indata[pass_filter]

        if not data.empty:
            # data transformations
            data_unit = count_per_unit[pass_filter]
            data_feature = (data.sum(axis=1, level='feature').astype(np.int32))
            counts = data.values.reshape(
                data_unit.shape[0],
                data_unit.shape[1],
                data_feature.shape[1])
            # JSD terms
            mix_entropy = entropy(data_feature.values)
            avg_entropy = np.average(
                entropy(counts, axis=2), weights=data_unit.fillna(0), axis=1)
            # output
            div = data_feature
            div.insert(0, 'JSD_bit_', LOG2E * (mix_entropy - avg_entropy))
            div.insert(1, 'sample size', sample_size[pass_filter])
            div.insert(2, 'HMIX_bit_', LOG2E * mix_entropy)
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
            continue

    return None


def groupname(by=None, name=None, fname=None):
    """Return filename of the subgroup.
    """

    # ensure that name is a tuple of strings
    name = tuple(str(key) for key in name)

    group = '_and_'.join('_'.join(items) for items in zip(by, name))

    old_suffix = fname.split('.')[-1]

    new_suffix = '.'.join([group, old_suffix])

    return fname.replace(old_suffix, new_suffix)


def chrom_supremum(tabixfiles, chrom):
    """Return the least upper bound for the chrom end coordinate.
    """

    end_coordinate = list()

    for f in tabixfiles:
        tabix = subprocess.Popen(["tabix", f, chrom], stdout=subprocess.PIPE)
        tail = subprocess.Popen(
            ["tail", "-1"], stdin=tabix.stdout, stdout=subprocess.PIPE)
        cut = subprocess.Popen(
            ["cut", "-f3"], stdin=tail.stdout, stdout=subprocess.PIPE)
        # Allow first process to receive a SIGPIPE if process 2 exits.
        tabix.stdout.close()
        try:
            base_position = int(cut.communicate()[0])
            end_coordinate.append(base_position)
        except ValueError:
            continue

    try:
        out = np.max(end_coordinate)
    except ValueError:
        out = None

    return out


def nsites_supremum(tabixfiles, chrom):
    '''Return the least upper bound for the number of covered sites.
    '''

    sites = list()

    for f in tabixfiles:
        tabix = subprocess.Popen(["tabix", f, chrom], stdout=subprocess.PIPE)
        wcl = subprocess.Popen(
            ["wc", "-l"], stdin=tabix.stdout, stdout=subprocess.PIPE)
        tabix.stdout.close()  # Allow tabix to receive a SIGPIPE if wcl exits.
        try:
            site_count = int(wcl.communicate()[0])
            sites.append(site_count)
        except ValueError:
            continue

    try:
        out = np.max(sites)
    except ValueError:
        out = None

    return out


def entropy(counts, axis=1, method='mle'):
    """Entropy (in nat) of feature frequency profile.
    """

    expression = "sum(where(pr > 0, -pr * log(pr), 0), axis={})".format(axis)

    if method == 'mle':
        '''Maximum-likelihood/plug-in/non-parametric estimator
        '''
        sumCounts = counts.sum(axis)[..., np.newaxis]
        pr = counts / sumCounts
        result = ne.evaluate(expression)

    return result


def impute(data, method='pseudocount'):
    if method == 'pseudocount':
        # given by the parameters of a uninformative Dirichlet prior on the
        # probabilities
        value = 1
    return value
