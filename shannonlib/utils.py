# utils.py
"""Utility functions.
"""

import subprocess

import numpy as np


def supremum_numsites(tabixfiles, chrom):
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


def supremum_position(tabixfiles, chrom):
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


def groupname(by=None, name=None, fname=None):
    """Return filename of the subgroup.
    """

    # ensure that name is a tuple of strings
    name = tuple(str(key) for key in name)

    group = '_and_'.join('_'.join(items) for items in zip(by, name))

    old_suffix = fname.split('.')[-1]

    new_suffix = '.'.join([group, old_suffix])

    return fname.replace(old_suffix, new_suffix)
