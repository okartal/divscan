# -*- coding:utf-8 -*-
# preprocessing.py

"""Functions for preprocessing data.
"""


def partition(dframe, levels=None):
    """Map IDs to subsets at each consecutive level of hierarchy.
    """

    maps = []

    for i, k in enumerate(levels):
        maps.append(dframe[levels[:i + 1]]
                    .apply(lambda x: ' > '.join(x.astype(str)), axis=1)
                    .to_dict())

    return maps

    
def impute(data, method='pseudocount'):
    """Imputes missing values of a data frame.
    """
    if method == 'pseudocount':
        # given by the parameters of a uninformative Dirichlet prior on the
        # probabilities
        return 1

def groupname(by=None, name=None, fname=None):
    """Return filename of the subgroup.
    """

    # ensure that name is a tuple of strings
    name = tuple(str(key) for key in name)

    group = '_and_'.join('_'.join(items) for items in zip(by, name))

    old_suffix = fname.split('.')[-1]

    new_suffix = '.'.join([group, old_suffix])

    return fname.replace(old_suffix, new_suffix)