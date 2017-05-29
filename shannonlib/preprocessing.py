# preprocessing.py
"""Functions for preprocessing data.
"""


def impute(data, method='pseudocount'):
    """Imputes missing values of a data frame.
    """
    if method == 'pseudocount':
        # given by the parameters of a uninformative Dirichlet prior on the
        # probabilities
        return 1