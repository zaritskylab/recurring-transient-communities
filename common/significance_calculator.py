import itertools
from enum import Enum
from typing import List
import numpy as np
import pandas as pd
import scipy.stats


class SignificanceMethod(Enum):
    ZSCORE = 'ZSCORE'
    BOOTSTRAP = 'BOOTSTRAP'

class StatisticalTest(Enum):
    MANN_WHITNEY_U = 'MANN_WHITNEY_U'
    FISHERS_EXACT_TEST = 'FISHERS_EXACT_TEST'
    KRUSKAL_WALLIS = 'KRUSKAL_WALLIS'

def statistical_test_wrapper(df: pd.DataFrame, field: str, test_type: StatisticalTest, static_experiment_type: str = None):
    """Wrapper for mann whitney u test.
        If static_experiment_type is provided, then the test is performed between the static_experiment_type and all other experiment types.
        If static_experiment_type is not provided, then the test is performed between all experiment types."""
    results = {}
    if test_type == StatisticalTest.MANN_WHITNEY_U:
        stat_func = mann_whitney_u_test
    elif test_type == StatisticalTest.KRUSKAL_WALLIS:
        stat_func = kruskal_wallis_test
    else:
        stat_func = fishers_exact_test

    if static_experiment_type:
        experiment_type_1 = static_experiment_type
        experiment_types = df[df['experiment_type'] != experiment_type_1]['experiment_type'].unique()
        for experiment_type_2 in experiment_types:
            p_value = stat_func(df, experiment_type_1, experiment_type_2, field)
            results[(experiment_type_1, experiment_type_2)] = p_value

    else:
        experiment_types = df['experiment_type'].unique()
        for experiment_type_1, experiment_type_2 in itertools.combinations(experiment_types, 2):
            p_value = stat_func(df, experiment_type_1, experiment_type_2, field)
            results[(experiment_type_1, experiment_type_2)] = p_value

    return {k:round(v,5) for k,v in results.items()} #{k:pvalue_to_asterisks(v) for k,v in results.items()}

def fishers_exact_test(df: pd.DataFrame, experiment_type_1: str, experiment_type_2: str, field: str):
    experiment_1_data = df[df['experiment_type'] == experiment_type_1][field]
    experiment_2_data = df[df['experiment_type'] == experiment_type_2][field]

    exp_1_trues = (experiment_1_data <= 0.05).sum()
    exp_1_falses = (experiment_1_data > 0.05).sum()

    exp_2_trues = (experiment_2_data <= 0.05).sum()
    exp_2_falses = (experiment_2_data > 0.05).sum()

    contingency_table = np.array([[exp_1_trues, exp_2_trues], [exp_1_falses, exp_2_falses]])
    oddsratio, pvalue = scipy.stats.fisher_exact(contingency_table)

    return pvalue

def kruskal_wallis_test(df: pd.DataFrame, experiment_type_1: str, experiment_type_2: str, field: str):
    import scipy.stats as stats

    # define two groups with uneven sample sizes
    group1 = df[df['experiment_type'] == experiment_type_1][field]
    group2 = df[df['experiment_type'] == experiment_type_2][field]

    # perform Mann-Whitney U test
    statistic, p_value = stats.kruskal(group1, group2)
    # print the results
    return p_value


def mann_whitney_u_test(df: pd.DataFrame, experiment_type_1: str, experiment_type_2: str, field: str):
    # The Mann-Whitney U test is based on the ranks of the data, rather than the actual values. The test determines whether the two groups come from the same population or not by comparing the sum of ranks between the two groups. The test involves the following steps:
    #
    # Rank all the data from both groups together, from smallest to largest.
    # Calculate the sum of ranks for each group.
    # Calculate the U statistic using the following formula:
    # U = n1 x n2 + (n1 x (n1+1))/2 - R1
    #
    # where:
    #
    # U is the Mann-Whitney U statistic
    # n1 is the sample size of the first group
    # n2 is the sample size of the second group
    # R1 is the sum of ranks for the first group
    # Compare the calculated U value to a critical value in the Mann-Whitney U distribution table. If the calculated U value is less than or equal to the critical value, we can conclude that there is no significant difference between the two groups. Otherwise, we reject the null hypothesis and conclude that there is a significant difference between the two groups.
    # The Mann-Whitney U test is commonly used in various fields such as medicine, biology, and social sciences to compare two independent groups, particularly when the data is non-normal or when the sample size is small.
    import scipy.stats as stats

    # define two groups with uneven sample sizes
    group1 = df[df['experiment_type'] == experiment_type_1][field]
    group2 = df[df['experiment_type'] == experiment_type_2][field]

    # perform Mann-Whitney U test
    statistic, p_value = stats.mannwhitneyu(group1, group2)
    # print the results
    return p_value

class SignificanceCalculator:

    def _calc_zscore(self, actual_value, shuffled_values):
        mean = np.mean(shuffled_values)
        std = np.std(shuffled_values)
        actual_zscore = (actual_value - mean) / std
        p_value = scipy.stats.norm.sf(abs(actual_zscore))
        return p_value

    def _calc_bootstrap(self, actual_value, shuffled_values, less_is_more):

        shuffled_values.sort(reverse=(not less_is_more))
        denominator = (len(shuffled_values) + 1)

        if actual_value in shuffled_values:
            bootstrap_score = (shuffled_values.index(actual_value) +
                               shuffled_values.count(actual_value) +
                               1) / denominator
        else:
            bootstrap_score = 1
            for ix, shuffled_val in enumerate(shuffled_values):
                if (not less_is_more) and actual_value > shuffled_val:
                    bootstrap_score = (ix + 1) / denominator
                    break
                elif less_is_more and actual_value < shuffled_val:
                    bootstrap_score = (ix + 1) / denominator
                    break
        return bootstrap_score

    def calculate(self, method: SignificanceMethod, actual_value, shuffled_values: List,
                  less_is_more=False):
        if method == SignificanceMethod.ZSCORE:
            return self._calc_zscore(actual_value, shuffled_values)
        elif method == SignificanceMethod.BOOTSTRAP:
            return self._calc_bootstrap(actual_value, shuffled_values, less_is_more)
        else:
            raise Exception('Unknown SignificanceMethod value')


def pvalue_to_asterisks(pvalue):
    if pvalue < 0.0001:
        return '****'
    if pvalue < 0.001:
        return '***'
    elif pvalue < 0.01:
        return '**'
    elif pvalue < 0.05:
        return '*'
    else:
        return 'ns'

def asterisks_to_pvalue(asterisks):
    if asterisks == '****':
        return 0.0001
    if asterisks == '***':
        return 0.001
    elif asterisks == '**':
        return 0.01
    elif asterisks == '*':
        return 0.05
    else:
        return 1