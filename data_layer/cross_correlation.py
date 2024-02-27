import math
from typing import List

import numpy as np
import pandas as pd
import scipy.stats

from common.euclidean_distance import EuclideanMatrixSingleton


# def xcorr_w_lags(self, cells_data: List[pd.DataFrame], cell_num: int, time_series_len: int, lags: int,
#                  cutoff: float = -1.0, p_val_cutoff=0.05):
#     """
#     Calculate CrossCorrelation of all unique pairs in cells_data.
#     Return:
#         1. corrs_passed_cutoff_by_lag - dictionary of {lag: [corrs of pairs passed the cutoff]}
#         2. out - num of pairs didn't pass the cutoff
#         3. pairs_edges - list of pairs with correlation score > cutoff (used for graph analysis)
#         4. total_pairs_passed_cutoff_by_rounded_distance - for each natural number representing an euclidian distance,
#                                                            returns the total num of pairs passed the cutoff.
#         5. pairs_passed_cutoff_max_corr_and_distance_and_lags - for each pair, returns the following structure:
#                                                              {'dist': pair-distance, 'lag': best-lag, 'score': max-corr}
#     :param cells_data:
#     :param cell_num:
#     :param time_series_len:
#     :param lags:
#     :param cutoff:
#     :return:
#     """
#     out = 0
#     cell_pair = 0
#
#     distance_matrix = self.euclidean_matrix.get_distance_matrix()
#     rounded_distance_matrix: np.ndarray = distance_matrix.round()
#     distinct_rounded_distances = np.unique(rounded_distance_matrix)  # [1:] # drop 0 distance value
#
#     max_corrs_by_lag = {i: 0 for i in range(-lags, lags + 1, 1)}
#     corrs_passed_cutoff_by_lag = {i: [] for i in range(-lags, lags + 1, 1)}
#     total_pairs_passed_cutoff_by_rounded_distance = {dist: 0 for dist in distinct_rounded_distances}
#     pairs_passed_cutoff_max_corr_and_distance_and_lags = {}
#     pairs_corr_at_fixed_lag_above_p_val = {}
#     pairs_edges = []
#
#     for i in range(cell_num):
#         for j in range(i + 1, cell_num):
#             cell_i = cells_data[i]['normalized_intensity']
#             cell_j = cells_data[j]['normalized_intensity']
#
#             y1 = cell_i.values
#             y2 = cell_j.values
#
#             p_val_by_lag = {}
#             pearson_by_lag = {}
#             for lag in max_corrs_by_lag.keys():
#                 if lag < 0:
#                     tmp_y1 = y1[int(math.fabs(lag)):]
#                     tmp_y2 = y2[0:time_series_len + lag]
#                 elif lag > 0:
#                     tmp_y1 = y1[0: time_series_len - lag]
#                     tmp_y2 = y2[int(math.fabs(lag)):]
#                 else:
#                     tmp_y1 = y1
#                     tmp_y2 = y2
#                 pearson_by_lag[lag] = scipy.stats.pearsonr(tmp_y1, tmp_y2)[0]
#                 p_val_by_lag[lag] = scipy.stats.pearsonr(tmp_y1, tmp_y2)[1]
#
#             # Filtering out lags with correlation scores with insignificant p-values
#             if min(p_val_by_lag.values()) < 0.05:
#                 filtered_lags = dict(filter(lambda tup: tup[1] < 0.05,
#                                             p_val_by_lag.items()))  # TODO: change to p_val_by_lag to include also insignificant pairs
#             else:
#                 filtered_lags = dict(filter(lambda tup: tup[1] < p_val_cutoff,
#                                             p_val_by_lag.items()))  # TODO: change to p_val_by_lag to include also insignificant pairs
#             filtered_pearson_by_lag = dict(filter(lambda tup: tup[0] in (filtered_lags.keys()), pearson_by_lag.items()))
#
#             # For fixed-lags
#             if len(filtered_pearson_by_lag) > 0:
#                 best_significant_score_at_fixed_lag = max(filtered_pearson_by_lag.get(-lags, -9),
#                                                           filtered_pearson_by_lag.get(lags, -9))
#                 if best_significant_score_at_fixed_lag > -9:
#                     pairs_corr_at_fixed_lag_above_p_val[(i, j)] = {'dist': distance_matrix[i][j],
#                                                                    'lag': lags,
#                                                                    'score': best_significant_score_at_fixed_lag}
#
#             # For optimal-lags
#             if len(filtered_pearson_by_lag) > 0:
#                 limited_ccor = list(pearson_by_lag.values())
#                 curr_max_corr = max(filtered_pearson_by_lag.values())
#                 lag_with_curr_max_corr = list(limited_ccor).index(curr_max_corr) - lags
#                 if curr_max_corr > cutoff:
#                     max_corr = curr_max_corr
#                     lag_with_max_corr = lag_with_curr_max_corr
#
#                     current_pair_original_distance = distance_matrix[i][j]
#                     current_pair_rounded_distance = rounded_distance_matrix[i][j]
#                     total_pairs_passed_cutoff_by_rounded_distance[current_pair_rounded_distance] += 1
#                     pairs_passed_cutoff_max_corr_and_distance_and_lags[(i, j)] = {
#                         'dist': current_pair_original_distance, 'lag': lag_with_max_corr, 'score': max_corr,
#                         'pval': p_val_by_lag[lag_with_max_corr]}
#
#                     corrs_passed_cutoff_by_lag[lag_with_max_corr].append(max_corr)
#                     max_corrs_by_lag[lag_with_max_corr] = max_corrs_by_lag[lag_with_max_corr] + 1
#                     pairs_edges.append((i, j, curr_max_corr))
#                 else:
#                     out += 1
#                 cell_pair += 1
#
#     return corrs_passed_cutoff_by_lag, \
#         out, \
#         pairs_edges, \
#         total_pairs_passed_cutoff_by_rounded_distance, \
#         pairs_passed_cutoff_max_corr_and_distance_and_lags, \
#         pairs_corr_at_fixed_lag_above_p_val



def xcorr_w_lags(cells_data: List[pd.DataFrame], cell_num: int, time_series_len: int, lags: int,
                 correlation_cutoff: float = -1.0, p_val_cutoff=0.05):
    """ Calculate CrossCorrelation of all unique pairs in cells_data."""
    cell_pair = 0

    distance_matrix = EuclideanMatrixSingleton(cells_data).get_distance_matrix()
    max_corrs_by_lag = {i: 0 for i in range(-lags, lags + 1, 1)}
    pairs_above_correlation_cutoff_and_pval_cutoff = {}

    for i in range(cell_num):
        for j in range(i + 1, cell_num):
            cell_i = cells_data[i]['normalized_intensity']
            cell_j = cells_data[j]['normalized_intensity']

            y1 = cell_i.values
            y2 = cell_j.values

            # Calculate Pearson Correlation Coefficient for each lag
            p_val_by_lag = {}
            pearson_by_lag = {}
            for lag in max_corrs_by_lag.keys():
                if lag < 0:
                    tmp_y1 = y1[int(math.fabs(lag)):]
                    tmp_y2 = y2[0:time_series_len + lag]
                elif lag > 0:
                    tmp_y1 = y1[0: time_series_len - lag]
                    tmp_y2 = y2[int(math.fabs(lag)):]
                else:
                    tmp_y1 = y1
                    tmp_y2 = y2
                pearson_by_lag[lag] = scipy.stats.pearsonr(tmp_y1, tmp_y2)[0]
                p_val_by_lag[lag] = scipy.stats.pearsonr(tmp_y1, tmp_y2)[1]

            # Filtering out lags with correlation scores with insignificant p-values
            filtered_lags = dict(filter(lambda tup: tup[1] < p_val_cutoff, p_val_by_lag.items()))
            filtered_pearson_by_lag = dict(filter(lambda tup: tup[0] in (filtered_lags.keys()), pearson_by_lag.items()))

            # Find the lag with the highest correlation score that passed the cutoff
            if len(filtered_pearson_by_lag) > 0:
                limited_ccor = list(pearson_by_lag.values())
                curr_max_corr = max(filtered_pearson_by_lag.values())
                lag_with_curr_max_corr = list(limited_ccor).index(curr_max_corr) - lags
                if curr_max_corr > correlation_cutoff:
                    max_corr = curr_max_corr
                    lag_with_max_corr = lag_with_curr_max_corr

                    current_pair_original_distance = distance_matrix[i][j]
                    pairs_above_correlation_cutoff_and_pval_cutoff[(i, j)] = {
                        'dist': current_pair_original_distance, 'lag': lag_with_max_corr, 'score': max_corr,
                        'pval': p_val_by_lag[lag_with_max_corr]}

                    max_corrs_by_lag[lag_with_max_corr] = max_corrs_by_lag[lag_with_max_corr] + 1
                cell_pair += 1

    return pairs_above_correlation_cutoff_and_pval_cutoff