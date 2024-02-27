import math
from typing import List

import numpy as np
import pandas as pd
import scipy.stats

from common.euclidean_distance import EuclideanMatrixSingleton


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