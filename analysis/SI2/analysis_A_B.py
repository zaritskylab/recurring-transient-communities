import math
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from common.constants import SI_FIGURES_LOCATION, PROCESSED_TIMESERIES_LOCATION
from data_layer.cross_correlation import xcorr_w_lags
from data_layer.preprocess_data import _get_cells_data_from_csv


def plot_wt_pairs_corr_by_dist(distances, corrs):
    def log_func(x, a, b):
        return a + b * np.log(x)

    def _fit_goodness(y, y_fit):
        # residual sum of squares
        ss_res = np.sum((y - y_fit) ** 2)

        # total sum of squares
        ss_tot = np.sum((y - np.mean(y)) ** 2)

        # r-squared
        r2 = 1 - (ss_res / ss_tot)
        return r2

    fig, ax = plt.subplots(figsize=(14, 10))

    x_start = int(math.floor(min(distances)))
    x_end = int(min(math.ceil(max(distances)), 155))

    x_ticks_major = np.arange(0, 151, 50)
    plt.scatter(distances, corrs, s=20)

    popt_log, pcov_log = curve_fit(log_func, distances, corrs, maxfev=1000)
    log_r2 = _fit_goodness(np.array(corrs), log_func(np.array(distances), *popt_log))
    plt.plot(np.linspace(x_start, x_end, x_end - x_start),
             log_func(np.linspace(x_start, x_end, x_end - x_start), *popt_log), c='red',
             label=f'Log (R2 = {round(log_r2, 2)})')

    plt.xticks(fontsize=16)
    plt.yticks([0, 0.5, 1], fontsize=16, rotation=90)
    plt.xlabel("Pair's Distance (\u03BCm)", fontsize=22)
    plt.ylabel("Pair's Correlation Score", fontsize=22)
    plt.ylim(0, 1)
    plt.xlim(0, x_end)
    ax.axvline(x=14, color='grey', linestyle='--', linewidth=2)
    plt.text(15, max(corrs) / 1.3, 'Close cell pair (\u2264 14 \u03BCm)', rotation=90, fontsize=15)

    ax.set_xticks(x_ticks_major, major=True)
    ax.tick_params(width=2, length=10, which='major')
    plt.tight_layout()
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI2_A.png')
    plt.close(fig=fig)

    fig, ax = plt.subplots(figsize=(14, 10))

    plt.plot(np.linspace(x_start, x_end, x_end - x_start),
             np.gradient(log_func(np.linspace(x_start, x_end, x_end - x_start), *popt_log)) / np.gradient(
                 np.array(np.linspace(x_start, x_end, x_end - x_start))), c='green')

    plt.xticks(fontsize=16)
    plt.yticks([-0.1, -0.05, 0], fontsize=16, rotation=90)
    plt.xlabel("Pair's Distance (\u03BCm)", fontsize=22)
    plt.ylabel("Slope of Correlation Function Fit (Derivative)", fontsize=22)
    plt.ylim(-0.1, 0.01)
    plt.xlim(0, x_end)
    ax.axvline(x=7, color='grey', linestyle='--', linewidth=2)
    plt.text(8, -0.05, '7 \u03BCm', rotation=90, fontsize=15)
    ax.axvline(x=14, color='grey', linestyle='--', linewidth=2)
    plt.text(15, -0.05, '14 \u03BCm', rotation=90, fontsize=15)

    ax.set_xticks(x_ticks_major, major=True)
    ax.tick_params(width=2, length=10, which='major')


    plt.tight_layout()
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI2_B.png')
    plt.close(fig=fig)

def main():
    all_scores = []
    all_dists = []

    folders = os.listdir(f'{PROCESSED_TIMESERIES_LOCATION}/mid_third_wild_type/')

    for folder in folders:
        num_of_cells = int(len(os.listdir(f'{PROCESSED_TIMESERIES_LOCATION}/mid_third_wild_type/{folder}/'))) - 1
        cells_data = _get_cells_data_from_csv(f'{PROCESSED_TIMESERIES_LOCATION}/mid_third_wild_type/{folder}/', num_of_cells)
        time_series_len = len(cells_data[0])

        pairs_above_correlation_cutoff_and_pval_cutoff = \
            xcorr_w_lags(cells_data, num_of_cells, time_series_len, 0, correlation_cutoff=0)
        dist = [pair_dict['dist'] for _, pair_dict in
                pairs_above_correlation_cutoff_and_pval_cutoff.items()]
        all_dists.extend(dist)
        score = [pair_dict['score'] for _, pair_dict in
                 pairs_above_correlation_cutoff_and_pval_cutoff.items()]
        all_scores.extend(score)

    plot_wt_pairs_corr_by_dist(all_dists, all_scores)
