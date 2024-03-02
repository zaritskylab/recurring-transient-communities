import os
import numpy as np
from matplotlib import pyplot as plt
from common.constants import FIGURES_LOCATION, EXPERIMENTS_LOCATION
from data_layer.cross_correlation import xcorr_w_lags
from data_layer.preprocess_data import _get_cells_data_from_csv


def plot(distances, corrs):
    fig, ax = plt.subplots(figsize=(14, 10))
    x_ticks_major = np.arange(0, 101, 50)

    plt.scatter(distances, corrs, s=20)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16, rotation=90)
    plt.xlabel("Pair's Distance (\u03BCm)", fontsize=22)
    plt.ylabel("Pair's Correlation Score", fontsize=22)
    plt.ylim(0, 1)

    ax.axvline(x=14, color='grey', linestyle='--', linewidth=2)
    plt.text(15, max(corrs) / 1.3, 'Close cell pair (\u2264 14 \u03BCm)', rotation=90, fontsize=15)

    ax.set_xticks(x_ticks_major, major=True)
    ax.tick_params(width=2, length=10, which='major')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}Fig_1_B.png')

    plt.close(fig=fig)



def pairs_correlation_by_distance(cells_data, total_cells, time_series_len):
    """ Plots the correlation score of pairs of cells by their distance. """
    pairs_above_correlation_cutoff_and_pval_cutoff = xcorr_w_lags(cells_data, total_cells, time_series_len, 0, correlation_cutoff=0)
    distances = [pair_dict['dist'] for _, pair_dict in pairs_above_correlation_cutoff_and_pval_cutoff.items()]
    scores = [pair_dict['score'] for _, pair_dict in pairs_above_correlation_cutoff_and_pval_cutoff.items()]
    plot(distances, scores)


def main():
    relevant_experiment_folder = f'{EXPERIMENTS_LOCATION}processed/wild_type/sample_7/'
    total_cells = len(list(filter(lambda file_name: file_name.startswith('cell_') and file_name.endswith('.csv'),
                                  os.listdir(relevant_experiment_folder)
                                  )))
    cells_data = _get_cells_data_from_csv(relevant_experiment_folder, total_cells)
    pairs_correlation_by_distance(cells_data, total_cells, len(cells_data[0]))
