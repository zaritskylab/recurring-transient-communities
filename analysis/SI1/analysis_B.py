import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from common.config import PAIR_DISTANCE_THRESHOLD
from common.constants import SI_FIGURES_LOCATION, PROCESSED_TIMESERIES_LOCATION
from data_layer.cross_correlation import xcorr_w_lags
from data_layer.preprocess_data import _get_cells_data_from_csv


def plot(series_a: pd.Series, series_a_name: str, series_b: pd.Series, series_b_name: str):
    fig, ax = plt.subplots(figsize=(10, 4))

    # getting data of the histogram
    count_a, bins_count_a = np.histogram(series_a, bins=1000)
    count_b, bins_count_b = np.histogram(series_b, bins=1000)

    # finding the PDF of the histogram using count values
    pdf_a = count_a / sum(count_a)
    pdf_b = count_b / sum(count_b)

    # using numpy np.cumsum to calculate the CDF
    # We can also find using the PDF values by looping and adding
    cdf_a = 1 - np.cumsum(pdf_a)
    cdf_b = 1 - np.cumsum(pdf_b)

    # plotting PDF and CDF
    plt.plot(bins_count_a[1:], cdf_a, label=series_a_name)
    plt.plot(bins_count_b[1:], cdf_b, label=series_b_name)

    # plt.grid(True, alpha=0.7, linestyle='--')
    x_min = ax.get_xticks()[0]
    x_max = ax.get_xticks()[-1]
    plt.xticks([round(i, 1) for i in np.arange(x_min, x_max + 0.1, 0.2)])

    # plt.title('KDE - ' + title, fontsize=18)
    plt.xlim(-0.5, 1)
    plt.ylim(0, 1)
    plt.xlabel('Cell pair correlation', fontsize=22)
    plt.ylabel('Probability', fontsize=22)

    plt.xticks(fontsize=16)
    plt.yticks([0, 0.5, 1], fontsize=16, rotation=90)

    plt.legend(fontsize=18)

    plt.tight_layout()
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI1_B.png')

    plt.close(fig=fig)


def correlation_density_function_by_distance_group(cells_data, total_cells, time_series_len, distance_threshold):
    pairs_above_correlation_cutoff_and_pval_cutoff = xcorr_w_lags(cells_data, total_cells, time_series_len, 0, correlation_cutoff=-1, p_val_cutoff=1)
    distances = [pair_dict['dist'] for _, pair_dict in pairs_above_correlation_cutoff_and_pval_cutoff.items()]
    scores = [pair_dict['score'] for _, pair_dict in pairs_above_correlation_cutoff_and_pval_cutoff.items()]
    distance_groups = ['Near' if distance <= distance_threshold else 'Far' for distance in distances]

    df = pd.DataFrame(data={'distance': distances, 'score': scores, 'distance_group': distance_groups})

    near_scores = df[(df['distance_group'] == 'Near')]['score']
    far_scores = df[(df['distance_group'] == 'Far')]['score']

    plot(far_scores, 'Far Pairs', near_scores, 'Close Pairs')

def main():
    relevant_experiment_folder = f'{PROCESSED_TIMESERIES_LOCATION}/mid_third_wild_type/sample_8/'
    total_cells = len(list(filter(lambda file_name: file_name.startswith('cell_') and file_name.endswith('.csv'),
                                  os.listdir(relevant_experiment_folder)
                                  )))
    cells_data = _get_cells_data_from_csv(relevant_experiment_folder, total_cells)
    correlation_density_function_by_distance_group(cells_data, total_cells, len(cells_data[0]), distance_threshold=PAIR_DISTANCE_THRESHOLD)
