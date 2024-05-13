import collections
import math
import os
from typing import List

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from common.constants import SI_FIGURES_LOCATION, PROCESSED_TIMESERIES_LOCATION, ARCOS_OUTPUT_LOCATION
from data_layer.cross_correlation import xcorr_w_lags
from data_layer.preprocess_data import _get_cells_data_from_csv
from shapely.geometry import Point, Polygon


def _calc_mean_density_in_pair_region(arcos_df, pair, r, all_points):
    def _get_circle(x, y, r):
        return Point(x, y).buffer(r)

    def _get_relevant_points(points: List[Point], poly: Polygon):
        return [ix for ix, point in enumerate(points) if poly.contains(point)]

    cell_1 = pair[0]
    cell_1_x = arcos_df[arcos_df['cell_id'] == cell_1]['x_microns'].mean()
    cell_1_y = arcos_df[arcos_df['cell_id'] == cell_1]['y_microns'].mean()
    cell_1_circle = _get_circle(cell_1_x, cell_1_y, r)
    cell_1_density = len(_get_relevant_points(all_points, cell_1_circle))

    cell_2 = pair[1]
    cell_2_x = arcos_df[arcos_df['cell_id'] == cell_2]['x_microns'].mean()
    cell_2_y = arcos_df[arcos_df['cell_id'] == cell_2]['y_microns'].mean()
    cell_2_circle = _get_circle(cell_2_x, cell_2_y, r)
    cell_2_density = len(_get_relevant_points(all_points, cell_2_circle))

    return (cell_1_density + cell_2_density) / 2


def plot(distances, densities, corrs):
    density_bins = [3.25, 5.5]

    binned_distance = collections.defaultdict(list)
    binned_score = collections.defaultdict(list)
    binned_density = collections.defaultdict(list)

    for i in range(len(distances)):
        for bin in density_bins:
            if distances[i] <= bin:
                binned_distance[bin].append(distances[i])
                binned_score[bin].append(distances[i])
                binned_density[bin].append(distances[i])
                break


    colors = ['tab:blue' if density <= 3.25 else 'tab:red' for density in densities]

    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(14, 10))


    x_ticks_major = np.arange(0, 101, 50)

    plt.scatter(distances, corrs, s=20, c=colors)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16, rotation=90)
    plt.xlabel("Pair's Distance (\u03BCm)", fontsize=22)
    plt.ylabel("Pair's Correlation Score", fontsize=22)
    plt.ylim(-0.2, 1)
    ax.axhline(y=0, color='k')

    ax.axvline(x=14, color='grey', linestyle='--', linewidth=2)
    plt.text(15, max(corrs) / 1.3, 'Close cell pair (\u2264 14 \u03BCm)', rotation=90, fontsize=15)

    ax.set_xticks(x_ticks_major, major=True)
    ax.tick_params(width=2, length=10, which='major')

    custom = [Line2D([], [], marker='.', color='tab:blue', linestyle='None', markersize=14),
              Line2D([], [], marker='.', color='tab:red', linestyle='None', markersize=14)]

    plt.legend(handles=custom, labels=['Low density pairs', 'High density pairs'], loc="upper right", fontsize=14)

    plt.tight_layout()
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI1_A.png')
    plt.close(fig=fig)


def pairs_correlation_by_distance(cells_data, total_cells, time_series_len):
    """ Plots the correlation score of pairs of cells by their distance. """
    arcos_df = pd.read_csv(f'{ARCOS_OUTPUT_LOCATION}/mid_third_wild_type/sample_8.csv')
    pairs_above_correlation_cutoff_and_pval_cutoff = xcorr_w_lags(cells_data, total_cells, time_series_len, 0, correlation_cutoff=-1, p_val_cutoff=1)
    distances = [pair_dict['dist'] for _, pair_dict in pairs_above_correlation_cutoff_and_pval_cutoff.items()]
    scores = [pair_dict['score'] for _, pair_dict in pairs_above_correlation_cutoff_and_pval_cutoff.items()]

    cell_ids = list(range(arcos_df['cell_id'].max()+1))
    relevant_pairs = {pair: pair_dict for pair, pair_dict in
                      pairs_above_correlation_cutoff_and_pval_cutoff.items()}

    all_points = [Point(arcos_df[arcos_df['cell_id'] == cid]['x_microns'].mean(),
                        arcos_df[arcos_df['cell_id'] == cid]['y_microns'].mean()) for cid in cell_ids]
    densities = [_calc_mean_density_in_pair_region(arcos_df, pair, r=7, all_points=all_points) for pair in relevant_pairs.keys()]
    plot(distances, densities, scores)


def main():
    relevant_experiment_folder = f'{PROCESSED_TIMESERIES_LOCATION}/mid_third_wild_type/sample_8/'
    total_cells = len(list(filter(lambda file_name: file_name.startswith('cell_') and file_name.endswith('.csv'),
                                  os.listdir(relevant_experiment_folder)
                                  )))
    cells_data = _get_cells_data_from_csv(relevant_experiment_folder, total_cells)
    pairs_correlation_by_distance(cells_data, total_cells, len(cells_data[0]))
