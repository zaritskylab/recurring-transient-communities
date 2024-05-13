import os
from typing import List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common.config import LAG_IN_SECONDS
from common.constants import SI_FIGURES_LOCATION, PROCESSED_TIMESERIES_LOCATION, ARCOS_COMMUNITIES_STATS_LOCATION
from common.frame_to_frequency import FrequencyTranslator
from data_layer.arcos.hotspots_mapper import EXPERIMENTS_HOT_SPOTS
from data_layer.cross_correlation import xcorr_w_lags
from data_layer.preprocess_data import _get_cells_data_from_csv


def plot(df, title, group_a, group_b, bins, figname):
    relevant_df = df[df['bin'].isin(bins)]

    # Here we plot the scores of pairs that are in the same community vs. the scores of pairs that are not in the same community,
    # as a scatter plot (using sns package) with distance bins (+random float) as x axis
    fig = plt.figure()
    sns.violinplot(x='bin', y='score', hue='type', data=relevant_df, split=False, inner="quart", linewidth=2)
    plt.xlabel('Pairs Distance (\u03BCm)')
    plt.ylabel('Pairs Correlation Score')
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI6.png')
    plt.close(fig=fig)

def main():
    experiment_type = 'mid_third_wild_type'

    distance_bins = [7, 14, 21, 28, 35, 42]
    in_community_scores_by_bins = {bin: [] for bin in distance_bins}
    out_community_scores_by_bins = {bin: [] for bin in distance_bins}

    ######################################
    ########## COMMUNITY DATA ############
    ######################################

    data_folder = f"{ARCOS_COMMUNITIES_STATS_LOCATION}/{experiment_type}/"
    folders = os.listdir(f'{PROCESSED_TIMESERIES_LOCATION}/{experiment_type}')
    hotspot_scores_by_bins = {bin: [] for bin in distance_bins}
    community_scores_by_bins = {bin: [] for bin in distance_bins}

    for folder in folders:
        community_filename = folder + '_communities_statistics.csv'
        community_df = pd.read_csv(f"{data_folder}{community_filename}")

        num_of_cells = len(os.listdir(f'{PROCESSED_TIMESERIES_LOCATION}/{experiment_type}/{folder}')) - 1
        cells_data = _get_cells_data_from_csv(f'{PROCESSED_TIMESERIES_LOCATION}/{experiment_type}/{folder}', num_of_cells)
        time_series_len = len(cells_data[0])
        experiment_lag = FrequencyTranslator().temporal_lag_to_frames(experiment_type, folder, LAG_IN_SECONDS)

        lagged_pairs_passed_cutoff_max_corr_and_distance_and_lags = xcorr_w_lags(cells_data,
                                                                                 num_of_cells,
                                                                                 time_series_len,
                                                                                 experiment_lag,
                                                                                 correlation_cutoff=-1)

        distances = [pair_dict['dist'] for _, pair_dict in
                     lagged_pairs_passed_cutoff_max_corr_and_distance_and_lags.items()]
        scores = [pair_dict['score'] for _, pair_dict in
                  lagged_pairs_passed_cutoff_max_corr_and_distance_and_lags.items()]

        def _find_closest_ceil_bin(distance: int, bins: List[int]):
            for bin in bins:
                if distance <= bin:
                    return bin
            return None

        validated_hotspots = ['sample_6__1', 'sample_6__2', 'sample_6__3', 'sample_10__0', 'sample_10__1',
                              'sample_7__0', 'sample_9__0', 'sample_11__0']
        relevant_hotspots = [hotspot for hotspot in EXPERIMENTS_HOT_SPOTS
                             if experiment_type == hotspot.experiment_type and folder == hotspot.experiment_name
                             and f'{hotspot.experiment_name}__{hotspot.id}' in validated_hotspots]

        for ix, distance in enumerate(distances):
            relevant_bin = _find_closest_ceil_bin(distance, distance_bins)
            if not relevant_bin:
                continue
            pair = list(lagged_pairs_passed_cutoff_max_corr_and_distance_and_lags.keys())[ix]

            is_community_pair = False
            is_A_appear_in_community = False
            is_B_appear_in_community = False
            #################################
            ########## Communities ##########
            #################################

            for _, community_row in community_df.iterrows():
                community_cells = eval(community_row['cells'])
                if pair[0] in community_cells and pair[1] in community_cells:
                    is_community_pair = True
                    break
                elif pair[0] in community_cells:
                    is_A_appear_in_community = True
                elif pair[1] in community_cells:
                    is_B_appear_in_community = True

            if is_community_pair:
                in_community_scores_by_bins[relevant_bin].append(scores[ix])
            elif is_A_appear_in_community and is_B_appear_in_community:
                out_community_scores_by_bins[relevant_bin].append(scores[ix])

            #################################
            ########## Hotspots #############
            #################################

        for ix, distance in enumerate(distances):
            relevant_bin = _find_closest_ceil_bin(distance, distance_bins)
            if not relevant_bin:
                continue
            pair = list(lagged_pairs_passed_cutoff_max_corr_and_distance_and_lags.keys())[ix]

            is_community_pair = False
            is_hotspot_cell = False

            for _, community_row in community_df.iterrows():
                community_cells = eval(community_row['cells'])
                if pair[0] in community_cells and pair[1] in community_cells:
                    is_community_pair = True
                    break

            if is_community_pair:
                for hotspot in relevant_hotspots:
                    hotspot_cells = hotspot.hotspot_cells
                    if pair[0] in hotspot_cells and pair[1] in hotspot_cells:
                        is_hotspot_cell = True
                        break

            if is_hotspot_cell:
                hotspot_scores_by_bins[relevant_bin].append(scores[ix])
            if is_community_pair:
                community_scores_by_bins[relevant_bin].append(scores[ix])

    in_df = pd.DataFrame(data=[{'type': 'intra community', 'bin': curr_bin, 'score': score}
                               for curr_bin, scores in in_community_scores_by_bins.items() for score in scores])
    out_df = pd.DataFrame(data=[{'type': 'inter community', 'bin': curr_bin, 'score': score}
                                for curr_bin, scores in out_community_scores_by_bins.items() for score in scores])
    unify_df = pd.concat([in_df, out_df])

    plot(unify_df, 'community_vs_not_community', group_a='intra community', group_b='inter community', bins=distance_bins,
         figname='reviewer_2_comment_3_1')
