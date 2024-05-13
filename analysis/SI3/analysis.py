import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common.constants import SI_FIGURES_LOCATION, PROCESSED_TIMESERIES_LOCATION, ARCOS_OUTPUT_LOCATION
from data_layer.arcos.hotspots_mapper import EXPERIMENTS_HOT_SPOTS


def plot_spike_magnitude(df):
    fig = plt.figure()
    sns.violinplot(x=df.type, y=df.zscore, inner="quart", linewidth=1)
    plt.ylabel('Spike Magnitude (z-score)')
    plt.xlabel('')
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI3.png')
    plt.close(fig=fig)


def main():
    wt_folders = os.listdir(f'{PROCESSED_TIMESERIES_LOCATION}/mid_third_wild_type/')
    zpg_folders = os.listdir(f'{PROCESSED_TIMESERIES_LOCATION}/mid_third_zpg_RNAi/')

    wt_noise_model_by_experiment = {}
    wt_active_model_by_experiment = {}
    zpg_noise_model_by_experiment = {}
    zpg_active_model_by_experiment = {}

    non_community_zscores = []
    in_community_zscores = []

    for folder_type in ['mid_third_wild_type', 'mid_third_zpg_RNAi']:
        if folder_type == 'mid_third_wild_type':
            folders = wt_folders
            noise_model_by_experiment = wt_noise_model_by_experiment
            active_model_by_experiment = wt_active_model_by_experiment
        else:
            folders = zpg_folders
            noise_model_by_experiment = zpg_noise_model_by_experiment
            active_model_by_experiment = zpg_active_model_by_experiment

        for folder in folders:
            arcos_df = pd.read_csv(f'{ARCOS_OUTPUT_LOCATION}/{folder_type}/{folder}.csv')
            num_of_cells = arcos_df['cell_id'].nunique()

            experiment_per_cell_noise = {}
            experiment_per_cell_active = {}

            ########################
            ###### HOTSPOTS ########
            ########################

            validated_hotspots = ['sample_6__1', 'sample_6__2', 'sample_6__3', 'sample_10__0', 'sample_10__1',
                                  'sample_7__0', 'sample_9__0', 'sample_11__0']
            relevant_hotspots = [hotspot for hotspot in EXPERIMENTS_HOT_SPOTS
                                 if folder_type == hotspot.experiment_type and folder == hotspot.experiment_name
                                 and f'{hotspot.experiment_name}__{hotspot.id}' in validated_hotspots]
            ########################
            ###### HOTSPOTS ########
            ########################

            for cell_id in range(num_of_cells):
                cell_df = arcos_df[arcos_df['cell_id'] == cell_id]
                active_frames = np.where(cell_df['normalized_intensity.bin'] == 1)[0]

                # remove frames that are around active frames from noise calculation
                around_active_frames = np.where(cell_df['normalized_intensity.bin'].rolling(7, center=True).sum() > 0)[
                    0]

                inactive_frames = np.where(cell_df['normalized_intensity.bin'] == 0)[0]
                inactive_frames = np.setdiff1d(inactive_frames, around_active_frames)

                if len(active_frames) == 0:
                    continue

                cell_noise_mean = np.mean(cell_df['mean_intensity'].iloc[inactive_frames])
                cell_noise_std = np.std(cell_df['mean_intensity'].iloc[inactive_frames])
                cell_peaks_zscores = (cell_df['mean_intensity'].iloc[active_frames] - cell_noise_mean) / cell_noise_std
                cell_mean_zscore = np.mean(cell_peaks_zscores)

                experiment_per_cell_noise[cell_id] = {'noise_mean': cell_noise_mean, 'noise_std': cell_noise_std}
                experiment_per_cell_active[cell_id] = {'cell_mean_zscore': cell_mean_zscore}

                if cell_id in [hotspot_cell for hotspot in relevant_hotspots for hotspot_cell in hotspot.hotspot_cells]:
                    in_community_zscores.append(cell_mean_zscore)
                elif cell_df['collid'].nunique() > 1:
                    in_community_zscores.append(cell_mean_zscore)
                else:
                    non_community_zscores.append(cell_mean_zscore)

            noise_model_by_experiment[folder] = experiment_per_cell_noise
            active_model_by_experiment[folder] = experiment_per_cell_active


    df = pd.DataFrame(data=
                      [{'type': 'Community cells', 'zscore': zscore} for zscore in in_community_zscores] +
                      [{'type': 'Non-community cells', 'zscore': zscore} for zscore in non_community_zscores])
    plot_spike_magnitude(df)
