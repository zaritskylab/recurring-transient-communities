from typing import List
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from common.constants import FIGURES_LOCATION
from common.frame_to_frequency import FrequencyTranslator
from common.significance_calculator import statistical_test_wrapper, StatisticalTest


def plot(df: pd.DataFrame, experiments_order: List[str]):
    fig, ax = plt.subplots(figsize=(8, 8))
    field = 'info_spread_rate (avg between 2 adjacent cells) - cells per second'

    pvalues = statistical_test_wrapper(
        df[df['experiment_type'].isin(experiments_order)][['experiment_type', field]], field,
        test_type=StatisticalTest.KRUSKAL_WALLIS,
        static_experiment_type='wild_type_data')

    axis_map = {}
    for i, d in enumerate(experiments_order):
        pivot_df = df[df['experiment_type'] == d]
        y = pivot_df[field]
        y_mean = y.mean()

        x = np.random.normal(i + 1, 0.1, len(y))
        axis, = plt.plot(x, y, mec='k', ms=6, marker="o", linestyle="None")
        plt.plot(i + 1, y_mean, mec='k', ms=12, marker="o", linestyle="None", color='red')
        plt.text(i + 0.7, y_mean * 1.05, f"{y_mean:.2f}", horizontalalignment='center', fontsize=18,
                 fontdict={'weight': 'bold'})
        axis_map[d] = axis

    plt.xlim(0, len(experiments_order) + 0.7)
    plt.ylabel('Info Spread Rate (\u03BCm per second)', fontsize=18)
    x_ticks = [i + 1 for i in range(len(experiments_order))]
    plt.xticks(x_ticks, rotation=30, fontsize=18, labels=experiments_order)
    plt.yticks(fontsize=18, ticks=[0, 3, 6], rotation=90)
    plt.tight_layout()

    plt.savefig(f'{FIGURES_LOCATION}/Fig_1_I.png')
    plt.close(fig)



def main():
    communities_stats_df = pd.read_csv('./data/cross_experiment_communities_statistics.csv')

    # Limit to experiments with FPS of > 0.25 (frames per second)
    FPS = 0.25
    freqs_df = FrequencyTranslator().get_records_below_fps(FPS)[['experiment_type', 'experiment_name']]

    narrowed_df = pd.DataFrame()
    for _, row in freqs_df.iterrows():
        _df = communities_stats_df[(communities_stats_df['experiment_type'] == row['experiment_type']) &
                                      (communities_stats_df['experiment_name'] == row['experiment_name'])]
        narrowed_df = pd.concat((narrowed_df, _df))

    experiments_order = ['mid_third_wild_type','mid_third_zpg_RNAi','cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout']
    plot(narrowed_df, experiments_order)
