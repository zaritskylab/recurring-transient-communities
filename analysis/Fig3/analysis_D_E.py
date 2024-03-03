import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from common.constants import FIGURES_LOCATION
from common.significance_calculator import statistical_test_wrapper, StatisticalTest

def plot(df: pd.DataFrame, experiments_order: list, fig_name: str):

    fig, ax = plt.subplots(figsize=(8, 8))
    x_labels = experiments_order

    field = 'MEC_magnitude'
    y_label = 'Mean MEC Magnitude'
    ymax = df[field].max() if df[field].max() >= 5 else 5
    plt.ylim(-0.3, ymax)

    pvalues = statistical_test_wrapper(df[['experiment_type', field]], field,
                                       test_type=StatisticalTest.KRUSKAL_WALLIS,
                                       static_experiment_type='wild_type_data')

    axis_map = {}
    pivot_df = df[['experiment_type', field]].pivot(columns='experiment_type', values=field)

    for i, d in enumerate(experiments_order):
        y = pivot_df[d].dropna()
        x = np.random.normal(i + 1, 0.1, len(y))
        axis, = plt.plot(x, y, mec='k', ms=8, marker="o", linestyle="None")
        axis_map[d] = axis

    plt.suptitle('')
    plt.axhline(y=1, color='k', linestyle='--', linewidth=2)
    plt.xticks([i + 1 for i in range(len(x_labels))], rotation=30, fontsize=16, labels=x_labels)
    plt.yticks([0, 3, 6], fontsize=16, rotation=90)

    plt.xlabel('Experiment Type', fontsize=22)
    plt.ylabel("Magnitude of Change in #ARCOS-Events", fontsize=22)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}/Fig_3_{fig_name}.png')
    plt.close(fig=fig)

def main():
    selected_experiments = ['late_second_wild_type', 'early_third_wild_type', 'mid_third_wild_type']
    df = pd.read_csv(f'./data/cross_experiment_arcos_agg_data.csv')
    df = df[df['experiment_type'].isin(selected_experiments)]
    plot(df, experiments_order=selected_experiments, fig_name='D')

    selected_experiments = ['late_second_zpg_RNAi', 'early_third_zpg_RNAi', 'mid_third_zpg_RNAi']
    df = pd.read_csv(f'./data/cross_experiment_arcos_agg_data.csv')
    df = df[df['experiment_type'].isin(selected_experiments)]
    plot(df, experiments_order=selected_experiments, fig_name='E')