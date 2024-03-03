import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from common.constants import FIGURES_LOCATION
from common.significance_calculator import statistical_test_wrapper, StatisticalTest


def plot(df: pd.DataFrame, experiments_order: list):
    field = 'p_value'
    x_labels = experiments_order

    pvalues = statistical_test_wrapper(df[['experiment_type', field]], field,
                                       test_type=StatisticalTest.FISHERS_EXACT_TEST,
                                       static_experiment_type='early_third_wild_type')

    fig, ax = plt.subplots(figsize=(8, 10))
    y = df[field]

    pivot_df = df[['experiment_type', field]].pivot(columns='experiment_type', values=field)
    axis_map = {}
    for i, d in enumerate(experiments_order):
        y_sign = sum(pivot_df[d].dropna() < 0.05)
        y_total = len(pivot_df[d].dropna())
        x = np.repeat(i + 1, len(y))
        axis = plt.bar(x, y_sign, width=0.8)
        _ = plt.bar(x, y_total - y_sign, width=0.8, bottom=y_sign, color='lightgrey')
        axis_map[d] = axis
    plt.suptitle('')
    x_ticks = [i + 1 for i in range(len(x_labels))]
    plt.xticks(x_ticks, rotation=30, fontsize=16, labels=x_labels)
    plt.yticks([0, 4, 8, 12], fontsize=16, rotation=90)
    plt.ylim(0, 14)

    plt.xlabel('Experiment Type', fontsize=24)
    plt.ylabel('#Spatially Significant Experiments', fontsize=24)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}/Fig_1_G.png')
    plt.close(fig=fig)

def main():
    selected_experiments = ['mid_third_wild_type', 'mid_third_zpg_RNAi', 'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout']
    df = pd.read_csv(f'./data/cross_experiment_arcos_agg_data.csv')
    df = df[df['experiment_type'].isin(selected_experiments)]
    plot(df, experiments_order=selected_experiments)