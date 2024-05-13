from typing import List
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from common.constants import SI_FIGURES_LOCATION

def plot(df: pd.DataFrame, x_col: str, y_col: str, title: str,
                              selected_experiments: List[str] = None):
    fig, ax = plt.subplots(figsize=(8, 8))
    prefix = title

    if selected_experiments:
        df = df[df['experiment_type'].isin(selected_experiments)]

    experiments_order = ['mid_third_wild_type','mid_third_zpg_RNAi',
                         'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout']

    m, b = np.polyfit(df[df['experiment_type'] == 'mid_third_wild_type'][x_col].to_list(),
                      df[df['experiment_type'] == 'mid_third_wild_type'][y_col].to_list(), 1)

    # add linear regression line to scatterplot
    plt.plot(df[df['experiment_type'] == 'mid_third_wild_type'][x_col].to_list(),
             m * np.array(df[df['experiment_type'] == 'mid_third_wild_type'][x_col].to_list()) + b, linestyle=':',
             label=f'Wild-Type Fit')

    for experiment_type in experiments_order:
        if experiment_type not in df['experiment_type'].unique():
            continue
        experiment_df = df[df['experiment_type'] == experiment_type]
        densities = experiment_df[x_col].to_list()
        activations = experiment_df[y_col].to_list()
        plt.scatter(densities, activations, label=experiment_type)

    xlabel = "Mean Local Density (#)" if x_col == "mean_density" else "Cell Activation Probability (%)"
    ylabel = "Cell Probability for Event (%)" if y_col == "MEC_per_minute" else "Cell Activation Probability (%)"
    plt.xlabel(xlabel, fontsize=22)
    plt.ylabel(ylabel, fontsize=22)

    x_ticks = [0, 1, 2, 3] if x_col == 'mean_density' else [0, 0.15, 0.3]
    y_ticks = [0, 0.06, 0.12] if y_col == 'MEC_per_minute' else [0, 0.15, 0.3]

    plt.xticks(x_ticks, fontsize=16)
    plt.yticks(y_ticks, fontsize=16, rotation=90)
    if len(selected_experiments) > 1:
        plt.legend(fontsize=14)

    plt.tight_layout()

    if prefix == 'WT':
        if y_col == 'MEC_per_minute' and x_col == 'mean_density':
            plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI4_A.png')
        elif y_col == 'MEC_per_minute' and x_col == 'mean_activations_per_minute':
            plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI4_B.png')
        elif y_col == 'mean_activations_per_minute' and x_col == 'mean_density':
            plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI4_C.png')
    elif prefix == 'ALL':
        if y_col == 'MEC_per_minute' and x_col == 'mean_density':
            plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI7_A.png')
        elif y_col == 'MEC_per_minute' and x_col == 'mean_activations_per_minute':
            plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI7_B.png')
        elif y_col == 'mean_activations_per_minute' and x_col == 'mean_density':
            plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI7_C.png')
    elif prefix == 'WT_ZPG' and y_col == 'MEC_per_minute' and x_col == 'mean_activations_per_minute':
        plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI7_D.png')

    plt.close(fig=fig)


def main():
    df = pd.read_csv('./data/cross_experiment_arcos_agg_data.csv')
    plot(df, 'mean_density', 'MEC_per_minute', title='WT', selected_experiments=['mid_third_wild_type'])
    plot(df, 'mean_density', 'MEC_per_minute', title='ALL', selected_experiments=['mid_third_wild_type', 'mid_third_zpg_RNAi', 'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout'])
    plot(df, 'mean_density', 'mean_activations_per_minute', title='WT', selected_experiments=['mid_third_wild_type'])
    plot(df, 'mean_density', 'mean_activations_per_minute', title='ALL', selected_experiments=['mid_third_wild_type', 'mid_third_zpg_RNAi', 'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout'])
    plot(df, 'mean_activations_per_minute', 'MEC_per_minute', title='WT', selected_experiments=['mid_third_wild_type'])
    plot(df, 'mean_activations_per_minute', 'MEC_per_minute', title='WT_ZPG', selected_experiments=['mid_third_wild_type', 'mid_third_zpg_RNAi'])
    plot(df, 'mean_activations_per_minute', 'MEC_per_minute', title='ALL', selected_experiments=['mid_third_wild_type', 'mid_third_zpg_RNAi', 'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout'])
