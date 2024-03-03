import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from common.constants import FIGURES_LOCATION
from matplotlib.lines import Line2D


def plot(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(5, 8))
    max_y = 0
    y = df['MEC_magnitude']
    pvalues = df['p_value']
    x = np.random.normal(1, 0.2, len(y))
    max_y = y.max() if y.max() > max_y else max_y
    plt.scatter(x, y, s=80, marker="o", linestyle="None", label='mid_third_wild_type',
                facecolors=['#1f77b4' if p < 0.05 else 'none' for p in pvalues], edgecolors='black',
                color=['#1f77b4' if p < 0.05 else '#ff7f0e' for p in pvalues])
    plt.xticks([1], fontsize=16, labels=['mid_third_wild_type'])
    plt.yticks([0, 2, 4], fontsize=16, rotation=90)
    plt.xlabel('Experiment Type', fontsize=22)
    plt.ylabel('MEC Magnitude', fontsize=22)

    line1 = Line2D([], [], color="white", marker='o', markerfacecolor='#1f77b4', markeredgecolor='black',
                   markersize=10)
    line2 = Line2D([], [], color="white", marker='o', markerfacecolor="none", markeredgecolor='black',
                   markersize=10)
    plt.legend((line1, line2), ('Significant Experiment', 'Insignificant Experiment'), numpoints=1, loc=1, fontsize=14)
    plt.axhline(y=1, color='black', linestyle='--', linewidth=2)
    plt.ylim(ymin=0, ymax=4)
    plt.xlim(xmin=0, xmax=2)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}Fig_1_F.png')
    plt.close(fig=fig)

def main():
    df = pd.read_csv('./data/cross_experiment_arcos_agg_data.csv')
    plot(df[df['experiment_type'] == 'mid_third_wild_type'])
