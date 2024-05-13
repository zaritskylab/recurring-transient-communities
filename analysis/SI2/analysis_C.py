import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from common.constants import SI_FIGURES_LOCATION


def plot(cross_experiment_stats: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(8, 8))
    wt_data = cross_experiment_stats[cross_experiment_stats['experiment_type'] == 'mid_third_wild_type']['community_size']
    bins = np.arange(int(wt_data.max()) + 1) - 0.5
    plt.hist(wt_data, bins, density=True)
    plt.ylabel('Probability', fontsize=24)
    plt.xlabel('Community size', fontsize=24)
    plt.xlim(2.5,15.5)
    plt.xticks([3,7,11,15], fontsize=18)
    plt.yticks(rotation=90, fontsize=18)
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI2_C.png')
    plt.close()


def main():
    cross_experiment_stats = pd.read_csv(f'./data/cross_experiment_communities_statistics.csv')
    plot(cross_experiment_stats)

