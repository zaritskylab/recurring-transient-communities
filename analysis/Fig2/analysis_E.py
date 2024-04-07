import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from common.constants import FIGURES_LOCATION, HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION
from common.significance_calculator import statistical_test_wrapper, StatisticalTest


def plot(df: pd.DataFrame, x_labels_and_order):

    x = [i + 1 for i in range(len(x_labels_and_order))]

    fig, ax = plt.subplots(figsize=(20, 10))
    cmap = plt.cm.Greys
    axis_map = {}

    y_hotspots = df['hotspots']
    y_enough_statistics = df['enough_statistics']
    y_significant = df['is_significant']

    axis_gray = plt.bar(x, y_enough_statistics, width=0.8, color='lightgrey')

    axis = plt.bar(x, y_significant, width=0.6)
    x_start = np.array([plt.getp(item, 'x') - 0.1 for item in axis])
    x_end = x_start + [plt.getp(item, 'width') + 0.2 for item in axis]

    plt.hlines(y_hotspots, x_start, x_end, linestyle='--', color='black', linewidth=4)

    for i, start in enumerate(x_start):
        plt.text(x=start, y=y_significant[i] + 0.1, s=f'{y_significant[i]}', fontdict={'fontsize': 16}) if \
            y_significant[i] > 0 else None

        plt.text(x=start, y=y_enough_statistics[i] + 0.1, s=f'{y_enough_statistics[i]}',
                 fontdict={'fontsize': 16}) if y_enough_statistics[i] > 0 else None
        plt.text(x=start, y=y_hotspots[i] + 0.1, s=f'{y_hotspots[i]}', fontdict={'fontsize': 16})

    x_ticks = [i + 1 for i in range(len(x_labels_and_order))]
    plt.xticks(x_ticks, rotation=30, fontsize=18, labels=x_labels_and_order)
    plt.yticks([i for i in range(0, 20)], fontsize=18)

    plt.xlabel('Experiment Type', fontsize=24)
    main_y_axis = "#Validated Hotspots"
    plt.ylabel(main_y_axis, fontsize=24)

    # create manual symbols for legend
    patch1 = mpatches.Patch(facecolor=axis[0].get_facecolor(), label=main_y_axis)
    handles, _ = plt.gca().get_legend_handles_labels()

    patch2 = mpatches.Patch(facecolor='lightgrey', label='#of Hotspots with >100 permutations')
    line1 = Line2D([0], [0], label='Total Hotspots', color='black', linestyle='--', linewidth=4)
    handles.extend([patch1, patch2, line1])

    plt.legend(handles=handles, fontsize=18, loc='upper right')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}Fig_2_E.png')
    plt.close(fig=fig)



def analyze_aggregated_data():
    df = pd.read_csv(f'{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/aggregated_results.csv')
    df['hotspots'] = 1

    x_labels_and_order = ['mid_third_wild_type', 'mid_third_zpg_RNAi', 'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout']
    df = df[df['experiment_type'].isin(x_labels_and_order)]

    # Calculating the p-values based on the is_significant field
    p_val_field = 'is_significant'
    df[p_val_field] = df.apply(lambda row: 1 if row['pvalue'] < 0.05 and row['permutations_count'] > 100 else 0, axis=1)
    df['enough_statistics'] = df.apply(lambda row: 1 if row['permutations_count'] > 100 else 0, axis=1)

    pvalues = statistical_test_wrapper(df[['experiment_type', p_val_field]], p_val_field,
                                       test_type=StatisticalTest.FISHERS_EXACT_TEST,
                                       static_experiment_type='mid_third_wild_type')
    df = df[['experiment_type', 'is_significant', 'enough_statistics', 'hotspots']].groupby('experiment_type',
                                                                as_index=False).sum()

    # Adding rows for experiments that don't exist in the data
    for label in x_labels_and_order:
        if label not in df['experiment_type'].unique():
            df = df.append({'experiment_type': label,
                            'is_significant': 0,
                            'enough_statistics': 0,
                            'hotspots': 0}, ignore_index=True)

    df['order'] = df['experiment_type'].apply(lambda x: x_labels_and_order.index(x))
    df.sort_values(by='order', inplace=True)
    df.drop(columns=['order'], inplace=True)
    df.reset_index(inplace=True, drop=True)

    plot(df, x_labels_and_order)


if __name__ == '__main__':
    analyze_aggregated_data()