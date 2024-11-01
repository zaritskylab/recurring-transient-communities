import matplotlib.pyplot as plt
import pandas as pd

from common.constants import SI_FIGURES_LOCATION


def plot():
    fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(15, 10))
    si_df = pd.read_csv('./data/cross_experiment_arcos_agg_si_data.csv')

    for ix, row in si_df.iterrows():
        curr_ax = axs[ix // 4, ix % 4]
        experiment_name = row['experiment_name']

        sim_mec = [float(ele) for ele in row['simulated_MEC'][1:-1].split()]
        min_x = round(min(sim_mec) ,2)
        max_x = round(max(sim_mec) ,2)

        if max_x <= 1:
            text_ratio = 0.9
        elif max_x <= 2:
            text_ratio = 0.95
        else:
            text_ratio = 0.97

        bins = curr_ax.hist(sim_mec, bins=10)
        curr_ax.axvline(x=row['MEC'], ymin=0, linewidth=2, color='k', linestyle='--')
        curr_ax.text(x=row['MEC'] * text_ratio, y=max(bins[0]) / 2, s=f'Actual MEC',
                     rotation=90, fontsize=9)
        curr_ax.set_xticks([min_x, round(min_x + ((max_x -min_x ) /2), 2), max_x], labels=[min_x, round(min_x + ((max_x -min_x ) /2) ,2), max_x])
        curr_ax.set_ylim(0, 300)
        curr_ax.set_yticks([0 ,150 ,300], labels=[0 ,150 ,300])

        if experiment_name == 'sample_12':
            color = 'firebrick'
            for axis in ['top', 'bottom', 'left', 'right']:
                curr_ax.spines[axis].set_linewidth(2.5)  # change width
                curr_ax.spines[axis].set_color(color)

    fig.text(0.5, 0.04, 'Mean Events per Cell (#events/cell)', ha='center', fontsize=16)
    fig.text(0.04, 0.5, '#Simulated Experiments', va='center', rotation='vertical', fontsize=16)
    plt.savefig(f'{SI_FIGURES_LOCATION}Fig_SI5.png')
    plt.close()

def main():
    plot()