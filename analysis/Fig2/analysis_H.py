import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from common.constants import FIGURES_LOCATION, HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION, ARCOS_OUTPUT_LOCATION
from common.frame_to_frequency import FrequencyTranslator
from data_layer.arcos.hotspots_mapper import EXPERIMENTS_HOT_SPOTS


def plot(heat_map):
    fig, ax = plt.subplots(figsize=(8, 8))

    max_value = int(heat_map.max() + 1)
    cmap = plt.cm.get_cmap('hot')

    # normalize chosen colormap
    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_value, clip=True)

    plt.imshow(heat_map.T, interpolation='nearest', cmap=cmap, norm=norm, origin='lower')
    plt.xlabel("Hotspot cells in a community (#)", fontsize=18)
    plt.ylabel("Non-hotspot cells in a community (#)", fontsize=18)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, color='white', alpha=0.75, zorder=0)

    above = 0
    below = 0
    for i in range(len(heat_map)):
        for j in range(len(heat_map)):
            if i > j:
                below += heat_map[i][j]
            if j > i:
                above += heat_map[i][j]

    # Create a ScalarMappable object
    from matplotlib.cm import ScalarMappable
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Setting an empty array to make it work with colorbar

    # Add color bar
    plt.colorbar(sm, label='# of Communities')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}Fig_2_H.png')
    plt.close(fig=fig)


def analyze():

    agg_df = pd.read_csv(f'{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/aggregated_results.csv')

    # Calculating the p-values based on the is_significant field
    p_val_field = 'is_significant'
    agg_df[p_val_field] = agg_df.apply(lambda row: 1 if row['pvalue'] < 0.05 and row['permutations_count'] > 100 else 0, axis=1)
    agg_df['enough_statistics'] = agg_df.apply(lambda row: 1 if row['permutations_count'] > 100 else 0, axis=1)
    agg_df = pd.merge(agg_df, FrequencyTranslator().get_names_df(), on=['experiment_type', 'experiment_name'], how='inner')

    global_hs_prop = []
    global_non_hs_prop = []

    for experiment_hotspot in EXPERIMENTS_HOT_SPOTS:
        if experiment_hotspot.experiment_type != 'mid_third_wild_type':
            continue

        arcos_df = pd.read_csv(f'{ARCOS_OUTPUT_LOCATION}/{experiment_hotspot.experiment_type}/{experiment_hotspot.experiment_name}.csv')

        hotspot_cells = set(experiment_hotspot.hotspot_cells)
        arcos_df = arcos_df[~arcos_df['collid'].isna()]

        all_coll_events = arcos_df['collid'].drop_duplicates().values.tolist()
        relevant_coll_events = []

        for coll_event in all_coll_events:
            coll_event_cells = set(arcos_df[arcos_df['collid'] == coll_event]['cell_id'].drop_duplicates().values.tolist())
            if len(coll_event_cells.intersection(hotspot_cells)) > 0:
                relevant_coll_events.append(coll_event)

        # hs_proportion = []
        # hs_total_cells = len(hotspot_cells)
        # non_hs_total_cells = []
        # first_cells_in_hotspot = []
        agg_df_record = agg_df[(agg_df['experiment_type'] == experiment_hotspot.experiment_type) & (agg_df['experiment_name'] == experiment_hotspot.experiment_name) & (agg_df['hotspot_id'] == experiment_hotspot.id)]
        if len(agg_df_record) == 1 and int(agg_df_record['is_significant']) == 1:
            for coll_event in relevant_coll_events:
                coll_event_df = arcos_df[arcos_df['collid'] == coll_event]
                coll_event_cells = set(coll_event_df['cell_id'].drop_duplicates().values.tolist())

                # first_frame = coll_event_df['frame'].min()
                # first_cells = set(coll_event_df[coll_event_df['frame'] == first_frame]['cell_id'].drop_duplicates().values.tolist())
                # first_cells_in_hotspot.extend(list(first_cells))
                #
                # if len(first_cells.intersection(hotspot_cells)) == 0:
                #     inbound += 1
                # elif 0 < len(first_cells.intersection(hotspot_cells)) < len(first_cells):
                #     mixed_start += 1
                #     # start_in_hs_and_propagated += 1
                # elif 0 < len(first_cells.intersection(hotspot_cells)) and len(coll_event_cells.intersection(hotspot_cells)) == len(coll_event_cells):
                #     start_in_hs_and_stayed += 1
                # else:
                #     start_in_hs_and_propagated += 1

                # import itertools
                # first_frame_per_cell = coll_event_df[['cell_id', 'frame']].groupby('cell_id', as_index=True).min()
                coll_event_hs_cells = coll_event_cells.intersection(hotspot_cells)
                # hs_proportion.append(len(coll_event_hs_cells)/ len(coll_event_cells))
                # non_hs_total_cells.append(coll_event_cells.difference(coll_event_hs_cells))
                global_hs_prop.append(len(coll_event_hs_cells))
                global_non_hs_prop.append(len(coll_event_cells.difference(coll_event_hs_cells)))

                # hs_cells_with_non_hs_cells = list(itertools.product(coll_event_hs_cells, coll_event_cells.difference(coll_event_hs_cells)))
                # hs_cells_with_non_hs_cells = list(filter(lambda cells: _get_cells_dist(cells[0], cells[1], coll_event_df) <= 14, hs_cells_with_non_hs_cells))
                # hs_frames_with_non_hs_frames = [(first_frame_per_cell.loc[hsc]['frame'], first_frame_per_cell.loc[non_hsc]['frame']) for
                #         (hsc, non_hsc) in hs_cells_with_non_hs_cells]
                #
                # hs_before += sum([1 if hsc < non_hsc else 0 for (hsc, non_hsc) in hs_frames_with_non_hs_frames])
                # non_hs_before += sum([1 if non_hsc < hsc else 0 for (hsc, non_hsc) in hs_frames_with_non_hs_frames])

            # reduced_non_hs_total_cells = set()
            # for s in non_hs_total_cells:
            #     reduced_non_hs_total_cells = reduced_non_hs_total_cells.union(s)


    global_max = max(max(global_non_hs_prop), max(global_hs_prop)) + 1
    heatmap = np.zeros(shape=(global_max, global_max))
    for i, j in zip(global_hs_prop, global_non_hs_prop):
        heatmap[i][j] += 1

    plot(heatmap)

def main():
    analyze()

