from copy import deepcopy
import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import math
from shapely.geometry import Point, Polygon
from common.constants import FIGURES_LOCATION, ARCOS_OUTPUT_LOCATION


def plot(all_event_heatmaps, all_stacked_heatmaps, all_active_cells_indices, cells_xs, cells_ys, max_limit=None):
    def _agg_list(original_list, agg_size):
        sublists = []
        for i in range(0, len(original_list), agg_size):
            sublist = original_list[i:i + agg_size]
            sublists.append(sublist)
        return sublists

    n_groups = 9
    dot_size = 1

    heatmaps_per_agg_group = int(math.ceil(len(all_event_heatmaps) / n_groups))
    tmp_all_event_heatmaps = _agg_list(all_event_heatmaps, heatmaps_per_agg_group)
    tmp_all_stacked_heatmaps = _agg_list(all_stacked_heatmaps, heatmaps_per_agg_group)
    tmp_all_active_cells_indices = _agg_list(all_active_cells_indices, heatmaps_per_agg_group)

    all_event_heatmaps = [np.max(array_list, axis=0) for array_list in tmp_all_event_heatmaps]
    all_stacked_heatmaps = [np.max(array_list, axis=0) for array_list in tmp_all_stacked_heatmaps]
    all_active_cells_indices = [sum(sublist, []) for sublist in tmp_all_active_cells_indices]

    fig, axes = plt.subplots(2, n_groups, figsize=(32, 8))

    for ix in range(len(all_stacked_heatmaps)):
        if ix >= 9:
            continue
        event_heatmap = all_event_heatmaps[ix]
        event_cells = all_active_cells_indices[ix]
        ax_event = axes[0, ix]
        stacked_heatmap = all_stacked_heatmaps[ix]
        ax_stacked = axes[1, ix]

        if max_limit:
            stacked_max_value = event_max_value = max_limit
        else:
            stacked_max_value = int(np.max(all_stacked_heatmaps) + 1)
            event_max_value = int(np.max(all_event_heatmaps)) + 0.1

        cmap_stacked = plt.cm.get_cmap('hot')
        cmap_event = plt.cm.get_cmap('bone')

        # normalize chosen colormap
        stacked_norm = matplotlib.colors.Normalize(vmin=0, vmax=stacked_max_value, clip=True)
        event_norm = matplotlib.colors.Normalize(vmin=0, vmax=event_max_value, clip=True)

        ax_stacked.imshow(stacked_heatmap.T, interpolation='nearest', cmap=cmap_stacked, norm=stacked_norm)
        ax_stacked.scatter(cells_xs, cells_ys, s=dot_size, facecolors='white', edgecolors='white', linewidth=2)

        ax_event.imshow(event_heatmap.T, interpolation='nearest', cmap=cmap_event, norm=event_norm)
        ax_event.scatter([x for x_id, x in enumerate(cells_xs) if x_id not in event_cells],
                         [y for y_id, y in enumerate(cells_ys) if y_id not in event_cells], s=dot_size,
                         facecolors='white', edgecolors='white', linewidth=2)

        ax_event.scatter([x for x_id, x in enumerate(cells_xs) if x_id in event_cells],
                         [y for y_id, y in enumerate(cells_ys) if y_id in event_cells], s=dot_size, facecolors='red',
                         edgecolors='red', linewidth=2)

        min_x = min(cells_xs)
        min_y = min(cells_ys)
        max_x = max(cells_xs)
        max_y = max(cells_ys)

        ax_stacked.set_ylim(min_y - 10, max_y + 10)
        ax_stacked.set_xlim(min_x - 10, max_x + 10)
        ax_stacked.set_xticks([])
        ax_stacked.set_yticks([])

        ax_event.set_ylim(min_y - 10, max_y + 10)
        ax_event.set_xlim(min_x - 10, max_x + 10)
        ax_event.set_xticks([])
        ax_event.set_yticks([])

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}Fig_2_A.png')
    plt.close(fig=fig)

def _analyze_hotspot_formation_in_time(arcos_df):

    cells_xs, cells_ys, polygons = [], [], []
    distinct_cells = set()

    for i in range(arcos_df['cell_id'].max()+1):
        current_cell = arcos_df[arcos_df['cell_id']==i]
        cells_xs.append(current_cell['x'].mean())
        cells_ys.append(current_cell['y'].mean())

    ## Create heatmap, each cell is 1micron X 1micron
    all_stacked_heat_map = []
    all_event_heat_map = []
    all_event_cell_ids = []
    stacked_heat_map = np.zeros(shape=(math.ceil(max(arcos_df['x']) + 10), math.ceil(max(arcos_df['y']) + 10)))

    all_events_df = arcos_df[~arcos_df['collid'].isna()]
    event_first_frame_by_event_id = all_events_df.groupby('collid')['frame'].min().to_dict()
    all_events_id = all_events_df['collid'].drop_duplicates().to_list()
    all_events_id = sorted(all_events_id, key=lambda x: event_first_frame_by_event_id[x])

    for ix, event_id in enumerate(all_events_id):
        event_heat_map = np.zeros(shape=(math.ceil(max(arcos_df['x']) + 10), math.ceil(max(arcos_df['y']) + 10)))

        event_df = all_events_df[all_events_df['collid'] == event_id]

        distinct_participating_cells = event_df['cell_id'].drop_duplicates().to_list()

        cells_coordinates = []
        for cell_id in distinct_participating_cells:
            distinct_cells.add(cell_id)
            cell_df = event_df[event_df['cell_id'] == cell_id]
            cell_x = cell_df['x'].mean()
            cell_y = cell_df['y'].mean()
            cells_coordinates.append(
                Point(math.floor(cell_x), math.floor(cell_y))
            )

        try:
            if len(cells_coordinates) >= 3:
                convex_poly = Polygon([[p.x, p.y] for p in cells_coordinates]).convex_hull
                polygons.append(convex_poly)

                min_x, min_y, max_x, max_y = convex_poly.bounds
                for i in range(int(min_x), int(max_x+1)):
                    for j in range(int(min_y), int(max_y+1)):
                        tmp_point = Point(i,j)
                        if tmp_point.distance(convex_poly) < 0.75:  # Visual threshold
                            event_heat_map[i][j] = 1
                        if tmp_point.within(convex_poly):
                            stacked_heat_map[i][j] += 1
        except Exception as e:
            print("exception")

        all_stacked_heat_map.append(deepcopy(stacked_heat_map))
        all_event_heat_map.append(event_heat_map)
        all_event_cell_ids.append(distinct_participating_cells)

    plot(all_event_heat_map, all_stacked_heat_map, all_event_cell_ids, cells_xs, cells_ys)


def main():
    relevant_experiment_file = f'{ARCOS_OUTPUT_LOCATION}/mid_third_wild_type/sample_11.csv'
    arcos_df = pd.read_csv(relevant_experiment_file)
    _analyze_hotspot_formation_in_time(arcos_df)


