import json
from copy import deepcopy
import matplotlib
import pandas as pd
from matplotlib import pyplot as plt
from common.config import LAG_IN_SECONDS
from common.constants import FIGURES_LOCATION, ARCOS_OUTPUT_LOCATION, HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION, \
    ARCOS_INPUT_LOCATION
from common.frame_to_frequency import FrequencyTranslator
from data_layer.arcos.arcos_hotspots_detection import find_arcos_hot_spots
from data_layer.arcos.arcos_params import DEFAULT_PARAMS
from data_layer.arcos.arcos_wrapper import arcos_wrapper
from data_layer.arcos_hotspots_shuffle_analysis import _swap_cells


def plot(heat_map, cells_xs, cells_ys, is_swapped: bool = False, max_limit=None,
         is_diff: bool = False, min_x=None, min_y=None, max_x=None, max_y=None):

    fig, ax = plt.subplots(figsize=(8, 8))

    if max_limit:
        max_value = max_limit
        min_value = -max_limit if is_diff else 0
    else:
        max_value = int(heat_map.max() + 1)
        min_value = -max_limit if is_diff else 0

    # normalize chosen colormap
    cmap = plt.cm.get_cmap('hot') if not is_diff else plt.cm.get_cmap('seismic')
    norm = matplotlib.colors.Normalize(vmin=min_value, vmax=max_value, clip=True)

    plt.imshow(heat_map.T, interpolation='nearest', cmap=cmap, norm=norm)
    plt.scatter(cells_xs, cells_ys, s=40, facecolors='none', edgecolors='green', linewidth=2)

    min_x = min(cells_xs) if not min_x else min_x
    min_y = min(cells_ys) if not min_y else min_y
    max_x = max(cells_xs) if not max_x else max_x
    max_y = max(cells_ys) if not max_y else max_y
    plt.ylim(min_y - 10, max_y + 10)
    plt.xlim(min_x - 10, max_x + 10)

    plt.xlabel("Cell's X-Coordinates (Micrometer)", fontsize=18)
    plt.ylabel("Cell's Y-Coordinates (Micrometer)", fontsize=18)

    # Create a ScalarMappable object
    from matplotlib.cm import ScalarMappable
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Setting an empty array to make it work with colorbar

    # Add color bar
    plt.colorbar(sm, label='#Communities', ticks=[min_value, (min_value + max_value) / 2, max_value]) if not is_diff \
        else plt.colorbar(sm, label='Difference in #Communities', ticks=[min_value, (min_value + max_value) / 2, max_value])

    plt.tight_layout()
    if is_swapped:
        plt.savefig(f'{FIGURES_LOCATION}Fig_2_D.png')
    else:
        plt.savefig(f'{FIGURES_LOCATION}Fig_2_C.png')
    plt.close(fig=fig)

def main():
    relevant_experiment_file = f'{ARCOS_OUTPUT_LOCATION}/mid_third_wild_type/sample_11.csv'
    arcos_df = pd.read_csv(relevant_experiment_file)

    # Figure 2C
    _, heat_map_orig, _, _, _, cells_xs, cells_ys, _ = find_arcos_hot_spots(arcos_df)
    max_limit = heat_map_orig.max()+1
    plot(heat_map_orig, cells_xs, cells_ys, is_swapped=False, max_limit=max_limit)


    # Figure 2D
    relevant_experiment_best_swap_file = f"{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/mid_third_wild_type/sample_11/permutations/0_best_perm.json"
    arcos_input_df = pd.read_csv(f'{ARCOS_INPUT_LOCATION}/mid_third_wild_type/sample_11.csv')
    with open(relevant_experiment_best_swap_file, 'r') as f:
        best_swap = json.load(f)

    old_cells = [swap['old_cell'] for swap in best_swap['swaps']]
    new_cells = [swap['new_cell'] for swap in best_swap['swaps']]

    for ix in range(len(old_cells)):
        arcos_input_df = _swap_cells(arcos_input_df, old_cells[ix], new_cells[ix])

    arcos_params = deepcopy(DEFAULT_PARAMS)
    arcos_params.nPrev = FrequencyTranslator().temporal_lag_to_frames('mid_third_wild_type', 'sample_11', LAG_IN_SECONDS)
    arcos_output_df = arcos_wrapper("", arcos_params, arcos_input_df)

    _, heat_map_shuffled, _, _, _, cells_xs, cells_ys, _ = find_arcos_hot_spots(arcos_output_df)
    plot(heat_map_shuffled, cells_xs, cells_ys, is_swapped=True, max_limit=max_limit)

