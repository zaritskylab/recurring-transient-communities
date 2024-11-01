import os
from copy import deepcopy
from typing import List
from pathlib import Path
import pandas as pd
import random
import json
from common.config import LAG_IN_SECONDS
from common.constants import HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION, ARCOS_INPUT_LOCATION, ARCOS_OUTPUT_LOCATION
from common.density import get_activations_per_cell
from common.frame_to_frequency import FrequencyTranslator
from data_layer.arcos.arcos_hotspots_detection import find_arcos_hot_spots
from data_layer.arcos.arcos_params import DEFAULT_PARAMS
from data_layer.arcos.arcos_wrapper import arcos_wrapper
from data_layer.arcos.hotspots_mapper import EXPERIMENTS_HOT_SPOTS, create_swap_candidates, SwapPermutation, Swap


def _swap_cells(df: pd.DataFrame, old_c: int, new_c: int) -> pd.DataFrame:
    """ Swaps the cells in the given DataFrame """

    pd.options.mode.chained_assignment = None  # default='warn'
    untouched_df = df[~df['cell_id'].isin([old_c, new_c])]
    old_cell_df = df[df['cell_id'] == old_c]
    old_cell_df.reset_index(inplace=True, drop=True)
    new_cell_df = df[df['cell_id'] == new_c]
    new_cell_df.reset_index(inplace=True, drop=True)

    old_cell_x = old_cell_df['x_microns']
    old_cell_y = old_cell_df['y_microns']
    new_cell_x = new_cell_df['x_microns']
    new_cell_y = new_cell_df['y_microns']
    old_cell_df['x_microns'] = new_cell_x
    old_cell_df['y_microns'] = new_cell_y
    new_cell_df['x_microns'] = old_cell_x
    new_cell_df['y_microns'] = old_cell_y

    old_cell_x = old_cell_df['x']
    old_cell_y = old_cell_df['y']
    new_cell_x = new_cell_df['x']
    new_cell_y = new_cell_df['y']
    old_cell_df['x'] = new_cell_x
    old_cell_df['y'] = new_cell_y
    new_cell_df['x'] = old_cell_x
    new_cell_df['y'] = old_cell_y

    swapped_df = pd.concat((untouched_df, old_cell_df, new_cell_df))
    swapped_df.reset_index(inplace=True, drop=True)
    pd.options.mode.chained_assignment = 'warn'
    return swapped_df

def _get_events_among_hs_cells(arcos_df: pd.DataFrame, hs_cells: List[int]) -> int:
    """ Returns the amount of events that happened in the given cells """
    events_df = arcos_df[~arcos_df['collid'].isna()]
    hs_cells_set = set(hs_cells)
    res = 0

    for event_id, event_df in events_df.groupby('collid'):
        event_cells = set(event_df['cell_id'].drop_duplicates().to_list())
        if len(event_cells.intersection(hs_cells_set)) >= 2:
            res += 1
    return res

def swapping_same_activity_analysis(arcos_input_df: pd.DataFrame, arcos_output_df: pd.DataFrame, experiment_name: str, experiment_type: str, n_samples: int, lag_in_seconds=None):
    """ Performs the swapping analysis for the hotspots in a given experiment.
    The hotspot cells are swapped based on their activation similarity, so the cells with greater or equal activation are swapped. """

    results = []

    relevant_hotspots = [hotspot for hotspot in EXPERIMENTS_HOT_SPOTS
                      if experiment_type == hotspot.experiment_type and experiment_name == hotspot.experiment_name]

    params = deepcopy(DEFAULT_PARAMS)
    if lag_in_seconds:
        params.nPrev = FrequencyTranslator().temporal_lag_to_frames(experiment_type, experiment_name, lag_in_seconds)

    for experiment_hotspot in relevant_hotspots:
        hotspot_id = experiment_hotspot.id
        activations = get_activations_per_cell(arcos_output_df)
        old_cells = experiment_hotspot.hotspot_cells
        _, heat_map, _, _, _, _, _, orig_mean_events_in_hs = find_arcos_hot_spots(arcos_output_df, check_mean_in_area=old_cells)
        max_limit_for_experiment = int(heat_map.max() + 1)

        # # Creating swap-permutations by activation similarity
        # all_valid_swap_perms, delta_acts_per_swap = create_swap_candidates(old_cells, activations, n_samples)
        #
        # sampled_perms = random.sample(all_valid_swap_perms, n_samples) if len(all_valid_swap_perms) >= n_samples else all_valid_swap_perms
        sampled_perms = {}
        for file in os.listdir(f'data/hotspots_shuffle_analysis/{experiment_type}/{experiment_name}/permutations/'):
            if file.endswith('.json') and file.startswith(f'{hotspot_id}_perm'):
                with open(f'data/hotspots_shuffle_analysis/{experiment_type}/{experiment_name}/permutations/{file}', 'r') as f:
                    import json
                    current_p = json.load(f)
                    sampled_perms[current_p['id']] = SwapPermutation(swaps=[Swap(old_cell=s['old_cell'], new_cell=s['new_cell']) for s in  current_p['swaps']])

        best_perm = None
        best_perm_id = None
        best_perm_decrease = 0
        if len(sampled_perms) == 0:
            results.append({
                'experiment_type': experiment_type,
                'experiment_name': experiment_name,
                'experiment_size': len(activations),
                'hotspot_id': hotspot_id,
                'hotspot_size': len(old_cells),
                'max_events_in_hotspot_area': max_limit_for_experiment-1,
                'permutation_id': 0,
                'total_permutations_tested': len(sampled_perms),
                'original_mean_events_in_hotspot_area': orig_mean_events_in_hs,
                'swapped_mean_events_in_hotspot_area': -1,
                'decrease_in_hotspot_mean_events': -1,
            })
        else:
            for perm_id, permutation in sampled_perms.items():
                new_cells = [swap.new_cell for swap in permutation.swaps]
                swapped_df = deepcopy(arcos_input_df)

                # Swapping the cells according to the permutation
                for swap in permutation.swaps:
                    swapped_df = _swap_cells(swapped_df, swap.old_cell, swap.new_cell)

                swapped_arcos_df = arcos_wrapper("", params, swapped_df)
                _, _, _, _, _, _, _, swap_mean_events_in_hs = find_arcos_hot_spots(swapped_arcos_df, check_mean_in_area=new_cells)

                decrease_in_MEC = round((swap_mean_events_in_hs - orig_mean_events_in_hs) /orig_mean_events_in_hs, 2)
                results.append({
                    'experiment_type': experiment_type,
                    'experiment_name': experiment_name,
                    'experiment_size': len(activations),
                    'hotspot_id': hotspot_id,
                    'hotspot_size': len(old_cells),
                    'max_events_in_hotspot_area': max_limit_for_experiment-1,
                    'permutation_id': perm_id,
                    'total_permutations_tested': len(sampled_perms),
                    'original_mean_events_in_hotspot_area': orig_mean_events_in_hs,
                    'swapped_mean_events_in_hotspot_area': swap_mean_events_in_hs,
                    'decrease_in_hotspot_mean_events': decrease_in_MEC,
                })

                permutation_loc = f"{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/{experiment_type}/{experiment_name}/permutations/"
                # Path(permutation_loc).mkdir(parents=True, exist_ok=True)
                # current_p = {'id': perm_id,
                #              'swaps': [{'old_cell': swap.old_cell, 'new_cell': swap.new_cell} for swap in permutation.swaps]}
                # with open(f'{permutation_loc}/{hotspot_id}_perm_{perm_id}.json', 'w') as f:
                #     import json
                #     json.dump(current_p, f)

                if decrease_in_MEC < best_perm_decrease:
                    best_perm = permutation
                    best_perm_id = perm_id
                    best_perm_decrease = decrease_in_MEC
                if perm_id % 100 == 0:
                    print(f'Finished {perm_id} / {len(sampled_perms)}')

            if best_perm:
                best_json = {'id': best_perm_id,
                             'swaps': [{'old_cell': swap.old_cell, 'new_cell': swap.new_cell} for swap in best_perm.swaps]}
                with open(f'{permutation_loc}/{hotspot_id}_best_perm.json', 'w') as f:
                    import json
                    json.dump(best_json, f)
    return results

def aggregate_data():
    hotspots_df = pd.read_csv(f'{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/all_results.csv')
    hotspots_df['pvalue'] = hotspots_df['swapped_mean_events_in_hotspot_area'] >= hotspots_df[
        'original_mean_events_in_hotspot_area']
    cnt = hotspots_df.groupby(['experiment_type', 'experiment_name', 'hotspot_id'])['experiment_name'].count()
    hotspots_df = hotspots_df.groupby(['experiment_type', 'experiment_name', 'hotspot_id'])[['decrease_in_hotspot_mean_events', 'pvalue']].agg(
        {'decrease_in_hotspot_mean_events': 'mean', 'pvalue': 'sum'}).rename(
        columns={'decrease_in_hotspot_mean_events': '%delta (hotspot swap)'})
    hotspots_df['permutations_count'] = cnt
    hotspots_df['pvalue'] = hotspots_df['pvalue']/hotspots_df['permutations_count']
    hotspots_df.to_csv(f'{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/aggregated_results.csv')

if __name__ == '__main__':
    experiment_types = ['late_second_wild_type', 'early_third_wild_type', 'mid_third_wild_type',
                        'late_second_zpg_RNAi', 'early_third_zpg_RNAi', 'mid_third_zpg_RNAi',
                        'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout']

    swapping_results = []
    for experiment_type in experiment_types:
        for filename in os.listdir(f'{ARCOS_INPUT_LOCATION}/{experiment_type}'):
            if filename.endswith('.csv'):
                experiment_name = filename.split('.csv')[0]
                arcos_input_df = pd.read_csv(f'{ARCOS_INPUT_LOCATION}/{experiment_type}/{filename}')
                arcos_output_df = pd.read_csv(f'{ARCOS_OUTPUT_LOCATION}/{experiment_type}/{filename}')

                print(f"{experiment_name} cells = {arcos_input_df.cell_id.nunique()}/{arcos_output_df.cell_id.nunique()}")
                try:
                    res = swapping_same_activity_analysis(arcos_input_df, arcos_output_df, experiment_name, experiment_type, 1000, lag_in_seconds=LAG_IN_SECONDS)
                except Exception:
                    continue
                if res:
                    swapping_results.append(res)

    flattened_swapping_results = [item for sublist in swapping_results for item in sublist]
    swapping_results_df = pd.DataFrame(data=flattened_swapping_results)
    swapping_results_df.to_csv(f'{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/all_results.csv')
    aggregate_data()
