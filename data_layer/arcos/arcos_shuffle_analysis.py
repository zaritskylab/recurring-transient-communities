import json
import os
from pathlib import Path
from typing import Dict, List
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon
from common.config import LAG_IN_SECONDS
from common.frame_to_frequency import FrequencyTranslator
from data_layer.arcos.arcos_hotspots_detection import find_arcos_hot_spots
from data_layer.arcos.arcos_params import DEFAULT_PARAMS
from common.constants import ARCOS_INPUT_LOCATION, ARCOS_OUTPUT_LOCATION, ARCOS_SHUFFLE_ANALYSIS_LOCATION
from data_layer.arcos.arcos_wrapper import arcos_wrapper


def _resolve_target_file(index: int = None):
    if index != None:
        return f'shuffled_data_{index}.json'
    else:
        return f'actual_data.json'

def _dump_to_json(file_path: str, experiment_type: str, experiment: str, data: Dict, index: int = None):
    """ Dumps the given data to the given file path """
    experiment_data_folder = '/'.join(file_path.split('/')[:-1])
    if not os.path.exists(experiment_data_folder):
        Path(f'{experiment_data_folder}').mkdir(parents=True, exist_ok=True)

    if not os.path.exists(file_path):
        with open(file_path, 'w+') as f:
            if index != None:
                data = {experiment_type: {experiment: {index: data}}}
            else:
                data = {experiment_type: {experiment: data}}
            json.dump(data, f)

    else:
        with open(file_path, 'r+') as f:
            new_data = json.load(f)
            if experiment_type not in new_data:
                new_data[experiment_type] = {}

            if index != None:
                new_data[experiment_type].setdefault(experiment, {})[index] = data
            else:
                new_data[experiment_type][experiment] = data

            f.seek(0)
            json.dump(new_data, f)
            f.truncate()


def _write_hotspot_results(experiment_type, experiment, polygons, heat_map, all_events_df, num_of_cells, distinct_cells,
                             cells_xs, cells_ys, index=None):
    """ Writes the hot spots results to a JSON file"""
    areas = [poly.area for poly in polygons]
    amnt_polys_with_max_events = 0
    for poly in polygons:
        for x, y in zip(np.where(heat_map == heat_map.max())[0], np.where(heat_map == heat_map.max())[1]):
            if poly.contains(Point(x, y)):
                amnt_polys_with_max_events += 1
                break

    events_sizes = all_events_df[['cell_id', 'collid']].drop_duplicates() \
        .groupby('collid').count() \
        .sort_values('cell_id', ascending=False)['cell_id'].to_list()
    events_per_cell_id = all_events_df[['cell_id', 'collid']].drop_duplicates().groupby('cell_id')['collid'].count()

    data = {
        '#cells': num_of_cells,
        '#active_cells': len(distinct_cells),
        'experiment_mm_size': Polygon([Point(cells_xs[i], cells_ys[i]) for i in range(num_of_cells - 1)]).area,
        '#polys': len(polygons),
        'min_poly_area': min(areas) if len(areas) > 0 else 0,
        'max_poly_area': max(areas) if len(areas) > 0 else 0,
        'mean_poly_area': np.mean(areas) if len(areas) > 0 else 0,
        'std_poly_area': np.std(areas) if len(areas) > 0 else 0.00000001,
        'events_sizes': events_sizes,
        'events_per_cell': [int(events_per_cell_id[cell_id]) if cell_id in events_per_cell_id else 0 for cell_id in
                            range(num_of_cells)],
        'maximum_events_in_area': heat_map.max(),
        '#polys_reached_max_events': amnt_polys_with_max_events,
    }
    file_path = f'{ARCOS_SHUFFLE_ANALYSIS_LOCATION}/{experiment_type}/{experiment}/{_resolve_target_file(index)}'
    _dump_to_json(file_path, experiment_type, experiment, data, index)



def _generate_distinct_permutations(l: List, n: int = 1000):
    """ Generates N distinct permutations of the given list """
    mylist = np.array(l)
    perms = set()
    iters = 0
    for i in range(n):  # (1) Draw N samples from permutations Universe U (#U = k!)
        while True:  # (2) Endless loop
            iters += 1
            perm = np.random.permutation(mylist)  # (3) Generate a random permutation form U
            is_switch_places = sum([x!=ix for ix,x in enumerate(perm)]) == len(perm)
            key = tuple(perm)
            if key not in perms and is_switch_places:  # (4) Check if permutation already has been drawn (hash table)
                perms.add(key)  # (5) Insert into set
                break  # (6) Break the endless loop
    return [list(p) for p in perms]


def generate_significance_data_for_experiment(experiment_type, experiment_name, gen_amnt, lag_in_seconds):
    """ Generates the significance data for the given experiment by shuffling the data. The data is shuffled N times,
    and for each shuffle, the hot spots are detected and written to a JSON file."""
    # Get data from ACTUAL experiment ARCOS output, and write to JSON
    original_arcos_input_df = pd.read_csv(f'{ARCOS_INPUT_LOCATION}{experiment_type}/{experiment_name}.csv')
    original_arcos_output_df = pd.read_csv(f'{ARCOS_OUTPUT_LOCATION}{experiment_type}/{experiment_name}.csv')

    polygons, heat_map, all_events_df, num_of_cells, distinct_cells, cells_xs, cells_ys = \
        find_arcos_hot_spots(original_arcos_output_df)

    _write_hotspot_results(experiment_type, experiment_name, polygons,
                             heat_map, all_events_df, num_of_cells, distinct_cells,
                             cells_xs, cells_ys)

    # Get data from SHUFFLED experiments ARCOS output, and write to JSON
    cell_ids = list(original_arcos_input_df['cell_id'].unique())
    permutations = _generate_distinct_permutations(cell_ids, n=gen_amnt)

    params = DEFAULT_PARAMS
    if lag_in_seconds:
        params.nPrev = FrequencyTranslator().temporal_lag_to_frames(experiment_type, experiment_name, lag_in_seconds)

    for i, perm in enumerate(permutations):
        shuffled_arcos_input_df = pd.DataFrame()

        for ix, shuffled_ix in enumerate(perm):
            ix_cols = ['frame', 'mean_intensity', 'normalized_intensity', 'cell_id']
            shuffled_ix_cols = ['x', 'y', 'x_microns', 'y_microns', 'x_pixels', 'y_pixels']

            tmp_df = pd.concat((original_arcos_input_df[original_arcos_input_df['cell_id']==ix][ix_cols].reset_index(drop=True),
                                original_arcos_input_df[original_arcos_input_df['cell_id']==shuffled_ix][shuffled_ix_cols].reset_index(drop=True)),
                               axis=1)
            shuffled_arcos_input_df = pd.concat((shuffled_arcos_input_df, tmp_df))
        shuffled_arcos_input_df.reset_index(inplace=True, drop=True)

        shuffled_arcos_output_df = arcos_wrapper("", params, shuffled_arcos_input_df)

        polygons, heat_map, all_events_df, num_of_cells, distinct_cells, cells_xs, cells_ys = \
            find_arcos_hot_spots(shuffled_arcos_output_df)


        _write_hotspot_results(experiment_type, experiment_name, polygons, heat_map,
                                 all_events_df, num_of_cells, distinct_cells,
                                 cells_xs, cells_ys, index=i)

        if i % 100 == 0:
            print(f"Finished index {i}")


def unify_significance_data():

    unified_actual_data = {}
    unified_shuffled_data = {}

    unified_actual_data_filename = 'unified_observed_data.json'
    unified_shuffled_data_filename = 'unified_in_silico_data.json'

    for dir in os.listdir(ARCOS_SHUFFLE_ANALYSIS_LOCATION):
        experiment_type = dir
        path = f"{ARCOS_SHUFFLE_ANALYSIS_LOCATION}/{dir}"
        if os.path.isdir(path):
            sub_dirs = os.listdir(path)
            for sub_dir in sub_dirs:
                if os.path.isdir(f'{path}/{sub_dir}'):

                    curr_actual_data = f'{path}/{sub_dir}/actual_data.json'
                    with open(curr_actual_data, 'r') as f:
                        curr_actual = json.load(f)
                    if experiment_type not in unified_actual_data.keys():
                        unified_actual_data.update(curr_actual)
                    else:
                        unified_actual_data[experiment_type].update(curr_actual[experiment_type])


                    all_shuffle_files = [f for f in os.listdir(f'{path}/{sub_dir}') if f.startswith('shuffled_data_')]
                    for shuffle_file in all_shuffle_files:
                        with open(f'{path}/{sub_dir}/{shuffle_file}', 'r') as f:
                            curr_shuff = json.load(f)
                        if experiment_type not in unified_shuffled_data.keys():
                            unified_shuffled_data.update(curr_shuff)
                        elif sub_dir not in unified_shuffled_data[experiment_type].keys():
                            unified_shuffled_data[experiment_type][sub_dir] = curr_shuff[experiment_type][sub_dir]
                        else:
                            unified_shuffled_data[experiment_type][sub_dir].update(curr_shuff[experiment_type][sub_dir])

    with open(f'{ARCOS_SHUFFLE_ANALYSIS_LOCATION}/{unified_actual_data_filename}', 'w') as f:
        json.dump(unified_actual_data, f)
    with open(f'{ARCOS_SHUFFLE_ANALYSIS_LOCATION}/{unified_shuffled_data_filename}', 'w') as f:
        json.dump(unified_shuffled_data, f)


def main():
    for experiment_type in os.listdir(ARCOS_INPUT_LOCATION):
        for experiment_name_file in os.listdir(f'{ARCOS_INPUT_LOCATION}/{experiment_type}'):
            experiment_name = experiment_name_file.replace('.csv', '')
            generate_significance_data_for_experiment(experiment_type, experiment_name, 1000, LAG_IN_SECONDS)

    unify_significance_data()