import os
from copy import deepcopy
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

from common.constants import ARCOS_OUTPUT_LOCATION, ARCOS_COMMUNITIES_STATS_LOCATION
from common.data_classes import Community
from common.frame_to_frequency import FrequencyTranslator

def _get_closest_cell_id(cell_x, cell_y, current_active_cells_data):
    """ Returns the closest cell id to the given cell coordinates """
    closest_cell_id = np.inf
    for cell_id, cell_data in current_active_cells_data.items():
        dist = np.sqrt((cell_x - cell_data['x_microns']) ** 2 + (cell_y - cell_data['y_microns']) ** 2)
        if dist < closest_cell_id:
            closest_cell_id = cell_id
    return closest_cell_id, dist

def _calc_avg_adjacent_spread_time(df: pd.DataFrame, frequence_in_seconds):
    """ Calculates the average spread time between adjacent cells in the given dataframe"""
    microns_per_second_ratio = []

    df.sort_values(by=['frame'], ascending=True, inplace=True)
    first_frame = df['frame'].iloc[0]

    current_active_cells = set(df[df['frame'] == first_frame]['cell_id'].values.tolist())
    current_active_cells_data = {
        cell_id: {'x_microns': df[(df['frame'] == first_frame) & (df['cell_id'] == cell_id)]['x_microns'].iloc[0],
                  'y_microns': df[(df['frame'] == first_frame) & (df['cell_id'] == cell_id)]['y_microns'].iloc[0],
                  'start_frame': first_frame}
        for cell_id in current_active_cells
    }

    for frame_id, frame_df in df.groupby('frame'):
        if frame_id == first_frame:
            continue

        frame_cells = set(frame_df['cell_id'].values.tolist())
        new_cells = frame_cells.difference(current_active_cells)
        existing_cells = frame_cells.intersection(current_active_cells)

        new_cells_data = {}
        for cell_id in new_cells:
            cell_x = frame_df[frame_df['cell_id'] == cell_id]['x_microns'].iloc[0]
            cell_y = frame_df[frame_df['cell_id'] == cell_id]['y_microns'].iloc[0]

            closest_cell_id, dist = _get_closest_cell_id(cell_x, cell_y, current_active_cells_data)
            new_cells_data[cell_id] = {'x_microns': cell_x,
                                       'y_microns': cell_y,
                                       'start_frame': frame_id}
            microns_per_second_ratio.append(dist / (frequence_in_seconds * (frame_id - current_active_cells_data[closest_cell_id]['start_frame'])))

        # Updating the current active cells
        current_active_cells = set([cell_id for cell_id in current_active_cells if cell_id in existing_cells])
        current_active_cells.update(new_cells)
        current_active_cells_data = {k:v for k,v in current_active_cells_data.items() if k in existing_cells}
        current_active_cells_data = {**current_active_cells_data, **new_cells_data}

    if microns_per_second_ratio == []:
        return 0  # all cells shine at once (first frame)
    return np.mean(microns_per_second_ratio)

def _calc_cells_activated_amount(collid_df: pd.DataFrame):
    current_cells = set()
    amount = 0
    for frame, frame_df in collid_df.groupby('frame'):
        cells_in_frame = frame_df['cell_id'].unique()
        new_cells = set(cells_in_frame).difference(current_cells)
        amount += len(new_cells)
        current_cells = cells_in_frame
    return amount
def _get_community_by_id(community_id: int, all_communities: List[Community]) -> Community:
    for community in all_communities:
        if community.community_id == community_id:
            return community
    return None

def _get_collid_to_community_mapping(df: pd.DataFrame):
    all_communities = []
    i = 0
    for collid, collid_df in df.groupby('collid'):
        all_communities.append(Community(
            cells=set(collid_df['cell_id'].unique()),
            community_id=i,
            community_size=int(len(collid_df['cell_id'].unique())),
            original_coll_ids=[int(collid)],
            coll_ids=[int(collid)],
        ))
        i += 1

    all_communities = sorted(all_communities, key=lambda x: x.community_id, reverse=False)

    community_to_collids = {community.community_id: {'original': community.original_coll_ids,
                                                     'additional': [i for i in community.coll_ids if i not in community.original_coll_ids]}
                            for community in all_communities}
    collid_to_communities = {}
    for collid in df[~df['collid'].isna()]['collid'].unique():
        collid_int = int(collid)
        collid_to_communities[collid_int] = {'original': -1, 'additional': []}
        for community in all_communities:
            if collid_int in community.original_coll_ids:
                collid_to_communities[collid_int]['original'] = community.community_id
            elif collid_int in community.coll_ids:
                collid_to_communities[collid_int]['additional'].append(community.community_id)

    return community_to_collids, collid_to_communities, all_communities


def write_experiment_communities_statistics_to_csv(all_communities, experiment_type, experiment_name):
    df = pd.DataFrame()
    for community in all_communities:
        df = df.append(community.to_dict(), ignore_index=True)

    target_folder = f'{ARCOS_COMMUNITIES_STATS_LOCATION}/{experiment_type}/'
    if not os.path.exists(target_folder):
        Path(f'{target_folder}').mkdir(parents=True, exist_ok=True)

    df.to_csv(target_folder + f'{experiment_name}_communities_statistics.csv')


def generate_per_sample_community_statistics(experiment_types: List[str]):
    freq_translator = FrequencyTranslator()
    all_experiments_df = pd.DataFrame()

    for experiment_type in experiment_types:
        experiment_outputs = os.listdir(f"{ARCOS_OUTPUT_LOCATION}/{experiment_type}")

        for arcos_output in experiment_outputs:
            statistics_df = pd.DataFrame()
            df = pd.read_csv(f"{ARCOS_OUTPUT_LOCATION}/{experiment_type}/{arcos_output}")
            experiment_name = arcos_output.replace('.csv', '')

            experiment_stats = {
                'experiment_type': experiment_type,
                'experiment_name': experiment_name,
                'experiment_size': df['cell_id'].nunique(),
                'experiment_length_in_minutes': freq_translator.get_movie_in_minutes(experiment_type, experiment_name),
                'num_of_frames': df['frame'].nunique(),
                'frames_frequency (sec)': freq_translator.get_fps(experiment_type, experiment_name),
                'num_of_arcos_events': df['collid'].nunique(),
            }

            community_to_collids, collid_to_communities, all_communities = _get_collid_to_community_mapping(df)

            if len(all_communities) == 0:
                no_communities_record = {
                    'community_id': -1,
                    'additional_communities': [],
                    'community_size': 0,
                    'event_cells': [],
                    'event_start_cell': [],
                    'event_end_cell': [],
                    'event_length_in_seconds': 0,
                    'activated_cells': 0,
                    'info_spread_rate (avg between 2 adjacent cells) - cells per second': 0,
                }
                no_communities_record.update(deepcopy(experiment_stats))
                statistics_df = pd.concat((statistics_df, pd.DataFrame([no_communities_record])), ignore_index=True)

            for collid, collid_df in df.groupby('collid'):
                collid_int = int(collid)
                experiment_stats_copy = deepcopy(experiment_stats)
                community_id, additional_communities_ids = collid_to_communities[collid_int]['original'], collid_to_communities[collid_int]['additional']
                community = _get_community_by_id(community_id, all_communities)

                experiment_stats_copy['community_id'] = community.community_id
                experiment_stats_copy['additional_communities'] = additional_communities_ids
                experiment_stats_copy['community_size'] = community.community_size

                experiment_stats_copy['event_cells'] = collid_df['cell_id'].unique()
                experiment_stats_copy['event_start_cell'] = collid_df[collid_df['frame'] == collid_df['frame'].min()]['cell_id'].unique()
                experiment_stats_copy['event_end_cell'] = collid_df[collid_df['frame'] == collid_df['frame'].max()]['cell_id'].unique()
                experiment_stats_copy['event_length_in_seconds'] = collid_df['frame'].nunique() * experiment_stats['frames_frequency (sec)']
                experiment_stats_copy['activated_cells'] = _calc_cells_activated_amount(collid_df)
                experiment_stats_copy['info_spread_rate (avg between 2 adjacent cells) - cells per second'] = _calc_avg_adjacent_spread_time(collid_df, experiment_stats['frames_frequency (sec)'])

                statistics_df = pd.concat([statistics_df, pd.DataFrame([experiment_stats_copy])], ignore_index=True)

            write_experiment_communities_statistics_to_csv(all_communities, experiment_type, experiment_name)
            all_experiments_df = pd.concat((all_experiments_df, statistics_df), ignore_index=True)

    all_experiments_df.to_csv('./data/cross_experiment_communities_statistics.csv', index=False)

