import json

import numpy as np
import pandas as pd

from common.constants import EXPERIMENTS_LOCATION, ARCOS_OUTPUT_LOCATION
from common.density import get_cell_density, get_activations_per_cell
from common.frame_to_frequency import FrequencyTranslator
from common.significance_calculator import SignificanceMethod, SignificanceCalculator


def generate_cross_experiment_spatial_significance_data():
    significance_method = SignificanceMethod.BOOTSTRAP

    observed_data_loc = f'{EXPERIMENTS_LOCATION}/arcos_shuffle_analysis/unified_observed_data.json'
    in_silico_data_loc = f'{EXPERIMENTS_LOCATION}/arcos_shuffle_analysis/unified_in_silico_data.json'

    with open(observed_data_loc,'r') as f:
        observed_data = json.load(f)
        observed_data = observed_data
    with open(in_silico_data_loc,'r') as f:
        in_silico_data = json.load(f)
        in_silico_data = in_silico_data

    frequency_translator = FrequencyTranslator()
    experiment_types = ['mid_third_wild_type', 'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout',
                        'mid_third_zpg_RNAi',
                        'late_second_wild_type', 'early_third_wild_type',
                        'late_second_zpg_RNAi','early_third_zpg_RNAi'
                        ]
    data = {
        'experiment_type': [],
        'p_value': [],
        'mean_density': [],
        'mean_activations': [],
        'mean_activations_per_minute': [],
        'mean_events_per_cell': [],
        'MEC_magnitude': [],
        'MEC_per_minute': [],
    }

    si_data = {
        'experiment_name': [],
        'MEC': [],
        'simulated_MEC': [],
        'pvalue': [],
        'magnitude': []
    }

    for experiment_type in experiment_types:
        experiment_names = list(observed_data[experiment_type].keys())

        for experiment_name in experiment_names:
            # Calculates Mean #Events per Cell
            actual_val = np.mean(observed_data[experiment_type][experiment_name]['events_per_cell'])
            shuffled_vals = [v['events_per_cell'] for i,v in in_silico_data[experiment_type][experiment_name].items()]
            shuffled_vals = [np.mean(val) for val in shuffled_vals]

            p_value = SignificanceCalculator().calculate(significance_method, actual_val, shuffled_vals)

            arcos_df = pd.read_csv(f"{ARCOS_OUTPUT_LOCATION}{experiment_type}/{experiment_name}.csv")
            mean_density = np.mean(get_cell_density(arcos_df))
            mean_activations = np.mean(get_activations_per_cell(arcos_df))
            MEC = np.mean(observed_data[experiment_type][experiment_name]['events_per_cell'])
            MEC_per_minute = MEC / frequency_translator.get_movie_in_minutes(experiment_type, experiment_name)

            simulated_MEC = np.array([np.mean(sim['events_per_cell']) for sim in in_silico_data[experiment_type][experiment_name].values()])
            simulated_MEC_per_minute = simulated_MEC / frequency_translator.get_movie_in_minutes(experiment_type, experiment_name)
            magnitude = MEC / np.mean(simulated_MEC)


            data['experiment_type'].append(experiment_type)
            data['p_value'].append(p_value)
            data['mean_density'].append(mean_density)
            data['mean_activations'].append(mean_activations)
            data['mean_activations_per_minute'].append(mean_activations / frequency_translator.get_movie_in_minutes(experiment_type, experiment_name))
            data['mean_events_per_cell'].append(MEC)
            data['MEC_magnitude'].append(magnitude)
            data['MEC_per_minute'].append(MEC_per_minute)

            if experiment_type == 'mid_third_wild_type':
                si_data['experiment_name'].append(experiment_name)
                si_data['MEC'].append(MEC)
                si_data['simulated_MEC'].append(simulated_MEC)
                si_data['pvalue'].append(p_value)
                si_data['magnitude'].append(magnitude)


    df = pd.DataFrame(data)
    df.to_csv('./data/cross_experiment_arcos_agg_data.csv', index=False)
