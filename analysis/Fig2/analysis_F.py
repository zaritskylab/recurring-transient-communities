import pandas as pd
import matplotlib.pyplot as plt
from common.constants import HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION, ARCOS_OUTPUT_LOCATION, FIGURES_LOCATION
from common.frame_to_frequency import FrequencyTranslator
from data_layer.arcos.hotspots_mapper import EXPERIMENTS_HOT_SPOTS


def plot(df: pd.DataFrame):

    experiment_type = 'mid_third_wild_type'

    field_required = 'interaction_with_external_cells_ratio'
    fig, ax = plt.subplots(figsize=(8, 4))

    means = []
    experiment_df = df[df['experiment_type'] == experiment_type]
    if len(experiment_df) > 0:
        y = experiment_df['hotspot_size']
        x = experiment_df[field_required]
        means.append(y.mean())
        ax.scatter(x, y, label=experiment_type)

        ax.set_xlim(0, 1.1)

    plt.xticks([0, 0.5, 1], fontsize=14)
    plt.xlabel('External Interaction Probability (%)', fontsize=18)


    plt.ylabel('Hotspot Size (#Cells)', fontsize=18)
    plt.ylim(0, 16)
    plt.yticks([0, 8, 16], fontsize=14, rotation=90)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}Fig_2_F.png')
    plt.close(fig)


def analyze():

    data = []
    agg_df = pd.read_csv(f'{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}/aggregated_results.csv')

    # Calculating the p-values based on the is_significant field
    p_val_field = 'is_significant'
    agg_df[p_val_field] = agg_df.apply(lambda row: 1 if row['pvalue'] < 0.05 and row['permutations_count'] > 100 else 0, axis=1)
    agg_df['enough_statistics'] = agg_df.apply(lambda row: 1 if row['permutations_count'] > 100 else 0, axis=1)
    agg_df = pd.merge(agg_df, FrequencyTranslator().get_names_df(), on=['experiment_type','experiment_name'], how='inner')

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

        start_in_hs_and_stayed = 0
        start_in_hs_and_propagated = 0
        inbound = 0
        mixed_start = 0
        hs_proportion = []
        non_hs_total_cells = []
        first_cells_in_hotspot = []

        agg_df_record = agg_df[(agg_df['experiment_name'] == experiment_hotspot.experiment_name) & (agg_df['hotspot_id'] == experiment_hotspot.id)]
        if len(agg_df_record) == 1 and int(agg_df_record['is_significant']) == 1:
            for coll_event in relevant_coll_events:
                coll_event_df = arcos_df[arcos_df['collid'] == coll_event]
                coll_event_cells = set(coll_event_df['cell_id'].drop_duplicates().values.tolist())

                first_frame = coll_event_df['frame'].min()
                first_cells = set(coll_event_df[coll_event_df['frame'] == first_frame]['cell_id'].drop_duplicates().values.tolist())
                first_cells_in_hotspot.extend(list(first_cells))

                if len(first_cells.intersection(hotspot_cells)) == 0:
                    inbound += 1
                elif 0 < len(first_cells.intersection(hotspot_cells)) < len(first_cells):
                    mixed_start += 1
                elif 0 < len(first_cells.intersection(hotspot_cells)) and len(coll_event_cells.intersection(hotspot_cells)) == len(coll_event_cells):
                    start_in_hs_and_stayed += 1
                else:
                    start_in_hs_and_propagated += 1

                coll_event_hs_cells = coll_event_cells.intersection(hotspot_cells)
                hs_proportion.append(len(coll_event_hs_cells)/ len(coll_event_cells))
                non_hs_total_cells.append(coll_event_cells.difference(coll_event_hs_cells))
                global_hs_prop.append(len(coll_event_hs_cells))
                global_non_hs_prop.append(len(coll_event_cells.difference(coll_event_hs_cells)))

            env_interaction_ratio = (mixed_start + inbound + start_in_hs_and_propagated) / (mixed_start + inbound + start_in_hs_and_propagated + start_in_hs_and_stayed)
            reduced_non_hs_total_cells = set()
            for s in non_hs_total_cells:
                reduced_non_hs_total_cells = reduced_non_hs_total_cells.union(s)

            data.append({
                'experiment_type': experiment_hotspot.experiment_type,
                'experiment_name': experiment_hotspot.experiment_name,
                'hotspot_id': experiment_hotspot.id,
                'hotspot_size': len(experiment_hotspot.hotspot_cells),
                'interaction_with_external_cells_ratio': env_interaction_ratio,
            })

    df = pd.DataFrame(data)
    plot(df)

def main():
    analyze()

