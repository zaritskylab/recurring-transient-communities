from typing import List
import pandas as pd
import os
from common.constants import EXPERIMENTS_LOCATION
from common.euclidean_distance import EuclideanMatrixSingleton
from pathlib import Path

def _get_cells_data_from_csv(folder: str, num_of_cells: int):
    """ Returns a list of dataframes, each dataframe represents a cell time series"""
    cells_data : List[pd.DataFrame] = []

    for i in range(num_of_cells):
        output_location = f"""{folder}/cell_{i}.csv"""
        cell = pd.read_csv(output_location)
        cell.drop(columns=['Unnamed: 0'], inplace=True)
        cell["cell_num"] = i
        cells_data.append(cell)
    return cells_data

def get_num_of_cells(df: pd.DataFrame):
    return round(len(df.columns)/3)

def _normalize(Ft, Fmin, Fmax):
    """ Min-Max normalization of the intensity of a cell at time T """
    return (Ft - Fmin)/(Fmax - Fmin)


def _extract(raw_data, num_of_cells: int):
    """ Extracts the data of each cell from the raw data """
    cells_data: List[pd.DataFrame] = []
    for i in range(num_of_cells):
        curr_cell_df: pd.DataFrame = raw_data.iloc[:, (1 + (3 * i)):(4 + (3 * i))]
        curr_cell_cols = list(curr_cell_df.columns.values)
        mean_intensity_col = curr_cell_cols[0]
        x_col = curr_cell_cols[1]
        y_col = curr_cell_cols[2]
        curr_cell_df.rename(columns={mean_intensity_col: "mean_intensity", x_col: "x", y_col: "y"},
                            inplace=True)
        cells_data.append(curr_cell_df)
    return cells_data

def _handle_missing_data(cell_df: pd.DataFrame):
    """ Handles missing data by filling it with the mean of the previous and next values """
    cell_df[cell_df["mean_intensity"] == 0] = pd.NA

    back_filled = cell_df.bfill()
    forward_filled = cell_df.ffill()

    back_filled[back_filled.isnull().any(axis=1)] = forward_filled[back_filled.isnull().any(axis=1)]
    forward_filled[forward_filled.isnull().any(axis=1)] = back_filled[forward_filled.isnull().any(axis=1)]

    back_filled["mean_intensity"] = pd.to_numeric(back_filled["mean_intensity"])
    back_filled["x"] = pd.to_numeric(back_filled["x"])
    back_filled["y"] = pd.to_numeric(back_filled["y"])
    forward_filled["mean_intensity"] = pd.to_numeric(forward_filled["mean_intensity"])
    forward_filled["x"] = pd.to_numeric(forward_filled["x"])
    forward_filled["y"] = pd.to_numeric(forward_filled["y"])
    return pd.concat([forward_filled, back_filled]).groupby(level=0).mean()

def _get_key_feature(cells_data: List[pd.DataFrame]):
    """ Calculates the key features of all cells """
    features_cols = ["min_intensity", "max_intensity", "mean_intensity", "normalized_mean", "median_intensity",
                     "normalized_median", "0.9_quantile", "normalized_0.9_quantile"]
    cells_features = pd.DataFrame(columns=features_cols)
    for cell in cells_data:
        min = cell[["mean_intensity"]].min()
        max = cell[["mean_intensity"]].max()
        mean = cell[["mean_intensity"]].mean()
        normalized_mean = _normalize(mean, min, max)
        median = cell[["mean_intensity"]].median()
        normalized_median = _normalize(median, min, max)
        quantile = pd.Series(index=["mean_intensity"], data=cell["mean_intensity"].quantile(q=0.9))
        normalized_quantile = _normalize(quantile, min, max)

        cell[['normalized_intensity']] = cell['mean_intensity'].apply(lambda Ft: _normalize(Ft, min, max))

        cell_features_df: pd.DataFrame = pd.concat(
            [min, max, mean, normalized_mean, median, normalized_median, quantile, normalized_quantile], axis=1)
        curr_cols = list(cell_features_df.columns.values)
        mapper = dict(zip(curr_cols, features_cols))
        cell_features_df.rename(columns=mapper, inplace=True)
        cells_features = pd.concat((cells_features, cell_features_df))

    cells_features.reset_index(inplace=True)
    return cells_features

def extract_clean_and_normalize(path: str):
    dirs = list(filter(lambda dir: dir.endswith('.csv'), os.listdir(path)))

    for filename in dirs:
        # if re.search("\.csv$", filename):
        raw_data: pd.DataFrame = pd.read_csv(f"{path}/{filename}")
        num_of_cells = get_num_of_cells(raw_data)

        print(f"Extracting each cell's data from: {path}/{filename}")
        cells_data = _extract(raw_data, num_of_cells)

        print("Handling missing data")
        for i in range(num_of_cells):
            cell_df = cells_data[i]
            cells_data[i] = _handle_missing_data(cell_df)

        print("Remove duplicate cells readouts (2 cells that represent one)")
        euclidean_matrix = EuclideanMatrixSingleton(cells_data).get_distance_matrix()
        cell_indices_to_remove = []
        for i in range(num_of_cells):
            for j in range(i+1, num_of_cells):
                if euclidean_matrix[i][j] == 0.0:
                    cell_indices_to_remove.append(j)

        num_of_cells -= len(cell_indices_to_remove)
        for index in sorted(cell_indices_to_remove, reverse=True):
            del cells_data[index]

        print("Calculating Key Features Per Cell")
        cells_features = _get_key_feature(cells_data)

        print(f"Saving Data To {num_of_cells} CSVs + 1 for all Features")
        filename_without_extension = filename.replace('.csv', '')
        output_location = path.replace('raw', 'processed')+f'/{filename_without_extension}'
        Path(f'{output_location}').mkdir(parents=True, exist_ok=True)

        for i in range(num_of_cells):
            name = f"cell_{i}"
            cells_data[i].to_csv(f'{output_location}/{name}.csv')
        cells_features.to_csv(f"{output_location}/All_Features.csv")



# def write_append_to_csv(analysis_folder, df, file_name, should_create_new_file):
#     if analysis_folder is None: #for unified CSVs
#         loc = f'{PROJECT_LOCATION}{file_name}.csv'
#     else:
#         loc = f'{PROJECT_LOCATION}{analysis_folder}/{file_name}.csv'
#
#     if should_create_new_file or (not os.path.isfile(loc)):
#         df.to_csv(loc, index=False, header='column_names')
#     else:  # else it exists so append without writing the header
#         df.to_csv(loc, index=False, mode='a', header=False)


# def generate_corrs_data_for_experiment_type(experiment_type, lags, cutoffs):
#     loc = f'{EXPERIMENTS_LOCATION}{experiment_type}/'
#     data_folder_subdirs = [name for name in os.listdir(loc) if os.path.isdir(os.path.join(loc, name))]
#
#     for data_folder_subdir in data_folder_subdirs:
#         experiment_data_folder = f'{loc}{data_folder_subdir}/data/'
#         num_of_cells = int(data_folder_subdir.split('_')[0])
#
#         print(f"\n>>>>>> START WORKING ON DIR: {data_folder_subdir} >>>>>>")
#         cells_data = _get_cells_data_from_csv(experiment_data_folder, num_of_cells)
#         all_data = pd.concat(cells_data, axis=0)
#         all_data.set_index(['cell_num'], append=True, inplace=True)
#         all_data.index.rename(names=["sample", "cell_num"], inplace=True)
#         time_series_len = len(cells_data[0])
#
#         logic = BaseLogic(cells_data, lags, cutoffs, experiment_data_folder)
#         for lag in lags:
#             for cutoff in cutoffs:
#                 corrs_passed_cutoff_by_lag, out, pairs_edges, total_pairs_passed_cutoff_by_rounded_distance, \
#                 pairs_passed_cutoff_max_corr_and_distance_and_lags, pairs_corr_at_fixed_lag_above_p_val = logic.xcorr_w_lags(cells_data, num_of_cells,
#                                                                                         time_series_len, lag,
#                                                                                         cutoff=cutoff)
#
#                 xcorr_data_tuple = (corrs_passed_cutoff_by_lag, out, pairs_edges, total_pairs_passed_cutoff_by_rounded_distance,
#                                     pairs_passed_cutoff_max_corr_and_distance_and_lags, pairs_corr_at_fixed_lag_above_p_val)
#
#                 with open(f"{experiment_data_folder}xcorr_lag_{lag}_cutoff_{cutoff}.pkl", 'wb+') as pkl:
#                     pickle.dump(xcorr_data_tuple, pkl, pickle.HIGHEST_PROTOCOL)
#
#         print(f"<<<<<< FINISHED WORKING ON DIR: {data_folder_subdir} <<<<<<")
#

if __name__ == '__main__':
    raw_data_loc = '../data/raw/'
    for experiment_folder in os.listdir(raw_data_loc):
        print(f"Extracting data from {experiment_folder}")
        extract_clean_and_normalize(raw_data_loc + experiment_folder)
        print(f"Finished extracting data from {experiment_folder}")
        print("------------------------------------------------------")