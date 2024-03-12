from common.config import LAG_IN_SECONDS
from common.constants import ARCOS_INPUT_LOCATION, ARCOS_OUTPUT_LOCATION
from common.frame_to_frequency import FrequencyTranslator
from data_layer.arcos.arcos_params import DEFAULT_PARAMS
from data_layer.arcos.arcos_wrapper import arcos_wrapper

def apply_arcos(experiment_type: str, experiment_name: str):
    frequency_translator = FrequencyTranslator()
    nPrev = frequency_translator.temporal_lag_to_frames(experiment_type, experiment_name, LAG_IN_SECONDS)

    arcos_params = DEFAULT_PARAMS
    arcos_params.nPrev = nPrev

    arcos_res_df = arcos_wrapper(f'{ARCOS_INPUT_LOCATION}{experiment_type}/{experiment_name}.csv', arcos_params)
    arcos_res_df.to_csv(f"{ARCOS_OUTPUT_LOCATION}{experiment_type}/{experiment_name}.csv")
