from dataclasses import dataclass
from typing import Optional
from arcos4py import ARCOS
from arcos4py.tools.filter_events import filterCollev
import pandas as pd


@dataclass
class ArcosParams():
    smoothK: int
    binThr: float
    biasMet: str

    neighborhoodSize: int
    minClsz: int
    nPrev: int

    minDuration: int
    minTotalEventSize: int

    x_col: str = 'x_microns'
    y_col: str = 'y_microns'
    frame_column: str = "frame"
    id_column: str = "cell_id"
    measurement_column: str = "normalized_intensity"
    clid_column: str = "collid"

    peakThr: Optional[float] = None
    biasK: Optional[int] = None


def _arcos_wrapper(input_df: pd.DataFrame, params: ArcosParams, should_binarize=True):
    """ Wrapper for the ARCOS framework. """
    arcos = ARCOS(input_df,
                  posCols=[params.x_col, params.y_col],
                  frame_column=params.frame_column,
                  id_column=params.id_column,
                  measurement_column=params.measurement_column,
                  clid_column=params.clid_column)

    if not should_binarize:
        binarized_df = input_df
    else:
        binarized_df = arcos.bin_measurements(smoothK=params.smoothK,
                                              biasK=params.biasK,
                                              peakThr=params.peakThr,
                                              binThr=params.binThr,
                                              polyDeg=None,
                                              biasMet=params.biasMet)

    coll_events_df = arcos.trackCollev(eps=params.neighborhoodSize,
                                       minClsz=params.minClsz,
                                       nPrev=params.nPrev)

    arcos_filter = filterCollev(data=coll_events_df,
                                frame_column=params.frame_column,
                                obj_id_column=params.id_column)
    filtered_coll_events_df = arcos_filter.filter(coll_duration=params.minDuration,
                                                  coll_total_size=params.minTotalEventSize)

    arcos_final_df = pd.merge(left=binarized_df,right=filtered_coll_events_df, how='left', on=list(binarized_df.columns))

    return arcos_final_df

def arcos_wrapper(input_data_loc: str, params: ArcosParams, input_df: pd.DataFrame=None, should_binarize=True):
    if not (type(input_df) == pd.DataFrame and len(input_df) > 0):
        input_df = pd.read_csv(input_data_loc)
    return _arcos_wrapper(input_df, params, should_binarize=should_binarize)
