import math
import os
from pathlib import Path
from typing import List
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

from common.constants import HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION


def _get_mean_events_in_area(arcos_df: pd.DataFrame, heat_map, area_of_interest: List[int]):
    from shapely.geometry import Point, Polygon
    area_points = []
    for cell_id in area_of_interest:
        cell_df = arcos_df[arcos_df['cell_id'] == cell_id]
        cell_x = cell_df['x'].mean()
        cell_y = cell_df['y'].mean()
        area_points.append(
            Point(math.floor(cell_x), math.floor(cell_y))
        )
    area_polly = Polygon([[p.x, p.y] for p in area_points]).convex_hull

    area_events = []
    min_x, min_y, max_x, max_y = area_polly.bounds
    for i in range(int(min_x), int(max_x + 1)):
        for j in range(int(min_y), int(max_y + 1)):
            tmp_point = Point(i, j)
            if tmp_point.within(area_polly):
                area_events.append(heat_map[i][j])

    return np.mean(area_events)

def find_arcos_hot_spots(arcos_df: pd.DataFrame, check_mean_in_area: List[int] = None, dump_loc: str = None):

    num_of_cells = arcos_df.cell_id.nunique()
    cells_xs = []
    cells_ys = []
    polygons = []
    distinct_cells = set()

    for i in range(arcos_df['cell_id'].max() + 1):
        current_cell = arcos_df[arcos_df['cell_id'] == i]
        cells_xs.append(current_cell['x'].mean())
        cells_ys.append(current_cell['y'].mean())

    ## Create padded heatmap, each cell is 1micron X 1micron
    heat_map = np.zeros(shape=(math.ceil(max(arcos_df['x']) + 10), math.ceil(max(arcos_df['y']) + 10)))

    all_events_df = arcos_df[~arcos_df['collid'].isna()]
    all_events_id = all_events_df['collid'].drop_duplicates().to_list()
    for event_id in all_events_id:
        event_df = all_events_df[all_events_df['collid'] == event_id]

        distinct_participating_cells = event_df['cell_id'].drop_duplicates().to_list()

        cells_coordinates = []
        for cell_id in distinct_participating_cells:
            distinct_cells.add(cell_id)
            cell_df = event_df[event_df['cell_id'] == cell_id]
            cell_x = cell_df['x'].mean()
            cell_y = cell_df['y'].mean()
            cells_coordinates.append(Point(math.floor(cell_x), math.floor(cell_y)))

        try:
            if len(cells_coordinates) < 3:
                continue

            else:
                convex_poly = Polygon([[p.x, p.y] for p in cells_coordinates]).convex_hull
                polygons.append(convex_poly)

                min_x, min_y, max_x, max_y = convex_poly.bounds
                for i in range(int(min_x), int(max_x + 1)):
                    for j in range(int(min_y), int(max_y + 1)):
                        tmp_point = Point(i, j)
                        if tmp_point.within(convex_poly):
                            heat_map[i][j] += 1

        except Exception as e:
            print("exception")

    # Contours
    max_events = max(heat_map.flatten())
    min_events_requirement = 5
    if max_events >= min_events_requirement:
        events_threshold = max(max_events // 2, min_events_requirement)
        bin_heat_map = heat_map >= events_threshold
        contour_heat_map = np.zeros(bin_heat_map.shape, dtype=bool)

        for i in range(1, len(bin_heat_map) - 1):
            for j in range(1, len(bin_heat_map.T) - 1):
                if bin_heat_map[i][j] == False:
                    continue

                else:
                    contour_heat_map[i][j] = False in [
                        bin_heat_map[i - 1][j - 1],
                        bin_heat_map[i - 1][j],
                        bin_heat_map[i - 1][j + 1],
                        bin_heat_map[i][j + 1],
                        bin_heat_map[i + 1][j + 1],
                        bin_heat_map[i + 1][j],
                        bin_heat_map[i + 1][j - 1],
                        bin_heat_map[i][j - 1]
                    ]
    area_mean = None
    if check_mean_in_area:
        area_mean = _get_mean_events_in_area(arcos_df, heat_map, check_mean_in_area)
    if dump_loc:
        target_folder = f'{HOTSPOTS_SHUFFLE_ANALYSIS_LOCATION}{dump_loc}'
        if not os.path.exists(target_folder):
            Path(f'{target_folder}').mkdir(parents=True, exist_ok=True)
        np.save(target_folder, heat_map)

    return polygons, heat_map, all_events_df, num_of_cells, distinct_cells, cells_xs, cells_ys, area_mean


