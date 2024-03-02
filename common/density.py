import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

def _generate_rectangle_polys(p: Point, r=7):
    x, y = p.x, p.y

    p1 = (x+r, y+r)
    p2 = (x-r, y+r)
    p3 = (x-r, y-r)
    p4 = (x+r, y-r)
    return Polygon([p1,p2,p3,p4])


def get_cell_density(arcos_df: pd.DataFrame):
    coordinates_df = arcos_df[['cell_id', 'x_microns', 'y_microns']].groupby('cell_id').mean()
    coordinates = list(coordinates_df.T.to_dict().values())
    points = [Point((c['x_microns'], c['y_microns'])) for c in coordinates]
    point_areas = [_generate_rectangle_polys(p) for p in points]

    densities = []
    for i in range(len(points)):
        neighbors = 0
        for j in range(len(points)):
            if i == j:
                continue

            if points[j].within(point_areas[i]):
                neighbors += 1

        densities.append(neighbors)

    return densities


def get_activations_per_cell(arcos_df: pd.DataFrame):
    activations = []

    for cell_id, cell_df in arcos_df.groupby('cell_id'):
        activation_series = cell_df['normalized_intensity.bin'].tolist()
        activation_series = ''.join([str(i) for i in activation_series])

        n_activations = activation_series.count('10')
        if activation_series[-1] == '1':
            n_activations += 1

        activations.append(n_activations)

    return activations
