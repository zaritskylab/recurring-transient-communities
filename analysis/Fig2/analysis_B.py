from typing import Dict, List

import pandas as pd
from libpysal.cg.shapes import Point
import matplotlib
import numpy as np
from matplotlib import pyplot as plt, cm
from matplotlib.patches import Polygon
from scipy.spatial import Voronoi
from common.constants import FIGURES_LOCATION, ARCOS_OUTPUT_LOCATION
from common.density import get_cell_density, get_activations_per_cell


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

def plot_voronoi(vor: Voronoi, points: Dict, region_colors: List[int], points_colors: List[int], cbar_label: str):

    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot(111)

    points_colors_minima = min(points_colors)
    points_colors_maxima = max(points_colors)

    # normalize chosen colormap
    max_val = max(region_colors)
    region_colors_norm = matplotlib.colors.Normalize(vmin=0, vmax=max_val, clip=True)
    region_colors_mapper = cm.ScalarMappable(norm=region_colors_norm, cmap=plt.cm.get_cmap('Reds', max_val))
    region_colors_mapper.set_array([])
    cbar = plt.colorbar(region_colors_mapper, ticks=np.linspace(0, max_val, max_val + 1))
    cbar.set_label(cbar_label, size=20)
    cbar.ax.tick_params(labelsize=18)  # set your label size here

    points_colors_norm = matplotlib.colors.Normalize(vmin=points_colors_minima, vmax=points_colors_maxima, clip=True)
    points_colors_mapper = cm.ScalarMappable(norm=points_colors_norm, cmap=cm.Greys)

    # plot Voronoi diagram, and fill finite regions with color mapped from speed value
    regions, vertices = voronoi_finite_polygons_2d(vor)
    for i in range(len(regions)):
        region = regions[i]
        polygon = vertices[region]
        colored_cell = Polygon(polygon,
                               linewidth=0.7,
                               facecolor=region_colors_mapper.to_rgba(region_colors[i]),
                               edgecolor='black')
        ax.add_patch(colored_cell)
        ax.annotate(i, (points[i][0], points[i][1] * 0.995), fontsize=9, weight='bold',
                    color='black')

    colors = [points_colors_mapper.to_rgba(points_colors[i]) for i in range(len(points))]
    plt.scatter([p[0] for p in points.values()], [p[1] for p in points.values()], s=10, c=colors)

    min_x = min([p[0] for p in points.values()])
    min_y = min([p[1] for p in points.values()])
    max_x = max([p[0] for p in points.values()])
    max_y = max([p[1] for p in points.values()])
    plt.ylim(min_y - 10, max_y + 10)
    plt.xlim(min_x - 10, max_x + 10)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_LOCATION}Fig_2_B.png')
    plt.close(fig=fig)


def calculate_activations_vornoi(arcos_df):
    points = arcos_df[['cell_id', 'x_microns', 'y_microns']].groupby('cell_id').mean().T.to_dict()
    points = {i: (point_dict['x_microns'], -point_dict['y_microns']) for i, point_dict in points.items()}

    vor = Voronoi(list(points.values()))
    regions, vertices = voronoi_finite_polygons_2d(vor)

    polygons_by_region_id = {}
    cell_id_to_region_id = {}

    for region_id, region in enumerate(regions):
        polygon = Polygon([Point(p) for p in list(vertices[region])])
        polygons_by_region_id[region_id] = polygon

    for cell_id, coords in points.items():
        for region_id, polygon in polygons_by_region_id.items():
            if polygon.contains_point(Point(coords)):
                cell_id_to_region_id[cell_id] = region_id

    # find min/max values for normalization
    activations = get_activations_per_cell(arcos_df)
    densities = get_cell_density(arcos_df)

    plot_voronoi(vor, points, activations, densities, cbar_label='#Activations')

def main():
    relevant_experiment_file = f'{ARCOS_OUTPUT_LOCATION}/mid_third_wild_type/sample_11.csv'
    arcos_df = pd.read_csv(relevant_experiment_file)
    calculate_activations_vornoi(arcos_df)
