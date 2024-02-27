from typing import List
import pandas as pd


class EuclideanMatrixSingleton:
    def __init__(self, cells_data: List[pd.DataFrame]):
        self.matrix = None
        self.set_distance_matrix(cells_data)

    def set_distance_matrix(self, cells_data: List[pd.DataFrame]):
        def Euclidean_Dist(df1):
            from scipy.spatial.distance import pdist, squareform
            return squareform(pdist(df1))

        mean_coordinates_df = pd.DataFrame()
        for cell in cells_data:
            cell_mean_x = cell['x'].mean()
            cell_mean_y = cell['y'].mean()
            cell_mean_coordinates = pd.DataFrame(data=[[cell_mean_x, cell_mean_y]], columns=['x', 'y'])
            mean_coordinates_df = pd.concat((mean_coordinates_df, cell_mean_coordinates))
        mean_coordinates_df.reset_index(inplace=True)
        mean_coordinates_df.drop(columns='index', inplace=True)
        euclidean_dist = Euclidean_Dist(mean_coordinates_df)
        self.matrix = euclidean_dist

    def get_distance_matrix(self):
        return self.matrix

    def get_above_diagonal(self):
        vals = []
        num_of_cells = self.matrix.shape[0]
        for i in range(num_of_cells):
            for j in range(i+1, num_of_cells):
                vals.append(self.matrix[i,j])
        return vals

    def get_topological_neighbors(self, cell_id, distance):
        neighbors = []
        potential_neighbors = list(self.matrix[cell_id])
        for ix, neighbor_dist in enumerate(potential_neighbors):
            if 0 < neighbor_dist <= distance:
                if ix < cell_id:
                    neighbors.append((ix, cell_id))
                else:
                    neighbors.append((cell_id, ix))
        return neighbors

