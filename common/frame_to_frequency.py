import math
import pandas as pd

class FrequencyTranslator():
    def __init__(self):
        self.frequencies_df = pd.read_csv(f'./common/frame_to_seconds.csv')

    def get_movie_in_minutes(self, experiment_type: str, experiment_name: str):
        """ Returns the movie length in minutes """
        try:
            frequency_record = self.frequencies_df[(self.frequencies_df['experiment_name'] == experiment_name) &
                                       (self.frequencies_df['experiment_type'] == experiment_type)].to_dict(orient='records')
            movie_in_minutes = frequency_record[0]['frames_number'] * frequency_record[0]['frame_to_seconds'] / 60
            return movie_in_minutes
        except Exception as e:
            return None

    def temporal_lag_to_frames(self, experiment_type: str, experiment_name: str, temporal_lag: int):
        """ Returns the movie length in minutes """
        import math
        try:
            frequency_record = self.frequencies_df[(self.frequencies_df['experiment_name'] == experiment_name) &
                                       (self.frequencies_df['experiment_type'] == experiment_type)].to_dict(orient='records')
            frames = int(math.ceil(temporal_lag / frequency_record[0]['frame_to_seconds']))
            return frames
        except Exception as e:
            return None

    def frame_lag_to_seconds(self, experiment_type: str, experiment_name: str, frame_lag: int):
        """ Returns the movie length in minutes """
        try:
            frequency_record = self.frequencies_df[(self.frequencies_df['experiment_name'] == experiment_name) &
                                       (self.frequencies_df['experiment_type'] == experiment_type)].to_dict(orient='records')
            seconds = int(math.floor(frame_lag * frequency_record[0]['frame_to_seconds']))
            return seconds
        except Exception as e:
            return None

    def get_fps(self, experiment_type: str, experiment_name: str):
        """ Returns the movie length in minutes """
        try:
            frequency_record = self.frequencies_df[(self.frequencies_df['experiment_name'] == experiment_name) &
                                       (self.frequencies_df['experiment_type'] == experiment_type)].to_dict(orient='records')
            fps = frequency_record[0]['frame_to_seconds']
            return fps
        except Exception as e:
            return None

    def get_max_fps(self):
        """ Returns the movie length in minutes """
        try:
            max_fps = self.frequencies_df['frame_to_seconds'].max()
            return max_fps
        except Exception as e:
            return None

    def get_all_lengths(self) -> pd.DataFrame:
        """ Returns the movies length in minutes """
        self.frequencies_df['movie_length'] = self.frequencies_df['frames_number'] * self.frequencies_df['frame_to_seconds'] / 60
        return self.frequencies_df[['experiment_type', 'movie_length']]

    def get_names_df(self) -> pd.DataFrame:
        """ Returns the movies length in minutes """
        return self.frequencies_df[['experiment_type', 'experiment_name']]

    def get_records_below_fps(self, fps: float) -> pd.DataFrame:
        """ Returns the movies length in minutes """
        return self.frequencies_df[self.frequencies_df['frame_to_seconds'] < (1/fps)]