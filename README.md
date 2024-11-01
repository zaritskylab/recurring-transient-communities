
# main module
[main.py](main.py) module is the main code that runs all analyses.
It generates aggregated data based on the preprocessed and shuffled data.
### generate_cross_experiment_spatial_significance_data 
This function generates the data for the spatial significance analysis, as in Figures 1-3.
### generate_per_sample_community_statistics
This function generates the data for the community statistics analysis, as in Figures 1-3.

# Other executables
The following modules are preprocessing steps that are not executed by the main code.
## arcos_shuffle_analysis
[arcos_shuffle_analysis.py](data_layer%2Farcos%2Farcos_shuffle_analysis.py) module is responsible for creating the data for the ARCOS shuffle analysis, as in Figure1.
The main code doesn't run it directly, as this function is randomized and may result in different outputs from one execution to another.

## arcos_hotspots_shuffle_analysis
[arcos_hotspots_shuffle_analysis.py](data_layer%2Farcos_hotspots_shuffle_analysis.py) module is responsible for creating the data for the ARCOS hotspots shuffle analysis, as in Figure2.
The main code doesn't run it directly, as this function is randomized and may result in different outputs from one execution to another.

## preprocess_data
[preprocess_data.py](data_layer%2Fpreprocess_data.py) module is responsible for preprocessing the raw data received from the imaging.
The main code doesn't run it directly, as we are using the preprocessed data.
Tha module is included in the repository for reference.
