HT_Chemotaxis_Toolbox

This repository contains a set of MATLAB functions that can be used for analysis of high-throughput live-cell chemotaxis experiments.
It relies on the Matlab's Image Processing Toolbox.

Usage:
The chemotaxis_load_new_data_script provides a template script for loading and processing data.
This script will create HTChemotaxisExperiment objects which will contain processed data, and has associated functions for processing and viewing data.
The script also contains comments to provide guidance on how to load and process data.

Setting processing parameters:
Default parameters for analyzing images and cell movement data are contained in the .m files in the default_parameters folder.
These parameters should be checked and adjusted as needed for a given microscope system.

Data organization - folder structure, filenames, etc.
The repository is designed with a particular folder structure for image data in mind, as well as naming conventions for image files, and image file formats.
It expects .tif files for individual images.
It expects the name of the imaging channel (e.g. RFP) to be present in the filenames.
It expects images to be organized with a subfolder for each well inside a master folder for each 96-well plate experiment.
The functions can be adjusted for other image formats and file conventions.
In particular, the functions in the "read" folder are the ones that would need to be adjusted.