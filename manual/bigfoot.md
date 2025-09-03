# Bigfoot Data Processing Example

This example demonstrates how to transform Bigfoot data files into HDF5 format and visualize the averaged result using the `mdsam` library.

## Files

All files need to be stored in a file structure as saved by default  by Bigfoot. 
Each scan is stored in a single file with the sample name extracted from the measurement metadata as well as date and run number as replacement for unique timestamp.

## Example Script

```python

import sys
from mdsam.BF import Bigfoot as bf
from mdsam.BF_plots import BFPlots as bfplot 


# %% Load data 
data_folder = r'\data\BIGFOOT\2025-06-13\Run031' 

h5_folder_and_filename = bf.transform_to_HDF5(data_folder)
data,header = bf.read_HDF5_file(h5_folder_and_filename)

# %%plot data
plot_range = [1545,1565] # plot range in meV
bfplot.plot_Reph_Re_Im_Abs(data,header,plot_range) # plot 2D rephasing spectra 
