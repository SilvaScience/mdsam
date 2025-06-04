# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 16:53:19 2025

@author: marfe
"""
# Import required modules from mdsam

import time
t_start = time.time()
import sys
import os
from os import listdir
import numpy as np # for array and matrix calculations
import matplotlib as mpl
import h5py
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PT_data_workup_functions import mlList
from PT_data_workup_functions import backgroundCorrect
from PT_data_workup_functions import hammingWindow
from PT_data_workup_functions import rowWiseLoop
from DATA_PROCESS import import_data
from DATA_PROCESS_2 import process_data
from PLOTTING import plotting_data

polyFitRange = list()
#from mdsam import PT_data_workup_functions

# DIRECTORY AND FILE SELECTION
fileDirectory = r"C:\Users\marfe\Downloads\example_data_processing_code_python_and_matlab" #Write the correct directory
dataSet = ["20250312", "33"] # WRITE the Prefix(NAME) and Experiment Number from QuickControl (usually 0)
waitingTimeNum =4 # SPECIFY the DELAY selected
numScans = 0 # number of scan files to average. If '0' all scans are averaged
SampleName = dataSet[0] 

# Save to HDF5 file
output_h5_file = os.path.join(fileDirectory, f"{SampleName}_saved.h5")

print(f"dataSet saved to: {output_h5_file}")

with h5py.File(output_h5_file, "w") as f:
    ascii_data = np.array(dataSet, dtype='S')  # Convert to byte strings for HDF5
    f.create_dataset("dataSet", data=ascii_data)
    
# Load the saved HDF5 data    
    input_h5_file = os.path.join(fileDirectory, f"{SampleName}_saved.h5")

with h5py.File(input_h5_file, "r") as f:
    raw_data = f["dataSet"][:]
    dataSet_loaded = [s.decode('utf-8') for s in raw_data]

print("Loaded dataSet from file:", dataSet_loaded)


# 2D PLOT
indexRangePump = [0, 500]  # SELECT range for the pump axis without Calibration
indexRangeProbe = [0, 32]   # SELECT range for the probe axis without Calibration
symmetricContours = True     # plots with contours and a colormap symmetric around 0
nContours = 20 # number of contours to plot
manualContourRange = False * [-0.4, 0.05]

# BACKGROUND CORRECTION
bkgdCorrect = False
polyFitRange += mlList(1, 10)
polyFitRange += mlList(100, 127)
polyFitOrder = 1

# DIAGONAL
plotDiagonal = True # plots a diagonal line (probe freq = pump freq)
diagSlope = 1
diagIntercept = 233


# CALIBRATION
calibrate = True # FREQUENCY CALIBRATION OF AXIS (Using Absorption spectra, read manual for more info.)
calibPixels = [25, 17]
calibFreqs =  [1965, 1980] # in wavemnumbers
freqRangePump = [1900, 2000] # SELECT range for the pump axis
freqRangeProbe = [1900,2000] # SELECT range for the probe axis
showProjections = True # adds projections along each frequency axis


# FIGURE SAVING
flag_save_fig = True # flag that turns on and off saving the figure
figure_savepath = os.path.join(fileDirectory, f"{SampleName}__{waitingTimeNum}.png") #Filepath



FTdata, dataTemp, nT1 =import_data(fileDirectory, SampleName, waitingTimeNum, dataSet, numScans, start=0)
axisList, FF, avgTF, data =process_data(FTdata, bkgdCorrect, polyFitRange, polyFitOrder, nT1, numScans)