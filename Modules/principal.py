# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 10:25:30 2025

@author: marfe
"""
"From fileName import NameClass"
from Data_process import DataImporter
from Data_process import DataProcess
from Data_plot import AxisCalib
from Data_plot import Plot

# Create an instance of DataImporter
fileDirectory=r"C:\Users\marfe\Downloads\example_data_processing_code_python_and_matlab"
waitingTimeNum= 2

importer = DataImporter(
    fileDirectory=r"C:\Users\marfe\Downloads\example_data_processing_code_python_and_matlab",
    dataSet=["20250312", "33"],
    waitingTimeNum= 2
)

# Call the import_data method on the instance
FTdata, dataTemp, nT1, numScans, nPixels,filePrefix, output_h5_file = importer.import_data(numScans=0, start=0)

process = DataProcess(
    FTdata, nT1, numScans, fftLength=0, fftAxis=0, apodizeData= True, bkgdCorrect=False
)

axisList, FF, avgTF, data,fftLength = process.data_process(FTdata)

Calib = AxisCalib (nPixels, fftLength, numScans, calibrate = True, plotDiagonal = True, diagSlope = 1, diagIntercept = 233, specDisplay = True, cropUncalibratedPlot = True)

probeIdx, pumpIdx, probeFreqs, pumpFreqs,diagonal = Calib.Calibration_data(diagMethod = int(0))

Plotting = Plot (probeIdx, pumpIdx, filePrefix, fileDirectory, diagonal, numScans, waitingTimeNum, start = 0, calibrate = True, specDisplay = True,  symmetricContours = True, 
             manualContourRange = False * [-0.4, 0.05], nContours = 20, manualAxisAspect = False, showProjections = True, plotDiagonal = True, flag_save_fig = True)

Plotting.Plot_data(probeFreqs, pumpFreqs, FF, diagMethod = int(0))