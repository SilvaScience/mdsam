# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 15:26:43 2025

@author: marfe
"""

from Data_process import DataImporter
from Data_process import DataProcess
from Data_plot import AxisCalib
from Data_plot import Plot


class Main_:
    def __init__(self,fileDirectory,dataSet,waitingTimeNum):
        
        self.fileDirectory = fileDirectory
        self.dataSet = dataSet
        self.waitingTimeNum =waitingTimeNum 
        
    def Analysis(self):
            importer = DataImporter(
                fileDirectory=self.fileDirectory,
                dataSet=self.dataSet,
                waitingTimeNum= self.waitingTimeNum
            )
            
            # Call the import_data method on the instance
            FTdata, dataTemp, nT1, numScans, nPixels,filePrefix = importer.import_data(numScans=0, start=0)
            
            process = DataProcess(
                FTdata, nT1, numScans, fftLength=0, fftAxis=0, apodizeData= True, bkgdCorrect=False
            )
            
            axisList, FF, avgTF, data,fftLength = process.data_process(FTdata)
            
            Calib = AxisCalib (nPixels, fftLength, numScans, calibrate = True, plotDiagonal = True, diagSlope = 1, diagIntercept = 233, specDisplay = True, cropUncalibratedPlot = True)
            
            probeIdx, pumpIdx, probeFreqs, pumpFreqs,diagonal = Calib.Calibration_data(diagMethod = int(0))
            
            Plotting = Plot (probeIdx, pumpIdx, filePrefix, self.fileDirectory, diagonal, numScans, self.waitingTimeNum, start = 0, calibrate = True, specDisplay = True,  symmetricContours = True, 
                         manualContourRange = False * [-0.4, 0.05], nContours = 20, manualAxisAspect = False, showProjections = True, plotDiagonal = True, flag_save_fig = True)
            
            return Plotting.Plot_data(probeFreqs, pumpFreqs, FF, diagMethod = int(0))
