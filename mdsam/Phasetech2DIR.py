# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 13:02:15 2025

@author: marfe
"""

import os
import sys
import numpy as np
from os import listdir
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mdsam.PT_data_workup_functions import mlList
from mdsam.PT_data_workup_functions import backgroundCorrect
from mdsam.PT_data_workup_functions import hammingWindow
from mdsam.PT_data_workup_functions import rowWiseLoop

class ImportData:
    '''
        Module for 2DIR spectra from 2DQuickIR PhaseTech spectrometer
    '''
    def __init__(self, fileDirectory, dataSet, waitingTimeNum=0):
        self.fileDirectory = os.path.abspath(fileDirectory)
        self.dataSet = dataSet
        self.waitingTimeNum = waitingTimeNum

    def import_data(self, numScans, start):
        print("Directory:", self.fileDirectory)
        print("\tDetermining number of files to average.")

        fileList = listdir(self.fileDirectory)  # list of all files in the directory
        filePrefix = "#".join(self.dataSet) + (("_T" + "{:02}".format(self.waitingTimeNum)) if self.waitingTimeNum > 0 else "")
        fileExtension = ".scan"
        filteredList = [i for i in fileList if filePrefix in i and fileExtension in i]

        totalScans = len(filteredList)
        print("\t", totalScans, "files found")

        if totalScans == 0:
            print("No data files found matching prefix",
                  filePrefix, "in directory", self.fileDirectory, "\n")
            return
            # sys.exit("Stopping script execution.")
        
        output_h5_file = os.path.join(self.fileDirectory, f"{filePrefix}_saved.h5")

        print(f"dataSet saved to: {output_h5_file}")

        with h5py.File(output_h5_file, "w") as f:
            ascii_data = np.array(self.dataSet, dtype='S')  # Convert to ASCII string type
            f.create_dataset("self.dataSet", data=ascii_data)
        
        temp = np.loadtxt(os.path.join(self.fileDirectory, filteredList[0]))
        nT1, nPixels = temp.shape
        nT1 -= 1  # time points (ending row is not a time point)
        nPixels -= 1  # number of pixels (column 0 is time, not a pixel)

        A = numScans
        if A == 0:
            A = totalScans - start
        if A > totalScans:
            A = totalScans
        if A < 0:
            if start == 0:
                start = totalScans + A
            else:
                start = start + A
            if start < 0:
                start = 0
                A = -totalScans
            A = -A

        numScans = A
        scans2Avg = np.arange(numScans)
        FTdata = np.zeros([nT1, nPixels, numScans])

        for i in scans2Avg:
            dataTemp = np.loadtxt(os.path.join(self.fileDirectory, filteredList[start + i]))
            print("\tLoading file", i + 1, "of", numScans)
            FTdata[..., i] = dataTemp[:-1, 1:].copy()

        return FTdata, dataTemp, nT1, numScans, nPixels,filePrefix,output_h5_file
    
    
class DataProcess:
        def __init__(self, FTdata, nT1, numScans, fftLength=0, fftAxis=0, bkgdCorrect = False, apodizeData = True):
            self.FTdata = FTdata
            self.nT1 = nT1
            self.numScans = numScans
            self.fftLength = fftLength
            self.fftAxis = fftAxis
            self.bkgdCorrect = False
            self.apodizeData = True
      
    # Background correction
        def data_process(self, data):
            polyFitRange = list()
            polyFitRange += mlList(1, 10)
            polyFitRange += mlList(100, 127)
            polyFitOrder = 1
            self.data= self.FTdata
            if self.bkgdCorrect:
                print("#### Background Correction ####")
                self.data = rowWiseLoop(self.data, backgroundCorrect,self.polyFitRange)
    
    # Apodization
            if self.apodizeData:
                print("### Apodization ###")
                print("\tCalculating windowing function")
                H = hammingWindow(self.nT1, self.data)
                print("\tApodizing data")
                self.data = np.multiply(H,data)
    
    #  averaging
            print("### Averaging Scans ###")
            if len(self.data.shape) == 3:
                avgTF = np.mean(self.data[...,np.arange(self.numScans)], axis = 2)
            elif len(data.shape) == 2:
                print("\tData has only 2 dimensions. Not averaging data.")
            else:
                print("\tUnexpected shape of data:",data.shape)
                print("\tUnable to take average")
    
    # Fourier transform
    # calculate the FFT length, if requested
            if self.fftLength == 0:
                self.fftLength = 2**(np.ceil(np.log2(self.nT1))+1)
                axisList = ['column','row','page']
                avgTF[0,...] *= 0.5 # see Hamm and Zanni section 9.5.3
                print("### Fourier Transform ###")
                print("\tFFT Length:", self.fftLength)
                print("\tZero-padding by an additional factor of 2")
    # See Hamm & Zanni section 9.5.4 regarding FFT Length
            FF = np.real(np.fft.rfft(avgTF, n = int(2*self.fftLength), axis = self.fftAxis))
            FF = FF[:-1,...]
    
            return axisList, FF, avgTF, data, self.fftLength

class AxisCalib:
    def __init__(self,  nPixels, fftLength, numScans, calibrate = False, plotDiagonal = True, diagSlope = 1, diagIntercept = 233, specDisplay = True, cropUncalibratedPlot = True, start = 0):
        self.nPixels = nPixels
        self.fftLength = fftLength
        # self.filePrefix = filePrefix
        self.numScans = numScans
        self.calibrate = calibrate
        self.plotDiagonal = plotDiagonal
        self.diagSlope = diagSlope
        self.diagIntercept = diagIntercept
        self.specDisplay = specDisplay
        self.cropUncalibratedPlot = cropUncalibratedPlot
        self.start = start

    def Calibration_data(self, diagMethod = int(0)):
        self.diagMethod = diagMethod
# CALIBRATION
    # frequency calibration of axis
        calibPixels = [25, 17]
        calibFreqs =  [1965, 1980] # in wavemnumbers
        freqRangePump = [1900, 2000] # SELECT range for the pump axis
        freqRangeProbe = [1900,2000] # SELECT range for the probe axis
        showProjections = True # adds projections along each frequency axis

        diagMethod = int(0)  # don't change, only implemented method
        
        # 2D PLOT
        indexRangePump = [0, 500]  # SELECT range for the pump axis without Calibration
        indexRangeProbe = [0, 32]   # SELECT range for the probe axis without Calibration
        symmetricContours = True     # plots with contours and a colormap symmetric around 0
        nContours = 20 # number of contours to plot
        manualContourRange = False * [-0.4, 0.05]

# FREQUENCY AXES #
# calculate the probe and pump indices
        self.probeIdx = np.arange(self.nPixels)
        self.pumpIdx = np.arange(self.fftLength)
        probeFreqs = self.probeIdx # useful if not calibrating data
        pumpFreqs = self.pumpIdx # useful if not calibrating data

        if self.diagMethod == 0:
            self.diagonal = self.diagSlope*(self.nPixels-probeFreqs) + self.diagIntercept
            diagIdx = np.arange(len(self.diagonal))

        if self.calibrate:
            print("### Calibrating pump and probe axes ###")
            calibFreqs = np.array(calibFreqs)
            nm2cm1 = 1E7 # conversion between nm and cm-1
            calibParams = np.polyfit(calibPixels,nm2cm1/calibFreqs,len(calibPixels)-1)
            probeFreqs = nm2cm1/np.polyval(calibParams,probeFreqs)
            pumpFreqs = self.nPixels - (1/self.diagSlope)*(pumpFreqs-self.diagIntercept)
            pumpFreqs = nm2cm1/np.polyval(calibParams,pumpFreqs)
            diagonal = probeFreqs
    
            freqRangePump.sort()
            freqRangeProbe.sort()
    # Find the appropriate indices
            self.pumpIdx = ((pumpFreqs >= freqRangePump[0])*
                            (pumpFreqs<=freqRangePump[1]))
            self.probeIdx = ((probeFreqs >= freqRangeProbe[0])*
                             (probeFreqs<=freqRangeProbe[1]))
            diagIdx = ((diagonal >= freqRangePump[0])*
                       (diagonal <= freqRangePump[1]))
        elif not self.calibrate and self.cropUncalibratedPlot:
                self.pumpIdx = ((pumpFreqs >= indexRangePump[0])*
               (pumpFreqs<=indexRangePump[1]))
                self.probeIdx = ((probeFreqs >= indexRangeProbe[0])*
                (probeFreqs<= indexRangeProbe[1]))
                diagIdx = ((self.diagonal >= indexRangePump[0])*
               (self.diagonal <= indexRangePump[1]))
        else:
                    self.pumpIdx = np.arange(0,len(pumpFreqs))
                    self.probeIdx = np.arange(0,len(probeFreqs))
        return self.probeIdx, self.pumpIdx, probeFreqs, pumpFreqs,diagonal
    
class Plot:
    def __init__(self, probeIdx, pumpIdx,filePrefix, fileDirectory, diagonal, numScans, waitingTimeNum, calibrate=False, start = 0, specDisplay = True,  symmetricContours = True, 
                 manualContourRange = False * [-0.4, 0.05], nContours = 20, manualAxisAspect = False, showProjections = True, plotDiagonal = True, flag_save_fig = True): 
        self.probeIdx = probeIdx
        self.pumpIdx = pumpIdx
        self.filePrefix = filePrefix
        self.numScans = numScans
        self.fileDirectory = os.path.abspath(fileDirectory)
        self.waitingTimeNum = waitingTimeNum
        self.specDisplay = specDisplay
        self.calibrate = calibrate
        self.diagonal = diagonal
        self.symmetricContours = symmetricContours
        self.manualContourRange = manualContourRange
        self.nContours = nContours
        self.manualAxisAspect = manualAxisAspect
        self.showProjections = showProjections
        self.plotDiagonal = plotDiagonal
        self.flag_save_fig = flag_save_fig
        self.start = start
        
    def Plot_data(self,probeFreqs, pumpFreqs, FF, diagMethod = int(0)):
        self.diagMethod = diagMethod   
        
        polyFitRange = list()
        calibUnits = 'wn' # for displaying later

        # contour plot colormap options
        plotColorMap = 'bwr' # default color map for the contour plot
        symmetricContoursColorMap = 'bwr' # color map for the contour plot with symmetricContours = True
        colorbar = True # display a colorbar
        swapAxes = False # places prob freq on vertical axis
    
        if self.specDisplay:
            print('### Plotting ###')
            xStr = "probe"
            yStr = "pump"
            xInd = self.probeIdx
            yInd = self.pumpIdx
            x = probeFreqs[xInd].copy()
            y = pumpFreqs[yInd].copy()
            z = FF[yInd,:][:,xInd].copy()
        
        if not self.calibrate:
            xStr += " (pixels)"
            yStr += " (freq index)"
        else:
            xStr += " freq / " + calibUnits
            yStr += " freq / " + calibUnits
        
        if swapAxes:
            x, y, z  = y, x, z.T
            xStr, yStr = yStr, xStr
            self.diagonal, probeFreqs = probeFreqs, self.diagonal

        if not self.symmetricContours:
            if self.manualContourRange:
                lowZ = np.min(self.manualContourRange)
                highZ = np.max(self.manualContourRange)
            else:
                lowZ = np.min(FF)
                highZ = np.max(FF)
            mapString = plotColorMap
        else:
            if self.manualContourRange:
                highZ = np.min(np.abs(self.manualContourRange))
                if np.abs(self.manualContourRange[0]) != np.abs(self.manualContourRange[1]):
                    print("\tWarning: symmetric z contours requested, but manual "
                          "bounds supplied not symmetric.")
                    print("\tUsing the bound provided"
                          f" with the smallest absolute value ({highZ}) as max.")
            else:
                highZ = np.max(np.abs(FF))
            lowZ = -highZ
            mapString = symmetricContoursColorMap
        cints = np.linspace(lowZ-0.05*np.abs(lowZ),
                            highZ+0.05*np.abs(highZ),self.nContours)


        fig = plt.figure(figsize=[3.25,3.25],dpi= 300)
        
        titlestring = (self.filePrefix,'#',str(self.start+1).zfill(3),'-',
                        str(self.start+self.numScans).zfill(3),'.scan');
        titlestring = ''.join(titlestring)
        
        if self.manualAxisAspect:
            axisAspect = self.manualAxisAspect # alternate: 'auto'
        elif self.calibrate and not self.manualAxisAspect:
            axisAspect = 'equal'
        elif not self.calibrate and not self.manualAxisAspect:
            axisAspect = 'auto'

        mpl.rcParams['contour.negative_linestyle'] = 'solid'


        if not self.showProjections:
            contourPlotF = plt.contourf(x, y, z, cints, alpha = 0.75, 
                                        cmap = mapString)
            contourPlot = plt.contour(x, y, z, cints, colors = 'black', 
                                      linewidths=0.5)
            ax = fig.gca()
            if self.plotDiagonal:
                diag = plt.plot(probeFreqs, self.diagonal, color = 'black', 
                                linewidth = 0.5, linestyle = 'solid')
            ax.set_xlim(np.min(x), np.max(x))
            ax.set_ylim(np.min(y), np.max(y))
        else:
            fig.subplots_adjust(top = 0.9, left = 0.25)
            ax = fig.add_subplot(111)
        
            contourPlotF = plt.contourf(x, y, z, cints, alpha = 0.75, 
                                        cmap=mapString)
            contourPlot = plt.contour(x, y, z, cints, colors='black', 
                                      linewidths=0.5)
            if self.plotDiagonal:
                diag = plt.plot(probeFreqs, self.diagonal, color = 'black', 
                                linewidth = 0.25, linestyle = 'solid')
            ax.set_xlim(np.min(x), np.max(x))
            ax.set_ylim(np.min(y), np.max(y))
        
            
            sumX = np.sum(z, axis = 0)
            sumY = np.sum(z, axis = 1)
            sumAbsX = np.sum(np.abs(z), axis = 0)
            sumAbsY = np.sum(np.abs(z), axis = 1)
            if swapAxes:
                sumX = sumAbsX
            else:
                sumY = sumAbsY
            
            div = make_axes_locatable(ax)
            axSumX = div.append_axes("top", 0.5, pad = 0., sharex = ax)
            axSumY = div.append_axes("right", 0.5, pad = 0., sharey = ax) 
            axSumX.plot(x, sumX, linewidth = 0.75, color = "black", alpha = 0.75)
            axSumY.plot(sumY, y, linewidth = 0.75, color = "black", alpha=0.75)
            axSumY.tick_params(bottom = True, top = False, left = False, right = False,
                                   labelbottom = False, labelleft = False)
            axSumX.tick_params(bottom = False, top = False, left = True, right = False,
                                   labelbottom = False, labelleft = False)
            axSumY.spines['right'].set_visible(False)
            axSumY.spines['top'].set_visible(False)
            axSumX.spines['right'].set_visible(False)
            axSumX.spines['top'].set_visible(False)
            axSumX.set_title(titlestring, fontsize = 12)
        if colorbar:
            if not self.showProjections:
                div = make_axes_locatable(ax)
            axColorbar = div.append_axes("right", 0.05, pad = 0.05)
            axColorbar.tick_params(bottom = False, top = False, left = False, right = True,
                                   labelbottom = False, labelright = True, 
                                   labelleft = False)
            cbar = fig.colorbar(contourPlotF, cax = axColorbar)
            cbar.ax.tick_params(labelsize = 8)
            rounds = np.ceil(-np.log10(np.abs(cints)))
            newCints = np.zeros(rounds.shape)
            for i in enumerate(cints):
                r = round(i[1], int(rounds[i[0]]))
                newCints[i[0]] = r
    #        cbarInts = np.round(cints, 1)
            cbarTicks = [newCints[1], 0, newCints[-2]]
            cbar.set_ticks(cbarTicks)
        ax.set_aspect(axisAspect)
        ax.set_xlabel(xStr)
        ax.set_ylabel(yStr)
        if not self.calibrate and swapAxes:
            ax.invert_yaxis()
        elif not self.calibrate and not swapAxes:
            ax.invert_xaxis()
       


        self.flag_save_fig = True # flag that turns on and off saving the figure
        if self.flag_save_fig:
            figure_savepath = os.path.join(self.fileDirectory, f"{self.filePrefix}__{self.waitingTimeNum}.png") #Filepath
            plt.savefig(figure_savepath)
        plt.show()
        return plt.show()
    
class PhaseTech2DIRData:
    def __init__(self,fileDirectory,dataSet):
        
        self.fileDirectory = fileDirectory
        self.dataSet = dataSet
        
    def Analysis(self):
            #with open("DataAnalysis.txt", "w") as f:
             #   f.write (f"{nPixels} {fftLength} {numScans} {filePrefix} {FF} \n")
            file_path = "Data_analysis.txt"
            with open(file_path,"w") as f:
                for a in range(11):
                # Call the import_data method on the instance
                  
                        self.waitingTimeNum = a
                        try:
                            importer = ImportData(
                                fileDirectory=self.fileDirectory,
                                dataSet=self.dataSet,
                                waitingTimeNum = self.waitingTimeNum
                                )
                            FTdata, dataTemp, nT1, numScans, nPixels,filePrefix, output_h5_file = importer.import_data(numScans=0, start=0)
                            
                            process = DataProcess(FTdata, nT1, numScans, fftLength=0, fftAxis=0, apodizeData= True, bkgdCorrect=False)
                    
                            axisList, FF, avgTF, data,fftLength = process.data_process(FTdata) 
                            # nPixels = np.array(nPixels)
                            # fftLength = np.array(fftLength)
                            # numScans = np.array(numScans)
                            # filePrefix = np.array(filePrefix)
                            # FF = np.array(FF)
                            f.write(f"{nPixels}\n")
                            f.write(f"{fftLength}\n")
                            f.write(f"{numScans}\n")
                            f.write(f"{filePrefix}\n")
                         #   f.write(f"{FF}\n")
                        except Exception as e:
                                               print(f'Maximum number of data files is: {e}')


            with open("array3D.txt","w") as f:
                 for a in range(11):
                     self.waitingTimeNum = a
                     try:
                           importer = ImportData(
                               fileDirectory=self.fileDirectory,
                               dataSet=self.dataSet,
                               waitingTimeNum = self.waitingTimeNum
                               )
                           
                           FTdata, dataTemp, nT1, numScans, nPixels,filePrefix, output_h5_file = importer.import_data(numScans=0, start=0)
                           for i, slice_2d in enumerate(FTdata):
                                       f.write(f"Slice {i}:\n")
                                       np.savetxt(f, slice_2d, fmt="%.6f", delimiter="\t")
                                       f.write("\n") 
                                     
                     except Exception as e:
                                      print()                       
            with open("FF.txt","w") as f:
                  for a in range(11):
                      self.waitingTimeNum = a
                      try:
                            importer = ImportData(
                                fileDirectory=self.fileDirectory,
                                dataSet=self.dataSet,
                                waitingTimeNum = self.waitingTimeNum
                                )
                            FTdata, dataTemp, nT1,   numScans, nPixels,filePrefix, output_h5_file = importer.import_data(numScans=0, start=0)
                            
                            process = DataProcess(FTdata, nT1, numScans, fftLength=0, fftAxis=0, apodizeData= True, bkgdCorrect=False)
                            axisList, FF, avgTF, data,fftLength = process.data_process(FTdata) 

                           
                            np.savetxt(f, FF, delimiter="\t")
                                           
                      except Exception as e:
                                        print()

    def Graphs(self):
                   
                    nPixels_ = []
                    fftLength_ = []
                    numScans_ = []
                    filePrefix_ = []
                    with open("Data_analysis.txt", "r", encoding="utf-8") as f:
                        lines = f.readlines()
                    k=0
                    for line in lines:
                        parts = line.strip().split()  # split by spaces
                        if k == 0:
                            nPixels_.append(int(parts[0]))
                        if k == 1:
                            fftLength_.append(float(parts[0]))
                        if k == 2:
                            numScans_.append(int(parts[0]))
                        if k == 3:
                            filePrefix_.append(parts[0])
                        k= k +1
                        if k ==4:
                            k=0
                    with open("FF.txt", "r", encoding="utf-8") as f:
                                lines = f.readlines()
                               
                    k = 0  
                    FF_ = []
                    FF_d = []
                    for line in lines:
                                    if k  == 0:
                                        FF_d = []
                                    parts = line.strip().split()
                                    parts_ = list(map(float, parts))
                                    FF_d.append(parts_)
                                    k= k + 1
                                    if k==512:
                                        FF_.append(np.array(FF_d))
                                        k = 0
 
                 
        
                    for a in range(len(FF_)):
                        if a!=0:
                            nPixels = nPixels_[a]
                            fftLength = fftLength_[a]
                            numScans = numScans_[a]
                            filePrefix = filePrefix_[a]
                            FF = FF_[a]
                            
                            Calib = AxisCalib (nPixels, fftLength, numScans, calibrate = True, plotDiagonal = True, diagSlope = 1, diagIntercept = 233, specDisplay = True, cropUncalibratedPlot = True)
                                  
                            probeIdx, pumpIdx, probeFreqs, pumpFreqs,diagonal = Calib.Calibration_data(diagMethod = int(0))
                                  
                            Plotting = Plot (probeIdx, pumpIdx, filePrefix, self.fileDirectory, diagonal, numScans, self.waitingTimeNum, start = 0, calibrate = True, specDisplay = True,  symmetricContours = True, 
                                                 manualContourRange = False * [-0.4, 0.05], nContours = 20, manualAxisAspect = False, showProjections = True, plotDiagonal = True, flag_save_fig = True)
                                    
                            Plotting.Plot_data(probeFreqs, pumpFreqs, FF, diagMethod = int(0))     
                                       


                    
                   
                
