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
from pathlib import Path
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

        fileList = listdir(self.fileDirectory)  # list of all files in the directory
        filePrefix = "#".join(self.dataSet) + (("_T" + "{:02}".format(self.waitingTimeNum)) if self.waitingTimeNum > 0 else "")
        fileExtension = ".scan"
        filteredList = [i for i in fileList if filePrefix in i and fileExtension in i]

        totalScans = len(filteredList)

        if totalScans == 0:
            raise FileNotFoundError("No data files found matching prefix",
                  filePrefix, "in directory", self.fileDirectory, "\n")
        
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
            FTdata[..., i] = dataTemp[:-1, 1:].copy()

        return FTdata, dataTemp, nT1, numScans, nPixels,filePrefix,fileList

class DataProcess:
        def __init__(self, FTdata, nT1, numScans, fftLength=0, fftAxis=0, bkgdCorrect = False, apodizeData = True):
            self.data = FTdata
            self.nT1 = nT1
            self.numScans = numScans
            self.fftLength = fftLength
            self.fftAxis = fftAxis
            self.bkgdCorrect = False
            self.apodizeData = True
    # Background correction
        def data_process(self):
            polyFitRange = list()
            polyFitRange += mlList(1, 10)
            polyFitRange += mlList(100, 127)
            polyFitOrder = 1
            if self.bkgdCorrect:
                self.data = rowWiseLoop(self.data, backgroundCorrect,self.polyFitRange)
    # Apodization
            if self.apodizeData:
                ### Apodization ###
                #Calculating windowing function
                H = hammingWindow(self.nT1, self.data)
                #Apodizing data
                self.data = np.multiply(H,self.data)
    #  averaging
            ### Averaging Scans if more than 2 dimensions###
            if len(self.data.shape) == 3:
                avgTF = np.mean(self.data[...,np.arange(self.numScans)], axis = 2)
            else:
                raise UserWarning("Unexpected shape of data %s unable to take average"%self.data.shape)
    
    # Fourier transform
    # calculate the FFT length, if requested
            if self.fftLength == 0:
                self.fftLength = 2**(np.ceil(np.log2(self.nT1))+1)
                axisList = ['column','row','page']
                avgTF[0,...] *= 0.5 # see Hamm and Zanni section 9.5.3
                ### Fourier Transform ###
                #Zero-padding by an additional factor of 2
    # See Hamm & Zanni section 9.5.4 regarding FFT Length
            FF = np.real(np.fft.rfft(avgTF, n = int(2*self.fftLength), axis = self.fftAxis))
            FF = FF[:-1,...]
            return axisList, FF, avgTF, self.data, self.fftLength

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
            ### Calibrating pump and probe axes ###
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
        self.pumpFreqs=pumpFreqs
        self.probeFreqs=probeFreqs
        self.diagonal=diagonal
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
        """
            Plots the data held in the 
        """
        self.diagMethod = diagMethod
        polyFitRange = list()
        calibUnits = 'wn' # for displaying later

        # contour plot colormap options
        plotColorMap = 'bwr' # default color map for the contour plot
        symmetricContoursColorMap = 'bwr' # color map for the contour plot with symmetricContours = True
        colorbar = True # display a colorbar
        swapAxes = False # places prob freq on vertical axis
    
        if self.specDisplay:
            ### Plotting ###
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
                    raise UserWarning("\tWarning: symmetric z contours requested, but manual bounds supplied not symmetric. Using the bound provided with the smallest absolute value ({highZ}) as max.")
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

        if self.flag_save_fig:
            figure_savepath = os.path.join(self.fileDirectory, f"{self.filePrefix}__{self.waitingTimeNum}.png") #Filepath
            plt.savefig(figure_savepath)
        plt.show()
        return fig,ax
    
class PhaseTech2DIRData:
    def __init__(self,fileDirectory,dataSet=None):
        '''
         Loads the dataset contained in the specified file directory and analyzes it.
         input:
            - fileDirectory (str or Path): Path to the directory holding the dataset
            - dataSet (list of 2 str): A list holding two strings specifying the date tag (e.g. 20250312) and the experiment number (e.g 33).
                These assume the filenames are labelled according to YYYYMMDD#XX_TYY#ZZZ.scan where YYYYMMDD is the date and XX the experiment number
        ''' 
        self.fileDirectory = fileDirectory
        if os.path.isfile(fileDirectory):
            self.import_from_h5(fileDirectory)
        else:
            self.dataSet = dataSet
            self.Analysis()
        
    def Analysis(self):
        """
            Analyzes the 2DIR data held at self.fileDirectory and self.dataSet provided during initialization
        """
        for a in range(11):
        # Call the import_data method on the instance
            self.waitingTimeNum = a
            try:
                importer = ImportData(
                    fileDirectory=self.fileDirectory,
                    dataSet=self.dataSet,
                    waitingTimeNum = self.waitingTimeNum
                    )
                self.FTdata, self.dataTemp, self.nT1, self.numScans, self.nPixels,self.filePrefix,self.fileList= importer.import_data(numScans=0, start=0)
                self.process = DataProcess(self.FTdata, self.nT1, self.numScans, fftLength=0, fftAxis=0, apodizeData= True, bkgdCorrect=False)
                self.axisList,self.FF,self.avgTF,self.data,self.fftLength = self.process.data_process()
            except Exception as e:
                print(f'Maximum number of data files is: {e}')

    def Graphs(self):
        self.calib = AxisCalib (self.nPixels, self.fftLength, self.numScans, calibrate = True, plotDiagonal = True, diagSlope = 1, diagIntercept = 233, specDisplay = True, cropUncalibratedPlot = True)
        self.calib.Calibration_data(diagMethod = int(0))
        Plotting = Plot (self.calib.probeIdx, self.calib.pumpIdx, self.filePrefix, self.fileDirectory, self.calib.diagonal, self.numScans, self.waitingTimeNum, start = 0, calibrate = True, specDisplay = True,  symmetricContours = True, 
                                manualContourRange = False * [-0.4, 0.05], nContours = 20, manualAxisAspect = False, showProjections = True, plotDiagonal = True, flag_save_fig = False)
        Plotting.Plot_data(self.calib.probeFreqs, self.calib.pumpFreqs, self.FF, diagMethod = int(0))

    def export_to_h5(self, h5_filename):
        """
            Exports the data contained in the object to a H5 file
            input:
            - h5_filename: relative or absolute path where to save the H5 file
        """
        # check if provided data_path is absolute
        p = Path(h5_filename)
        if not p.is_absolute():
            cwd = Path(os.getcwd())
            output_path= cwd / p
        else:
            output_path=h5_filename
        print(output_path)
        data={}
        data['FTdata']=self.FTdata
        data['FF']=self.FF
        data['avgTF']=self.avgTF
        data['data']=self.data
        data['pumpFreqs']=self.calib.pumpFreqs
        data['probeFreqs']=self.calib.probeFreqs
        header={}
        header['filePrefix']=self.filePrefix
        header['fileList']=self.fileList
        header['dataSet']=self.dataSet
        header['nT1']=self.nT1
        header['numScans']=self.numScans
        header['nPixels']=self.nPixels
        header['fftLength']=self.fftLength
        header['waitingTimeNum']=self.waitingTimeNum
        #### save to HDF5 ####
        with h5py.File(output_path, 'w') as hdf:
            # Save all numerical arrays
            for key, value in data.items():
                if isinstance(value, np.ndarray):
                    hdf.create_dataset(f"data/{key}", data=value)
                elif isinstance(value, dict):
                    grp = hdf.create_group(f"data/{key}")
                    for subkey, subval in value.items():
                        if isinstance(subval, np.ndarray):
                            grp.create_dataset(subkey, data=subval)
                        elif isinstance(subval, dict):
                            subgrp = hdf.create_group(f"data/{key}/{subkey}")
                            for subsubkey, subsubval in subval.items():
                                subgrp.create_dataset(subsubkey, data=subsubval)
            # Save header information
            header_grp = hdf.create_group("header")
            for key, val in header.items():
                if isinstance(val, (str, float, int)):
                    header_grp.attrs[key] = val
                elif isinstance(val, (list, np.ndarray)):
                    arr = np.array(val)
                    if arr.dtype.kind in {'U', 'S'}:  # it's a string or unicode string array
                        header_grp.attrs[key] = [str(x) for x in arr]
                    else:
                        header_grp.create_dataset(key, data=arr)
                elif isinstance(val, dict):
                    subgrp = header_grp.create_group(key)
                    for sk, sv in val.items():
                        try:
                            subgrp.attrs[sk] = sv
                        except TypeError:
                            pass  # Skip complex/unserializable types

    def import_from_h5(self,file_path):
        """Read data and header from a 2DIR HDF5 file.
        
        Args:
            file_path (str): Path to the .h5 file
        
        Returns:
            data (dict): Nested dictionary of numerical datasets
            header (dict): Metadata dictionary
        """
        data = {}
        header = {}
        
        with h5py.File(file_path, 'r') as hdf:
            # Load data. Check at every level of nested dict whether item is data or subgroup
            data_group = hdf['data']
            for key in data_group:
                item = data_group[key]
                if isinstance(item, h5py.Dataset):
                    data[key] = item[()]
                elif isinstance(item, h5py.Group):
                    data[key] = {}
                    for subkey in item:
                        if isinstance(item[subkey], h5py.Group):
                            data[key][subkey] = {}
                            for subsubkey in item[subkey].keys():
                                data[key][subkey][subsubkey] = item[subkey][subsubkey][()]
                        else:
                            data[key][subkey] = item[subkey][()]
        
            # Load header
            header_group = hdf['header']
            for key in header_group.attrs:
                header[key] = header_group.attrs[key]
        
            for key in header_group:
                item = header_group[key]
                if isinstance(item, h5py.Dataset):
                    header[key] = item[()]
                elif isinstance(item, h5py.Group):
                    header[key] = {}
                    for subkey in item.attrs:
                        header[key][subkey] = item.attrs[subkey]
                    for subkey in item:
                        header[key][subkey] = item[subkey][()]
        # Now import the dicts inside attributes of the class
        for key in data:
            setattr(self, key, data[key])
        for key in header:
            setattr(self, key, header[key])


