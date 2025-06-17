# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 12:16:19 2025

@author: marfe
"""
import os
import sys
import numpy as np
from os import listdir
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PT_data_workup_functions import mlList
from PT_data_workup_functions import backgroundCorrect
from PT_data_workup_functions import hammingWindow
from PT_data_workup_functions import rowWiseLoop

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