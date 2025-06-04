# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 08:17:26 2025

@author: marfe
"""
import time
t_start = time.time()
import sys
import os
from os import listdir
#import json
#import pickle
import numpy as np # for array and matrix calculations
import matplotlib as mpl
import h5py
#mpl.use("pdf")
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PT_data_workup_functions import mlList
from PT_data_workup_functions import backgroundCorrect
from PT_data_workup_functions import hammingWindow
from PT_data_workup_functions import rowWiseLoop
polyFitRange = list()
diagMethod = int(0)  # don't change, only implemented method


def plotting_data()
# FREQUENCY AXES #
# calculate the probe and pump indices
probeIdx = np.arange(nPixels)
pumpIdx = np.arange(fftLength)
probeFreqs = probeIdx # useful if not calibrating data
pumpFreqs = pumpIdx # useful if not calibrating data

if diagMethod == 0:
    diagonal = diagSlope*(nPixels-probeFreqs) + diagIntercept
    diagIdx = np.arange(len(diagonal))

if calibrate:
    print("### Calibrating pump and probe axes ###")
    calibFreqs = np.array(calibFreqs)
    nm2cm1 = 1E7 # conversion between nm and cm-1
    calibParams = np.polyfit(calibPixels,nm2cm1/calibFreqs,len(calibPixels)-1)
    probeFreqs = nm2cm1/np.polyval(calibParams,probeFreqs)
    pumpFreqs = nPixels - (1/diagSlope)*(pumpFreqs-diagIntercept)
    pumpFreqs = nm2cm1/np.polyval(calibParams,pumpFreqs)
    diagonal = probeFreqs
    
    freqRangePump.sort()
    freqRangeProbe.sort()
    # Find the appropriate indices
    pumpIdx = ((pumpFreqs >= freqRangePump[0])*
               (pumpFreqs<=freqRangePump[1]))
    probeIdx = ((probeFreqs >= freqRangeProbe[0])*
                (probeFreqs<=freqRangeProbe[1]))
    diagIdx = ((diagonal >= freqRangePump[0])*
               (diagonal <= freqRangePump[1]))
elif not calibrate and cropUncalibratedPlot:
    pumpIdx = ((pumpFreqs >= indexRangePump[0])*
               (pumpFreqs<=indexRangePump[1]))
    probeIdx = ((probeFreqs >= indexRangeProbe[0])*
                (probeFreqs<= indexRangeProbe[1]))
    diagIdx = ((diagonal >= indexRangePump[0])*
               (diagonal <= indexRangePump[1]))
else:
    pumpIdx = np.arange(0,len(pumpFreqs))
    probeIdx = np.arange(0,len(probeFreqs))
    
    if specDisplay:
    print('### Plotting ###')
    xStr = "probe"
    yStr = "pump"
    xInd = probeIdx
    yInd = pumpIdx
    x = probeFreqs[xInd].copy()
    y = pumpFreqs[yInd].copy()
    z = FF[yInd,:][:,xInd].copy()
    
    if not calibrate:
        xStr += " (pixels)"
        yStr += " (freq index)"
    else:
        xStr += " freq / " + calibUnits
        yStr += " freq / " + calibUnits
    
    if swapAxes:
        x, y, z  = y, x, z.T
        xStr, yStr = yStr, xStr
        diagonal, probeFreqs = probeFreqs, diagonal

    if not symmetricContours:
        if manualContourRange:
            lowZ = np.min(manualContourRange)
            highZ = np.max(manualContourRange)
        else:
            lowZ = np.min(FF)
            highZ = np.max(FF)
        mapString = plotColorMap
    else:
        if manualContourRange:
            highZ = np.min(np.abs(manualContourRange))
            if np.abs(manualContourRange[0]) != np.abs(manualContourRange[1]):
                print("\tWarning: symmetric z contours requested, but manual "
                      "bounds supplied not symmetric.")
                print("\tUsing the bound provided"
                      f" with the smallest absolute value ({highZ}) as max.")
        else:
            highZ = np.max(np.abs(FF))
        lowZ = -highZ
        mapString = symmetricContoursColorMap
    cints = np.linspace(lowZ-0.05*np.abs(lowZ),
                        highZ+0.05*np.abs(highZ),nContours)


    fig = plt.figure(figsize=[3.25,3.25],dpi= 300)
    
    titlestring = (filePrefix,'#',str(start+1).zfill(3),'-',
                   str(start+numScans).zfill(3),'.scan');
    titlestring = ''.join(titlestring)
    
    if manualAxisAspect:
        axisAspect = manualAxisAspect # alternate: 'auto'
    elif calibrate and not manualAxisAspect:
        axisAspect = 'equal'
    elif not calibrate and not manualAxisAspect:
        axisAspect = 'auto'

    mpl.rcParams['contour.negative_linestyle'] = 'solid'


    if not showProjections:
        contourPlotF = plt.contourf(x, y, z, cints, alpha = 0.75, 
                                    cmap = mapString)
        contourPlot = plt.contour(x, y, z, cints, colors = 'black', 
                                  linewidths=0.5)
        ax = fig.gca()
        if plotDiagonal:
            diag = plt.plot(probeFreqs, diagonal, color = 'black', 
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
        if plotDiagonal:
            diag = plt.plot(probeFreqs, diagonal, color = 'black', 
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
        if not showProjections:
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
    if not calibrate and swapAxes:
        ax.invert_yaxis()
    elif not calibrate and not swapAxes:
        ax.invert_xaxis()
   


if flag_save_fig:
 #fn = os.path.join(os.path.abspath(figure_savepath),figure_filename)
 plt.savefig(figure_filename)

# h5_filename = os.path.join(figure_savepath, file_basename + ".h5")

 #with h5py.File(h5_filename, 'w') as hf:
  #  hf.create_dataset("x_data", data=nPixels)  
   # hf.create_dataset("y_data", data=nT1) 
    #hf.create_dataset("z_data", data=numScans) 
 #plt.savefig(h5_filename)
 print("### Saving Figure ###")
 
plt.show()
    
      return  plt.savefig(figure_filename), plt.show()