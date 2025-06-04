# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 14:02:25 2025

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
from PT_data_workup_functions import mlList
from PT_data_workup_functions import backgroundCorrect
from PT_data_workup_functions import hammingWindow
from PT_data_workup_functions import rowWiseLoop
polyFitRange = list()

def process_data(FTdata, bkgdCorrect, polyFitRange, polyFitOrder, nT1, numScans, fftLength=0, fftAxis=0, apodizeData= True):
    
    
    # Time domain calculations
    data = FTdata
    
    # Background correction
    if bkgdCorrect:
        print("#### Background Correction ####")
        data = rowWiseLoop(data,backgroundCorrect,polyFitRange = polyFitRange,
                              polyFitOrder = polyFitOrder)
    
    # Apodization
    if apodizeData:
        print("### Apodization ###")
        print("\tCalculating windowing function")
        H = hammingWindow(nT1, data)
        print("\tApodizing data")
        data = np.multiply(H,data)
    
    #  averaging
    print("### Averaging Scans ###")
    if len(data.shape) == 3:
        avgTF = np.mean(data[...,np.arange(numScans)], axis = 2)
    elif len(data.shape) == 2:
        print("\tData has only 2 dimensions. Not averaging data.")
    else:
        print("\tUnexpected shape of data:",data.shape)
        print("\tUnable to take average")
    
    # Fourier transform
    # calculate the FFT length, if requested
    if fftLength == 0:
        fftLength = 2**(np.ceil(np.log2(nT1))+1)
    axisList = ['column','row','page']
    avgTF[0,...] *= 0.5 # see Hamm and Zanni section 9.5.3
    print("### Fourier Transform ###")
    print("\tFFT Length:", fftLength)
    print("\tZero-padding by an additional factor of 2")
    # See Hamm & Zanni section 9.5.4 regarding FFT Length
    FF = np.real(np.fft.rfft(avgTF, n = int(2*fftLength), axis = fftAxis))
    FF = FF[:-1,...]
    
    return axisList, FF, avgTF, data