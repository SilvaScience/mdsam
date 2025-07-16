# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 15:29:17 2020

company: PhaseTech Spectroscopy, Inc. <http://phasetechspectroscopy.com/>
author: Thomas Brinzer <support@phasetechspectroscopy.com>
description: helper code with a few functions

"""

import numpy as np # for array and matrix calculations

def mlList(x1,x2):
    '''Makes a numerical list in MatLab array style
    
    Parameters
    ----------
    x1 int
        lower bound (inclusive)
    x2 int
        upper bound (inclusive)
    
    Returns
    -------
    a list
    '''
    a = list(range(x1, x2 + 1))
    return a

# FUNCTION DEFINITIONS. SCRIPT BEGINS AT LINE 123
def backgroundCorrect(bkgdCorrectArray,polyFitRange,polyFitOrder):
    ''' Performs a polynomial background correction on a row vector
    
    Performs a numpy.polyfit on the selected regions of an input row vector.
    Regions are selected by passing chosen indices (polyFitRange) to
    backgroundCorrect(). The resulting polynomial is then subtracted from the
    input vector, and the (background-corrected) result is returned.
    
    Parameters
    ----------
    bkgdCorrectArray: numpy.ndarray
    polyFitRange: list
        desired indices to use for background correction
    polyFitOrder: int
        order of desired polynomial fit
        
    Returns
    -------
    correctedArray: numpy.ndarray
    '''
    # selecting a subset of the input array for background correction
    polyFitArray = bkgdCorrectArray[polyFitRange]
    # polynomial fit
    polyCoeffs = np.polyfit(polyFitRange,polyFitArray,polyFitOrder) 
    correctedArray = (bkgdCorrectArray - 
                      np.polyval(polyCoeffs, np.arange(len(bkgdCorrectArray))))
    return correctedArray


def hammingWindow(nT1, dataIn):
    ''' Creates a Hamming window for apodizing the 2D data R(t1,w3;t2)
    
    Parameters
    ----------
    nT1: int
        number of elements of the desired time window (t >= 0)
        
    Returns
    -------
    windowArray: numpy.ndarray
    '''
    # calculate full Hamming window (-1 accounts for only having one 0)
    fullHammingWindow = np.hamming(2*len(np.linspace(0,1,nT1)) - 1) 
    zeroIndex = int(np.floor(len(fullHammingWindow)/2)) # Find t=0 index 
    h = fullHammingWindow[zeroIndex:] # select t>=0
    
    if len(dataIn.shape) == 3:
        [m,n,p] = dataIn.shape
        H = np.repeat(h[:,np.newaxis], n, axis = 1)
        H = np.repeat(H[:,:,np.newaxis], p, axis = 2)
    elif len(dataIn.shape) == 2:
        [m,n] = dataIn.shape
        H = np.repeat(h[:,np.newaxis], n, axis = 1)
    else:
        return
    return H

def rowWiseLoop(loopDataIn,functionIn,**kwargs):
    ''' Apply a 1d array function to each row of a 3D matrix
    
    Applies a function that operates on a 1d array to each row of an 2- or 3d
    input matrix. The resulting 1d array is stored in an (m*n'*p) output matrix
    
    Parameters
    ----------
    loopDataIn: numpy.ndarray
    functionHandle: function
    **kwargs: see below
    
    Returns
    -------
    loopDataOut: numpy.ndarray
    
    Notes
    -----
    **kwargs (keyword arguments) are the additional input arguments that
        are required for each 1d array function to operate.
        
    '''    
    if len(loopDataIn.shape) == 2: # handles a case with a 2D matrix
        loopDataIn = np.expand_dims(loopDataIn, axis=2)
    [m,n,p] = loopDataIn.shape
    firstIteration = True
    for j in range(p):
        print('\tProcessing scan',j + 1,'of',p)
        for i in range(m):
            resultArray = functionIn(loopDataIn[i,:,j],**kwargs)
            if firstIteration:
                loopDataOut = np.zeros((m,len(resultArray),p))
                firstIteration = False
            loopDataOut[i,:,j] = resultArray
    loopDataOut = np.squeeze(loopDataOut)
    return loopDataOut