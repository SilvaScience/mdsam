import numpy as np
import h5py
import os

class PhaseTech2DIRData:
    '''
        Insert description of class here
    '''
    def __init__(self,folder_path,filePrefix):
        '''
            describe what this function does (example in https://github.com/SilvaScience/colberto/blob/main/src/compute/beams.py)
        '''
        self.Data={}
        self.Data['FourierTransformed'],self.Data['TimeDomain'],self.Data['ProbeTransmission']=self.import_data(folder_path,filePrefix=filePrefix)

    def import_data(fileDirectory, filePrefix, waitingTimeNum, dataSet, numScans, start):
        """
            Add a description of how this 
        """
        fileDirectory = os.path.abspath(fileDirectory)
        fileList = listdir(fileDirectory) # list of all files in the directory
        filePrefix = "#".join(dataSet) + (waitingTimeNum > 0)*("_T" + "{:02}".format(waitingTimeNum))
        #filePrefix = "#".join(dataSet) + (waitingTimeNum > 0)*("_T" + "{:02}".format(waitingTimeNum)) + '#'
        fileExtension = ".scan"
        filteredList = [i for i in fileList if filePrefix in i]
        filteredList = [i for i in filteredList if fileExtension in i]
        totalScans = len(filteredList)
        
        # LOAD FILES, STORE IN MEMORY, AND DETERMINE THEIR SIZES
        if totalScans != 0:
            temp = np.loadtxt(os.path.join(fileDirectory,filteredList[0]))
            [nT1,nPixels] = temp.shape
            nT1 -= 1 # time points (ending row is not a time point)
            nPixels -= 1 # number of pixels (column 0 is time, not a pixel)
        else:
            raise FileNotFoundError("No data files found matching found matching prefix",
                filePrefix,"in directory",fileDirectory,"\n")
        
        # DETERMINING THE NUMBER OF FILES TO AVERAGE
        A = numScans
        if A == 0:
            A = totalScans - start
        if A > totalScans: # don't allow more scans than the total number of scans
            A = totalScans
        if A < 0: # negative numbers allow for averaging from end
            if start == 0:
                start = totalScans + A
            else:
                start = start + A
            if start < 0:
                start = 0
                A = -totalScans
            A = -A
        numScans = A
        scans2Avg = np.arange(numScans) # 0:numScans - 1
        
        # initialize data matrices
        FTdata = np.zeros([nT1,nPixels,numScans])
        
        
        
        # Data loading
        print("### Loading Scans ###")
        for i in scans2Avg:
            dataTemp = np.loadtxt(os.path.join(fileDirectory,filteredList[start + i]))
            print("\tLoading file", i + 1, "of", numScans)
            FTdata[...,i] = dataTemp[:-1,1:].copy()
        return FTdata, dataTemp, nT1

    def save(self):