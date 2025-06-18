import numpy as np
import h5py
import os

class PhaseTech2DIRData:
    '''
        Insert description of class here
    '''
    def __init__(self,dataSetName,fileDirectory):
        '''
            describe what this function does (example in https://github.com/SilvaScience/colberto/blob/main/src/compute/beams.py)
        '''
        self.fileDirectory = fileDirectory
        self.dataSet = dataSetName
        # function that scans the directory and finds all the waiting time numbers
        DataImporter(fileDirectory=self.fileDirectory,dataSet=self.dataSet)
        #self.waitingTimeNum =waitingTimeNum 
        
class DataImporter:
    def __init__(self, fileDirectory, dataSet, waitingTimeNum=None):
        self.fileDirectory = os.path.abspath(fileDirectory)
        self.dataSet = dataSet
        if waitingTimeNum is None:
            self.find_waiting_times()
        else:
            self.waitingTimeNums = [waitingTimeNum]
        
    def find_waiting_times(self):
        self.waitingTimeNums=[]
        for filename in os.listdir(path):
            #if re.search(pattern, filename):
            if self.dataSet in filename:
                m = re.search('T(.+?)#', filename)
                if m:
                    found = int(m.group(1))
                    if not found in waitingTimeNums:
                        self.waitingTimeNums.append(found)
                    
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
            sys.exit("Stopping script execution.")
        
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