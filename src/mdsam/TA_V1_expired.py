# Author : Simon Daneau
# Date : 2025/04/09

import os                       # type: ignore
import numpy as np              # type: ignore
import matplotlib.pyplot as plt # type: ignore
import h5py                     # type: ignore

from pathlib import Path

class Database:
    # Initialize the attributes with dimensions where necessary
    def __init__(self, NbScan, NbDelay, NbWavelength):
        self.Background = self.Background(NbScan, NbWavelength)
        self.Measurement = self.Measurement(NbScan, NbDelay, NbWavelength)
        self.Statistics = self.Statistics(NbScan, NbDelay, NbWavelength)

    class Background:
        # Initialize the attributes with dimensions where necessary
        def __init__(self, NbScan, NbWavelength):
            self.DelayTA = np.zeros(shape=(NbScan,), dtype=np.float32)                                  # ps
            self.DelayTB = np.zeros(shape=(NbScan,), dtype=np.float32)                                  # ps
            self.GeneralInformations = np.zeros(shape=(), dtype=np.str_)                                # ---
            self.NotPumpedSignalCount = np.zeros(shape=(NbScan,), dtype=np.int32)                       # cnt
            self.NotPumpedSignalDeviation = np.zeros(shape=(NbScan,), dtype=np.float32)                 # V ?
            self.NotPumpedSpectra = np.zeros(shape=(NbScan, NbWavelength), dtype=np.float32)            # V ?
            self.PumpPower = np.zeros(shape=(NbScan,), dtype=np.float32)                                # mW ?
            self.PumpedSignalCount = np.zeros(shape=(NbScan,), dtype=np.int32)                          # cnt
            self.PumpedSignalDeviation = np.zeros(shape=(NbScan,), dtype=np.float32)                    # V ?
            self.PumpedSpectra = np.zeros(shape=(NbScan, NbWavelength), dtype=np.float32)               # V ? 
            self.TotalSignalCount = np.zeros(shape=(NbScan,), dtype=np.int32)                           # cnt
            self.Wavelength = np.zeros(shape=(NbWavelength,), dtype=np.float32)                         # nm

    class Measurement:
        # Initialize the attributes with dimensions where necessary
        def __init__(self, NbScan, NbDelay, NbWavelength):
            self.AveragedDelay = np.zeros(shape=(NbDelay,), dtype=np.float32)                           # ps
            self.AveragedSignal = np.zeros(shape=(NbDelay, NbWavelength), dtype=np.float32)             # mOD ?
            self.DelayTA = np.zeros(shape=(NbScan, NbDelay), dtype=np.float32)                          # ps
            self.DelayTB = np.zeros(shape=(NbScan, NbDelay), dtype=np.float32)                          # ps
            self.GeneralInformations = np.zeros(shape=(), dtype=np.str_)                                # ---
            self.NotPumpedSignalCount = np.zeros(shape=(NbScan, NbDelay), dtype=np.int32)               # cnt
            self.NotPumpedSignalDeviation = np.zeros(shape=(NbScan, NbDelay), dtype=np.float32)         # V ?
            self.NotPumpedSpectra = np.zeros(shape=(NbScan, NbDelay, NbWavelength), dtype=np.float32)   # V ?
            self.PumpPower = np.zeros(shape=(NbScan, NbDelay), dtype=np.float32)                        # mW ?
            self.PumpedSignalCount = np.zeros(shape=(NbScan, NbDelay), dtype=np.int32)                  # cnt
            self.PumpedSignalDeviation = np.zeros(shape=(NbScan, NbDelay), dtype=np.float32)            # V ?
            self.PumpedSpectra = np.zeros(shape=(NbScan, NbDelay, NbWavelength), dtype=np.float32)      # V ?
            self.TotalSignalCount = np.zeros(shape=(NbScan, NbDelay), dtype=np.int32)                   # cnt
            self.Wavelength = np.zeros(shape=(NbWavelength,), dtype=np.float32)                         # nm

    class Statistics:
        # Initialize the attributes with dimensions where necessary
        def __init__(self, NbScan, NbDelay, NbWavelength):
            self.ChopperActualFrequency = np.zeros(shape=(), dtype=np.float32)                          # Hz
            self.ChopperActualPhaseError = np.zeros(shape=(), dtype=np.float32)                         # rad
            self.ChopperClockSource = np.zeros(shape=(), dtype=np.str_)                                 # ---
            self.ChopperLocked = np.zeros(shape=(), dtype=np.str_)                                      # ---
            self.CryostatTemperature = np.zeros(shape=(), dtype=np.float32)                             # K
            self.DetectorMotorizedVNDFTransmittance = np.zeros(shape=(), dtype=np.float32)              # V ?
            self.GeneralInformations = np.zeros(shape=(), dtype=np.str_)                                # ---
            self.GratingBlazeWavelength = np.zeros(shape=(), dtype=np.int32)                            # nm
            self.GratingGrooves = np.zeros(shape=(), dtype=np.float32)                                  # gr/nm
            self.HarpiaMotorizedVNDFTransmittance = np.zeros(shape=(), dtype=np.float32)                # V ?
            self.LaserTriggerClockDivider = np.zeros(shape=(), dtype=np.int32)                          # ---
            self.MeasurementFinished = np.zeros(shape=(), dtype=np.str_)                                # ---
            self.MeasurementStarted = np.zeros(shape=(), dtype=np.str_)                                 # ---
            self.NbDelay = np.zeros(shape=(), dtype=np.int32); self.NbDelay = NbDelay                   # cnt
            self.NbScan = np.zeros(shape=(), dtype=np.int32); self.NbScan = NbScan                      # cnt
            self.NbWavelength = np.zeros(shape=(), dtype=np.int32); self.NbWavelength = NbWavelength    # cnt
            self.NumberRun = np.zeros(shape=(), dtype=np.int32)                                         # cnt
            self.NumberSpectraBackground = np.zeros(shape=(), dtype=np.int32)                           # cnt
            self.ProbeMotorizedVNDFTransmittance = np.zeros(shape=(), dtype=np.float32)                 # V ?
            self.PumpLaserFrequency = np.zeros(shape=(), dtype=np.int32)                                # Hz
            self.PumpMotorizedVNDFTransmittance = np.zeros(shape=(), dtype=np.float32)                  # V ?
            self.SampleMoverActive = np.zeros(shape=(), dtype=np.str_)                                  # ---
            self.SampleName = np.zeros(shape=(), dtype=np.str_)                                         # ---
            self.TargetFrequency = np.zeros(shape=(), dtype=np.int32)                                   # Hz

class TransientAbsorption:

    def __init__(self):
        pass

    def get_data(self, FilesPath, MeasurementDataFileName, MeasurementStatisticsFileName, AveragedOutputFileName, SampleName, Temperature):
        
        # Open the files
        with open(os.path.join(FilesPath, MeasurementDataFileName)) as MeasurementDataFile:
            MeasurementDataLines = MeasurementDataFile.readlines()

        with open(os.path.join(FilesPath, MeasurementStatisticsFileName)) as MeasurementStatisticsFile:
            MeasurementStatisticsLines = MeasurementStatisticsFile.readlines()
        
        with open(os.path.join(FilesPath, AveragedOutputFileName)) as AveragedOutputFile:
            AveragedOutputLines = AveragedOutputFile.readlines()

        # Get some dimensions
        NbScan = int(MeasurementStatisticsLines[7][16:-1])-1
        NbDelay = len(AveragedOutputLines)-3
        NbWavelength = len(MeasurementStatisticsLines[19][12:-1].split('\t'))

        # Define the classes
        Data = Database(NbScan, NbDelay, NbWavelength)

        # Get the data from Measurement Statistics
        Data.Statistics.GeneralInformations = MeasurementStatisticsLines[0][:-1]
        Data.Statistics.SampleName = SampleName
        Data.Statistics.CryostatTemperature = Temperature
        Data.Statistics.MeasurementStarted = MeasurementStatisticsLines[1][21:-1]
        Data.Statistics.LaserTriggerClockDivider = MeasurementStatisticsLines[2][29:-1]
        Data.Statistics.ChopperLocked = MeasurementStatisticsLines[3][19:-1]
        Data.Statistics.GratingBlazeWavelength = MeasurementStatisticsLines[4][30:-1]
        Data.Statistics.GratingGrooves = MeasurementStatisticsLines[5][23:-1]
        Data.Statistics.NumberSpectraBackground = MeasurementStatisticsLines[6][49:-1]
        Data.Statistics.NumberRun = MeasurementStatisticsLines[7][16:-1]
        Data.Statistics.PumpLaserFrequency = MeasurementStatisticsLines[8][22:-4]
        Data.Statistics.PumpMotorizedVNDFTransmittance = MeasurementStatisticsLines[9][35:-1]
        Data.Statistics.ProbeMotorizedVNDFTransmittance = MeasurementStatisticsLines[10][36:-1]
        Data.Statistics.DetectorMotorizedVNDFTransmittance = MeasurementStatisticsLines[11][39:-1]
        Data.Statistics.HarpiaMotorizedVNDFTransmittance = MeasurementStatisticsLines[12][40:-1]
        LineStr = MeasurementStatisticsLines[13][:-1].split(', ')
        Data.Statistics.ChopperClockSource = LineStr[0][22:]
        Data.Statistics.TargetFrequency = LineStr[1][17:-3]
        Data.Statistics.ChopperActualFrequency = MeasurementStatisticsLines[14][26:-4]
        Data.Statistics.ChopperActualPhaseError = MeasurementStatisticsLines[15][28:-5]
        Data.Statistics.SampleMoverActive = MeasurementStatisticsLines[16][21:-1]
        Data.Statistics.MeasurementFinished = MeasurementStatisticsLines[-1][22:-1]
        
        Data.Background.GeneralInformations = MeasurementStatisticsLines[17][12:-1]
        Data.Background.Wavelength = [np.float32(i) for i in MeasurementStatisticsLines[19][12:-1].split('\t')]

        Data.Measurement.GeneralInformations = MeasurementStatisticsLines[18][13:-1]
        Data.Measurement.Wavelength = [np.float32(i) for i in MeasurementStatisticsLines[19][12:-1].split('\t')]

        Delay = 0
        for i in range(20, len(MeasurementStatisticsLines)-6, 6):
            LineStr = MeasurementStatisticsLines[i][:-1].split(', ')
            Type = LineStr[0]
            Scan = np.int32(LineStr[1][-1])-1
            if Type == "Background":
                Data.Background.DelayTA[Scan] = MeasurementStatisticsLines[i][26:42]
                Data.Background.DelayTB[Scan] = MeasurementStatisticsLines[i][43:-1]
                Data.Background.PumpPower[Scan] = MeasurementStatisticsLines[i+1][5:-1]
                LineStr = MeasurementStatisticsLines[i+2][:-1].split(', ')
                Data.Background.PumpedSignalCount[Scan] = LineStr[0][18:]
                Data.Background.PumpedSignalDeviation[Scan] = LineStr[2][22:]
                LineStr = MeasurementStatisticsLines[i+4][:-1].split(', ')
                Data.Background.NotPumpedSignalCount[Scan] = LineStr[0][21:]
                Data.Background.NotPumpedSignalDeviation[Scan] = LineStr[2][25:]
                Data.Background.TotalSignalCount[Scan] = LineStr[1][17:]
            elif Type == "Measurement":
                Data.Measurement.DelayTA[Scan][Delay] = MeasurementStatisticsLines[i][33:49]
                Data.Measurement.DelayTB[Scan][Delay] = MeasurementStatisticsLines[i][50:-1]
                Data.Measurement.PumpPower[Scan][Delay] = MeasurementStatisticsLines[i+1][5:-1]
                LineStr = MeasurementStatisticsLines[i+2][:-1].split(', ')
                Data.Measurement.PumpedSignalCount[Scan][Delay] = LineStr[0][18:]
                Data.Measurement.PumpedSignalDeviation[Scan][Delay] = LineStr[2][22:]
                LineStr = MeasurementStatisticsLines[i+4][:-1].split(', ')
                Data.Measurement.NotPumpedSignalCount[Scan][Delay] = LineStr[0][21:]
                Data.Measurement.NotPumpedSignalDeviation[Scan][Delay] = LineStr[2][25:]
                Data.Measurement.TotalSignalCount[Scan][Delay] = LineStr[1][17:]
                Delay = Delay+1
                if Delay == NbDelay: Delay = 0

        # Get the data from Measurement Data
        Delay = 0
        for i in range(4, len(MeasurementDataLines)-6, 6):
            LineStr = MeasurementDataLines[i][:-1].split(', ')
            Type = LineStr[0]
            Scan = np.int32(LineStr[1][-1])-1
            if Type == "Background":
                Data.Background.PumpedSpectra[Scan][:] = [np.float32(j) for j in MeasurementDataLines[i+2][:].split('\t')]
                Data.Background.NotPumpedSpectra[Scan][:] = [np.float32(j) for j in MeasurementDataLines[i+4][:].split('\t')]
            elif Type == "Measurement":
                Data.Measurement.PumpedSpectra[Scan][Delay][:] = [np.float32(j) for j in MeasurementDataLines[i+2][:].split('\t')]
                Data.Measurement.NotPumpedSpectra[Scan][Delay][:] = [np.float32(j) for j in MeasurementDataLines[i+4][:].split('\t')]
                Delay = Delay+1
                if Delay == NbDelay: Delay = 0

        # Get the data from Averaged Output
        Delay = 0
        for i in range(3, len(AveragedOutputLines), 1):
            LineStr = AveragedOutputLines[i][:-1].split('    ') # It's four or five spaces...
            Data.Measurement.AveragedDelay[Delay] = LineStr[1]
            Data.Measurement.AveragedSignal[Delay][:] = [np.float32(j) for j in LineStr[2:]]
            Delay = Delay+1

        return Data
    
    def plot_averaged(self, Data):
        
        # Define parameters
        Wavelength = Data.Measurement.Wavelength
        AveragedDelay = Data.Measurement.AveragedDelay
        AveragedSignal = Data.Measurement.AveragedSignal
        SampleName = Data.Statistics.SampleName
        Temperature = Data.Statistics.CryostatTemperature

        # Plot the 2D map 
        fig, ax = plt.subplots()
        pcolor = ax.pcolormesh(Wavelength, AveragedDelay, AveragedSignal, cmap = "RdBu")
        cbar = fig.colorbar(pcolor, ax = ax)
        cbar.set_label('Difference absorption', rotation = 90)
        ax.set_title("Transient absorption: "+str(SampleName)+", T = "+str(Temperature)+" K", fontweight = "bold")
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel("Delay [ps]")
        plt.show()

    def get_classes(self, HDF5FilePath, HDF5FileName):

        # Read HDF5
        HDF5File = h5py.File(HDF5FilePath+HDF5FileName, "r")

        # Get dimensions
        NbScan = HDF5File['Statistics']['NbScan'][()]
        NbDelay = HDF5File['Statistics']['NbDelay'][()]
        NbWavelength = HDF5File['Statistics']['NbWavelength'][()]

        # Get the classes
        Data = Database(NbScan, NbDelay, NbWavelength)

        return Data

    def HDF5File_name(self, Data):
        
        # Define parameters
        SampleName = Data.Statistics.SampleName
        Hour = Data.Statistics.MeasurementStarted[11:-4].replace(":", "-")
        Temperature = round(Data.Statistics.CryostatTemperature, 1)
        PumpPower = round(np.mean(Data.Measurement.PumpPower[Data.Measurement.PumpPower != 0]), 3)

        # Define HDF5File name
        HDF5FileName = str(SampleName)+"_H-M-S_"+str(Hour)+"_TransientAbsorption_Temperature_"+str(Temperature)+"K_PumpPower_"+str(PumpPower)+"mW.hdf5"

        return HDF5FileName

class GeneralFunctions:

    def create_HDF5(self, Data, HDF5FilePath):
        
        # Generate HDF5File name
        HDF5FileName = TransientAbsorption.HDF5File_name(self, Data)

        # Check if the file already exists
        if Path(HDF5FilePath+HDF5FileName).exists():
            answer = input("\nThe file already exists, do you want to overwrite it? [Yes/No] ")
            if answer == "No": return
            elif answer == "Yes": os.remove(HDF5FilePath+HDF5FileName)
            else: return

        # Create HDF5
        HDF5File = h5py.File(HDF5FilePath+HDF5FileName, "w")

        # Create groups and create dataset
        for i in range(0, len(dir(Data))-27, 1):
            group_str = dir(Data)[i]
            group_obj = getattr(Data, group_str)
            group = HDF5File.create_group(group_str)
            for j in range(0, len(dir(group_obj))-27, 1):
                dataset_str = dir(group_obj)[j]
                group.create_dataset(dir(group_obj)[j], data = getattr(group_obj, dataset_str))    

        HDF5File.close()

    def read_HDF5(self, HDF5FilePath, HDF5FileName):
        
        # Read HDF5
        HDF5File = h5py.File(HDF5FilePath+HDF5FileName, "r")

        # Define the classes
        Data = TransientAbsorption.get_classes(self, HDF5FilePath, HDF5FileName)

        # Recover the database
        for i in range(0, len(HDF5File), 1):
            group_str = list(HDF5File)[i]
            group_obj = getattr(Data, group_str)
            for j in range(0, len(HDF5File[group_str]), 1):
                dataset_str = list(HDF5File[group_str])[j]
                setattr(group_obj, list(HDF5File[group_str])[j], HDF5File[group_str][dataset_str][()])

        return Data