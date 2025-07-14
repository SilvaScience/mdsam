# Author : Simon Daneau
# Date : 2025/05/04

import os                       # type: ignore
import numpy as np              # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib               # type: ignore
import h5py                     # type: ignore
import re                       # type: ignore
import tkinter as tk            # type: ignore

from tkinter import filedialog
from pathlib import Path

class TransientAbsorption:

    @staticmethod
    def create_database_dict(self):
        return {
            "NumberDelay": self.NbDelay,
            "NumberScan": self.NbScan,
            "NumberWavelength": self.NbWavelength,
            "AveragedDelay": np.zeros((self.NbDelay,), dtype=np.float32),
            "AveragedSignal": np.zeros((self.NbDelay, self.NbWavelength), dtype=np.float32),
            "Wavelength": np.zeros((self.NbWavelength,), dtype=np.float32),
            "Background": {
                "DelayTA": np.zeros((self.NbScan,), dtype=np.float32),
                "DelayTB": np.zeros((self.NbScan,), dtype=np.float32),
                "GeneralInformations": np.zeros((), dtype=np.str_),
                "NotPumpedSignalCount": np.zeros((self.NbScan,), dtype=np.int32),
                "NotPumpedSignalDeviation": np.zeros((self.NbScan,), dtype=np.float32),
                "NotPumpedSpectra": np.zeros((self.NbScan, self.NbWavelength), dtype=np.float32),
                "PumpPower": np.zeros((self.NbScan,), dtype=np.float32),
                "PumpedSignalCount": np.zeros((self.NbScan,), dtype=np.int32),
                "PumpedSignalDeviation": np.zeros((self.NbScan,), dtype=np.float32),
                "PumpedSpectra": np.zeros((self.NbScan, self.NbWavelength), dtype=np.float32),
                "TotalSignalCount": np.zeros((self.NbScan,), dtype=np.int32),},
            "Measurement": {
                "DelayTA": np.zeros((self.NbScan, self.NbDelay), dtype=np.float32),
                "DelayTB": np.zeros((self.NbScan, self.NbDelay), dtype=np.float32),
                "GeneralInformations": np.zeros((), dtype=np.str_),
                "NotPumpedSignalCount": np.zeros((self.NbScan, self.NbDelay), dtype=np.int32),
                "NotPumpedSignalDeviation": np.zeros((self.NbScan, self.NbDelay), dtype=np.float32),
                "NotPumpedSpectra": np.zeros((self.NbScan, self.NbDelay, self.NbWavelength), dtype=np.float32),
                "PumpPower": np.zeros((self.NbScan, self.NbDelay), dtype=np.float32),
                "PumpedSignalCount": np.zeros((self.NbScan, self.NbDelay), dtype=np.int32),
                "PumpedSignalDeviation": np.zeros((self.NbScan, self.NbDelay), dtype=np.float32),
                "PumpedSpectra": np.zeros((self.NbScan, self.NbDelay, self.NbWavelength), dtype=np.float32),
                "TotalSignalCount": np.zeros((self.NbScan, self.NbDelay), dtype=np.int32),},
            "Statistics": {
                "ChopperActualFrequency": np.zeros((), dtype=np.float32),
                "ChopperActualPhaseError": np.zeros((), dtype=np.float32),
                "ChopperClockSource": np.zeros((), dtype=np.str_),
                "ChopperLocked": np.zeros((), dtype=np.str_),
                "DetectorMotorizedVNDFTransmittance": np.zeros((), dtype=np.float32),
                "GeneralInformations": np.zeros((), dtype=np.str_),
                "GratingBlazeWavelength": np.zeros((), dtype=np.int32),
                "GratingGrooves": np.zeros((), dtype=np.float32),
                "HarpiaMotorizedVNDFTransmittance": np.zeros((), dtype=np.float32),
                "LaserTriggerClockDivider": np.zeros((), dtype=np.int32),
                "MeasurementFinished": np.zeros((), dtype=np.str_),
                "MeasurementStarted": np.zeros((), dtype=np.str_),
                "NumberRun": np.zeros((), dtype=np.int32),
                "NumberSpectraBackground": np.zeros((), dtype=np.int32),
                "ProbeMotorizedVNDFTransmittance": np.zeros((), dtype=np.float32),
                "PumpLaserFrequency": np.zeros((), dtype=np.int32),
                "PumpMotorizedVNDFTransmittance": np.zeros((), dtype=np.float32),
                "SampleMoverActive": np.zeros((), dtype=np.str_),
                "SampleName": np.zeros((), dtype=np.str_),
                "TargetFrequency": np.zeros((), dtype=np.int32),}}
    
    def __init__(self, FilesPath, MeasurementDataFileName, SampleName, extra_metadata=None):

        self.Data = {}
        self.FilesPath = FilesPath
        self.BaseName = MeasurementDataFileName
        self.SampleName = SampleName
        self.extra_metadata = extra_metadata
    
    def transform_to_HDF5(self):

        # --- Read data files ---
        def read_lines_from_file(path):
            with open(path) as f:
                return f.readlines()
            
        BaseName = os.path.splitext(self.BaseName)[0] # Removes .dat
        MeasurementStatisticsFileName = f'{BaseName}_stats.dat'
        AveragedOutputFileName = f'{BaseName}_matrix.dat'

        MeasurementDataLines = read_lines_from_file(os.path.join(self.FilesPath, self.BaseName))
        MeasurementStatisticsLines = read_lines_from_file(os.path.join(self.FilesPath, MeasurementStatisticsFileName))
        AveragedOutputLines = read_lines_from_file(os.path.join(self.FilesPath, AveragedOutputFileName))

        # --- Initialize storage parameters ---
        self.NbDelay = len(AveragedOutputLines) - 2 - 1  # Number of delay points (skip 2 header lines)
        self.NbScan = int(MeasurementStatisticsLines[7][16:-1]) - 1  # Number of scans (0-indexed)
        WavelengthStr = MeasurementStatisticsLines[19][12:-1].split('\t')
        self.NbWavelength = len(WavelengthStr) # Number of wavelength channels

        # --- Create data structure ---
        self.Data = self.create_database_dict(self)
        self.Data["Wavelength"] = [np.float32(w) for w in WavelengthStr]  # Store wavelength values

        # --- Parse averaged output ---
        for delay_index in range(self.NbDelay):
            LineStr = re.split(r'\s+', AveragedOutputLines[delay_index + 3].strip())  # Tokenize line
            self.Data["AveragedDelay"][delay_index] = LineStr[0]  # Delay value
            self.Data["AveragedSignal"][delay_index] = [np.float32(value) for value in LineStr[1:]]  # Signal values as float32

        # --- Parse measurement data file ---
        # Extract general metadata for background and measurement
        self.Data["Background"]["GeneralInformations"] = MeasurementStatisticsLines[17][12:-1]
        self.Data["Measurement"]["GeneralInformations"] = MeasurementStatisticsLines[18][13:-1]

        delay_index = 0  # Keeps track of the delay step for measurement entries

        # Iterate through measurement data blocks (6 lines per block)
        for i in range(20, len(MeasurementStatisticsLines) - 6, 6):
            line = MeasurementStatisticsLines[i][:-1].split(', ')
            measurement_type = line[0]                          # "Background" or "Measurement"
            scan_index = int(line[1][-1]) - 1                   # Extract scan number (1-based to 0-based)

            # Helper function to extract and store parsed values into Data dict
            def store_measurement_data(base, is_background):
                line2 = MeasurementStatisticsLines[i + 2][:-1].split(', ')
                line4 = MeasurementStatisticsLines[i + 4][:-1].split(', ')

                # Indexing shortcut
                data = self.Data[base]

                # Parse shared fields
                if is_background:
                    data["DelayTA"][scan_index]               = MeasurementStatisticsLines[i][26:42]
                    data["DelayTB"][scan_index]               = MeasurementStatisticsLines[i][43:-1]
                    data["PumpPower"][scan_index]             = MeasurementStatisticsLines[i + 1][5:-1]
                    data["PumpedSignalCount"][scan_index]     = line2[0][18:]
                    data["PumpedSignalDeviation"][scan_index] = line2[2][22:]
                    data["NotPumpedSignalCount"][scan_index]  = line4[0][21:]
                    data["NotPumpedSignalDeviation"][scan_index] = line4[2][25:]
                    data["TotalSignalCount"][scan_index]      = line4[1][17:]
                else:
                    data["DelayTA"][scan_index][delay_index]               = MeasurementStatisticsLines[i][33:49]
                    data["DelayTB"][scan_index][delay_index]               = MeasurementStatisticsLines[i][50:-1]
                    data["PumpPower"][scan_index][delay_index]             = MeasurementStatisticsLines[i + 1][5:-1]
                    data["PumpedSignalCount"][scan_index][delay_index]     = line2[0][18:]
                    data["PumpedSignalDeviation"][scan_index][delay_index] = line2[2][22:]
                    data["NotPumpedSignalCount"][scan_index][delay_index]  = line4[0][21:]
                    data["NotPumpedSignalDeviation"][scan_index][delay_index] = line4[2][25:]
                    data["TotalSignalCount"][scan_index][delay_index]      = line4[1][17:]

            # Dispatch based on measurement type
            if measurement_type == "Background":
                store_measurement_data("Background", is_background=True)
            elif measurement_type == "Measurement":
                store_measurement_data("Measurement", is_background=False)
                delay_index += 1
                if delay_index == self.NbDelay:
                    delay_index = 0  # Reset delay index after one full cycle

        # --- Parse measurement data lines (spectra) ---
        delay_index = 0  # Tracks current delay step for measurement entries

        # Process every 6-line block starting from line 4
        for i in range(4, len(MeasurementDataLines) - 6, 6):
            line = MeasurementDataLines[i][:-1].split(', ')
            measurement_type = line[0]                        # "Background" or "Measurement"
            scan_index = int(line[1][-1]) - 1                 # Convert scan number to 0-based index

            # Convert a tab-separated string of numbers to a float32 numpy list
            def parse_spectrum(line_str):
                return [np.float32(value) for value in line_str.split('\t')]

            if measurement_type == "Background":
                self.Data["Background"]["PumpedSpectra"][scan_index][:] = parse_spectrum(MeasurementDataLines[i + 2])
                self.Data["Background"]["NotPumpedSpectra"][scan_index][:] = parse_spectrum(MeasurementDataLines[i + 4])

            elif measurement_type == "Measurement":
                self.Data["Measurement"]["PumpedSpectra"][scan_index][delay_index][:] = parse_spectrum(MeasurementDataLines[i + 2])
                self.Data["Measurement"]["NotPumpedSpectra"][scan_index][delay_index][:] = parse_spectrum(MeasurementDataLines[i + 4])
                delay_index += 1
                if delay_index == self.NbDelay:
                    delay_index = 0  # Reset after completing all delay points

        # Get the data from Measurement Statistics
        self.Data["Statistics"]["GeneralInformations"] = MeasurementStatisticsLines[0][:-1]
        self.Data["Statistics"]["SampleName"] = self.SampleName
        self.Data["Statistics"]["MeasurementStarted"] = MeasurementStatisticsLines[1][21:-1]
        self.Data["Statistics"]["LaserTriggerClockDivider"] = MeasurementStatisticsLines[2][29:-1]
        self.Data["Statistics"]["ChopperLocked"] = MeasurementStatisticsLines[3][19:-1]
        self.Data["Statistics"]["GratingBlazeWavelength"] = MeasurementStatisticsLines[4][30:-1]
        self.Data["Statistics"]["GratingGrooves"] = MeasurementStatisticsLines[5][23:-1]
        self.Data["Statistics"]["NumberSpectraBackground"] = MeasurementStatisticsLines[6][49:-1]
        self.Data["Statistics"]["NumberRun"] = MeasurementStatisticsLines[7][16:-1]
        self.Data["Statistics"]["PumpLaserFrequency"] = MeasurementStatisticsLines[8][22:-4]
        self.Data["Statistics"]["PumpMotorizedVNDFTransmittance"] = MeasurementStatisticsLines[9][35:-1]
        self.Data["Statistics"]["ProbeMotorizedVNDFTransmittance"] = MeasurementStatisticsLines[10][36:-1]
        self.Data["Statistics"]["DetectorMotorizedVNDFTransmittance"] = MeasurementStatisticsLines[11][39:-1]
        self.Data["Statistics"]["HarpiaMotorizedVNDFTransmittance"] = MeasurementStatisticsLines[12][40:-1]

        LineStr = MeasurementStatisticsLines[13][:-1].split(', ')
        self.Data["Statistics"]["ChopperClockSource"] = LineStr[0][22:]
        self.Data["Statistics"]["TargetFrequency"] = LineStr[1][17:-3]

        self.Data["Statistics"]["ChopperActualFrequency"] = MeasurementStatisticsLines[14][26:-4]
        self.Data["Statistics"]["ChopperActualPhaseError"] = MeasurementStatisticsLines[15][28:-5]
        self.Data["Statistics"]["SampleMoverActive"] = MeasurementStatisticsLines[16][21:-1]
        self.Data["Statistics"]["MeasurementFinished"] = MeasurementStatisticsLines[-1][22:-1]

        # Inject extra metadata if provided
        if self.extra_metadata:
            for key, value in self.extra_metadata.items():
                self.Data["Statistics"][key] = value

        # Generate HDF5 file name using f-string for better readability
        Hour = self.Data["Statistics"]["MeasurementStarted"][11:-4].replace(":", "-")
        HDF5FileName = f"{self.SampleName}_H-M-S_{Hour}_TransientAbsorption.hdf5"

        # Create HDF5
        #HDF5Helper.save_to_hdf5_with_prompt(Data, HDF5FileName)
        HDF5Helper.save_to_hdf5(self.Data, self.FilesPath, HDF5FileName)

    def plot_averaged(self, save_png=False):
        """
        Plots the 2D map of averaged signal with respect to wavelength and delay.

        Parameters:
        - Data: Dictionary containing the measurement and statistics data
        - save_png (bool): If True, open dialog to save plot as PNG
        """
        
        # Extracting necessary data
        Wavelength = self.Data['Wavelength']
        AveragedDelay = self.Data['AveragedDelay']
        AveragedSignal = self.Data['AveragedSignal']
        SampleName = self.Data['Statistics']['SampleName']

        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 6))
        pcolor = ax.pcolormesh(Wavelength, AveragedDelay, AveragedSignal, cmap="RdBu", shading='auto')

        cbar = fig.colorbar(pcolor, ax=ax)
        cbar.set_label('Differential absorption', rotation=90)

        ax.set_title(f"Transient absorption: {SampleName}", fontweight="bold")
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel("Delay [ps]")

        plt.tight_layout()

        # Save to PNG without triggering the GUI backend
        if save_png:
            # Switch backend to Agg for non-interactive saving
            matplotlib.use("Agg")
            root = tk.Tk()
            root.withdraw()

            filepath = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
                title="Save plot as PNG"
            )
            if filepath:
                fig.savefig(filepath, dpi=300)
                print(f"Plot saved to {filepath}")
            else:
                print("Save cancelled.")

            plt.close(fig)  # Close the figure after saving to avoid GUI issues
        else:
            plt.show()

class HDF5Helper:

    @staticmethod
    def save_to_hdf5_with_prompt(data, default_filename="data.h5"):
        """
        Open a dialog to choose filename and save HDF5 file.

        Parameters:
        - data (dict): The data to save
        - default_filename (str): Suggested default file name
        """
       
        # Initialize and hide the Tkinter root window
        root = tk.Tk()
        root.withdraw()

        # Ask user to choose filename and location
        filepath = filedialog.asksaveasfilename(
            title="Save HDF5 file as",
            defaultextension=".hdf5",
            initialfile=default_filename,
            filetypes=[("HDF5 files", "*.hdf5 *.h5"), ("All files", "*.*")])

        if not filepath:
            print("Save cancelled.")
            return

        # If file exists, remove it
        if os.path.exists(filepath):
            os.remove(filepath)

        # Write the data to the HDF5 file
        with h5py.File(filepath, 'w') as h5f:
            HDF5Helper._recursively_save(h5f, '', data)

        print(f"Data saved to {filepath}")

    @staticmethod
    def save_to_hdf5(data, filepath, filename):
        """
        Save the provided data to an HDF5 file.

        Parameters:
        - data (dict): The data to save.
        - filepath (str): Directory where the file will be saved.
        - filename (str): File name (with or without extension).
        """

        # Ensure the directory exists
        os.makedirs(filepath, exist_ok=True)

        # Add default extension if missing
        base, ext = os.path.splitext(filename)
        if ext == '':
            ext = '.h5'
        full_path = os.path.join(filepath, base + ext)

        # If file exists, remove it
        if os.path.exists(full_path):
            os.remove(full_path)

        # Write the data to the HDF5 file
        with h5py.File(full_path, 'w') as h5f:
            HDF5Helper._recursively_save(h5f, '', data)

        print(f"Data saved to {full_path}")

    @staticmethod
    def _recursively_save(h5file, path, dic):
        """Recursively save a nested dictionary into an HDF5 file."""
        
        for key, item in dic.items():
            # Build the dataset path within the HDF5 structure
            key_path = f"{path}/{key}" if path else key
            if isinstance(item, dict):
                # Recurse if the value is another dictionary
                HDF5Helper._recursively_save(h5file, key_path, item)
            else:
                # Create intermediate group if needed
                group_path = os.path.dirname(key_path)
                if group_path and group_path not in h5file:
                    h5file.require_group(group_path)
                # Save the actual dataset
                h5file.create_dataset(key_path, data=item)

    @staticmethod
    def load_from_hdf5_prompt():
        """Open file dialog to load HDF5 data and return it as a nested dictionary."""
        
        # Initialize and hide the Tkinter root window
        root = tk.Tk()
        root.withdraw()

        # Ask user to select an HDF5 file
        full_path = filedialog.askopenfilename(
            title="Choose HDF5 file to open",
            filetypes=[("HDF5 files", "*.hdf5 *.h5"), ("All files", "*.*")])
        if not full_path:
            print("Load cancelled.")
            return None
        
        # Split into directory and filename
        filepath, filename = os.path.split(full_path)

        # Load and return the data
        return HDF5Helper.load_from_hdf5(filepath, filename)

    @staticmethod
    def load_from_hdf5(filepath, filename):
        """
        Load HDF5 data from the specified directory and filename, and return it as a nested dictionary.

        Parameters:
        - filepath (str): Directory where the file is located.
        - filename (str): Name of the HDF5 file (with or without extension).

        Returns:
        - dict: Nested dictionary of loaded data, or None if the file doesn't exist.
        """

        # Add default extension if missing
        base, ext = os.path.splitext(filename)
        if ext == '':
            ext = '.h5'
        full_path = os.path.join(filepath, base + ext)

        # Check file existence
        if not os.path.exists(full_path):
            print(f"File not found: {full_path}")
            return None

        # Load and return the data
        with h5py.File(full_path, 'r') as h5f:
            return HDF5Helper._recursively_load(h5f)

    @staticmethod
    def _recursively_load(h5group):
        """Recursively load data from an HDF5 group into a nested dictionary."""
        
        result = {}
        for key, item in h5group.items():
            if isinstance(item, h5py.Group):
                # Recurse into sub-groups
                result[key] = HDF5Helper._recursively_load(item)
            elif isinstance(item, h5py.Dataset):
                # Read the dataset
                data = item[()]
                # Decode bytes to string if necessary
                if isinstance(data, bytes):
                    data = data.decode('utf-8')
                result[key] = data
        return result
