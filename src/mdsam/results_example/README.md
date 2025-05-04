# Transient Absorption Data Processing Example

This example demonstrates how to transform transient absorption data files into HDF5 format and visualize the averaged result using the `mdsam` library.

## Files

Make sure the following data files exist in your specified directory (`FilesPath`):

- `fgl530-fp-test_TA_0.dat` — raw measurement data  
- `fgl530-fp-test_TA_0_stats.dat` — measurement metadata/statistics  
- `fgl530-fp-test_TA_0_matrix.dat` — averaged signal output  

## Example Script

```python
# Define the file paths and sample parameters
FilesPath = '/Users/simondaneau/Desktop/'
MeasurementDataFileName = 'fgl530-fp-test_TA_0.dat'
MeasurementStatisticsFileName = 'fgl530-fp-test_TA_0_stats.dat'
AveragedOutputFileName = 'fgl530-fp-test_TA_0_matrix.dat'
SampleName = 'CuGeO3'
Temperature = 5.0

# Import required modules from mdsam
from mdsam.TA_V2 import TransientAbsorption, HDF5Helper

# Convert the data into an HDF5 file
TransientAbsorption.transform_to_HDF5(
    FilesPath,
    MeasurementDataFileName,
    MeasurementStatisticsFileName,
    AveragedOutputFileName,
    SampleName,
    Temperature
)

# Load the data from an HDF5 file (with a file selection dialog)
Data = HDF5Helper.load_from_hdf5_prompt()

# Plot the averaged signal (and save to PNG if desired)
TransientAbsorption.plot_averaged(Data, save_to_file=True)
```

## Notes

- When running this script, a dialog will appear to let you select the HDF5 file to load.
- The final plot will display the 2D transient absorption map and optionally save it to a PNG file.
