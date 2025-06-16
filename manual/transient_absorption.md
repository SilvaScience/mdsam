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
# MeasurementStatisticsFileName = 'fgl530-fp-test_TA_0_stats.dat'
# AveragedOutputFileName = 'fgl530-fp-test_TA_0_matrix.dat'
SampleName = 'fgl530'

# Optional: add any custom metadata as key-value pairs
extra_metadata = {
    "User": "Simon Daneau",
    "Comment": "Test scan at low temperature with fine alignment",
    "PumpPower": "50uW",
    "Temperature": "5K"}

# Import required modules from mdsam
from mdsam import Harpia

# Convert the data into an HDF5 file, i.e. run the __init__ function and plot.
TA = Harpia.TransientAbsorption(
    FilesPath,
    MeasurementDataFileName,
    SampleName,
    extra_metadata)