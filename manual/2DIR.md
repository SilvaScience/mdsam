**2DIR Data Processing Example**

This example demonstrates how to manage data files and visualize the graphs using the mdsam library.

**Files**

Make sure to download all data files in the route: *SilvaScience/mdsam/2DIR_Mareny/Test files* and save them  in your specified directory (FilesPath):

**Example Script**

from mdsam.Phasetech2DIR import PhaseTech2DIRData

"Steps: 1. Enter the appropiate filepath in fileDirectory, 2. Enter the name of your file in dataSet"
"Ex. if your complete file name is: 20250312#33_T01#001 write dataSet=[20250312,33]"

main = PhaseTech2DIRData (fileDirectory=r"./Test_files", dataSet=["20250312", "33"])
    
main.Analysis()

main.Graphs()
