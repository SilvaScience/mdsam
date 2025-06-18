# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 09:45:08 2025

@author: marfe
"""

from mdsam.phasetech2dir import PhaseTech2DIRData

fileDirectory=r"C:\Users\marfe\Downloads\example_data_processing_code_python_and_matlab",
dataSet="20250312"

data2dir=PhaseTech2DIRData(dataSetName,fileDirectory)
data2dir.plot_data()
data2dir.save_hdf5(savepath)



    
