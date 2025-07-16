# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 10:20:56 2025

@author: marfe
"""

#from mdsam.Phasetech2DIR import PhaseTech2DIRData

from Phasetech2DIR import PhaseTech2DIRData

"Steps: 1. Enter the appropiate filepath in fileDirectory 2. Enter the name of your file in dataSet"
"Ex. if your complete file name is: 20250312#33_T01#001 write dataSet=[20250312,33]"
main = PhaseTech2DIRData (fileDirectory=r"C:\Users\marfe\Downloads\example_data_processing_code_python_and_matlab",
    dataSet=["20250312", "33"])
main.Analysis()
