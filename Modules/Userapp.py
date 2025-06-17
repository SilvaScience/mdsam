# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 15:33:55 2025

@author: marfe
"""

from Main import Main_

main = Main_(fileDirectory=r"C:\Users\marfe\Downloads\example_data_processing_code_python_and_matlab",
    dataSet=["20250312", "33"],
    waitingTimeNum= 9)
main.Analysis()