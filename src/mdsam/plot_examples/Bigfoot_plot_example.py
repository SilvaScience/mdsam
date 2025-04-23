# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 23:50:11 2025

@author: David Tiede
"""

import os                       
import numpy as np              
import matplotlib.pyplot as plt                 
import matplotlib
import sys
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
from helper_functions import utils as ut
from matplotlib.ticker import FormatStrFormatter
from BF import Bigfoot as bf


# %% Load data 
# Get current directory 
cwd = os.getcwd() # get current working directory
project_folder = 'src'
project_data_idx = cwd.find(project_folder)
data_folder = cwd[:project_data_idx+len(project_folder) + 1] + 'Data\\Bigfoot_GaAs\\Run020' +  '\\'
output_path = cwd[:project_data_idx+len(project_folder) + 1] + 'Data\\Bigfoot_GaAs\\test_file.h5'

bf.transform_to_HDF5(data_folder,output_path)
data,header = bf.read_HDF5_file(output_path)

em_axis = data['em_axis']
ex_axis = data['ex_axis']

mat = data['2Dabs']
matreal = data['2Dreal']
matimag = data['2Dimag']

# %% Plot data 
ims = {} # Dict for colorbars 
plot_axis = [1545,1565]
f, plts = plt.subplots(1,3,figsize=(16*ut.cm,6*ut.cm))
#for i in range(3): plts[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

### Plot Data ###
# plot Real 
norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(matreal), vcenter=0, vmax=np.max(matreal))
ims[0] = plts[0].contourf(em_axis, ex_axis, np.transpose(matreal), levels=30, cmap="RdBu_r", norm=norm) #/np.min(matreal)
plts[0].contour(em_axis, ex_axis, np.transpose(matreal)/np.min(matreal), levels=10, colors='k', linewidths=0.2, norm=norm)
plts[0].set_ylabel("Absorption Energy (meV)")

plts[0].plot(plot_axis,[-plot_axis[0],-plot_axis[1]], ex_axis, linestyle="-.", color="grey")
plts[0].set_title('Real(A) Reph. ')


# plot Imag 
norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(matimag), vcenter=0, vmax=np.max(matimag))
ims[1] = plts[1].contourf(em_axis, ex_axis, np.transpose(matimag), levels=30, cmap="RdBu_r", norm=norm) #/np.min(matimag)
plts[1].contour(em_axis, ex_axis, np.transpose(matimag)/np.min(matimag), levels=10, colors='k', linewidths=0.2, norm=norm)

plts[1].plot(plot_axis,[-plot_axis[0],-plot_axis[1]], ex_axis, linestyle="-.", color="grey")
plts[1].set_title('Imag(A) Reph. ')



# plot Abs 
norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(mat), vcenter = np.max(mat)/2, vmax=np.max(mat))
ims[2] = plts[2].contourf(em_axis, ex_axis, np.transpose(mat), levels=30, cmap="Reds", norm=norm) #/np.max(mat)
plts[2].contour(em_axis, ex_axis, np.transpose(mat)/np.max(mat), levels=10, colors='k', linewidths=0.2, norm=norm)


plts[2].plot(plot_axis,[-plot_axis[0],-plot_axis[1]], ex_axis, linestyle="-.", color="grey")
plts[2].set_title('$|A|^2$ Reph. ')


f.suptitle('GaAs QW', y=1.04)
f.subplots_adjust(wspace =0)
for i in range(3): 
    plts[i].set_xlim(plot_axis)
    plts[i].set_ylim(-plot_axis[1],-plot_axis[0])
    plts[i].set_xlabel("Emission Energy (meV)")
    f.colorbar(ims[i], orientation='horizontal', pad= 0.2,shrink=0.9, format = '%.1f')
    plts[i].minorticks_on()
    plts[i].tick_params(axis='both', which='major', length=6, width=1, direction='in', bottom=True, top=True, left=True, right=True)
    plts[i].tick_params(axis='both', which='minor', length=3, width=1, direction='in', bottom=True, top=True, left=True, right=True)
for i in range(2): plts[i+1].yaxis.set_ticklabels([])

plt.show()
