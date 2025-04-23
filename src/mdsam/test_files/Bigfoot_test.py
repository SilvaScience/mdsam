# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 16:52:53 2025

@author: David Tiede
"""

import os                       
import numpy as np              
import matplotlib.pyplot as plt 
import h5py                    
import matplotlib
import sys
import inspect
from helper_functions import utils as ut
from matplotlib.ticker import FormatStrFormatter
from BF import Bigfoot as bf
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
#from TA import * 

def parse_mdcs_header_to_dict(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    result = {
        "headers": [],
        "emission_energy": [],
        "stepped_axis_energy": [],
        "population_scan_steps": [],
        "value": None,
        "scan_params": {},
        "scan_notes": []
    }

    # Step 1: Extract headers
    result["headers"] = lines[0].strip().replace("Rows:", "").split("\t")

    # Step 2: Read data sections
    mode = 1
    for idx, line in enumerate(lines[1:], start=1):
        #print(line)
        stripped = line.strip()
        if not stripped:
            continue

        if mode == 1:
            # Emission and stepped axis energy (same block)
            if stripped.startswith('-'):  # start of population scan steps
                result["stepped_axis_energy"].extend(map(float, stripped.split()))
                mode =2
            else:
                result["emission_energy"].extend(map(float, stripped.split()))

        if mode == 2:
            # Population scan steps
            if stripped.startswith("Scan type") or '\t' not in stripped:
                try:
                    result["population_scan_steps"].extend(map(float, stripped.split()))
                except ValueError:
                    # Possibly the single value
                    try:
                        result["value"] = float(stripped)
                    except ValueError:
                        pass
                    mode = 3
            else:
                mode = 3

        if mode == 3:
            # Scan parameters
            if stripped.startswith("Scan type"):
                keys = stripped.split("\t")
                values = lines[idx + 1].strip().split("\t")
                result["scan_params"] = dict(zip(keys, values))
                mode = 4

        if mode == 4:
            # Scan notes
            if stripped.startswith("Scan Notes:"):
                mode = 5

        if mode == 5:
            result["scan_notes"].append(stripped)

    return result


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

"""
# get .hdf5 files in directory
#files = [file for file in os.listdir(data_folder) if file.find('.hdf5')!= -1]
files_raw = [file for file in os.listdir(data_folder) if file.find('.tsv')!= -1]
header_file = [file for file in os.listdir(data_folder) if file.find('Header.txt')!= -1][0]
subfolder = [file for file in os.listdir(data_folder) if file.find('_Averages')!= -1][0]
subfolder_files = os.listdir(data_folder + subfolder)


header_data = parse_mdcs_header_to_dict(data_folder + header_file)

filenames = ['2DSpec_amp','2DSpec_phase','Linear_amp','TimeSpec_amp','TimeSpec_phase']
subfolder_filenames = ['2DSpec_amp','2DSpec_phase','Linear_amp','RefAB_pos0']
data = {}
for file in files_raw:
    for f in filenames: 
        if file.find(f)!= -1:
            data[f] = np.loadtxt(data_folder + file)
avg = list(range(len(subfolder_files)))
for i,file in enumerate(subfolder_files):
    print(file)
    avg[i] = file[file.find('avg'):file.find('.tsv')]
avg = np.unique(avg)
for av in avg:
    data[av]= {}
for file in subfolder_files: 
    for f in subfolder_filenames: 
        if file.find(f)!= -1:
            av = file[file.find('avg'):file.find('.tsv')]
            data[av][f] = np.loadtxt(data_folder + subfolder +'\\' + file)
            
amp_2D_file = [file for file in files_raw if file.find('2DSpec_amp_T0')!= -1][0]
phase_2D_file = [file for file in files_raw if file.find('2DSpec_phase_T0')!= -1][0]
amp_2D_data = np.loadtxt(data_folder + amp_2D_file)
phase_2D_data = np.loadtxt(data_folder + phase_2D_file)

range_min = float(header_data['scan_params']['Plot range min (units)'])
range_max = float(header_data['scan_params']['Plot range max (units)'])

em_axis = np.array(header_data['emission_energy'])
ex_axis = np.array(header_data['stepped_axis_energy'])

em_range = np.r_[ut.find_idx(em_axis,range_min):ut.find_idx(em_axis,range_max)]
ex_range = np.r_[ut.find_idx(ex_axis,-range_max):ut.find_idx(ex_axis,-range_min)]
max_idx = ut.find_idx(em_axis,range_max)-9 


amp_2D_data = amp_2D_data[em_range[0]:em_range[-1],ex_range[0]:ex_range[-1]]
phase_2D_data = phase_2D_data[em_range[0]:em_range[-1],ex_range[0]:ex_range[-1]]
em_axis = em_axis[em_range[:-1]]
ex_axis = ex_axis[ex_range[:-1]]

mat = amp_2D_data
matreal = amp_2D_data*np.sin(phase_2D_data)
matimag = amp_2D_data*np.cos(phase_2D_data)

"""

# %% Plot data 
ims = {} # Dict for colorbars 
plot_axis = [1545,1565]
f, plts = plt.subplots(1,3,figsize=(16*ut.cm,6*ut.cm))
#for i in range(3): plts[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

### Plot Data ###
# plot Real 
norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(matreal), vcenter=0, vmax=np.max(matreal))
ims[0] = plts[0].contourf(em_axis, ex_axis, np.transpose(matreal), levels=30, linewidths=0.5, cmap="RdBu_r", norm=norm) #/np.min(matreal)
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


f.suptitle(f'GaAs QW', y=1.02)
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

"""
# %%
f = plt.figure(figsize = (3,3))
norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(amp_2D_data),vcenter= 5, vmax=np.max(amp_2D_data))
plt.contourf(em_axis, ex_axis, np.transpose(amp_2D_data), levels=30, linewidths=0.5, cmap="Reds", norm=norm ) #/np.min(matreal)
#plt.pcolor(em_axis, ex_axis, np.transpose(amp_2D_data), cmap="viridis") #/np.min(matreal)
plt.colorbar()
plt.contour(em_axis, ex_axis, np.transpose(amp_2D_data)/np.min(amp_2D_data), levels=10, colors='k', linewidths=0.2, norm=norm)
lims = [1545,1565]
plt.xlim(lims)
plt.ylim(-lims[1],-lims[0])
plt.plot(lims,[-lims[0],-lims[1]], ex_axis, linestyle="-.", color="grey")
plt.xlabel('Emission Energy (eV)')
plt.ylabel('Excitation Energy (eV)')

"""