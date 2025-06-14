# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 16:47:16 2025

@author: David Tiede
"""

import os                       
import numpy as np              
import matplotlib.pyplot as plt 
import h5py        
import re             

from pathlib import Path
from mdsam.helper_functions import HelperFunctions  as ut

class Bigfoot:
                       
    def transform_to_HDF5(data_folder):
        """ Transform Bigfoot data to .hdf5 file 
        
        Args:
            data_folder (str): relative or full path to the data as saved from Bigfoot software
        
        Returns: .hdf5 file 
        """
        # check if provided data_path is absolute
        p = Path(data_folder)
        if not p.is_absolute():
            cwd = os.getcwd()
            data_folder = cwd + data_folder
        
        # append to data_folder to make access files insider folder 
        data_folder = data_folder +  '\\'
        
        #### load data ####
        # find data files in data folder
        files_raw = [file for file in os.listdir(data_folder) if file.find('.tsv')!= -1]
        if files_raw == []:
            print('WARNING: No files found in data folder!')
        header_file = [file for file in os.listdir(data_folder) if file.find('Header.txt')!= -1][0]
        subfolder = [file for file in os.listdir(data_folder) if file.find('_Averages')!= -1][0]
        subfolder_files = os.listdir(data_folder + subfolder)
        
        #prepare output filename
        num_run = header_file[-14:-11]
        sample_name = header_file[:-14]
        date = Path(data_folder).parent.name
        output_filenname = sample_name + '_' + date + '_' + num_run + '_Bigfoot.h5'
        output_path = str(Path(data_folder).parent) + '\\' + output_filenname

        # process header data in separate function
        header = Bigfoot_helper_functions.parse_mdcs_header_to_dict(data_folder + header_file)
        header['file_id'] = output_filenname[:-3]

        # order files in data dictionary 
        filenames = ['2DSpec_amp','2DSpec_phase','Linear_amp','TimeSpec_amp','TimeSpec_phase']
        subfolder_filenames = ['2DSpec_amp','2DSpec_phase','Linear_amp','RefAB_pos0']
        # check if 3D or 2D dataset
        if len(header['population_scan_steps'])>1: # 3D dataset 
            
            # preallocate data structure 
            data = {}
            for val in header['population_scan_steps']:
                T_pop = round(val,2)
                data[T_pop] = {}
                data[T_pop]['raw'] = {}
                for i in range(int(float(header['scan_params']['# averages']))):
                    data[T_pop]['avg' + str(i)]= {}
                
            # get T_pop times 
            T_pop_array = np.unique(list(data.keys()))
            
            # store raw data
            for file in files_raw: # treat averaged data
                numbers_in_filename = re.findall(r'\d+', file)
                T_pop = T_pop_array[int(numbers_in_filename[-1])]
                for f in filenames: 
                    if file.find(f)!= -1:
                        data[T_pop]['raw'][f] = np.loadtxt(data_folder + file)
            
            # store averages
            for i,file in enumerate(subfolder_files):
                numbers_in_filename = re.findall(r'\d+', file)
                T_pop = T_pop_array[int(numbers_in_filename[-2])]
                avg_num = int(numbers_in_filename[-1])
                for f in subfolder_filenames: 
                    if file.find(f)!= -1:
                        data[T_pop]['avg' + str(avg_num)][f] = np.loadtxt(data_folder + subfolder +'\\' + file)

                        
            #### process data ####
            # extract valid range from header
            range_min = float(header['scan_params']['Plot range min (units)'])
            range_max = float(header['scan_params']['Plot range max (units)'])
    
            # get 2D axis from header
            em_axis = np.array(header['emission_energy'])
            ex_axis = np.array(header['stepped_axis_energy'])
    
            # find indices for valid range 
            em_range = np.r_[ut.find_idx(em_axis,range_min):ut.find_idx(em_axis,range_max)]
            ex_range = np.r_[ut.find_idx(ex_axis,-range_max):ut.find_idx(ex_axis,-range_min)] 
            
            # calculate 2D spec for every T_pop
            for T_pop  in T_pop_array:
                # crop 2D data to valid range
                amp_2D_data = data[T_pop]['raw']['2DSpec_amp'][em_range[0]:em_range[-1],ex_range[0]:ex_range[-1]]
                phase_2D_data = data[T_pop]['raw']['2DSpec_phase'][em_range[0]:em_range[-1],ex_range[0]:ex_range[-1]]
                
                # store processed data
                data[T_pop]['em_axis'] = em_axis[em_range[:-1]]
                data[T_pop]['ex_axis'] = ex_axis[ex_range[:-1]]
                data[T_pop]['2Dabs'] = amp_2D_data
                data[T_pop]['2Dreal'] = amp_2D_data*np.sin(phase_2D_data)
                data[T_pop]['2Dimag'] = amp_2D_data*np.cos(phase_2D_data)
        else:
            data = {}
            data['raw'] = {}
            for file in files_raw: # treat averaged data
                for f in filenames: 
                    if file.find(f)!= -1:
                        data['raw'][f] = np.loadtxt(data_folder + file)
            avg = list(range(len(subfolder_files))) # order individual scans in avg sub-directories
            for i,file in enumerate(subfolder_files):
                avg[i] = file[file.find('avg'):file.find('.tsv')]
            avg = np.unique(avg)
            for av in avg:
                data[av]= {}
            for file in subfolder_files: 
                for f in subfolder_filenames: 
                    if file.find(f)!= -1:
                        av = file[file.find('avg'):file.find('.tsv')]
                        data[av][f] = np.loadtxt(data_folder + subfolder +'\\' + file)
                        
            #### process data ####
            # extract valid range from header
            range_min = float(header['scan_params']['Plot range min (units)'])
            range_max = float(header['scan_params']['Plot range max (units)'])
    
            # get 2D axis from header
            em_axis = np.array(header['emission_energy'])
            ex_axis = np.array(header['stepped_axis_energy'])
    
            # find indices for valid range 
            em_range = np.r_[ut.find_idx(em_axis,range_min):ut.find_idx(em_axis,range_max)]
            ex_range = np.r_[ut.find_idx(ex_axis,-range_max):ut.find_idx(ex_axis,-range_min)] 
    
            # crop 2D data to valid range
            amp_2D_data = data['raw']['2DSpec_amp'][em_range[0]:em_range[-1],ex_range[0]:ex_range[-1]]
            phase_2D_data = data['raw']['2DSpec_phase'][em_range[0]:em_range[-1],ex_range[0]:ex_range[-1]]
            
            # store processed data
            data['em_axis'] = em_axis[em_range[:-1]]
            data['ex_axis'] = ex_axis[ex_range[:-1]]
            data['2Dabs'] = amp_2D_data
            data['2Dreal'] = amp_2D_data*np.sin(phase_2D_data)
            data['2Dimag'] = amp_2D_data*np.cos(phase_2D_data)
        
        #### save to HDF5 ####
        with h5py.File(output_path, 'w') as hdf:
            # Save all numerical arrays
            for key, value in data.items():
                if isinstance(value, np.ndarray):
                    hdf.create_dataset(f"data/{key}", data=value)
                elif isinstance(value, dict):
                    grp = hdf.create_group(f"data/{key}")
                    for subkey, subval in value.items():
                        if isinstance(subval, np.ndarray):
                            grp.create_dataset(subkey, data=subval)
                        elif isinstance(subval, dict):
                            subgrp = hdf.create_group(f"data/{key}/{subkey}")
                            for subsubkey, subsubval in subval.items():
                                subgrp.create_dataset(subsubkey, data=subsubval)
                                
            # Save header information
            header_grp = hdf.create_group("header")
            for key, val in header.items():
                if isinstance(val, (str, float, int)):
                    header_grp.attrs[key] = val
                elif isinstance(val, (list, np.ndarray)):
                    arr = np.array(val)
                    if arr.dtype.kind in {'U', 'S'}:  # it's a string or unicode string array
                        header_grp.attrs[key] = [str(x) for x in arr]
                    else:
                        header_grp.create_dataset(key, data=arr)
                elif isinstance(val, dict):
                    subgrp = header_grp.create_group(key)
                    for sk, sv in val.items():
                        try:
                            subgrp.attrs[sk] = sv
                        except TypeError:
                            pass  # Skip complex/unserializable types
        
        print(f"Data and header saved to {output_path}")
        return output_path
        
        
    def read_HDF5_file(file_path):
        """Read data and header from a Bigfoot HDF5 file.
        
        Args:
            file_path (str): Path to the .h5 file
        
        Returns:
            data (dict): Nested dictionary of numerical datasets
            header (dict): Metadata dictionary
        """
        data = {}
        header = {}
        
        with h5py.File(file_path, 'r') as hdf:
            # Load data. Check at every level of nested dict whether item is data or subgroup
            data_group = hdf['data']
            for key in data_group:
                item = data_group[key]
                if isinstance(item, h5py.Dataset):
                    data[key] = item[()]
                elif isinstance(item, h5py.Group):
                    data[key] = {}
                    for subkey in item:
                        if isinstance(item[subkey], h5py.Group):
                            data[key][subkey] = {}
                            for subsubkey in item[subkey].keys():
                                data[key][subkey][subsubkey] = item[subkey][subsubkey][()]
                        else:
                            data[key][subkey] = item[subkey][()]
        
            # Load header
            header_group = hdf['header']
            for key in header_group.attrs:
                header[key] = header_group.attrs[key]
        
            for key in header_group:
                item = header_group[key]
                if isinstance(item, h5py.Dataset):
                    header[key] = item[()]
                elif isinstance(item, h5py.Group):
                    header[key] = {}
                    for subkey in item.attrs:
                        header[key][subkey] = item.attrs[subkey]
                    for subkey in item:
                        header[key][subkey] = item[subkey][()]
        
        return data, header
        
                    
                    
class Bigfoot_helper_functions: 
    
   def parse_mdcs_header_to_dict(filepath):
       """Parses all content of the header file.
       
       Args:
           file_path (str): Path to the header.txt file
       
       Returns:
           header (dict): Metadata dictionary
       """
       with open(filepath, 'r') as f:
           lines = f.readlines()

       result = {
           "headers": [],
           "emission_energy": [],
           "stepped_axis_energy": [],
           "population_scan_steps": [],
           "scan_params": {},
           "scan_notes": []
       }

       # Step 1: Extract headers
       result["headers"] = lines[0].strip().replace("Rows:", "").split("\t")[1:] # [1:] removes first empty entry

       # Step 2: Read data sections line by line and append accordingly to dict
       for idx, line in enumerate(lines[1:], start=1):
           stripped = line.strip()
           if not stripped:
               continue

           if idx == 1:
               result["emission_energy"].extend(map(float, stripped.split()))
           elif idx == 2:
               result["stepped_axis_energy"].extend(map(float, stripped.split()))
           elif idx == 3:
               result["population_scan_steps"].extend(map(float, stripped.split()))
           elif idx == 4:
               keys = stripped.split("\t")
               values = lines[idx + 1].strip().split("\t")
               result["scan_params"] = dict(zip(keys, values))
           elif idx == 5:
               pass # already inclued in scan parameters 
           elif idx > 5: 
               if stripped.startswith("Scan Notes:"):
                   stripped = stripped[12:]  # [12:] removes "scan notes: " comment
               result["scan_notes"].append(stripped)
              
       return result
    

