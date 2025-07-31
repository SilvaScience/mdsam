# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 12:08:57 2025

@author: David Tiede
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mdsam.helper_functions import HelperFunctions  as ut
from matplotlib.ticker import MaxNLocator

class BFPlots:

    def plot_Reph_Re_Im_Abs(data,header,plot_range):
        """ Plots 2D rephasing spectra with real, imag and absolute 2D maps. 
            Axis limits are specified through range. 
        
        Args:
            data (dict): dict with data loaded from a BF .h5 file
            header (dict): dict with metadata loaded from a BF .h5 file
            plot range (list): list with two entries for min and max values of axis (in meV)
            
        """
        
        #get data from data dict
        try: 
            em_axis = data['em_axis']
            ex_axis = data['ex_axis']
        except KeyError as e:
            raise KeyError(f"Axis not found in the highest h5 hieracrcy. \
                           Check for 3D or 2D data. Missing expected key: {e}") 
            
            
    
        mat = data['2Dabs']
        matreal = data['2Dreal']
        matimag = data['2Dimag']
        
        # Initialize plot
        ims = {} # Dict for colorbars 
        plot_axis = plot_range
        f, plts = plt.subplots(1,3,figsize=(16*ut.cm,6*ut.cm))
        
        ### Plot Data ###
        # plot Real 
        norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(matreal), vcenter=0, vmax=np.max(matreal))
        ims[0] = plts[0].contourf(em_axis, ex_axis, np.transpose(matreal), levels=30, cmap="RdBu_r", norm=norm) 
        plts[0].contour(em_axis, ex_axis, np.transpose(matreal)/np.min(matreal), levels=10, colors='k', linewidths=0.2, norm=norm)
        plts[0].set_ylabel("Absorption Energy (meV)")
        plts[0].set_title('Real(A)')
        
        
        # plot Imag 
        norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(matimag), vcenter=0, vmax=np.max(matimag))
        ims[1] = plts[1].contourf(em_axis, ex_axis, np.transpose(matimag), levels=30, cmap="RdBu_r", norm=norm) #/np.min(matimag)
        plts[1].contour(em_axis, ex_axis, np.transpose(matimag)/np.min(matimag), levels=10, colors='k', linewidths=0.2, norm=norm)
        plts[1].set_title('Imag(A)')
    
    
    
        # plot Abs 
        norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(mat), vcenter = np.max(mat)/2, vmax=np.max(mat))
        ims[2] = plts[2].contourf(em_axis, ex_axis, np.transpose(mat), levels=30, cmap="Reds", norm=norm) #/np.max(mat)
        plts[2].contour(em_axis, ex_axis, np.transpose(mat)/np.max(mat), levels=10, colors='k', linewidths=0.2, norm=norm)
        plts[2].set_title('$|A|^2$')
        
        title = header['file_id'][:-8]
        f.suptitle(title, y=1.06)
        f.subplots_adjust(wspace =0)
        for i in range(3): 
            plts[i].plot(plot_axis,[-plot_axis[0],-plot_axis[1]], ex_axis, linestyle="-.", color="grey")
            plts[i].set_xlim(plot_axis)
            plts[i].set_ylim(-plot_axis[1],-plot_axis[0])
            plts[i].set_xlabel("Emission Energy (meV)")
            cbar = f.colorbar(ims[i], orientation='horizontal', pad= 0.2,shrink=0.9, format = '%.1f')
            cbar.locator = MaxNLocator(nbins=4)
            cbar.update_ticks()
            plts[i].minorticks_on()
            plts[i].tick_params(axis='both', which='major', length=6, width=1, direction='in', bottom=True, top=True, left=True, right=True)
            plts[i].tick_params(axis='both', which='minor', length=3, width=1, direction='in', bottom=True, top=True, left=True, right=True)
        for i in range(2): plts[i+1].yaxis.set_ticklabels([])
        plt.show()