# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:27:18 2024

@author: David Tiede
"""
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

class utils:

    cm = 1/2.54 
    
    def create_custom_colormap(colors, positions=None):
        if positions is None:
            positions = np.linspace(0, 1, len(colors))
        return LinearSegmentedColormap.from_list("", list(zip(positions, colors)))
    
    def create_brightness_colormap(base_color, name='brightness_colormap'):
        """
        Create a colormap that transitions from a given base color to white (increasing brightness).
        
        Parameters:
        - base_color: tuple of (R, G, B) values normalized between 0 and 1 (e.g., (0.2, 0.6, 0.8)).
        - name: str, the name of the colormap.
        
        Returns:
        - colormap: A custom LinearSegmentedColormap.
        """
        # Define the color points (start at the base color, end at white)
        r = base_color[0]
        g = base_color[1]
        b = base_color[2] 
        c2 = 0.6 
        #cdict = {
        #    'red':   [(0.0, 1, r+(1-r)*c2), (0.5, r+(1-r)*c, r+(1-r)*c), (1.0, r, 0.1)],
        #    'green': [(0.0, 1, g+(1-g)*c2), (0.5, g+(1-g)*c, g+(1-g)*c), (1.0, g, 0.1)],
        #    'blue':  [(0.0, 1, b+(1-b)*c2), (0.5, b+(1-b)*c, b+(1-b)*c), (1.0, b, 0.1)],
        #}
        cdict = {
            'red':   [(0.0, 1, r+(1-r)*c2), (1.0, r, 0.1)],
            'green': [(0.0, 1, g+(1-g)*c2), (1.0, g, 0.1)],
            'blue':  [(0.0, 1, b+(1-b)*c2), (1.0, b, 0.1)],
        }
        return LinearSegmentedColormap(name, cdict)
    
    
    def find_idx(array,value):
        return np.argmin(abs(array-value))
    
    def get_sec(time_str):
        """Get seconds from time."""
        h, m, s = time_str.split(':')
        return int(h) * 3600 + int(m) * 60 + int(s)
    
    def lin_fit(x,a,b):
        return a*x+b
    
    def nexpdec(x,*args):
        return_function = 0
        num_gauss = int(len(args)/2)
        for i in range(num_gauss):
            return_function += args[i]*np.exp(-args[num_gauss + i]*x)  
        return return_function
    
    def nexpdec_log(x,*args):
        return_function = 0
        num_gauss = int(len(args)/2)
        for i in range(num_gauss):
            return_function += (10**args[i])*np.exp(-(10**args[num_gauss + i])*x)  
        return return_function
    
    def nexpdec_deriv(x,*args):
        return_function = 0
        num_gauss = int(len(args)/2)
        for i in range(num_gauss):
            return_function -= args[num_gauss + i]*args[i]*np.exp(-args[num_gauss + i]*x)  
        return return_function
    
    def nexpdec_deriv_log(x,*args):
        return_function = 0
        num_gauss = int(len(args)/2)
        for i in range(num_gauss):
            return_function -= (10**args[num_gauss + i])*(10**args[i])*np.exp(-(10**args[num_gauss + i])*x)  
        return return_function
    
    def polyfits(x,*args):
        return_function = np.zeros(len(x))
        num_poly = len(args)
        for i in range(num_poly):
            return_function += args[i]*x**i
        return return_function
    
    def sci_print(number: float) -> str:
        sci_notation = f"{number:.1e}"
        base, exponent = sci_notation.split('e')
        return f"${base} \\cdot 10^{{{int(exponent)}}}$"

    def time_to_seconds(time_str):
        """
        Converts a time string in hh:mm:ss format to total seconds.
    
        Args:
        time_str (str): The time string in the format 'hh_mm_ss'.
    
        Returns:
        int: The total time in seconds.
        """
        try:
            # Split the time string into hours, minutes, and seconds
            hours, minutes, seconds = map(int, time_str.split("_"))
            # Convert to total seconds
            total_seconds = hours * 3600 + minutes * 60 + seconds
            return total_seconds
        except ValueError:
            raise ValueError("Invalid time format. Please use 'hh_mm_ss'.")

