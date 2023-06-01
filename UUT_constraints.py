# -*- coding: utf-8 -*-
"""
Created on Wed May 31 20:35:51 2023

@author: Ryan Larson
"""

import numpy as np
# import pygad
# import logging
from shapely.geometry import Point, Polygon
from shapely.prepared import prep
import matplotlib.pyplot as plt
import random


def probe_type(type_selection):
    """
    Select a probe type and get the probe parameters.

    Parameters
    ----------
    type_selection : string
        String input with name of probe type from the following list:
            100mil
            75mil
            50mil
            39mil

    Raises
    ------
    ValueError
        Prints text specifying an invalid type was selected.

    Returns
    -------
    probe_diameter : float
        Probe diameter in inches.
    probe_min_ctc_spc : float
        Minimum probe center to center spacing, in inches.
    probe_force : float
        Probe force in ounces.

    """
    if type_selection == "100mil":
        probe_diameter = 0.036
        probe_min_ctc_spc = 0.1
        probe_force = 5.5
    elif type_selection == "75mil":
        probe_diameter = 0.025
        probe_min_ctc_spc = 0.075
        probe_force = 5.5
    elif type_selection == "50mil":
        probe_diameter = 0.021
        probe_min_ctc_spc = 0.050
        probe_force = 5.5
    elif type_selection == "39mil":
        probe_diameter = 0.015
        probe_min_ctc_spc = 0.039
        probe_force = 3.6
    else:
        raise ValueError("ERROR: Invalid probe type selection")
    
    return probe_diameter, probe_min_ctc_spc, probe_force


def pressure_rod_type(type_selection):
    if type_selection == "0.25":
        rod_diameter_top = 0.25
        rod_diameter_tip = 0.0625
    else:
        raise ValueError("ERROR: Invalid probe type selection")
        
    return rod_diameter_top, rod_diameter_tip


# pressure_rod_ctc = 0.375    # Pressure rod center to center spacing, inches

def max_pressure_rods(UUT_poly, pressure_rod_ctc):
    # Get the bounding box of the UUT polygon exterior
    xmin, ymin, xmax, ymax = UUT_poly.bounds
    
    offset_pcts = [0.0, 0.33, 0.5, 0.67]
    offset_vals = [offset*pressure_rod_ctc for offset in offset_pcts]
    
    max_rods = 0
    for x_offset in offset_vals:
        for y_offset in offset_vals:
            rod_pts = []
            for x in np.arange(xmin + x_offset, xmax, pressure_rod_ctc):
                for y in np.arange(ymin + y_offset, ymax, pressure_rod_ctc):
                    rod_pts.append(Point(x, y))
                
            # Validate that each point falls inside shape
            valid_pts = []
            valid_pts.extend(filter(UUT_poly.contains, rod_pts))
            n_rods = len(valid_pts)
            if n_rods > max_rods:
                max_rods = n_rods
                xoff = x_offset
                yoff = y_offset
    
    return max_rods, xoff, yoff
    
    
            
        