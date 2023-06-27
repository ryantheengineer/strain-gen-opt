# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 23:28:58 2023

@author: Ryan Larson

Functions for interpreting initial XML design file and constructing a shapely
polygon around it.

THIS WILL BE THE FINAL LIB FILE FOR GENERATING CONSTRAINTS
FROM THE XML DESIGN INPUT.
"""

from shapely.geometry import Point, Polygon, MultiPolygon
import xmltodict
import pandas as pd
import numpy as np
import random
from tkinter import Tk
from tkinter import filedialog as fd

def xml_to_dict(initialdir):
    filetypes = (("XML", ["*.xml"]),) 
    root = Tk()
    root.wm_attributes('-topmost', 1)
    inputfile = fd.askopenfilename(
            title="Select FEA input file",
            initialdir=initialdir,
            filetypes=filetypes
            )
    root.destroy()
    
    with open(inputfile, 'r', encoding='utf-8') as file:
        UUT_xml = file.read()    
    UUT_dict = xmltodict.parse(UUT_xml)
    return UUT_dict

