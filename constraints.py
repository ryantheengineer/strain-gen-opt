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
import xml.etree.ElementTree as ET
import pickle
import pathlib

# def xml_to_dict(initialdir):
#     filetypes = (("XML", ["*.xml"]),) 
#     root = Tk()
#     root.wm_attributes('-topmost', 1)
#     inputfile = fd.askopenfilename(
#             title="Select FEA input file",
#             initialdir=initialdir,
#             filetypes=filetypes
#             )
#     root.destroy()
    
#     with open(inputfile, 'r', encoding='utf-8') as file:
#         UUT_xml = file.read()    
#     UUT_dict = xmltodict.parse(UUT_xml)
#     return UUT_dict

def parse_input_xml(initialdir):
    filetypes = (("XML", ["*.xml"]),) 
    root = Tk()
    root.wm_attributes('-topmost', 1)
    inputfile = fd.askopenfilename(
            title="Select FEA input file",
            initialdir=initialdir,
            filetypes=filetypes
            )
    root.destroy()
    
    
    # Parse the XML file
    tree = ET.parse(inputfile)
    root = tree.getroot()
    
    # Define empty lists to store region and hole vertices
    region_vertices = []
    hole_vertices = []
    
    # Iterate over the Panel structs
    data = []
    structs = root.findall('.//struct')
    # for struct in structs:
    #     if struct.attrib.get("identifier") == "Panel":
    
    polyshapes = root.findall('.//polyshape')
    polyshape_ids = []
    for polyshape in polyshapes:
        polyshape_ids.append(polyshape.attrib.get("identifier"))
        
    # Gather data into usable form
    if "pBoards" in polyshape_ids:
        pBoards = root.findall('.//polyshape[@identifier="pBoards"]')
        regions = pBoards[0].findall('.//region')
        pBoards_regions = []
        for region in regions:
            # Extract vertices for the region
            vertices = []
            for vertex in region.findall('.//vertex'):
                x, y = map(float, vertex.text.split('|'))
                vertices.append([x, y])
            pBoards_regions.append(vertices)

        holes = pBoards[0].findall('.//hole')
        pBoards_holes = []
        for hole in holes:
            # Extract vertices for the hole
            vertices = []
            for vertex in hole.findall('.//vertex'):
                x, y = map(float, vertex.text.split('|'))
                vertices.append([x, y])
            pBoards_holes.append(vertices)
            
        # Create pandas DataFrames
        if len(pBoards_regions) == 0:
            pBoards_regions_df = None
        else:
            pBoards_regions_df = pd.DataFrame(pBoards_regions, columns=['Vertices'])
        if len(pBoards_holes) == 0:
            pBoards_holes_df = None
        else:
            pBoards_holes_df = pd.DataFrame(pBoards_holes, columns=['Vertices'])
    else:
        # pBoards_regions = None
        pBoards_regions_df = None
        pBoards_holes_df = None
        
        
    return pBoards_regions_df, pBoards_holes_df
            
            
    #     # Iterate over the pBoards in the panel
    #     for pboard in panel.findall('.//pBoards'):
    #         # Iterate over the regions in the pBoard
    #         for region in pboard.findall('.//Region'):
    #             # Extract vertices for the region
    #             vertices = []
    #             for vertex in region.findall('.//Vertex'):
    #                 x, y = map(float, vertex.text.split('|'))
    #                 vertices.append((x, y))
    #             region_vertices.append(vertices)
    
    #         # Iterate over the holes in the pBoard
    #         for hole in pboard.findall('.//Hole'):
    #             # Extract vertices for the hole
    #             vertices = []
    #             for vertex in hole.findall('.//Vertex'):
    #                 x, y = map(float, vertex.text.split('|'))
    #                 vertices.append((x, y))
    #             hole_vertices.append(vertices)
    
    # # Create pandas DataFrames
    # regions_df = pd.DataFrame(region_vertices, columns=['Vertices'])
    # holes_df = pd.DataFrame(hole_vertices, columns=['Vertices'])
    
    # # Print the DataFrames
    # print("Regions DataFrame:")
    # print(regions_df)
    
    # print("\nHoles DataFrame:")
    # print(holes_df)


if __name__ == "__main__":
    # Load previously chosen FEA path here
    filename = 'FEApath.pk'
    with open(filename, 'rb') as fi:
        FEApath = pickle.load(fi)
    
    initialdir = str(pathlib.Path(FEApath).parent) + "Examples"
    pBoards_regions_df, pBoards_holes_df = parse_input_xml(initialdir)