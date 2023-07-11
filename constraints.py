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
import matplotlib.pyplot as plt

# Get the root of the XML tree that will be used for all other input parsing
def get_XML_tree(initialdir):
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
    treeroot = tree.getroot()
    return treeroot

########## Obtain relevant Panel info ##########
# pBoards
def get_pBoards(root):
    polyshapes = root.findall('.//polyshape')
    pBoards = None
    region_polys = None
    hole_polys = None
    for polyshape in polyshapes:
        if polyshape.attrib.get("identifier") == "pBoards":
            pBoards = polyshape
    
    def verts_to_poly(polys, elem):
        vertices_elem = elem.findall('.//vertex')
        vertices = [vertex.text.split("|") for vertex in vertices_elem]
        polys.append(Polygon(vertices))
        
        
    if pBoards:
        numRegions = int(pBoards.find('.//numRegions').text)
        numHoles = int(pBoards.find('.//numHoles').text)
        if numRegions:
            # Create the Polygons of each region
            regions = pBoards.findall('.//region')
            region_polys = []
            for region in regions:
                verts_to_poly(region_polys, region)
                # vertices_elem = region.findall('.//vertex')
                # vertices = [vertex.text.split("|") for vertex in vertices_elem]
                # region_polys.append(Polygon(vertices))
        if numHoles:
            # Create the Polygons of each hole
            holes = pBoards.findall('.//hole')
            hole_polys = []
            for hole in holes:
                verts_to_poly(hole_polys, hole)
                # vertices_elem = hole.findall('.//vertex')
                # vertices = [vertex.text.split("|") for vertex in vertices_elem]
                # hole_polys.append(Polygon(vertices))
    else:
        raise ValueError("No pBoards element found")
        
    # If there are holes, determine which region they belong to and subtract
    # them to make a list of complete polygons
    # NOTE: This assumes that all of the polygons produced by get_pBoards are
    # unique and will not overlap each other.
    if numHoles:
        holes_grouped = [[] for poly in region_polys]
        for hole_poly in hole_polys:
            for i,region_poly in enumerate(region_polys):
                if region_poly.contains(hole_poly):
                    holes_grouped[i].append(hole_poly)
                    continue
                
        # Subtract the holes from their regions
        for i in range(len(region_polys)):
            for hole_poly in holes_grouped[i]:
                region_polys[i] = region_polys[i].difference(hole_poly)
    
    return region_polys, numRegions, numHoles


# Utility functions
def plot_poly_list_w_holes(poly_list):
    fig, ax = plt.subplots(figsize=(10,8),dpi=300)
    ax.set_aspect('equal')
    for poly in poly_list:
        xe, ye = poly.exterior.xy
        plt.plot(xe, ye, color="blue")
        for inner in poly.interiors:
            xi, yi = zip(*inner.coords[:])
            plt.plot(xi, yi, color="blue")


def get_region_list_centroids(region_polys):
    region_centroids = [poly.centroid for poly in region_polys]
    return region_centroids



def get_table_rows(root, identifier):
    tables = root.findall('.//table')
    rows = None
    for table in tables:
        if table.attrib.get("identifier") == identifier:
            # Extract rows
            rows = []
            for row in table.findall('.//row'):
                rows.append(row.text.split('|'))
    return rows
                        


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


# def parse_input_xml(initialdir):
#     filetypes = (("XML", ["*.xml"]),) 
#     root = Tk()
#     root.wm_attributes('-topmost', 1)
#     inputfile = fd.askopenfilename(
#             title="Select FEA input file",
#             initialdir=initialdir,
#             filetypes=filetypes
#             )
#     root.destroy()
    
    
#     # Parse the XML file
#     tree = ET.parse(inputfile)
#     root = tree.getroot()
    
#     # Define empty lists to store region and hole vertices
#     region_vertices = []
#     hole_vertices = []
    
#     # Iterate over the Panel structs
#     data = []
#     structs = root.findall('.//struct')
#     # for struct in structs:
#     #     if struct.attrib.get("identifier") == "Panel":
    
#     polyshapes = root.findall('.//polyshape')
#     polyshape_ids = []
#     for polyshape in polyshapes:
#         polyshape_ids.append(polyshape.attrib.get("identifier"))
        
#     # tables = root.findall('.//table')
#     # for table in tables:
#     #     if table.attrib.get("identifier") == "PressureRods":
#     pressurerods_rows = get_table_rows(root, "PressureRods")
#     standoffs_rows = get_table_rows(root, "Standoffs")
        
#     # Gather data into usable form
#     if "pBoards" in polyshape_ids:
#         pBoards = root.findall('.//polyshape[@identifier="pBoards"]')
#         regions = pBoards[0].findall('.//region')
#         pBoards_regions = []
#         for region in regions:
#             # Extract vertices for the region
#             vertices = []
#             for vertex in region.findall('.//vertex'):
#                 x, y = map(float, vertex.text.split('|'))
#                 vertices.append([x, y])
#             pBoards_regions.append(vertices)

#         holes = pBoards[0].findall('.//hole')
#         pBoards_holes = []
#         for hole in holes:
#             # Extract vertices for the hole
#             vertices = []
#             for vertex in hole.findall('.//vertex'):
#                 x, y = map(float, vertex.text.split('|'))
#                 vertices.append([x, y])
#             pBoards_holes.append(vertices)
            
#         # Create pandas DataFrames
#         if len(pBoards_regions) == 0:
#             pBoards_regions_df = None
#         else:
#             pBoards_regions_df = pd.DataFrame(pBoards_regions, columns=['Vertices'])
#         if len(pBoards_holes) == 0:
#             pBoards_holes_df = None
#         else:
#             pBoards_holes_df = pd.DataFrame(pBoards_holes, columns=['Vertices'])
#     else:
#         # pBoards_regions = None
#         pBoards_regions_df = None
#         pBoards_holes_df = None
        
        
#     return pBoards_regions_df, pBoards_holes_df
            
            
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
    # pBoards_regions_df, pBoards_holes_df = parse_input_xml(initialdir)
    root = get_XML_tree(initialdir)
    pBoards, numRegions, numHoles = get_pBoards(root)
    plot_poly_list_w_holes(pBoards)
    print(f"Regions:\t{numRegions}\nHoles:\t{numHoles}")