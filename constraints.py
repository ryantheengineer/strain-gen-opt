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
from shapely.affinity import rotate, translate
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
            break
    
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
        region_polys = None
        numRegions = None
        numHoles = None
        
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


def get_region_shifts(region_polys):
    polys_copy = region_polys.copy()
    refpoly = polys_copy.pop(0)
    angles = [0, 45, 90, 135, 180, 225, 270, 315]
    
    xshifts = [0.0]
    yshifts = [0.0]
    rotations = [0.0]
    
    for poly in polys_copy:
        for angle in angles:
            origin = poly.centroid
            rotpoly = rotate(poly, angle=angle, origin=origin)
            shifted, xshift, yshift = check_shifted(refpoly, rotpoly)
            if shifted==True:
                xshifts.append(xshift)
                yshifts.append(yshift)
                rotations.append(angle)
                break
            elif shifted==False and angle==angles[-1]:
                xshifts.append(xshift)
                yshifts.append(xshift)
                rotations.append(None)
    
    return xshifts, yshifts, rotations
                
            
def check_shifted(refpoly, poly):
    # Only checks whether poly is identical to refpoly, just translated
    refcoords = refpoly.exterior.coords
    coords = poly.exterior.coords
    xdiffs = []
    ydiffs = []
    sensitivity = 0.0001
    for i in range(len(coords)):
        xdiff = coords[i][0] - refcoords[i][0]
        ydiff = coords[i][1] - refcoords[i][1]
        xdiffs.append(xdiff)
        ydiffs.append(ydiff)
    if len(set(xdiffs))==1 and len(set(ydiffs))==1:
        shifted = True
        xshift = xdiffs[0]
        yshift = ydiffs[0]
    elif np.std(list(set(xdiffs)))<sensitivity and np.std(list(set(ydiffs)))<sensitivity:
        shifted = True
        xshift = np.median(xdiffs)
        yshift = np.median(ydiffs)
    else:
        shifted = False
        xshift = None
        yshift = None
        
    return shifted, xshift, yshift
    


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
    
    xshifts, yshifts, rotations = get_region_shifts(pBoards)
    df = pd.DataFrame({"xshift":xshifts, "yshift":yshifts, "rotation":rotations})
    print(df)