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
from shapely import buffer
from shapely.ops import unary_union
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
import runFEA
import os

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
    return treeroot, inputfile

########## Obtain relevant Fixture info ##########
# Panel, Plates
def get_fixture_geometry(root, identifier):
    polyshapes = root.findall('.//polyshape')
    fixture = None
    region_polys = None
    hole_polys = None
    for polyshape in polyshapes:
        if polyshape.attrib.get("identifier") == identifier:
            fixture = polyshape
            break
    
    def verts_to_poly(polys, elem):
        vertices_elem = elem.findall('.//vertex')
        vertices = [vertex.text.split("|") for vertex in vertices_elem]
        polys.append(Polygon(vertices))
        
        
    if fixture:
        numRegions = int(fixture.find('.//numRegions').text)
        numHoles = int(fixture.find('.//numHoles').text)
        if numRegions:
            # Create the Polygons of each region
            regions = fixture.findall('.//region')
            region_polys = []
            for region in regions:
                verts_to_poly(region_polys, region)
                # vertices_elem = region.findall('.//vertex')
                # vertices = [vertex.text.split("|") for vertex in vertices_elem]
                # region_polys.append(Polygon(vertices))
        if numHoles:
            # Create the Polygons of each hole
            holes = fixture.findall('.//hole')
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
                
    dfshifts = pd.DataFrame({"xshift":xshifts, "yshift":yshifts, "rotation":rotations})
    
    return dfshifts
                
            
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
    

# FEA functions
def runFEA_valid_circles(valid_circles, df_PressureRods, root, inputfile):
    # Ensure correct type is used for integer columns
    df_PressureRods['unimplemented1'] = df_PressureRods['unimplemented1'].astype(int)
    df_PressureRods['unimplemented2'] = df_PressureRods['unimplemented2'].astype(int)
    df_PressureRods['unimplemented3'] = df_PressureRods['unimplemented3'].astype(int)
    df_PressureRods['unimplemented4'] = df_PressureRods['unimplemented4'].astype(int)
    df_PressureRods['unimplemented5'] = df_PressureRods['unimplemented5'].astype(str)
    df_PressureRods['unimplemented5'] = df_PressureRods['unimplemented5'].apply(str.lower)
    df_PressureRods['unimplemented6'] = df_PressureRods['unimplemented6'].astype(int)
    df_PressureRods['unimplemented7'] = df_PressureRods['unimplemented7'].astype(int)
    df_PressureRods['unimplemented8'] = df_PressureRods['unimplemented8'].astype(int)
    df_PressureRods['unimplemented10'] = df_PressureRods['unimplemented10'].astype(int)
    df_PressureRods['unimplemented11'] = df_PressureRods['unimplemented11'].astype(int)
    df_PressureRods['type'] = df_PressureRods['type'].astype(str)
    df_PressureRods['unimplemented12'] = df_PressureRods['unimplemented12'].astype(str)
    df_PressureRods['unimplemented12'] = df_PressureRods['unimplemented12'].apply(str.lower)
    df_PressureRods['unimplemented13'] = df_PressureRods['unimplemented13'].astype(int)
    
    df = df_PressureRods.loc[0,:]
    basic_vals = {}
    for col in df.index:
        basic_vals[col] = df[col]
    # new_vals = [[] for col in df.index]
    
    xnew = []
    ynew = []
    for valid_circle in valid_circles:
        xnew.append(valid_circle.centroid.x)
        ynew.append(valid_circle.centroid.y)
    
    newdata = {}
    for col in df.index:
        if col == "x":
            newdata[col] = xnew
        elif col == "y":
            newdata[col] = ynew
        else:
            newdata[col] = [basic_vals[col] for valid_circle in valid_circles]
            
    df_PressureRods_update = pd.DataFrame(newdata)
    
    # Get the rows of df_PressureRods_update as a list of strings
    vals = df_PressureRods_update.to_string(header=False,
                  index=False,
                  index_names=False).split('\n')
    for j,row in enumerate(vals):
        splitrow = row.split()
        newrow = []
        for i,ele in enumerate(splitrow):
            if i==13:
                newrow.append(' '.join(splitrow[13:15]))
            elif i==14:
                continue
            elif i==18:
                newrow.append(''.join(splitrow[18:]))
            elif i==19:
                continue
            elif i==20:
                continue
            else:
                newrow.append(ele)
        newrow = '|'.join(newrow)
        vals[j] = newrow
    # vals = ['|'.join(ele.split()) for ele in x
    
    # Delete the XML rows that need to be replaced
    rows = root.find('.//table[@identifier="PressureRods"].//rows')
    del rows[:]
    
    rowlist = [ET.SubElement(rows,'row') for val in vals]

    for i,val in enumerate(vals):
        rowlist[i].text = val
    # for i,val in enumerate(vals):
    #     rowlist[i].text = val
        
    # for val in vals:
    #     row = ET.SubElement(rows, 'row')
    #     row.text = val
    #     # # val_elem = ET.SubElement(rows, 'row')
    #     # val_elem = ET.Element('row')
    #     # val_elem.text = val
    #     # rows.append(val_elem)
    
    # tables = root.findall('.//table')
    # for table in tables:
    #     if table.attrib.get("identifier") == "PressureRods":            
    #         for row in table.findall('.//row'):
    #             table.remove(row)
            
    #         for val in vals:
    #             table.insert(val)
    
    # Write the new element tree
    tree = ET.ElementTree(root)
    ET.indent(tree, '  ')
    # new_inputfile = 
    # path = pathlib.Path("/path/to/file.txt")
    new_filename = "FEA_random_pressure_rods.xml"
    path, filename = os.path.split(inputfile)
    new_path = path + "/" + new_filename
    
    # new_path = os.inputfile.join(os.inputfile.dirname(inputfile), new_filename)
    # inputfile.rename(new_filename)
    tree.write(new_path)
    # with open('FEA_random_pressure_rods.xml', 'w') as f:
    #     root.write(f, encoding='unicode')
    
    FEApath = runFEA.loadFEApath('FEApath.pk')
    runFEA.runFEA(FEApath, new_path)
    
    dfmesh = runFEA.resultsToDataframe(inputfile)
    
    strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max = runFEA.getFitness(dfmesh)
    
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max
    

# Utility functions
def plot_poly_list_w_holes(poly_list, fig, ax, color, linestyle, label):
    # fig, ax = plt.subplots(figsize=(10,8),dpi=300)
    # ax.set_aspect('equal')
    first = True
    for poly in poly_list:
        # Check if poly is a Polygon or MultiPolygon
        if poly.geom_type == "Polygon":
            xe, ye = poly.exterior.xy
            if first == True:
                ax.plot(xe, ye, color=color, label=label, linestyle=linestyle, linewidth=0.5)
                first = False
            else:
                ax.plot(xe, ye, color=color, linestyle=linestyle, linewidth=0.5)
            for inner in poly.interiors:
                xi, yi = zip(*inner.coords[:])
                ax.plot(xi, yi, color=color, linestyle=linestyle, linewidth=0.5)
        elif poly.geom_type == "MultiPolygon":
            for geom in poly.geoms:
                xe, ye = geom.exterior.xy
                if first == True:
                    ax.plot(xe, ye, color=color, label=label, linestyle=linestyle, linewidth=0.5)
                    first = False
                else:
                    ax.plot(xe, ye, color=color, linestyle=linestyle, linewidth=0.5)
                for inner in geom.interiors:
                    xi, yi = zip(*inner.coords[:])
                    ax.plot(xi, yi, color=color, linestyle=linestyle, linewidth=0.5)            
        else:
            raise IOError("Shape is not a polygon")
            
def plot_probes_guidepins(df, fig, ax, linestyle, identifier):
    first = True
    for row in range(len(df)):
        x = df.loc[row,"x"]
        y = df.loc[row,"y"]
        radius = df.loc[row,"diameter"]/2
        circle = Point(x, y).buffer(radius)
        xe, ye = circle.exterior.xy
        color = tuple(df.loc[row,"color"])
        if first == True:
            ax.plot(xe, ye, color=color, label=identifier, linestyle=linestyle, linewidth=0.5)
            first = False
        else:
            ax.plot(xe, ye, color=color, linestyle=linestyle, linewidth=0.5)
            
def plot_pressurerods_standoffs(df, fig, ax, linestyle, identifier):
    if identifier == "PressureRods":
        radius = 0.0625/2
    elif identifier == "Standoffs":
        radius = 0.15/2
    else:
        raise ValueError("Invalid identifier fed to plot_pressurerods_standoffs")
        
    first = True
    for row in range(len(df)):
        x = df.loc[row,"x"]
        y = df.loc[row,"y"]
        circle = Point(x, y).buffer(radius)
        xe, ye = circle.exterior.xy
        color = tuple(df.loc[row,"color"])
        if first == True:
            ax.plot(xe, ye, color=color, label=identifier, linestyle=linestyle, linewidth=0.5)
            first = False
        else:
            ax.plot(xe, ye, color=color, linestyle=linestyle, linewidth=0.5)
    
    

def get_region_list_centroids(region_polys):
    region_centroids = [poly.centroid for poly in region_polys]
    return region_centroids



# def get_table_rows(root, identifier):
#     tables = root.findall('.//table')
#     rows = None
#     for table in tables:
#         if table.attrib.get("identifier") == identifier:
#             # Extract rows
#             rows = []
#             for row in table.findall('.//row'):
#                 rows.append(row.text.split('|'))
#     return rows

def get_point_geometry(root, identifier):
    tables = root.findall('.//table')
    column_names = None
    types = None
    rows = None
    df_output = None
    for table in tables:
        if table.attrib.get("identifier") == identifier:
            column_names = []
            for col in table.findall('.//col'):
                column_names.append(col.text)
            
            types = []
            for typ in table.findall('.//type'):
                types.append(typ.text)
            
            rows = []
            for row in table.findall('.//row'):
                rows.append(row.text.split('|'))
            
            if len(column_names) != len(types) or len(column_names) != len(rows[0]):
                raise ValueError(f"Mismatch between number of columns, types, or row elements for {identifier}")
            
            columns = []
            for i in range(len(column_names)):
                columns.append([])
                for j in range(len(rows)):
                    columns[i].append(interpret_type_description(types[i],rows[j][i]))
                    
            output_dict = {column_names[i]:columns[i] for i in range(len(column_names))}
            df_output = pd.DataFrame(output_dict)
             
    return df_output
            

def interpret_type_description(type_str, elem):
    if type_str == "double":
        if elem[0]=='[' and elem[-1]==']':
            res = elem.strip('][').split(',')
            reslist = [float(x) for x in res]
            return reslist
        else:
            return float(elem)
    elif type_str == "string":
        return elem
    elif type_str == "logical":
        if elem.lower() == "false":
            return False
        else:
            return True
    

def place_circle(x,y,radius):
    circle = Point(x, y).buffer(radius)
    return circle

def random_perturb(poly, maxmag):
    while True:
        xp = random.uniform(-maxmag, maxmag)
        yp = random.uniform(-maxmag, maxmag)
        magnitude = np.sqrt(xp**2 + yp**2)
        if magnitude > maxmag:
            continue
        else:
            break
    
    # Adjust polygon by perturbation
    xe, ye = poly.exterior.xy
    xe = list(xe)
    ye = list(ye)
    xe_perturb = [x+xp for x in xe]
    ye_perturb = [y+yp for y in ye]
    xpyp = zip(xe_perturb, ye_perturb)
    xpyp = tuple(xpyp)
    poly_perturb = Polygon(xpyp)
    return poly_perturb

def random_valid_perturb(poly, outer_poly, maxmag):
    count = 0
    lim = 100
    if outer_poly.area > poly.area:
        while count < lim:
            # while True:
            #     xp = random.uniform(-maxmag, maxmag)
            #     yp = random.uniform(-maxmag, maxmag)
            #     magnitude = np.sqrt(xp**2 + yp**2)
            #     if magnitude > maxmag:
            #         continue
            #     else:
            #         break
                
            # Approach: Choose a random magnitude, then choose the first
            # component randomly that would give that magnitude, then solve
            # for the remaining component
            mag = random.uniform(0, maxmag)
            xp = random.uniform(-mag, mag)
            yp = np.sqrt(mag**2 - xp**2)
            yp = random.choice([yp, -yp])
                
            # Adjust polygon by perturbation
            xe, ye = poly.exterior.xy
            xe = list(xe)
            ye = list(ye)
            xe_perturb = [x+xp for x in xe]
            ye_perturb = [y+yp for y in ye]
            xpyp = zip(xe_perturb, ye_perturb)
            xpyp = tuple(xpyp)
            poly_perturb = Polygon(xpyp)
            
            if outer_poly.contains(poly_perturb):
                break
            else:
                count += 1
        return poly_perturb
    else:
        return poly


# %% Load data and read in the XML definition
if __name__ == "__main__":
    # Load previously chosen FEA path here
    filename = 'FEApath.pk'
    with open(filename, 'rb') as fi:
        FEApath = pickle.load(fi)
    
    initialdir = str(pathlib.Path(FEApath).parent) + "Examples"
    # pBoards_regions_df, pBoards_holes_df = parse_input_xml(initialdir)
    root, inputfile = get_XML_tree(initialdir)
    
    # Panel
    pBoards, _, _ = get_fixture_geometry(root, "pBoards")
    pOutline, _, _ = get_fixture_geometry(root, "pOutline")
    pShape, _, _ = get_fixture_geometry(root, "pShape")
    pComponentsTop, _, _ = get_fixture_geometry(root, "pComponentsTop")
    pComponentsBot, _, _ = get_fixture_geometry(root, "pComponentsBot")
    # fig1, ax1 = plt.subplots(figsize=(10,8),dpi=300)
    # ax1.set_aspect('equal')
    # labels = []
    # line = "-"
    # if pBoards:
    #     plot_poly_list_w_holes(pBoards, fig1, ax1, "blue", line, "pBoards")
    # if pOutline:
    #     plot_poly_list_w_holes(pOutline, fig1, ax1, "green", line, "pOutline")
    # if pShape:
    #     plot_poly_list_w_holes(pShape, fig1, ax1, "black", line, "pShape")
    # if pComponentsTop:
    #     plot_poly_list_w_holes(pComponentsTop, fig1, ax1, "purple", line, "pComponentsTop")
    # if pComponentsBot:
    #     plot_poly_list_w_holes(pComponentsBot, fig1, ax1, "orange", line, "pComponentsBot")
    # ax1.legend()
    
    # Plates
    Pressure, _, _ = get_fixture_geometry(root, "Pressure")
    I_Plate, _, _ = get_fixture_geometry(root, "I_Plate")
    Stripper, _, _ = get_fixture_geometry(root, "Stripper")
    Probe, _, _ = get_fixture_geometry(root, "Probe")
    Countersink, _, _ = get_fixture_geometry(root, "Countersink")
    df_Probes = get_point_geometry(root, "Probes")
    df_GuidePins = get_point_geometry(root, "GuidePins")
    df_PressureRods = get_point_geometry(root, "PressureRods")
    df_Standoffs = get_point_geometry(root, "Standoffs")
    
    # #########################################################################
    # ############# Plot Plates and Features ##############################
    # #########################################################################
    # # Top Plate/Pressure Plate (doesn't hold pressure rods if I_plate exists)
    # if I_Plate:
    #     fig1, ax1 = plt.subplots(figsize=(10,8),dpi=500)
    #     line="-"
    #     plot_poly_list_w_holes(I_Plate, fig1, ax1, "blue", line, "I_Plate")
    #     plot_poly_list_w_holes(Pressure, fig1, ax1, "red", line, "Pressure")
    #     if pBoards:
    #         plot_poly_list_w_holes(pBoards, fig1, ax1, "black", line, "pBoards")
    #     if pComponentsTop:
    #         plot_poly_list_w_holes(pComponentsTop, fig1, ax1, "purple", line, "pComponentsTop")
    #     if df_PressureRods is not None:
    #         plot_pressurerods_standoffs(df_PressureRods, fig1, ax1, line, "PressureRods")
    #     ax1.set_title("Pressure and Intermediate Plates")
    #     ax1.legend()
    # else:
    #     print("\nExample has no intermediate plate\n")
    #     fig1, ax1 = plt.subplots(figsize=(10,8),dpi=500)
    #     line="-"
    #     plot_poly_list_w_holes(Pressure, fig1, ax1, "blue", line, "Pressure")
    #     if pBoards:
    #         plot_poly_list_w_holes(pBoards, fig1, ax1, "black", line, "pBoards")
    #     if pComponentsTop:
    #         plot_poly_list_w_holes(pComponentsTop, fig1, ax1, "purple", line, "pComponentsTop")
    #     if df_PressureRods is not None:
    #         plot_pressurerods_standoffs(df_PressureRods, fig1, ax1, line, "PressureRods")
    #     ax1.set_title("Pressure Plate")
    #     ax1.legend()
        
    # # Probe Protector Plate (Stripper Plate)
    # if Stripper:
    #     fig2, ax2 = plt.subplots(figsize=(10,8),dpi=500)
    #     line="-"
    #     plot_poly_list_w_holes(Stripper, fig2, ax2, "blue", line, "Stripper")
    #     plot_poly_list_w_holes(Probe, fig2, ax2, "red", line, "Probe")
    #     if pBoards:
    #         plot_poly_list_w_holes(pBoards, fig2, ax2, "black", line, "pBoards")
    #     if pComponentsBot:
    #         plot_poly_list_w_holes(pComponentsBot, fig2, ax2, "purple", line, "pComponentsBot")
    #     if df_Probes is not None:
    #         plot_probes_guidepins(df_Probes, fig2, ax2, line, "Probes")
    #     if df_GuidePins is not None:
    #         plot_probes_guidepins(df_GuidePins, fig2, ax2, line, "GuidePins")
    #     if df_Standoffs is not None:
    #         plot_pressurerods_standoffs(df_Standoffs, fig2, ax2, line, "Standoffs")
    #     ax2.set_title("Probe and Probe Protector (Stripper)Plates")
    #     ax2.legend()
    # else:
    #     print("\nExample has no probe protector plate (stripper plate)\n")
    #     fig2, ax2 = plt.subplots(figsize=(10,8),dpi=500)
    #     line="-"
    #     plot_poly_list_w_holes(Probe, fig2, ax2, "red", line, "Probe")
    #     if pBoards:
    #         plot_poly_list_w_holes(pBoards, fig2, ax2, "black", line, "pBoards")
    #     if pComponentsBot:
    #         plot_poly_list_w_holes(pComponentsBot, fig2, ax2, "purple", line, "pComponentsBot")
    #     if df_Probes is not None:
    #         plot_probes_guidepins(df_Probes, fig2, ax2, line, "Probes")
    #     if df_GuidePins is not None:
    #         plot_probes_guidepins(df_GuidePins, fig2, ax2, line, "GuidePins")
    #     if df_Standoffs is not None:
    #         plot_pressurerods_standoffs(df_Standoffs, fig2, ax2, line, "Standoffs")
    #     ax2.set_title("Probe and Probe Protector (Stripper)Plates")
    #     ax2.legend()
        
    # # plt.show()
    
    # %% Place pressure rods
    ### Randomly place pressure rod circles on top side of board ###
    # Add buffers around components and edge of boards
    buffer_dist = 0.025
    pBoards_dilated = []
    for board in pBoards:
        UUT_poly_ext = Polygon(board.exterior.coords)
        UUT_poly_dilated_ext = buffer(UUT_poly_ext,-buffer_dist)
        pBoards_dilated.append(UUT_poly_dilated_ext)
    
    topcomponents_dilated = []
    for inner in pComponentsTop:
        topcomponents_dilated.append(buffer(inner, buffer_dist))
    
    topcomponents_dilated = unary_union(topcomponents_dilated)
    if topcomponents_dilated.geom_type == "MultiPolygon":
        topcomponents_dilated = list(topcomponents_dilated.geoms)
    
    # Subtract dilated interiors from the exterior polygon or multipolygon
    for i,board in enumerate(pBoards_dilated):
        for inner in topcomponents_dilated:
            if inner.intersects(board):
                pBoards_dilated[i] = pBoards_dilated[i].difference(inner)
                
    # Plot the original board with components and then the dilated version
    fig3, ax3 = plt.subplots(figsize=(10,8),dpi=500)
    ax3.set_aspect('equal')
    line="-"
    plot_poly_list_w_holes(pBoards, fig3, ax3, "black", line, "pBoards")
    if pComponentsTop:
        plot_poly_list_w_holes(pComponentsTop, fig3, ax3, "purple", line, "pComponentsTop")
    plot_poly_list_w_holes(pBoards_dilated, fig3, ax3, "blue", line, "pBoards_dilated")
    
    # Place circles on grid within UUT
    valid_circles = []
    first = True
    for board in pBoards_dilated:
        # FIXME: Here this should iterate on the polygons that define the 
        # available space, not the outline of the PCBs. Each PCB could have
        # several possible areas that could have a new grid of circles in it.
        # Additionally, every area that could have a circle, should have a
        # circle to make sure the optimization fully explores the design space.
        # The grid resolution concept may need to be improved on to make this
        # happen.
        
        xmin, ymin, xmax, ymax = board.bounds
        circles = []
        radius = 0.0625/2
        resolution = radius*2 + 0.2
        perturb = 0.05/2
        for x in np.arange(xmin+radius, xmax, resolution):
            for y in np.arange(ymin+radius, ymax, resolution):
                c = place_circle(x, y, radius)
                circles.append(c)
        
        # Randomly perturb each circle
        for i,circle in enumerate(circles):
            circles[i] = random_valid_perturb(circle, board, perturb)
            # circles[i] = random_perturb(circle, perturb)
            # NOTE: Before accepting randomly perturbed circle, make sure it is 
            # within spec. Force rods can't be too close together.
        
        # Validate that each circle falls inside shape
        valid_circles.extend(filter(board.contains, circles)) # FIX: This is an old version that doesn't check for multipolygons. Some circle overlaps exist
        
        # Plot circles
        for circle in valid_circles:
            xe, ye = circle.exterior.xy
            if first==True:
                ax3.plot(xe, ye, color="red", label="random pressure rods")
                first = False
            else:
                ax3.plot(xe, ye, color="red")
                
    
    ax3.set_title("Dilated boundaries around board and components")
    ax3.legend()
        
    plt.show()
    
    # labels = []
    # line = ":"
    # if Pressure:
    #     plot_poly_list_w_holes(Pressure, fig1, ax1, "blue", line, "Pressure")
    # if I_Plate:
    #     plot_poly_list_w_holes(I_Plate, fig1, ax1, "green", line, "I_Plate")
    # if Stripper:
    #     plot_poly_list_w_holes(Stripper, fig1, ax1, "black", line, "Stripper")
    # if Probe:
    #     plot_poly_list_w_holes(Probe, fig1, ax1, "purple", line, "Probe")
    # if Countersink:
    #     plot_poly_list_w_holes(Countersink, fig1, ax1, "orange", line, "Countersink")
    
    
    # line = "-"
    # if df_Probes is not None:
    #     plot_probes_guidepins(df_Probes, fig1, ax1, line, "Probes")
    # if df_GuidePins is not None:
    #     plot_probes_guidepins(df_GuidePins, fig1, ax1, line, "GuidePins")
    # if df_PressureRods is not None:
    #     plot_pressurerods_standoffs(df_PressureRods, fig1, ax1, line, "PressureRods")
    # if df_Standoffs is not None:
    #     plot_pressurerods_standoffs(df_Standoffs, fig1, ax1, line, "Standoffs")
    
    # ax1.legend()
    
    # plot_poly_list_w_holes(pBoards)
    # print(f"Regions:\t{numRegions}\nHoles:\t{numHoles}")
    
    
    
    
    
    dfshifts = get_region_shifts(pBoards)
    
    # %% Try to run FEA using the randomly placed pressure rods
    # Run FEA with valid_circles defining the pressure rod points
    strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max = runFEA_valid_circles(valid_circles, df_PressureRods, root, inputfile)
    