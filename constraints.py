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
# from shapely import buffer
from shapely.ops import unary_union
# import xmltodict
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
import time
import multiprocessing
import queue

# %% XML Functions
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


# %% Fixture Reading Functions
# Super function for getting fixture data
def get_constraint_geometry():
    # Load previously chosen FEA path here
    filename = 'FEApath.pk'
    with open(filename, 'rb') as fi:
        FEApath = pickle.load(fi)
    
    initialdir = str(pathlib.Path(FEApath).parent) + "Examples"
    root, inputfile = get_XML_tree(initialdir)
    
    print("--- Reading in UUT geometry ---")
    # Panel
    pBoards, _, _ = get_fixture_geometry(root, "pBoards")
    pOutline, _, _ = get_fixture_geometry(root, "pOutline")
    pShape, _, _ = get_fixture_geometry(root, "pShape")
    pComponentsTop, _, _ = get_fixture_geometry(root, "pComponentsTop")
    pComponentsBot, _, _ = get_fixture_geometry(root, "pComponentsBot")
    
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
    
    results = (root, inputfile, pBoards, pOutline, pShape, pComponentsTop,
               pComponentsBot, Pressure, I_Plate, Stripper, Probe, Countersink,
               df_Probes, df_GuidePins, df_PressureRods, df_Standoffs)
    return results


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
    

#%% Circle placing functions
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

class PressureRod():
    def __init__(self, x, y, rod_type, on):
        self.x = x
        self.y = y
        self.rod_type = rod_type
        self.on = on
        self.ctc = 0.375
        self.tip_from_UUT_edge = 0.0
        self.update_pressure_rod(self.x, self.y, self.rod_type, self.on)        
        
    def select_rod_type(self, rod_type):
        if rod_type == 'Press-Fit Tapered':
            dtop = 0.25
            dtip = 0.09
        elif rod_type == 'Press-Fit Flat':
            dtop = 0.25
            dtip = 0.19
        elif rod_type == '3.325" Tapered':
            dtop = 0.315
            dtip = 0.10
        elif rod_type == '3.325" Flat':
            dtop = 0.315
            dtip = 0.315
        else:
            raise ValueError("ERROR: Invalid rod_type string")
        self.rod_type = rod_type
        self.dtop = dtop
        self.dtip = dtip
        
    def update_pressure_rod(self, new_x, new_y, new_rod_type, new_on):
        self.x = new_x
        self.y = new_y
        self.center = Point(new_x, new_y)
        self.select_rod_type(new_rod_type)
        self.rtip = self.dtip/2.0
        self.rtop = self.dtop/2.0
        self.tip = place_circle(new_x, new_y, self.rtip)
        self.top = place_circle(new_x, new_y, self.rtop)
        self.tip_from_components = self.rtip + 0.035
        self.top_from_components = self.rtop + 0.0625
        self.tip_from_top_probes = self.rtip + 0.125
        self.top_from_top_probes = self.rtop + 0.04
        self.tip_component_buffer = place_circle(new_x, new_y, self.rtip + self.tip_from_components)
        self.top_component_buffer = place_circle(new_x, new_y, self.rtop + self.top_from_components)
        self.tip_UUT_buffer = place_circle(new_x, new_y, self.rtip + self.tip_from_UUT_edge)
        self.tip_from_top_probe_buffer = place_circle(new_x, new_y, self.rtip + self.tip_from_top_probes)
        self.top_from_top_probe_buffer = place_circle(new_x, new_y, self.rtop + self.top_from_top_probes)
        self.prod_to_prod_buffer = place_circle(new_x, new_y, self.rtop + 0.025)
        self.on = new_on
        
        
        
def grid_nprods(pBoards, pComponentsTop):
    ### Randomly place pressure rod circles on top side of board ###
    # Add buffers around components and edge of boards
    # buffer_dist = 0.025
    pBoards_diff = []
    for board in pBoards:
        UUT_poly_ext = Polygon(board.exterior.coords)
        # UUT_poly_dilated_ext = buffer(UUT_poly_ext,-buffer_dist)
        pBoards_diff.append(UUT_poly_ext)
    
    topcomponents = []
    for inner in pComponentsTop:
        topcomponents.append(inner)
    
    topcomponents = unary_union(topcomponents)
    if topcomponents.geom_type == "MultiPolygon":
        topcomponents = list(topcomponents.geoms)
    
    # Subtract interiors from the exterior polygon or multipolygon
    for i,board in enumerate(pBoards_diff):
        for inner in topcomponents:
            if inner.intersects(board):
                pBoards_diff[i] = pBoards_diff[i].difference(inner)
    
    # Place circles on grid within UUT
    valid_prods_small = []
    valid_prods_large = []
    for board in pBoards_diff:
        # FIXME: Here this should iterate on the polygons that define the 
        # available space, not the outline of the PCBs. Each PCB could have
        # several possible areas that could have a new grid of circles in it.
        # Additionally, every area that could have a circle, should have a
        # circle to make sure the optimization fully explores the design space.
        # The grid resolution concept may need to be improved on to make this
        # happen.
        
        xmin, ymin, xmax, ymax = board.bounds
        prods_sm = []
        prods_lg = []
        radius_sm = 0.09/2
        radius_lg = 0.315/2
        resolution = 0.375
        # perturb = 0.05/2
        for x in np.arange(xmin+radius_sm, xmax, resolution):
            for y in np.arange(ymin+radius_sm, ymax, resolution):
                c_sm = place_circle(x, y, radius_sm)
                prods_sm.append(c_sm)
        for x in np.arange(xmin+radius_lg, xmax, resolution):
            for y in np.arange(ymin+radius_lg, ymax, resolution):
                c_lg = place_circle(x, y, radius_lg)
                prods_lg.append(c_lg)
                
        
        
        # # Randomly perturb each circle
        # for i,circle in enumerate(circles):
        #     circles[i] = random_valid_perturb(circle, board, perturb)
        #     # circles[i] = random_perturb(circle, perturb)
        #     # NOTE: Before accepting randomly perturbed circle, make sure it is 
        #     # within spec. Force rods can't be too close together.
        
        # Validate that each circle falls inside shape
        valid_prods_small.extend(filter(board.contains, prods_sm)) # FIXME: This is an old version that doesn't check for multipolygons. Some circle overlaps exist
        valid_prods_large.extend(filter(board.contains, prods_lg)) # FIXME: This is an old version that doesn't check for multipolygons. Some circle overlaps exist
        
    nprods_small = len(valid_prods_small)
    nprods_large = len(valid_prods_large)
    
    return nprods_small, nprods_large, pBoards_diff


def create_chromosome(nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff, all_on=False, on_prob=0.5, rod_type="All"):
    # Initialize random chromosome, where the first ncircles entries are
    # the x coordinates, the next ncircles entries are y coordinates, then
    # radii, and on/off binary values
    
    # Get pBoards as a MultiPolygon so pressure rods can be placed
    # anywhere within the UUT grid
    pBoards_multi = MultiPolygon(pBoards_diff)
    
    df_Probes_top = df_Probes[df_Probes["side"]==1]
    df_Probes_top.reset_index(inplace=True)
    top_probes = []
    for i,row in df_Probes_top.iterrows():
        top_probes.append(place_circle(row.x, row.y, row.diameter/2))
        
    top_probes = unary_union(top_probes)
    if top_probes.geom_type != "MultiPolygon":
        top_probes = MultiPolygon(top_probes)
    
    topcomponents = []
    for inner in pComponentsTop:
        topcomponents.append(inner)
    
    topcomponents = unary_union(topcomponents)
    if topcomponents.geom_type != "MultiPolygon":
        topcomponents = MultiPolygon(topcomponents)
    
    xmin, ymin, xmax, ymax = pBoards_multi.bounds
    rod_types = ['Press-Fit Tapered',
                 'Press-Fit Flat',
                 '3.325" Tapered',
                 '3.325" Flat']
    
    chromosome_x = []
    chromosome_y = []
    chromosome_rod_type = []
    chromosome_on = []
    prods_chosen = []
    tries = 10
    
    while len(chromosome_x) < nprods:
        valid = True
        x = random.uniform(xmin,xmax)
        y = random.uniform(ymin,ymax)
        if rod_type not in rod_types:
            rod_type_i = random.randint(0,3)
        else:
            rod_type_i = rod_types.index(rod_type)
        if all_on:
            on = 1
        else:
            # on_prob = 0.5   # Percentage chance of a pressure rod being on initially
            on_chance = random.uniform(0,1)
            if on_prob >= on_chance:
                on = 1
            else:
                on = 0
            # on = random.randint(0,1)
        prod = PressureRod(x,y,rod_types[rod_type_i],on)
        intersects_topcomponents = False
        intersects_topprobes = False
        centroids = []
        perturbing = False
        
        # Make sure pressure rod is within the UUT and make sure it doesn't intersect any components, using the appropriate buffer sizes
        if not prod.center.intersects(pBoards_multi):
            continue
        # if not prod.center.within(pBoards_multi):
        #     continue
        
        ### TRYING MOVING AWAY FROM INTERSECTION VIOLATIONS TO SALVAGE DESIGN ###
        # Gather status of currently placed prod
        if prod.tip_component_buffer.intersects(topcomponents):
            intersects_topcomponents = True
            intersection_topcomponents = prod.tip_component_buffer.intersection(topcomponents)
            
        if prod.tip_from_top_probe_buffer.intersects(top_probes):
            intersects_topprobes = True
            intersection_topprobes = prod.tip_from_top_probe_buffer.intersection(top_probes)
        
        # Parse the status of intersections
        if not intersects_topcomponents and not intersects_topprobes:
            pass
        elif intersects_topcomponents and not intersects_topprobes:
            perturbing = True
            if intersection_topcomponents.geom_type == "Polygon":
                centroids.append(intersection_topcomponents.centroid)
            else:
                for poly in intersection_topcomponents.geoms:
                    centroids.append(poly.centroid)
        elif not intersects_topcomponents and intersects_topprobes:
            perturbing = True
            if intersection_topprobes.geom_type == "Polygon":
                centroids.append(intersection_topprobes.centroid)
            else:
                for poly in intersection_topprobes.geoms:
                    centroids.append(poly.centroid)                    
        else:
            perturbing = True
            # both top components and top probes are intersected by the current prod
            if intersection_topcomponents.geom_type == "Polygon":
                centroids.append(intersection_topcomponents.centroid)
            else:
                for poly in intersection_topcomponents.geoms:
                    centroids.append(poly.centroid)
                    
            if intersection_topprobes.geom_type == "Polygon":
                centroids.append(intersection_topprobes.centroid)
            else:
                for poly in intersection_topprobes.geoms:
                    centroids.append(poly.centroid)
        
        if perturbing:                    
            centroids_x = [centroid.x for centroid in centroids]
            centroids_y = [centroid.y for centroid in centroids]
            
            avg_x = np.mean(centroids_x)
            avg_y = np.mean(centroids_y)
            
            diff_x = prod.x - avg_x
            diff_y = prod.y - avg_y
            
            mag_diff = np.sqrt(diff_x**2 + diff_y**2)
            
            unit_x = diff_x / mag_diff
            unit_y = diff_y / mag_diff
            
            # print("")
            # print("-"*60)
            # print("PERTURBING PROD AWAY FROM INTERSECTION")
            
            for i in range(tries):
                valid = True
                prodxmin, prodymin, prodxmax, prodymax = prod.tip_component_buffer.bounds
                stepsize = (prodxmax - prodxmin)/tries
                new_x = prod.x + unit_x*stepsize
                new_y = prod.y + unit_y*stepsize
                
                prod.update_pressure_rod(new_x, new_y, prod.rod_type, prod.on)
                
                topcomponent_intersection_area = prod.tip_component_buffer.intersection(topcomponents).area
                top_probe_intersection_area = prod.tip_from_top_probe_buffer.intersection(top_probes).area
                
                # print(f"\n{topcomponent_intersection_area} intersection with top components")
                # print(f"{top_probe_intersection_area} intersection with top probes")
                
                if not prod.center.intersects(pBoards_multi):
                    valid = False
                    continue                
                if prod.tip_component_buffer.intersects(topcomponents):
                    valid = False
                    continue
                if prod.tip_from_top_probe_buffer.intersects(top_probes):
                    valid = False
                    continue
                
                if valid == True:
                    break
                
        if not prod.center.intersects(pBoards_multi):
            valid = False
            continue
            
        
        
        # ## END NEW CODE ###
        
        # # ## ORIGINAL CODE ###
        # if prod.tip_component_buffer.intersects(topcomponents):
        #     continue
        
        # # Make sure the pressure rod isn't too close to any top probes
        # if prod.tip_from_top_probe_buffer.intersects(top_probes):
        #     continue
            
        # Make sure pressure rod doesn't conflict with any previously-placed pressure rods
        if len(prods_chosen) > 0:
            for prod_chosen in prods_chosen:
                if prod_chosen.top.intersects(prod.top_from_top_probe_buffer):
                    valid = False
                    break
                # dist = centroid_distance(prod.tip,prod_chosen.tip)
                # if dist < prod.ctc:
                #     valid = False
                #     break
            if valid == False:
                continue
        ### END ORIGINAL CODE ###
        
        if valid == True:
            prods_chosen.append(prod)
            chromosome_x.append(prod.x)
            chromosome_y.append(prod.y)
            chromosome_rod_type.append(rod_type_i)
            chromosome_on.append(on)
            
    chromosome = []
    chromosome.extend(chromosome_x)
    chromosome.extend(chromosome_y)
    chromosome.extend(chromosome_rod_type)
    chromosome.extend(chromosome_on)

    return chromosome

# def validate_chromosome(chromosome, nprods, top_constraints):
#     prods = interpret_chromosome_to_prods(chromosome, nprods)
#     pBoards_multi = top_constraints[0]
#     top_probes = top_constraints[1]
#     topcomponents = top_constraints[2]
#     xmin, ymin, xmax, ymax = pBoards_multi.bounds

def get_top_constraints(pBoards, pComponentsTop, df_Probes, pBoards_diff):
    pBoards_multi = MultiPolygon(pBoards_diff)
    
    df_Probes_top = df_Probes[df_Probes["side"]==1]
    df_Probes_top.reset_index(inplace=True)
    top_probes = []
    for i,row in df_Probes_top.iterrows():
        top_probes.append(place_circle(row.x, row.y, row.diameter/2))
        
    top_probes = unary_union(top_probes)
    if top_probes.geom_type != "MultiPolygon":
        top_probes = MultiPolygon(top_probes)
    
    topcomponents = []
    for inner in pComponentsTop:
        topcomponents.append(inner)
    
    topcomponents = unary_union(topcomponents)
    if topcomponents.geom_type != "MultiPolygon":
        topcomponents = MultiPolygon(topcomponents)
    
    top_constraints = (pBoards_multi, top_probes, topcomponents)
    return top_constraints

def plot_prods_top_constraints(prods, top_constraints, title):
    # Quick plot to visually check designs as they are produced by crossover
    # or other relevant functions
    fig, ax = plt.subplots(dpi=300, figsize=(10,8))
    ax.set_aspect('equal')
    
    # UUT boundaries
    color = 'k'
    linestyle = '-'
    label = "UUT"
    plot_multipolygon_w_holes(top_constraints[0], fig, ax, color, linestyle, label)
    
    # Top probes
    color = 'b'
    label = "Top Probes"
    plot_multipolygon_w_holes(top_constraints[1], fig, ax, color, linestyle, label)
    
    # Top components
    color = 'purple'
    label = "Top Components"
    plot_multipolygon_w_holes(top_constraints[2], fig, ax, color, linestyle, label)
    
    # "On" pressure rods with offsets
    prods_poly_list_on = [prod.tip for prod in prods if prod.on]
    prods_tip_component_buffer_on = [prod.tip_component_buffer for prod in prods if prod.on]
    prods_tip_UUT_buffer_on = [prod.tip_UUT_buffer for prod in prods if prod.on]
    prods_tip_from_top_probe_buffer_on = [prod.tip_from_top_probe_buffer for prod in prods if prod.on]
    linestyle = '-'
    label = "Pressure Rods - On"
    plot_poly_list_w_holes(prods_poly_list_on, fig, ax, 'g', linestyle, label)
    plot_poly_list_w_holes(prods_tip_UUT_buffer_on, fig, ax, 'lime', linestyle, label)
    plot_poly_list_w_holes(prods_tip_component_buffer_on, fig, ax, 'turquoise', linestyle, label)
    plot_poly_list_w_holes(prods_tip_from_top_probe_buffer_on, fig, ax, 'teal', linestyle, label)
    
    # "Off pressure rods with offsets
    prods_poly_list_off = [prod.tip for prod in prods if not prod.on]
    prods_tip_component_buffer_off = [prod.tip_component_buffer for prod in prods if not prod.on]
    prods_tip_UUT_buffer_off = [prod.tip_UUT_buffer for prod in prods if not prod.on]
    prods_tip_from_top_probe_buffer_off = [prod.tip_from_top_probe_buffer for prod in prods if not prod.on]
    
    linestyle = '-'
    label = "Pressure Rods - Off"
    plot_poly_list_w_holes(prods_poly_list_off, fig, ax, 'r', linestyle, label)
    plot_poly_list_w_holes(prods_tip_UUT_buffer_off, fig, ax, 'yellow', linestyle, label)
    plot_poly_list_w_holes(prods_tip_component_buffer_off, fig, ax, 'gold', linestyle, label)
    plot_poly_list_w_holes(prods_tip_from_top_probe_buffer_off, fig, ax, 'darkorange', linestyle, label)
    
    fig.suptitle(title)
    
    
    
def validate_prod(prod, prods_chosen, top_constraints):
    # Validate a single pressure rod for the top side
    pBoards_multi = top_constraints[0]
    top_probes = top_constraints[1]
    topcomponents = top_constraints[2]
    xmin, ymin, xmax, ymax = pBoards_multi.bounds
    
    valid = True
    while True:        
        # Make sure pressure rod is within the UUT and make sure it doesn't intersect any components, using the appropriate buffer sizes
        if not prod.center.intersects(pBoards_multi):
            valid = False
            break
        if prod.tip_component_buffer.intersects(topcomponents):
            valid = False
            break
        
        # Make sure the pressure rod isn't too close to any top probes
        if prod.tip_from_top_probe_buffer.intersects(top_probes):
            valid = False
            break

        # If the pressure rod is turned off, it's automatically valid
        if prod.on == 0:
            valid = True
            break
            
        # Make sure pressure rod doesn't conflict with any previously-placed pressure rods
        for prod_chosen in prods_chosen:
            if prod_chosen.top.intersects(prod.top_from_top_probe_buffer):
                valid = False
                break
            # dist = centroid_distance(prod.tip,prod_chosen.tip)
            # if dist < prod.ctc:
            #     valid = False
            #     break
        if valid == False:
            break
        
        if valid == True:
            break
            
    return valid

def validate_prods(prods_chosen, top_constraints):
    for prod in prods_chosen:
        valid = validate_prod(prod, prods_chosen, top_constraints)
        if valid == False:
            break
    return valid


# def worker_function(result_queue,nprods,pBoards, pComponentsTop, df_Probes, pBoards_diff):
#     while True:
#         chromosome = create_chromosome(nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff)
#         result_queue.put(chromosome)

        
# def initialize_population_multiprocessing(nchromosomes, nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff):
#     num_processes = multiprocessing.cpu_count()  # Number of available CPU cores
#     result_queue = multiprocessing.Queue()
    
#     # Create and start worker processes
#     processes = []
#     for _ in range(num_processes):
#         process = multiprocessing.Process(target=worker_function, args=(result_queue, nprods,pBoards, pComponentsTop, df_Probes, pBoards_diff))
#         process.start()
#         processes.append(process)
        
#     # Wait for all worker processes to finish
#     for process in processes:
#         process.join()
        
#     # Retrieve results from the result queue
#     initial_population = []
#     while not result_queue.empty():
#         chromosome = result_queue.get()
#         initial_population.append(chromosome)
        
#     return initial_population

def initialize_population_simple(npop, nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff, all_on, on_prob, rod_type):
    initial_population = []
    for i in range(npop):
        initial_population.append(create_chromosome(nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff, all_on, on_prob, rod_type))
        print(f"Chromosome {i} of {npop} created")
    return initial_population

def interpret_chromosome_to_prods(chromosome, nprods):
    
    def interpret_rod_type_int(rod_type_int):
        rod_types = ['Press-Fit Tapered',
                     'Press-Fit Flat',
                     '3.325" Tapered',
                     '3.325" Flat']
        return rod_types[rod_type_int]
    
    prods = []
    for i in range(nprods):
        x = chromosome[i]
        y = chromosome[i+nprods]
        rod_type_int = int(chromosome[i+2*nprods])
        rod_type = interpret_rod_type_int(rod_type_int)
        on = chromosome[i+3*nprods]
        prod = PressureRod(x, y, rod_type, on)
        prods.append(prod)
    
    return prods

def interpret_prods_to_chromosome(prods):
    
    def interpret_rod_type(rod_type):
        rod_types = ['Press-Fit Tapered',
                     'Press-Fit Flat',
                     '3.325" Tapered',
                     '3.325" Flat']
        
        return rod_types.index(rod_type)
        
    chromosome_x = []
    chromosome_y = []
    chromosome_rod_type = []
    chromosome_on = []
    for prod in prods:
        chromosome_x.append(prod.x)
        chromosome_y.append(prod.y)
        chromosome_rod_type.append(interpret_rod_type(prod.rod_type))
        chromosome_on.append(prod.on)
            
    chromosome = []
    chromosome.extend(chromosome_x)
    chromosome.extend(chromosome_y)
    chromosome.extend(chromosome_rod_type)
    chromosome.extend(chromosome_on)
    
    return chromosome
    

def prods_to_valid_circles(prods):
    valid_circles = [prod.tip for prod in prods if prod.on]
    return valid_circles


# %% FEA functions
def runFEA_valid_circles(valid_circles, df_PressureRods, root, inputfile, gen, iteration):
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
    
    # Set the salesOrder field (used in naming the output folder) to include
    # generation and iteration numbers
    salesOrder = root.find('.//salesOrder')
    salesOrder.text = f"GEN{gen}_ITER{iteration}"
    
    # Write the new element tree
    tree = ET.ElementTree(root)
    ET.indent(tree, '  ')
    new_filename = f"FEA_GEN{gen}_ITER{iteration}.xml"
    path, filename = os.path.split(inputfile)
    new_path = path + "/" + new_filename
    tree.write(new_path)
    
    FEApath = runFEA.loadFEApath('FEApath.pk')
    exit_code = runFEA.runFEA(FEApath, new_path)
    
    # results = (gen, iteration)
    
    return (exit_code)
    

def read_FEA_results(root, inputfile, gen, iteration):
    # Write the new element tree
    tree = ET.ElementTree(root)
    ET.indent(tree, '  ')
    new_filename = f"FEA_GEN{gen}_ITER{iteration}.xml"
    path, filename = os.path.split(inputfile)
    new_path = path + "/" + new_filename
    
    dfmesh = runFEA.resultsToDataframe(new_path)
    
    strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max = runFEA.getFitness(dfmesh)
    
    results = (strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max)
    
    return results
    

def design_to_xml(valid_circles, df_PressureRods, root, inputfile, gen, iteration):
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
    
    # Set the salesOrder field (used in naming the output folder) to include
    # generation and iteration numbers
    salesOrder = root.find('.//salesOrder')
    salesOrder.text = f"GEN{gen}_ITER{iteration}"
    
    # Write the new element tree
    tree = ET.ElementTree(root)
    ET.indent(tree, '  ')
    new_filename = f"FEA_GEN{gen}_ITER{iteration}.xml"
    path, filename = os.path.split(inputfile)
    new_path = path + "/" + new_filename
    tree.write(new_path)
    
    return new_path


def runFEA_new_path(new_path_xml):
    # Simplified and separated method that should run off a separately created XML file
    FEApath = runFEA.loadFEApath('FEApath.pk')
    runFEA.runFEA(FEApath, new_path_xml)
    
    dfmesh = runFEA.resultsToDataframe(new_path_xml)
    
    strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max = runFEA.getFitness(dfmesh)
    
    results = (strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max)
    
    return results
    


# %% Utility functions
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
            
def plot_multipolygon_w_holes(multipoly, fig, ax, color, linestyle, label):
    first = True
    if multipoly.geom_type == "MultiPolygon":
        for geom in multipoly.geoms:
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
        raise IOError("Shape is not a multipolygon")
            
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

def centroid_distance(poly1,poly2):
    x1 = poly1.centroid.x
    y1 = poly1.centroid.y
    x2 = poly2.centroid.x
    y2 = poly2.centroid.y
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

# %% Load data and read in the XML definition
if __name__ == "__main__":
    results = get_constraint_geometry()
    root = results[0]
    inputfile = results[1]
    pBoards = results[2]
    pOutline = results[3]
    pShape = results[4]
    pComponentsTop = results[5]
    pComponentsBot = results[6]
    Pressure = results[7]
    I_Plate = results[8]
    Stripper = results[9]
    Probe = results[10]
    Countersink = results[11]
    df_Probes = results[12]
    df_GuidePins = results[13]
    df_PressureRods = results[14]
    df_Standoffs = results[15]
    
    # Estimate a number of pressure rods for the top side that would make sense
    print("--- Estimating possible pressure rods ---")
    nprods_small, nprods_large, pBoards_diff = grid_nprods(pBoards,pComponentsTop)
    
    # Choose the largest number of variables between len(df_PressureRods),
    # nprods_small, and nprods_large
    print("--- Selecting pressure rod quantity ---")
    nprods = 20
    # nprods = np.max([len(df_PressureRods), nprods_small, nprods_large])
    
    # Generate random population of pressure rod designs
    npop = 10
    start_time = time.time()
    print(f"--- Generating population of {npop} chromosomes with {nprods*4} variables each---")
    # chromosomes = [create_chromosome(nprods,pBoards,pComponentsTop,df_Probes,pBoards_diff) for _ in range(npop)] # Could this be modified to use multiprocessing? This will become very time intensive with larger populations
    initial_population = initialize_population_simple(npop, nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff)
    print(f"--- Population generated in {(time.time()-start_time)/60} minutes ---")
    # print("--- COMPLETE ---")
    
    # start_time = time.time()
    # print(f"--- Generating population of {npop} chromosomes with {nprods*4} variables each---")
    # initial_population = initialize_population_multiprocessing(npop, nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff)
    # print(f"--- Population generated in {(time.time()-start_time)/60} minutes USING MULTIPROCESSING ---")
    
    # print("--- COMPLETE ---")
    
    
