# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.prepared import prep
from shapely import buffer
from shapely.ops import unary_union
import matplotlib.pyplot as plt
# import pandas as pd
import random
import numpy as np


class ForceRod:
    def __init__(self,x,y,radius,offset):
        self.x = x
        self.y = y
        self.radius = radius
        self.offset = offset
        self.circle = place_circle(x,y,radius)
        self.bufferzone = buffer(self.circle,self.offset)
        
    def perturb_rod(self,maxperturb):
        while True:
            xp = random.uniform(-maxperturb, maxperturb)
            yp = random.uniform(-maxperturb, maxperturb)
            magnitude = np.sqrt(xp**2 + yp**2)
            if magnitude > maxperturb:
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
        self.circle = poly_perturb
        self.bufferzone = buffer(self.circle,self.offset)
    

def UUT_square(maxbnd):
    coords = ((0.,0.),(0.,maxbnd),(maxbnd,maxbnd),(maxbnd,0.),(0.,0.))
    poly = Polygon(coords)
    return poly

def random_rectangle(mindim,maxdim,maxbnd):
    w = random.uniform(mindim,maxdim)
    h = random.uniform(mindim,maxdim)
    
    xlc = random.uniform(0.0, maxbnd-w)
    ylc = random.uniform(h, maxbnd)
    
    coords = ((xlc,ylc),(xlc+w,ylc),(xlc+w,ylc-h),(xlc,ylc-h),(xlc,ylc))
    poly = Polygon(coords)
    return poly

def random_circle(radius,maxbnd):
    x = random.uniform(radius,maxbnd-radius)
    y = random.uniform(radius,maxbnd-radius)
    circle = Point(x, y).buffer(radius)
    return circle

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
        

if __name__ == "__main__":
    
    # Create a random UUT board with some rectangular features that should be
    # avoided by circles. Plot these in blue. Create a series of random circles
    # to represent pins that do not conflict with any features in the UUT.
    maxbnd = 12.
    n_rectangles = 30
    min_rect = 0.2
    max_rect = 2.0
    n_circles = 10
    radius = 0.125
    resolution = 0.375
    perturb = 0.1
    
    UUT_poly = UUT_square(maxbnd)
    
    poly_int = random_rectangle(min_rect, max_rect, maxbnd)
    
    interiors = [poly_int]
    
    for i in range(n_rectangles-1):
        while True:
            p = random_rectangle(min_rect, max_rect, maxbnd)
            conflict = False
            for interior in interiors:
                if interior.intersects(p):
                    conflict = True
                    break
            if conflict == False:
                interiors.append(p)
                break
            else:
                continue
                
    for interior in interiors:
        UUT_poly = UUT_poly.difference(interior)
        
    # Dilate UUT so no variable components can be placed in incompatible
    # locations
    bufferdist_edge = radius
    bufferdist_component = 0.0625
    UUT_poly_ext = Polygon(UUT_poly.exterior.coords)
    UUT_poly_dilated_ext = buffer(UUT_poly_ext,-bufferdist_edge)
    interiors_dilated = []
    for inner in UUT_poly.interiors:
        UUT_poly_int = Polygon(inner.coords)
        interiors_dilated.append(buffer(inner,bufferdist_component))
    
    print(f"Before:\t{len(interiors_dilated)}")
    # Unary union of dilated interiors
    interiors_dilated = unary_union(interiors_dilated)
    interiors_dilated = list(interiors_dilated.geoms)
    print(f"After:\t{len(interiors_dilated)}")
    
    # Subtract dilated interiors from dilated exterior
    UUT_poly_dilated = UUT_poly_dilated_ext
    for inner in interiors_dilated:
        inner = Polygon(inner.exterior.coords)
        UUT_poly_dilated = UUT_poly_dilated.difference(inner)
    
    
    # Place circles on grid within UUT
    xmin, ymin, xmax, ymax = UUT_poly.bounds
    circles = []
    # NOTE: 
    for x in np.arange(xmin, xmax, resolution):
        for y in np.arange(ymin, ymax, resolution):
            c = place_circle(x, y, radius)
            circles.append(c)
    
    
    # valid_circles.extend(filter(UUT_poly.contains, circles))
    
    # Randomly perturb each circle
    for i,circle in enumerate(circles):
        circles[i] = random_perturb(circle, perturb)
        # NOTE: Before accepting randomly perturbed circle, make sure it is 
        # within spec. Force rods can't be too close together.
    
    # Validate that each circle falls inside shape
    valid_circles = []
    valid_circles.extend(filter(UUT_poly_dilated.contains, circles))
    # NOTE: Add check here to make sure no circles overlap with other circles
    
    # valid_circles = list(filter(UUT_poly.contains, valid_circles))
    
    # # Randomly place circles
    # circles = []
    # for i in range(n_circles):
    #     while True:
    #         p = random_circle(radius, maxbnd)
    #         conflict = False
    #         for circle in circles:
    #             if circle.intersects(p):
    #                 conflict = True
    #                 break
    #             elif UUT_poly.contains(circle) == False:
    #                 conflict = True
    #                 break
    #         if conflict == False:
    #             circles.append(p)
    #             break
    #         else:
    #             continue
    
    
    fig, ax = plt.subplots(dpi=300)
    ax.axis('equal')
    
    # Plot Polygon
    xe, ye = UUT_poly.exterior.xy
    for inner in UUT_poly.interiors:
        xi, yi = zip(*inner.coords[:])
        ax.plot(xi, yi, color="blue")
     
    ax.plot(xe, ye, color="blue", label="UUT")
    
    # Plot Dilated Polygon
    if UUT_poly_dilated.geom_type == 'MultiPolygon':
        UUT_poly_dilated_polys = list(UUT_poly_dilated.geoms)
        for poly in UUT_poly_dilated_polys:
            xe, ye = poly.exterior.xy
            ax.plot(xe, ye, color="green", label="UUT Offset")
            for inner in poly.interiors:
                xi, yi = zip(*inner.coords[:])
                ax.plot(xi, yi, color="green")
    elif UUT_poly_dilated.geom_type == 'Polygon':
        xe, ye = UUT_poly_dilated.exterior.xy
        ax.plot(xe, ye, color="green", label="UUT Offset")
        for inner in UUT_poly_dilated.interiors:
            xi, yi = zip(*inner.coords[:])
            ax.plot(xi, yi, color="green")
    else:
        raise IOError('Dilated UUT boundary is not a polygon.')
    
    # xe, ye = UUT_poly_dilated.exterior.xy
    # for inner in UUT_poly_dilated.interiors:
    #     xi, yi = zip(*inner.coords[:])
    #     ax.plot(xi, yi, color="green")
     
    
    # xe, ye = UUT_poly_dilated_ext.exterior.xy
    # for inner in interiors_dilated:
    #     xi, yi = inner.exterior.xy
    #     ax.plot(xi, yi, color="green")
     
    # ax.plot(xe, ye, color="green", label="UUT Offset")
    
    
    # Plot circles
    for circle in valid_circles:
        xe, ye = circle.exterior.xy
        ax.plot(xe, ye, color="red")
    plt.legend()
    plt.show()
    
     
    # # If there are any Interiors
    # # Retrieve coordinates for all interiors
    # for inner in UUT_poly.interiors:
    #     xi, yi = zip(*inner.coords[:])
    #     ax.plot(xi, yi, color="blue")
     
    # ax.plot(xe, ye, color="blue")
    # # ax.axis([0, 100, 0, 100])
    # plt.show()