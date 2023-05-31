# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from shapely.geometry import Point, Polygon
from shapely.prepared import prep
import matplotlib.pyplot as plt
# import pandas as pd
import random
import numpy as np

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


if __name__ == "__main__":
    
    # Create a random UUT board with some rectangular features that should be
    # avoided by circles. Plot these in blue. Create a series of random circles
    # to represent pins that do not conflict with any features in the UUT.
    maxbnd = 100.
    n_rectangles = 20
    min_rect = 2.
    max_rect = 20.
    n_circles = 10
    radius = 2.
    resolution = 5.
    
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
        
    
    # Place circles on grid within UUT
    xmin, ymin, xmax, ymax = UUT_poly.bounds
    
    prep_UUT_poly = prep(UUT_poly)
    
    circles = []
    for x in np.arange(xmin, xmax, resolution):
        for y in np.arange(ymin, ymax, resolution):
            c = place_circle(x, y, radius)
            circles.append(c)
    
    # Validate that each circle falls inside shape
    valid_circles = []
    valid_circles.extend(filter(UUT_poly.contains, circles))
    
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
    
    
    fig, ax = plt.subplots()
    
    # Plot Polygon
    xe, ye = UUT_poly.exterior.xy
    for inner in UUT_poly.interiors:
        xi, yi = zip(*inner.coords[:])
        ax.plot(xi, yi, color="blue")
     
    ax.plot(xe, ye, color="blue")
    
    # Plot circles
    for circle in valid_circles:
        xe, ye = circle.exterior.xy
        ax.plot(xe, ye, color="red")
    plt.show()
    
     
    # # If there are any Interiors
    # # Retrieve coordinates for all interiors
    # for inner in UUT_poly.interiors:
    #     xi, yi = zip(*inner.coords[:])
    #     ax.plot(xi, yi, color="blue")
     
    # ax.plot(xe, ye, color="blue")
    # # ax.axis([0, 100, 0, 100])
    # plt.show()