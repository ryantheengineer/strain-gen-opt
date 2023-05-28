# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
import pandas as pd
import random

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




if __name__ == "__main__":
    
    # Create a random UUT board with some rectangular features that should be
    # avoided by circles. Plot these in blue. Create a series of random circles
    # to represent pins that do not conflict with any features in the UUT.
    maxbnd = 100.
    n_rectangles = 40
    min_rect = 2.
    max_rect = 20.
    
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
    
    fig, ax = plt.subplots()
    
    # Plot Polygon
    xe, ye = UUT_poly.exterior.xy
    for inner in UUT_poly.interiors:
        xi, yi = zip(*inner.coords[:])
        ax.plot(xi, yi, color="blue")
     
    ax.plot(xe, ye, color="blue")
    # ax.axis([0, 100, 0, 100])
    plt.show()
    
     
    # # If there are any Interiors
    # # Retrieve coordinates for all interiors
    # for inner in UUT_poly.interiors:
    #     xi, yi = zip(*inner.coords[:])
    #     ax.plot(xi, yi, color="blue")
     
    # ax.plot(xe, ye, color="blue")
    # # ax.axis([0, 100, 0, 100])
    # plt.show()