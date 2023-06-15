# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 23:14:11 2023

@author: Ryan Larson
"""

from shapely.geometry import Point, Polygon, MultiPolygon
from shapely import buffer
import matplotlib.pyplot as plt
import random


def get_circle(x,y,radius):
    circle = Point(x, y).buffer(radius)
    return circle

def interpret_chromosome(ncircles,chromosome):
    # Interpret the chromosome and create circles (Polygons) at the proper
    # positions and sizes
    circles = []
    for i in range(ncircles):
        feature_list = []
        for j in range(len(chromosome)):
            if j % ncircles == i:
                feature_list.append(chromosome[j])
        x = feature_list[0]
        y = feature_list[1]
        radius = feature_list[2]
        on_off = feature_list[3]
        print(f"\nCircle {i}:")
        print(f"x = {x}")
        print(f"y = {y}")
        print(f"radius = {radius}")
        if on_off == 0:
            print("Status = OFF")
        else:
            print("Status = ON")
        circle = get_circle(x, y, radius)
        circles.append(circle)
        
    return circles


if __name__ == "__main__":
    xmin = 0.
    xmax = 10.
    ymin = 10.
    ymax = 20.
    radius_sizes = [0.1, 0.2, 0.3]
    ncircles = 10
    
    # Initialize random chromosome, where the first ncircles entries are
    # the x coordinates, the next ncircles entries are y coordinates, then
    # radii, and on/off binary values
    chromosome = []
    chromosome.extend([random.uniform(xmin,xmax) for _ in range(ncircles)])
    chromosome.extend([random.uniform(ymin,ymax) for _ in range(ncircles)])
    chromosome.extend([random.choice(radius_sizes) for _ in range(ncircles)])
    chromosome.extend([random.randint(0,1) for _ in range(ncircles)])
    
    # Interpret the chromosome into the individual circles it defines
    circles = interpret_chromosome(ncircles,chromosome)
        
    # Plot circles
    fig, ax = plt.subplots(dpi=300)
    ax.axis('equal')
    for circle in circles:
        xe, ye = circle.exterior.xy
        ax.plot(xe, ye, color="red")
    plt.legend()
    plt.show()