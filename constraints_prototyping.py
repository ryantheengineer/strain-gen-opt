# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
import pandas as pd
import random

def random_rectangle(mindim,maxdim,maxbnd):
    pass




if __name__ == "__main__":
    maxbnd = 40.
    coords1 = ((0.,0.),(0.,maxbnd),(maxbnd,maxbnd),(maxbnd,0.),(0.,0.))
    coords2 = ((1.,1.),(1.,2.),(2.,4.),(1.,1.)) # interior points
    
    poly1 = Polygon(coords1, [coords2])
    
    fig,ax = plt.subplot()
    