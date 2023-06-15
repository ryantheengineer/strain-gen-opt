# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 22:27:20 2023

@author: Ryan Larson (ChatGPT code assist)
"""

import random
from shapely.geometry import Point, Polygon

def place_circle(x,y,radius):
    circle = Point(x, y).buffer(radius)
    return circle

def UUT_square(maxbnd):
    coords = ((0.,0.),(0.,maxbnd),(maxbnd,maxbnd),(maxbnd,0.),(0.,0.))
    poly = Polygon(coords)
    return poly

def generate_random_circle(main_polygon, circles, circle_size, min_dist):
    while True:
        # Generate a random point within the main polygon
        x = random.uniform(main_polygon.bounds[0], main_polygon.bounds[2])
        y = random.uniform(main_polygon.bounds[1], main_polygon.bounds[3])

        # Create the new circle
        new_circle = place_circle(x, y, circle_size)

        # Check if the new circle overlaps any previously placed circles
        if any(existing_circle.intersects(new_circle) for existing_circle in circles):
            continue

        # Check if the new circle satisfies the minimum distance (offset) requirement
        if any(existing_circle.distance(new_circle) < min_dist for existing_circle in circles):
            continue

        return new_circle

def generate_maximum_circles_in_polygon(circle_size, min_dist, main_polygon, num_holes=0):
    # Generate the main polygon
    maxbnd = 12.
    main_polygon = UUT_square(maxbnd)

    # # Generate random holes within the main polygon
    # for _ in range(num_holes):
    #     hole_coords = [...]  # Specify the coordinates of each hole's vertices
    #     hole_polygon = Polygon(hole_coords)
    #     main_polygon = main_polygon.difference(hole_polygon)

    # Generate circles within the main polygon
    circles = []
    while True:
        new_circle = generate_random_circle(main_polygon, circles, circle_size, min_dist)
        if new_circle is None:
            break

        circles.append(new_circle)

    return circles

# Example usage
circle_size = 1.0
min_dist = 0.5
num_holes = 0

circles = generate_maximum_circles_in_polygon(circle_size, min_dist, num_holes)
print("Number of circles generated:", len(circles))
