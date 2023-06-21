# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 23:20:26 2023

@author: Ryan Larson
"""
import xmltodict
from shapely.geometry import Polygon
import matplotlib.pyplot as plt


if __name__ == "__main__":
    with open('FEA_Example2.xml', 'r', encoding='utf-8') as file:
        FEA_xml = file.read()
        
    FEA_dict = xmltodict.parse(FEA_xml)
    
    FEA_parameters = FEA_dict["FEA"]["struct"][0]
    
    fixture_dict = FEA_dict["FEA"]["struct"][1]
    
    nsprings = int(fixture_dict["SpringCount"])
    
    # Get relevant information about panels
    panels = fixture_dict["struct"][0]["polyshape"]
    panel_types = [panel["@identifier"] for panel in panels]
    t_panel = fixture_dict["struct"][0]["thickness"]
    xMax_panel = fixture_dict["struct"][0]["xMax"]
    xMin_panel = fixture_dict["struct"][0]["xMin"]
    yMax_panel = fixture_dict["struct"][0]["yMax"]
    yMin_panel = fixture_dict["struct"][0]["yMin"]
    
    for panel in panels:
        if 'pBoards' in panel.values():
            pBoards_vertices = panel['region']['vertex']
            break
        else:
            pBoards_vertices = None
    
    # Get relevant information about plates
    plates = fixture_dict["struct"][1]["polyshape"]
    plates_present = [plate["@identifier"] for plate in plates]
    plate_types = ["Pressure", "I_Plate", "Stripper", "Probe"]
    
    plate_polygons = []
    for platenum,plate in enumerate(plates):
        # Get the exterior boundary of the plate
        plate_region = plate["region"]["vertex"]
        for i,vertex in enumerate(plate_region):
            plate_region[i] = vertex.split("|")
        
        # Get the hole boundaries, if there are any
        plate_holes = plate["hole"]
        for i,hole in enumerate(plate_holes):
            plate_holes[i] = hole["vertex"]
            for j,vertex in enumerate(plate_holes[i]):
                plate_holes[i][j] = vertex.split("|")
                
        # Create a shapely Polygon from the exterior and the holes
        plate_poly = Polygon(plate_region)
        for hole in plate_holes:
            hole_poly = Polygon(hole)
            plate_poly = plate_poly.difference(hole_poly)
        plate_polygons.append(plate_poly)
        
        # Plot the plate
        fig, ax = plt.subplots(dpi=300)
        ax.axis('equal')
        
        # Plot Polygon
        xe, ye = plate_poly.exterior.xy
        for inner in plate_poly.interiors:
            xi, yi = zip(*inner.coords[:])
            ax.plot(xi, yi, color="blue")
         
        ax.plot(xe, ye, color="blue", label=f"{plates_present[platenum]}")
        plt.legend()
        plt.show()
    
    # NOTE: Need to respond differently depending on which plates are part of the design






# import pandas as pd
# import xml.etree.ElementTree as ET

# # Load the XML file
# tree = ET.parse('FEA_Example2.xml')

# # Get the root element
# root = tree.getroot()

# # Create a DataFrame
# df = pd.DataFrame()

# # Iterate over the child elements of the root element
# for child in root:
#     # Add the child element's text to the DataFrame
#     df[child.tag] = child.text

# # Print the DataFrame
# print(df)


