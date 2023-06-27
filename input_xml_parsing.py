# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 23:20:26 2023

@author: Ryan Larson
"""
import xmltodict
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Use for inputs that are an XML description of a data table. These include
# the inputs for probes, guide pins, pressure rods on top and bottom, and
# standoffs (probes_dict, guide_pins_dict, pressure_rods_dict,
# bottom_pressure_rods_dict, and standoffs_dict).
def feature_dict_to_dataframe(feature_dict):
    cols = feature_dict["columns"]["col"]
    rows = feature_dict["rows"]["row"]
    new_rows = [elem.split("|") for elem in rows]
    identifier = feature_dict["@identifier"]
    
    if identifier == "Probes":
        int_indices = [0,1,4,6,8,10,16]
        float_indices = [2,3,5,13,14]
        bool_indices = [12]
        list_indices = [15]
        nan_indices = [11]
    elif identifier == "GuidePins":
        int_indices = [0,1,2,3,5,9,11,12,15,19]
        float_indices = [6,7,8,10,16]
        bool_indices = [4,14,18]
        list_indices = [17]
        nan_indices = []
    elif identifier == "PressureRods":
        int_indices = [0,1,2,3,5,8,9,11,12,15]
        float_indices = [6,7,10,16]
        bool_indices = [4,14]
        list_indices = [17]
        nan_indices = []
    elif identifier == "BottomPressureRods":
        int_indices = []
        float_indices = []
        bool_indices = []
        list_indices = []
        nan_indices = []
    elif identifier == "Standoffs":
        int_indices = [0,1,2,3,5,8,9,11,12,15]
        float_indices = [6,7,10,16]
        bool_indices = [4,14]
        list_indices = [17]
        nan_indices = []
    else:
        raise ValueError("Identifier not recognized")
    
    # Convert elements of new_rows to the appropriate data type before
    # converting to DataFrame
    for row in new_rows:
        for i,elem in enumerate(row):
            if i in int_indices:
                row[i] = int(elem)
            elif i in float_indices:
                row[i] = float(elem)
            elif i in bool_indices:
                if elem == "false":
                    row[i] = False
                else:
                    row[i] = True
            elif i in list_indices:
                row[i] = elem.strip('][').split(', ')
            elif i in nan_indices:
                row[i] = np.nan
    
    df = pd.DataFrame(new_rows, columns=cols)
    return df

if __name__ == "__main__":
    with open('FEA_Example2.xml', 'r', encoding='utf-8') as file:
        FEA_xml = file.read()
        
    FEA_dict = xmltodict.parse(FEA_xml)
    
    FEA_parameters = FEA_dict["FEA"]["struct"][0]
    
    fixture_dict = FEA_dict["FEA"]["struct"][1]
    
    nsprings = int(fixture_dict["SpringCount"])
    
    # Get dicts of features of interest (design variables)
    probes_dict = fixture_dict["table"][0]
    guide_pins_dict = fixture_dict["table"][1]
    pressure_rods_dict = fixture_dict["table"][2]
    bottom_pressure_rods_dict = fixture_dict["table"][3]
    standoffs_dict = fixture_dict["table"][4]
    
    # Convert the dicts to dataframes
    probes_df = feature_dict_to_dataframe(probes_dict)
    guide_pins_df = feature_dict_to_dataframe(guide_pins_dict)
    pressure_rods_df = feature_dict_to_dataframe(pressure_rods_dict)
    # bottom_pressure_rods_df = feature_dict_to_dataframe(bottom_pressure_rods_dict)
    standoffs_df = feature_dict_to_dataframe(standoffs_dict)
    
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

    # NOTE: Down the line it would be useful to have several different plotting
    # functions:
        # Plot a specific chromosome (individual design)
        # Plot the initial design
        # Show all plate designs simultaneously or in a way that they can be compared
        # 




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


