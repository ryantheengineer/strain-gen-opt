# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 23:20:26 2023

@author: Ryan Larson
"""
import xmltodict



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
    plate_types = [plate["@identifier"] for plate in plates]








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


