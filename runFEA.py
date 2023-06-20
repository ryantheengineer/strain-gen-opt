# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 22:53:50 2023

@author: Ryan Larson
"""

import subprocess
import os, glob
import pathlib
import pandas as pd
from tkinter import Tk
from tkinter import filedialog as fd



def runFEA(FEApath, inputfile):
    args = [FEApath, inputfile]
    subprocess.call(args, shell=False)

def resultsToDataframe(inputfile):
    directory = pathlib.Path(inputfile)
    directory = str(directory.parent) + "\Output"
    
    # Get the most recently modified subdirectory
    latest_subdir = max(glob.glob(os.path.join(directory, '*/')), key=os.path.getmtime)
    
    meshfile = latest_subdir + "FEA_MeshNodes.csv"
    
    df = pd.read_csv(meshfile)
    return df

def getFitness(dfmesh):
    abssums = dfmesh.abs().sum()
    strain_xx = abssums["strain_xx"]
    strain_yy = abssums["strain_yy"]
    strain_xy = abssums["strain_xy"]
    principalStrain_min = abssums["principalStrain_min"]
    principalStrain_max = abssums["principalStrain_max"]
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max
    

if __name__ == "__main__":
    
    FEApath = "C:/Users/Ryan Larson/github/strain-gen-opt/FEA/FEA.exe"
    
    initialdir = str(pathlib.Path(FEApath).parent) + "Examples"
    
    filetypes = (("XML", ["*.xml"]),)
    
    root = Tk()
    root.wm_attributes('-topmost', 1)
    inputfile = fd.askopenfilename(
            title="Select FEA input file",
            initialdir=initialdir,
            filetypes=filetypes
            )
    root.destroy()
    
    # inputfile = "C:/Users/Ryan Larson/github/strain-gen-opt/FEA/Examples/Example 1/FEA.xml"
    
    runFEA(FEApath, inputfile)
    
    dfmesh = resultsToDataframe(inputfile)
    
    strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max = getFitness(dfmesh)