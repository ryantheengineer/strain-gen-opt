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
import pickle

def chooseFEApath(initialdir):
    filetypes = (("Executable", ["*.exe"]),)
    
    root = Tk()
    root.wm_attributes('-topmost', 1)
    FEApath = fd.askopenfilename(
            title="Select FEA executable",
            initialdir=initialdir,
            filetypes=filetypes
            )
    root.destroy()
    
    filename = 'FEApath.pk'
    
    with open(filename, 'wb') as fi:
        # dump data into file
        pickle.dump(FEApath, fi)
    return FEApath

def loadFEApath(filename):
    with open(filename, 'rb') as fi:
        FEApath = pickle.load(fi)
    return FEApath

def runFEA(FEApath, inputfile):
    args = [FEApath, inputfile]
    subprocess.call(args, shell=False)

def resultsToDataframe(inputfile):
    directory = pathlib.Path(inputfile)
    directory = str(directory.parent) + "\Output"
    path, filename = os.path.split(inputfile)
    filename = os.path.splitext(filename)[0]
    
    # Get the most recently modified subdirectory that matches the needed substring from the inputfile
    latest_subdir = max(glob.glob(os.path.join(directory, f'{filename}*/')), key=os.path.getmtime) # FIXME: Can't use this method with multiprocessing - gives multiple fitnesses that are identical
    # latest_subdir = max(glob.glob(os.path.join(directory, '*/')), key=os.path.getmtime) # FIXME: Can't use this method with multiprocessing - gives multiple fitnesses that are identical
    
    meshfile = latest_subdir + "FEA_MeshNodes.csv"
    
    df = pd.read_csv(meshfile)
    return df

def getFitness(dfmesh):
    # abssums = dfmesh.abs().sum()
    absmeans = dfmesh.abs().mean()
    # stdevs = dfmesh.std()
    strain_xx = absmeans["strain_xx"]
    strain_yy = absmeans["strain_yy"]
    strain_xy = absmeans["strain_xy"]
    principalStrain_min = absmeans["principalStrain_min"]
    principalStrain_max = absmeans["principalStrain_max"]
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max
    

if __name__ == "__main__":
    
    # FEApath = "C:/Users/Ryan Larson/github/strain-gen-opt/FEA/FEA.exe"
    
    # # Choose FEA path here
    # initialdirFEA = "C:/Users/Ryan Larson/github/strain-gen-opt/FEA"
    # FEApath = chooseFEApath(initialdirFEA)
    
    # Load previously chosen FEA path here
    FEApath = loadFEApath('FEApath.pk')
    
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