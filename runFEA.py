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
import re
from datetime import datetime
import numpy as np

def chooseFEApath(initialdir):
    """
    Saves the path to FEA.exe for convenient use later. This should be run
    during setup of the optimization codebase, before the first optimization
    run.

    Parameters
    ----------
    initialdir : str
        Path to the initial directory where the FEA.exe file is located (for
        convenience).

    Returns
    -------
    FEApath : path
        Full path to FEA executable.

    """
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
    """
    Loads the path to FEA.exe from a pickle file name.

    Parameters
    ----------
    filename : str
        Pickle file where FEA.exe path data is stored. Use 'FEApath.pk' as this
        input.

    Returns
    -------
    FEApath : path
        Full path to FEA.exe

    """
    with open(filename, 'rb') as fi:
        FEApath = pickle.load(fi)
    return FEApath

def runFEA(FEApath, inputfile):
    """
    Run FEA on a chosen inputfile.

    Parameters
    ----------
    FEApath : path
        Path to FEA.exe, loaded by loadFEApath().
    inputfile : XML filepath
        Filepath to XML design file that is to be optimized. This is chosen via
        Tkinter file dialog in the main optimization. The file can have
        pressure rods and standoffs in it already; those will be removed and
        replaced by the ones created by the optimization.

    Returns
    -------
    int
        Return code for FEA.exe. If it is anything other than 0, an error
        occurred during solving.

    """
    args = f'/input "{inputfile}" /noprogressbar'
    command = f'"{FEApath}" {args}'
    result = subprocess.run(command)
    if result.returncode == 0:
        # print(f"{inputfile} ran successfully")
        pass
    else:
        print(f"{inputfile} failed with code {result.returncode}")
    return result.returncode

def find_latest_folder_with_substring(base_dir, substring):
    """
    Checks for the most recently-created folder with a given substring. Used to
    choose the latest version of a given output folder, particularly if an
    optimization run has started without deleting data from past runs in
    base_dir.

    Parameters
    ----------
    base_dir : str
        Output folder where FEA solves are sent.
    substring : str
        String name for the current solve, including generation and iteration
        information

    Returns
    -------
    latest_folder : path
        Path to the most recently-created folder that has substring in its
        name.

    """
    latest_folder = None
    latest_timestamp = None

    # Iterate through the folders in the base directory
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)

        # Check if the folder name contains the specified substring
        if substring in folder_name:
            # Extract the timestamp from the folder name using regex
            timestamp_match = re.search(r'(\d{8}-\d{4})', folder_name)
            if timestamp_match:
                timestamp_str = timestamp_match.group(1)
                timestamp = datetime.strptime(timestamp_str, "%Y%m%d-%H%M")

                # Compare the timestamp with the latest one found
                if latest_timestamp is None or timestamp > latest_timestamp:
                    latest_timestamp = timestamp
                    latest_folder = folder_path

    return latest_folder

def resultsToDataframe_mesh(inputfile):
    """
    Get the output MeshNodes.csv data associated with the chosen input file
    and convert it to a Pandas dataframe for processing.

    Parameters
    ----------
    inputfile : str
        Path to the chosen input file.

    Returns
    -------
    df : Pandas dataframe
        Dataframe of the FEA_MeshNodes.csv file associated with the chosen
        input file.

    """
    directory = pathlib.Path(inputfile)
    directory = str(directory.parent) + "/Output"
    path, filename = os.path.split(inputfile)
    filename = os.path.splitext(filename)[0]
    
    # Get the most recently modified subdirectory that matches the needed substring from the inputfile
    latest_subdir = find_latest_folder_with_substring(directory, filename)
    
    # meshfile = latest_subdir + "\\FEA_MeshNodes.csv"
    meshfile = latest_subdir + "\\FEA_AllMeshNodes.csv"
    
    df = pd.read_csv(meshfile)
    return df


def resultsToDataframe_report(inputfile):
    """
    Get the output FEAReport.csv data associated with the chosen input file and
    convert it to a Pandas dataframe for processing.

    Parameters
    ----------
    inputfile : str
        Path to the chosen input file.

    Returns
    -------
    df : Pandas dataframe
        Dataframe of the FEAReport.csv file associated with the chosen input
        file.

    """
    directory = pathlib.Path(inputfile)
    directory = str(directory.parent) + "/Output"
    path, filename = os.path.split(inputfile)
    filename = os.path.splitext(filename)[0]
    
    # Get the most recently modified subdirectory that matches the needed substring from the inputfile
    latest_subdir = find_latest_folder_with_substring(directory, filename)
    
    meshfile = latest_subdir + "\\FEAFilteredReport.csv"
    # meshfile = latest_subdir + "\\FEAReport.csv"
    
    df = pd.read_csv(meshfile)
    return df



def getFitness(dfmesh):
    """
    Get the design fitness parameters using the maximum absolute value of the
    mesh strain. MUST use data from FEA_MeshNodes.csv.

    Parameters
    ----------
    dfmesh : Pandas dataframe
        Dataframe of the FEA_Meshnodes.csv output.

    Returns
    -------
    strain_xx : float
        Maximum magnitude of strain in the x direction.
    strain_yy : float
        Maximum magnitude of strain in the y direction.
    strain_xy : float
        Maximum magnitude of shear strain.
    principalStrain_min : float
        Maximum magnitude of minimum principal strain.
    principalStrain_max : float
        Maximum magnitude of maximum principal strain.

    """
    absmax = dfmesh.abs().max()
    strain_xx = absmax["strain_xx"]
    strain_yy = absmax["strain_yy"]
    strain_xy = absmax["strain_xy"]
    
    principalStrain_min = absmax["principalStrain_min"]
    principalStrain_max = absmax["principalStrain_max"]
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max

def getFitness_meshmax(dfmesh):
    """
    Get the design fitness parameters using the maximum absolute value of the
    mesh strain. MUST use data from FEA_MeshNodes.csv.

    Parameters
    ----------
    dfmesh : Pandas dataframe
        Dataframe of the FEA_Meshnodes.csv output.

    Returns
    -------
    strain_xx : float
        Maximum magnitude of strain in the x direction.
    strain_yy : float
        Maximum magnitude of strain in the y direction.
    strain_xy : float
        Maximum magnitude of shear strain.
    principalStrain_min : float
        Maximum magnitude of minimum principal strain.
    principalStrain_max : float
        Maximum magnitude of maximum principal strain.

    """
    absmax = dfmesh.abs().max()
    strain_xx = absmax["strain_xx"]
    strain_yy = absmax["strain_yy"]
    strain_xy = absmax["strain_xy"]
    
    principalStrain_min = absmax["principalStrain_min"]
    principalStrain_max = absmax["principalStrain_max"]
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max

def getFitness_meshsum(dfmesh, npts):
    """
    Get the design fitness parameters using the maximum absolute value of the
    mesh strain. MUST use data from FEA_MeshNodes.csv.

    Parameters
    ----------
    dfmesh : Pandas dataframe
        Dataframe of the FEA_Meshnodes.csv output.

    Returns
    -------
    strain_xx : float
        Maximum magnitude of strain in the x direction.
    strain_yy : float
        Maximum magnitude of strain in the y direction.
    strain_xy : float
        Maximum magnitude of shear strain.
    principalStrain_min : float
        Maximum magnitude of minimum principal strain.
    principalStrain_max : float
        Maximum magnitude of maximum principal strain.

    """
    dfmesh_copy = dfmesh.copy()
    dfmesh_copy = dfmesh_copy.abs()
    
    strain_xx = dfmesh_copy.nlargest(npts, columns='strain_xx').sum()['strain_xx']
    strain_yy = dfmesh_copy.nlargest(npts, columns='strain_yy').sum()['strain_yy']
    strain_xy = dfmesh_copy.nlargest(npts, columns='strain_xy').sum()['strain_xy']
    principalStrain_min = dfmesh_copy.nlargest(npts, columns='principalStrain_min').sum()['principalStrain_min']
    principalStrain_max = dfmesh_copy.nlargest(npts, columns='principalStrain_max').sum()['principalStrain_max']
    # strain_xx = abssum["strain_xx"]
    # strain_yy = abssum["strain_yy"]
    # strain_xy = abssum["strain_xy"]
    
    # principalStrain_min = abssum["principalStrain_min"]
    # principalStrain_max = abssum["principalStrain_max"]
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max

def getFitness_meshmean(dfmesh):
    """
    Get the design fitness parameters using the mean absolute value of the
    mesh strain. MUST use data from FEA_MeshNodes.csv.

    Parameters
    ----------
    dfmesh : Pandas dataframe
        Dataframe of the FEA_Meshnodes.csv output.

    Returns
    -------
    strain_xx : float
        Mean magnitude of strain in the x direction.
    strain_yy : float
        Mean magnitude of strain in the y direction.
    strain_xy : float
        Mean magnitude of shear strain.
    principalStrain_min : float
        Mean magnitude of minimum principal strain.
    principalStrain_max : float
        Mean magnitude of maximum principal strain.

    """
    absmean = dfmesh.abs().mean()
    strain_xx = absmean["strain_xx"]
    strain_yy = absmean["strain_yy"]
    strain_xy = absmean["strain_xy"]
    
    principalStrain_min = absmean["principalStrain_min"]
    principalStrain_max = absmean["principalStrain_max"]
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max

def getFitness_meshmedian(dfmesh):
    """
    Get the design fitness parameters using the median absolute value of the
    mesh strain. MUST use data from FEA_MeshNodes.csv.

    Parameters
    ----------
    dfmesh : Pandas dataframe
        Dataframe of the FEA_Meshnodes.csv output.

    Returns
    -------
    strain_xx : float
        Mean magnitude of strain in the x direction.
    strain_yy : float
        Mean magnitude of strain in the y direction.
    strain_xy : float
        Mean magnitude of shear strain.
    principalStrain_min : float
        Mean magnitude of minimum principal strain.
    principalStrain_max : float
        Mean magnitude of maximum principal strain.

    """
    absmedian = dfmesh.abs().median()
    strain_xx = absmedian["strain_xx"]
    strain_yy = absmedian["strain_yy"]
    strain_xy = absmedian["strain_xy"]
    
    principalStrain_min = absmedian["principalStrain_min"]
    principalStrain_max = absmedian["principalStrain_max"]
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max


def getFitness_v2(dfreport):
    """
    Get the design fitness parameters using the maximum magnitude value found
    in FEAReport.csv.

    Parameters
    ----------
    dfreport : Pandas dataframe
        Dataframe of the FEAReport.csv output.

    Returns
    -------
    strain_xx : float
        Maximum magnitude of strain in the x direction.
    strain_yy : float
        Maximum magnitude of strain in the y direction.
    strain_xy : float
        Maximum magnitude of shear strain.
    principalStrain_min : float
        Maximum magnitude of minimum principal strain.
    principalStrain_max : float
        Maximum magnitude of maximum principal strain.

    """
    dfreport.set_index('Row', inplace=True)
    strain_xx = np.max([np.abs(dfreport.loc['horizontalStrain_max','Value']),
                        np.abs(dfreport.loc['horizontalStrain_min','Value'])])
    strain_yy = np.max([np.abs(dfreport.loc['verticalStrain_max','Value']),
                        np.abs(dfreport.loc['verticalStrain_min','Value'])])
    strain_xy = np.max([np.abs(dfreport.loc['shearStrain_max','Value']),
                        np.abs(dfreport.loc['shearStrain_min','Value'])])
    principalStrain_min = np.abs(dfreport.loc['principalStrain_max','Value'])
    principalStrain_max = np.abs(dfreport.loc['principalStrain_min','Value'])
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max

def getFitness_reportmax(dfreport):
    """
    Get the design fitness parameters using the maximum magnitude value found
    in FEAReport.csv.

    Parameters
    ----------
    dfreport : Pandas dataframe
        Dataframe of the FEAReport.csv output.

    Returns
    -------
    strain_xx : float
        Maximum magnitude of strain in the x direction.
    strain_yy : float
        Maximum magnitude of strain in the y direction.
    strain_xy : float
        Maximum magnitude of shear strain.
    principalStrain_min : float
        Maximum magnitude of minimum principal strain.
    principalStrain_max : float
        Maximum magnitude of maximum principal strain.

    """
    dfreport.set_index('Row', inplace=True)
    strain_xx = np.max([np.abs(dfreport.loc['horizontalStrain_max','Value']),
                        np.abs(dfreport.loc['horizontalStrain_min','Value'])])
    strain_yy = np.max([np.abs(dfreport.loc['verticalStrain_max','Value']),
                        np.abs(dfreport.loc['verticalStrain_min','Value'])])
    strain_xy = np.max([np.abs(dfreport.loc['shearStrain_max','Value']),
                        np.abs(dfreport.loc['shearStrain_min','Value'])])
    principalStrain_min = np.abs(dfreport.loc['principalStrain_max','Value'])
    principalStrain_max = np.abs(dfreport.loc['principalStrain_min','Value'])
    return strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max
    

if __name__ == "__main__":    
    # # Choose FEA path here
    # initialdirFEA = "E:/github/strain-gen-opt/FEA"
    # FEApath = chooseFEApath(initialdirFEA)
    
    # Load previously chosen FEA path here
    FEApath = loadFEApath('FEApath.pk')
    
    # initialdir = str(pathlib.Path(FEApath).parent) + "Examples"
    
    # filetypes = (("XML", ["*.xml"]),)
    
    # root = Tk()
    # root.wm_attributes('-topmost', 1)
    # inputfile = fd.askopenfilename(
    #         title="Select FEA input file",
    #         initialdir=initialdir,
    #         filetypes=filetypes
    #         )
    # root.destroy()
    
    # runFEA(FEApath, inputfile)
    
    # dfmesh = resultsToDataframe(inputfile)
    
    # strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max = getFitness(dfmesh)