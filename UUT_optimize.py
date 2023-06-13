# -*- coding: utf-8 -*-
"""
Created on Wed May 31 20:03:21 2023

@author: Ryan Larson
"""

import numpy as np
import pygad
import logging
from shapely.geometry import Point, Polygon
from shapely.prepared import prep
import matplotlib.pyplot as plt
import random

def fitness_func(ga_instance, solution, solution_idx):
    # Read in parameter values from parameter list
    
    model = translate_params_to_MATLAB(solution)
    
    strain_values = run_MATLAB_FEA(model)
    
    fitness = get_fitness(strain_values)
    
    return fitness

def translate_params_to_MATLAB(solution):
    # Turn generated solution values into inputs that can be fed to MATLAB FEA
    pass

def run_MATLAB_FEA(model):
    # return strain_values
    pass
    
def get_fitness(strain_values):
    # return fitness given strain (could be maximin fitness)
    pass
    
def on_generation(ga_instance):
    ga_instance.logger.info("Generation = {generation}".format(generation=ga_instance.generations_completed))
    ga_instance.logger.info("Fitness    = {fitness}".format(fitness=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]))
    




        
        
        
    
    
    
if __name__ == "__main__":
    print("No main function defined yet")