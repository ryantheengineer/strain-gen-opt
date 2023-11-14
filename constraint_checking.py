# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 20:48:29 2023

@author: Ryan Larson
"""

import constraints
import optimize
from shapely.geometry import Point, Polygon, MultiPolygon
import numpy as np
import pickle

if __name__ == "__main__":
    reload = True
    filename = 'Example_1_40prods_pop20_5children_per.pkl'

    # Parameters
    print("Setting genetic algorithm parameters")
    pop_size = 20              # initial number of chromosomes
    rate_crossover = 5         # number of chromosomes that we apply crossover to
    rate_mutation = 5          # number of chromosomes that we apply mutation to
    chance_mutation = 0.3       # normalized percent chance that an individual pressure rod will be mutated
    n_searched = 5      # number of chromosomes that we apply local_search to
    chance_localsearch = 0.5
    fliprate = 0.3
    perturbrate = 1.0
    maxmag = 0.1             # coordinate displacement during local_search
    typerate = 0.3
    maximum_generation = 15    # number of iterations
    nobjs = 4
    
    nprods = 40
    
    if not reload:
        print("Initial setup - reading in constraints")
        constraint_geom = constraints.get_constraint_geometry()
        root = constraint_geom[0]
        inputfile = constraint_geom[1]
        pBoards = constraint_geom[2]
        pOutline = constraint_geom[3]
        pShape = constraint_geom[4]
        pComponentsTop = constraint_geom[5]
        pComponentsBot = constraint_geom[6]
        Pressure = constraint_geom[7]
        I_Plate = constraint_geom[8]
        Stripper = constraint_geom[9]
        Probe = constraint_geom[10]
        Countersink = constraint_geom[11]
        df_Probes = constraint_geom[12]
        df_GuidePins = constraint_geom[13]
        df_PressureRods = constraint_geom[14]
        df_Standoffs = constraint_geom[15]
        
        nprods_small, nprods_large, pBoards_diff = constraints.grid_nprods(pBoards,pComponentsTop)
        
        top_constraints = constraints.get_top_constraints(pBoards, pComponentsTop, df_Probes, pBoards_diff)
        
        # Generate population
        pop = constraints.initialize_population_simple(pop_size, nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff)    # initial parents population P
        pop = np.asarray(pop)
        
        # Pickle the variables
        with open(filename, 'wb') as f:
            pickle.dump([pop, top_constraints], f)
    
    else:
        # Load the pickled variables
        with open(filename, 'rb') as f:
            pop, top_constraints = pickle.load(f)
        
    pBoards_multi = top_constraints[0]
    top_probes = top_constraints[1]
    topcomponents = top_constraints[2]
    
    # Test offspring using different methods
    offspring_from_crossover = optimize.crossover_prods(pop, rate_crossover, nprods, top_constraints)
    offspring_from_mutation = optimize.mutation(pop, rate_mutation, chance_mutation, nprods, top_constraints)
    offspring_from_local_search = optimize.local_search(pop, n_searched, chance_localsearch, fliprate, perturbrate, maxmag, typerate, nprods, top_constraints)