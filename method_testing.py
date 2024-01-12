# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 21:13:50 2024

@author: Ryan Larson

Method testing for optimization to make it easier to see where things could
be going wrong

"""
import numpy as np
import constraints
import runFEA
import optimize
import itertools

def crossover_test(pop, crossover_rate, nprods_top, nprods_bot, top_constraints, bot_constraints):
    offspring = optimize.crossover_prods(pop, crossover_rate, nprods_top, nprods_bot, top_constraints, bot_constraints)
    
    offspring_shape = offspring.shape
    offspring_prods = []
    for i in range(offspring.shape[0]):
        print(f"Checking offspring {i}")
        prods = constraints.interpret_chromosome_to_prods_v2(offspring[i], nprods_top, nprods_bot)
        offspring_prods.append(prods)
        
        # Compare all combinations of pressure rods to see if there are 
        # any distances less than the .ctc parameter
        combos = []
        for r in range(len(prods)):
            combos.extend(itertools.combinations(prods,2))
        
        printflag = False
        for combo in combos:
            dist = constraints.centroid_distance(combo[0].tip, combo[1].tip)
            if dist < combo[0].ctc or dist < combo[1].ctc:
                print(f"Pressure rods found that are too close together - offspring {i}")
                printflag = True
                break
        
        if printflag is True:
            title = f"Offspring {i}"
            constraints.plot_prods_top_constraints(prods, top_constraints, title)
    
    
    # Plot all offspring and save
    
    return offspring_prods

def initial_pop_test(pop, nprods_top, nprods_bot, top_constraints, bot_constraints):
    # Get prod versions of pop
    pop_prods = []
    for i in range(pop.shape[0]):
        print(f"Checking initial design {i}")
        prods = constraints.interpret_chromosome_to_prods_v2(pop[i], nprods_top, nprods_bot)
        pop_prods.append(prods)
        
        # Compare all combinations of pressure rods to see if there are 
        # any distances less than the .ctc parameter
        combos = []
        for r in range(len(prods)):
            combos.extend(itertools.combinations(prods,2))
        
        printflag = False
        for combo in combos:
            dist = constraints.centroid_distance(combo[0].tip, combo[1].tip)
            if dist < combo[0].ctc or dist < combo[1].ctc:
                print(f"Pressure rods found that are too close together - initial design {i}")
                printflag = True
                break
        
        if printflag is True:
            title = f"Design {i}"
            constraints.plot_prods_top_constraints(prods, top_constraints, title)
        
    
    
if __name__ == "__main__":
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
    
    # Estimate a number of pressure rods for the top side that would make sense
    # print("--- Estimating possible pressure rods ---")
    nprods_small, nprods_large, pBoards_diff = constraints.grid_nprods(pBoards,pComponentsTop) # FIXME: Need to make this consider bottom side pressure rods
    
    # top_constraints = constraints.get_top_constraints(pBoards, pComponentsTop, df_Probes, pBoards_diff)
    top_constraints = constraints.get_board_constraints_single_side(pBoards, pComponentsTop, 1, df_Probes, pBoards_diff)
    bot_constraints = constraints.get_board_constraints_single_side(pBoards, pComponentsBot, 2, df_Probes, pBoards_diff)
    
    # Parameters
    print("Setting genetic algorithm parameters")
    pop_size = 20              # initial number of chromosomes
    rate_crossover = 150         # number of chromosomes that we apply crossover to
    rate_mutation = 9         # number of chromosomes that we apply mutation to
    chance_mutation = 0.2       # normalized percent chance that an individual pressure rod will be mutated
    n_searched = 9              # number of chromosomes that we apply local_search to
    chance_localsearch = 0.2
    on_prob_initial = 0.5   # Initial percentage chance that a pressure rod will be on (only in the initial population)
    on_prob = 0.8           # Likelihood an "off" pressure rod will be switched on
    perturbrate = 1.0
    maxmag = 0.1             # coordinate displacement during local_search
    typerate = 0.1
    maximum_generation = 12    # number of iterations
    nobjs = 5
    
    end_early = True
    # FIXME: Add ability to pickle the variables needed to continue an optimization later
    
    nprods = 64
    nprods_top = 64
    nprods_bot = 0
    print(f"nprods_small = {nprods_small}")
    print(f"nprods_large = {nprods_large}")
    nprods_top_input = input(f"Current nprods_top: {nprods_top}\n If this quantity is adequate press enter. Otherwise choose an integer value and press enter.\n")
    if len(nprods_top_input) == 0:
        print(f"\nCurrent value of {nprods_top} accepted.")
        pass
    else:
        nprods_top = int(nprods_top_input)
        print(f"New value of {nprods_top} accepted.")
    # nprods = 40
    # nprods = nprods_small
    # nprods = len(df_PressureRods)
    # nprods = np.max([len(df_PressureRods), nprods_small, nprods_large])
    
    design_accepted = False     # Flag for deciding whether to end optimization early if criteria are met
    
    print("Creating initial random population")
    all_on = True
    rod_type = 'Press-Fit Tapered'
    # rod_types = ['Press-Fit Tapered',
    #              'Press-Fit Flat',
    #              '3.325" Tapered',
    #              '3.325" Flat']
    pop = constraints.initialize_population_simple_v2(pop_size, nprods_top, nprods_bot, top_constraints, bot_constraints, all_on, on_prob, rod_type)    # initial parents population P
    # pop = constraints.initialize_population_simple(pop_size, nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff, all_on, on_prob_initial, rod_type)    # initial parents population P
    pop = np.asarray(pop)
    
    
    
    ##### TESTING BEGINS HERE #####
    # Testing initial population for pressure rods that are getting placed
    # too closely together
    print("\n\nTesting initial population")
    initial_pop_test(pop, nprods_top, nprods_bot, top_constraints, bot_constraints)
    
    
    # Testing crossover for pressure rods that end up too close to a component
    # or to another pressure rod
    print("\n\nTesting crossover")
    crossover_test(pop, rate_crossover, nprods_top, nprods_bot, top_constraints, bot_constraints)