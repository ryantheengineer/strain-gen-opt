# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 22:14:37 2023

@author: Ryan Larson
"""

import random
import numpy as np
import matplotlib.pyplot as plt
import math
import constraints
import runFEA
import time
import multiprocessing
import copy

# MINIMIZATION

# On each iteration, out of 2 randomly selected parents we create 2 offsprings
# by taking fraction of genes from one parent and remaining fraction from other parent 
def crossover_prods(pop, crossover_rate, nprods, top_constraints):
    print(f"Performing crossover to create {crossover_rate} child designs")
    offspring = np.zeros((crossover_rate, pop.shape[1]))
    for i in range(crossover_rate):
        # Crossover with complete pressure rods.
        complete = False
        while complete is False:
            # Get parents
            r1 = np.random.randint(0, pop.shape[0])
            r2 = np.random.randint(0, pop.shape[0])
            while r1 == r2:
                r1 = np.random.randint(0, pop.shape[0])
                r2 = np.random.randint(0, pop.shape[0])
            parent1 = pop[r1]
            parent2 = pop[r2]
            
            # Interpret each parent into PressureRod representation
            parent1_prods = constraints.interpret_chromosome_to_prods(parent1, nprods)
            parent2_prods = constraints.interpret_chromosome_to_prods(parent2, nprods)
            
            # Perform crossover on PressureRod representation, with the understanding
            # that both parents were previously validated against all constraints,
            # so the only way crossed genes in a child can be invalid is if pressure
            # rods are too close together.
            child_prods = copy.deepcopy(parent1_prods)
            
            crossover_point = np.random.randint(1, nprods)
            
            for j in range(0, crossover_point):
                # Get gene at position j from parent 2
                child_prods[j] = copy.deepcopy(parent2_prods[j])
                
                # Check if the current potentially crossed gene violates any
                # other genes, and if so, replace that gene instead of gene j
                n_violated = 0
                idx_violated = []
                for k,prod in enumerate(child_prods):
                    if k == j:
                        continue
                    dist = constraints.centroid_distance(child_prods[j].tip, prod.tip)
                    if dist < prod.ctc:
                        n_violated += 1
                        idx_violated.append(k)
                if n_violated == 0:
                    continue
                elif n_violated == 1:
                    child_prods[j] = parent1_prods[j]
                    child_prods[idx_violated[0]] = parent2_prods[j]
                else:
                    continue
                    
            # Validate that at least one gene has been changed
            for j in range(len(child_prods)):
                if child_prods[j] != parent1_prods[j]:
                    complete = True
                    break
            
        # Interpret the child back to chromosome form
        child_chromosome = constraints.interpret_prods_to_chromosome(child_prods)
        offspring[i, :] = child_chromosome
        
        # # Plot the parent and child designs for examination
        # constraints.plot_prods_top_constraints(parent1_prods, top_constraints, f"Offspring {i}: Parent 1")
        # constraints.plot_prods_top_constraints(parent2_prods, top_constraints, f"Offspring {i}: Parent 2")
        # constraints.plot_prods_top_constraints(child_prods, top_constraints, f"Offspring {i}: Child")
        # plt.show()
        
    return offspring


# Perform crossover with mutation
def mutation(pop, n_mutated, mutation_rate, nprods, top_constraints):
    print("Entering mutation phase...creating children from crossover")
    offspring = np.zeros((n_mutated, pop.shape[1]))
    for i in range(n_mutated):
        # Crossover with complete pressure rods.
        complete = False
        while complete is False:
            # Get parents
            r1 = np.random.randint(0, pop.shape[0])
            r2 = np.random.randint(0, pop.shape[0])
            while r1 == r2:
                r1 = np.random.randint(0, pop.shape[0])
                r2 = np.random.randint(0, pop.shape[0])
            parent1 = pop[r1]
            parent2 = pop[r2]
            
            # Interpret each parent into PressureRod representation
            parent1_prods = constraints.interpret_chromosome_to_prods(parent1, nprods)
            parent2_prods = constraints.interpret_chromosome_to_prods(parent2, nprods)
            
            # Perform crossover on PressureRod representation, with the understanding
            # that both parents were previously validated against all constraints,
            # so the only way crossed genes in a child can be invalid is if pressure
            # rods are too close together.
            child_prods = copy.deepcopy(parent1_prods)
            
            crossover_point = np.random.randint(1, nprods)
            
            for j in range(0, crossover_point):
                # Get gene at position j from parent 2
                child_prods[j] = copy.deepcopy(parent2_prods[j])
                
                # Check if the current potentially crossed gene violates any
                # other genes, and if so, replace that gene instead of gene j
                n_violated = 0
                idx_violated = []
                for k,prod in enumerate(child_prods):
                    if k == j:
                        continue
                    dist = constraints.centroid_distance(child_prods[j].tip, prod.tip)
                    if dist < prod.ctc:
                        n_violated += 1
                        idx_violated.append(k)
                if n_violated == 0:
                    continue
                elif n_violated == 1:
                    child_prods[j] = parent1_prods[j]
                    child_prods[idx_violated[0]] = parent2_prods[j]
                else:
                    continue
            
            # Mutate prod form of child
            rod_types = ['Press-Fit Tapered',
                         'Press-Fit Flat',
                         '3.325" Tapered',
                         '3.325" Flat']
            pBoards_multi = top_constraints[0]
            top_probes = top_constraints[1]
            topcomponents = top_constraints[2]
            xmin, ymin, xmax, ymax = pBoards_multi.bounds
            # print("")
            for j in range(len(child_prods)):
                chance = random.uniform(0,1)
                if mutation_rate > chance:
                    # print("Mutating pressure rod gene")
                    while True:
                        valid = True
                        while True:
                            # Replace the pressure rod with a random new one (flipping on/off and changing sizes will be handled by local search)
                            x = random.uniform(xmin,xmax)
                            y = random.uniform(ymin,ymax)
                            rod_type_i = random.randint(0,3)
                            on = random.randint(0,1)
                            prod_new = constraints.PressureRod(x,y,rod_types[rod_type_i],on)
                            # Make sure pressure rod is within the UUT and make sure it doesn't intersect any components, using the appropriate buffer sizes
                            if not prod_new.tip_UUT_buffer.within(pBoards_multi):
                                valid = False
                                break
                            if prod_new.tip_component_buffer.intersects(topcomponents):
                                valid = False
                                break
                            
                            # Make sure the pressure rod isn't too close to any top probes
                            if prod_new.tip_from_top_probe_buffer.intersects(top_probes):
                                valid = False
                                break
    
                            # If the pressure rod is turned off, it's automatically valid
                            if prod_new.on == 0:
                                valid = True
                                break
                                
                            # Make sure pressure rod doesn't conflict with any previously-placed pressure rods
                            for k,prod_chosen in enumerate(child_prods):
                                if k == j:
                                    # Don't compare against the current pressure rod
                                    continue
                                dist = constraints.centroid_distance(prod_new.tip, prod_chosen.tip)
                                if dist < prod_new.ctc:
                                    valid = False
                                    break
                            if valid == False:
                                break
                            
                            if valid == True:
                                break
                            
                        if valid == True:
                            child_prods[j] == prod_new
                            # print("Mutation produced a new pressure rod")
                            break
                    
            
            # Validate that at least one gene has been changed
            for j in range(len(child_prods)):
                if child_prods[j] != parent1_prods[j]:
                    complete = True
                    break
            
        # Interpret the child back to chromosome form
        child_chromosome = constraints.interpret_prods_to_chromosome(child_prods)
        offspring[i, :] = child_chromosome
        
        # # Plot the parent and child designs for examination
        # constraints.plot_prods_top_constraints(parent1_prods, top_constraints, f"Offspring {i}: Parent 1")
        # constraints.plot_prods_top_constraints(parent2_prods, top_constraints, f"Offspring {i}: Parent 2")
        # constraints.plot_prods_top_constraints(child_prods, top_constraints, f"Offspring {i}: Child (mutated)")
        # plt.show()

    return offspring    # arr(mutation_size x n_var)

# Create some amount of offspring Q by adding fixed coordinate displacement to some 
# randomly selected parent's genes/coordinates
def local_search(pop, n_searched, localsearch_rate, fliprate, perturbrate, maxmag, typerate, nprods, top_constraints):
    print("Entering local search phase...creating altered versions of other chromosomes")
    offspring = np.zeros((n_searched, pop.shape[1]))
    for i in range(n_searched):
        # print(f"Local search {i} of {n_searched}")
        complete = False
        while complete is False:
            r1 = np.random.randint(0, pop.shape[0])
            parent1 = pop[r1]
            
            # Interpret each parent into PressureRod representation
            parent1_prods = constraints.interpret_chromosome_to_prods(parent1, nprods)
            child_prods = copy.deepcopy(parent1_prods)
            
            
            rod_types = ['Press-Fit Tapered',
                         'Press-Fit Flat',
                         '3.325" Tapered',
                         '3.325" Flat']
            pBoards_multi = top_constraints[0]
            top_probes = top_constraints[1]
            topcomponents = top_constraints[2]
            xmin, ymin, xmax, ymax = pBoards_multi.bounds
            for j in range(len(child_prods)):
                chance = random.uniform(0,1)
                if localsearch_rate > chance:
                    while True:
                        valid = True
                        while True:
                            # Chance to flip on/off parameter
                            old_on = child_prods[j].on
                            flip_chance = random.uniform(0,1)
                            if fliprate > flip_chance:
                                if old_on == 1:
                                    new_on = 0
                                else:
                                    new_on = 1
                                # print(f"Flipped prod {j} from {old_on} to {new_on}")
                                # print(f"\nChild prod {j} before:\t{child_prods[j].on}")
                                # print(f"Parent prod {j} before:\t{parent1_prods[j].on}")
                                child_prods[j].update_pressure_rod(child_prods[j].x, child_prods[j].y, child_prods[j].rod_type, new_on)
                                # print(f"Child prod {j} after:\t{child_prods[j].on}")
                                # print(f"Parent prod {j} after:\t{parent1_prods[j].on}")
                                
                                # Validate that the pressure rod doesn't conflict with other pressure rods
                                if child_prods[j].on == 0.0:
                                    valid = True
                                    break
                                # Make sure pressure rod doesn't conflict with any previously-placed pressure rods
                                for k, prod_chosen in enumerate(child_prods):
                                    if k == j:
                                        # Don't compare against the current pressure rod
                                        continue
                                    dist = constraints.centroid_distance(child_prods[j].tip, prod_chosen.tip)
                                    if dist < child_prods[j].ctc:
                                        # print(f"prod {j} was too close to another prod")
                                        valid = False
                                        break
                            else:
                                break
                                
                        if valid == False:
                            break
                        
                        if valid == True:
                            break
                        
                    # If the pressure rod is on, then it may be perturbed in other
                    # ways that will affect the design (position, type)
                    if child_prods[j].on:
                        chance = random.uniform(0,1)
                        if perturbrate > chance:
                            while True:
                                valid = True
                                while True:
                                    # Position perturb first
                                    mag = random.uniform(0, maxmag)
                                    xp = random.uniform(-mag, mag)
                                    yp = np.sqrt(mag**2 - xp**2)
                                    yp = random.choice([yp, -yp])
                                    
                                    # Adjust x and y by perturbation
                                    xnew = child_prods[j].x + xp
                                    ynew = child_prods[j].y + yp
                                    
                                    # print(f"\nChild prod {j} before:\t{child_prods[j].x}, {child_prods[j].y}")
                                    # print(f"Parent prod {j} before:\t{parent1_prods[j].x}, {parent1_prods[j].y}")                                    
                                    child_prods[j].update_pressure_rod(xnew, ynew, child_prods[j].rod_type, child_prods[j].on)
                                    # print(f"\nChild prod {j} after:\t{child_prods[j].x}, {child_prods[j].y}")
                                    # print(f"Parent prod {j} after:\t{parent1_prods[j].x}, {parent1_prods[j].y}")                                    
                                    
                                    
                                    # Validate the perturbation here
                                    # Make sure pressure rod is within the UUT and make sure it doesn't intersect any components, using the appropriate buffer sizes
                                    if not child_prods[j].tip_UUT_buffer.within(pBoards_multi):
                                        valid = False
                                        break
                                    if child_prods[j].tip_component_buffer.intersects(topcomponents):
                                        valid = False
                                        break
                                    
                                    # Make sure the pressure rod isn't too close to any top probes
                                    if child_prods[j].tip_from_top_probe_buffer.intersects(top_probes):
                                        valid = False
                                        break
                                        
                                    # Make sure pressure rod doesn't conflict with any previously-placed pressure rods
                                    for k,prod_chosen in enumerate(child_prods):
                                        if k == j:
                                            # Don't compare against the current pressure rod
                                            continue
                                        dist = constraints.centroid_distance(child_prods[j].tip, prod_chosen.tip)
                                        if dist < child_prods[j].ctc:
                                            valid = False
                                            break
                                        
                                else:
                                    break
                                    
                                if valid == False:
                                    break
                                
                                if valid == True:
                                    break
                        
                        
                        chance = random.uniform(0,1)
                        if typerate > chance:
                            # Get the index of the current rod_type in the rod_types list
                            current_type_idx = rod_types.index(child_prods[j].rod_type)
                            # Randomly shuffle the indices of the rod types that are not used
                            available_indices = [0,1,2,3]
                            available_indices.remove(current_type_idx)
                            random.shuffle(available_indices)
                            
                            for idx in available_indices:
                                valid = True
                                # print(f"\nChild prod {j} before:\t{child_prods[j].rod_type}")
                                # print(f"Child prod {j} before:\t{parent1_prods[j].rod_type}")
                                child_prods[j].update_pressure_rod(child_prods[j].x, child_prods[j].y, rod_types[idx], child_prods[j].on)
                                # print(f"Child prod {j} after:\t{child_prods[j].rod_type}")
                                # print(f"Child prod {j} after:\t{parent1_prods[j].rod_type}")
                                
                                # Validate the new prod type. If it is valid, break.
                                # Make sure pressure rod is within the UUT and make sure it doesn't intersect any components, using the appropriate buffer sizes
                                if not child_prods[j].tip_UUT_buffer.within(pBoards_multi):
                                    valid = False
                                    continue
                                if child_prods[j].tip_component_buffer.intersects(topcomponents):
                                    valid = False
                                    continue
                                
                                # Make sure the pressure rod isn't too close to any top probes
                                if child_prods[j].tip_from_top_probe_buffer.intersects(top_probes):
                                    valid = False
                                    continue
                                    
                                # Make sure pressure rod doesn't conflict with any previously-placed pressure rods
                                for k,prod_chosen in enumerate(child_prods):
                                    if k == j:
                                        # Don't compare against the current pressure rod
                                        continue
                                    dist = constraints.centroid_distance(child_prods[j].tip, prod_chosen.tip)
                                    if dist < child_prods[j].ctc:
                                        valid = False
                                        break
                            
                                if valid:
                                    break
                            
                            # If no other rod type is valid, put it back the way it was
                            if valid == False:
                                # print("Had to use original prod type")
                                child_prods[j].update_pressure_rod(child_prods[j].x, child_prods[j].y, rod_types[current_type_idx], child_prods[j].on)
                                
            # Validate that at least one gene has been changed
            for j in range(len(child_prods)):
                if child_prods[j] != parent1_prods[j]:
                    complete = True
                    # print(f"Prod {j} has been changed")
                    # print("At least one gene changed via local search")
                    break
                else:
                    pass
                    # print("\nNO PRESSURE RODS WERE CHANGED COMPARED TO THE PARENT DESIGN\n")

        # Interpret the child back to chromosome form
        child_chromosome = constraints.interpret_prods_to_chromosome(child_prods)
        offspring[i, :] = child_chromosome
    return offspring    # arr(loc_search_size x n_var)

# Calculate fitness (obj function) values for each chromosome/solution
def evaluation(pop, nobjs, gen, nprods, inputfile, constraint_geom):
    # fitness_values = np.zeros((pop.shape[0], nobjs))
    
    # Read in constraint_geom (output of constraints.get_constraint_geometry())
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
    
    ncpus = multiprocessing.cpu_count()
    
    prods = []
    valid_circles = []
    for i, chromosome in enumerate(pop):
        prods.append(constraints.interpret_chromosome_to_prods(chromosome, nprods))
        valid_circles.append(constraints.prods_to_valid_circles(prods[i]))
        
    # Generate the XML files sequentially
    xml_filenames = []
    for i, valid_design in enumerate(valid_circles):
        xml_filenames.append(constraints.design_to_xml(valid_design, df_PressureRods, root, inputfile, gen, i))
        
    pool = multiprocessing.Pool(processes=ncpus)
    # arg_tuples = [(xml_filenames[i]) for i in range(pop.shape[0])]
    arg_tuples = [(valid_circles[i], df_PressureRods, root, inputfile, gen, i) for i in range(pop.shape[0])]
    # arg_tuples = [(valid_circles[i], df_PressureRods, root, xml_filenames[i], gen, i) for i in range(pop.shape[0])]
    # results_mp = pool.starmap(constraints.runFEA_new_path, arg_tuples)  # FIXME: runFEA_valid_circles is designed to reuse FEA.xml, not write a new version to include the new parameters. This probably causes issues with running multiple at the same time.
    results_mp = pool.starmap(constraints.runFEA_valid_circles, arg_tuples)  # FIXME: runFEA_valid_circles is designed to reuse FEA.xml, not write a new version to include the new parameters. This probably causes issues with running multiple at the same time.
    pool.close()
    pool.join()
    
    fitness_values = np.array(results_mp)
    # fitness_values = np.asarray(list(zip(*results_mp)))
        
        # # FIXME: Need code to translate chromosomes into designs, then run FEA
        # prods = constraints.interpret_chromosome_to_prods(chromosome, nprods)
        # valid_circles = constraints.prods_to_valid_circles(prods)
        # strain_xx, strain_yy, strain_xy, principalStrain_min, principalStrain_max = constraints.runFEA_valid_circles(valid_circles, df_PressureRods, root, inputfile, gen, i)
        
        # fitness_values[i,0] = strain_xx
        # fitness_values[i,1] = strain_yy
        # fitness_values[i,2] = strain_xy
        # fitness_values[i,3] = principalStrain_min + principalStrain_max

    return fitness_values


# Estimate how tightly clumped fitness values are on Pareto front. 
def crowding_calculation(fitness_values):
    pop_size = len(fitness_values[:, 0])
    fitness_value_number = len(fitness_values[0, :])                    # == n of objective functions
    matrix_for_crowding = np.zeros((pop_size, fitness_value_number))    # arr(pop_size x 2) 
    normalized_fitness_values = (fitness_values - fitness_values.min(0))/fitness_values.ptp(0)  # arr.ptp(0) array of max elem in each col
    
    for i in range(fitness_value_number):
        crowding_results = np.zeros(pop_size)
        crowding_results[0] = 1 # extreme point has the max crowding distance
        crowding_results[pop_size - 1] = 1 # extreme point has the max crowding distance
        sorted_normalized_fitness_values = np.sort(normalized_fitness_values[:,i])
        sorted_normalized_values_index = np.argsort(normalized_fitness_values[:,i])
        # crowding distance calculation. Say for fitness1[i], crowding = fitness1[i+1] - fitness1[i-1]
        crowding_results[1:pop_size - 1] = (sorted_normalized_fitness_values[2:pop_size] - sorted_normalized_fitness_values[0:pop_size - 2])
        re_sorting = np.argsort(sorted_normalized_values_index)
        matrix_for_crowding[:, i] = crowding_results[re_sorting]
    
    crowding_distance = np.sum(matrix_for_crowding, axis=1) # on fitness1 - fitness2 plot, each point on pareto front has crowding distance number

    return crowding_distance    # arr(pop_size,)

# Crowding distance is used to maintain diversity of solutions on Pareto front. 
# Remove some amount of solutions that are clumped together to much
def remove_using_crowding(fitness_values, number_solutions_needed):
    pop_index = np.arange(fitness_values.shape[0])
    crowding_distance = crowding_calculation(fitness_values)
    selected_pop_index = np.zeros(number_solutions_needed)
    selected_fitness_values = np.zeros((number_solutions_needed, len(fitness_values[0, :])))    # arr(num_sol_needed x 2)
    for i in range(number_solutions_needed):
        pop_size = pop_index.shape[0]
        solution_1 = random.randint(0, pop_size - 1)
        solution_2 = random.randint(0, pop_size - 1)
        if crowding_distance[solution_1] >= crowding_distance[solution_2]:
            # solution 1 is better than solution 2
            selected_pop_index[i] = pop_index[solution_1]
            selected_fitness_values[i, :] = fitness_values[solution_1, :]
            pop_index = np.delete(pop_index, (solution_1), axis=0)
            fitness_values = np.delete(fitness_values, (solution_1), axis=0)
            crowding_distance = np.delete(crowding_distance, (solution_1), axis=0)
        else:
            # solution 2 is better than solution 1
            selected_pop_index[i] = pop_index[solution_2]
            selected_fitness_values[i, :] = fitness_values[solution_2, :]
            pop_index = np.delete(pop_index, (solution_2), axis=0)
            fitness_values = np.delete(fitness_values, (solution_2), axis=0)
            crowding_distance = np.delete(crowding_distance, (solution_2), axis=0)
    
    selected_pop_index = np.asarray(selected_pop_index, dtype=int)

    return selected_pop_index   # arr(n_sol_needed,)

# find indices of solutions that dominate others
def pareto_front_finding(fitness_values, pop_index):
    pop_size = fitness_values.shape[0]
    pareto_front = np.ones(pop_size, dtype=bool)    # all True initially
    for i in range(pop_size):
        for j in range(pop_size):
            if all(fitness_values[j] <= fitness_values[i]) and any(fitness_values[j] < fitness_values[i]):
                pareto_front[i] = 0 # i is not in pareto front becouse j dominates i
                break

    return pop_index[pareto_front]  # arr(len_pareto_front,)

# repeat Pareto front selection to build a population within defined size limits
def selection(pop, fitness_values, pop_size):
    
    pop_index_0 = np.arange(pop.shape[0])   # unselected pop ids
    pop_index = np.arange(pop.shape[0])     # all pop ids. len = len(pop_size)
    pareto_front_index = []
    
    while len(pareto_front_index) < pop_size:   # pop_size = initial_pop_size
        new_pareto_front = pareto_front_finding(fitness_values[pop_index_0, :], pop_index_0)
        total_pareto_size = len(pareto_front_index) + len(new_pareto_front)

        # check the size of pareto_front, if larger than pop_size then remove some
        if total_pareto_size > pop_size:
            number_solutions_needed = pop_size - len(pareto_front_index)
            selected_solutions = remove_using_crowding(fitness_values[new_pareto_front], number_solutions_needed)
            new_pareto_front = new_pareto_front[selected_solutions]
        
        pareto_front_index = np.hstack((pareto_front_index, new_pareto_front))
        remaining_index = set(pop_index) - set(pareto_front_index)
        pop_index_0 = np.array(list(remaining_index))
        
    selected_pop = pop[pareto_front_index.astype(int)]

    return selected_pop     # arr(pop_size x n_var)

def main_optimization():
    start_time = time.time()
    # Initial setup
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
    
    # Estimate a number of pressure rods for the top side that would make sense
    # print("--- Estimating possible pressure rods ---")
    nprods_small, nprods_large, pBoards_diff = constraints.grid_nprods(pBoards,pComponentsTop)
    
    top_constraints = constraints.get_top_constraints(pBoards, pComponentsTop, df_Probes, pBoards_diff)
    
    # Parameters
    # n_var = 3                   # chromosome has 3 coordinates/genes
    # lb = [-5, -5, -5]
    # ub = [5, 5, 5]
    print("Setting genetic algorithm parameters")
    pop_size = 40              # initial number of chromosomes
    rate_crossover = 10         # number of chromosomes that we apply crossover to
    rate_mutation = 10          # number of chromosomes that we apply mutation to
    chance_mutation = 0.3       # normalized percent chance that an individual pressure rod will be mutated
    n_searched = 10      # number of chromosomes that we apply local_search to
    chance_localsearch = 0.5
    fliprate = 0.3
    perturbrate = 1.0
    maxmag = 0.1             # coordinate displacement during local_search
    typerate = 0.3
    maximum_generation = 10    # number of iterations
    nobjs = 4
    
    # nprods = 64
    nprods = 30
    # nprods = np.max([len(df_PressureRods), nprods_small, nprods_large])
    print("Creating initial random population")
    pop = constraints.initialize_population_simple(pop_size, nprods, pBoards, pComponentsTop, df_Probes, pBoards_diff)    # initial parents population P
    pop = np.asarray(pop)
    end_setup_time = time.time()
    
    best_fitnesses_1 = []
    best_fitnesses_2 = []
    best_fitnesses_3 = []
    best_fitnesses_4 = []
    # NSGA-II main loop
    for i in range(maximum_generation):
        offspring_from_crossover = crossover_prods(pop, rate_crossover, nprods, top_constraints)
        offspring_from_mutation = mutation(pop, rate_mutation, chance_mutation, nprods, top_constraints)
        offspring_from_local_search = local_search(pop, n_searched, chance_localsearch, fliprate, perturbrate, maxmag, typerate, nprods, top_constraints)
        
        # Append children (crossover, mutation, local search) to parents
        pop = np.append(pop, offspring_from_crossover, axis=0)
        pop = np.append(pop, offspring_from_mutation, axis=0)
        pop = np.append(pop, offspring_from_local_search, axis=0)
        # print(pop.shape)
        print("Evaluating fitnesses...")
        fitness_values = evaluation(pop, nobjs, i, nprods, inputfile, constraint_geom)
        j = fitness_values[:,0].argmin()
        best_fitnesses_1.append(fitness_values[j,:])
        j = fitness_values[:,1].argmin()
        best_fitnesses_2.append(fitness_values[j,:])
        j = fitness_values[:,2].argmin()
        best_fitnesses_3.append(fitness_values[j,:])
        j = fitness_values[:,3].argmin()
        best_fitnesses_4.append(fitness_values[j,:])
        pop = selection(pop, fitness_values, pop_size)  # we arbitrary set desired pereto front size = pop_size
        print('Generation:', i)
        fig,ax = plt.subplots(dpi=300)
        for j in range(len(pop)):
            x1 = fitness_values[j][0]
            x2 = fitness_values[j][1]
            # x1 = pop[j][0]
            # x2 = pop[j][1]
            ax.scatter(x1,x2,marker='o',color='b')
        ax.set_xlabel('x1')
        ax.set_ylabel('x2')
        ax.set_title(f"Generation: {i}")
        plt.show()
        # fig = plt.figure(dpi=300)
        # ax = fig.add_subplot(projection='3d')
        # for j in range(len(pop)):
        #     x1 = pop[j][0]
        #     x2 = pop[j][1]
        #     x3 = pop[j][2]
        #     ax.scatter(x1,x2,x3, marker='o', color='b')
        # ax.set_xlabel('x1')
        # ax.set_ylabel('x2')
        # ax.set_zlabel('x3')
        # ax.set_xlim3d(-5, 5)
        # ax.set_ylim3d(-5, 5)
        # ax.set_zlim3d(-5, 5)
        # ax.set_title(f"Generation: {i}")
    
    # Pareto front visualization
    fitness_values = evaluation(pop, nobjs, i, nprods, inputfile, constraint_geom)
    index = np.arange(pop.shape[0]).astype(int)
    pareto_front_index = pareto_front_finding(fitness_values, index)
    pop = pop[pareto_front_index, :]
    # print("_________________")
    # print("Optimal solutions:")
    # print("       x1               x2                 x3")
    # print(pop) # show optimal solutions
    fitness_values = fitness_values[pareto_front_index]
    # print("______________")
    # print("Fitness values:")
    # print("  objective 1    objective 2")
    # print(fitness_values)
    best_fitnesses_1 = np.asarray(best_fitnesses_1)
    best_fitnesses_2 = np.asarray(best_fitnesses_2)
    best_fitnesses_3 = np.asarray(best_fitnesses_3)
    best_fitnesses_4 = np.asarray(best_fitnesses_4)
    plt.figure(dpi=300)
    plt.scatter(fitness_values[:, 0],fitness_values[:, 1], label='Pareto optimal front')
    plt.scatter(best_fitnesses_1[:,0],best_fitnesses_1[:,1], label="Optimal objective 1")
    plt.scatter(best_fitnesses_2[:,0],best_fitnesses_2[:,1], label="Optimal objective 2")
    plt.legend(loc='best')
    plt.xlabel('Objective function F1')
    plt.ylabel('Objective function F2')
    plt.title('Optimal designs over all generations')
    # plt.grid(b=1)
    plt.show()
    
    end_time = time.time()
    
    print(f"\n\nSetup time:\t{end_setup_time-start_time}")
    print(f"Total elapsed time:\t{end_time-start_time}")
    
    plt.figure(dpi=300)
    plt.plot(best_fitnesses_1[:,0], label="Avg strain xx")
    plt.plot(best_fitnesses_2[:,1], label="Avg strain yy")
    plt.plot(best_fitnesses_3[:,2], label="Avg strain xy")
    plt.plot(best_fitnesses_4[:,3], label="Sum principal strains")
    plt.title("Fitnesses by objective")
    plt.legend()
    
    
    best_fitnesses = np.concatenate((best_fitnesses_1, best_fitnesses_2, best_fitnesses_3, best_fitnesses_4), axis=1)
    
    return fitness_values, best_fitnesses, pop

if __name__ == "__main__":
    fitness_values, best_fitnesses, pop = main_optimization()