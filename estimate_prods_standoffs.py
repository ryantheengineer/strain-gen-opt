# -*- coding: utf-8 -*-
"""
Created on Mon May 20 21:34:03 2024

@author: ryanl
"""
import constraints
import optimize
import numpy as np
import itertools
import multiprocessing
from datetime import datetime
import time
import matplotlib.pyplot as plt
import copy

def objective(xs, top_constraints, bot_constraints, all_on, on_prob, rod_type):
    # Create nsims random chromosomes and evaluate them. Sum the first three
    # objective values to get the overall fitness, then average the results
    
    nprods = xs[0]
    nstandoffs = xs[1]
    nsims = 3
    objval = 0
    gen = 'setup'
    maxgen = 100
    
    chromosomes = constraints.initialize_population_simple_v3(nsims, nprods, nstandoffs, top_constraints, bot_constraints, all_on, on_prob, rod_type)
    pop = np.asarray(chromosomes)
    
    fitness_values, fitness_report, fitness_mesh = optimize.evaluation(pop, gen, maxgen, nprods, nstandoffs, inputfile, constraint_geom, all_on, on_prob, rod_type)
    
    fitness_values = np.asarray(fitness_report)
    fitness_values = fitness_values[:,0:3]
    objval += np.sum(fitness_values)
    objval /= nsims
    
    return objval
    
def Boltzmann(dE, dEavg, T):
    P = np.exp(-dE/(dEavg*T))
    return P

def SimAnneal(nprods, nprods_bounds, nstandoffs, nstandoffs_bounds, constraint_geom, all_on, on_prob, rod_type):
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
    pBoards_diff_top = constraint_geom[16]
    pBoards_diff_bot = constraint_geom[17]
    top_constraints = constraint_geom[18]
    bot_constraints = constraint_geom[19]
    
    # Starting design
    xs = np.array([nprods, nstandoffs])
    xub = np.array([nprods_bounds[0], nstandoffs_bounds[0]])
    xlb = np.array([nprods_bounds[1], nstandoffs_bounds[1]])
    for i in range(len(xub)):
        if xs[i] > xub[i]:
            raise ValueError("Out of bounds")
        elif xs[i] < xlb[i]:
            raise ValueError("Out of bounds")
    
    
    fs = objective(xs, top_constraints, bot_constraints, all_on, on_prob, rod_type)
    xsearch = [xs]

    # Select Ps, Pf, N, and calculate Ts, Tf, and F
    Ps = 0.3                # Probability of acceptance at start
    Pf = 0.00001             # Probability of acceptance at finish
    N = 7                 # Number of cycles

    Ts = -1/np.log(Ps)      # Temperature at start
    Tf = -1/np.log(Pf)      # Temperature at finish
    F = (Tf/Ts)**(1/(N-1))  # Temperature reduction factor each cycle

    # Perturbation information
    delta = 20               # Max perturbation
    n = 4                   # Starting number of perturbations per cycle

    # Holding variables
    dE = 0.0
    dEavg = 0.0
    perturbations = list(range(N))
    objvals = [0] * N

    # Set starting values
    xc = xs
    fc = fs
    T = Ts

    # Step through the cycles
    for i in range(N):
        print(f'\nCycle {i}')
        # Add the current objective value to the objective vector for plotting
        print(f'Calculating objective at current position:\t{xc}')
        objvals[i] = objective(xc, top_constraints, bot_constraints, all_on, on_prob, rod_type)

        # Step through the perturbations
        for j in range(n):
            # Perturb xc by some random value within delta. If any perturbed
            # value falls outside the specified bounds, retry the perturbation.
            while True:
                p1 = np.random.uniform(-delta, delta)
                p2 = np.random.uniform(-delta, delta)
                # Sp = np.random.uniform(-delta, delta)
                # xpp = np.random.uniform(-delta, delta)
                # ypp = np.random.uniform(-delta, delta)
                perturb = np.array([p1,p2])
                xp = xc + perturb
                xp = xp.astype(int)

                if xp[0] > xub[0]:
                    continue
                elif xp[0] < xlb[0]:
                    continue
                else:
                    break

            # print(xp)

            # Get the objective value at the perturbed point
            print(f'Calculating objective at perturbed position:\t{xp}')
            fp = objective(xp, top_constraints, bot_constraints, all_on, on_prob, rod_type)

            # Calculate values for Boltzmann function in case they're needed
            dE = np.abs(fp - fc)
            if i == 1 and j == 1:
                dEavg = dE
            else:
                dEavg = (dEavg + dE)/2

            P = Boltzmann(dE, dEavg, T)

            # Check if the new design is better than the old design
            if fp < fc:
                xc = copy.deepcopy(xp)     # Accept as current design if better
                fc = copy.deepcopy(fp)
            else:
                # If the new design is worse, generate a random number and
                # compare to the Boltzmann probability. If the random number is
                # lower than the Boltzmann probability, accept the worse design
                # as the current design
                randnum = np.random.uniform(0,1)

                if randnum < P:
                    xc = copy.deepcopy(xp)
                    fc = copy.deepcopy(fp)

        # Decrease the temperature by factor F
        T = F*T

        # Increase the number of perturbations every few cycles
        if (i % 3) == 1:
            n += 1

        # Save the new search position at the end of each cycle
        xsearch.append(xc)

    return perturbations, objvals, xsearch

if __name__ == '__main__':
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
    pBoards_diff_top = constraint_geom[16]
    pBoards_diff_bot = constraint_geom[17]
    top_constraints = constraint_geom[18]
    bot_constraints = constraint_geom[19]
    
    all_on = True
    rod_type = 'Press-Fit Tapered'
    on_prob = 0.8
    
    start_time = datetime.now()
    
    # # Estimate a number of pressure rods for the top side that would make sense
    adjust_factor = 0.5
    nprods_small, nprods_large = constraints.grid_nprods_v2(pBoards_diff_top)
    nstandoffs_small, nstandoffs_large = constraints.grid_nprods_v2(pBoards_diff_bot)
    nprods_bounds = (int(adjust_factor*nprods_small), 1)
    # nstandoffs_bounds = (int(adjust_factor*nstandoffs_small), 1)
    nstandoffs_bounds = copy.deepcopy(nprods_bounds)
    
    nprods = int(np.mean(nprods_bounds))
    nstandoffs = int(np.mean(nstandoffs_bounds))
    
    print(f'Pressure rods:\t{nprods}')
    print(f'Standoffs:\t{nstandoffs}')
    
    perturbations, objvals, xsearch = SimAnneal(nprods, nprods_bounds, nstandoffs, nstandoffs_bounds, constraint_geom, all_on, on_prob, rod_type)
    
    
    # perturbations, objvals, xsearch = SimAnneal(2000.0)

    x_end = xsearch[-1][:]

    best_ind = objvals.index(min(objvals))
    x_best = xsearch[best_ind][:]
    # S_best = xsearch[best_ind][1]
    # wp_best = xsearch[best_ind][2]
    # hp_best = xsearch[best_ind][3]

    # Plot the cooling curve
    fig1 = plt.figure(1, figsize=(12,8))
    plt.plot(perturbations, objvals)
    plt.xlabel("Cycles")
    plt.ylabel("Objective")
    plt.title("Cooling Curve")
    end_annotation = f"Final design:\nObj: {np.around(objvals[-1],1)}\nProds: {x_end[0]}\nStandoffs: {x_end[1]}"
    plt.annotate(end_annotation, (perturbations[-1], objvals[-1]), (0.8*perturbations[-1], 0.8*(np.max(objvals)-np.min(objvals))+np.min(objvals)), arrowprops=dict(facecolor='black', shrink = 0.01, width=0.5))
    best_annotation = f"Best design:\nObj: {np.around(objvals[best_ind],1)}\nProds: {x_best[0]}\nStandoffs: {x_best[1]}"
    mid_ind = int(len(perturbations)/3)
    plt.annotate(best_annotation, (perturbations[best_ind], objvals[best_ind]), (mid_ind, 0.75*(np.max(objvals)-np.min(objvals))+np.min(objvals)), arrowprops=dict(facecolor='black', shrink = 0.01, width=0.5))
    plt.show()







# nsteps = 5
# stepsize_prods = int((nprods_small/2 - nprods_large)/nsteps)
# # stepsize_prods = int((nprods_small - nprods_large)/nsteps)
# stepsize_standoffs = int((nstandoffs_small - nstandoffs_large)/nsteps)

# prods_range = np.arange(nprods_large, int(nprods_small/2), stepsize_prods)
# # standoffs_range = np.arange(nstandoffs_large, nstandoffs_small, stepsize_standoffs)
# standoffs_range = np.arange(2,20,int((20-2)/nsteps))

# # if nprods_small not in prods_range:
# #     prods_range = np.append(prods_range, nprods_small)
    
# # if nstandoffs_small not in standoffs_range:
# #     standoffs_range = np.append(standoffs_range, nstandoffs_small)
# len_prods_range = len(prods_range)
# len_standoffs_range = len(standoffs_range)
    
# pop = []
# vals = []

# for nprods in prods_range:
#     row = []
#     rowvals = []
#     for nstandoffs in standoffs_range:
#         chromosome = constraints.create_chromosome_v3(nprods, nstandoffs, top_constraints, bot_constraints, all_on, on_prob, rod_type)
#         print(f'Chromosome created with {nprods} prods and {nstandoffs} standoffs')
#         row.append(chromosome)
#         rowvals.append((nprods,nstandoffs))
#     pop.append(row)
#     vals.append(rowvals)
    
# popstack = list(itertools.chain.from_iterable(pop))
# valstack = list(itertools.chain.from_iterable(vals))

# prods = []
# valid_circles = []
# tip_radii = []
# drill_radii = []
# top_radii = []
# for i, chromosome in enumerate(popstack):
#     nprods = valstack[i][0]
#     nstandoffs = valstack[i][1]
#     prods.append(constraints.interpret_chromosome_to_prods_v2(chromosome, nprods, nstandoffs))
#     # prods.append(constraints.interpret_chromosome_to_prods(chromosome, nprods))
#     valid_circles.append(constraints.prods_to_valid_circles(prods[i]))
#     tip_radii.append(constraints.prods_to_tip_radii(prods[i]))
#     drill_radii.append(constraints.prods_to_drill_radii(prods[i]))
#     top_radii.append(constraints.prods_to_top_radii(prods[i]))

# print("Running FEA")

# # # Straight calculation version
# # results_mp = [constraints.runFEA_valid_circles_v2(valid_circles[i], tip_radii[i], drill_radii[i], top_radii[i], nprods_top, df_PressureRods, df_Standoffs, root, inputfile, gen, i) for i in range(pop.shape[0])]

# # Multiprocessing version
# ncpus = multiprocessing.cpu_count()
# gen = 0
# pool = multiprocessing.Pool(processes=ncpus)
# arg_tuples = [(valid_circles[i], tip_radii[i], drill_radii[i], top_radii[i], nprods, df_PressureRods, df_Standoffs, root, inputfile, gen, i) for i in range(pop.shape[0])]
# results_mp = pool.starmap(constraints.runFEA_valid_circles_v2, arg_tuples)
# pool.close()
# pool.join()

# # Once all the FEA cases have been accounted for, retrieve results
# results = []
# results_report_all = []
# results_mesh_all = []
# for i in range(pop.shape[0]):
#     # results.append(constraints.read_FEA_results(root, inputfile, gen, i))
#     results_blend, results_report, results_mesh = constraints.read_FEA_results_blend(root, inputfile, gen, i, 10)
#     results.append(results_blend)
#     results_report_all.append(results_report)
#     results_mesh_all.append(results_mesh)

# fitness_values = np.array(results)
# fitness_report = np.array(results_report_all)
# fitness_mesh = np.array(results_mesh_all)

# fitness_report_sum = np.sum(fitness_report, axis=1)
# min_fit_index = np.argmin(fitness_report_sum)

# print('Best combination:')
# print(f'Pressure rods: {valstack[min_fit_index][0]}')
# print(f'Standoffs: {valstack[min_fit_index][1]}')

# end_time = time.time()

# print(f"Total elapsed time:\t{end_time-start_time}")