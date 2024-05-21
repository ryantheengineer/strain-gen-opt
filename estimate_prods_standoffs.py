# -*- coding: utf-8 -*-
"""
Created on Mon May 20 21:34:03 2024

@author: ryanl
"""
import constraints
import numpy as np
import itertools
import multiprocessing
from datetime import datetime
import time

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
nprods_small, nprods_large = constraints.grid_nprods_v2(pBoards_diff_top)
nstandoffs_small, nstandoffs_large = constraints.grid_nprods_v2(pBoards_diff_bot)

nsteps = 5
stepsize_prods = int((nprods_small/2 - nprods_large)/nsteps)
# stepsize_prods = int((nprods_small - nprods_large)/nsteps)
stepsize_standoffs = int((nstandoffs_small - nstandoffs_large)/nsteps)

prods_range = np.arange(nprods_large, int(nprods_small/2), stepsize_prods)
# standoffs_range = np.arange(nstandoffs_large, nstandoffs_small, stepsize_standoffs)
standoffs_range = np.arange(2,20,int((20-2)/nsteps))

# if nprods_small not in prods_range:
#     prods_range = np.append(prods_range, nprods_small)
    
# if nstandoffs_small not in standoffs_range:
#     standoffs_range = np.append(standoffs_range, nstandoffs_small)
len_prods_range = len(prods_range)
len_standoffs_range = len(standoffs_range)
    
pop = []
vals = []

for nprods in prods_range:
    row = []
    rowvals = []
    for nstandoffs in standoffs_range:
        chromosome = constraints.create_chromosome_v3(nprods, nstandoffs, top_constraints, bot_constraints, all_on, on_prob, rod_type)
        print(f'Chromosome created with {nprods} prods and {nstandoffs} standoffs')
        row.append(chromosome)
        rowvals.append((nprods,nstandoffs))
    pop.append(row)
    vals.append(rowvals)
    
popstack = list(itertools.chain.from_iterable(pop))
valstack = list(itertools.chain.from_iterable(vals))

prods = []
valid_circles = []
tip_radii = []
drill_radii = []
top_radii = []
for i, chromosome in enumerate(popstack):
    nprods = valstack[i][0]
    nstandoffs = valstack[i][1]
    prods.append(constraints.interpret_chromosome_to_prods_v2(chromosome, nprods, nstandoffs))
    # prods.append(constraints.interpret_chromosome_to_prods(chromosome, nprods))
    valid_circles.append(constraints.prods_to_valid_circles(prods[i]))
    tip_radii.append(constraints.prods_to_tip_radii(prods[i]))
    drill_radii.append(constraints.prods_to_drill_radii(prods[i]))
    top_radii.append(constraints.prods_to_top_radii(prods[i]))

print("Running FEA")

# # Straight calculation version
# results_mp = [constraints.runFEA_valid_circles_v2(valid_circles[i], tip_radii[i], drill_radii[i], top_radii[i], nprods_top, df_PressureRods, df_Standoffs, root, inputfile, gen, i) for i in range(pop.shape[0])]

# Multiprocessing version
ncpus = multiprocessing.cpu_count()
gen = 0
pool = multiprocessing.Pool(processes=ncpus)
arg_tuples = [(valid_circles[i], tip_radii[i], drill_radii[i], top_radii[i], nprods, df_PressureRods, df_Standoffs, root, inputfile, gen, i) for i in range(pop.shape[0])]
results_mp = pool.starmap(constraints.runFEA_valid_circles_v2, arg_tuples)
pool.close()
pool.join()

# Once all the FEA cases have been accounted for, retrieve results
results = []
results_report_all = []
results_mesh_all = []
for i in range(pop.shape[0]):
    # results.append(constraints.read_FEA_results(root, inputfile, gen, i))
    results_blend, results_report, results_mesh = constraints.read_FEA_results_blend(root, inputfile, gen, i, 10)
    results.append(results_blend)
    results_report_all.append(results_report)
    results_mesh_all.append(results_mesh)

fitness_values = np.array(results)
fitness_report = np.array(results_report_all)
fitness_mesh = np.array(results_mesh_all)

fitness_report_sum = np.sum(fitness_report, axis=1)
min_fit_index = np.argmin(fitness_report_sum)

print('Best combination:')
print(f'Pressure rods: {valstack[min_fit_index][0]}')
print(f'Standoffs: {valstack[min_fit_index][1]}')

end_time = time.time()

print(f"Total elapsed time:\t{end_time-start_time}")