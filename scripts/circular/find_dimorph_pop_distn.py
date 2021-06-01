# for each dimorphic steady state, find the proportion of each morph on the large island
# (large island only bc one morph has zero dispersal)

import pandas as pd
import os
import csv

import sys
sys.path.insert(0,'../../functions') # so I can import the functions
from cycle_funcs import sample_wsV_dimorph

# parameters
# ---

# the low-dispersal morph has zero dispersal
m_res1 = 0 # I know this from the TEP

# which results to use
suffix = '_1'   # which parameter set to use

# which IDs to do
IDV = [0]                                               # default values
IDV += list(range(1,19)) + list(range(95,104))          # vary dispersal cost, c
IDV += [36, 37, 38, 39, 40] + list(range(72, 95))       # vary intrinsic growth rate, r
IDV += [22, 23, 24, 25, 26] + list(range(104, 114))     # vary probability of catastrophe

# where results will be stored
dir_results = '../../results/circular/'

# algorithm parameter values
params = {
        'tol_res': 1e-10,           # the euclidean distance in final resident structure before say it's equilibriated
        'tol_res_dimorph': 1e-6,    # the euclidean distance in final dimorphic resident structure before say it's equilibriated
        'tol_mut': 1e-10,           # the euclidean distance in final mutant structure before say it's equilibriated
        'delta_m': 1e-9,            # step size to estimate gradients; default heuristic sqrt(machine error) * approx value
        }


# get the results
# ---

# open csv file and read
fname = dir_results + 'dimorph_steady_state' + suffix + '.csv'
df = pd.read_csv(fname)

# find all the rows corresponding to the IDV
df2 = df.loc[df['ID'].isin(IDV)]


# if the csv file doesn't exist yet, write it and write headers
# ---

# column headers for the parameter values
par_names = ['suffix', 'layout', 'f', 'c', 'r', 'KL', 'h', 'KS_mean', 'p_cat']

# if the csv file isn't there, write the headers

fname = dir_results + 'dimorph_pop_distn' + suffix + '.csv'
if not os.path.isfile(fname):

    # append column headers for results from this script
    col_names = ['ID'] + par_names + ['m_res2_ss', 'lo_disp_popsize', 'hi_disp_popsize']

    with open(fname, "w", newline="") as ftarget:

        writer = csv.writer(ftarget)
        writer.writerow(col_names)




# find the population size of each morph on the mainland at steady state
# ---

for index, row in df2.iterrows():

    ID = row['ID']
    print('doing ID = ' + str(ID))

    # update params
    for par_name in par_names:
        params[par_name] = row[par_name]

    # run the population at the dimorphic evolutionary steady to population-dynamic steady state
    m_res2 = row['m_res2_ss']
    wsV, nL = sample_wsV_dimorph(m_res1, m_res2, params, return_nT = True)

    # calculate the proportion of zero-dispersal morphs on mainland
    lo_disp_popsize = nL[0][0] # zero-dispersal morphs on mainland
    hi_disp_popsize = nL[1][0] # high-dispersal morphs on mainland

    # write result row to file

    with open(fname, "a", newline="") as ftarget:

        writer = csv.writer(ftarget)
        res = [ID] + [ params[par_name] for par_name in par_names ] + [m_res2, lo_disp_popsize, hi_disp_popsize]
        writer.writerow(res)

