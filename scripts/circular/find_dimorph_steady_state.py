# This should only be applied to parameter values like the default shown in the main text (high f)
# where we know that one isocline is an attractor on zero dispersal for one morph and the other isocline intersects it.

import numpy as np
import os
import pandas as pd
import csv

import sys
sys.path.insert(0,'../../functions') # so I can import the functions
from cycle_funcs import calc_dfitness, sample_wsV_dimorph


# model parameters
# ---

# which parameter set to find the dimorphic evolutionary steady state for
suffix = '_1'

# default
# IDV = [0]

'''
# vary dispersal cost, c
IDV = list(range(1,19)) + list(range(95,104))

# vary intrinsic growth rate, r
IDV = [36, 37, 38, 39, 40] + list(range(72, 95))

# vary probability of catastrophe
IDV = [22, 23, 24, 25, 26] + list(range(104, 114))
'''


# vary dispersal cost, c


# vary probability of catastrophe

# where results will be stored
dir_results = '../../results/circular/'

# maximum derivative (for the resident 2 mutant invasion fitness) that we'll accept as being close enough to 0 
tol_dfit = 1e-7

# default algorithm parameters
params = {
        'tol_res': 1e-10,       # the euclidean distance in final resident structure before say it's equilibriated
        'tol_res_dimorph': 1e-6,# the euclidean distance in final dimorphic resident structure before say it's equilibriated
        'tol_mut': 1e-10,       # the euclidean distance in final mutant structure before say it's equilibriated
        'delta_m': 1e-9,        # step size to estimate gradients; default heuristic sqrt(machine error) * approx value
        }


# write the output file and headers if it doesn't exist already
# ---

# column headers for the parameter values
par_names = ['suffix', 'layout', 'f', 'c', 'r', 'KL', 'h', 'KS_mean', 'p_cat']

fname_out = dir_results + 'dimorph_steady_state' + suffix + '.csv'

if not os.path.isfile(fname_out):

    # append column headers for results from this script
    col_names = ['ID'] + par_names + ['m_res2_ss', 'dfit2']

    with open(fname_out, "w", newline="") as ftarget:

        writer = csv.writer(ftarget)
        writer.writerow(col_names)


# read in the info from the sing_strat_x.csv file 
# ---

fname_in = dir_results + 'sing_strat' + suffix + '.csv'
df = pd.read_csv(fname_in)

for ID in IDV:

    print('---')
    print('finding dimorphic ss for ID = ' + str(ID))

    # repopulate params dictionary
    # ---

    ss_res = df.iloc[ID] # the particular row we want

    for par_name in par_names:

        params[par_name] = ss_res[par_name]


    # find the resident-2 strategy where the mutant invasion fitness gradient is zero
    # ---

    # this is the upper bound on the resident-2 strategy
    # where the mutant invasion fitness gradient is negative
    y_intercept = ss_res['y_intercept']

    # the dimorphic steady state occurs when one of the morphs has zero dispersal
    m_res1 = 0

    # initialise search
    m_res2_lo = y_intercept
    dfit2_lo = -1   # at this point, the resident-2 fitness gradient is negative
    nL = None

    # find where resident 2 invasion fitness gradient goes positive (lower bound)
    # ---

    print('find where gradient goes positive...')

    while dfit2_lo < 0:

        print(m_res2_lo)

        # update
        m_res2_hi = m_res2_lo
        dfit2_hi = dfit2_lo
        nL_hi = nL

        # halve previous bound for the new low estimate
        m_res2_lo = m_res2_hi / 2

        # find the fitness gradient here
        wsV, nL = sample_wsV_dimorph(m_res1, m_res2_lo, params, return_nT=True, nL=nL)
        dfit2_lo = calc_dfitness(m_res2_lo, wsV, params)

    nL_lo = nL

    print('m_res2_lo = ' + str(m_res2_lo))


    # use bisection method to find the root
    # ---

    # using bisection because it allows me to use previously found 
    # the resident-species population structure nL, which speeds things up

    print('find root...')

    dfit2_mid = dfit2_lo
    while abs(dfit2_mid) > tol_dfit:

        m_res2_mid = (m_res2_lo + m_res2_hi) / 2
        nL_mid = [ (nL_lo[0]+nL_hi[0])/2 , (nL_lo[1]+nL_hi[1])/2 ]
        wsV_mid, nL_mid = sample_wsV_dimorph(m_res1, m_res2_mid, params, return_nT=True, nL=nL_mid)
        dfit2_mid = calc_dfitness(m_res2_mid, wsV_mid, params) 

        if np.sign(dfit2_lo) == np.sign(dfit2_mid):
            nL_lo = nL_mid
            wsV_lo = wsV_mid
            dfit2_lo = dfit2_mid
            m_res2_lo = m_res2_mid
        else:
            nL_hi = nL_mid
            wsV_hi = wsV_mid
            dfit2_hi = dfit2_mid
            m_res2_hi = m_res2_mid

        print(m_res2_mid)

    print('m_res2_mid = ' + str(m_res2_mid))


    # write this result to the csv
    # ---

    with open(fname_out, "a", newline="") as ftarget:

        writer = csv.writer(ftarget)

        # the current ID
        res = [ID]

        # get the parameter values
        res += [ params[par_name] for par_name in par_names ]

        # append results from this loop step
        res += [ m_res2_mid, dfit2_mid ]

        # write the result row
        writer.writerow(res)

