# find for various parameter values: where the mutant-fitness isocline intersects the mutant-strategy axis, the singular strategy, 
# and the derivative of the fitness gradient at the singular strategy

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import os
import csv

import sys
sys.path.insert(0,'../../functions') # so I can import the functions
from funcs import get_M
from cycle_funcs import *


# which parameter we vary and its values
# ---

# comment out those I've done

vary_par_name = 'c'
vary_par_vals = [0.01] #, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5]


# model parameters
# ---

# where results will be stored
dir_results = '../../results/circular/'

# default parameter values
params = {
        # file suffix for keeping results together
        'suffix': '_2',
        # algorithm parameters
        'tol_res': 1e-10,       # the euclidean distance in final resident structure before say it's equilibriated
        'tol_mut': 1e-10,       # the euclidean distance in final mutant structure before say it's equilibriated
        'delta_m': 1e-9,        # step size to estimate gradients; default heuristic sqrt(machine error) * approx value
        # model parameters
        'layout': 'circular',   # islands arranged in a circle around mainland
        'f': 0.2 ,              # the max proportion of the mainland that disperses (chosen here so net flow mainland -> small islands)
        'c': 0.01,              # dispersal cost, i.e., probability to die during dispersal
        'r': 0.1,               # recovery rate of populations
        'KL': 700,              # carrying capacity of single large island
        'h': 7,                 # how many small islands
        'KS_mean': 30,          # carrying capacity on small islands
        'p_cat': 0.0005,        # the probability of a catastrophe on any small island (1 every 2000 years)
        }


# if the csv file doesn't exist yet, write it and write headers
# ---

# column headers for the parameter values
par_names = ['suffix', 'layout', 'f', 'c', 'r', 'KL', 'h', 'KS_mean', 'p_cat']

# if the csv file isn't there, write the headers

suffix = params['suffix']
fname = dir_results + 'sing_strat' + suffix + '.csv'
if not os.path.isfile(fname):

    # append column headers for results from this script
    col_names = ['ID'] + par_names + ['y_intercept', 'm_ss', 'ddfit_ss']

    with open(fname, "w", newline="") as ftarget:

        writer = csv.writer(ftarget)
        writer.writerow(col_names)

    # set ID to initial value
    ID = -1

else:

    # get the ID from the last row

    with open(fname, "r", errors="ignore") as ftarget:

        reader = csv.reader(ftarget)

        for row in reader:
            pass
    
        ID = int(row[0])

# loop through parameter values, calculating y_intercept, sing strat, and ddfit for each
# ---

for vary_par_val in vary_par_vals:


    # update the parameter value
    # ---

    params[vary_par_name] = vary_par_val
    print( vary_par_name + ' = ' + str(vary_par_val) )


    # find intercept of isocline with mutant-strategy axis (y_intercept)
    # ---

    print('find where mutant fitness isocline intersects mutant strategy axis')

    m_res = 0
    wsV = sample_wsV(m_res, params)

    # start by checking if the mutant with strategy m = 1 can invade
    m_mut = 1
    fit_fnc = lambda m: calc_fitness(m, wsV, params)

    fit = fit_fnc(m_mut)

    if fit > 0:

        # can invade, then the isocline doesn't intersect the axis
        y_intercept = 1

    else:

        # find where the isocline intersects

        # find where fitness becomes positive
        m_hi = 1   # fit < 0
        m_lo = 0.5 # first guess for fit > 0

        while fit_fnc(m_lo) < 0:

            m_hi = m_lo
            m_lo = m_lo / 2

        # use brentq to find good estimate
        y_intercept = brentq(fit_fnc, m_lo, m_hi)


    # find the singular strategy (m_ss)
    # ---

    print('find singular strategy')

    dfit_fnc = lambda m: calc_dfitness(m, sample_wsV(m, params), params)

    # using brentq, so need an upper bound (where dfit < 0) and lower bound (where dfit > 0)
    m_hi = 1
    #m_lo = 1e-6 # first guess
    m_lo = 0.01 # debugging

    while dfit_fnc(m_lo) < 0:

        m_hi = m_lo
        m_lo = m_lo / 2

    # now find the singular strategy, where dfit = 0
    m_ss = brentq(dfit_fnc, m_lo, m_hi)


    # find the second derivative of the fitness gradient at the singular strategy (dfitgrad_ss)
    # ---

    wsV_ss = sample_wsV(m_ss, params)
    ddfit_ss = calc_ddfitness(m_ss, wsV_ss, params)


    # write result to file
    # ---

    with open(fname, "a", newline="") as ftarget:

        writer = csv.writer(ftarget)

        # the current ID
        ID += 1
        res = [ID]

        # get the parameter values
        res += [ params[par_name] for par_name in par_names ]

        # append results from this loop step
        res += [ y_intercept, m_ss, ddfit_ss ]

        # write the result row
        writer.writerow(res)


