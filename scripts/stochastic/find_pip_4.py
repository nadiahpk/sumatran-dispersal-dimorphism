# create the data needed to plot the PIP for the stochastic invasion fitness of the mutant
# NOTE same as find_pip_1 but I've strengthened the requirement for mutant equilibrium
# 'tol_mut': 1e-6 in params below

import numpy as np
import os
import pickle
import csv

import sys
sys.path.insert(0,'../../functions') # so I can import the functions

from funcs import generate_catastrophes, get_M
from funcs import sample_wsV2 as sample_wsV_stoch
from funcs import calc_fitness2 as calc_fitness

from cycle_funcs import equilib_mutant
from cycle_funcs import sample_wsV as sample_wsV_cycle


# default model parameter values
# ---

'''
run_number = 3  # an ID number so we can do multiple runs of the same parameter set
tlimit = 10000  # length of stochastic sequence used to estimate stochastic invasion fitness
ngrid = 50      # number of grid points to plot on a log scale

run_number = 2  # an ID number so we can do multiple runs of the same parameter set
tlimit = 100000  # length of stochastic sequence used to estimate stochastic invasion fitness
ngrid = 50      # number of grid points to plot on a log scale

run_number = 1  # an ID number so we can do multiple runs of the same parameter set
tlimit = 10000  # length of stochastic sequence used to estimate stochastic invasion fitness
ngrid = 50      # number of grid points to plot on a log scale

run_number = 4  # an ID number so we can do multiple runs of the same parameter set
tlimit = 100000  # length of stochastic sequence used to estimate stochastic invasion fitness
ngrid = 50      # number of grid points to plot on a log scale

run_number = 5  # an ID number so we can do multiple runs of the same parameter set
tlimit = 100000  # length of stochastic sequence used to estimate stochastic invasion fitness
ngrid = 50      # number of grid points to plot on a log scale
'''

run_number = 6  # an ID number so we can do multiple runs of the same parameter set
tlimit = 100000  # length of stochastic sequence used to estimate stochastic invasion fitness
ngrid = 50      # number of grid points to plot on a log scale

pwr_min = -6  # plot will range from [10**pwr_min, 1]
tburnin = 10000

# where results will be stored
dir_results = '../../results/stochastic/'

# default parameter values
params = {
        # file suffix for keeping results together
        'suffix': '_4',
        # algorithm parameters
        'tol_res': 1e-10,       # the euclidean distance in final resident structure before say it's equilibriated
        'tol_mut': 1e-9,       # the euclidean distance in final mutant structure before say it's equilibriated
        'delta_m': 1e-9,        # step size to estimate gradients; default heuristic sqrt(machine error) * approx value
        # model parameters
        'layout': 'circular',   # islands arranged in a circle around mainland
        'f': 0.2 ,              # the max proportion of the mainland that disperses (chosen here so net flow mainland -> small islands)
        'c': 0.01,              # dispersal cost, i.e., probability to die during dispersal
        'r': 0.1,               # recovery rate of populations
        'KL': 500,              # carrying capacity of single large island
        'h': 5,                 # how many small islands
        'KS_mean': 30,          # carrying capacity on small islands
        'p_cat': 0.0005,        # the probability of a catastrophe on any small island (1 every 2000 years)
        }


# get a random catastrophe sequence
# ---

# we might have already generated one, or need to create one
suffix = params['suffix']
fname = dir_results + 'catastrophes' + suffix + '_' + str(run_number) + '.pkl'

if not os.path.isfile(fname):

    # generate a new one
    catastropheD = generate_catastrophes(params, tlimit, fname)

else:

    # read in old one
    f = open(fname, 'rb')
    ss = pickle.load( f )           # explanatory string
    catastropheD = pickle.load( f ) # dictionary of catastrophes { time: (islands having catast)}
    f.close()


# if the output csv file doesn't exist yet write it
# ---

# column headers
col_names = ['m_res', 'm_mut', 'fitness'] # resident strategy, mutant strategy, mutant invasion fitness

# if the csv file isn't there, write the headers

suffix = params['suffix']
fname = dir_results + 'pip' + suffix + '_' + str(run_number) + '.csv'
if not os.path.isfile(fname):

    with open(fname, "w", newline="") as ftarget:

        writer = csv.writer(ftarget)
        writer.writerow(col_names)


# calculate the stochastic invasion fitness of the mutant on a grid
# ---

pwrV = np.linspace(pwr_min, 0, ngrid)
m_resV = [ 10**pwr for pwr in pwrV ]
m_mutV = [ 10**pwr for pwr in pwrV ]

# fitV = list()
for idx, m_res in enumerate(m_resV):

    print(' -- resident index ' + str(idx) )
    print(m_res)

    # sample the catastrophes' effects on resident species and growth rates
    # ---

    # get a sample of the cyclic to use for burnins
    wsV_cycle, nT0 = sample_wsV_cycle(m_res, params, return_nT=True)

    # create the sample of the stochastic to use for calculations
    # - initialised with the steady state from the cyclic
    wsV_stoch = sample_wsV_stoch( m_res, params, catastropheD, tlimit, initial_nT = nT0, tburnin=tburnin )


    for m_mut in m_mutV:

        print(m_mut)
        if m_res == m_mut:

            fit = 0 # on diagonal

        else:

            # calculate the invasion fitness of the mutant given the catastrophes and growth rates
            # ---

            # first, equilibriate the population structure of the mutant using the cyclic result
            # - doing this because it's a faster way to obtain particularly when m_mut close to m_res
            initial_nT_mut = equilib_mutant(m_mut, wsV_cycle, params)

            # then use that equilibrium as the initial structure to calculate the fitness
            fit = calc_fitness(m_mut, wsV_stoch, params, initial_nT = initial_nT_mut, tburnin=tburnin)


        # write line to csv
        # ---

        with open(fname, "a", newline="") as ftarget:

            # write the result row
            res = [m_res, m_mut, fit]
            writer = csv.writer(ftarget)
            writer.writerow(res)



