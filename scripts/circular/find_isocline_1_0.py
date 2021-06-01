
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import pandas as pd
import csv
import os

import sys
sys.path.insert(0,'../../functions') # so I can import the functions
from funcs import get_M
from cycle_funcs import *


# model parameters
# ---

# choose grid

eps = 1e-8 # a small epsilon to keep the search off the diagonal
'''

# general grid for linear scale
m_mutV = np.linspace(0, 1, 50)

# get the curve
m_mutV = np.linspace(0.02, 0.06, 30)
m_mutV = np.linspace(0.0009, 0.0198, 10)
'''
m_mutV = np.linspace(0.001, 0.0028, 10)

'''
# get the area below the singular strategy

eps = 1e-9 # a small epsilon to keep the search off the diagonal

pwrV = np.linspace(-6, -4, 20)
m_mutV = [ 10**pwr for pwr in pwrV ]
pwrV = np.linspace(-4, -3.067, 10)[1:-1] # up to near where the sing strat is
m_mutV = [ 10**pwr for pwr in pwrV ]
'''

# which parameter set to do the pip for
suffix = '_1'
ID = 0


# where results will be stored
dir_results = '../../results/circular/'

# default algorithm parameters
params = {
        'tol_res': 1e-10,       # the euclidean distance in final resident structure before say it's equilibriated
        'tol_mut': 1e-10,       # the euclidean distance in final mutant structure before say it's equilibriated
        'delta_m': 1e-9,        # step size to estimate gradients; default heuristic sqrt(machine error) * approx value
        }

# names of missing parameter values
par_names = ['suffix', 'layout', 'f', 'c', 'r', 'KL', 'h', 'KS_mean', 'p_cat']


# construct file names for getting and writing results
# ---

fname_in = dir_results + 'sing_strat' + suffix + '.csv'
fname_out = dir_results + 'isocline' + suffix + '_' + str(ID) + '.csv'


# read in the information from the sing_strat_x.csv file
# ---

df = pd.read_csv(fname_in)

# the particular row we want
ss_res = df.iloc[ID]

# useful information about this parameter-value set
y_intercept = ss_res['y_intercept'] # where the isocline intersects the mutant-strategy axis (1 if it doesn't)
m_ss = ss_res['m_ss']               # singular strategy
ddfit_ss = ss_res['ddfit_ss']       # derivative of the fitness gradient (how strong divergent selection is)


# put the parameter values back in the dictionary
# ---

for par_name in par_names:

    params[par_name] = ss_res[par_name]

'''
# double-check that I've got the parameter values right and still finding the same m_ss
# ---

wsV_ss = sample_wsV(m_ss, params)
dfit = calc_dfitness(m_ss, wsV_ss, params)      # -- check this is close to zero
ddfit = calc_ddfitness(m_ss, wsV_ss, params)    # -- check this is close to ss_res['ddfit_ss']
'''

# find isocline's mutant-strategy values for grid of resident-strategy values
# ---



# initalise results matrix

resM = list() # [m_res, m_mut]

# append values that are already known if needed
if not os.path.isfile(fname_out):
    resM.append( [m_ss, m_ss] )         # isocline passes through the singular strategy
    if y_intercept < 1:
        resM.append( [y_intercept, 0] ) # where isocline intersects mutant-strategy axis


# loop through grid

for m_mut in m_mutV:

    print(m_mut)

    # fitness of the invading mutant as a function of the resident strategy
    fitF = lambda m_res: calc_fitness(m_mut, sample_wsV(m_res, params), params)

    # the range we search depends on whether we're above or below the m_ss
    if m_mut < m_ss:
        m_mut_lo = m_mut + eps # little eps is needed because the fitness at the diagonal is approx 0
        m_mut_hi = 1
    else:
        m_mut_lo = 0
        m_mut_hi = m_mut - eps

    # check the range we've picked has opposite-signed mutant invasion fitnesses
    fit_lo = fitF(m_mut_lo)
    fit_hi = fitF(m_mut_hi)

    if np.sign(fit_lo) == np.sign(fit_hi):

        m_res_iso = np.nan # skip this one

    else:

        # find the resident strategy at which the invading mutant fitness is zero
        m_res_iso = brentq(fitF, m_mut_lo, m_mut_hi)

    resM.append( [m_res_iso, m_mut] ) # store new value


# write a new isocline file, or append if it already exists
# ---

df = pd.DataFrame(resM, columns = ['m_res', 'm_mut'])

if not os.path.isfile(fname_out):

    df.to_csv(fname_out, header=True, index=False)

else:

    df.to_csv(fname_out, mode='a', header=False, index=False)

