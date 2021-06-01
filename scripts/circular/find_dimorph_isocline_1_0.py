# get the data needed to plot the TEP with a quiver plot on top of it

import numpy as np
#import pickle
import os
import pandas as pd
import csv
from scipy.interpolate import interp1d

import sys
sys.path.insert(0,'../../functions') # so I can import the functions
from cycle_funcs import calc_fitness, calc_dfitness, sample_wsV_dimorph

# =====================================================================================================================

# model parameters
# ---

# which parameter set to find the isoclines for (comment out ones I've done)
suffix = '_1'

'''
# default parameter values
ID = 0
ngrid = 31 # number of grid points in log space along the resident1 axis
'''
# high dispersal cost parameter values
ID = 17
ngrid = 31 # number of grid points in log space along the resident1 axis

# where results will be stored
dir_results = '../../results/circular/'

# default algorithm parameters
params = {
        'tol_res': 1e-10,       # the euclidean distance in final resident structure before say it's equilibriated
        'tol_res_dimorph': 1e-6,# the euclidean distance in final dimorphic resident structure before say it's equilibriated
        'tol_mut': 1e-10,       # the euclidean distance in final mutant structure before say it's equilibriated
        'delta_m': 1e-9,        # step size to estimate gradients; default heuristic sqrt(machine error) * approx value
        }


# read in the info from the sing_strat_x.csv file and repopulate params dictionary
# ---

par_names = ['suffix', 'layout', 'f', 'c', 'r', 'KL', 'h', 'KS_mean', 'p_cat'] # names of missing parameter values

fname = dir_results + 'sing_strat' + suffix + '.csv'
df = pd.read_csv(fname)
ss_res = df.iloc[ID] # the particular row we want

for par_name in par_names:

    params[par_name] = ss_res[par_name]


# read in the isoclines to define the TEP region
# ---

fname = dir_results + 'isocline' + suffix + '_' + str(ID) + '.csv'
df = pd.read_csv(fname)
df = df.sort_values(by=['m_mut'])

# grab the isoclines
m_iso1v = df['m_res'].values
m_iso2v = df['m_mut'].values

# flip all points below the diagonal so that we define the TEP region
region = [ (m_mut, m_res) if m_res > m_mut else (m_res, m_mut) for m_res, m_mut in zip(m_iso1v, m_iso2v) ]

# make sure the first entry is 0,0
if region[0] != (0,0):
    region += [(0,0)] + region

# make sure the last entry is 0,y_intercept
if region[-1][1] == 1:
    region[-1] = (0, ss_res['y_intercept'])
else:
    region += [(0, ss_res['y_intercept'])]


# split the region boundary into an upper and lower bound on resident 2
# ---

# find the extremal point of the region in the resident 1 dimension, the bulging out to the right of the PIP graph
m_iso1V, m_iso2V = zip(*region)
extremal_res1 = max(m_iso1V)
extremal_idx = m_iso1V.index( extremal_res1 )
extremal_res2 = m_iso2V[extremal_idx]

# split the region boundary into two lines, one for a lower bound on res2, and one for an upper bound
line_lo = [ (m_res1, m_res2) for m_res1, m_res2 in region if m_res2 <= extremal_res2 ]
line_hi = [ (m_res1, m_res2) for m_res1, m_res2 in region if m_res2 >= extremal_res2 ]

# create functions that will return the lower and upper bound on res2 for a given res1
m_res1V, m_res2V = zip(*line_lo)
f_lo = interp1d(m_res1V, m_res2V)
m_res1V, m_res2V = zip(*line_hi)
f_hi = interp1d(m_res1V, m_res2V)


# create a grid along the resident 1 dimension, and find the isocline at each point along that grid
# ---

# do it in log space
pwrV = np.linspace(-6, np.log10(extremal_res1), ngrid)[:-1]
m_res1V = [ 10**pwr for pwr in pwrV ]

# if the csv file doesn't exist yet, create it, and include the resident 1 = 0 point in the grid
# ---

fname = dir_results + 'dimorph_isocline' + suffix + '_' + str(ID) + '.csv'
if not os.path.isfile(fname):

    # add the zero point to our search
    m_res1V = [0] + m_res1V 

    # write the column headers
    with open(fname, "w", newline="") as ftarget:
        writer = csv.writer(ftarget)
        writer.writerow( ['m_res1', 'm_res2', 'dfit2'] )



# for each resident 1 strategy, find where the resident 2 mutant invasion fitness gradient equals 0
# (I know that this is an attractor)
# ---

for m_res1 in m_res1V:


    m_res2_hi_bnd = f_hi([m_res1])[0]
    m_res2_lo_bnd = f_lo([m_res1])[0]

    # find where resident 2 invasion fitness gradient goes positive
    # ---

    # initialise
    m_res2_lo = m_res2_hi_bnd
    dfit2_lo = -1   # at this point, the resident-2 fitness gradient is negative
    nL = None

    print('find where gradient goes positive')

    while dfit2_lo < 0:

        print(m_res2_lo)

        # update
        m_res2_hi = m_res2_lo
        dfit2_hi = dfit2_lo
        nL_hi = nL

        # halve distance to lower bound for the new low estimate
        m_res2_lo = (m_res2_lo_bnd + m_res2_hi) / 2

        # find the fitness gradient here
        wsV, nL = sample_wsV_dimorph(m_res1, m_res2_lo, params, return_nT=True, nL=nL)
        dfit2_lo = calc_dfitness(m_res2_lo, wsV, params)

    nL_lo = nL

    print('m_res2_lo = ' + str(m_res2_lo))


    # use bisection method to find the root
    # ---

    print('find root')

    tol_dfit = 1e-7 # maximum derivative that we'll accept as being close enough to 0 
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


    # write this isocline point to the csv
    # ---

    with open(fname, "a", newline="") as ftarget:

        writer = csv.writer(ftarget)
        writer.writerow([m_res1, m_res2_mid, dfit2_mid])

