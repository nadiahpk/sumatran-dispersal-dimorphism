# get the data needed to plot the TEP with a quiver plot on top of it

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

import sys
sys.path.insert(0,'../../functions') # so I can import the functions
from cycle_funcs import calc_fitness, calc_dfitness, sample_wsV_dimorph

# =====================================================================================================================

# model parameters
# ---

# which parameter set to use
suffix = '_1'
ID = 0
#ngrid = 31 # number of grid points in log space along the resident1 axis
ngrid = 20 # number of grid points in log space along the resident1 axis

# where results are kept
dir_results = '../../results/circular/'

# algorithm parameter values
params = {
        'tol_res': 1e-10,           # the euclidean distance in final resident structure before say it's equilibriated
        'tol_res_dimorph': 1e-6,    # the euclidean distance in final dimorphic resident structure before say it's equilibriated
        'tol_mut': 1e-10,           # the euclidean distance in final mutant structure before say it's equilibriated
        'delta_m': 1e-9,            # step size to estimate gradients; default heuristic sqrt(machine error) * approx value
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

tep_polygon = Polygon(region)


# define the space to grid over
# ---

m_iso1V, m_iso2V = zip(*region)

# find the extents of the region
min_m_iso1 = min(m_iso1V); max_m_iso1 = max(m_iso1V)
min_m_iso2 = min(m_iso2V); max_m_iso2 = max(m_iso2V)

# put the grid on a log scale, but:
#   start it at 1e-6 instead of 0
#   skip the maximum in the resident 1 direction because that is the extremal point on the curve
#pwr_res1V = np.linspace( -6, np.log10(max_m_iso1), ngrid+1 )[:-1]
#pwr_res2V = np.linspace( -6, np.log10(max_m_iso2), ngrid )
pwr_res1V = np.linspace( -6, 0, ngrid ) # make evenly spaced instead
pwr_res2V = np.linspace( -6, 0, ngrid )
res1_gridV = [ 10**pwr for pwr in pwr_res1V ]
res2_gridV = [ 10**pwr for pwr in pwr_res2V ]

lenlen = ngrid*ngrid

# calculate the derivatives for each grid point in the dimorphism region
# ---

m_res1V = list(); m_res2V = list()  # store which res1, res2 pairs were within the region
dfit1V = list(); dfit2V = list()    # store derivates found at each res1, res2 pair
nL_prev = None
cnt = 0
for m_res2 in reversed(res2_gridV):

    # going in reverse direction because res1 far from res2 is easiest to equilibriate,
    # and we can use the nL found at the previous point as the first estimate for the next point
    nL = nL_prev
    store_nL_flag = True

    # nest resident 1 so we can stop loop and save when we get too close to the bottom left corner
    # NOTE if you reverse this one too, might be easier
    for m_res1 in res1_gridV:

        point = Point(m_res1, m_res2)

        if tep_polygon.contains(point): # if this point is within the region, find the derivatives

            cnt += 1
            print(str(cnt) + ' of ' + str(lenlen))

            # append point to our list
            m_res1V.append(m_res1)
            m_res2V.append(m_res2)

            wsV, nL = sample_wsV_dimorph(m_res1, m_res2, params, return_nT = True, nL = nL) # calculate derivatives

            if store_nL_flag:
                # if this is the first point in m_res2, store nL for next time
                nL_prev = nL
                store_nL_flag = False

            # calculate derivatives and append
            dfit1 = calc_dfitness(m_res1, wsV, params) 
            dfit2 = calc_dfitness(m_res2, wsV, params) 
            # fake data for checking
            #dfit1 = -1 + 2*np.random.rand()
            #dfit2 = -1 + 2*np.random.rand()
            dfit1V.append( dfit1 )
            dfit2V.append( dfit2 )


# write output to csv file
# ---

df = pd.DataFrame(list(zip(m_res1V, m_res2V, dfit1V, dfit2V)), columns =['m_res1', 'm_res2', 'dfit1', 'dfit2'])
fname = dir_results + 'dimorph_derivs_ngrid' + str(ngrid) + suffix + '_' + str(ID) + '.csv'
df.to_csv(fname, index=False)

'''
In [265]: m_res1
Out[265]: 1e-06

In [266]: m_res2
Out[266]: 2.06913808111479e-06
'''


