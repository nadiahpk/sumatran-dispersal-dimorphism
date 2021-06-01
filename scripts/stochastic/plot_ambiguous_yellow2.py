# the point of this is to illustrate that the stochastic simulation is a bit ambiguous about exactly where the isocline is
# I'm going to mark the ambiguity by finding sign changes and put them in yellow
# differs from plot_ambiguous_yellow.py in that I don't plot the periodic-model isoclines on top

# the point of this is to illustrate that the stochastic simulation is a bit ambiguous about exactly where the isocline is
# I'm going to mark the ambiguity by finding sign changes and put them in yellow

import copy
import numpy as np
import pandas as pd
import matplotlib.colors as colors
import matplotlib 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pandas as pd

# from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# parameter values
# ---

vmax = 1e-1

# where results will be stored
dir_results = '../../results/stochastic/'

# new grid number
# ngrid = 100
ngrid = 50

'''
suffixV = ['_4'] # _4 has the more stringent requirement for mutant population structure equilibrium
run_numberV = [2, 4] # run 2 and runs >= 4 have the more stringent tlimit = 100,000
'''
suffixV = ['_1', '_4'] # _4 has the more stringent requirement for mutant population structure equilibrium
run_numberV = [1, 2, 3] # run 2 and runs >= 4 have the more stringent tlimit = 100,000

fit_store_sign = None
fitVV_avg = np.zeros( (ngrid, ngrid) )
cnt = 0
csV = list()
for suffix in suffixV:

    for run_number in run_numberV:

        # read data for stochastic run
        # ---

        fname = dir_results + 'pip' + suffix + '_' + str(run_number) + '.csv'
        df = pd.read_csv(fname)

        m_resV = df['m_res'].values
        m_mutV = df['m_mut'].values
        fitV = df['fitness'].values


        # create new grid
        # ---

        pwrV = np.linspace(-6, 0, ngrid)
        m_resVV = [ 10**pwr for pwr in pwrV ]
        m_mutVV = [ 10**pwr for pwr in pwrV ]

        XX, YY = np.meshgrid(m_resVV, m_mutVV)
        points = (m_resV, m_mutV)
        fitVV = griddata(points, fitV, (XX, YY), method='linear')


        # make a note of which locations have a sign that contradicts previous signs
        if fit_store_sign is None:

            # if this is the first matrix, make this the fitness sign matrix
            fit_store_sign = np.sign(fitVV)

        else:

            # compare the signs we have now, and mark contradictions with a nan
            fit_store_sign[ fit_store_sign != np.sign(fitVV) ] = np.nan


        # add to averaging
        # ---

        fitVV_avg = fitVV_avg + fitVV
        cnt += 1


# center points, for use with shading=flat
# ---

del_pwr = (pwrV[1] - pwrV[0])/2
pwrV2 = pwrV - del_pwr
pwrV2 = list(pwrV2) + [pwrV[-1]+del_pwr]
m_resVV2 = [ 10**pwr for pwr in pwrV2 ]
m_mutVV2 = [ 10**pwr for pwr in pwrV2 ]


# color plot of average of stochastic runs
# ---

# average the fitnesses
fitVV_avg = fitVV_avg/cnt

# modify it by the nans
fitVV_avg[ np.isnan(fit_store_sign) ] = np.nan

# mark every diagonal element as ambiguous because I know they are
np.fill_diagonal(fitVV_avg, np.nan)

current_cmap = copy.copy(matplotlib.cm.get_cmap("seismic_r"))
current_cmap.set_bad(color='yellow')

# fitsurf = plt.pcolormesh(m_resVV2, m_mutVV2, fitVV_avg, vmin=-vmax, vmax=vmax, cmap='seismic_r', shading='flat')
fitsurf = plt.pcolormesh(m_resVV2, m_mutVV2, fitVV_avg, norm=colors.SymLogNorm(linthresh=1e-5, linscale=1, base=10, vmin=-vmax, vmax=vmax), cmap=current_cmap, shading='flat')

cbar = plt.colorbar(fitsurf)
cbar.set_label('mutant stochastic invasion fitness', fontsize='x-large')

# plot decorations
# ---

plt.yscale('log')
plt.xscale('log')

plt.axis([1e-6, 1, 1e-6, 1])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel(r'resident dispersal trait, $m$', fontsize='x-large')
plt.ylabel(r"mutant dispersal trait, $m'$", fontsize='x-large')

legend_elements = [
        Line2D([0], [0], marker='s', color='w', label='mutant can invade', markerfacecolor='blue', markersize=10),
        Line2D([0], [0], marker='s', color='w', label='mutant cannot invade', markerfacecolor='red', markersize=10),
        Line2D([0], [0], marker='s', color='w', label='invasibility varies between\nstochastic estimates', markerfacecolor='yellow', markersize=10)]

# Create the figure
plt.legend(handles=legend_elements, loc='lower right', framealpha=1, fontsize='medium')

# plt.show()

fname = dir_results + 'ambiguous_yellow2.pdf'
plt.tight_layout()
plt.savefig(fname)
plt.close()
