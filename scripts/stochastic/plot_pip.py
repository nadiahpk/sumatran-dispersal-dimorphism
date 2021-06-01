# plot the PIP using interpolation as needed

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


# parameter values
# ---

type_vmax = '_lo_v'
vmax = 1e-7

'''
suffix = '_1' 
run_number = 1
run_number = 2
run_number = 3
'''

'''
suffix = '_4'
run_number = 1
run_number = 2
run_number = 3
'''

# where results will be stored
dir_results = '../../results/stochastic/'

# new grid number
# ngrid = 100
ngrid = 50


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


# center points, for use with shading=flat
# ---

del_pwr = (pwrV[1] - pwrV[0])/2
pwrV2 = pwrV - del_pwr
pwrV2 = list(pwrV2) + [pwrV[-1]+del_pwr]
m_resVV2 = [ 10**pwr for pwr in pwrV2 ]
m_mutVV2 = [ 10**pwr for pwr in pwrV2 ]


# plot
# ---

fname = dir_results + 'pip' + suffix + '_' + str(run_number) + type_vmax + '.pdf'

plt.pcolormesh(m_mutVV2, m_resVV2, fitVV, vmin=-vmax, vmax=vmax, cmap='seismic_r', shading='flat')

plt.yscale('log')
plt.xscale('log')

plt.axis([1e-6, 1, 1e-6, 1])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel(r'resident strategy, $m$')
plt.ylabel(r"mutant strategy, $m'$")

plt.colorbar()
plt.tight_layout()
plt.savefig(fname)
plt.close()



