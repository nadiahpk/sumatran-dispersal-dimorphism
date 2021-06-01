# plot the PIP using interpolation as needed

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


# parameter values
# ---

vmax = 1e-2

# where results will be stored
dir_results = '../../results/stochastic/'

# new grid number
# ngrid = 100
ngrid = 50

#suffixV = ['_1', '_4']
#run_numberV = [1, 3]
suffixV = ['_1', '_4']
run_numberV = [2]

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


        # add to averaging
        # ---

        fitVV_avg = fitVV_avg + fitVV
        cnt += 1


        # plot
        # ---

        # plt.contour(m_resVV, m_mutVV, fitVV, levels=[0], linewidths=[0.5], colors=['black'])
        cs = plt.contour(m_resVV, m_mutVV, fitVV, levels=[0], linewidths=[0.5], colors=['black'])
        csV.append(cs)

plt.close()

for cs in csV:
    ps = cs.collections[0].get_paths()

    # create one long path omitting diagonal elements
    long_p = list()
    for p in ps:

        v = p.vertices
        x = v[:,0]
        y = v[:,1]

        # omit diagonals
        xs = list()
        ys = list()

        for xx, yy in zip(x,y):
            if abs(xx - yy) > xx/1e3:
                long_p.append( (xx,yy) )
                #xs.append(xx)
                #ys.append(yy)

    # sort by y axis
    long_p_s = sorted(long_p, key=lambda v: v[1])

    # plot
    xs, ys = zip(*long_p_s)
    plt.plot(xs, ys, lw=0.5, color='black')


# center points, for use with shading=flat
# ---

del_pwr = (pwrV[1] - pwrV[0])/2
pwrV2 = pwrV - del_pwr
pwrV2 = list(pwrV2) + [pwrV[-1]+del_pwr]
m_resVV2 = [ 10**pwr for pwr in pwrV2 ]
m_mutVV2 = [ 10**pwr for pwr in pwrV2 ]


# color plot of average of stochastic runs
# ---
fitVV_avg = fitVV_avg/cnt

fitsurf = plt.pcolormesh(m_resVV2, m_mutVV2, fitVV_avg, vmin=-vmax, vmax=vmax, cmap='seismic_r', shading='flat')

cbar = plt.colorbar(fitsurf)
cbar.set_label('mutant stochastic invasion fitness')


# put on top the isoclines found using the cyclical disturbance
# ---

plt.yscale('log')
plt.xscale('log')

plt.axis([1e-6, 1, 1e-6, 1])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel(r'resident strategy, $m$')
plt.ylabel(r"mutant strategy, $m'$")

plt.show()

'''
fname = dir_results + 'pip' + suffix + '_' + str(run_number) + type_vmax + '.pdf'
plt.tight_layout()
plt.savefig(fname)
plt.close()
'''



