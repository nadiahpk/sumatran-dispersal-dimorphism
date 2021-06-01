# Redo a PIP the default values along with some annotations to explain how to interpret the thing
# I'm going to mark the ambiguity by finding sign changes and put them in yellow

import numpy as np
import pandas as pd
#import matplotlib 
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D


# parameter values
# ---

# where results will be stored
dir_results = '../../results/circular/'


suffix_iso = '_2'
ID_iso = 0 # default parameter values I am also using here for the stochastic run

# plot isocline 
# ---

fname = dir_results + 'isocline' + suffix_iso + '_' + str(ID_iso) + '.csv'
df = pd.read_csv(fname)

df = df.sort_values(by=['m_mut'])   # sort so we can plot the line nicely

# grab the isoclines
m_iso1v = df['m_res'].values
m_iso2v = df['m_mut'].values

# filter out the nans
m_iso1_2v = list(zip(m_iso1v, m_iso2v))
m_iso1_2v_filt = [ v for v in m_iso1_2v if not np.isnan(v[0]) ]
m_iso1v, m_iso2v = zip(*m_iso1_2v_filt)

# plot log scale
diagx = [0, 1]; diagy = [0, 1]
plt.plot(diagx, diagy, color='black')
plt.plot(m_iso1v, m_iso2v, color='black')


# plot a point at the singular strategy
# ---

plt.yscale('log')
plt.xscale('log')

# read in the information from the sing_strat_x.csv file
fname_in = dir_results + 'sing_strat' + suffix_iso + '.csv'
df = pd.read_csv(fname_in)
ss_res = df.iloc[ID_iso]                # the particular row we want
m_ss = ss_res['m_ss']               # singular strategy

plt.scatter([m_ss], [m_ss], color='black', zorder=9999)

plt.text(5e-4, 8e-4, r'singular strategy, $m^*$', fontsize='medium', ha='right')

# fill colour regions
# ---

# split m_iso2v up
m_iso2v_lo, m_iso1v_lo = zip(*[ (y, x) for y, x in zip(m_iso2v, m_iso1v) if y <= m_ss])
m_iso2v_hi, m_iso1v_hi = zip(*[ (y, x) for y, x in zip(m_iso2v, m_iso1v) if y >= m_ss])

plt.fill_between([0,1], [0,0], 1, facecolor='#ffcdcd', interpolate=True)

plt.fill_between(m_iso1v_hi, m_iso2v_hi, m_iso1v_hi, facecolor='#d9d9ff', interpolate=True)
plt.fill_between(list(m_iso1v_lo) + [1], list(m_iso2v_lo) + [1], 0, facecolor='#d9d9ff', interpolate=True)


'''

# draw mutations
# ---

# mutations below m_ss

pwrV = np.linspace(np.log10(2e-6), np.log10(2e-4), 12)
vV = [ 10**pwr for pwr in pwrV ]

xV = list()
for v in vV[:-1]:
    xV.append(v)
    xV.append(v)
xV.append(vV[-1])

yV = [vV[0]]
for v in vV[1:]:
    yV.append(v)
    yV.append(v)

plt.plot(xV, yV, color='black', lw=0.5)

plt.plot(xV[-1]-xV[-1]/10, yV[-1], marker=5, color='black')

plt.text(1.2e-2, 1.5e-3, 'mutations drive\ntrait towards\nsingular strategy', fontsize='medium', ha='left')

# mutations above m_ss

pwrV = np.linspace(np.log10(7e-3), np.log10(7e-1), 12)
vV = [ 10**pwr for pwr in reversed(pwrV) ]

xV = list()
for v in vV[:-1]:
    xV.append(v)
    xV.append(v)
xV.append(vV[-1])

yV = [vV[0]]
for v in vV[1:]:
    yV.append(v)
    yV.append(v)

plt.plot(xV, yV, color='black', lw=0.5)

plt.plot(xV[-1]+xV[-1]/10, yV[-1], marker=4, color='black')


# draw the divergent selection
# ----

pwr_mss = np.log10(m_ss)

# up arrow
yV = [ 10**(pwr_mss+0.2), 10**(pwr_mss+1) ]
plt.plot([m_ss, m_ss], yV, color='black')
plt.scatter([m_ss], [yV[-1]], marker=6, color='black')

# down arrow
yV = [ 10**(pwr_mss-0.2), 10**(pwr_mss-1) ]
plt.plot([m_ss, m_ss], yV, color='black')
plt.scatter([m_ss], [yV[-1]], marker=7, color='black')

plt.text(3e-6, 1e-2, 'both mutants with higher\nand lower dispersal can\ninvade singular strategy\ni.e., divergent selection', fontsize='medium', ha='left')

'''

# plot decorations
# ---

# legend
legend_elements = [
        Line2D([0], [0], marker='s', color='w', label='mutant can invade', markerfacecolor='#d9d9ff', markersize=10),
        Line2D([0], [0], marker='s', color='w', label='mutant cannot invade', markerfacecolor='#ffcdcd', markersize=10)]

plt.legend(handles=legend_elements, loc='lower right', framealpha=1, fontsize='medium')


plt.axis([1e-6, 1, 1e-6, 1])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel(r'resident dispersal trait, $m$', fontsize='x-large')
plt.ylabel(r"mutant dispersal trait, $m'$", fontsize='x-large')

# plt.show()

fname = dir_results + 'pip_pretty_2.pdf'
plt.tight_layout()
plt.savefig(fname)
plt.close()

