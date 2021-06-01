# TEP with quiver

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# model parameters
# ---


# which parameter set to find the isoclines for (comment out ones I've done)
suffix = '_1'
min_val = 1e-6 # for plot
ngrid = 20
# ID = 0
ID = 17

# where results will be stored
dir_results = '../../results/circular/'


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

m_iso1V, m_iso2V = zip(*region)


# read in the dimorphic isocline to show the maximum for resident 2
# ---

fname = dir_results + 'dimorph_isocline' + suffix + '_' + str(ID) + '.csv'
df = pd.read_csv(fname)
df = df.sort_values(by=['m_res1'])

m_dimiso1V = df['m_res1'].values
m_dimiso2V = df['m_res2'].values


# read in the dimorphic derivatives to create the quiver plot
# ---

fname = dir_results + 'dimorph_derivs_ngrid' + str(ngrid) + suffix + '_' + str(ID) + '.csv'
df = pd.read_csv(fname)
m_res1V = df['m_res1'].values
m_res2V = df['m_res2'].values
dfit1V = df['dfit1'].values
dfit2V = df['dfit2'].values
speedV = [ np.sqrt( dfit1**2 + dfit2**2 ) for dfit1, dfit2 in zip(dfit1V, dfit2V) ]

# maximum speeds
# ID = 0,  max = 0.5875829787576762
# ID = 17, max = 0.6691351767205614

# plot them together
# ---

# log scale
plt.axis('scaled')
plt.xscale('log')
plt.yscale('log')
plt.plot([0, 1], [0, 1], color='black', alpha=0.5)
plt.xlim( (min_val, 1) )
plt.ylim( (min_val, 1) )

# isoclines
plt.plot(m_iso1V, m_iso2V, color='black')

# dimorphic isocline for resident 2
plt.plot(m_dimiso1V, m_dimiso2V, lw=3, color='black')

# quiver plot of derivs in dimorphic space
qvr = plt.quiver(m_res1V, m_res2V, dfit1V, dfit2V, speedV, cmap='plasma')

cbar = plt.colorbar(qvr)
plt.clim(0.1,0.65)
cbar.set_label('selection gradient', fontsize='x-large')


# plot a point at the singular strategy
# ---

plt.yscale('log')
plt.xscale('log')

# read in the information from the sing_strat_x.csv file
fname_in = dir_results + 'sing_strat' + suffix + '.csv'
df = pd.read_csv(fname_in)
ss_res = df.iloc[ID]                # the particular row we want
m_ss = ss_res['m_ss']               # singular strategy

plt.scatter([m_ss], [m_ss], color='black', zorder=9999)


# plot decorations and save
# ---

plt.xlabel(r'low-dispersal morph trait, $m_L$', fontsize='x-large')
plt.ylabel(r'high-dispersal morph trait, $m_H$', fontsize='x-large')
#plt.show()
plt.tight_layout()
plt.savefig(dir_results + 'tep_quiver_pretty_logscale' + suffix + '_' + str(ID) + '.pdf')
plt.close()

