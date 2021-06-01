# TEP with quiver

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# model parameters
# ---


# which parameter set to find the isoclines for (comment out ones I've done)
suffix = '_1'
ngrid = 60
ID = 62
min_val = 1e-4 # choose higher so we can see the dynamics more clearly
max_val = 1e-0

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


# read in the dimorphic derivatives to create the quiver plot
# ---

fname = dir_results + 'dimorph_derivs_ngrid' + str(ngrid) + suffix + '_' + str(ID) + '.csv'
df = pd.read_csv(fname)
m_res1V = df['m_res1'].values
m_res2V = df['m_res2'].values
dfit1V = df['dfit1'].values
dfit2V = df['dfit2'].values
speedV = [ np.sqrt( dfit1**2 + dfit2**2 ) for dfit1, dfit2 in zip(dfit1V, dfit2V) ]

# plot them together
# ---

# log scale
plt.axis('scaled')
plt.xscale('log')
plt.yscale('log')
plt.plot([0, 1], [0, 1], color='black', alpha=0.5)
plt.xlim( (2e-4, 1e-2) )
plt.ylim( (5e-3, 3e-1) )

# isoclines
plt.plot(m_iso1V, m_iso2V, color='black')

# quiver plot of derivs in dimorphic space
#qvr = plt.quiver(m_res1V, m_res2V, dfit1V/speedV, dfit2V/speedV, speedV, scale=80, width=0.002, cmap='plasma')
qvr = plt.quiver(m_res1V, m_res2V, dfit1V/speedV, dfit2V/speedV, speedV, width=0.004, cmap='plasma')

cbar = plt.colorbar(qvr)
plt.clim(0.1,0.65)
cbar.set_label('selection gradient', fontsize='x-large')

plt.xlabel(r'low-dispersal morph trait, $m_L$', fontsize='x-large')
plt.ylabel(r'high-dispersal morph trait, $m_H$', fontsize='x-large')
#plt.show()
plt.tight_layout()
plt.savefig(dir_results + 'tep_quiver_logscale' + suffix + '_' + str(ID) + '.pdf')
plt.close()

