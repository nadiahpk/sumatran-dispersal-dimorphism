import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# model parameters
# ---

# which parameter set to find the isoclines for (comment out ones I've done)
suffix = '_1'

'''
# default
min_val = 1e-6 # lower bound for plotting
ID = 0
include_dimorphic_isocline = True
'''
# low f but divergent selection
min_val = 1e-4
ID = 62
include_dimorphic_isocline = False

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

if include_dimorphic_isocline:

    fname = dir_results + 'dimorph_isocline' + suffix + '_' + str(ID) + '.csv'
    df = pd.read_csv(fname)
    df = df.sort_values(by=['m_res1'])

    m_res1V = df['m_res1'].values
    m_res2V = df['m_res2'].values


# plot them together
# ---

# log scale
plt.axis('scaled')
plt.xscale('log')
plt.yscale('log')
plt.plot([0, 1], [0, 1], color='black', alpha=0.5)
plt.xlim( (min_val, 1) )
plt.ylim( (min_val, 1) )

# isoclines that define the TEP region
plt.plot(m_iso1V, m_iso2V, color='black')

# isocline for resident 2 (high-dispersal morph)
if include_dimorphic_isocline:
    plt.plot(m_res1V, m_res2V, lw=3, color='black')

plt.xlabel('resident 1 strategy')
plt.ylabel('resident 2 strategy')
#plt.show()
plt.tight_layout()
plt.savefig(dir_results + 'tep_logscale' + suffix + '_' + str(ID) + '.pdf')
plt.close()
