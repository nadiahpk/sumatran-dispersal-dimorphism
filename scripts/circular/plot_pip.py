import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# parameters
# ---

# which parameter set to do the pip for (comment out ones I've done)
'''
suffix = '_1'
# ID = 0      # default
# ID = 31     # high f = 0.25
# ID = 17     # high cost of dispersal c = 0.26
# ID = 26     # high catastrophe c = 0.001 
# ID = 42     # low f = 0.1
# ID = 41     # low f = 0.05
# ID = 43     # low f = 0.01
# ID = 44       # low f = 0.07
# ID = 45       # low f = 0.12 NOTE - not done
'''
suffix = '_2'
ID = 0


# where results will be stored
dir_results = '../../results/circular/'


# get isocline
# ---

fname = dir_results + 'isocline' + suffix + '_' + str(ID) + '.csv'
df = pd.read_csv(fname)

df = df.sort_values(by=['m_mut'])   # sort so we can plot the line nicely

# grab the isoclines
m_iso1v = df['m_res'].values
m_iso2v = df['m_mut'].values

# filter out the nans
m_iso1_2v = list(zip(m_iso1v, m_iso2v))
m_iso1_2v_filt = [ v for v in m_iso1_2v if not np.isnan(v[0]) ]
m_iso1v, m_iso2v = zip(*m_iso1_2v_filt)

# plot
# ---

# linear scale
plt.axis('auto')
plt.plot(m_iso1v, m_iso2v, 'black')
plt.xlabel('resident dispersal strategy')
plt.ylabel('mutant dispersal strategy')
plt.tight_layout()
plt.savefig(dir_results + 'pip' + suffix + '_' + str(ID) + '.pdf')
plt.close()

# log scale
min_val = min(m_iso2v)
plt.axis('scaled')
plt.xscale('log')
plt.yscale('log')
plt.plot([0, 1], [0, 1], color='black')
plt.xlim( (min_val, 1) )
plt.ylim( (min_val, 1) )
plt.plot(m_iso1v, m_iso2v, 'black')
#plt.plot(m_iso1v, m_iso2v, '-o')
plt.xlabel(r'resident dispersal trait, $m$', fontsize='x-large')
plt.ylabel(r'mutant dispersal trait, $m^\prime$', fontsize='x-large')
#plt.show()
plt.tight_layout()
plt.savefig(dir_results + 'pip_logscale' + suffix + '_' + str(ID) + '.pdf')
plt.close()
