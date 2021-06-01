# plot how the singular strategy and derivative of the fitness gradient varies with a parameter

import matplotlib.pyplot as plt
import pandas as pd


# parameters
# ---

'''
par_name = 'f'
idx_include = list(range(37))

par_name = 'c'
idx_include = list(range(28))

par_name = 'r'
idx_include = list(range(28))
'''

par_name = 'p_cat'
idx_include = None

# xlabels
xlabels = {
        'f': r'proportion mainland disperse, $f$',
        'c': r'dispersal cost, $c$',
        'r': r'intrinsic growth rate, $r$',
        'p_cat': r'probability of catastrophe, $p_c$',
        }

# which results to use
suffix = '_1'   # which parameter set to use
ID_default = 0  # where is the default value

# where results will be stored
dir_results = '../../results/circular/'


# get the results
# ---

# open csv file and read
fname = dir_results + 'sing_strat' + suffix + '.csv'
df = pd.read_csv(fname)

# find all the rows where par_name is not at its default value
res_default = df.iloc[ID_default]
par_val_default = res_default[par_name]
df2 = df[df[par_name] != par_val_default]

# append default and sort
df2 = df2.append(df.iloc[ID_default])
df2 = df2.sort_values(by=[par_name])   # sort so we can plot the line nicely

if idx_include is None:

    fV = df2[par_name].values
    m_ssV = df2['m_ss'].values
    ddfit_ssV = df2['ddfit_ss'].values

else:
    fV = df2[par_name].values[idx_include]
    m_ssV = df2['m_ss'].values[idx_include]
    ddfit_ssV = df2['ddfit_ss'].values[idx_include]


# plot stacked
# ---

fig, (ax1, ax2) = plt.subplots(2, sharex=True)

ax1.plot(fV, m_ssV, color='black')
ax2.plot(fV, ddfit_ssV, color='black')
#ax2.set_yscale('symlog')
#ax2.axhline(0, color='black', alpha=0.5)

ax1.set_ylabel('singular strategy')
ax2.set_ylabel('divergent selection strength')
ax2.set_xlabel(xlabels[par_name])

plt.tight_layout()
plt.savefig(dir_results + 'vary_' + par_name + '.pdf')
plt.close()
