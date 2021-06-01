# plot how the dimorphic steady state varies with a parameter

import matplotlib.pyplot as plt
import pandas as pd


# parameters
# ---

# comment out parameters done

par_name = 'c'
idx_include = None

'''
par_name = 'r'
idx_include = None

par_name = 'p_cat'
idx_include = None
'''

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
fname = dir_results + 'dimorph_steady_state' + suffix + '.csv'
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
    m_res2V = df2['m_res2_ss'].values

else:
    fV = df2[par_name].values[idx_include]
    m_res2V = df2['res2_ss'].values[idx_include]


# plot 
# ---

plt.plot(fV, m_res2V, color='black')
plt.ylabel('high-dispersal morph strategy')
plt.xlabel(xlabels[par_name])
plt.tight_layout()
plt.savefig(dir_results + 'dimorph_vary_' + par_name + '.pdf')
plt.close()
