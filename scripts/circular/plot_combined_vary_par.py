# plot how the singular strategy and derivative of the fitness gradient varies with a parameter

import matplotlib.pyplot as plt
import pandas as pd


# parameters
# ---

'''
par_name = 'f'
idx_include = list(range(37))

par_name = 'r'
idx_include = list(range(28))

par_name = 'c'
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


# get the singular strategy + divergence strength
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


# get the high-dispersal strategy
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

f2V = df2[par_name].values
m_res2V = df2['m_res2_ss'].values


# plot stacked
# ---

ylabel_coords = (-0.14, 0.5)
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(7,7))

ax1.plot(fV, m_ssV, color='black')
ax2.plot(fV, ddfit_ssV, color='black')
ax3.plot(f2V, m_res2V, color='black')
#ax2.set_yscale('symlog')
#ax2.axhline(0, color='black', alpha=0.5)

from matplotlib import ticker

ax1.set_ylabel('singular strategy,\n' + r'$m^*$', fontsize='large')
ax1.get_yaxis().set_label_coords(*ylabel_coords)
#ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2e}"))
ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#ax2.set_ylabel('divergent selection strength,\n' + r'$\left. \frac{\partial^2 w(m^\prime,m)}{\partial {m^\prime}^2} \right|_{m^\prime = m^*}$')
ax2.set_ylabel('divergent selection\nstrength', fontsize='large')
ax2.get_yaxis().set_label_coords(*ylabel_coords)
ax3.set_ylabel('high-dispersal\nmorph, ' + r'$m_H^*$', fontsize='large')
ax3.get_yaxis().set_label_coords(*ylabel_coords)

ax3.set_xlabel(xlabels[par_name], fontsize='x-large')

plt.tight_layout()
plt.savefig(dir_results + 'combined_vary_' + par_name + '.pdf')
plt.close()
