import matplotlib.pyplot as plt
import pandas as pd


# parameters
# ---

# which results to use
suffix = '_1'   # which parameter set to use
ID_default = 0  # where is the default value

# where results will be stored
dir_results = '../../results/circular/'


# get the results
# ---

# open csv file and read
fname = dir_results + 'dimorph_pop_distn' + suffix + '.csv'
df = pd.read_csv(fname)


# for starters, what's the relationship with simply m_res2?
# ---

res_default = df.iloc[ID_default]
vary_par_names = ['c', 'r', 'p_cat']

for vary_par_name in vary_par_names:

    par_val_default = res_default[vary_par_name]
    df2 = df[df[vary_par_name] != par_val_default]

    m_res2V = df2['m_res2_ss'].values
    lo_dispV = df2['lo_disp_popsize'].values
    hi_dispV = df2['hi_disp_popsize'].values

    prop_lo = lo_dispV / (lo_dispV + hi_dispV)

    plt.scatter(m_res2V, prop_lo, label=vary_par_name, alpha=0.5)

plt.ylabel('proportion of mainland population with zero dispersal')
plt.xlabel('dispersal probability of morph 2')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(dir_results + 'dimorph_pop_distn' + suffix + '_v_mres2.pdf')
plt.close()
