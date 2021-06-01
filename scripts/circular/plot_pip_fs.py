# plot the PIPs for different f values on top of each other so I can compare them

import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import seaborn as sns
# sns.set_palette("colorblind")
import tol_colors as tc
cset = tc.tol_cset('muted')


# parameters
# ---

# IDs containing variation in f


#       0.01    0.02    0.04    0.06    0.08    0.09    0.095   0.10
IDV = [ 43,     53,     54,     46,     47,     48,     55,     42]
IDM = [ IDV ]

#       0.095   0.10    0.105   0.11    0.12    0.14    0.16    0.2
IDV = [ 55,     42,     56,     49,     45,     50,     51,     0]
IDM.append(IDV)

# NOTE something very weird going on with 52

# which parameter set to do the pip for (comment out ones I've done)
suffix = '_1'

# where results will be stored
dir_results = '../../results/circular/'

for idx, IDV in enumerate(IDM):

    # get the f-values corresponding to each ID
    # ---

    fname_in = dir_results + 'sing_strat' + suffix + '.csv'
    df = pd.read_csv(fname_in)
    fV = [ df.iloc[ID]['f'] for ID in IDV ]


    # plot each isocline on top of each other
    # ---

    plt.axis('scaled')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot([0, 1], [0, 1], color='black')
    plt.rc('axes', prop_cycle=plt.cycler('color', list(cset)))

    min_val = 1
    for f, ID in zip(fV, IDV):

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

        # plot log scale
        # ---

        # log scale
        plt.plot(m_iso1v, m_iso2v, label=str(f), lw=3, alpha=0.9 )


    # update minimum for the axis ranges
    min_val = 5e-4
    max_val = 0.5

    plt.xlim( (min_val, max_val) )
    plt.ylim( (min_val, max_val) )
    #plt.plot(m_iso1v, m_iso2v, '-o')
    plt.xlabel(r'resident dispersal trait, $m$', fontsize='x-large')
    plt.ylabel(r'mutant dispersal trait, $m^\prime$', fontsize='x-large')
    plt.legend(loc='best', title='f:', fontsize='large')
    #plt.show()
    plt.tight_layout()
    plt.savefig(dir_results + 'pip_fs_logscale_part' + str(idx) + '.pdf')
    plt.close()
