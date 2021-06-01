# functions used to generate stochastic environments, simulate the environment imposed by the
# resident in that stochastic environment, and calculate a mutant strategy's stochastic invasion fitness
# and its derivatives

import pickle
import numpy as np


def get_M(m, params):
    '''
    Constructs the dispersal matrix.

    Inputs:
    ---

    m, float:
        Dispersal strategy, the probability of leaving an island (or mainland) to go to a new island or mainland.

    params, dict:
        Parameter values that define the model. Reqs: h, c. Optional: f, layout

    Outputs:
    ---

    M, np matrix:
        The dispersal probabilities matrix between the mainland (idx 0) and islands (idx 1..h)

    '''

    # get useful parameters
    # ---

    h = params['h'] # number of islands
    c = params['c'] # dispersal cost

    # the proportion of the mainland that can disperse at all (i.e., only those living near the coast)
    if 'f' in params: 
        f = params['f']
    else:
        f = 1

    # edges, list of tuples of integers: List of edges between islands and mainland where mainland has index 0.
    if 'layout' in params:

        layout = params['layout']

        if layout == 'circular':

            #         mainland  --  every isle           1 -- last     isle -- isle+1
            edges = [ (0, v) for v in range(1, h+1) ] + [ (1, h) ] + [ (u, u+1) for u in range(1, h) ]

        else:
        
            edges = [ (u, v) for u in range(h+1) for v in range(u+1,h+1) ]

    else: # all_connected is the default

            edges = [ (u, v) for u in range(h+1) for v in range(u+1,h+1) ]


    # construct the migration matrix M
    # ---

    # use edges to encode the network as a matrix
    M = np.zeros( (h+1, h+1) )
    for u, v in edges:
        M[u,v] = 1; M[v,u] = 1

    degs = np.sum(M, axis=0)    # find the degree of each island

    # multiply matrix by migration probability
    M = m*M

    # the first column is receipts from the mainland, so multiply the first column by the proportion that can leave the mainland, f
    M[:,0] = f*M[:,0]

    # destination chosen at random, so divide each column through by its degree
    for col in range(h+1):
        M[:,col] = M[:,col]/degs[col]

    # populate the diagonals
    M[0,0] = 1-m*f              # proportion who leave mainland = m*f
    for idx in range(1,h+1):    # proportion who leave island = m
        M[idx,idx] = 1-m


    # put in dispersal costs
    # ---

    # dispersal cost applies everywhere except the diagonal (staying home)
    C = (1-c) * ( np.ones((h+1,h+1), dtype=int) - np.identity(h+1, dtype=int) ) + np.identity(h+1, dtype=int)
    M = C * M

    return(M)


def generate_catastrophes(params, tlimit, fname=None):
    '''
    Generate a list of timesteps and the islands that had a catastrophe in
    that timestep

    Inputs:
    ---

    params, dict:
        Parameter values that define the model. Reqs: h, p_cat.

    tlimit, int:
        Number of timesteps to simulate catastrophes for

    fname, string:
        If specified, the name of the pickle file where they list of catastrophes will be stored

    Outputs:
    ---

    catastrophD, dictionary of tuples:
        The timestep (key) and list of small-island indexes (vals) that had a catastrophe

    Creates:
    ---

    fname
        A pickle file containing catastropheD. Saves us redoing the simulation for parameter-value exploration.
    '''

    h = params['h']         # number of small islands
    p_cat = params['p_cat'] # probability of a catastrophe on each island

    catastropheD = dict()   # catastrophes stored as a dictionary: {time: (list of islands having a catastrophe)}

    for t in range(tlimit):

        # how many islands have a catastrophe in year t?
        no_cat = h
        while no_cat == h: # disallow all small islands to have a catastrophe in the same year
            no_cat = np.random.binomial(h, p_cat)

        # if any islands had a catastrophe, store them
        if no_cat > 0:

            # find out which ones they were
            cat_idxs = tuple(np.random.choice(range(1,h+1), size=no_cat, replace=False))

            # append them to our dictionary
            catastropheD[t] = cat_idxs

    # if a file name was specified, write to a pickle file
    if not (fname is None):

        f = open(fname, 'wb')
        # a string explaining the pickle file
        ss  = 'A simulation of the stochastic environment as varying carrying capacities K.\n'
        ss += 'Created by generate_Kts().\n'
        ss += 'This pickle file contains the following:\n'
        ss += '0. ss, string: this string you are reading now.\n'
        ss += '1. catastropheD, dictionary of tuples: the timesteps (key) and island indexes (val) that had a catastrophe.\n'
        ss += '2. params, dictionary of floats: the parameter values used for this environment simulation.\n'
        pickle.dump( ss, f )
        pickle.dump( catastropheD, f )
        pickle.dump( params, f )
        f.close()

    return(catastropheD)


def sample_wsV2(m_res, params, catastropheD, tlimit, initial_nT = None, tburnin=None, return_nT = False):
    '''
    Given a time-series of catastrophes,
    simulate the resident dynamics to sample the growth rate matrices W.


    Inputs:
    ---

    m_res, float:
        Resident dispersal strategy.

    params, dict:
        Parameter values that define the model. Reqs: r, h, KL, KS_mean.

    catastrophD, dictionary of tuples:
        The timestep (key) and list of small-island indexes (vals) that had a catastrophe

    tlimit, int:
        How many timesteps to simulate the resident dynamics for.

    initial_nT, np array length h+1
        An initial population structure vector that can be used instead of the burnin.
        Default is carrying capacities

    tburnin, int:
        Number of timesteps to burn in the resident dynamics before recording growth rate matrices.
        Default = tlimit

    return_nT, boolean:
        Flag for if we want the final nT returned


    Outputs:
    ---

    wsV, list of lists:
        The diagonal of the W matrix at each timestep (rows)
        
    '''

    # get needed parameter values
    r = params['r']             # intrinsic growth rate
    h = params['h']             # number of small islands
    KL = params['KL']           # carrying capacity on the large island
    KS_mean = params['KS_mean'] # carrying capacity on the small islands

    # default carrying capacities
    Ks = [KL] + [KS_mean]*h

    # set the burnin parameter to default if needed
    if tburnin is None:
        tburnin = tlimit

    # get resident species' dispersal matrix
    M = get_M(m_res, params)


    # burn in resident dynamics
    # ---

    # choose the initial population sizes on each island
    if initial_nT is None:
        nT = np.array([KL] + [KS_mean]*h )
    else:
        nT = initial_nT

    for t in range(tburnin):

        # if there were catastrophes, set those carrying capacities to zero
        if t in catastropheD:
            cat_idxs = catastropheD[t]
            Kts = [ 0 if i in cat_idxs else Ks[i] for i in range(h+1) ]
        else:
            Kts = Ks

        # construct the matrix A and project it
        ws = [ 0 if Kts[i] == 0 else np.exp( r*(1 - ( nT[i] / Kts[i] )) ) for i in range(h+1) ]
        W = np.diag( ws )
        A = M @ W
        nT = A @ nT.transpose()

    # now have nT near the steady state
    

    # sample W matrix in the face of stochasticity and resident dynamics
    # ---

    wsV = list()

    for t in range(tlimit):

        # if there were catastrophes, set those carrying capacities to zero
        if t in catastropheD:
            cat_idxs = catastropheD[t]
            Kts = [ 0 if i in cat_idxs else Ks[i] for i in range(h+1) ]
        else:
            Kts = Ks

        # construct the matrix A and project it
        ws = [ 0 if Kts[i] == 0 else np.exp( r*(1 - ( nT[i] / Kts[i] )) ) for i in range(h+1) ]
        wsV.append(ws) # store for later
        W = np.diag( ws )
        A = M @ W
        nT = A @ nT.transpose()

    if return_nT:
        res = (wsV, nT)
    else:
        res = wsV

    return(res)


def calc_fitness2(m, wsV, params, initial_nT = None, tburnin = None):
    '''
    Given a stochastic resident dynamics, find the invasion fitness of a mutant dispersal strategy.

    Inputs:
    ---

    m, float:
        Mutant dispersal strategy.

    wsV, list of lists:
        The diagonal of the W matrix at each timestep (rows) imposed by the stochastic resident dynamics

    params, dict:
        Parameter values that define the model.

    initial_nT, np array length h+1
        An initial normalised population structure vector that can be used instead of (or in addition to) the burnin.
        Default is carrying capacities.

    tburnin, int:
        Number of timesteps to burn in the mutant dynamics before recording growth rates.
        Default is the length of the wsV vector.

    Outputs:
    ---

    loglambda, float:
        An estimate of the stochastic invasion fitness of the mutant.

    '''

    # get needed parameter values
    # ---

    # mutant's dispersal matrix
    M = get_M(m, params)

    # set the burnin parameter to default if needed
    if tburnin is None:
        tburnin = len(wsV)


    # set the initial population structure if needed
    # ---

    if initial_nT is None:
        h = params['h']             # number of small islands
        KL = params['KL']           # carrying capacity on the large island
        KS_mean = params['KS_mean'] # carrying capacity on the small islands
        n = np.array([KL] + [KS_mean]*h )
        n = n/np.sum(n)
    else:
        n = initial_nT


    # burn in the mutant dynamics
    # ---

    for ws in wsV[:tburnin]:

        # construct matrix A and project
        W = np.diag( ws )
        A = M @ W
        n = A @ n.transpose()
        N = np.sum(n)
        n = n / N # normalise for next timestep

    # now we have the mutant population structure somewhere near the steady state


    # estimate stochastic invasion fitness
    # ---

    grwthV = list() # a place to store mutant growth rates
    for ws in wsV:

        # construct matrix A and project
        W = np.diag( ws )
        A = M @ W
        n = A @ n.transpose()
        N = np.sum(n)
        n = n / N # normalise for next timestep

        # store mutant growthrates
        grwth = np.log(N)
        grwthV.append(grwth)


    # estimate longterm stochastic invasion finess loglambda using Tuljaparkur
    loglambda = np.mean(grwthV)

    return(loglambda)



def calc_dfitness(m, wsV, params, delta_m=1e-9, tburnin=1000):
    '''
    Given a stochastic resident dynamics, find the fitness gradient of an invading mutant's dispersal strategy.

    Inputs:
    ---

    m, float:
        The mutant's dispersal strategy the probability of leaving an island (or mainland) to go to a new island or mainland.

    wsV, list of lists:
        The diagonal of the W matrix at each timestep (rows) imposed by the stochastic resident dynamics

    params, dict:
        Parameter values that define the model.

    delta_m, float:
        The step size used to estimate the gradient. The default was chosen using the heuristic
        sqrt(machine error) * approx value

    tburnin, int:
        Number of timesteps to burn in the mutant dynamics before recording growth rates,
        which are averaged to find the stochastic invasion fitness.

    Outputs:
    ---

    fitgrad, float:
        The mutant fitness gradient at m.
    '''

    # mutant fitness gradient given resident sample
    if m < delta_m:
        md_lo = m
        md_hi = m+2*delta_m
    elif m > 1-delta_m:
        md_hi = m
        md_lo = m-2*delta_m
    else:
        md_lo = m - delta_m
        md_hi = m + delta_m
    

    Md_lo = get_M(md_lo, params)
    fit_lo = calc_fitness(Md_lo, wsV, params, tburnin=tburnin)

    Md_hi = get_M(md_hi, params)
    fit_hi = calc_fitness(Md_hi, wsV, params, tburnin=tburnin)

    fitgrad = (fit_hi - fit_lo) / (2*delta_m)

    return(fitgrad)

def calc_ddfitness(m, wsV, params, delta_m=1e-9, tburnin=1000):
    '''
    Given a stochastic resident dynamics, find the derivative of the mutant fitness gradient.
    This is used at a singular strategy to evaluate whether or not the strategy is evolutionarily stable
    or unstable. For unstable singular strategies, it is a measure of the strength of divergent selection.

    Inputs:
    ---

    m, float:
        The mutant's dispersal strategy the probability of leaving an island (or mainland) to go to a new island or mainland.

    wsV, list of lists:
        The diagonal of the W matrix at each timestep (rows) imposed by the stochastic resident dynamics

    params, dict:
        Parameter values that define the model.

    delta_m, float:
        The step size used to estimate the gradient. The default was chosen using the heuristic
        sqrt(machine error) * approx value

    tburnin, int:
        Number of timesteps to burn in the mutant dynamics before recording growth rates,
        which are averaged to find the stochastic invasion fitness.

    Outputs:
    ---

    dfitgrad, float:
        The derivative of the mutant fitness gradient at m.
    '''

    # mutant fitness gradient given resident sample
    if m < delta_m:
        md_lo = m
        md_hi = m+2*delta_m
    elif m > 1-delta_m:
        md_hi = m
        md_lo = m-2*delta_m
    else:
        md_lo = m - delta_m
        md_hi = m + delta_m
    

    dfit_lo = calc_dfitness(md_lo, wsV, params, tburnin=tburnin)
    dfit_hi = calc_dfitness(md_hi, wsV, params, tburnin=tburnin)

    dfitgrad = (dfit_hi - dfit_lo) / (2*delta_m)

    return(dfitgrad)

# ============================================================================

# older functions not used for final results (keep a copy for older scripts just in case)

# replaced by generate_catastrophes()
def generate_Kts(params, tlimit, fname=None):
    '''
    Generate a time-series of the stochastic environment, the carrying capacities K_i on each island i.
    This is a simulation where the mainland (i=0) has constant K, but the small islands suffer catastrophes
    with probability p such that K=0. Note that the possibility of all small islands suffering a catastrophe
    in the same timestep is excluded.

    # old function

    Inputs:
    ---

    params, dict:
        Parameter values that define the model.

    tlimit, int:
        Number of timesteps to simulate catastrophes for

    fname, string:
        If specified, the name of the pickle file where they list of catastrophes will be stored

    Outputs:
    ---

    Kts, list of lists:
        Carrying capacities on each island (cols) at each timestep (rows)

    Creates:
    ---

    Kts[suffix].pkl
        A pickle file containing the Kts. Saves us redoing the simulation for parameter-value exploration.
    '''

    h = params['h']
    p_cat = params['p_cat']
    KL = params['KL']
    KS_mean = params['KS_mean']

    Kts = list()
    for t in range(tlimit):

        # how many islands had a catastrophe in year t?
        no_cat = h
        while no_cat == h: # we're disallowing all h small islands to have a catastrophe in the same year
            no_cat = np.random.binomial(h, p_cat)

        # choose which islands to have the catastrophes
        cat_idxs = np.random.choice(range(h), size=no_cat, replace=False)

        # construct environments
        Kt = [ KL ] + [ 0 if i in cat_idxs else KS_mean for i in range(h) ]
        Kts.append(Kt)

    if not (fname is None):

        f = open(fname, 'wb')
        # a string explaining the pickle file
        ss  = 'A simulation of the stochastic environment as varying carrying capacities K.\n'
        ss += 'Created by generate_Kts().\n'
        ss += 'This pickle file contains the following:\n'
        ss += '0. ss, string: this string you are reading now.\n'
        ss += '1. Kts, list of lists: the carrying capacities on each island (cols) at each timestep (rows).\n'
        ss += '2. params, dictionary of floats: the parameter values used for this environment simulation.\n'
        pickle.dump( ss, f )
        pickle.dump( Kts, f )
        pickle.dump( params, f )
        f.close()

    return(Kts)


def sample_wsV(M, Kts, params, tburnin=1000, tlimit=10000):
    '''
    Given a time-series of the stochastic environment, 
    simulate the resident dynamics to sample the growth rate matrices W.

    Inputs:
    ---

    M, np matrix:
        Resident dispersal matrix.

    Kts, list of lists:
        Carrying capacities on each island (cols) at each timestep (rows)

    tburnin, int:
        Number of timesteps to burn in the resident dynamics before recording growth rate matrices.

    tlimit, int:
        How many timesteps to simulate the resident dynamics for.

    Outputs:
    ---

    wsV, list of lists:
        The diagonal of the W matrix at each timestep (rows)
        
    '''

    r = params['r']
    h = params['h']

    # construct initial population vector
    KL = params['KL']
    KS_mean = params['KS_mean']
    nT = np.array([KL] + [KS_mean]*h )

    # burn in resident dynamics

    for t in range(tburnin):

        # construct matrix A and project
        ws = [ 0 if Kts[t][i] == 0 else np.exp( r*(1 - ( nT[i] / Kts[t][i] )) ) for i in range(h+1) ]
        W = np.diag( ws )
        A = M @ W
        nT = A @ nT.transpose()

    # now have nT near the steady state
    # project the resident population

    wsV = list()
    for t in range(tlimit):

        # construct matrix A and project
        ws = [ 0 if Kts[t][i] == 0 else np.exp( r*(1 - ( nT[i] / Kts[t][i] )) ) for i in range(h+1) ]
        wsV.append(ws) # store for later
        W = np.diag( ws )
        A = M @ W
        nT = A @ nT.transpose()

    return(wsV)


def calc_fitness(M, wsV, params, tburnin=100):
    '''
    Given a stochastic resident dynamics, find the invasion fitness of a mutant dispersal strategy.

    Inputs:
    ---

    M, np matrix:
        The dispersal probabilities matrix of the mutant between the mainland (idx 0) and islands (idx 1..size)
        This defines the mutant dispersal strategy.

    wsV, list of lists:
        The diagonal of the W matrix at each timestep (rows) imposed by the stochastic resident dynamics

    params, dict:
        Parameter values that define the model.

    tburnin, int:
        Number of timesteps to burn in the mutant dynamics before recording growth rates.

    Outputs:
    ---

    loglambda, float:
        An estimate of the stochastic invasion fitness of the mutant.
    '''


    # initial vector
    h = params['h']
    KL = params['KL']
    KS_mean = params['KS_mean']
    n = np.array([KL] + [KS_mean]*h )
    n = n/np.sum(n)

    for ws in wsV[:tburnin]:

        # construct matrix A and project
        W = np.diag( ws )
        A = M @ W
        n = A @ n.transpose()
        N = np.sum(n)
        n = n / N # normalise for next timestep

    # now have n near steady state
    # sample mutant fitness

    grwthV = list() # a place to store mutant growth rates

    for ws in wsV: # for each stochastic environment-determined year

        # construct matrix A and project
        W = np.diag( ws )
        A = M @ W
        n = A @ n.transpose()
        N = np.sum(n)
        n = n / N # normalise for next timestep

        # store mutant growth rates
        grwth = np.log(N)
        grwthV.append(grwth)

    # estimate longterm stochastic invasion fitness loglambda using Tuljaparkur
    loglambda = np.mean(grwthV)

    return(loglambda)
