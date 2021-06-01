import numpy as np

import sys
#sys.path.insert(0,'../../functions') # so I can import the functions
from funcs import get_M

def equilib_mutant(m, wsV, params):

    # extract the needful
    # ---

    h = params['h']
    tol_mut = params['tol_mut']
    KL = params['KL']
    KS_mean = params['KS_mean']


    # create the dispersal matrix
    # ---

    M = get_M(m, params)


    # initialise the mutant population structure
    # ---

    if m == 0:

        # I know that this one equilibriates with only individuals on the mainland
        n = np.array([1] + [0]*h)

    else:

        # initialise at the carrying capacity
        n = np.array([KL] + [KS_mean]*h )
        n = n/np.sum(n)


    # equilibriate the population structure of the mutant
    # ---

    prev_n = n
    dist = 1+tol_mut
    while dist > tol_mut:

        # one cycle of catastrophes
        for ws in wsV:

            # construct matrix A and project
            W = np.diag( ws )
            A = M @ W
            n = A @ n.transpose()
            N = np.sum(n)
            n = n / N # normalise for next timestep

        dist = np.linalg.norm(n-prev_n)
        prev_n = n

    return(n)

def calc_fitness(m, wsV, params):


    # extract the needful
    # ---

    h = params['h']
    tol_mut = params['tol_mut']
    KL = params['KL']
    KS_mean = params['KS_mean']


    # create the dispersal matrix
    # ---

    M = get_M(m, params)


    # initialise the mutant population structure
    # ---

    if m == 0:

        # I know that this one equilibriates with only individuals on the mainland
        n = np.array([1] + [0]*h)

    else:

        # initialise at the carrying capacity
        n = np.array([KL] + [KS_mean]*h )
        n = n/np.sum(n)


    # equilibriate the population structure of the mutant
    # ---

    prev_n = n
    dist = 1+tol_mut
    while dist > tol_mut:

        # one cycle of catastrophes
        for ws in wsV:

            # construct matrix A and project
            W = np.diag( ws )
            A = M @ W
            n = A @ n.transpose()
            N = np.sum(n)
            n = n / N # normalise for next timestep

        dist = np.linalg.norm(n-prev_n)
        prev_n = n

    # now equilibriated, calculate growth
    # ---

    grwthV = list() # a place to store mutant growth rates

    # one cycle of catastrophes is all that's needed because it loops
    for ws in wsV:

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

def sample_wsV(m, params, return_nT = False):


    # extract some needed parameters
    # ---

    r = params['r']
    h = params['h']
    tol_res = params['tol_res']
    KL = params['KL']
    KS_mean = params['KS_mean']
    p_cat = params['p_cat'] 


    # calculate time between catastrophes
    # ---

    tau = 1/p_cat       # expected time to the next catastrophe on one island (geometric distn)
    tt = int(tau / h)   # assume islands go in a cycle, the time between any island catastrophe


    # construct initial population vector
    # ---

    if m == 0:

        # I know that this one equilibriates with only individuals on the mainland
        nT = np.array([KL] + [0]*h)

    else:

        # initialise at the carrying capacity
        nT = np.array([KL] + [KS_mean]*h)


    # other useful things
    # ---

    Kt = [KL] + [KS_mean]*h     # carrying capacity when there are no catastrophes
    M = get_M(m, params)        # dispersal matrix


    # equilibriate resident dynamics
    # ---

    prev_nT = nT
    dist = 1+tol_res
    while dist > tol_res:

        # nV = list() 
        # one cycle of catastrophes
        for hi in range(h):

            # first generation drops carrying capacity of island hi to 0
            ws = [ np.exp( r*(1 - ( nT[i] / Kt[i] )) ) for i in range(h+1) ]
            ws[hi+1] = 0

            # step
            W = np.diag( ws )
            A = M @ W
            nT = A @ nT.transpose()
            #nV.append(list(nT))

            # remaining tt-1 generations use full carrying capacity
            for t in range(tt-1):

                ws = [ np.exp( r*(1 - ( nT[i] / Kt[i] )) ) for i in range(h+1) ]
                W = np.diag( ws )
                A = M @ W
                nT = A @ nT.transpose()
                #nV.append(list(nT))

        dist = np.linalg.norm(nT-prev_nT)
        prev_nT = nT


    # gather one cycle of resident dynamics at the equilibrium
    # ---

    #nV = list() 
    wsV = list()

    # one cycle of catastrophes
    for hi in range(h):

        # first generation drops carrying capacity of island hi to 0
        ws = [ np.exp( r*(1 - ( nT[i] / Kt[i] )) ) for i in range(h+1) ]
        ws[hi+1] = 0
        wsV.append(ws) # store for later

        # step
        W = np.diag( ws )
        A = M @ W
        nT = A @ nT.transpose()
        #nV.append(list(nT))

        # remaining tt-1 generations use full carrying capacity
        for t in range(tt-1):

            ws = [ np.exp( r*(1 - ( nT[i] / Kt[i] )) ) for i in range(h+1) ]
            W = np.diag( ws )
            wsV.append(ws) # store for later
            A = M @ W
            nT = A @ nT.transpose()
            #nV.append(list(nT))

    '''
    # plot to check
    # ---

    nM = np.array(nV)
    for i in range(h+1):
        plt.plot(nM[:,i])
    plt.show() # looks good

    wsM = np.array(wsV)
    for i in range(h+1):
        plt.plot(wsM[:,i])
    plt.show() # looks good
    '''

    if return_nT:
        res = (wsV, nT)
    else:
        res = wsV

    return(res)


def calc_dfitness(m, wsV, params):
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

    Outputs:
    ---

    fitgrad, float:
        The mutant fitness gradient at m.
    '''

    # the step I take depends on if there's room for it
    # but generally I'm using the middle method

    delta_m = params['delta_m']

    if m < delta_m:
        md_lo = m
        md_hi = m+2*delta_m
    elif m > 1-delta_m:
        md_hi = m
        md_lo = m-2*delta_m
    else:
        md_lo = m - delta_m
        md_hi = m + delta_m
    

    # gradient approximated 

    fit_lo = calc_fitness(md_lo, wsV, params)
    fit_hi = calc_fitness(md_hi, wsV, params)
    fitgrad = (fit_hi - fit_lo) / (2*delta_m)

    return(fitgrad)

def calc_ddfitness(m, wsV, params):
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

    # the step I take depends on if there's room for it
    # but generally I'm using the middle method

    delta_m = params['delta_m']*1e3 # NOTE: somewhat heuristic choice

    if m < delta_m:
        md_lo = m
        md_hi = m+2*delta_m
    elif m > 1-delta_m:
        md_hi = m
        md_lo = m-2*delta_m
    else:
        md_lo = m - delta_m
        md_hi = m + delta_m
    

    dfit_lo = calc_dfitness(md_lo, wsV, params)
    dfit_hi = calc_dfitness(md_hi, wsV, params)
    dfitgrad = (dfit_hi - dfit_lo) / (2*delta_m)

    return(dfitgrad)

def sample_wsV_dimorph(m_res1, m_res2, params, return_nT = False, nL = None):
    '''
    nL, initial guess for population structure
    '''

    # extract some needed parameters
    # ---

    r = params['r']
    h = params['h']
    tol_res_dimorph = params['tol_res_dimorph']
    KL = params['KL']
    KS_mean = params['KS_mean']
    p_cat = params['p_cat'] 


    # useful things
    # ---

    Kt = [KL] + [KS_mean]*h     # carrying capacity when there are no catastrophes

    M1 = get_M(m_res1, params)  # each morph's dispersal matrix
    M2 = get_M(m_res2, params)
    ML = [M1, M2]

    # calculate time between catastrophes
    tau = 1/p_cat               # expected time to the next catastrophe on one island (geometric distn)
    tt = int(tau / h)       # assume islands go in a cycle, the time between any island catastrophe


    # initialise population sizes on mainland and islands
    # ---

    if nL is None:

        # initialise with the smaller morph as being at carrying capacity on the mainland
        n_lo_m = np.array([KL] + [0]*h)
        n_hi_m = np.array([0] + [KS_mean]*h)
        nL = [n_lo_m, n_hi_m] if m_res1 < m_res2 else [n_hi_m, n_lo_m] 


    # equilibriate resident dynamics
    # ---

    prev_N = np.array(list(nL[0]) + list(nL[1]))
    prev_N = prev_N / sum(prev_N)
    dist = 1+tol_res_dimorph

    while dist > tol_res_dimorph:

        # one cycle of catastrophes

        for hi in range(h):

            # one step with a 0 carrying capacity for island h_i

            ws = [ np.exp( r*(1 - ( sum( n[i] for n in nL ) / Kt[i] )) ) for i in range(h+1) ]
            ws[hi+1] = 0

            nL_new = list()
            for n, M in zip(nL, ML):

                W = np.diag( ws )
                A = M @ W
                n_new = A @ n.transpose()
                nL_new.append(n_new)

            nL = nL_new

            # remaining tt-1 generations use full carrying capacity

            for t in range(tt-1):

                ws = [ np.exp( r*(1 - ( sum( n[i] for n in nL ) / Kt[i] )) ) for i in range(h+1) ]

                # step
                nL_new = list()
                for n, M in zip(nL, ML):

                    W = np.diag( ws )
                    A = M @ W
                    n_new = A @ n.transpose()
                    nL_new.append(n_new)

                nL = nL_new

        # check distance of population from previous (how far from eqm)

        N = np.array(list(nL[0]) + list(nL[1]))
        N = N / sum(N)
        dist = np.linalg.norm(N - prev_N)
        prev_N = N


    # gather one cycle of resident dynamics at the equilibrium
    # ---

    wsV = list()
    for hi in range(h):

        # one step with a 0 carrying capacity for island h_i

        ws = [ np.exp( r*(1 - ( sum( n[i] for n in nL ) / Kt[i] )) ) for i in range(h+1) ]
        ws[hi+1] = 0
        wsV.append(ws)

        nL_new = list()
        for n, M in zip(nL, ML):

            W = np.diag( ws )
            A = M @ W
            n_new = A @ n.transpose()
            nL_new.append(n_new)

        nL = nL_new

        # remaining tt-1 generations use full carrying capacity

        for t in range(tt-1):

            ws = [ np.exp( r*(1 - ( sum( n[i] for n in nL ) / Kt[i] )) ) for i in range(h+1) ]
            wsV.append(ws)

            # step
            nL_new = list()
            for n, M in zip(nL, ML):

                W = np.diag( ws )
                A = M @ W
                n_new = A @ n.transpose()
                nL_new.append(n_new)

            nL = nL_new


    # return results
    # ---

    if return_nT:
        res = (wsV, nL)
    else:
        res = wsV

    return(res)

