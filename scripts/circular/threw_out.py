# from find_dimorph_isocline.py

# looking at mutations in resident 1 (low dispersal morph)
if False:

    # check the 0,1 point
    # ---

    m_res1 = 0
    m_res2 = 1

    wsV, nL = sample_wsV_dimorph(m_res1, m_res2, params, return_nT=True, nL=None)

    dfit1 = calc_dfitness(m_res1, wsV, params) # --- takes too long

    # check the 0,0.26 point, where the isocline is some reasonable distance out from the axis
    # ---

    m_res2 = 0.26
    wsV, nL = sample_wsV_dimorph(m_res1, m_res2, params, return_nT=True, nL=None)

    dfit1 = calc_dfitness(m_res1, wsV, params) # --- it's negative

    # now try some distance out, where previously I found it would go positive
    m_res1 = 0.00005 # half the way to the isocline
    wsV, nL = sample_wsV_dimorph(m_res1, m_res2, params, return_nT=True, nL=None)
    # [array([433.32779939,  22.10235444,  22.10235455,  22.10235455, 22.10235444,  22.10235432]),
    #  array([65.83291608,  7.79066622,  7.7906662 ,  7.7906662 ,  7.79066622, 7.79066623])]


    dfit1 = calc_dfitness(m_res1, wsV, params) # --- it's positive, so there is some push away from the isocline there

    # try on the isocline
    m_res1, m_res2 = (0.00012075549556944323, 0.26530612244897955)
    wsV, nL = sample_wsV_dimorph(m_res1, m_res2, params, return_nT=True, nL=None)
    # [array([499.70725917,  29.98420696,  29.98420696,  29.98420696, 29.98420696,  29.98420696]),
    #  array([0.22868724, 0.02702231, 0.02702231, 0.02702231, 0.02702231, 0.02702231])]

    dfit1 = calc_dfitness(m_res1, wsV, params) # --- it's also positive, so that implies it will tip it ouside of the region? 0.06460391318279933

    # is the issue that I am not running the dimorphic population to steady state properly?
    # because one of them should have gone extinct for that combination
    params['tol_res_dimorph'] = 1e-10
    wsV, nL = sample_wsV_dimorph(m_res1, m_res2, params, return_nT=True, nL=nL) # takes too long



'''


# define the space to grid over
# ---

# flip over the part below the diagonal, so we're defining the dimorphism region
region = [ (m_mut, m_res) if m_res > m_mut else (m_res, m_mut) for m_res, m_mut in zip(m_resV, m_mutV) ]
m_iso1V, m_iso2V = zip(*region)

# find the extents
min_m_iso1 = min(m_iso1V); max_m_iso1 = max(m_iso1V)
min_m_iso2 = min(m_iso2V); max_m_iso2 = max(m_iso2V)

# define a grid

eps = 1e-5 # to get it off the edges a bit
ngrid = 50

if False:

    # on log scale
    pwr_res1V = np.linspace( np.log10(min_m_iso1+eps), np.log10(max_m_iso1-eps), ngrid )
    pwr_res2V = np.linspace( np.log10(min_m_iso2+eps), np.log10(max_m_iso2-eps), ngrid )

    res1_gridV = [ 10**pwr for pwr in pwr_res1V ]
    res2_gridV = [ 10**pwr for pwr in pwr_res2V ]

else:

    res1_gridV = np.linspace( min_m_iso1+eps, max_m_iso1-eps, ngrid )
    res2_gridV = np.linspace( min_m_iso2+eps, max_m_iso2-eps, ngrid )

# find which points in the grid are within the dimorphism space

tep_polygon = Polygon(region)

m_res1_res2V = list()
for res1 in res1_gridV:
    for res2 in res2_gridV:

        point = Point(res1, res2)

        if tep_polygon.contains(point):

            m_res1_res2V.append( (res1, res2) )


if False:
    # plot to check
    plt.plot(m_iso1V, m_iso2V);
    m_res1V, m_res2V = zip(*m_res1_res2V)
    plt.scatter(m_res1V, m_res2V)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.show()


# find the derivatives at each point
# ---

dfit1V = list()
dfit2V = list()
speedV = list()
cnt = 0
for m_res1, m_res2 in m_res1_res2V:

    print(cnt)
    cnt += 1
    # fake data
    #dfit1 = -1 + 2*np.random.rand()
    #dfit2 = -1 + 2*np.random.rand()
    wsV = sample_wsV_dimorph(m_res1, m_res2, params, False)
    dfit1 = calc_dfitness(m_res1, wsV, params) 
    dfit2 = calc_dfitness(m_res2, wsV, params) 

    speed = np.sqrt(dfit1**2 + dfit2**2)

    dfit1V.append( dfit1 )
    dfit2V.append( dfit2 )
    speedV.append( speed )


# save to file
# ---

fname = dir_results + 'tep.pkl'
f = open(fname, 'wb')

# a string explaining the pickle file
ss  = 'Mutant fitness derivatives in dimorphic region.\n'
ss  += '1. ss, this string.\n'
ss  += '2. region, polygon defining the dimorphic region.\n'
ss  += '3. res1_gridV, full grid in the region.\n'
ss  += '4. res1_gridV, full grid in the region.\n'
ss  += '5. m_res1_res2V, grid within the region.\n'
ss  += '6. dfit1V, fitness gradient resident 1 direction.\n'
ss  += '7. dfit2V, fitness gradient resident 2 direction.\n'

# write contents
pickle.dump( ss, f )
pickle.dump( region, f )
pickle.dump( res1_gridV, f )
pickle.dump( res2_gridV, f )
pickle.dump( m_res1_res2V, f )
pickle.dump( dfit1V, f )
pickle.dump( dfit2V, f )
f.close()
'''
