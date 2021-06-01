import sys
sys.path.insert(0,'../../functions') # so I can import the functions
from funcs import get_M

import numpy as np

# parameters
# ---

h = 4
c = 0
f = 0.25
m = 0.3

# circular
# [(0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 2), (2, 3), (3, 4)]

# all connected
edges = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]


# try to create M based on edges
# ---

# turn edges into a matrix
M = np.zeros( (h+1, h+1) )
for u, v in edges:
    M[u,v] = 1; M[v,u] = 1


# multiply the whole thing by migration probability
M = m*M

# multiply the first column by the proportion that can migrate from the mainland, f
M[:,0] = f*M[:,0]

# divide each column through by its degree
degs = np.sum(M, axis=0) # find the degree of each island
for col in range(h+1):
    M[:,col] = M[:,col]/degs[col]

# populate the diagonals
M[0,0] = 1-m*f
for idx in range(1,h+1):
    M[idx,idx] = 1-m


# double check against our function
# ---

params = {
        'h': h,
        'c': c,
        'f': f,
        }
M_orig = get_M(m, params)

# good, they match
