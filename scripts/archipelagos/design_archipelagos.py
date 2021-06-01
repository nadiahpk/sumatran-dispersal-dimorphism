# define some archipelagos 
# all have one mainland, which is indexed "0", and all others are small islands
# plot them to check

import networkx as nx
import matplotlib.pyplot as plt


# one mainland and 4 islands
dir_results = '../../results/archipelagos/'

h = 4
color_map = ['blue'] + h*['red']

# every island and mainland connected to every other

G = nx.Graph()
edges = [ (u, v) for u in range(h+1) for v in range(u+1,h+1) ]
# [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
G.add_edges_from(edges)

# write to file
nx.draw_circular(G, node_color=color_map)
plt.savefig(dir_results + 'all_connected.pdf')
plt.close()

# islands connected in a circle around the mainland

G = nx.Graph()

#         mainland connected every isle      1 -- last     isles connected 1 above
edges = [ (0, v) for v in range(1, h+1) ] + [ (1, h) ] + [ (u, u+1) for u in range(1, h) ]
# [(0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 2), (2, 3), (3, 4)]
G.add_edges_from(edges)

# write to file
nx.draw_spring(G, node_color=color_map)
plt.savefig(dir_results + 'circle_connected.pdf')
plt.close()
