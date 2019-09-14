import networkx as nx
import random as r
import math as math

def compute_tau(G,T):
	tau = 0
	weight = nx.get_edge_attributes(G,'weight')
	for a in G.edges():
		# print tuple(sorted(a))
		# print T.edges()
		if not a in T.edges():
			R_a = compute_path_resistance(a,G,T)
			tau += R_a / weight[(a[0],a[1])]
	return tau


def compute_path_resistance(a,G,T):
	T_undirected = T.to_undirected()
	path = next(nx.all_simple_paths(T_undirected, source=a[1], target=a[0]))
	weight = nx.get_edge_attributes(G,'weight')

	R_a = weight[(a[0],a[1])]
	tree_arcs = T.edges()
	# print weight
	for i in range(1,len(path)):
		if (path[i-1],path[i]) in tree_arcs:
			R_a += weight[(path[i-1],path[i])]
		else:
			R_a += weight[(path[i],path[i-1])]
	return R_a

def compute_directed_mst(G):
	T = nx.minimum_spanning_tree(G.to_undirected())
	weight = nx.get_edge_attributes(G,'weight')
	T1 = nx.DiGraph()
	
	for a in T.edges():
		if not a in G.edges():
			T1.add_edge(a[1],a[0],weight=weight[(a[1],a[0])])
			# print 'changed direction of ', a
		else:
			T1.add_edge(a[0],a[1],weight=weight[(a[0],a[1])])
			# print 'kept direction of ', a
	return T1

def create_weighted_connected_graph(n,p):
	G = nx.fast_gnp_random_graph(n,p,directed=True)
	while not nx.is_connected(G.to_undirected()):
		G = nx.fast_gnp_random_graph(n,p,directed=True)
	# generate random weights	
	weight = {}
	for edge in G.edges():
		weight[edge] = r.random()
	nx.set_edge_attributes(G,'weight',weight)
	return G

###########################################################################
###########################################################################
###########################################################################
## main
###########################################################################
###########################################################################
###########################################################################

G = create_weighted_connected_graph(20,0.2)

T = compute_directed_mst(G)
tau = compute_tau(G,T)

for a in T.edges():
	T12 = T.remove_edge(a[0],a[1])
	

















