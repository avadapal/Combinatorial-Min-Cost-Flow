import subprocess
import os
import networkx as nx
import pkg_resources
import heapq
import pprint
import os
import time
import random as r

def create_small_feasible_grid_instance(x,y,C,U):
	J = nx.grid_2d_graph(x,y,False)
	G = convert_ids(J)
	H = nx.MultiDiGraph()

	for edge in G.edges():
		H.add_edge(edge[0],edge[1],0)
		H.add_edge(edge[1],edge[0],0)

	demands 		= {}

	for node in H.nodes():
		demands[node] = 0
	nx.set_node_attributes(H,'demands'			, demands)

	for edge in H.edges(keys=True):
		# print edge
		cost = r.randint(-C,C)
		capacity = r.randint(1,U)
		H[edge[0]][edge[1]][edge[2]]['costs']       = cost
		H[edge[0]][edge[1]][edge[2]]['capacities'] 	= capacity

		if cost < 0:
			H.remove_edge(edge[0],edge[1])
			H.add_edge(edge[1],edge[0],costs=-cost,capacities=capacity)
			H.node[edge[0]]['demands'] += capacity
			H.node[edge[1]]['demands'] -= capacity
	return H


def convert_ids(G):
	mapping = {}
	i = 1
	for node in G.nodes():
		mapping[node] = i
		i += 1
	H=nx.relabel_nodes(G,mapping)
	return H


k = 2
while k<=512:
	l = 2
	while l<=512:
		for count in range(5):
			G = create_small_feasible_grid_instance(250,250,k,l)
			nx.write_graphml(G, '/KM/disopt/archive00/sar/generated_grids/grid250x250_' + str(k).zfill(3) + '_' + str(l).zfill(3) + '_' + str(count) + '.graphml', prettyprint=True)
		l *= 2
	k *= 2
