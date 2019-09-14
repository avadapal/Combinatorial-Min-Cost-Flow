import subprocess
import os
import networkx as nx
import pkg_resources
import heapq
import pprint
import os
import time
import random as r

def create_small_random_feasible_instance(n,p,C,U,B):
	while True:
		try:
			G = nx.fast_gnp_random_graph(n, p, seed=None, directed=True)
			
			# remove selfloops
			for arc in G.edges():
				# print "arc", arc[0], " ", arc[1]
				if arc[0]==arc[1]:
					G.remove_edge(arc[0],arc[1])
				for arc2 in G.edges():
					if arc2[0]==arc[1] and arc2[1]==arc[0]:
						G.remove_edge(arc[0],arc[1])
			
			demands 		= {}
			capacities	= {}
			costs 			= {}
			
			supply_sum = 0
			all_zero = True
			for node in G.nodes():
				if (node<G.number_of_nodes()-1):
					rint = r.randint(-B,B)
					if rint!=0:
						all_zero = False
					demands[node]  = rint
					supply_sum 		+= rint
				else:
					demands[node]  = -supply_sum
			if all_zero:
				continue
			for edge in G.edges():
				costs[edge] 			= r.randint(0,C)
				capacities[edge] 	= r.randint(1,U);
			
			nx.set_node_attributes(G,'demands'			, demands)
			nx.set_edge_attributes(G,'capacities'		, capacities)
			nx.set_edge_attributes(G,'costs'				, costs)
			
			# pprint.pprint(G.nodes(data=True))
			# pprint.pprint(G.edges(data=True))

			if not nx.is_connected(G.to_undirected()):
				# print "not connected"
				continue
			
			nx.network_simplex(G, demand='demands', capacity='capacities', weight='costs')
		
			return G
		except nx.exception.NetworkXUnfeasible:
			string = "Do nothing!"
			# print string



# shifts the names from 0,...,n-1 to 1,...,n
def convert_ids(G):
	mapping = {node : int(node)+1 for node in G.nodes() }
	H=nx.relabel_nodes(G,mapping)
	return H




ID = 1

while True:	
	name = 'fails_34.graphml'
	graph_name = "../../data/instances/debugging_graphs/" + name
	print "new graph..."
	G = create_small_random_feasible_instance(20,0.3,10,50,5)
	G = convert_ids(G)
	print 'ID: ', ID
	# nx.write_graphml(G, name[:-8] + str(ID) + '.graphml', prettyprint=True)
	nx.write_graphml(G, graph_name, prettyprint=True)
	# folderpath = "../../data/instances/small_generated_graphs/"
	# liste = os.listdir(folderpath)
	# liste = [graph for graph in liste if graph.endswith('.graphml')]
	
	
	# for graph in liste:
	cmd = 'time echo ' + graph_name + ' | ../min_cost_flow_generic_newopt'
	ret = os.system(cmd)
	if ret != 0:
		break
	ID += 1
