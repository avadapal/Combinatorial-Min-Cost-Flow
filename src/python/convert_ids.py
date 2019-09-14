import networkx as nx
import os

# shifts the names from 0,...,n-1 to 1,...,n
def convert_ids(G):
	mapping = {node : int(node)+1 for node in G.nodes() }
	H=nx.relabel_nodes(G,mapping)
	return H

# main
folderpath = "../../data/instances/small_generated_graphs/"
liste = os.listdir(folderpath)
liste = [graph for graph in liste if graph.endswith('.graphml')]

for graphfile in liste:
  G = nx.read_graphml(folderpath + graphfile)
  H = convert_ids(G)
  nx.write_graphml(H, folderpath + graphfile)