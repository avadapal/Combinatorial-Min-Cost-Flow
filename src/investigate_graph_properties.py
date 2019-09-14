#######################################################################
#######################################################################
import networkx as nx
import numpy as np
import os
import shutil
#######################################################################
#######################################################################

# parse input file and output networkx graph
def produce_nx_graph(filename):
  f = open(filename,'r')
  lines = f.readlines()
  f.close()
  lines = lines[1:]

  G = nx.DiGraph()
  
  for s in lines:
    curarc = [int(v) for v in s.rsplit("\t")]
    if int(curarc[2])>0:
      G.add_edge(int(curarc[0]),int(curarc[1]))
    else:
      G.add_edge(int(curarc[1]),int(curarc[0]))
  return G    


# This creates a node arc incidence matrix, i couldn t find that for directed graphs?
def create_node_arc_inc_matrix(G):
	nodes = G.nodes()
	edges = G.edges()
	adj_mat = np.zeros( (len(nodes) , len(edges)) ) 
	edge_index = 0
	for uv in edges:
		u = int(uv[0])-1
		v = int(uv[1])-1
		adj_mat[u][edge_index] = -1
		adj_mat[v][edge_index] = 1
		edge_index += 1
	return adj_mat	


#######################################################################
#######################################################################
## main

# read graphs
path_to_files = "../data/"
folders = [v for v in os.listdir(os.path.abspath(path_to_files)) if (v.startswith('Exp_') and not v.endswith('tar.gz'))]

if not os.path.exists(path_to_files + 'python_output/'):
	os.makedirs(os.path.abspath(path_to_files + 'python_output/'))
else:
	shutil.rmtree(os.path.abspath(path_to_files + 'python_output'))
	os.makedirs(os.path.abspath(path_to_files + 'python_output/'))

for folder in folders:
	abs_folder = path_to_files+folder+'/'
	files = []
	for v in os.listdir(os.path.abspath(abs_folder)):
		if v.endswith('.txt') and not v.endswith('info.txt'):
			files.append(v)

	for filename in files:

		## produce graph from .txt file in data/Exp...
		G = produce_nx_graph(abs_folder+filename)
		
		## write graph down in graphml format
		nx.write_graphml(G, abs_folder+filename.split('.')[0] + ".graphml")

		## create one file for each graph size and collect info there
		if not os.path.isfile(path_to_files + 'python_output/' + str(len(G)) + 'nodes' + '.txt'):
			graphsize_file = open(path_to_files + 'python_output/' + str(len(G)) + 'nodes' + '.txt', 'w')
			graphsize_file.write('seed \t number of nodes \t node_connectivity \t edge_connectivity \t algebraic connectivity \t average degree \t clique number \t average node betweeness')
			graphsize_file.write('\n')
		else:
			graphsize_file = open(path_to_files + 'python_output/' + str(len(G)) + 'nodes' + '.txt', 'a+')
		

		# seed of graph
		graphsize_file.write(filename.split('_')[2].split('.')[0])
		# number of nodes
		graphsize_file.write('\t' + filename.split('_')[1])
		# node,edge connectivity
		ncon = nx.node_connectivity(nx.Graph(G))
		econ = nx.edge_connectivity(nx.Graph(G))
		graphsize_file.write('\t' + str(ncon))
		graphsize_file.write('\t' + str(econ))

		# compute Laplacian matrix and its eigenvalues
		adj_mat = create_node_arc_inc_matrix(G)
		L = np.dot(adj_mat, np.transpose(adj_mat))
		eig = [ round(value, 10 ) for value in sorted(np.linalg.eigvalsh(L))]
		
		# algebraic connectivity
		graphsize_file.write('\t' + str(eig[1]))
		
		# average degree
		graphsize_file.write('\t' + str(2.0*len(G.edges())/len(G.nodes())))

		# max clique
		clnum = nx.graph_clique_number(G.to_undirected())
		graphsize_file.write('\t' + str(clnum))

		nbetw = nx.betweenness_centrality(G)
		betw  = np.mean([nbetw[i+1] for i in range(len(G))])
		graphsize_file.write('\t' + str(betw))

		# close file
		graphsize_file.write('\n')
		graphsize_file.close()

		## I currently dont write the eigenvalues in the common file
		# # eigenvalues
		# graphsize_file.write('eigenvalues of laplacian: ')
		# for v in eig: 
		# 	graphsize_file.write(str(v)+', ')


		# ## write single file for every graph
		# f = open(abs_folder+filename.split('.')[0]+"_info.txt",'w')

		# f.write("Eigenvalues of Laplacian: ")
		# for v in eig: 
		# 	f.write(str(v)+', ')
		# f.write("\n")

		# f.write("Algebraic connectivity: ")
		# f.write(str(eig[1]))
		# f.write("\n")

		# f.write("Node connectivity: ")
		# f.write(str(ncon))
		# f.write("\n")

		# f.write("Edge connectivity: ")
		# f.write(str(econ))
		# f.write("\n")

		# f.write("Average degree: ")
		# f.write(str(2.0*len(G.edges())/len(G.nodes())))
		# f.write("\n")

		# f.close()