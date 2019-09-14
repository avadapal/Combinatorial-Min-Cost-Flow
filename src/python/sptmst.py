import networkx as nx
import random as r
import math as math

def compute_tau(G,T):
	tau = 0
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

def compute_shortest_path_tree(G,source):
	paths = nx.single_source_shortest_path(G.to_undirected(), source)
	weight = nx.get_edge_attributes(G,'weight')

	T = nx.DiGraph()
	for v in G.nodes():
		T.add_node(v)

	for target in paths:
		v = paths[target][0]
		for w in paths[target][1:]:
			if not (v,w) in G.edges():
				T.add_edge(w,v,weight=weight[(w,v)])
				# print 'changed direction of ', a
			else:
				T.add_edge(v,w,weight=weight[(v,w)])
				# print 'kept direction of ', a
			v = w
	return T

def read_resistances(G,filename,iteration_number):
	f = open(filename,'r')
	# lines = f.readlines()
	i = iteration_number
	while i>0:
		line = f.next()
		if line.startswith("resistances"):
			i -= 1

	line = f.next()
	
	assert(line.startswith("1"))

	resistances = {}
	# e_iterator = nx.edges_iter(G)
	no_edges = nx.number_of_edges(G)
	edge = 0
	for edge in nx.edges(G):
		edgenumber, resistance = line.split(':')
		line = f.next()
		resistances[edge] = resistance
		print edge, resistance	
	f.close()

# main

# read graph
filename = "../../data/netg/cap1.graphml"
G = nx.read_graphml(filename)
# for edge in G.edges():
# 	print edge
# 	print edge[0],edge[1]

resistances = read_resistances(G,"../cap1_log",2)

# nx.set_edge_attributes(G,"weight",resistances)

# G = nx.DiGraph()
# for s in lines:
#   curarc = [int(v) for v in s.rsplit("\t")]
#   if int(curarc[2])>0:
#     G.add_edge(int(curarc[0]),int(curarc[1]))
#   else:
#     G.add_edge(int(curarc[1]),int(curarc[0]))

# # this is for msts
# for i in range(30):
# 	G = nx.fast_gnp_random_graph(100,0.3,directed=True)
# 	if not nx.is_connected(G.to_undirected()):
# 		print "not connected"
# 		break
# 	# generate random weights
# 	weight = {}
# 	for edge in G.edges():
# 		weight[edge] = r.random()
# 	nx.set_edge_attributes(G,'weight',weight)
	
# 	# compute mst
# 	T = compute_directed_mst(G)
	
# 	# sanity check
# 	for a in T.edges():
# 		assert T[a[0]][a[1]]['weight']==G[a[0]][a[1]]['weight']
	
# 	# compute tree condition number of mst
# 	print compute_tau(G,T)
# 	if compute_tau(G,T)<0:
# 		print nx.info(G)

# node_numbers = [10,20,40,80,160,320]
# file = open('python_output.txt', 'w')
# for n in node_numbers:
# 	G = nx.fast_gnp_random_graph(n,0.1,directed=True)
# 	while not nx.is_connected(G.to_undirected()):
# 		G = nx.fast_gnp_random_graph(n,0.1,directed=True)
# 	# generate random weights
# 	weight = {}
# 	for edge in G.edges():
# 		weight[edge] = r.random()
# 	nx.set_edge_attributes(G,'weight',weight)
# 	# print "graph edges are ", G.edges()
# 	m = G.number_of_edges()
	
# 	print 'n is ' + str(n) + ' m is ' + str(m)
# 	file.write('\n')

# 	T = compute_directed_mst(G)
# 	tau = compute_tau(G,T)
# 	print 'msts tau ', tau
# 	print 'msts tau/mlog2m ' + str(tau/ m / math.log(m)/math.log(m))


# 	count = 0

# 	sum_taudurchmlog2m = 0

# 	for source in G.nodes():
# 		T = compute_shortest_path_tree(G,source)
# 		# print "tree edges are ", T.edges()
# 		assert T.number_of_nodes() ==T.number_of_edges() + 1
# 		tau = compute_tau(G,T)
# 		print 'tau ', tau
# 		# file.write('source is ' + str(source) + ' tau/mlogm is ' + str(tau/ m / math.log(m)))
# 		# file.write('\n')
# 		# print 'source is ' + str(source) + ' tau/mlog^2m is ' + str(tau/ m / math.log(m)/math.log(m))
# 		sum_taudurchmlog2m += tau/ m / math.log(m)/math.log(m)
# 		count += 1
# 		if count == 10:
# 			break
# 	print 'n ',n,'m~ ', m, 'tau/mlog2m ', sum_taudurchmlog2m/10