import networkx as nx
import numpy as np
import scipy as sp
import random as r
import pprint as pp
import time as t

def create_small_graph(L1, L2):
	G = nx.MultiDiGraph()
	G.add_node('v', demand = -1.0, y=0.0, pi=0.0)
	G.add_node('w', demand =  1.0, y=0.0, pi=0.0)
	G.add_edge('v', 'w', 1, length=L1, battery=-L1, s=L1 , x=1.0, q=0.0)
	G.add_edge('v', 'w', 2, length=L2, battery=-L2, s=L2 , x=1.0, q=0.0)
	return G

def create_small_random_feasible_instance(n,m,C):
	while True:
		try:
			G0 = nx.gnm_random_graph(n, m, seed=None, directed=True)
			G = nx.MultiDiGraph()
			for node in G0.nodes():
				G.add_node(node)
			for arc in G0.edges():
				G.add_edge(arc[0],arc[1],key=1)

			costs 			= {}
						
			for edge in G.edges(keys=True):
				costs[edge] 			= float(r.randint(1,C))
			
			batteries = { edge : -costs[edge] for edge in G.edges(keys=True)}
			s 				= { edge :  costs[edge] for edge in G.edges(keys=True)}
			x 				= { edge :  				1.0 for edge in G.edges(keys=True)}
			q 				= { edge :  				0.0 for edge in G.edges(keys=True)}
			y 				= { node :  				0.0 for node in G.nodes()}
			pi 				= { node :  				0.0 for node in G.nodes()}
			demands 	= {}
			for node in G.nodes():
				if node == 0:
					demands[node] = -1.0
				else:
					if node == G.number_of_nodes()-1:
						demands[node] =  1.0
					else:
						demands[node] =  0.0


			nx.set_edge_attributes(G,'length'	, costs)
			nx.set_edge_attributes(G,'battery', batteries)
			nx.set_edge_attributes(G,'s'			, s)
			nx.set_edge_attributes(G,'x'			, x)
			nx.set_edge_attributes(G,'q'			, q)

			nx.set_node_attributes(G,'demand'	, demands)
			nx.set_node_attributes(G,'y'			, y)
			nx.set_node_attributes(G,'pi'			, pi)

			
			if not nx.is_connected(G.to_undirected()):
				print "not connected"
				continue
			
			path = nx.shortest_path(G, source=0, target = G.number_of_nodes()-1, weight='length')
	
			print "Shortest path: ", path
			return G
		except nx.exception.NetworkXNoPath:
			print "network not feasible"

def compute_electrical_flow(G):
	l_1 = G['v']['w'][1]['length']
	l_2 = G['v']['w'][2]['length']

	r_1 = G['v']['w'][1]['s']/G['v']['w'][1]['x']
	r_2 = G['v']['w'][2]['s']/G['v']['w'][2]['x']

	g_1 = G['v']['w'][1]['battery']
	g_2 = G['v']['w'][2]['battery']

	G['v']['w'][1]['q'] =  (r_2 + l_2 - l_1)/ (r_1 + r_2)
	G['v']['w'][2]['q'] =  (r_1 + l_1 - l_2)/ (r_1 + r_2)

	G.node['v']['pi'] = 0.0
	G.node['w']['pi'] = G['v']['w'][1]['q'] * r_1 - g_1
	assert abs(G.node['w']['pi'] - (G['v']['w'][2]['q'] * r_2 - g_2)) <= 0.000001


def compute_electrical_flow_general(G,node_dict,arc_dict,g,b,A,R,L,rhs,pi,q):
	# solve les
	pi, info = sp.sparse.linalg.isolve.cg(L,rhs)
	assert info==0
	# set potentials
	for node in G.nodes():
		G.node[node]['pi'] = pi[node_dict[node]] - pi[node_dict[G.nodes()[0]]]
	# set electrical flows
	q = sp.sparse.linalg.inv(R).dot(A.T.dot(pi)+g)
	for arc in G.edges(keys=True):
		v,w,ID = arc
		G[v][w][ID]['q'] = q[arc_dict[arc]]

def dual_update(G, h):
	for v in G.nodes():
		G.node[v]['y'] = (1-h) * G.node[v]['y'] + h * G.node[v]['pi']
	for a in G.edges(keys=True):
		v, w, ID = a
		G[v][w][ID]['s'] = G.node[v]['y'] - G.node[w]['y'] + G[v][w][ID]['length']

def step(G, h, d_update):
	for a in G.edges(keys=True):
		v, w, ID = a
		G[v][w][ID]['x'] = (1-h) * G[v][w][ID]['x'] + h * G[v][w][ID]['q']
		if d_update:
			dual_update(G,h)


def run_experiment(G,h,step_count,eps,d_update):
	all_edges_fine = False
	A = sp.sparse.csc_matrix(nx.incidence_matrix(G,nodelist=G.nodes(),oriented=True))
	g = np.array([G[v][w][ID]['battery'] for (v, w, ID) in G.edges(keys=True)])
	b = np.array([G.node[v]['demand'] for v in G.nodes()])
	i = 0
	node_dict = {}
	for node in G.nodes():
		node_dict[node] = i
		i += 1
	i = 0
	arc_dict = {}
	for arc in G.edges(keys=True):
		arc_dict[arc] = i
		i += 1
	pi = np.zeros(G.number_of_nodes())
	q  = np.zeros(G.number_of_edges())

	while not all_edges_fine:
		R = sp.sparse.diags([[G[a[0]][a[1]][a[2]]['s']/G[a[0]][a[1]][a[2]]['x'] for a in G.edges(keys=True)]],[0])
		L 	= A.dot(sp.sparse.linalg.inv(R).dot(A.T))
		rhs = b - A.dot(sp.sparse.linalg.inv(R).dot(g))
		
		compute_electrical_flow_general(G, node_dict, arc_dict, g, b, A, R, L, rhs,pi, q)
		step(G,h,d_update)
		all_edges_fine = True
		for a in G.edges(keys=True):
			v, w, ID = a
			if not (G[v][w][ID]['x']<=eps or G[v][w][ID]['x']>=1-eps): 
				all_edges_fine = False
	
def print_path(G, eps):
	path_edges = []
	for a in G.edges(keys=True):
		v, w, ID = a
		if G[v][w][ID]['x']>=1-eps: 
			path_edges.append(a)
	print path_edges
	

	# for a in G.edges(keys=True):
	# 	v, w, ID = a
	# 	print "arc", a, ":",
	# 	print "x:", format(G[v][w][ID]['x'], '10f'), ", ",
	# 	print "s:", format(G[v][w][ID]['s'], '10f'), ", ",
	# 	# print "c:", format(G[v][w][ID]['length'], '10f'), ", ",
	# 	# print "g:", format(G[v][w][ID]['battery'], '10f'), ", ",
	# 	print "r:", format(G[v][w][ID]['s']/G[v][w][ID]['x'], '10f'), ", ",
	# 	print "q:", format(G[v][w][ID]['q'], '10f')
	# for v in G.nodes():
	# 	print "node", v, ":",
	# 	print "y:" , format(G.node[v]['y'], '10f'), ", ",
	# 	print "pi:", format(G.node[v]['pi'], '10f')
	# print


# n 	= 20
# m 	= 40
# C 	= 10
# G = create_small_random_feasible_instance(n, m, C)

# h 	= 1.0/G.number_of_nodes()
# eps = 1.0/G.number_of_edges()

# pp.pprint(G.edges(data=True))
# start = t.time()
# print "Starting ..."
# run_experiment(G, h, 13061990, eps, True)
# print_path(G, eps)
# print "With dual update: ", t.time()-start, "seconds needed."
# start = t.time()
# print "Starting ..."
# run_experiment(G, h, 13061990, eps, False)
# print_path(G, eps)
# print "Wout dual update: ", t.time()-start, "seconds needed."



G = nx.MultiDiGraph()

G.add_node('a', demand =  0.0, y=0.0, pi=0.0)
G.add_node('b', demand =  0.0, y=0.0, pi=0.0)
G.add_node('c', demand =  0.0, y=0.0, pi=0.0)
G.add_node('d', demand =  0.0, y=0.0, pi=0.0)

G.add_edge('a', 'b', 1, length=0, battery=2.0, s=0.0 , x=1.0, q=0.0)
G.add_edge('a', 'c', 1, length=0, battery=4.0, s=0.0 , x=2.0, q=0.0)
G.add_edge('a', 'd', 1, length=0, battery=2.0, s=0.0 , x=2.0, q=0.0)
G.add_edge('b', 'c', 1, length=0, battery=2.0, s=0.0 , x=2.0, q=0.0)
G.add_edge('c', 'd', 1, length=0, battery=2.0, s=0.0 , x=4.0, q=0.0)

A = sp.sparse.csc_matrix(nx.incidence_matrix(G,nodelist=G.nodes(),oriented=True))
g = np.array([G[v][w][ID]['battery'] for (v, w, ID) in G.edges(keys=True)])
b = np.array([G.node[v]['demand'] for v in G.nodes()])
i = 0
node_dict = {}
for node in G.nodes():
	node_dict[node] = i
	i += 1
i = 0
arc_dict = {}
for arc in G.edges(keys=True):
	arc_dict[arc] = i
	i += 1
pi = np.zeros(G.number_of_nodes())
q  = np.zeros(G.number_of_edges())

R = sp.sparse.diags([[G[a[0]][a[1]][a[2]]['x'] for a in G.edges(keys=True)]],[0])
L 	= A.dot(sp.sparse.linalg.inv(R).dot(A.T))
rhs = b - A.dot(sp.sparse.linalg.inv(R).dot(g))

compute_electrical_flow_general(G, node_dict, arc_dict, g, b, A, R, L, rhs,pi, q)

for a in G.edges(keys=True):
	v,w,ID = a
	print "q[",v,w,"]:", G[v][w][ID]['q']

for v in G.nodes():
	print "pi[",v,"]:", G.node[v]['pi']

H = nx.MultiDiGraph()

H.add_node('A', demand =  0.0, y=0.0, pi=0.0)
H.add_node('B', demand =  0.0, y=0.0, pi=0.0)
H.add_node('C', demand =  0.0, y=0.0, pi=0.0)
H.add_node('D', demand =  0.0, y=0.0, pi=0.0)
H.add_node('E', demand =  0.0, y=0.0, pi=0.0)

H.add_edge('B', 'C', 1, length=0, battery=1.0, s=0.0 , x=0.5, q=0.0)
H.add_edge('C', 'D', 1, length=0, battery=1.0/3.0, s=0.0 , x=0.25, q=0.0)
H.add_edge('D', 'E', 1, length=0, battery=0.2, s=0.0 , x=1.0/6.0, q=0.0)
H.add_edge('C', 'A', 1, length=0, battery=0.5, s=0.0 , x=1.0, q=0.0)
H.add_edge('D', 'A', 1, length=0, battery=0.25, s=0.0 , x=1.0/3.0, q=0.0)

A = sp.sparse.csc_matrix(nx.incidence_matrix(H,nodelist=H.nodes(),oriented=True))
g = np.array([H[v][w][ID]['battery'] for (v, w, ID) in H.edges(keys=True)])
b = np.array([H.node[v]['demand'] for v in H.nodes()])
i = 0
node_dict = {}
for node in H.nodes():
	node_dict[node] = i
	i += 1
i = 0
arc_dict = {}
for arc in H.edges(keys=True):
	arc_dict[arc] = i
	i += 1
pi = np.zeros(H.number_of_nodes())
q  = np.zeros(H.number_of_edges())

R = sp.sparse.diags([[H[a[0]][a[1]][a[2]]['x'] for a in H.edges(keys=True)]],[0])
L 	= A.dot(sp.sparse.linalg.inv(R).dot(A.T))
rhs = b - A.dot(sp.sparse.linalg.inv(R).dot(g))

compute_electrical_flow_general(H, node_dict, arc_dict, g, b, A, R, L, rhs,pi, q)

for a in H.edges(keys=True):
	v,w,ID = a
	print "q[",v,w,"]:", H[v][w][ID]['q']

for v in H.nodes():
	print "pi[",v,"]:", H.node[v]['pi']
