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
				# for arc2 in G.edges():
				# 	if arc2[0]==arc[1] and arc2[1]==arc[0]:
				# 		G.remove_edge(arc[0],arc[1])
			
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
				capacities[edge] 	= r.randint(0,U);
			
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
			# print "not feasible"

def print_step_deficit_norm(G, v, step_count, deficit):
  print 
  print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  print "Step", 
  print(format(step_count, '3d')),
  print "from ",
  print( v ),
  print ";               ||deficit||_1 =",
  def_norm = sum([abs(deficit[node]) for node in G.nodes()])
  print(format(def_norm,'12d'))
  return def_norm

def print_step_deficit_norm_resp_cap(G, (i,j,ID), step_count, deficit, arc_deficit):
	print 
	print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	print "Step", 
	print(format(step_count, '3d')),
	print "from ",
	print( i,j,ID ),
	print ";               ||deficit||_1 =",
	def_norm = sum([abs(deficit[node]) for node in G.nodes()]) + sum([abs(arc_deficit[arc]) for arc in G.edges(keys=True) ])
	print(format(def_norm,'12d'))
	return def_norm

def	print_function_values(G, x, y,
	demand='demands', capacity='capacities', weight='costs'):
	dual_value = sum([y[node]*G.node[node][demand] for node in G.nodes()])
	for arc in G.edges(keys=True):
		v, w, ID = arc
		if x[arc]>0:
			dual_value -= G[v][w][ID][capacity]*(y[w]-y[v]-G[v][w][ID][weight])
	print "dual function value:   ", dual_value
	print "primal function value: ", sum([x[arc]*G[arc[0]][arc[1]][arc[2]][weight] for arc in G.edges(keys=True)])

# returns graph that is in file located at filepath from dimacs format
def read_dimacs(filepath):
	G = nx.MultiDiGraph()

	f = open(filepath, 'r')
	lines = f.readlines()
	f.close()
	while not lines[0].startswith('p'):
		lines = lines[1:]
	
	generaldata = [int(stuff) for stuff in lines[0].split()[2:] ]
	n = generaldata[0]
	m = generaldata[1]
	
	for i in range(1,n+1):
		G.add_node(i,demands=0)

	while not lines[0].startswith('n'):
		lines = lines[1:]

	for s in lines:
		if s[0].startswith('n'):
			s = s[2:]
			nodedata = [ int(stuff) for stuff in s.split() ]
			G.add_node(nodedata[0],demands=-nodedata[1])
		if s[0].startswith('a'):
			s = s[2:]
			arcdata  = [ int(stuff) for stuff in s.split() ]
			G.add_edge(arcdata[0],arcdata[1],capacities=arcdata[3], costs=arcdata[4])
			assert(arcdata[4]>=0)
			assert(arcdata[2]==0)
	return G

# shifts the names from 0,...,n-1 to 1,...,n
def convert_ids(G):
	mapping = {node : int(node)+1 for node in G.nodes() }
	H=nx.relabel_nodes(G,mapping)
	return H

# returns graph that is in nodefilepath and arcfilepath in lasvegas format
def read_lasvegas(nodefilepath,arcfilepath):
	G 		= nx.MultiDiGraph()
	node_f 		= open(nodefilepath, 'r')
	lines = node_f.readlines()
	node_f.close()
	lines = lines[1:]
	for node_line in lines:
		nodedata = [int(string) for string in node_line.split(',')]
		G.add_node(nodedata[0],demands=nodedata[1])

	arc_f 		= open(arcfilepath, 'r')
	lines = arc_f.readlines()
	arc_f.close()
	lines = lines[1:]
	for arc_line in lines:
		arcdata = [int(string) for string in arc_line.split(',')]
		G.add_edge(arcdata[1], arcdata[0], costs=arcdata[2], capacities=200)

	sum_b = 0
	for node in G.nodes():
		sum_b += G.node[node]['demands']
	assert sum_b == 0

	return G

# for asserting primal feasibility and optimality in the end
def check_primal_feasibility_complementarity(G,x,y,
	demand='demands', capacity='capacities', weight='costs'):
	# primal feasibility
	# flow conservation
	print "Checking primal feasibility..."
	for v in G.nodes():
		sum_v = 0
		for e in G.out_edges([v],keys=True):
			sum_v -= x[e]
		for e in G.in_edges([v],keys=True):
			sum_v += x[e]
		# print "node",v,"sum_v", sum_v, "demand", G.node[v]['demand']
		assert( sum_v == G.node[v][demand] ) 
	# capacity constraints
	for a in G.edges(keys=True):
		i,j,ID = a
		# print "arc ",a,":", x[a],"<=", G[i][j]['capacity']
		assert( x[a] <=G[i][j][ID][capacity] )
	print "Checking complementary slackness..."
	for a in G.edges(keys=True):
		i,j,ID = a
		s_a = G[i][j][ID][weight] + y[i] - y[j]
		assert s_a <= 0 or x[a] == 0
		assert s_a >= 0 or x[a] == G[i][j][ID][capacity]
		assert (x[a]==0 or x[a] == G[i][j][ID][capacity]) or s_a == 0

def check_primal_feasibility_complementarity_resp_cap(G, xlower, xupper, y,virt_pot, demand='demands', capacity='capacities', weight='costs'):
	print "Checking primal feasibility..."
	for a in G.edges(keys=True):
		assert G[a[0]][a[1]][a[2]][capacity] == xlower[a] + xupper[a]
	for v in G.nodes():
		sum_v = 0
		for e in G.out_edges([v],keys=True):
			sum_v -= xlower[e]
		for e in G.in_edges([v],keys=True):
			sum_v += xlower[e]
		assert( sum_v == G.node[v][demand] ) 
	print "Checking complementary slackness..."
	for a in G.edges(keys=True):
		i,j,ID = a
		s_lower = G[i][j][ID][weight] + y[i] - virt_pot[a]
		s_upper = 									y[j] - virt_pot[a]
		# print "s_lower", s_lower, "s_upper", s_upper, "x_lower", xlower[a], "x_upper", xupper[a] 
		assert s_lower*xlower[a] == 0
		assert s_upper*xupper[a] == 0

def update_sin_sout(G, w, Sin, Sout, visited, x, y,
	demand='demands', capacity='capacities', weight='costs'):
	for w_in_edge in G.in_edges([w],keys=True):
		i,j,ID = w_in_edge
		if not visited[i]:
			if x[w_in_edge] < G[i][j][ID][capacity]:
				heapq.heappush( Sin, 
					((G[i][j][ID][weight] + y[i] - y[w], i, j, ID), (w_in_edge, True)))
			if x[w_in_edge] > 0:
				heapq.heappush( Sout, 
					((-G[i][j][ID][weight] - y[i] + y[w], i, j, ID), (w_in_edge, False)))
	for w_out_edge in G.out_edges([w],keys=True):
		j,i,ID = w_out_edge
		if not visited[i]:
			if x[w_out_edge] < G[j][i][ID][capacity]:
				heapq.heappush( Sout, 
					((G[j][i][ID][weight] + y[w] - y[i], j, i, ID), (w_out_edge, True)))
			if x[w_out_edge] > 0:
				heapq.heappush( Sin, 
					((-G[j][i][ID][weight] - y[w] + y[i], j, i, ID), (w_out_edge, False)))
	return Sin, Sout

def update_x_by_S_cut(G, x, deficit_S, a_hat, forward,
	demand='demands', capacity='capacities', weight='costs'):
	xold = x[a_hat]
	if deficit_S < 0:
		if forward: 
			x[a_hat] 	= min( G[a_hat[0]][a_hat[1]][capacity], xold + abs(deficit_S))
		else:  
			x[a_hat] 	= max( 0 , xold - abs(deficit_S))
	if deficit_S > 0:
		if forward:
			x[a_hat] 	=  min( G[a_hat[0]][a_hat[1]][capacity], xold + deficit_S )
		else:
			x[a_hat] 	= max( 0 , xold - deficit_S)
	flowincr	= x[a_hat] - xold
	return flowincr, x

def update_x_by_tree(G, T, tree_depth, x, deficit, update_what_is_there,
	demand='demands', capacity='capacities', weight='costs'):
	imbalance 	= {node : deficit[node] for node in G.nodes()}
	# imbalance2	= {node : deficit[node] for node in G.nodes()}

	while T:
		depth, arc 	= heapq.heappop(T)
		i,j,ID 				= arc
		xold 				= x[arc]
		
		if not update_what_is_there:
			if tree_depth[i] > tree_depth[j]:
				imbalance[j] += imbalance[i]
				x[arc] 				= min(max(0, xold - imbalance[i]), G[i][j][capacity])
			else:
				x[arc] 				= min(max(0, xold + imbalance[j]), G[i][j][capacity])
				imbalance[i] += imbalance[j]
		else:
			if tree_depth[i] > tree_depth[j]:
				x[arc] 				= min(max(0, xold - deficit[i]), G[i][j][ID][capacity])
			else:
				x[arc] 				= min(max(0, xold + deficit[j]), G[i][j][ID][capacity])

		flowincr = x[arc] - xold

		deficit[j] -= flowincr
		deficit[i] += flowincr
	return x, deficit

# currently only update what is there makes sth plausible
def one_step(
	G, x, y, deficit, v, update_by_S_cut, update_what_is_there,
	demand='demands', capacity='capacities', weight='costs'):

	delta       = - y[v]
	deficit_S 	= deficit[v]
	Sin  				= []
	Sout 				= []
	tree_depth	= {node : 0 for node in G.nodes()}
	visited			= {node : 0 for node in G.nodes()}
	visited[v]	= 1
	no_visited	= 1

	if not update_by_S_cut:
		T = []

	# initially fill Sin and Sout
	Sin, Sout = update_sin_sout(G, v, Sin, Sout, visited, x, y, demand, capacity, weight)
	######
	# grow the set S
	while no_visited < G.number_of_nodes() and (Sin or Sout) :	
		# print_Sin_Sout_deficitS_nu( Sin,Sout,deficit_S)

		# pop a_hat
		if deficit_S<0 or not Sin:
			(s_a_hat, node_id1, node_id2, ID), (a_hat, forward) = heapq.heappop(Sout)
			delta = s_a_hat
			if forward:
				v,w,ID = a_hat 
			else:	
				w,v,ID = a_hat
			out_it_goes = True
			# print "popping outgoing arc: ", a_hat, " with value: ", s_a_hat
		else:
			(s_a_hat, node_id1, node_id2, ID), (a_hat, forward) = heapq.heappop(Sin)
			delta = - s_a_hat
			if forward: 
				w,v,ID = a_hat 
			else: 
				v,w,ID = a_hat
			out_it_goes = False
			# print "popping ingoing arc: ", a_hat, " with value: ", s_a_hat

		# if the edge is not a cut edge (anymore)
		if visited[w]:
			continue

		# set y[w] sign is considered in delta
		y[w] = y[w] + delta 
		tree_depth[w] = tree_depth[v] + 1
		heapq.heappush(T,(-tree_depth[w], a_hat))

		# update Sin and Sout
		update_sin_sout(G, w, Sin, Sout, visited, x, y, demand, capacity, weight)

		# put w in S
		deficit_S 	+= deficit[w]
		visited[w] 	=  1
		no_visited 	+= 1

		# update deficit vector
		# if update_by_S_cut:
		# 	if out_it_goes and forward or not out_it_goes and not forward:
		# 		deficit[v] += flowincr
		# 		deficit[w] -= flowincr
		# 	else:
		# 		deficit[v] -= flowincr
		# 		deficit[w] += flowincr

		# print_x_y_deficit(G,x,y,deficit)

	# if not update_by_S_cut:
		x, deficit = update_x_by_tree(G, T, tree_depth, x, deficit, update_what_is_there, demand, capacity, weight)
	return x,y,deficit

# tie breaking is bull shit
def update_sin_sout_node(G, w, Sin, Sout, visited, arc_visited, xlower, xupper, y, virt_pot, demand='demands', capacity='capacities', weight='costs'):
	for a_hat in G.out_edges([w], keys=True):
		if not a_hat in arc_visited:
			i,j,ID = a_hat
			heapq.heappush(Sout, ((G[i][j][ID][weight] + y[i] - virt_pot[a_hat],i,j), ((i,-1,-1),(i,j,ID))))
			if xlower[a_hat] > 0:
				heapq.heappush(Sin, ((-G[i][j][ID][weight] - y[i] + virt_pot[a_hat],j,i), ((i,j,ID),(i,-1,-1))))
	for a_hat in G.in_edges([w], keys=True):
		if not a_hat in arc_visited:
			i,j,ID = a_hat
			heapq.heappush(Sout, ((y[j] - virt_pot[a_hat],i,j), ((j,-1,-1),(i,j,ID))))
			if xupper[a_hat] > 0:
				heapq.heappush(Sin, ((- y[j] + virt_pot[a_hat],j,i), ((i,j,ID),(j,-1,-1))))
	return Sin, Sout

def update_sin_sout_arc(G, a_hat, Sin, Sout, visited, arc_visited, xlower, xupper, y, virt_pot, demand='demands', capacity='capacities', weight='costs'):
	i,j,ID = a_hat
	if not visited[i]:
		heapq.heappush(Sin, ((G[i][j][ID][weight] + y[i] - virt_pot[a_hat],i,j), ((i,-1,-1),(i,j,ID))))
		if xlower[a_hat] > 0:
			heapq.heappush(Sout, ((-G[i][j][ID][weight] - y[i] + virt_pot[a_hat],i,j), ((i,j,ID),(i,-1,-1))))
	if not visited[j]:
		heapq.heappush(Sin, ((y[j] - virt_pot[a_hat],i,j), ((j,-1,-1),(i,j,ID))))
		if xupper[a_hat] > 0:
			heapq.heappush(Sout, ((- y[j] + virt_pot[a_hat],i,j), ((i,j,ID),(j,-1,-1))))
	return Sin, Sout

# todo
def update_x_by_tree_resp_cap(G, T, tree_depth, arc_tree_depth, xlower, xupper, deficit, arc_deficit, demand='demands', capacity='capacities', weight='costs'):
	An = []
	while T:
		depth,((i1,j1,ID1),(i2,j2,ID2))	= heapq.heappop(T)
		if j1==-1:
			if i1==i2:
				xold = xlower[(i2,j2,ID2)]
				if tree_depth[i1] > arc_tree_depth[(i2,j2,ID2)]:
					xlower[(i2,j2,ID2)] = max(0, xold - deficit[i1])
					flowincr = xlower[(i2,j2,ID2)] - xold
				else:
					xlower[(i2,j2,ID2)] = max(0, xold + arc_deficit[(i2,j2,ID2)])
					flowincr = xlower[(i2,j2,ID2)] - xold
			else:
				assert i1==j2
				xold = xupper[(i2,j2,ID2)]
				if tree_depth[i1] > arc_tree_depth[(i2,j2,ID2)]:
					xupper[(i2,j2,ID2)] = max(0, xold - deficit[i1])
					flowincr = xupper[(i2,j2,ID2)] - xold
				else:
					xupper[(i2,j2,ID2)] = max(0, xold + arc_deficit[(i2,j2,ID2)])
					flowincr = xupper[(i2,j2,ID2)] - xold
			deficit[i1] += flowincr
			arc_deficit[(i2,j2,ID2)] -= flowincr
			if deficit[i1] < 0:
				heapq.heappush(An, (deficit[i1], (i1,-1,-1)))
			if arc_deficit[(i2,j2,ID2)]<0:
				heapq.heappush(An, (arc_deficit[(i2,j2,ID2)], (i2,j2,ID2)))
		else:
			assert j2==-1
			if i1==i2:
				xold = xlower[(i1,j1,ID1)]
				if tree_depth[i2] > arc_tree_depth[(i1,j1,ID1)]:
					xlower[(i1,j1,ID1)] = max(0, xold - deficit[i2])
					flowincr = xlower[(i1,j1,ID1)] - xold
				else:
					xlower[(i1,j1,ID1)] = max(0, xold + arc_deficit[(i1,j1,ID1)])
					flowincr = xlower[(i1,j1,ID1)] - xold
			else:
				assert i2==j1
				xold = xupper[(i1,j1,ID1)]
				if tree_depth[i2] > arc_tree_depth[(i1,j1,ID1)]:
					xupper[(i1,j1,ID1)] = max(0, xold - deficit[i2])
					flowincr = xupper[(i1,j1,ID1)] - xold
				else:
					xupper[(i1,j1,ID1)] = max(0, xold + arc_deficit[(i1,j1,ID1)])
					flowincr = xupper[(i1,j1,ID1)] - xold
			deficit[i2] += flowincr
			arc_deficit[(i1,j1,ID1)] -= flowincr
			# knoten werden mehrfach in An gemacht. das ist quatsch!!
			if deficit[i2] < 0:
				heapq.heappush(An, (deficit[i2], (i2,-1,-1)))
			if arc_deficit[(i1,j1,ID1)]<0:
				heapq.heappush(An, (arc_deficit[(i1,j1,ID1)], (i1,j1,ID1)))
	return xlower, xupper, deficit, arc_deficit, An

def pop_min(S_n, S_a):
	if not S_n:
		return heapq.heappop(S_a)
	if not S_a:
		return heapq.heappop(S_n)
	if S_n[0]<S_a[0]:
		return heapq.heappop(S_n)
	else:
		return heapq.heappop(S_a)


def one_step_resp_cap(G, n, m, xlower, xupper, y, virt_pot, deficit, arc_deficit, i, j, ID, demand='demands', capacity='capacities', weight='costs'):
	
	no_visited		 = 1
	Sin_n  					 = []
	Sout_n 					 = []
	Sin_a  					 = []
	Sout_a 					 = []
	tree_depth		 = {node : 0 for node in G.nodes()}
	arc_tree_depth = {}
	visited				 = {node : False for node in G.nodes()}
	arc_visited 	 = set()
	T = []
	
	if j==-1:
		deficit_S 	= deficit[i]
		visited[i]	= True
		Sin_n, Sout_n 	= update_sin_sout_node(G, i, Sin_n, Sout_n, visited, arc_visited,  xlower, xupper, y, virt_pot, demand, capacity, weight)	
	else: 
		deficit_S = arc_deficit[(i,j,ID)]
		arc_visited.add((i,j,ID))
		arc_tree_depth[(i,j,ID)] = 0
		Sin_a, Sout_a = update_sin_sout_arc(G, (i,j,ID), Sin_a, Sout_a, visited, arc_visited, xlower, xupper, y, virt_pot, demand, capacity, weight)	

	iteration = 0
	
	while no_visited < n+m and (Sin_n or Sout_n or Sin_a or Sin_a):	
		iteration += 1
		# print 
		# print "Sin: ",
		# pprint.pprint(Sin)
		# print "Sout: ",
		# pprint.pprint(Sout)
		# print "visited: ",
		# pprint.pprint(visited)
		# print "arc_visited: ",
		# pprint.pprint(arc_visited)
		# if no_visited > 1 and deficit_S == 0:
		# 	break
		if deficit_S<0 or not (Sin_n or Sin_a):
			out = True
			((s_a_hat, i,j), ((i1,j1,ID1), (i2,j2,ID2))) = pop_min(Sout_n, Sout_a)
			# print "popped", ((s_a_hat, i,j), ((i1,j1,ID1), (i2,j2,ID2)))
			delta = s_a_hat
		else:
			out = False
			((s_a_hat, i,j), ((i1,j1,ID1), (i2,j2,ID2))) = pop_min(Sin_n, Sin_a)
			# print "popped", ((s_a_hat, i,j), ((i1,j1,ID1), (i2,j2,ID2)))
			delta = -s_a_hat
		if j1==-1:
			if visited[i1] and (i2,j2,ID2) in arc_visited:
				# print "Skipping arc"
				continue
			no_visited += 1
			if out:
				arc_visited.add((i2,j2,ID2))
				# print "Update from", (i2,j2) 
				virt_pot[(i2,j2,ID2)] += delta
				update_sin_sout_arc(G, (i2,j2,ID2), Sin_a, Sout_a, visited, arc_visited, xlower, xupper, y, virt_pot, demand, capacity, weight)
				deficit_S += arc_deficit[(i2,j2,ID2)]
				arc_tree_depth[(i2,j2,ID2)] = tree_depth[i1] + 1
				heapq.heappush(T,(-arc_tree_depth[(i2,j2,ID2)], ((i1,j1,ID1),(i2,j2,ID2))))
			else:
				visited[i1] =  True
				# print "Update from", i1
				y[i1] += delta
				update_sin_sout_node(G, i1, Sin_n, Sout_n, visited, arc_visited, xlower, xupper, y, virt_pot, demand, capacity, weight)
				deficit_S += deficit[i1]
				tree_depth[i1] = arc_tree_depth[(i2,j2,ID2)] + 1
				heapq.heappush(T,(-tree_depth[i1], ((i1,j1,ID1),(i2,j2,ID2))))
		else:
			assert j2==-1
			if visited[i2] and (i1,j1,ID1) in arc_visited:
				# print "Skipping arc"
				continue
			no_visited  += 1
			if out:
				visited[i2] = True
				# print "Update from", i2
				y[i2] += delta
				update_sin_sout_node(G, i2, Sin_n, Sout_n, visited, arc_visited, xlower, xupper, y, virt_pot, demand, capacity, weight)
				deficit_S +=  deficit[i2]
				tree_depth[i2] = arc_tree_depth[(i1,j1,ID1)] + 1
				heapq.heappush(T,(-tree_depth[i2], ((i1,j1,ID1),(i2,j2,ID2))))
			else:
				arc_visited.add((i1,j1,ID1))
				# print "Update from", (i1,j1)
				virt_pot[(i1,j1,ID1)] += delta
				update_sin_sout_arc(G, (i1,j1,ID1), Sin_a, Sout_a, visited, arc_visited, xlower, xupper, y, virt_pot, demand, capacity, weight)
				deficit_S += arc_deficit[(i1,j1,ID1)]
				arc_tree_depth[(i1,j1,ID1)] = tree_depth[i2] + 1
				heapq.heappush(T,(-arc_tree_depth[(i1,j1,ID1)], ((i1,j1,ID1),(i2,j2,ID2))))
			
		# print "y: ",
		# pprint.pprint(y)
		# print "virt_pot: ",
		# pprint.pprint(virt_pot)
	# print iteration, "iterations"
	xlower, xupper, deficit, arc_deficit, An = update_x_by_tree_resp_cap(G, T, tree_depth, arc_tree_depth, xlower, xupper, deficit, arc_deficit, demand, capacity, weight)
	# print "xlower: ",
	# pprint.pprint(xlower)
	# print "xupper: ",
	# pprint.pprint(xupper)
	# print

	return xlower, xupper, y, virt_pot, deficit, arc_deficit, An

def ssp_variant_resp_cap(G, n, m, demand='demands', capacity='capacities', weight='costs'):
	xlower 				= { edge : 0 for edge in G.edges(keys=True) }
	xupper 				= { arc : G[arc[0]][arc[1]][arc[2]][capacity] for arc in G.edges(keys=True) }
	y 						= { node : 0 for node in G.nodes() }
	virt_pot			= {arc : 0 for arc in G.edges(keys=True)}

	# xupper 				= { arc : 0 for arc in G.edges(keys=True) }
	# in_capacity = {node : 0 for node in G.nodes()}
	# for arc in G.edges(keys=True): 
	# 		in_capacity[arc[1]] += G[arc[0]][arc[1]][arc[2]][capacity]
	# deficit 			= { node : G.node[node][demand] - in_capacity[node] for node in G.nodes() }
	# arc_deficit 	= {arc : G[arc[0]][arc[1]][arc[2]][capacity] for arc in G.edges(keys=True)}
	deficit     = { node : G.node[node][demand] for node in G.nodes() }
	arc_deficit = {arc : 0 for arc in G.edges(keys=True)}

	An = []
	for node in G.nodes():
		if deficit[node]<0:
	 		heapq.heappush(An, (deficit[node], (node,-1,-1)))
	# for arc in G.edges():
	# 	if arc_deficit[arc]<0:
	# 		heapq.heappush(An, (arc_deficit[arc], arc))
	
	step_count = 1
	# while step_count <=6:
	while An:
		for arc in G.edges(keys=True):
			if arc_deficit[arc]:
				print arc, arc_deficit[arc]
		# print An
		deficit_v, (i, j, ID) = heapq.heappop(An)
		print_step_deficit_norm_resp_cap(G, (i, j, ID), step_count, deficit, arc_deficit)
		print_function_values(G, xlower, y, demand, capacity, weight)
		print "Start step ..."
		xlower, xupper, y, virt_pot, deficit, arc_deficit, An = one_step_resp_cap(G, n, m, xlower, xupper, y, virt_pot, deficit, arc_deficit, i, j, ID, demand, capacity, weight)
		print "... end step."
		# An = []
		# for node in G.nodes():
		# 	if deficit[node]<0:
		#  		heapq.heappush(An, (deficit[node], (node,-1,-1)))
		# for arc in G.edges(keys=True):
		# 	if arc_deficit[arc]<0:
		# 		heapq.heappush(An, (arc_deficit[arc], arc))
		step_count += 1
	
	print "----------------", 
	print(format(step_count-1, '11d')),
	print "steps needed     ----------------"
	return G, xlower, xupper, y, virt_pot, step_count-1

def ssp_variant(
	G, update_by_S_cut, update_what_is_there,
	demand='demands', capacity='capacities', weight='costs'):
	x = { edge : 0 for edge in G.edges(keys=True) }
	y = { node : 0 for node in G.nodes() }

	deficit 			= { node : G.node[node][demand] for node in G.nodes() }

	# An = []
	# for node in G.nodes():
	# 	if deficit[node]<0:
	#  		heapq.heappush(An, (deficit[node], node))

	step_count = 1

	def_norm = print_step_deficit_norm(G, 4096, step_count, deficit)
	print_function_values(G, x, y, demand, capacity, weight)
	while def_norm != 0:
		# deficit_v, v = heapq.heappop(An)
		x,y,deficit = one_step(G, x, y, deficit, 4096, update_by_S_cut, update_what_is_there, demand, capacity, weight)

		def_norm = print_step_deficit_norm(G, 4096, step_count, deficit)
		print_function_values(G, x, y, demand, capacity, weight)

		# An = []
		# for node in G.nodes():
	 # 		if deficit[node]<0:
	 # 			heapq.heappush(An, (deficit[node], node))
		
		step_count += 1
	
	print "----------------", 
	print(format(step_count-1, '11d')),
	print "steps needed     ----------------"
	return G,x,y, step_count-1


###################################################################
############################## "main" #############################
###################################################################

## for generating feasible instances of varius sizes
# for n in range(5,6,10):
# 	print "n: ", n
# 	for count in range(50,101):
# 		print "count: ", count
# while True:	
# G = create_small_random_feasible_instance(40,1,100,50,50)
# nx.write_graphml(G,'fails6.graphml', prettyprint=True)


## for reading a single graphml or dimacs file and running the algorithm on it
G = read_dimacs('../../data/instances/netgen/big2.net')
# G = nx.MultiDiGraph()
# G = nx.read_graphml('../../data/instances/small_generated_graphs/50nodes_1.graphml')
# G = convert_ids(G)
# print G.nodes(data=True)
# print G.edges(data=True)
# G = read_dimacs('../../data/images/lena512.dimacs')


# G =nx.DiGraph()
# G.add_node(1,demands=-2)
# G.add_node(2,demands=1)
# G.add_node(3,demands=1)

# G.add_edge(1,2,capacities=1,costs=0)
# G.add_edge(1,3,capacities=1,costs=100)
# G.add_edge(2,3,capacities=1,costs=0)
# pprint.pprint(G.nodes(data=True))
# pprint.pprint(G.edges(data=True))
# n = G.number_of_nodes() 
# m = G.number_of_edges()
# start = time.time()
# G, xlower, xupper, y, virt_pot, step_count = ssp_variant_resp_cap(G, n, m, demand='demands', capacity='capacities', weight='costs')
# end = time.time()
# print("--- %s seconds ---" % (end - start))
# check_primal_feasibility_complementarity_resp_cap(G,xlower,xupper,y,virt_pot)
	# break
### for reading a lasvegas file and running the algorithm on it

folderpath = "../../data/instances/lasvegas/dataPUPAx/1000/"
G=read_lasvegas(folderpath+"out_8_cut_upha_nodes.txt",folderpath+"out_8_cut_upha_arcs.txt")

n = G.number_of_nodes() 
m = G.number_of_edges()
start = time.time()
G, xlower, xupper, y, virt_pot, step_count = ssp_variant_resp_cap(G, n, m, demand='demands', capacity='capacities', weight='costs')
end = time.time()
print("--- %s seconds ---" % (end - start))
check_primal_feasibility_complementarity_resp_cap(G,xlower,xupper,y,virt_pot, demand='demands', capacity='capacities', weight='costs')

# start = time.time()
# G,x,y, step_count = ssp_variant(G, False, True, demand='demands', capacity='capacities', weight='costs')
# end = time.time()
# print("--- %s seconds ---" % (end - start))
# check_primal_feasibility_complementarity(G,x,y, 'demands', 'capacities', 'costs')


### for reading all dimacs files in a folder and running the algorithm on them
# folderpath = "../../data/instances/gte_bad/"
# print "---------------------------------------------------------------------"
# print "This executes the ssp variant implementation on all files in " + folderpath 
# print "---------------------------------------------------------------------"
# print 


# liste = os.listdir(folderpath)
# liste = [graph for graph in liste if graph.endswith('.txt')]

### if you want to write sth in a result file
# resultpath = "../../data/results/"
# resultfile = open(resultpath + 'netgen_graphs_ssp_variant.txt', 'w')
# resultfile.write(
# 			'{:<15}'.format("number of nodes: ") 	+ " "
# 		+ '{:<15}'.format("number of arcs: ")  	+ " "
# 		+ '{:<15}'.format("time monotone: ")		+ " "
# 		+ '{:<15}'.format("steps monotone: ")		+ " "
# 		+ '{:<15}'.format("time nmonotone: ")		+ " "
# 		+ '{:<15}'.format("steps nmonotone: ")	+ " "
# 		+ '\n')

# for graph_str in liste:
# 	print "---------------------------------------------------------------------"
# 	print "Graph: ", graph_str
# 	print "---------------------------------------------------------------------"
# 	G = read_dimacs(folderpath + graph_str)
	
# 	start = time.time()
# 	G,x,y,no_steps_monotone = ssp_variant(G, False, False)
# 	time_monotone = time.time() - start
# 	print("----------------      %s seconds          ----------------" % time_monotone)
#
# 	check_primal_feasibility_complementarity(G,x,y)
# 	resultfile.write(
# 			'{:>15}'.format(str(G.number_of_nodes())) + " "
# 		+ '{:>15}'.format(str(G.number_of_edges())) + " " 
# 		+ '{:>15}'.format(str(time_monotone))				+ " "
# 		+ '{:>15}'.format(str(no_steps_monotone))		+ " "
# 		+ '{:>15}'.format(str(time_nmonotone))			+ " "
# 		+ '{:>15}'.format(str(no_steps_nmontone))		+ " "
# 		+ '\n')
