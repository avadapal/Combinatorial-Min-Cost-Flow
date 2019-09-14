import networkx as nx
import os
import shutil



def read_and_write_one_file(input_file):
	G = nx.MultiDiGraph()
	f = open( input_file)
	lines = f.readlines()
	f.close()
	
	demands 		= {}
	costs   		= {}
	capacities 	= {}
	arccounter  = 0
	
	for line in lines:
		print line
		if line.startswith('p'):
			p_line = line.split(' ')
			# print p_line
			n = int(p_line[2])
			# initial_demands(int(n))
			m = int(p_line[3])
			# initial_costs(int(m))
			# initial_capacities(int(m))
		if line.startswith('n'):
			n_line = line.split(' ')
			# print n_line
			v = int(n_line[1])
			demands[v] = - int(n_line[2])
		if line.startswith('a'):
			arccounter += 1
	
			a_line = line.split(' ')
			print a_line
			v 		= int(a_line[1])
			w 		= int(a_line[2])
			cap  	= int(a_line[4])
			cost 	= int(a_line[5])
	
			# print cap,cost
			G.add_edge(v,w, capacities=cap, costs=cost)
			# edge = G.edges(v,w)
			# print edge
			# costs[(v,w)] 			= cost
			# capacities[(v,w)] 	= cap 
	
			assert(int(a_line[3])==0)
	
	for i in range(1,n+1):
		if not i in demands.viewkeys():
			demands[i]=0
	# for edge in G.edges():
	# 	print edge
	# print "Costs: ", G.edgescosts
	# print 
	# print demands
	# nx.set_edge_attributes(G,'capacities',capacities)
	# nx.set_edge_attributes(G,'costs',costs)
	nx.set_node_attributes(G,'demands',demands)
	
	print input_file + ".graphml"
	
	nx.write_graphml(G, input_file[0:-4] + ".graphml")
	print input_file[0:-4] + ".graphml"


###################
# main
f = raw_input("Enter name of file: ")
#dirlist = os.listdir(os.path.abspath(input_folder))

#for f in dirlist:
	#print f[0:-4]
read_and_write_one_file(f)