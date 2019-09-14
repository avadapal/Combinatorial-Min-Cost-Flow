import networkx as nx
import pprint as pp

def create_nodes(filepath, G):
	f = open(filepath, 'r')
	lines = f.readlines()
	f.close()
	while len(lines)>1:
		coords 			= lines.pop(0)
		node_label 	=  '(' + coords.split(',')[0][8:] + ',' + str(int(coords.split(',')[1])) + ')'
		potential = lines.pop(0)
		potential = int(potential.split(':')[1])
		G.add_node( node_label, potential=potential, xx=int(coords.split(',')[1]), yy = int(coords.split(',')[0][8:]) )
		while lines and not lines[0].startswith('coords'):
			cur_line = lines.pop(0)
			if cur_line.startswith('apex potential:'):
			 	apex_potential = int(cur_line.split(':')[1])
			 	G.add_node('(-1,1)',potential=apex_potential)
	return G


def create_edges(filepath, G):
	f = open(filepath, 'r')
	lines = f.readlines()
	f.close()
	while lines:
		cur_line = lines.pop(0)
		if cur_line.startswith('coords:'):
			coords = cur_line
			y 				= int(coords.split(',')[0][8:])
			x 				= int(coords.split(',')[1])
		if cur_line.startswith('down'):
			print 'DOWN'
			first_arc 	= lines.pop(0)
			arc_data		= first_arc.split(';')
			slack 			= int(arc_data[1].split(':')[1])
			cost 				= int(arc_data[2].split(':')[1])
			flow 				= int(arc_data[4].split(':')[1])
			G.add_edge('('+str(y)+','+str(x)+')','('+str(y+1)+','+str(x)+')',cost=cost,slack=slack,flow=flow)
			print "Added edge:", '('+str(y)+','+str(x)+')','('+str(y+1)+','+str(x)+')'			

			second_arc 	= lines.pop(0)
			arc_data		= second_arc.split(';')
			slack 			= int(arc_data[1].split(':')[1])
			cost 				= int(arc_data[2].split(':')[1])
			flow 				= int(arc_data[4].split(':')[1])
			G.add_edge('('+str(y)+','+str(x)+')','('+str(y+1)+','+str(x)+')',cost=cost,slack=slack,flow=flow)
			print "Added edge:", '('+str(y)+','+str(x)+')','('+str(y+1)+','+str(x)+')'
			

		if cur_line.startswith('right'):
			print 'RIGHT'
			first_arc 	= lines.pop(0)
			arc_data		= first_arc.split(';')
			slack 			= int(arc_data[1].split(':')[1])
			cost 				= int(arc_data[2].split(':')[1])
			flow 				= int(arc_data[4].split(':')[1])
			G.add_edge('('+str(y)+','+str(x)+')','('+str(y)+','+str(x+1)+')',cost=cost,slack=slack,flow=flow)
			print "Added edge:", '('+str(y)+','+str(x)+')','('+str(y)+','+str(x+1)+')'
			
			second_arc 	= lines.pop(0)
			arc_data		= second_arc.split(';')
			slack 			= int(arc_data[1].split(':')[1])
			cost 				= int(arc_data[2].split(':')[1])
			flow 				= int(arc_data[4].split(':')[1])
			G.add_edge('('+str(y)+','+str(x)+')','('+str(y)+','+str(x+1)+')',cost=cost,slack=slack,flow=flow)
			print "Added edge:", '('+str(y)+','+str(x)+')','('+str(y)+','+str(x+1)+')'

		if cur_line.startswith('flow to apex'):
			flowapex = int(cur_line.split(':')[1])
			if flowapex>=0:
				flow_to_apex = flowapex
				flow_from_apex = 0
			else:
				flow_from_apex = -flowapex
				flow_to_apex = 0
			slacks			 = lines.pop(0).split(':')[1].split(',')
			cost_apex		 = int(lines.pop(0).split(':')[1])

			G.add_edge('('+str(x)+','+str(y)+')','(-1,-1)',cost=cost_apex,slack=int(slacks[0]),flow=flow_to_apex)
			G.add_edge('(-1,-1)','('+str(x)+','+str(y)+')',cost=0,slack=int(slacks[1]),flow=flow_from_apex)
	return G


def read_output_write_file(filepath):
	G = nx.MultiGraph()

	G = create_nodes(filepath, G)
	pp.pprint(G.nodes(data=True))
	print G.number_of_nodes()

	G = create_edges(filepath, G)
	pp.pprint(G.edges(keys=True,data=True))

	fn = filepath.rsplit("/")[-1].rsplit(".")
	nx.write_graphml(G, fn[0]+".graphml")
	print "Output file is called "+ fn[0]+".graphml"

read_output_write_file('../../data/results/lena8_last_iter.txt')