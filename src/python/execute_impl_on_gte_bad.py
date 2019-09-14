import subprocess
import os

gte_bad_liste = os.listdir("../../data/dimacs_bad_instances/")
gte_bad_liste = [graph for graph in gte_bad_liste if graph.endswith('.graphml')]


for graph in gte_bad_liste:
	cmd = 'time echo ../../data/dimacs_bad_instances/' + graph + '| ../min_cost_flow_read_graph_path_following > log' + graph[0:-8]
	os.system(cmd)