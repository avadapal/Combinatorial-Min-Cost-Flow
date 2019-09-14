import networkx as nx

def produce_nx_graph(filename):
  f = open(filename,'r')
  lines = f.readlines()
  f.close()
  lines = lines[1:]

  G = nx.MultiDiGraph()
  
  for s in lines:
    curarc = [int(v) for v in s.rsplit("\t")]
    if int(curarc[2])>0:
      G.add_edge(int(curarc[0]),int(curarc[1]))
    else:
      G.add_edge(int(curarc[1]),int(curarc[0]))
  return G    


# main

liste = os.listdir("../../data/instances/gte_bad/")
liste = [graph for graph in gte_bad_liste if graph.endswith('.txt')]

for graphfile in liste[1:1]:
  G = produce_nx_graph(graphfile)
  
  fn = graphfile.rsplit("/")[-1].rsplit(".")
  nx.write_graphml(G, fn[0]+".graphml")
  print "Output file is called "+ fn[0]+".graphml"
