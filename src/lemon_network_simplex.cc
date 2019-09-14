#include<iostream>
#include<fstream>
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/dimacs.h>
#include "graph.h"
#include "potential_reduction.h"
#include "min_cost_flow.h"
using namespace std;
using namespace lemon;

int main() {

  // Create a directed graph
  //DIGRAPH_TYPEDEFS(SmartDigraph);
  //SmartDigraph g;
//   ListGraph g;
//   // Create data structures (i.e. maps) associating values to nodes and arcs of the graph
//   ListGraph::EdgeMap<int> lower(g), capacity(g), cost(g);
//   ListGraph::NodeMap<int> supply(g);


  DIGRAPH_TYPEDEFS(SmartDigraph); 
  SmartDigraph g;
  
  // Create data structures (i.e. maps) associating values to nodes and arcs of the graph
  IntArcMap lower(g), capacity(g), cost(g);
  IntNodeMap supply(g);

  // Read DIMACS input file
  ifstream input("graph.txt"); 
  readDimacsMin(input, g, lower, capacity, cost, supply);  
  input.close();

  // Initialize NetworkSimplex algorithm object
  NetworkSimplex<SmartDigraph> ns(g);
  ns.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);

  // Run NetworkSimplex
  ns.run();
  
  // Print total flow cost
  printf("Total flow cost: %d\n\n", ns.totalCost());
  
 // Print flow values on the arcs
  vector<long double> x (1200000, 0);
   
  printf("Flow values on arcs:\n");
  for (ArcIt a(g); a != INVALID; ++a) {
    printf("Arc %d: %d/%d\n", g.id(a), ns.flow(a), capacity[a]);
    x[g.id(a)+1] = ns.flow(a);
    
  }

    long double dual_soln = 0;
    for (NodeIt a(g); a != INVALID; ++a) {
    dual_soln += -supply[a]*ns.potential(a); 
  }
   
  cout << "DUAL = " << dual_soln << endl;
  
    cout << "Enter filepath to input graph.. " << endl;
  string filepath;
  cin >> filepath; 
  int n = read_n(filepath);
  int m = read_m(filepath);
  
  
   
  typedef long double RationalType;
  typedef unbounded_integer<long long> IntegerType;
  cout << endl;
  for(int i = 1; i <=m; i++)
  {
    cout << "x[" << i << "]" << x[i] << endl;
  }
  cout << endl;
  //typedef mpz_t IntegerType;
  Graph<IntegerType, RationalType> G0(n,m); 
  G0.read_graph(filepath,m); 
   
  cout << " A look at the Graph later: " << endl;
  for(unsigned int i = 1; i <= G0.no_of_verticies; i++)
  {
    cout << i << " : " ;
    for(auto a: G0.incident_edges[i])
    {
      cout << a << " , " ;
    }
    cout << endl;
  }
  primal_sanity_check(G0, x);
  cout << "primal_sanity_check succesful" << endl;
  return 0;
}