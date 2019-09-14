#include <iostream>
#include <fstream>
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/cycle_canceling.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/capacity_scaling.h>
#include <lemon/dimacs.h>
#include <chrono>

using namespace std;
using namespace lemon;

int main() { 
	cout << "Enter filepath to input graph.. " << endl;
	string filename;
	cin >> filename;

	// check if file exists
	if (!std::ifstream(filename)) {
		std::cerr << "File not found: " << filename << std::endl;
		exit(1);
	}

	// Create a directed graph
  DIGRAPH_TYPEDEFS(SmartDigraph);
  SmartDigraph g;
  
  // Create data structures (i.e. maps) associating values to nodes and arcs of the graph
  IntArcMap lower(g), capacity(g), cost(g);
  IntNodeMap supply(g);

  // Read DIMACS input file
  ifstream input(filename);
  readDimacsMin(input, g, lower, capacity, cost, supply);
  input.close();


  // // Initialize CostScaling algorithm object and run 
  // CostScaling<SmartDigraph> cs(g);
  // cs.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
  auto start = chrono::steady_clock::now();
  // cs.run();
  auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
  // cout << "total time of cost scaling: " << elapsed_time.count() << "s" << std::endl;

  // // Initialize CapacityScaling algorithm object and run 
  // CapacityScaling<SmartDigraph> cas(g);
  // cas.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
  // start = chrono::steady_clock::now();
  // cas.run();
  // elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
  // cout << "total time of capacity scaling: " << elapsed_time.count() << "s" << std::endl;
  // start = chrono::steady_clock::now();
  // cas.run(false);
  // elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
  // cout << "total time of successive shortest path: " << elapsed_time.count() << "s" << std::endl;

  // Initialize NetworkSimplex algorithm object and run
  NetworkSimplex<SmartDigraph> ns(g);
  ns.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
  start = chrono::steady_clock::now();
  ns.run();
  elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
  cout << "total time of network simplex: " << elapsed_time.count() << "s" << std::endl;

  // // Initialize CycleCanceling algorithm object and run
  // CycleCanceling<SmartDigraph> cc(g);
  // cc.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
  // start = chrono::steady_clock::now();
  // cc.run();
  // elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
  // cout << "total time of cycle cancelling: " << elapsed_time.count() << "s" << std::endl;

  // // Print total flow cost
  // printf("Total flow cost: %d\n\n", cc.totalCost());
 	// // Print flow values on the arcs
  // printf("Flow values on arcs:\n");
  // vector<long double> x (1200000, 0);
  // for (ArcIt a(g); a != INVALID; ++a) {
  //   printf("Arc %d: %d/%d\n", g.id(a), cc.flow(a), capacity[a]);
  //    x[g.id(a)+1] = cc.flow(a);
  // }
    
  return 0;
}