#include "graph.h"
#include "random.h"
#include "min_cost_flow_sspvariant.h"
#include "min_cost_flow_sspvariant_default.h"
#include "min_cost_flow_sspvariant_apex_grid.h"

#include <iostream>
#include <queue>
#include <stack>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <boost/rational.hpp>
#include <tuple>
#include <chrono>
#include <dirent.h>

#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/cycle_canceling.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/capacity_scaling.h>
#include <lemon/dimacs.h>

using namespace std;
using namespace lemon;


int getdir (string dir, vector<string> &files){
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

bool starts_with(const string& s1, const string& s2) {
    return s2.size() <= s1.size() && s1.compare(0, s2.size(), s2) == 0;
}

bool ends_with(const string& full_string, const string& ending) {
    if (full_string.length() >= ending.length()) {
        return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

int main(){
  cout << "Enter path to directory.. " << endl;
  string directory;
  cin >> directory;
  vector<string> files = vector<string>();

  getdir(directory,files);
  
  for(auto filename : files){        
    const string s(filename);
    if(ends_with(s, ".min")){      
      cout << filename << endl;
      filename = directory + filename;
      // check if file exists
      if (!ifstream(filename)) {
        std::cerr << "File not found: " << filename << std::endl;
          exit(1);
      }
#ifdef OUR_IMPL
      using RationalType = long long int;
      using IntegerType = long long int;

      Graph<IntegerType, RationalType> G(filename);
#ifdef RESP_CAP
      Network<decltype(G), IntegerType, RationalType, SSPVariantNodeData<IntegerType, RationalType>, SSPVariantRCArcData<IntegerType, RationalType>> N(G, filename);
#else
      Network<decltype(G), IntegerType, RationalType, SSPVariantNodeData<IntegerType, RationalType>, BasicArcData<RationalType>> N(G, filename);
#endif

      cout << "Number of Edges: " << N.G.no_of_edges << endl;
      cout << "Number of Nodes: " << N.G.no_of_vertices << endl;

      auto start = chrono::steady_clock::now();
      auto start_clock = clock();

      cout.setstate(std::ios_base::failbit);
#ifdef RESP_CAP
      successive_shortest_path_rc<decltype(N), DefaultNodeIterator, DefaultEdgeIterator, DefaultNodeAccessor, DefaultEdgeAccessor, DefaultIncidentEdgesIterator>(N, DefaultNodeIterator(1), DefaultNodeIterator(N.G.no_of_vertices + 1), DefaultEdgeIterator(1), DefaultEdgeIterator(N.G.no_of_edges + 1));
#else
      successive_shortest_path<decltype(N), DefaultNodeIterator, DefaultEdgeIterator, DefaultNodeAccessor, DefaultEdgeAccessor, DefaultIncidentEdgesIterator>(N, DefaultNodeIterator(1), DefaultNodeIterator(N.G.no_of_vertices + 1));
#endif
      cout.clear();

      // auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
      // cout << "our implementation:\t" << elapsed_time.count() << " s" << std::endl;

      const auto time_chrono = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_chrono).count();
      const auto time_clock = ((float) (clock() - start_clock)) / CLOCKS_PER_SEC;

      std::cout << "total time with chrono: " << time_chrono << "s" << std::endl;
      std::cout << "total time with clock: " << time_clock << "s" << std::endl;
#else
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
    
      // compute m and n for lemon implementations
      auto arcs = 0;
      for (SmartDigraph::ArcIt a(g); a != INVALID; ++a){
        arcs ++;
      }
      cout << "Number of Edges: " << arcs << std::endl;
      
      auto nodes = 0;
      for (SmartDigraph::NodeIt n(g); n != INVALID; ++n){
        nodes++;
      }
      cout << "Number of Nodes: " << nodes << std::endl;


#ifdef COST_SCALING
      // Initialize CostScaling algorithm object and run 
      cout << "Cost Scaling" << endl;
      auto i = 1;
      while(i<=10){
        CostScaling<SmartDigraph> cs(g);
        cs.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
        auto start = chrono::steady_clock::now();
        cs.run();
        auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
        cout << elapsed_time.count() << " s" << endl;
        i++;
      }
      cout << endl;
#endif
#ifdef NETWORK_SIMPLEX
        // Initialize NetworkSimplex algorithm object and run
      cout << "Network Simplex" << endl;
      auto i = 1;
      while(i<=10){
        NetworkSimplex<SmartDigraph> ns(g);
        ns.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
        auto start = chrono::steady_clock::now();
        ns.run();
        auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
        cout << elapsed_time.count() << " s" << endl;
        i++;
      }
      cout << endl;
#endif
#ifdef SUCC_SHORTEST_PATH
      // Initialize CapacityScaling algorithm object and run 
      cout << "Successive Shortest Path" << endl;
      auto i = 1;
      while(i<=10){
        CapacityScaling<SmartDigraph> cas(g);
        cas.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
        auto start = chrono::steady_clock::now();
        cas.run(false);
        auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
        cout << elapsed_time.count() << " s" << std::endl;
        i++;
      }
      cout << endl;
#endif
#ifdef CAP_SCALING
      cout << "Capacity Scaling" << endl;
      auto i = 1;
      while(i<=10){
        CapacityScaling<SmartDigraph> cas(g);
        cas.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
        auto start = chrono::steady_clock::now();
        cas.run();
        auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
        cout << elapsed_time.count() << " s" << std::endl;
        i++;
      }
      cout << endl;
#endif
#ifdef CYCLE_CANCELING
      // Initialize CycleCanceling algorithm object and run
      cout << "Cycle Canceling" << endl;
      auto i = 1;
      while(i<=10){
        CycleCanceling<SmartDigraph> cc(g);
        cc.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
        auto start = chrono::steady_clock::now();
        cc.run();
        auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
        cout << elapsed_time.count() << " s" << std::endl;
        i++;
      }
      cout << endl;
#endif
#endif
    }
  }


  return 0;
}

