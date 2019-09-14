#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include<fstream>
#include "graph.h"
#include "random.h"
// #include "potential_reduction_path_following.h"
#include "min_cost_flow.h"
#include<boost/rational.hpp>

using namespace std;


int main(){

  int n = 0;
  cout << "Enter desired number of nodes: " << endl;
  cin>>n;

  typedef long double RationalType;
  typedef long double IntegerType;
  Graph<IntegerType, RationalType> G(n);
  
  // creates edges, name is stupid
  G.create_graph();
 
  Network<Graph<IntegerType,RationalType>, IntegerType, RationalType> N(G);

  // does not assign resistances, name is stupid
  assign_costs_resistances_capacities(N);
  assign_demands(N);

  cout << "Number of Edges: " << G.no_of_edges << endl;
  cout << "Number of Nodes: " << n << endl;

//  unsigned int q;

//  G0.print_lp();
  write_the_graph_in_dimacs_format( N );  

  set_initial_x_s_y( N );


 while(true){
   successive_shortest_path_step( N );
 }   
  

  cout << "Primal, Dual Solutions: " << endl;
  for(unsigned int a = 1; a <= G.no_of_edges; a++){
  }
  cout << endl;
  
  cout << "Graph: " << endl;
  cout << "Vertex: Incident Edges";
  for(unsigned int v = 1; v <= G.no_of_verticies; v++){
    cout << v << " :  ";
    for(auto a: G.incident_edges[v]){
      cout << a << ", ";
    }
    cout << endl;
  }

  return 0;
}

