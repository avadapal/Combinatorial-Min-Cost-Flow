#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include<fstream>
#include "graph.h"
#include "random.h"
#include "potential_reduction.h"
#include "min_cost_flow.h"
#include<boost/rational.hpp>

using namespace std;


int main(){
#ifndef NDEBUG
  cout << "Enter the number of nodes in the graph.. " << endl;
#endif
  int n = 0;
  cin>>n;
  
//   mpz_class f = 1000;
//   mpz_class d = -22;
//   
//   cout << f+d << endl;
  
  //typedef long double T;
  //typedef boost::rational<unbounded_integer<long long> > T;
  //typedef boost::rational<long long> T;
  //typedef long double T;
  typedef long double RationalType;
  typedef unbounded_integer<long long> IntegerType;
  //typedef mpz_t IntegerType;
  Graph<IntegerType, RationalType> G0(n);
   G0.initialize_DS();
  G0.create_graph();
 
  assign_costs_resistances_capacities(G0);
#ifndef NDEBUG
  cout << "No of Edges: " << G0.no_of_edges << endl;
  cout << "No of Verticies: " << n << endl;
#endif
  Graph<IntegerType, RationalType> G(n + G0.no_of_edges,3 * G0.no_of_edges);
  Graph<IntegerType, RationalType> G_delta_wye(n + 2 * G0.no_of_edges,3 * G0.no_of_edges);
  vector<RationalType> x;//(G.no_of_edges+1,0);
  vector<RationalType> y;//(G.no_of_verticies+1,0);
  vector<RationalType> s;//(G.no_of_edges+1,0);
  unsigned int q;

  #if VERBOSE
  G0.print_lp();
  #endif
  
 #if VERBOSE
  G0.print_lp();
 #endif

  preprocessing(G0,G, G_delta_wye,x,y,s,q);

  write_the_graph_in_dimacs_format(G0);

#ifdef VERBOSE
  cout << "preprocessing done" << endl;
#endif
  
#ifndef NDEBUG
  //sanity_check(G,x,y,s);
#endif
  
  RationalType minimum_dg = RationalType(1) - RationalType(1) / RationalType(100);
  //potential_reduction_algorithm(G, G_delta_wye, x, y, s, q, minimum_dg);   
  potential_reduction_algorithm(G0, G, x, y, s, q, minimum_dg);   
//  sanity_check(G,x,y,s);
#ifdef VERBOSE
  cout << "x: " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++)
  {
    cout << "x[" << i << "] = " << x[i] << " , " << endl;
  }
  cout << endl;  
#endif 

#ifdef VERBOSE
  cout << "demands" << endl;
    for(unsigned int i = 1; i <= G.no_of_verticies; i++)
  {
    cout << "x[" << i << "] = " << G.demands[i] << " , " << endl;
  }
#endif

#ifndef NDEBUG
//  print_obj_values_and_assert(G,x,y,s);
#endif
  
  vector<RationalType> y_tilde(G.no_of_verticies + 1, RationalType(0) );
  vector<RationalType> s_tilde(G.no_of_edges + 1, RationalType(0) );
   cout << "y: " << endl;
  for(unsigned int i = 1; i<= G.no_of_verticies; i++){
  
    cout << i << ": " << y[i] << " " << y_tilde[i] << endl;
  }
  postprocessing(G, x, y, s, y_tilde, s_tilde, RationalType(0));
  
  cout << "y: " << endl;
  for(unsigned int i = 1; i<= G.no_of_verticies; i++){
  
    cout << i << ": " << y[i] << " " << y_tilde[i] << endl;
  }
  cout << endl;
         cout << "x: " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++)
  {
    cout << "x[" << i << "] = " << x[i] << " , " << "s["<< i << "] = " << s[i] << " , " << "s_tilde[" << i << "] = " << s_tilde[i] << endl;
  }
#ifndef NDEBUG 
  //sanity_check(G,x,y_tilde,s_tilde);
#endif
  
  vector<RationalType> x0(G0.no_of_edges+1,0);
  //x = find_rounding_scheme_tree_solution(G, 1);
   reconstruct_orig_solution(G0,G,x0,x);
   
    unsigned int no_of_supplies_and_demands = 0;
   for(unsigned int i = 1; i <= G0.no_of_verticies; i++)
   {
     if(G0.demands[i] != 0)
     {
       no_of_supplies_and_demands++;
     }
   }
 
    //Graph<IntegerType, RationalType> G2(G0.no_of_verticies+2, G0.no_of_edges + no_of_supplies_and_demands);
   //ford_fulkerson(G0, G2, x0, x0_ s_tilde);
#ifdef VERBOSE
  cout << endl << "re-constructed original solution " << endl;
  for(unsigned int i=1;i<= G0.no_of_edges;++i){
    cout << "x0[" << i << "]: " << x0[i] << ", " << endl;
  }
  cout << endl;
#endif
  
#ifndef NDEBUG
  //print_obj_values_and_assert(G,x,y_tilde,s_tilde);
#endif
  
  //print_x_y( x, y_tilde ,G);

  // check_slack_null(G,y_tilde);
  
  return 0;
}

