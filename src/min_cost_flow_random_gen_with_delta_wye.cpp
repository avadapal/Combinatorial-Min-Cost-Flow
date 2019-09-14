#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include<fstream>
#include "graph.h"
#include "random.h"
//#include "potential_reduction.h"
#include "potential_reduction_DS.h"
#include "min_cost_flow.h"
#include<boost/rational.hpp>

using namespace std;


int main(){
#ifdef VERBOSE
  cout << "Enter the number of nodes in the graph.. " << endl;
#endif
  unsigned int n;
  cin>>n;
  
  typedef long double RationalType;
  typedef unbounded_integer<long long> IntegerType;
  //typedef mpz_t IntegerType;
 
  Graph<IntegerType, RationalType> G0(n);
  G0.initialize_DS();
  G0.create_graph();
 
  assign_costs_resistances_capacities(G0);
#ifdef VERBOSE
  cout << "No of Edges: " << G0.no_of_edges << endl;
  cout << "No of Verticies: " << n << endl;
#endif
  Graph<IntegerType, RationalType> G(n + G0.no_of_edges,3 * G0.no_of_edges);
  Graph<IntegerType, RationalType> G_delta_wye(n + 2 * G0.no_of_edges,3 * G0.no_of_edges);
  vector<RationalType> x;
  vector<RationalType> y;
  vector<RationalType> s;
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
  sanity_check(G,x,y,s);
#endif
  
  RationalType minimum_dg = RationalType(1) - RationalType(1) / RationalType(100);
  //potential_reduction_algorithm(G0, G, G_delta_wye, x, y, s, q, minimum_dg);   
  potential_reduction_algorithm_DS(G0, G, x, y, s, q, minimum_dg); 
  unsigned int edge_index_of_G0 = 0;
  for(unsigned int i = 1; i<= G.no_of_edges; i+=3)
  {
    edge_index_of_G0++;
    x[i] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1;
    x[i+1] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2;
    x[i+2] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3;
    s[i] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s1;
    s[i+1] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s2;
    s[i+2] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s3;
  }
  
  for(unsigned int i =1; i <= G0.no_of_edges; i++)
  {
    y[G0.tails[i]] = G0.original_auxiliary_transformed_corrospondence[i].y1;
    y[G0.heads[i]] = G0.original_auxiliary_transformed_corrospondence[i].y2;
    y[G0.no_of_verticies + i] = G0.original_auxiliary_transformed_corrospondence[i].y3;
  }
  
  
//  potential_reduction_algorithm(G, x, y, s, q, minimum_dg);   
#ifndef NDEBUG
  sanity_check(G,x,y,s);
#endif
  #ifdef VERBOSE
  cout << "x: " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++)
  {
    cout << i << " = " << x[i] << " , " << endl;
  }
  cout << endl;  
#endif 
  
#ifndef NDEBUG  
  print_obj_values_and_assert(G,x,y,s);
#endif
/**   
  vector<RationalType> y_tilde(G.no_of_verticies + 1, RationalType(0) );
  vector<RationalType> s_tilde(G.no_of_edges + 1, RationalType(0) );
  postprocessing(G, x, y, s, y_tilde, s_tilde);
  
#ifndef NDEBUG 
  sanity_check(G,x,y_tilde,s_tilde);
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
  
  print_obj_values_and_assert(G,x,y_tilde,s_tilde);

  
  //print_x_y( x, y_tilde ,G);

  // check_slack_null(G,y_tilde);*/
  
  return 0;
}

