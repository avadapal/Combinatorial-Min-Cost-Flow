#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include<boost/rational.hpp>
#include "graph.h"
#include "random.h"
//#include "potential_reduction.h"
#include "potential_reduction_DS.h"
#include "min_cost_flow.h"
#include <gmpxx.h>
using namespace std;


int main(){
  cout << "Enter filepath to input graph.. " << endl;
  string filepath;
  cin >> filepath; 
  int n = read_n(filepath);
  int m = read_m(filepath);
  //mpz_t a, b;
  typedef long double T;

 // typedef boost::rational<boost::multiprecision::mpz_int> T;
 //    typedef mpq_t T;
 // typedef boost::rational<long long> T;
  //typedef boost::rational<mpz_class> T;
  
  typedef long double RationalType;
  typedef unbounded_integer<long long> IntegerType; 
  //typedef mpz_t IntegerType;
  Graph<IntegerType, RationalType> G0(n,m);
  G0.initialize_DS();
  
  G0.read_graph(filepath,m);
  G0.print_lp();
  

  Graph<IntegerType, RationalType> G(n + m,3 * m);

  vector<int> original_cost(G0.no_of_edges + 1, 0);
  for(unsigned int a =1; a <= G0.no_of_edges; a++)
  {
    original_cost[a] = G0.costs[a];
    if(G0.costs[a] < 0)
    {
      G0.demands[G0.tails[a]] -= G0.capacities[a];
      G0.demands[G0.heads[a]] += G0.capacities[a];
      G0.costs[a] = -G0.costs[a];
    }
  }
  Graph<IntegerType, RationalType> G_delta_wye(n + 2 * G0.no_of_edges,3 * G0.no_of_edges); 
  cout << "G.no_of_verticies: " << G.no_of_verticies << endl;
  cout << "G.no_of_edges " << G.no_of_edges << endl;
  vector<RationalType> x; 
  vector<RationalType> y;
  vector<RationalType> s;
  unsigned int q;

  preprocessing(G0,G, G_delta_wye, x,y,s,q);

  write_the_graph_in_dimacs_format(G0);
 
  cout << "preprocessing done" << endl;
   
   //G.print_lp();
  
  const RationalType minimum_dg = RationalType(1) - RationalType(1) / RationalType(100); 
  //potential_reduction_algorithm(G, x, y, s, q, minimum_dg);  
  //potential_reduction_algorithm(G0, G, G_delta_wye, x, y, s, q, minimum_dg);
  potential_reduction_algorithm_DS(G0, G, x, y, s, q, minimum_dg);  
  
  unsigned int edge_index_of_G0 = 0;
  for(unsigned int i = 1; i<= 3*G0.no_of_edges; i+=3)
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
  
  //sanity_check(G,x,y,s);
  //print_obj_values_and_assert(G,x,y,s);

  vector<RationalType> y_tilde(G.no_of_verticies+1,0);
  vector<RationalType> s_tilde(G.no_of_edges+1,0);
  postprocessing(G, x, y, s, y_tilde, s_tilde, RationalType(0));
    
   unsigned int no_of_supplies_and_demands = 0;
   for(unsigned int i = 1; i <= G0.no_of_verticies; i++)
   {
     if(G0.demands[i] != 0)
     {
       no_of_supplies_and_demands++;
     }
   }
  
    sanity_check(G,x,y_tilde,s_tilde);

   vector<T> x0(G0.no_of_edges+1,0);
   vector<RationalType> x0_integral(G0.no_of_edges + no_of_supplies_and_demands +1, 0);
   vector<RationalType> s_tilde0 (G0.no_of_edges+ no_of_supplies_and_demands + 1, 0);
   
   reconstruct_orig_solution(G0,G,x0,x);
   reconstruct_orig_solution_s(G0,G, s_tilde0, s_tilde);
   
#ifdef VERBOSE
   cout << endl << endl;
   for(unsigned int i=1;i<= G0.no_of_edges;++i){
    cout << "s_tilde0[" << i << "]: " << s_tilde0[i] << ", " << endl;
  }
#endif

   Graph<IntegerType, RationalType> G2(G0.no_of_verticies+2, G0.no_of_edges + no_of_supplies_and_demands);
   
   ford_fulkerson(G0, G2, x0, x0_integral, s_tilde0);
  
  print_obj_values_and_assert(G,x,y_tilde,s_tilde);
  
  int actual_obj_val = 0;
  for(unsigned int i=1;i<= G0.no_of_edges;++i){
    actual_obj_val += original_cost[i]*x0_integral[i];
  }
  primal_sanity_check(G0, x0);
  
  cout << endl;  
  
 
 
  return 0;
}

