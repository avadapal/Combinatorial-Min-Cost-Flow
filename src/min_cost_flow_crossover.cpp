#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include<fstream>
#include "graph.h"
#include "simple_solver.h"
#include "random.h"
#include "potential_reduction.h"
#include "min_cost_flow.h"
#include "min_cost_crossover.h"
#include<boost/rational.hpp>

using namespace std;


int main(){
  Random rg;
  cout << "Enter the number of nodes in the graph.. " << endl;
  int n = 0;
  cin>>n;
  typedef long double RationalType;
  typedef long long IntegerType;
  //typedef boost::rational<long long> T;
  Graph<IntegerType,RationalType> G0(n);
  G0.create_graph();
  assign_costs_resistances_capacities(G0);

  Graph<IntegerType,RationalType> G(n + G0.no_of_edges,3 * G0.no_of_edges);

  vector<RationalType> x;//(G.no_of_edges+1,0);
  vector<RationalType> y;//(G.no_of_verticies+1,0);
  vector<RationalType> s;//(G.no_of_edges+1,0);
  unsigned int q;
//   G0.dijkstra(1);
//   G0.get_layers_in_spanning_tree(1,0);
//   G0.form_spanning_tree_alon();
  
  preprocessing(G0,G,x,y,s,q);
  sanity_check(G,x,y,s);

#if VERBOSE
  G.print_lp();
#endif

  RationalType minimum_dg = 0.99;
  
  RationalType duality_gap = calculate_duality_gap(G, x,s);
  vector<RationalType> s_hat_prime(G.no_of_edges + 1, 0);
  vector<RationalType> g_hat_prime(G.no_of_edges + 1, 0);
  vector<RationalType> g_prime(G.no_of_edges + 1, 0);
  vector<RationalType> s_prime(G.no_of_edges + 1, 0);
  vector<RationalType> z_prime(G.no_of_edges + 1, 0);
  RationalType beta = 0.5;
  unsigned int count = 0;
  while(duality_gap >=minimum_dg)
  {
    
    compute_s_hat_prime_g_hat_prime_s_prime_g_prime(G, s_hat_prime,  g_hat_prime, s_prime, g_prime, x, s, q, duality_gap);
//     for(unsigned int i = 1; i <= G.no_of_edges; i++)
//     {
//       T D = G.voltage_drop(i);
//       s_hat_prime[i] = D*x[i];
//       s_prime[i] = s[i]*x[i];
//       g_prime[i] = q*s_prime[i]/duality_gap - 1;            
//     }
//  
//  T s_prime_sum = 0.0;
//       for(unsigned int i = 1; i <= G.no_of_edges; i++)
//       {
// 	s_prime_sum += s_prime[i];
//       }
//     
//      for(unsigned int i = 1; i <= G.no_of_edges; i++)
//      {
//        g_hat_prime[i] = (q/s_prime_sum)*s_hat_prime[i] - 1;
//      }
//     
    RationalType t_num = 0.0;
    
    for(unsigned int i = 1; i <= G.no_of_edges; i++)
    {
      t_num += (g_hat_prime[i]*(g_hat_prime[i] - g_prime[i]));
    }
    
    RationalType t_den = 0.0;
    
    for(unsigned int i = 1; i <= G.no_of_edges; i++)
    {
     t_den += (g_prime[i] - g_hat_prime[i])*(g_prime[i] - g_hat_prime[i]); 
    }
    
    
    RationalType t = t_num/t_den;
    
    node v = ceil(rg.rng() * G.no_of_verticies);
   
    cout << "starting node for the rounding scheme : " << v << endl;
    for(unsigned int i = 1; i <= G.no_of_edges; i++)
    {
      const RationalType D = G.voltage_drop(i);
      s_hat_prime[i] = D*x[i];
      z_prime[i] = g_prime[i] - s_hat_prime[i];
    }
    
    RationalType z_norm = 0.0;
    
    for(unsigned int i =1; i <= G.no_of_edges; i++)
    {
      z_norm += z_prime[i]*z_prime[i];
    }
    
    cout << endl << "z_norm_squared: " << z_norm << endl;
    
    RationalType old_potential = calculate_potential_function(G,x,s,q);
    

    RationalType min_lambda = 0;
    RationalType min_potential = 2*old_potential; 
    cout << "min_potential: " << min_potential << endl;
    cout << "FINDING THE minimized LAMBDA: " << endl;
    
//     for(unsigned int lamba = 1; lamba <= 1000; lamba++)
//     {
//       T lamba1 = lamba * 0.001;
//     
//     vector<T> y_bar1(G.no_of_verticies + 1, 0);
//     vector<T> s_bar1(G.no_of_edges + 1, 0);  
//     vector<T> y_tilde = rounding_scheme( G, v , x, y, s);
//     for(unsigned int i = 1; i <= G.no_of_verticies; i++)
//     {
//       y_bar1[i] = lamba1*y_tilde[i] + (1-lamba1)*y[i];
//       
//     }
//     
//     s_bar1 = compute_dual_slacks( G, y_bar1);
// 
//     T new_potential1 = calculate_potential_function(G, x, s_bar1, q);
//       
//     cout << lamba1 << " - >  potential difference: " << new_potential1 - old_potential << endl;
//     if( min_potential > new_potential1 ) {
//       min_potential = new_potential1;
//       min_lambda = lamba1;
//     }
//     }
    
    min_lambda = get_optimum_lambda(G, x, s, y, q, v);
    
#ifdef VERBOSE
    cout << "minimum potential of " << min_potential << " at " << min_lambda << endl;
#endif
    vector<RationalType> y_bar(G.no_of_verticies + 1, 0);
    vector<RationalType> x_bar(G.no_of_edges + 1, 0);
    vector<RationalType> s_bar(G.no_of_edges + 1, 0);
    vector<RationalType> y_tilde = rounding_scheme( G, v , x, y, s);
    vector<RationalType> s_tilde = compute_dual_slacks( G, y_tilde);
    
    t = min_lambda;
    
    for(unsigned int i = 1; i <= G.no_of_verticies; i++)
    {
      y_bar[i] = t*y_tilde[i] + (1-t)*y[i];
      
    }
    
    s_bar = compute_dual_slacks( G, y_bar);
    
    RationalType new_potential = calculate_potential_function(G, x, s_bar, q);
    
//     cout << "x: ";
//     for (unsigned int i = 1; i <= G.no_of_edges; i++) 
//     {
//       cout << x[i] << " , " ;
//     }
//     
//     cout << endl;
//     cout << "s_bar: ";
//     for (unsigned int i = 1; i <= G.no_of_edges; i++) 
//     {
//       cout << s_bar[i] << " , " ;
//     }
//     cout << endl;
    cout << "new_potential: " << new_potential << endl;
    cout << "old_potential: " << old_potential << endl; 
    if(new_potential < old_potential)
    {
      cout << "dual step " << endl;
      s = s_bar;
      y = y_bar;
      old_potential = new_potential;
    }

    
    
 cout << "FINDING the minimized BETA: " << endl; 
  min_potential = 2 * old_potential;
//   T min_beta = 0.0;
//   for(unsigned int beta = 1; beta <= 1000; beta++)
//   {
//     vector<T> x_bar1(G.no_of_edges + 1, 0);
//     T beta1 = beta * 0.001;
//     rounding_scheme( G, v , x, y, s);
//     vector<T> x_tilde = find_rounding_scheme_tree_solution(G, v);
//     for(unsigned int i = 1; i <= G.no_of_edges; i++)
//     {
//       x_bar1[i] = beta1*x_tilde[i] + (1-beta1)*x[i];
//     }
//     T new_potential1 = calculate_potential_function(G, x_bar1, s, q);
//      cout << beta1 << " - >  potential difference: " << new_potential1 - old_potential << endl;
//      if( min_potential > new_potential1 ) {
//       min_potential = new_potential1;
//       min_beta = beta1;
//     }
//   }
 
// T min_beta = find_beta(G, x, s, y, q, v,  min_potential, old_potential);
 RationalType min_beta = get_optimum_beta(G, x, s, y, q, v);
  
#ifdef VERBOSE
 cout << endl << "minimum potential of " << min_potential << " at beta = " << min_beta << endl;
#endif
 
 beta  = min_beta;
#ifdef VERBOSE
   cout << "find_rounding_scheme_tree_solution() " << endl;
#endif
   rounding_scheme( G, v , x, y, s);
   vector<RationalType> x_tilde = find_rounding_scheme_tree_solution(G, v);
   
#ifdef VERBOSE
   cout << "x_tilde: " ;
   
   for(auto a : x_tilde)
   {
     cout << a << " , " ;
   }
   cout << endl;
#endif
//   cout << "x_bar: ";
   for(unsigned int i = 1; i <= G.no_of_edges; i++)
    {
      x_bar[i] = beta*x_tilde[i] + (1-beta)*x[i];
//       cout << x_bar[i] << " , ";
    }

    new_potential = calculate_potential_function(G, x_bar, s, q);

//         cout << "x_bar: ";
//     for (unsigned int i = 1; i <= G.no_of_edges; i++) 
//     {
//       cout << x_bar[i] << " , " ;
//     }
//     
//     cout << endl;
//     cout << "s: ";
//     for (unsigned int i = 1; i <= G.no_of_edges; i++) 
//     {
//       cout << s[i] << " , " ;
//     }
//     cout << endl;
    cout << "new_potential: " << new_potential << endl;
    cout << "old_potential: " << old_potential << endl; 
    if(new_potential < old_potential)
    {
      cout << "primal step " << endl;
      x = x_bar;
      old_potential = new_potential;
    }

    RationalType new_duality_gap = update_duality_gap(G,x,s);
    
    if(new_duality_gap < duality_gap)
    {
      count  = 0;
    }
    else
    {
      count++;
    }
    if(count == 2*G.no_of_verticies){
    
      cout << "count: " << count << "	no_of_verticies" << G.no_of_verticies << endl; 
    }
    assert(count != 4*G.no_of_verticies);
    duality_gap = new_duality_gap;
    sanity_check(G,x,y,s);
    cout << "duality_gap: " << duality_gap << "	  , t = " << t << endl;
  }
  
  sanity_check(G,x,y,s);
  print_obj_values_and_assert(G,x,y,s);
  
  vector<RationalType> y_tilde(G.no_of_verticies + 1,0);
  vector<RationalType> s_tilde(G.no_of_edges + 1,0);
  postprocessing(G, x, y, s, y_tilde, s_tilde);
  sanity_check(G,x,y_tilde,s_tilde);
  print_obj_values_and_assert(G,x,y_tilde,s_tilde);

  //print_x_y( x, y_tilde ,G);

  // check_slack_null(G,y_tilde);
  
  return 0;
}

