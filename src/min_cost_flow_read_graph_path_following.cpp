#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include<fstream>
#include "graph.h"
#include "random.h"
#include "potential_reduction_path_following.h"
#include "min_cost_flow.h"
#include<boost/rational.hpp>
#include<boost/math/common_factor_rt.hpp>

using namespace std;


int main() {
  cout << "Enter filepath to input graph.. " << endl;
  string filename;
  cin >> filename;

  using RationalType = long double;

#ifdef WithoutRounding
  using IntegerType = long double;
#else
  using IntegerType = unbounded_integer<long long>;
#endif



  Graph<IntegerType, RationalType> G(filename);
  Network<decltype(G), IntegerType, RationalType, FullNodeData<IntegerType, RationalType>, FullArcData<IntegerType, RationalType>> N(G, filename);
  
  cout << "G.no_of_verticies: " << N.G.no_of_vertices << endl;
  cout << "G.no_of_edges " << N.G.no_of_edges << endl;

  IntegerType initial_gamma = find_gamma(N);
  IntegerType initial_beta =  find_beta(N);
  
  find_U(N, initial_beta);
  find_C(N, initial_gamma);
  
  cout << "initial_gamma = " << to_double(initial_gamma) << endl;
  cout << "initial_beta = " << to_double(initial_beta) << endl;
  
  unsigned int n      = N.G.no_of_vertices;
  unsigned int m      = 3 * N.G.no_of_edges;
  RationalType delta  = RationalType(1)/RationalType(10);
  RationalType tau    = delta / (sqrt(m));
  RationalType b = RationalType(1) / RationalType(32);
  RationalType c = RationalType(1) / RationalType(16);
  IntegerType U = find_U(N,initial_beta);
  IntegerType C = find_C(N,initial_gamma); 
  cout << "U = " << to_double(U) << endl;
  cout << "C = " << to_double(C) << endl;
  

  long long beta_scaling_factor_ = ceil( to_double(n * m * m) / (b * delta));
  IntegerType beta_scaling_factor  = beta_scaling_factor_;
  long long gamma_scaling_factor_ = ceil(RationalType(44) * m * m * m * m * to_double(U * C * (initial_beta  * beta_scaling_factor)) / ( RationalType(5) * c * delta )); 
  IntegerType gamma_scaling_factor = gamma_scaling_factor_;
  
   
  scale_b_u_c( N, beta_scaling_factor, gamma_scaling_factor, delta );
  IntegerType beta  = find_beta(N);
  IntegerType gamma = find_gamma(N);
  assert(to_double(beta) >= n * m * m / (b * delta));
  assert(to_double(gamma) >= 44 * m * m * m * m * to_double(beta * U * C) / (5 * c * delta));
  cout << "scaling done: beta = " << to_double(beta) << ", gamma = " << to_double(gamma) << endl;


  IntegerType t_path_following = compute_t(N, beta, gamma, delta);
  cout << "t_path_following = " << to_double(t_path_following) << endl;
  
  // TODO: get rid of this q in path following
  RationalType q;

  preprocessing(N, q, t_path_following);

  // Compute mu_0
  IntegerType mu_0 = (t_path_following + (beta * gamma * max(U, C)));

  cout << "mu_0 = " << to_double(mu_0) << endl;
  cout << "Upper Bound on mu_0 = " << RationalType(3) * to_double(beta * gamma *  m * U * C) << endl;
  //assert(mu_0 * delta <= RationalType(3) * beta * gamma * m * U * C);
  
 
  cout << "primal variables after preprocessing: " << endl;
  print_primal_variable(N);

  cout << "dual variables after preprocessing: " << endl;
  print_dual_variable(N);
  
  IntegerType minimum_dg = beta * gamma;
  path_following_algorithm(N, t_path_following, minimum_dg, beta, gamma, delta, mu_0, tau);

//   vector<RationalType> y_tilde(G.no_of_vertices + 1, RationalType(0) );
//   vector<RationalType> s_tilde(G.no_of_edges + 1, RationalType(0) );
//   postprocessing(G, x, y, s, y_tilde, s_tilde, THRESHOLD_X);
//    cout << "After Post Processing: " ;
//    print_obj_values_and_assert(G,x,y_tilde,s_tilde);

//    for(unsigned int i = 1; i <= G.no_of_edges; i++)
//   {
//     if(x[i] >= RationalType(1)/(3*G.no_of_edges)) continue;
//
//     cout << "i = " << i << endl;
//                 node source = G.tail(i);
//      cout << "source = " << source << " " ;
//       node target = G.head(i);
//       cout << "target = " << target << " " ;
//    // while(true)
//     {
//     vector<int> visited( G.no_of_verticies+1,0 );
//     vector<arc> bfstree( G.no_of_verticies+1, 0 );
//     deque<node> order;
//     deque<node> Q;
//     Q.push_back( source );
//     visited[source] = 1;
//
//     while( !Q.empty() ){
//       const node v = Q.front();
//       if( v == target){
//  cout << "comes here" << endl;
//  break;
//       }
//       //#ifdef VERBOSE
//       cout << "v: " << v << endl;
//       //#endif
//       for( auto a : G.incident_edges[v] ) {
//  //#ifdef VERBOSE
//  cout << "a: " << a << endl;
//  //#endif
//
//  if(abs(a) == i) continue;
//
//  if(s[abs(a)] > RationalType(1)) continue;
//
//  if(a < 0 && (x[abs(a)] < RationalType(1)/(3*G.no_of_edges))) continue;
//
//  node w;
//  if( a > 0 ) {
//    assert( G.tails[a] == v );
//    w = G.heads[a];
//  } else {
//    assert( G.heads[-a] == v );
//    w = G.tails[-a];
//  }
//  visited[target] = 0;
//  if( !visited[w] ) {
//    Q.push_back( w );
//    visited[w] = 1;
//    bfstree[w] = a;
//
//    #ifdef VERBOSE
//    cout << "w = "  << w << endl;
//    cout << "bfstree[w] = " << a << " = " << bfstree[w] << endl;
//    #endif
//  }
//       }
//       Q.pop_front();
//       order.push_front( v );
//     }
//
//     #ifdef VERBOSE
//     cout << "Q.front() = " << target << endl;
//     #endif
// //     if(Q.front() != target){
// //       break;
// //     }
//
//     node u = target;
//     vector<arc> path;
//     cout << "mm" << endl;
//     for(auto a: bfstree){
//
//       cout << a << " , " ;
//     }
//     cout << endl;
//
//     while(true){
//
//       arc a = bfstree[u];
//       cout << "u = " << u << endl;
//       cout << "a = " << a << endl;
//
//       path.push_back(a);
//       node w;
//       if(G.head(a) == u){
//  w = G.tail(a);
//       }
//       else{
//  w = G.head(a);
//       }
//       if(w == source){
//  break;
//       }
//       u = w;
//     }
//     cout << "path size = " << path.size() << endl;
//     for(auto a1: path){
//       cout << a1 << endl << endl;
//       cout << "m" << endl;
//       if(a1 > 0){
//  x[abs(a1)] -= x[i];
//       }
//       else{
//  x[abs(a1)] += x[i];
//       }
//     }
//     x[i] = 0;
//
//     }
//   }

//   cout << "DONE" << endl;
//       for(unsigned int i = 1; i <= G.no_of_edges; i++){
//
//   if(x[i] <= RationalType(1)/(3*G.no_of_edges)){
//
//     G.demands[G.tail(i)] += x[i];
//     G.demands[G.head(i)] -= x[i];
//   }
// }
// #ifndef NDEBUG
//   print_obj_values_and_assert(G,x,y,s);
// #endif
//
//      unsigned int no_of_supplies_and_demands = 0;
//    for(unsigned int i = 1; i <= G0.no_of_verticies; i++)
//    {
//      if(G0.demands[i] != 0)
//      {
//        no_of_supplies_and_demands++;
//      }
//    }
//
//     //sanity_check(G,x,y_tilde,s_tilde);
//
//    vector<RationalType> x0(G0.no_of_edges+1,0);
//    vector<RationalType> x0_integral(G0.no_of_edges + no_of_supplies_and_demands +1, 0);
//    vector<RationalType> s_tilde0 (G0.no_of_edges+ no_of_supplies_and_demands + 1, 0);
//
//    reconstruct_orig_solution(G0,G,x0,x);
//    reconstruct_orig_solution_s(G0,G, s_tilde0, s_tilde);
//
// #ifdef VERBOSE
//    cout << endl << endl;
//    for(unsigned int i=1;i<= G0.no_of_edges;++i){
//     cout << "s_tilde0[" << i << "]: " << s_tilde0[i] << ", " << endl;
//   }
// #endif
//
//    Graph<IntegerType, RationalType> G2(G0.no_of_verticies+2, G0.no_of_edges + no_of_supplies_and_demands);
//
//    ford_fulkerson(G0, G2, x0, x0_integral, s_tilde0);
//     for(unsigned int i = 1; i <= G.no_of_edges; i++){
//
//   if(x[i] <= RationalType(1)/(3*G.no_of_edges)){
//
//     G.demands[G.tail(i)] += x[i];
//     G.demands[G.head(i)] -= x[i];
//   }
// }
//   print_obj_values_and_assert(G,x,y_tilde,s_tilde);
//
//   //primal_sanity_check(G0, x0_integral);
//   cout << "we did primal_sanity_check" << endl;
//   cout << endl;


// #ifndef NDEBUG
//  // sanity_check(G,x,y_tilde,s_tilde);
// #endif
//
//   vector<RationalType> x0(G0.no_of_edges+1,0);
//   //x = find_rounding_scheme_tree_solution(G, 1);
//    reconstruct_orig_solution(G0,G,x0,x);
//
//     unsigned int no_of_supplies_and_demands = 0;
//    for(unsigned int i = 1; i <= G0.no_of_verticies; i++)
//    {
//      if(G0.demands[i] != 0)
//      {
//        no_of_supplies_and_demands++;
//      }
//    }
//
//     //Graph<IntegerType, RationalType> G2(G0.no_of_verticies+2, G0.no_of_edges + no_of_supplies_and_demands);
//    //ford_fulkerson(G0, G2, x0, x0_ s_tilde);
// #ifdef VERBOSE
//   cout << endl << "re-constructed original solution " << endl;
//   for(unsigned int i=1;i<= G0.no_of_edges;++i){
//     cout << "x0[" << i << "]: " << x0[i] << ", " << endl;
//   }
//   cout << endl;
// #endif
//
// #ifndef NDEBUG
//   print_obj_values_and_assert(G,x,y_tilde,s_tilde);
// #endif
//
//   //print_x_y( x, y_tilde ,G);
//
//   // check_slack_null(G,y_tilde);
//
  return 0;
}

