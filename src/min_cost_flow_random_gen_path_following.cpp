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

using namespace std;


int main(){
#ifndef NDEBUG
  cout << "Enter the number of nodes in the graph.. " << endl;
#endif
  int n = 0;
  cin>>n;

  typedef long double RationalType;
  //typedef unbounded_integer<long long> IntegerType;
  typedef long double IntegerType;
  Graph<IntegerType, RationalType> G0(n);
  G0.create_graph();
  G0.initialize_DS();
 
 
//  N.arcdata[0].root = 4;
  
  Network<Graph<IntegerType,RationalType>, IntegerType, RationalType> N(G0);
  assign_costs_resistances_capacities( N );
#ifndef NDEBUG
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
  
 //#if VERBOSE
  G0.print_lp();
 //#endif

  RationalType t = preprocessing( N,x,y,s,q);
  cout << "t = " << t << endl;
//   cout << "Demands Later: " << endl;
//   for(unsigned int i = 1; i <= G.no_of_verticies; i++){
//   
//     cout << i << " = " << G.demands[i] << endl;
//   }
  write_the_graph_in_dimacs_format( N );

  RationalType THRESHOLD_X =  RationalType(1)/(9*G.no_of_edges);
  RationalType THRESHOLD_S =  RationalType(1)/(9*G.no_of_edges);

  #ifdef VERBOSE
    cout << "preprocessing done" << endl;
  #endif
    
  #ifndef NDEBUG
   sanity_check( N,x,y,s,THRESHOLD_X);
  #endif
  RationalType minimum_dg = RationalType(1) - RationalType(1) / RationalType(100);
  path_following_algorithm(G0, G, N, x, y, s, q, t, minimum_dg, THRESHOLD_X, THRESHOLD_S);   
  
  //sanity_check(G,x,y,s);
  

//#ifdef VERBOSE
  cout << "x: " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++)
  {
    cout << "x[" << i << "] = " << x[i] << " , " << "s[" << i << "]" << ": " << s[i] <<endl;
  }
  cout << endl;  
//#endif 

#ifdef VERBOSE
  cout << "demands" << endl;
    for(unsigned int i = 1; i <= G.no_of_verticies; i++)
  {
    cout << "x[" << i << "] = " << G.demands[i] << " , " << endl;
  }
#endif

#ifndef NDEBUG
  print_obj_values_and_assert(G,x,y,s);
#endif
  
  vector<RationalType> y_tilde(G.no_of_verticies + 1, RationalType(0) );
  vector<RationalType> s_tilde(G.no_of_edges + 1, RationalType(0) );
    
  postprocessing( N, x, y, s, y_tilde, s_tilde, THRESHOLD_X );
  
  cout << "After Post Processing: " ;
   print_obj_values_and_assert(G,x,y_tilde,s_tilde);
   cout << "y and y_tilde " << endl;
  for(unsigned int i = 1; i<= G.no_of_verticies; i++){
  
    cout << i << ": " << y[i] << " " << y_tilde[i] << endl;
  }
  cout << endl;
  
#ifdef VERBOSE
         cout << "x: " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++)
  {
    cout << "x[" << i << "] = " << x[i] << " , " << "s["<< i << "] = " << s[i] << " , " << "s_tilde[" << i << "] = " << s_tilde[i] << endl;
  }
  cout << endl;
  
  ////////////////////////////////////
  cout << "The Graph: " << endl;
  
  for(unsigned int i = 1; i <= G.no_of_verticies; i++)
  {
    cout << i << " :  ";
    for(auto a: G.incident_edges[i])
    {
      cout << a << " , ";
    }
    cout << endl;
  }
#endif
//   vector<arc> arcs_below_threshold;
//   for(unsigned int i = 1; i <= G.no_of_edges; i++)
//   {
//     if((x[i] <= RationalType(1)/(3*G.no_of_edges)))
//     {
//       arcs_below_threshold.push_back(i);
//     }
//   }
//   
//   for(auto i: arcs_below_threshold)
//    //for(unsigned int i = 1; i <= G.no_of_edges; i++)
//   {
//    // if((x[i] >= RationalType(1)/(3*G.no_of_edges))) continue;
//     cout << "threshold = " << RationalType(1)/(3*G.no_of_edges) << endl;
//     cout << "i = " << i << ", " << "x[i] = " << x[i] << endl;
//                 node source = G.tail(i);
// 	    cout << "source = " << source << " " ; 
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
// 	#ifdef VERBOSE
// 	cout << "comes here" << endl;
// 	#endif
// 	break;
//       }
//       #ifdef VERBOSE
//       cout << "v: " << v << endl;
//       #endif
//       for( auto a : G.incident_edges[v] ) {
// 	#ifdef VERBOSE
// 	cout << "a: " << a << endl;
// 	#endif
// 	#ifdef VERBOSE
// 	cout << "s[a] = " << s[abs(a)] << endl;
// 	#endif
// 	if(s[abs(a)] >= 1 && ( a > 0 || ( a < 0 && fabs(x[i] - x[-a] > 1e-3)) )) continue;
// 	//cout << "lllllll" << endl;
//        
// 	//\ if(a < 0 && (x[abs(a)] - x[i] < -1e-4)) continue; 
// 	//cout << "lllllll" << endl;
// 	node w;
// 	if( a > 0 ) {
// 	  assert( G.tails[a] == v );
// 	  w = G.heads[a];
// 	} else {
// 	  assert( G.heads[-a] == v );
// 	  w = G.tails[-a];
// 	}
// 	#ifdef VERBOSE
// 	cout << v << " < - > " << source << endl;
// 	#endif
// 	#ifdef VERBOSE
// 	cout << w << " < - > " << target << endl;
// 	#endif
// 	if(v==source && w == target) continue;
// 	#ifdef VERBOSE
// 	cout << "w = " << w << visited[w] << endl;
// 	#endif
// 	if( !visited[w] ) {
// 	  Q.push_back( w );
// 	  visited[w] = 1;
// 	  bfstree[w] = a;
// 	  
// 	  #ifdef VERBOSE
// 	  cout << "w = "  << w << endl;
// 	  cout << "bfstree[w] = " << a << " = " << bfstree[w] << endl;
// 	  #endif
// 	}
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
//     #ifdef VERBOSE
//     cout << "mm" << endl;
//     for(auto a: bfstree){
//     
//       cout << a << " , " ;
//     }
//     cout << endl;
//     #endif
//     while(true){
//       
//       arc a = bfstree[u];
//       #ifdef VERBOSE
//       cout << "u = " << u << endl;
//       cout << "a = " << a << endl;
//       #endif
//       path.push_back(a);
//       node w;
//       if(G.head(a) == u){
// 	w = G.tail(a);
//       }
//       else{
// 	w = G.head(a);
//       }
//       if(w == source){
// 	break;
//       }
//       u = w;
//     }
// //     cout << "path size = " << path.size() << endl;
//     for(auto a1: path){
//        cout << "a1 = " << a1 << endl;
// //       cout << "m" << endl;
//       if(a1 > 0){
// 	cout << "x initially + " << x[abs(a1)] << endl;
// 	x[abs(a1)] += x[i];
// 	cout << "x later + " << x[abs(a1)] << endl;
//       }
//       else{
// 	cout << "x initially - " << x[abs(a1)] << endl;
// 	x[abs(a1)] -= x[i];
// 	cout << "x later - " << x[abs(a1)] << endl;
//       } 
//     }
//     x[i] = 0;
// 
//     }
//   }
//   cout << "DONE" << endl;
//   for(unsigned int i = 1; i <= G.no_of_edges; i++){
//   
//     cout << i << " = " << x[i] << endl; 
//  
//   }
#ifndef NDEBUG 
  //sanity_check(G,x,y_tilde,s_tilde);
#endif
  #ifndef NDEBUG
 // print_obj_values_and_assert(G,x,y,s);
#endif
//   vector<RationalType> x0(G0.no_of_edges+1,0);
//   //x = find_rounding_scheme_tree_solution(G, 1);
//    reconstruct_orig_solution(G0,G,x0,x);
  /* 
    unsigned int no_of_supplies_and_demands = 0;
   for(unsigned int i = 1; i <= G0.no_of_verticies; i++)
   {
     if(G0.demands[i] != 0)
     {
       no_of_supplies_and_demands++;
     }
   }
 */
    //Graph<IntegerType, RationalType> G2(G0.no_of_verticies+2, G0.no_of_edges + no_of_supplies_and_demands);
   //ford_fulkerson(G0, G2, x0, x0_ s_tilde);
#ifdef VERBOSE
//   cout << endl << "re-constructed original solution " << endl;
//   for(unsigned int i=1;i<= G0.no_of_edges;++i){
//     cout << "x0[" << i << "]: " << x0[i] << ", " << endl;
//   }
//   cout << endl;

  
  cout << "y_tilde" << endl << endl;
  
  for(unsigned int i = 1; i <= G.no_of_verticies; i++){
  
    cout << i << " : " << y_tilde[i] << endl;
  }
  #endif
#ifndef NDEBUG
  //print_obj_values_and_assert(G,x,y_tilde,s_tilde);
#endif
  
  //print_x_y( x, y_tilde ,G);

  // check_slack_null(G,y_tilde);
  
  ///////////////////
  
//    sanity_check(G,x,y_tilde,s_tilde);

   vector<RationalType> x0(G0.no_of_edges+1,0);
//   vector<RationalType> x0_integral(G0.no_of_edges + no_of_supplies_and_demands +1, 0);
//   vector<RationalType> s_tilde0 (G0.no_of_edges+ no_of_supplies_and_demands + 1, 0);
   
//    reconstruct_orig_solution(G0,G,x0,x);
//    reconstruct_orig_solution_s(G0,G, s_tilde0, s_tilde);
   
#ifdef VERBOSE
   cout << endl << endl;
   for(unsigned int i=1;i<= G0.no_of_edges;++i){
    cout << "s_tilde0[" << i << "]: " << s_tilde0[i] << ", " << endl;
  }
#endif

// cout << "here" << endl;
// 
//    Graph<IntegerType, RationalType> G2(G0.no_of_verticies+2, G0.no_of_edges + no_of_supplies_and_demands);
//    
//    ford_fulkerson(G0, G2, x0, x0_integral, s_tilde0);
#ifdef VERBOSE 
   cout << "x0_integral " << endl;
   for(unsigned int i = 1; i <= G0.no_of_edges; i++){
   
     cout << i << ": " << x0_integral[i] << endl;
  }
#endif
 //cout << "AFTER POSTPROCESSING:"; print_obj_values_and_assert(G0, G,x0_integral,y_tilde,s_tilde);
  
//   primal_sanity_check(G0, x0_integral);
//   cout << "we did primal_sanity_check" << endl;
//   cout << endl;  
  return 0;
}

