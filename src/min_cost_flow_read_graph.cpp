#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include<boost/rational.hpp>
#include "graph.h"
#include "random.h"
#include "potential_reduction_DS.h"
//#include "potential_reduction_path_following.h"
#include "min_cost_flow.h"
#include "verbose_prints.h"
#include <gmpxx.h>
using namespace std;


int main(){
  
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

	cout << "number of edges = " << N.G.no_of_edges << endl;
	cout << "number of nodes = " << N.G.no_of_vertices << endl;

	cout << "capacities of arcs: " << endl;
	for(unsigned int a  = 1; a <= N.G.no_of_edges; a++){
	
	  cout << a << ": " << N.arcdata[a].capacity << endl;
	}
       vector<RationalType> x; 
       vector<RationalType> y;
       vector<RationalType> s;
       RationalType q;
	
	preprocessing(N, x, y, s, q);

	cout << "q = " << q << endl;
//#ifdef VERBOSE
        cout << "output after preprocessing" << endl;
	for(unsigned int i = 1; i <= G.no_of_edges; i++){
	  cout << "lower << " << i << " = " << N.arcdata[i].xlower << endl;
	  cout << "upper << " << i << " = " << N.arcdata[i].capacity - N.arcdata[i].xlower << endl;
	  cout << "infeasibility << " << i << " = " << N.arcdata[i].infeasibility << endl;
	  cout << "direction << " << i << " = " << N.arcdata[i].direction << endl;
	  cout << "==========================================================" << endl << endl;
	}
      	 
//#endif
   	//primal_feasibility_check( N );
   	
      	potential_reduction_algorithm_new(N, q);
      	
//   cout << "Enter filepath to input graph.. " << endl;
//   string filepath;
//   cin >> filepath; 
//   int n = read_n(filepath);
//   int m = read_m(filepath); 
//   //mpz_t a, b;
//   typedef long double T;
// 
//   typedef long double RationalType;
//   typedef unbounded_integer<long long> IntegerType;
//   
//   Graph<IntegerType, RationalType> G0(n,m); 
//   G0.initialize_DS();
//   
//   G0.read_graph(filepath,m);
//   G0.print_lp();
// 
//   for(unsigned int a =1; a <= G0.no_of_edges; a++)
//   {
//     if(G0.costs[a] < 0)
//     {
//       G0.demands[G0.tails[a]] -= G0.capacities[a];
//       G0.demands[G0.heads[a]] += G0.capacities[a];
//       G0.costs[a] = -G0.costs[a];
//     }
//   }
//   Graph<IntegerType, RationalType> G(n + m,3 * m);
//   Graph<IntegerType, RationalType> G_delta_wye(n + 2 * G0.no_of_edges,3 * G0.no_of_edges); 
//   cout << "G.no_of_verticies: " << G.no_of_verticies << endl;
//   cout << "G.no_of_edges " << G.no_of_edges << endl;
//   vector<RationalType> x; 
//   vector<RationalType> y;
//   vector<RationalType> s;
//   unsigned int q;
// 
//   preprocessing(G0,G, G_delta_wye, x,y,s,q);
// 
//   write_the_graph_in_dimacs_format(G0);
//   
//  
//    //G.print_lp();
//   
//   const RationalType minimum_dg = RationalType(1) - RationalType(1) / RationalType(100); 
//   potential_reduction_algorithm(G0, G, x, y, s, q, minimum_dg);  
//   //potential_reduction_algorithm(G, G_delta_wye, x, y, s, q, minimum_dg);
//   sanity_check(G,x,y,s);
//   
//   vector<RationalType> y_tilde(G.no_of_verticies+1,0);
//   vector<RationalType> s_tilde(G.no_of_edges+1,0);
//   postprocessing(G, x, y, s, y_tilde, s_tilde, RationalType(0));
//        cout << "x: " << endl;
//   for(unsigned int i = 1; i <= G.no_of_edges; i++)
//   {
//     cout << "x[" << i << "] = " << x[i] << " , " << "s["<< i << "] = " << s[i] << " , " << "s_tilde[" << i << "] = " << s_tilde[i] << endl;
//   }
//   cout << endl;
//    unsigned int no_of_supplies_and_demands = 0;
//    for(unsigned int i = 1; i <= G0.no_of_verticies; i++)
//    {
//      if(G0.demands[i] != 0)
//      {
//        no_of_supplies_and_demands++;
//      }
//    }
//   
//     sanity_check(G,x,y_tilde,s_tilde);
// 
//    vector<T> x0(G0.no_of_edges+1,0);
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
//   
//   print_obj_values_and_assert(G,x,y_tilde,s_tilde);
//   
//   primal_sanity_check(G0, x0_integral);
//   cout << "we did primal_sanity_check" << endl;
//   cout << endl;  
  
 
 
  return 0;
}

