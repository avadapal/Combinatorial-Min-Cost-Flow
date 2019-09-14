#include<iostream>
#include<queue>
#include<vector>
#include<algorithm>
#include "graph.h"
#include "rational_functions.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include<boost/rational.hpp>
#include "random.h"

//Random rg;

using namespace std;

template <typename T> int sign(T val) {
  return (T(0) < val) - (val < T(0));
}


/** \brief does a cycle update of the current in a naive way
 * 
 * @param G The Graph
 * @param edge_index The index of the non-tree edge selected
 * @param pushed_current The amount of current being updated 
 * 
 */
template<typename IntegerType, typename RationalType>
void 
update_current(Graph<IntegerType, RationalType> &G, 
	       unsigned int edge_index, 
	       RationalType pushed_current
)
{
  #ifdef VERBOSE
    cout << "update current called" << endl;  
  #endif
  

  // push current around cycle
  std::vector<arc>& cyc = G.cycle[edge_index];
  
  for( auto y : cyc ){
    
    if(y<0){
      G.currents[-y] -= pushed_current;    // This is for the case when the edge is in the opposite direction - Therefore, the negative sign
    }
    else{	
      G.currents[y]  += pushed_current;
    }	
  }

  #ifdef VERBOSE
    cout << endl;
    cout << "currents after updating.." << endl;
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << i << " : " << G.currents[i] << endl;
    }
  #endif
}

template<typename IntegerType, typename RationalType>
IntegerType compute_alpha(
           SpanningTree<Graph<IntegerType, RationalType>,IntegerType, RationalType>& ST, 
	   Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N,
           unsigned int edge_index
)
{   

  Graph<IntegerType, RationalType>& G = ST.G;
  
    
  const arc a = G.non_tree_edges[edge_index];
  
  
  #ifdef VERBOSE
    cout << endl << "Edge Index Sampled: " << edge_index << endl;
    cout << "the non tree edge: " << a << endl;
  #endif
  

  // Get R_a, r_a, and v_u - v_v
  const IntegerType& R_a = G.resistances_accross_cycles[edge_index];
  const IntegerType& voltage_drop = ST.query(edge_index); // v_u - v_v;
  const IntegerType& r_a = N.arcdata[a].resistance;
  
  
  // Get f_a
  IntegerType f_a = G.currents[a];
  IntegerType f_ar_a = 0;
  if(f_a != 0 && r_a != 0){
    // f_a will now be f_a r_a
    multiply_by_bit_shift(f_a, r_a);
    f_ar_a = f_a;
  }
  
  const IntegerType& Delta = voltage_drop - f_ar_a;

  RationalType alpha_unrounded = - Delta / R_a ;
  
  // comment this in for rounding (check rounding procedure)
  //  IntegerType alpha = round_alpha(alpha_unrounded, Delta, R_a);
  IntegerType alpha = alpha_unrounded;

  return alpha;
}


/** \brief does a cycle update of the curernt using the datastructure
 * 
 * @param G0 The Original Graph
 * @param ST The Spanning Tree
 * @param edge_index
 * 
 */
template<typename IntegerType, typename RationalType>
void do_update_with_values(
		       SpanningTree<Graph<IntegerType, RationalType>,IntegerType, RationalType>& ST, 
		       unsigned int edge_index,
           IntegerType alpha
)
{   
  
  Graph<IntegerType, RationalType>& G = ST.G;


  // updates the data structure 
  ST.update( alpha, edge_index);
  
  // updates the non-tree current
  const arc a = G.non_tree_edges[edge_index];
  G.currents[a] -= alpha;  
  
  
}
/** \brief updates current using the LCA
 * 
 * @param ST 
 * @param edge_index 
 * @param pushed_current 
 * 
 */
template<typename IntegerType, typename RationalType>
void 
update_current_LCA( SpanningTree<Graph<IntegerType, RationalType>,IntegerType, RationalType>& ST, 
		    unsigned int edge_index, 
		    RationalType pushed_current
)
{
  Graph<IntegerType, RationalType>& G = ST.G;
  
  #ifdef VERBOSE
  cout << "update curernt" << endl;
  #endif
  
  const arc a = G.non_tree_edges[edge_index];
  const node r = ST.LCA[abs(a)];
  const node s = G.tail(a);
  const node t = G.head(a);
  #ifdef VERBOSE
  cout << r << " " << s << " " << t << " " << a << endl;
  #endif
  if( a > 0 ){
    G.currents[a] += pushed_current;
  } else {
    G.currents[-a] -= pushed_current;
  }
  for( node v = t; v != r; v = G.head(G.arcs_to_root[v]) ){
    const arc aa = G.arcs_to_root[v];
    #ifdef VERBOSE
    cout << G.tail(aa) << " " << G.head(aa) << " " << aa << endl;
    #endif
    assert( aa != 0 );
    if( aa > 0 ){
      G.currents[aa] += pushed_current;
    } else {
      G.currents[-aa] -= pushed_current;
    }
  }
  for( node v = s; v != r; v = G.head(G.arcs_to_root[v]) ){
    const arc aa = G.arcs_to_root[v];
    #ifdef VERBOSE
    cout << G.tail(aa) << " " << G.head(aa) << " " << aa << endl;
    #endif
    assert( aa != 0 );
    if( aa > 0 ){
      G.currents[aa] -= pushed_current;
    } else {
      G.currents[-aa] += pushed_current;
    }
  }
}

/** \brief Computes the gap pertaining to the electrical-flow-problem
 * 
 * @param G The graph on the electrical-flow-problem is being run
 * 
 * @return Returns the duality gap
 * 
 */ 
template<typename IntegerType, typename RationalType>
long double
compute_gap(vector<RationalType>&x, 
	    Graph<IntegerType, RationalType> &G,
	    Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N
)
{
  #ifdef VERBOSE
  cout << "compute gap" <<endl;
  #endif
  RationalType gap = 0;
  IntegerType orthocheck = 0;
  
  
  #ifdef VERBOSE
  cout << "arc\t\tflow\t\t\tresistance\t\tbatteries\t\tb-flow\t\tvoltage\t\tgap-contrib.\t\tKVL viol." << endl;
  #endif
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    
    
    
    const IntegerType& D =   G.voltage_drop(i);
    const IntegerType& r_a = N.arcdata[i].resistances;
    const IntegerType& f_a = G.currents[i];
    const IntegerType Delta = f_a*r_a - D;
    const IntegerType DeltaSquare = Delta*Delta;
    const long double contribution = to_double( DeltaSquare )/ to_double( r_a );
    
    const IntegerType& f_i = G.f[i];
    
    #ifdef VERBOSE    
    const RationalType g_a = G.batteries[i];
    cout << "g_a:  " << g_a << endl << "DeltaSquare: " << to_double(DeltaSquare) << endl << "Delta :" << to_double(Delta)
    << "f_a: " << to_double(f_a) << endl << "r_a: " << to_double(r_a) << endl << "D: " << to_double(D) << endl;   
    cout << "f_a*r_a : " << to_double(f_a*r_a) << endl << "D : " << to_double(D) << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "Delta: " << to_double(Delta)  << endl << "Delta_square: " << to_double(Delta*Delta) << endl << "R_a: " <<  to_double(r_a) << endl << "Contribution: " << contribution << endl; 
    #endif
    
    #ifdef VERBOSE 
    const RationalType KVL_violation = to_double(D) - to_double(f_a)*to_double(r_a);
    cout << i << fixed << ":\t\t" << f_a << "\t\t" << r_a << "\t\t" << g_a << "\t\t" << g_a/r_a << "\t\t" << D << "\t\t" << contribution << "\t\t" << KVL_violation << endl;
    #endif
    
    orthocheck += (f_i - f_a)*D;
    
    #ifdef VERBOSE
    cout << "ORTHOCHECK: " << to_double(orthocheck) << endl; 
    cout << "g_a: " << g_a << " , r_a: " << to_double(r_a) << " , f_a: " << to_double(f_a) << " , D: " << to_double(D) << endl; 
    #endif
    
    
    gap += contribution;
    
    if(std::isinf(gap) || std::isinf(-gap)){
      const RationalType g_a = G.batteries[i];
      cout << "g_a:  " << g_a << endl << "DeltaSquare: " << to_double(DeltaSquare) << endl << "Delta :" << to_double(Delta) << "f_a: " << to_double(f_a) << endl << "r_a: " << 
      to_double(r_a) << endl << "D: " << to_double(D) << endl;   
      cout << "f_a*r_a : " << to_double(f_a*r_a) << endl << "D : " << to_double(D) << endl;
      const IntegerType KVL_violation = D - f_a*r_a;
      //cout << i << fixed << ":\t\t" << f_a << "\t\t" << r_a << "\t\t" << g_a << "\t\t" << g_a/r_a << "\t\t" << D << "\t\t" << contribution << "\t\t" << KVL_violation << endl;
      cout << "ORTHOCHECK: " << to_double(orthocheck) << endl; 
      //cout << "g_a: " << g_a << " , r_a: " << r_a << " , f_a: " << f_a << " , D: " << D << endl; 
    }
    
    assert(!std::isinf(gap));
    assert(!std::isinf(-gap));
    
  }
  
  
  #ifdef VERBOSE
  cout << "gap = " << to_double( gap ) << endl;
  cout << "orthocheck = " << to_double( orthocheck ) << endl;
  #endif
  
  assert( orthocheck == IntegerType(0) );
  return gap; 
}

template<typename IntegerType, typename RationalType>
void 
cache_tree_induced_voltages( const SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST ) {
  for( node v = 1; v <= ST.G.no_of_verticies; ++ v ){
    ST.G.tree_induced_voltages[v] = ST.query(v);
  }
}

/** \brief Computes the currents in the tree edges using the updated currents of the non-tree edge
 * 
 * @param G0 The Original Graph
 * @param ST Spanning Tree
 * 
 */
template<typename IntegerType, typename RationalType>
void 
compute_tree_currents_by_non_tree_edge( 
				       SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST,
					Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N,
				       vector<RationalType> & x,
					RationalType THRESHOLD_X
				      )
{
  

  
  Graph<IntegerType, RationalType>& G = ST.G;
  
  #ifdef VERBOSE
    cout << "currents before finding non-tree currents"<< endl;
    for(unsigned int i = 1; i <= G.no_of_edges; i++){
      cout << i << " : " << G.currents[i] << endl;
    }
    cout << endl;
  #endif
  

  vector<IntegerType> imbalance( G.no_of_verticies+1,0 );
  vector<IntegerType> carry( G.no_of_verticies+1,0 );

  for( auto aa : G.non_tree_edges ){
    if(x[aa] < THRESHOLD_X) {
      continue;
    }
    assert( aa > 0 );

    const IntegerType alpha = G.currents[aa] - G.f_0[aa];
    imbalance[G.head(aa)] += alpha;
    imbalance[G.tail(aa)] -= alpha;
    
    #ifdef VERBOSE
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.currents[aa]) << " " << to_double(alpha) << endl;
    #endif
  }
  
  for( arc aa : boost::adaptors::reverse( G.tree_edges ) ){
    
    const node v = G.tail(aa);
    const node w = G.head(aa);
    
    assert( ST.depth[v] > ST.depth[w] );
    
    if( aa > 0 ) {
      const IntegerType& imb = imbalance[v];
      const IntegerType outflow = imb + G.f_0[aa];

      carry[w] += outflow;
      G.currents[aa] = outflow; 
      imbalance[w] += imb;
      
      #ifdef VERBOSE
        const IntegerType alpha = outflow - carry[v];
        cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.currents[aa]) << " " << to_double(alpha) << endl;
      #endif
    } else {
      const IntegerType& imb = imbalance[v];
      const IntegerType outflow = imb - G.f_0[-aa];
      
      carry[w] += outflow;
      G.currents[-aa] = -outflow; 
      imbalance[w] += imb;
      
      #ifdef VERBOSE
        const IntegerType alpha = outflow - carry[v];
        cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << -to_double(G.currents[-aa]) << " " << to_double(alpha) << endl;
      #endif
    }
  }
  #ifdef VERBOSE
    cout << "currents after finding non-tree currents"<< endl;
    for(unsigned int i = 1; i <= G.no_of_edges; i++){
      cout << i << " : " << G.currents[i] << endl;
    }
    cout << endl;
    cout << "imbalance of root: " << to_double( imbalance[ST.root] ) << endl;
    cout << "Tree Condition Number: " << ST.calculate_tree_condition_number();
  #endif
  
  assert( RationalType(1e3) * fabs( to_double( imbalance[ST.root] ) )  < RationalType(1) );
  
  compute_tree_induced_voltages_layered(ST, N);
}


/** \brief computes and the returns the gap pertaining to the electrical=flow-problem
 * 
 * @param G0 The Original Graph
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * 
 * @return Returns the gap
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType
compute_gap_by_non_tree_currents( 
				 SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> & ST,
				  Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N,
				 vector<RationalType>& x,
         vector<RationalType>& s,
         RationalType THRESHOLD_X,
	 RationalType THRESHOLD_S
)
{  

  Graph<IntegerType, RationalType>& G = ST.G;
  
  #ifdef VERBOSE
    cout << "compute gap by non-tree currents " <<endl;
    cout << "currents in non_tree_edges: " <<  endl;
    for(auto a: G.non_tree_edges){
      cout << a<< ", " << to_double(G.currents[a]) << ";" ;
    }
    cout << endl;
  #endif
  
  
  RationalType gap = 0;
  
  unsigned int edge_index = 0;
  
  IntegerType max_alpha(0);
  unsigned int edge_index_with_max_alpha = 0;
  

  // computes the maximal alpha (max_query) and the gap
  for( auto a : ST.G.non_tree_edges ){
    
    // if arc got removed or r_a=0 (because arc got contracted) - no contribution
    if( x[a] < THRESHOLD_X || s[a] < THRESHOLD_S ) {
      continue;
    }
    
    // this is the potential difference of i
    const IntegerType& D = ST.query(edge_index); 
    
    
    // get resistance and current
    const IntegerType& r_a = N.arcdata[a].resistance;
    IntegerType f_a = G.currents[a];
    
    // compute f_a*r_a
    IntegerType f_ar_a = IntegerType(0);
    if(f_a != 0 && r_a !=0){
      multiply_by_bit_shift(f_a, r_a);
      f_ar_a = f_a;
    }

    // non tree arc + query for the path
    const IntegerType Delta = f_a - D;
    const IntegerType DeltaSquare = Delta*Delta;
    
    // contribution to gap of this arc
    assert( r_a > 0 );
    const RationalType contribution = to_double( DeltaSquare )/ to_double( r_a );


    gap += contribution;

    // This is for storing the maximal query
    IntegerType R_a = G.resistances_accross_cycles[edge_index];
    RationalType new_alpha = ( f_ar_a - D ) / R_a ;
    if( new_alpha > max_alpha){
      max_alpha = new_alpha;
      edge_index_with_max_alpha = edge_index ;
    }

    edge_index++;
  }

  // now that we computed gap, we do one update with the arc that gave the largest alpha  
  assert(edge_index_with_max_alpha >= 0);
  do_update_with_values(ST, edge_index_with_max_alpha, max_alpha);
  
  return gap; 
}


/** \brief Computes the gap pertaining to the electrical-flow-problem
 * 
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * 
 * @return Returns the gap computed
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType 
compute_gap( SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> & ST )
{
  cache_tree_induced_voltages(ST);
  Graph<IntegerType, RationalType>& G = ST.G;
  #ifdef VERBOSE
  cout << "compute gap" <<endl;
  #endif
  RationalType gap = 0;
  for( auto i : ST.G.non_tree_edges ){
    
    const IntegerType& D =   G.voltage_drop(i);
    const IntegerType& r_a = G.resistances[i];
    const IntegerType& f_a = G.currents[i];
    const IntegerType Delta = f_a*r_a - D;
    const RationalType contribution = Delta*Delta/r_a;
    #ifdef VERBOSE   
    const RationalType& g_a = G.batteries[i];
    const IntegerType KVL_violation = D - f_a*r_a;
    cout << i << fixed << ":\t\t" << f_a << "\t\t" << r_a << "\t\t" << g_a << "\t\t" << D << "\t\t" << contribution << "\t\t" << KVL_violation << endl;
    #endif
    gap += contribution;
  }
  #ifdef VERBOSE
  cout << "gap = " << gap << endl;
  #endif
  return gap; 
}


/** \brief Computes the Gap given the initial current
 * 
 * @param G The Graph
 * @param initial_current The Initial Currents
 * 
 * @return Returns the gap
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType 
compute_gap_with_initial_current(Graph<IntegerType, RationalType>& G, 
				 std::vector<IntegerType>& initial_current
)
{
  for(unsigned int i=1; i<=G.no_of_verticies; i++){
    G.currents[i] = initial_current[i];
  }
  compute_tree_induced_voltages_layered(G);
  RationalType gap = compute_gap(G);
  return gap;
}


/** \brief computes and stores the the tree-induced-voltages in a layered manner
 * 
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * 
 */
template<typename IntegerType, typename RationalType>
void 
compute_tree_induced_voltages_layered( SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> & ST, 
				       Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N)

{
  
  Graph<IntegerType, RationalType>& G = ST.G;

  #ifdef VERBOSE
    cout<<"computing tree induced voltages..." << endl;
  #endif
  
  
  assert( G.head( G.tree_edges.front() ) == ST.root );
  

  // set voltage of root to zero
  G.tree_induced_voltages[ST.root] = 0;
  
  // 
  for( auto a : G.tree_edges ) {
    assert( a != 0 );
    assert( ST.depth[G.tail(a)] > ST.depth[G.head(a)] );
    
    #ifdef VERBOSE
      cout << a << " = (" << G.tail(a) << "," << G.head(a) << ")" << endl;
    #endif

    if( a > 0 ) {
      
      IntegerType f_a = G.currents[a];
      
      if(f_a != 0 && N.arcdata[a].resistance != 0){
        multiply_by_bit_shift(f_a, N.arcdata[a].resistance);
      } else {
        f_a = 0;
      }

      // set tree induced voltages of tail (lower node)
      G.tree_induced_voltages[G.tails[a]] = G.tree_induced_voltages[G.heads[a]] + f_a; 
      
    } else {
      
      const arc rev = -a;
      
      IntegerType f_a = G.currents[rev];
      
      if(f_a != 0 && N.arcdata[rev].resistances !=0){
        multiply_by_bit_shift(f_a, N.arcdata[rev].resistances);
      } else {
        f_a = 0;
      }
      
      // set tree induced voltage of head (lower node)
      G.tree_induced_voltages[G.heads[rev]] = G.tree_induced_voltages[G.tails[rev]] - f_a;   
    }

    #ifdef VERBOSE
      cout << to_double(G.tree_induced_voltages[G.tail(a)]) << " = " << to_double(G.tree_induced_voltages[G.head(a)])
      << " + " << to_double(G.resistances[abs(a)]) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)]) << endl;
    #endif
  }
  
  #if !defined(NDEBUG) || defined(VERBOSE)
  
  for( auto a : G.tree_edges ) {
    #ifdef VERBOSE
      cout << a << " = (" << G.tail(a) << "," << G.head(a) << ")" << " " << to_double( G.tree_induced_voltages[G.tail(a)] ) << " - " <<  to_double( G.tree_induced_voltages[G.head(a)] ) << " - " << to_double( G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) << endl;

      if(fabs( to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] )  ) >= 1e-3){
        cout << to_double( G.tree_induced_voltages[G.tail(a)] ) << " - " << to_double(G.tree_induced_voltages[G.head(a)]  ) 
        << " - " << to_double(G.resistances[abs(a)] ) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)] ) << endl;
        cout << to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) << endl;
        
        cout << "Assertion: " << fabs( to_double(  (G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)]) - ( G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) ) )  << " < - > " << 0.0 << endl;
        cout << to_double(G.tree_induced_voltages[G.tail(a)]) << endl;
        cout << to_double(G.tree_induced_voltages[G.head(a)]) << endl;
        cout << to_double(G.tree_induced_voltages[G.tail(a)] - G.tree_induced_voltages[G.head(a)]) << endl;
        cout << to_double(G.resistances[abs(a)]) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)]) << endl;
        cout << to_double(G.resistances[abs(a)] * sign(a)) << " * " << to_double(G.currents[abs(a)]) << endl;
        cout << to_double(G.resistances[abs(a)] * sign(a) * G.currents[abs(a)]) << endl;
      }
    #endif
    
    assert(  fabs( to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(N.arcdata[abs(a)].resistances*sign(a)*G.currents[abs(a)] )  ) < 1e-3 );
  }
  #endif
}

/** \brief Finds the Sum of the resistances across the cycle
 * 
 * @param ST Spanning Tree
 * @param a Non-Tree arcs_to_root
 * @param f The flow
 * 
 */
template< typename IntegerType, typename RationalType >
IntegerType
sum_over_cycle( const SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST, 
		arc a, const vector<IntegerType>& f, 
		int sign = 1 
) 
{
  const Graph<IntegerType, RationalType>& G = ST.G;
  const node r = ST.LCA[abs(a)];
  const node s = G.tail(a);
  const node t = G.head(a);
  IntegerType sum = a > 0 ? f[a] : sign*f[-a];
  #ifdef VERBOSE
  cout << r << " " << s << " " << t << " " << a << " " << (a > 0 ? f[a] : sign*f[-a]) << endl;
  #endif
  for( node v = t; v != r; v = G.head(G.arcs_to_root[v]) ){
    const arc aa = G.arcs_to_root[v];
    sum += aa > 0 ? f[aa] : sign*f[-aa];
    #ifdef VERBOSE
    cout << G.tail(aa) << " " << G.head(aa) << " " << aa << " " << (aa > 0 ? f[aa] : sign*f[-aa]) << endl;
    #endif
    assert( aa != 0 );
  }
  for( node v = s; v != r; v = G.head(G.arcs_to_root[v]) ){
    const arc aa = G.arcs_to_root[v];
    sum += aa > 0 ? sign*f[aa] : f[-aa];
    #ifdef VERBOSE
    cout << G.tail(aa) << " " << G.head(aa) << " " << aa << " " << -(aa > 0 ? f[aa] : sign*f[-aa]) << endl;
    #endif
    assert( aa != 0 );
  }
  return sum;
}

/** \brief the initial state to begin the update using datastructure is built
 * 
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * 
 */
template<typename IntegerType, typename RationalType> 
void 
build_initial_state(SpanningTree<Graph<IntegerType, RationalType> ,IntegerType, 
		    RationalType> & ST,
		    Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N
)
{
  
  Graph<IntegerType, RationalType>& G = ST.G;
  
  
  #ifdef VERBOSE
    cout << "build_initial_state() " << endl;
    cout << "Initial currents ... " << endl;
    for(unsigned int i=1; i <= G.no_of_edges; i++) {
      cout << i << " : " << to_double(G.currents[i]) << endl;
    }
    cout << "init datastructure with initial currents" << endl;
    cout << "tree edges first" << endl;
  #endif

  // clears d_drop and d_ext
  ST.clear_d();


  // carry of flow
  vector<IntegerType> carry( G.no_of_verticies+1, RationalType(0) );
  

  // calls update for the nodes from the leaves to the root
  for( arc a : boost::adaptors::reverse( G.tree_edges ) ) {
    
    const node v = G.tail(a);
    const node w = G.head(a);
    assert( ST.depth[v] > ST.depth[w] );
    
    if( a > 0 ) {
      const IntegerType& outflow = G.currents[a];
      const IntegerType alpha = outflow - carry[v];
      
      ST.update( v, alpha);
      carry[w] += outflow;

      #ifdef VERBOSE
        cout << "alpha = " << to_double(alpha) << endl;
      #endif
      
    } else {
      const IntegerType& inflow = G.currents[-a];
      const IntegerType alpha = - (inflow + carry[v]);
      
      ST.update( v, alpha );
      carry[w] -= inflow;

      #ifdef VERBOSE
        cout << "alpha = " << to_double(alpha) << endl;
      #endif   
    }
  }
  

  // compute tree induced voltages
  compute_tree_induced_voltages_layered( ST, N );
  
  // asserts tree induced voltages are correct by checking ohm violations
  #ifdef VERBOSE
    cout << "Ohm Violations" << endl;
  #endif
  #if !defined(NDEBUG) || defined(VERBOSE)
    for( arc a : G.tree_edges ) {
      const IntegerType& pi_head = G.tree_induced_voltages[G.head(a)];
      const IntegerType& pi_tail = G.tree_induced_voltages[G.tail(a)];
      const IntegerType& r_a = N.arcdata[abs(a)].resistances;
      assert( r_a >= 0 );
      const IntegerType voltage_drop =  pi_tail - pi_head;
      const IntegerType alpha = a > 0 ? G.currents[a] : -G.currents[-a];
      const IntegerType ohm_violation = alpha*r_a - voltage_drop;
      #ifdef VERBOSE  
        if(ohm_violation != 0){
          cout << "pi_head: " << to_double(pi_head) << endl << "pi_tail: " << to_double(pi_tail) << endl;
          cout << endl << "voltage_drop: " << to_double(voltage_drop) << endl;
          cout << "ohm assertion" << endl;
          cout << "alpha * r_a = " << to_double(alpha * r_a) << endl;
          cout << "alpha: " << to_double(alpha) << endl;
          cout << "r_a: " << to_double(r_a) << endl;
          cout << a << " = (" << G.tail(a) << "," << G.head(a) << "): " << to_double(alpha) << " * " << to_double(r_a) << " + " << to_double(pi_head) << " - " << to_double(pi_tail) << " = " << to_double(ohm_violation) << endl;
        }
      #endif
      assert( ohm_violation < 1e-3 );
    }
  #endif
  
  #ifdef VERBOSE
    cout << "non-tree edges" << endl;
    for( arc a : G.non_tree_edges ) {
      assert( a > 0 );

      const IntegerType alpha = G.currents[a];
      const IntegerType pi_head = ST.query( G.head(a) , 0 );
      const IntegerType pi_tail = ST.query( G.tail(a) , 0 );
      const IntegerType voltage_drop =  pi_tail - pi_head;
      
      const IntegerType ohm_violation = alpha*G.resistances[a] - voltage_drop;
      cout << a << " = (" << G.tail(a) << "," << G.head(a) << "): " << to_double(alpha) << " * " << to_double(G.resistances[a])
      << " - (" << to_double(pi_head) << " - " << to_double(pi_tail) << ") = " << to_double(ohm_violation) << endl;
    }
  #endif 
  
}

/** \brief A sanity check to see if there are imbalances at any node
 * 
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * 
 * @return Returns true if the imbalance check fails
 * 
 */
template<typename IntegerType, typename RationalType>
bool 
imbalances_check( const SpanningTree<Graph<IntegerType, 
		  RationalType> , IntegerType, RationalType> &ST,
		  Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N
		)
{
  #if VERBOSE
    cout << "imbalance check" << endl;
  #endif
  
  Graph<IntegerType, RationalType>& G = ST.G;
  vector<IntegerType> imbalance( G.no_of_verticies+1,0 );
  
  for( auto aa : G.tree_edges ){
    if( aa > 0 ) {
      const IntegerType alpha = IntegerType( static_cast<long long>(( to_double( G.voltage_drop(aa) ) / to_double( N.arcdata[aa].resistances ) ) ) ) - G.f_0[aa];
      imbalance[G.head(aa)] += alpha;
      imbalance[G.tail(aa)] -= alpha;
      
      #ifdef VERBOSE
        cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.voltage_drop(aa)) << " <-> " << to_double(G.resistances[aa]*G.currents[aa]) << endl;
      #endif
      
    } else{
      const IntegerType alpha = -IntegerType( static_cast<long long>( to_double( G.voltage_drop(-aa))/to_double(N.arcdata[-aa].resistances) ) ) + G.f_0[-aa];
      imbalance[G.head(aa)] += alpha;
      imbalance[G.tail(aa)] -=  alpha;
      
      #ifdef VERBOSE
        cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << -to_double(G.currents[-aa]) << " <-> " << - to_double(G.voltage_drop(-aa))/to_double(G.resistances[-aa]) << endl;
      #endif
      
    }
  }
  
  for( auto aa : G.non_tree_edges ){
    assert( aa > 0 );
    const IntegerType alpha = G.currents[aa] - G.f_0[aa];
    imbalance[G.head(aa)] +=  alpha;
    imbalance[G.tail(aa)] -= alpha;
    #ifdef VERBOSE
    cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.currents[aa]) << " " << to_double(alpha) << endl;
    #endif
  }
  #ifdef VERBOSE
  cout << "large imbalances: " << endl;
  #endif
  for( node v = 1; v <= G.no_of_verticies; ++ v ){
    #ifdef VERBOSE
    if( RationalType(1e1) * to_double(imbalance[v]) > RationalType(1) ){
      cout << "imbalance of " << v << " is " << to_double(imbalance[v]) << endl;
    }
    #endif
    
  }
  
  #ifdef VERBOSE
  cout << "imbalances_check done " ;
  #endif
  return true;
}

// This function is never called (currently) ...

/** \brief Updates the currents of the tree edges
 * 
 * @param ST Spanning Tree
 * 
 */
template<typename IntegerType, typename RationalType>
void 
update_current_tree_edges(SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST, Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N)
{
  Graph<IntegerType, RationalType>& G = ST.G;
  vector<IntegerType> imbalance( G.no_of_verticies+1,0 );
  for( auto aa : G.non_tree_edges ) {
    assert( aa > 0 );
    const IntegerType alpha = G.currents[aa] - G.f[aa];
    imbalance[G.head(aa)] += alpha;
    imbalance[G.tail(aa)] -= alpha;
    #ifdef VERBOSE
    cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.currents[aa]) << " " << to_double(alpha) << endl;
    #endif
  }
  
  for( arc aa : boost::adaptors::reverse( G.tree_edges ) ){
    const node v = G.tail(aa);
    const node w = G.head(aa);
    assert( ST.depth[v] > ST.depth[w] );
    if( aa > 0 ) {
      const IntegerType imb = imbalance[v];
      const IntegerType outflow = imb + G.f[aa];
      const RationalType alpha = outflow - G.voltage_drop(aa)/N.arcdata[aa].resistances;
      ST.update( v, alpha );
      G.currents[aa] = outflow;
      imbalance[w] += imb;
      
      
      #ifdef VERBOSE
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.currents[aa]) << " " << to_double(alpha) << endl;
      #endif
    } else{
      const IntegerType imb = imbalance[v];
      const IntegerType outflow = imb - G.f[-aa];
      const RationalType alpha = outflow + ST.voltage_drop(-aa)/N.arcdata[-aa].resistances;
      ST.update( v, alpha );
      G.currents[-aa] = -outflow;
      imbalance[w] += imb;
      
      #ifdef VERBOSE
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << -to_double(G.currents[-aa]) << " " << to_double(alpha) << endl;
      #endif
    }
  }
  #ifdef VERBOSE
  cout << "imbalance of root: " << to_double(imbalance[ST.root]) << endl;
  #endif
  assert( fabs( to_double(imbalance[ST.root]) ) < 1e-3 );
}

/** \brief Runs a sanity-check to see if the value returned by query is same as the tree induced voltages
 * 
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * 
 */
template<typename IntegerType, typename RationalType>
void assert_query_v_equals_tree_ind_vol(SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType>& ST){
  compute_tree_currents_by_non_tree_edge(ST);
  
  Graph<IntegerType, RationalType>& G = ST.G;
  
  #ifdef VERBOSE
  cout << endl << " Current later... " << endl;
  for(unsigned int i=1; i <= G.no_of_edges; i++){
    cout << i << " : " << to_double(G.currents[i]) << endl;
  }
  
  #ifdef VERBOSE
  cout << "current updated; " << endl ;
  
  cout << "assert_query_v_equals_tree_ind_vol" << endl;
  #endif
  
  for(unsigned int v=1; v <= G.no_of_verticies; v++){
    IntegerType voltage = ST.query(v, 0);
    #ifdef VERBOSE
    cout << "voltage of " << v << " is " << to_double(voltage) << endl;
    #endif
  }
  #endif
  
  for(unsigned int v=1; v <= G.no_of_verticies; v++){
    const IntegerType voltage = ST.query(v, 0);
    const RationalType err = fabs( to_double(voltage - G.tree_induced_voltages[v]));
    if( err >= 1e-3){
      #ifdef VERBOSE
      cout << "vertex going wrong: " << v << " ,  margin of error: " << fabs(to_double(voltage - G.tree_induced_voltages[v])) << endl; 
      #endif
    }
    assert( err < 1e-3);
  }
  

}

/** \brief Runs the electrical flow problem
 * 
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * @param delta 
 * 
 * @see build_initial_state(ST)
 * @see compute_gap(G)
 * 
 */
template<typename IntegerType, typename RationalType>
int electrical_flow_problem(
          vector<RationalType>& x, 
			    vector<RationalType>& s,
			    SpanningTree<Graph< IntegerType, RationalType>, IntegerType, RationalType >& ST,
			    Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N,
			    const RationalType& minimum_efg, 
			    RationalType& TCN,
			    RationalType THRESHOLD_X,
			    RationalType THRESHOLD_S
 			  )
{
  

  Graph<IntegerType, RationalType>& G = ST.G;

  #ifdef VERBOSE
    cout << "electrical_flow_problem is called with minimum electrical flow gap: " << minimum_efg <<endl;
    cout << "G.no_of_edges: " <<  G.no_of_edges << endl << "Number of Tree Edges: " << G.tree_edges.size() << endl << "Non Tree Edges: " << G.non_tree_edges.size() << endl;
  #endif
  
  // clears R_a for a notin T
  G.resistances_accross_cycles.clear();

  // compute resistances across cycles
  for( auto a : G.non_tree_edges ){    
    G.resistances_accross_cycles.push_back(ST.compute_resistance_accross( a,N ));
    #ifdef VERBOSE 
      cout << "Resistances accross " << a << " =  " << to_double( ST.compute_resistance_accross( a, N ) ) << endl;
    #endif
  }
  
  // initializes the data structure: d_ext, d_drop
  build_initial_state(ST, N);

  // setting of probabilities
  for( auto a : G.non_tree_edges ){
    // removed and contracted arcs get probability zero
    if( x[a] < THRESHOLD_X || s[a] < THRESHOLD_S ){
      G.probability[a] = RationalType(0); 
    }
  }
  
 #ifdef VERBOSE
  cout << "Removed / Contracted Arcs" << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++){
  
    if(x[i] < THRESHOLD_X || s[i] < THRESHOLD_S){
    
      cout << i << " , ";
    }
  }
  cout << endl;
  
  cout << "non_tree_edges" << endl;
  for( auto a : G.non_tree_edges ){
    
    cout << a << " , " ;
    
  }
 #endif
  // compute intitial electrical flow gap
  RationalType gap = compute_gap_by_non_tree_currents(ST, N, x, s, THRESHOLD_X, THRESHOLD_S);

  // initializing boost discrete sampler
  ArcSampler sample_arc( G );


  assert( imbalances_check(ST, N) );  


  // a lot of couts...
  #ifdef VERBOSE   
    cout << "arc sampler called.." << endl;    
    cout << "currents: " << endl;
    
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << i << " : " << to_double(G.currents[i]) << "  " << G.head(i) << " , " << G.tail(i) << endl; 
    }
    
    cout << "demands: " << endl;
    for(unsigned int i=1; i<=G.no_of_verticies; i++){
      cout << i << " : " << G.demands[i] << endl; 
    }
    
    cout << endl << "Node Tree Corrospondance" << endl << endl;
    
    for(unsigned int i=1; i<=G.no_of_verticies; i++){
      cout << i << " : "; 
      for(auto tree_height: ST.node_tree_corrospondance[i]){
        cout << " [ " << tree_height.tree_index <<" , " << tree_height.above_d << " , " << " , " << tree_height.lca << " ] " << " , " ;
      }
      cout << endl;
    }
    
    cout << endl << "Tree Decomposition Information" << endl << endl;
    
    unsigned int count = 0;
    for(auto TDI: ST.tree_decomposition){
      cout << count << " : " ;
      cout << " [ " << TDI.size <<" , " << TDI.s << " , " << TDI.d << " ] " << " , " ;
      count++;
      cout << endl; 
    }
  #endif

  // these are for the gap computation after every m iterations
  unsigned int iteration = 0;
  unsigned int nextcheck = G.no_of_edges;
  unsigned int offset = 0;
  
  #ifdef VERBOSE
    cout << endl << "while begins now" << endl;
    cout << "initial gap: " << gap << " Minimal Electrical Flow Gap: " << minimum_efg << endl;
  #endif
  
  while(gap > minimum_efg){
    
    iteration++;
    
    // sample the non-tree edge
    unsigned int non_tree_edge_index = sample_arc();

#if !defined(NDEBUG) || defined(VERBOSE)
    const arc aa = G.non_tree_edges[non_tree_edge_index];
#endif
    
   #ifdef VERBOSE
    cout << "arc sampled = " << aa << endl;
   #endif
    
    assert( x[aa] >= THRESHOLD_X && s[aa] >= THRESHOLD_S );

    
    #ifdef VERBOSE
      cout << endl << " Current before update... " << endl;
      for(unsigned int i=1; i <= G.no_of_edges; i++){
        cout << i << " : " << to_double(G.currents[i]) << endl;
      }
      cout << endl << "current updation will take place now with non_tree_edge_index: " << non_tree_edge_index << endl;
    #endif
    
    // update current
    IntegerType alpha = compute_alpha(ST, N, non_tree_edge_index);
    do_update_with_values(ST, non_tree_edge_index, alpha);
    
    #ifdef VERBOSE
      cout << endl << " Current after update... " << endl;
      for(unsigned int i=1; i <= G.no_of_edges; i++){
        cout << i << " : " << to_double(G.currents[i]) << endl;
      }
    #endif
    

    // compute gap after every m many iterations
    if( iteration > nextcheck + offset){
      gap = compute_gap_by_non_tree_currents(ST, N, x, s, THRESHOLD_X, THRESHOLD_S);
      cout << " electrical flow gap = " << gap << endl;
      nextcheck += G.no_of_edges;
    }
  
  // end of while loop
  }

  // compute currents on tree from currents on non tree arcs
  compute_tree_currents_by_non_tree_edge( ST,N, x, THRESHOLD_X );      

  assert( imbalances_check(ST, N) );

  return iteration;
}
    
