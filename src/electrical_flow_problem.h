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
  
  std::vector<arc>& cyc = G.cycle[edge_index];
  
  for( auto y : cyc ){
    #ifdef VERBOSE
    cout << y <<  " ";
    #endif
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
  for(unsigned int i=1; i<=G.no_of_edges; i++)
  {
    cout << i << " : " << G.currents[i] << endl;
  }
  #endif
}


/** \brief does a cycle update of the current with the datastructure given the voltage drop
 * 
 * @param G0 The original graph 
 * @param ST The Spanning Tree
 * @param edge_index Edge Index of the arc being updated
 * @param voltage_drop Voltage drop accross the arc
 */
template<typename IntegerType, typename RationalType>
void update_current_DS(Graph<IntegerType, RationalType>& G0, 
		       SpanningTree<Graph<IntegerType, RationalType>,IntegerType, RationalType>& ST, 
		       unsigned int edge_index, 
		       IntegerType& voltage_drop
)
{   
  
  Graph<IntegerType, RationalType>& G = ST.G;
  
  
  const arc a = G.non_tree_edges[edge_index];
  
  const IntegerType& R_a = G.resistances_accross_cycles[edge_index];
  
  #ifdef VERBOSE
  cout << "Voltage Drop: " << to_double(voltage_drop) << endl;
  #endif
  const IntegerType& r_a = G.resistances[a];
  
  
  
  #ifdef VERBOSE
  cout << endl << "r_a = " << to_double(r_a) << endl;
  #endif
  const IntegerType& Delta = voltage_drop - G.currents[a]*r_a;
  RationalType alpha_unrounded = -to_double( Delta )/ to_double( R_a );
  
  IntegerType alpha = round_alpha(alpha_unrounded, Delta, R_a);
  
  
  
  ST.update( alpha, edge_index);
  
  G.currents[a] -= alpha;   
  
  
}

/** \brief does a cycle update of the curernt using the datastructure
 * 
 * @param G0 The Original Graph
 * @param ST The Spanning Tree
 * @param edge_index
 * 
 */
template<typename IntegerType, typename RationalType>
void update_current_DS(Graph<IntegerType, RationalType>& G0, 
		       SpanningTree<Graph<IntegerType, RationalType>,IntegerType, RationalType>& ST, 
		       unsigned int edge_index
)
{   
  
  Graph<IntegerType, RationalType>& G = ST.G;
  
  
  #ifdef VERBOSE
  cout << "UPDATE CURRENT DS: " << endl << "---------" << endl;
  cout << "no_of_edges: " << G.no_of_edges << endl;
  cout << "tree_edges: " << G.tree_edges.size() << endl;
  cout << "non_tree_edges: " << G.non_tree_edges.size() << endl;
  cout << "---------------------" << endl;
  for(auto a: G.non_tree_edges){
    cout << a << " , " ;
  }
  cout << endl << endl;
  #endif
  
  #ifdef VERBOSE
  cout << "update current called: " << edge_index << endl;  
  cout << "currents and resistances before updating: " << endl;
  for(unsigned int i=1; i<= G.no_of_edges; i++){
    cout << i << " : " << to_double(G.currents[i]) << " " << to_double(G.resistances[i]) << endl;
  }
  
  for(unsigned int i = 1; i <= G0.no_of_edges; i++){
    cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde) << endl;
    cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde) << endl;
    cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow3_tilde) << endl;
  }
  #endif
  
  const arc a = G.non_tree_edges[edge_index];
  
  
  #ifdef VERBOSE
  cout << endl << "Edge Index Sampled: " << edge_index << endl;
  cout << "the non tree edge: " << G.non_tree_edges[edge_index] << endl;
  #endif
  
  const IntegerType& R_a = G.resistances_accross_cycles[edge_index];
  const IntegerType& voltage_drop = ST.query(edge_index); // v_u - v_w;
  #ifdef VERBOSE
  cout << "Voltage Drop: " << to_double(voltage_drop) << endl;
  #endif
  const IntegerType& r_a = G.resistances[a];
  
  
  #ifdef VERBOSE
  cout << endl << "r_a = " << to_double(r_a) << endl;
  #endif
  
  
  
  IntegerType f_a = G.currents[a];
  if(f_a != 0 && r_a != 0){
    multiply_by_bit_shift(f_a, r_a);
  }
  else{
    f_a = 0;
  }
  
  const IntegerType& Delta = voltage_drop - f_a;
  
  #ifdef VERBOSE
  cout << endl << "currents = " << to_double(G.currents[a]) << endl;
  cout << "Delta : " << to_double(Delta) << endl;
  cout << "R_a: " << to_double(R_a) << endl;
  #endif
  RationalType alpha_unrounded = -to_double( Delta )/ to_double( R_a );
  #ifdef VERBOSE
  cout << "unround alpha = " << alpha_unrounded << endl;
  #endif
  
  #ifdef VERBOSE
  cout << endl << "alpha before rounding = " << alpha_unrounded << endl;
  #endif
  
  
  IntegerType alpha = round_alpha(alpha_unrounded, Delta, R_a);
  
  
  #ifdef VERBOSE
  cout << "alpha after rounding (update_current_DS) = " << to_double(alpha) << endl;
  #endif
  
  #ifdef VERBOSE
  
  cout << endl << "calculating alpha for all non-tree edges" << endl;
  
  for( unsigned int i = 0; i < G.non_tree_edges.size(); i++ ){
    cout << endl << " ------ " << i << " ---------- " << endl;
    const IntegerType& D =   ST.query(i); //G.voltage_drop(i);
    cout << "D = "  << to_double(D) << endl;
    const IntegerType& r_a = G.resistances[i];
    cout << "r_a = " << to_double(r_a) << endl;
    const IntegerType& f_a = G.currents[i];
    cout << endl << "f_a = " << to_double(f_a) << endl;
    const IntegerType Delta = f_a*r_a - D;
    
    #ifdef VERBOSE  
    cout << endl << "Delta = " << to_double(Delta) << endl;
    cout << endl << " alpha = " << to_double(alpha) << endl;
    #endif    
  }
  
  #endif
  
  #ifdef VERBOSE
  cout << "alpha: " << to_double(alpha) << " = " << " ( "  << -to_double(voltage_drop) << " - " << " ( " << to_double(ST.G.currents[a])
  << " * " << to_double(ST.G.resistances[a]) << " ) " << " ) " << " / " << to_double(ST.G.resistances_accross_cycles[edge_index]) << endl;
  #endif
  
  
  #ifdef VERBOSE
  cout << "will update now" << endl;
  
  IntegerType DeltaSquare = Delta*Delta;
  const RationalType contribution = to_double( DeltaSquare )/ to_double( r_a );
  cout << "contribution  " << contribution <<  " * " << G.no_of_edges << endl;
  cout << contribution << " = " << to_double(Delta*Delta) << " div  " << to_double(r_a) << endl;
  cout << "Delta: " << to_double(Delta) << endl;
  cout << "Delta_square: " << to_double(Delta*Delta) << endl; 
  #endif
  
  
  {
    
    #ifdef VERBOSE
    cout << "called now!!" << endl;
    #endif
    ST.update( alpha, edge_index);
    
    G.currents[a] -= alpha;  
    
  }
  
  #ifdef VERBOSE
  cout << "currents and resistances after update: " << endl;
  for(unsigned int i=1; i<= G.no_of_edges; i++){
    cout << i << " : " << to_double(G.currents[i]) << " "  << to_double(G.resistances[i]) << endl;
  }
  #endif
  
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
	    Graph<IntegerType, RationalType> &G
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
    const IntegerType& r_a = G.resistances[i];
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
compute_tree_currents_by_non_tree_edge(Graph<IntegerType, RationalType>& G0, 
				       SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST )
{
  Graph<IntegerType, RationalType>& G = ST.G;
  vector<IntegerType> imbalance( G.no_of_verticies+1,0 );
  vector<IntegerType> carry( G.no_of_verticies+1,0 );
  for( auto aa : G.non_tree_edges ){
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
    unsigned int edge_index_of_G0 = abs(aa)/3 + 1;
    if(abs(aa) % 3 == 0) edge_index_of_G0--;
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
  #if VERBOSE
  cout << "imbalance of root: " << to_double( imbalance[ST.root] ) << endl;
  cout << "Tree Condition Number: " << ST.calculate_tree_condition_number();
  #endif
  
  assert( RationalType(1e3) * fabs( to_double( imbalance[ST.root] ) )  < RationalType(1) );
  
  compute_tree_induced_voltages_layered(ST);
}

/** \brief Computes the gap using using non-tree currents
 * 
 * @param x The Primal Solution
 * @param ST The Spanning Tree
 * @param ST_auxiallary Spanning Tree of Auxiliary Graph
 * 
 * @return Returns duality gap
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType
compute_gap_by_non_tree_currents( vector<RationalType>& x,
				  SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> & ST, 
				  SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType>& ST_auxiallary) 
{  
  
  Graph<IntegerType, RationalType>& G = ST_auxiallary.G;
  
  #ifdef VERBOSE
  cout << "compute gap by non-tree currents" <<endl;
  #endif
  
  #ifdef VERBOSE
  cout << "currents in non_tree_edges: " <<  endl;
  for(auto a: G.non_tree_edges) {
    cout << to_double(G.currents[a]) << " , " ;
  }
  cout << endl;
  #endif
  
  
  #ifdef VERBOSE
  cout << "arc\t\tflow\t\t\tresistance\t\tbatteries\t\tb-flow\t\tvoltage\t\tgap-contrib.\t\tKVL viol." << endl;
  #endif
  
  RationalType gap = 0;
  unsigned int edge_index = 0;
  for( auto i : G.non_tree_edges ) {
    if(x[i] != 0) {
      const IntegerType& D = ST.G.tree_induced_voltages[G.tail(i)] - ST.G.tree_induced_voltages[G.head(i)]; // G.voltage_drop(i);
      
      #ifdef VERBOSE
      cout << endl << "edge_index = " << edge_index << endl;
      #endif
      #ifdef VERBOSE
      
      cout << "D = " << to_double(D) << " " << "ST - " << to_double(ST.query(edge_index)) << endl;
      #endif
      edge_index++;
      
      
      const RationalType& r_a = G.rounded_resistances_auxiallary[i];
      
      const RationalType& f_a = G.unrounded_currents[i];
      
      #ifdef VERBOSE
      cout << "f_a = " << to_double(f_a) << endl;
      #endif
      
      const RationalType Delta = f_a*r_a - to_double(D);
      
      
      const RationalType DeltaSquare = Delta*Delta;
      const RationalType contribution = ( DeltaSquare )/ ( r_a );
      
      #ifdef VERBOSE   
      const RationalType& g_a = G.batteries[i];
      const IntegerType KVL_violation = D - f_a*r_a;
      cout << i << fixed << ":\t\t" << to_double(f_a) << "\t\t" << to_double(r_a) << "\t\t" << g_a << "\t\t" << to_double(D) << "\t\t" << contribution << "\t\t" << to_double(KVL_violation) << endl;
      #endif
      gap += contribution;
    }
  }
  
  #ifdef VERBOSE
  cout << "gap = " << gap << endl;
  #endif
  return gap; 
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
compute_gap_by_non_tree_currents(Graph<IntegerType, RationalType>& G0, 
				 SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> & ST 
)
{  
  Graph<IntegerType, RationalType>& G = ST.G;
  
  #ifdef VERBOSE
  cout << "compute gap by non-tree currents" <<endl;
  #endif
  
  #ifdef VERBOSE
  cout << "currents in non_tree_edges: " <<  endl;
  for(auto a: G.non_tree_edges){
    cout << to_double(G.currents[a]) << " , " ;
  }
  cout << endl;
  #endif
  
  #ifdef VERBOSE
  cout << "arc\t\tflow\t\t\tresistance\t\tbatteries\t\tb-flow\t\tvoltage\t\tgap-contrib.\t\tKVL viol." << endl;
  #endif
  
  RationalType gap = 0;
  unsigned int edge_index = 0;
  unsigned int edge_index_with_max_query = 0;
  IntegerType max_query(0);
  for( auto i : ST.G.non_tree_edges ){
    
    #ifdef VERBOSE
    cout << i << endl;
    #endif
    
    const IntegerType& D = ST.query(edge_index); 
    if(D > max_query){
      max_query = D;
      edge_index_with_max_query = edge_index;
    }
    
    #ifdef VERBOSE
    cout << "D = " << to_double(D) << " " << "ST - " << to_double(ST.query(edge_index)) << endl;
    #endif
    edge_index++;
    
    const IntegerType& r_a = G.resistances[i];
    
    IntegerType f_a = G.currents[i];
    
    #ifdef VERBOSE
    cout << "f_a = " << to_double(f_a) << endl;
    cout << "r_a = " << to_double(r_a) << endl;
    #endif
    
    if(f_a != 0 && r_a !=0){
      multiply_by_bit_shift(f_a, r_a);
    }
    else{
      f_a = 0;
    }
    const IntegerType Delta = f_a - D;
    
    #ifdef VERBOSE
    cout << endl << "Delta = " << to_double(Delta) << endl;
    #endif
    
    const IntegerType DeltaSquare = Delta*Delta;
    const RationalType contribution = to_double( DeltaSquare )/ to_double( r_a );
    #ifdef VERBOSE   
    const RationalType& g_a = G.batteries[i];
    const IntegerType KVL_violation = D - f_a*r_a;
    cout << i << fixed << ":\t\t" << to_double(f_a) << "\t\t" << to_double(r_a) << "\t\t" << g_a << "\t\t" << to_double(D) << "\t\t" << contribution << "\t\t" << to_double(KVL_violation) << endl;
    #endif
    
    gap += contribution;
    
  }
  
  #ifdef VERBOSE
  cout << "edge_index_with_max_query = " << edge_index_with_max_query << endl;
  #endif
  
  assert(edge_index_with_max_query >= 0);
  update_current_DS(G0, ST, edge_index_with_max_query, max_query);
  
  #ifdef VERBOSE
  cout << "gap = " << gap << endl;
  #endif
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
compute_tree_induced_voltages_layered( SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> & ST)

{
  Graph<IntegerType, RationalType>& G = ST.G;
  #ifdef VERBOSE
  cout<<"computing tree induced voltages..." << endl;
  #endif
  
  
  assert( G.head( G.tree_edges.front() ) == G.root );
  G.tree_induced_voltages[G.root] = 0;
  for( auto a : G.tree_edges ) {
    assert( a != 0 );
    assert( ST.depth[G.tail(a)] > ST.depth[G.head(a)] );
    
    #ifdef VERBOSE
    cout << a << " = (" << G.tail(a) << "," << G.head(a) << ")" << endl;
    #endif
    if( a > 0 ) {
      
      IntegerType f_a = G.currents[a];
      if(f_a != 0 && G.resistances[a] !=0){
	multiply_by_bit_shift(f_a, G.resistances[a]);
      }
      else{
	f_a = 0;
      }
      G.tree_induced_voltages[G.tails[a]] = G.tree_induced_voltages[G.heads[a]] + f_a; 
      
    } else {
      const arc rev = -a;
      
      IntegerType f_a = G.currents[rev];
      if(f_a != 0 && G.resistances[rev] !=0){
	multiply_by_bit_shift(f_a, G.resistances[rev]);
      }
      else{
	f_a = 0;
      }
      
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
    #endif
    
    
    #ifdef VERBOSE
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
    
    assert(  fabs( to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] )  ) < 1e-3 );
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
		    RationalType> & ST
)
{
  
  Graph<IntegerType, RationalType>& G = ST.G;
  ST.clear_d();
  
  #ifdef VERBOSE
  cout << "build_initial_state() " << endl;
  cout << "currents ... " << endl;
  for(unsigned int i=1; i <= G.no_of_edges; i++) {
    cout << i << " : " << to_double(G.currents[i]) << endl;
  }
  #endif
  
  #ifdef VERBOSE
  cout << "init datastructure with initial currents" << endl;
  cout << "tree edges first" << endl;
  #endif
  vector<IntegerType> carry( G.no_of_verticies+1,0 );
  for( arc a : boost::adaptors::reverse( G.tree_edges ) ) {
    
    const node v = G.tail(a);
    const node w = G.head(a);
    assert( ST.depth[v] > ST.depth[w] );
    
    if( a > 0 ) {
      const IntegerType& outflow = G.currents[a];
      const IntegerType alpha = outflow - carry[v];
      #ifdef VERBOSE
      cout << "alpha = " << to_double(alpha) << endl;
      #endif
      ST.update( v, alpha);
      carry[w] += outflow;
    } else {
      const IntegerType& inflow = G.currents[-a];
      const IntegerType alpha = -(inflow + carry[v]);
      #ifdef VERBOSE
      cout << "alpha = " << to_double(alpha) << endl;
      #endif
      ST.update( v, alpha );
      carry[w] -= inflow;
    }
  }
  
  //cache_tree_induced_voltages( ST );
  compute_tree_induced_voltages_layered( ST );
  #if !defined(NDEBUG) || defined(VERBOSE)
  for( arc a : G.tree_edges ) {
    const IntegerType& pi_head = G.tree_induced_voltages[G.head(a)];
    const IntegerType& pi_tail = G.tree_induced_voltages[G.tail(a)];
    const IntegerType& r_a = G.resistances[abs(a)];
    assert( r_a > 0 );
    const IntegerType voltage_drop =  pi_tail - pi_head;
    const IntegerType alpha = a > 0 ? G.currents[a] : -G.currents[-a];
    const IntegerType ohm_violation = alpha*r_a - voltage_drop;
    
    if(ohm_violation != 0){
      cout << "pi_head: " << to_double(pi_head) << endl << "pi_tail: " << to_double(pi_tail) << endl;
      cout << endl << "voltage_drop: " << to_double(voltage_drop) << endl;
      cout << "ohm assertion" << endl;
      cout << "alpha * r_a = " << to_double(alpha * r_a) << endl;
      cout << "alpha: " << to_double(alpha) << endl;
      cout << "r_a: " << to_double(r_a) << endl;
      cout << a << " = (" << G.tail(a) << "," << G.head(a) << "): " << to_double(alpha) << " * " << to_double(r_a) << " + " << to_double(pi_head) << " - " << to_double(pi_tail) << " = " << to_double(ohm_violation) << endl;
    }
    
    assert( ohm_violation == 0 );
  }
  #endif
  #ifdef VERBOSE
  cout << "non-tree edges" << endl;
  #endif
  
  #ifdef VERBOSE
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
		  RationalType> , IntegerType, RationalType> &ST)
{
  #if VERBOSE
  cout << "imbalance check" << endl;
  #endif
  Graph<IntegerType, RationalType>& G = ST.G;
  vector<IntegerType> imbalance( G.no_of_verticies+1,0 );
  for( auto aa : G.tree_edges ){
    if( aa > 0 ) {
      const IntegerType alpha = IntegerType( static_cast<long long>(( to_double( G.voltage_drop(aa) ) / to_double( G.resistances[aa] ) ) ) ) - G.f_0[aa];
      imbalance[G.head(aa)] += alpha;
      imbalance[G.tail(aa)] -= alpha;
      #ifdef VERBOSE
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.voltage_drop(aa)) << " <-> " << to_double(G.resistances[aa]*G.currents[aa]) << endl;
      #endif
      
    } else{
      const IntegerType alpha = -IntegerType( static_cast<long long>( to_double( G.voltage_drop(-aa))/to_double(G.resistances[-aa]) ) ) + G.f_0[-aa];
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

/** \brief Updates the currents of the tree edges
 * 
 * @param ST Spanning Tree
 * 
 */
template<typename IntegerType, typename RationalType>
void 
update_current_tree_edges(SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST)
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
      const RationalType alpha = outflow - G.voltage_drop(aa)/G.resistances[aa];
      ST.update( v, alpha );
      G.currents[aa] = outflow;
      imbalance[w] += imb;
      
      
      #ifdef VERBOSE
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.currents[aa]) << " " << to_double(alpha) << endl;
      #endif
    } else{
      const IntegerType imb = imbalance[v];
      const IntegerType outflow = imb - G.f[-aa];
      const RationalType alpha = outflow + ST.voltage_drop(-aa)/G.resistances[-aa];
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
int 
electrical_flow_problem(vector<RationalType>& x,  
			Graph<IntegerType, RationalType>& G0, 
			SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST, 
			const RationalType& delta, RationalType& TCN
		       )
{
  
  #ifdef VERBOSE
  cout << "electrical_flow_problem   " << "delta: " << delta <<endl;
  #endif
  
  Graph<IntegerType, RationalType>& G = ST.G;
  
  #ifdef VERBOSE
  cout << "G.no_of_edges: " <<  G.no_of_edges << endl << "Number of Tree Edges: " << G.tree_edges.size() << endl << "Non Tree Edges: " << G.non_tree_edges.size() << endl;
  #endif
  ST.clear_d();
  
  #ifdef VERBOSE
  cout << endl << "D_EXT and D_DROP " << endl;  
  for(unsigned int v = 0; v < ST.d_ext.size(); v++){
    cout << to_double(ST.d_ext[v]) << " = " << to_double(ST.d_drop[v]) << endl;
  }
  
  cout << endl << "electrical_flow_problem " << endl;     
  cout<< "Non tree edges size: " << G.non_tree_edges.size();
  cout<< endl;
  for(auto a: G.non_tree_edges){
    cout<< a << ", ";
  }
  cout << endl;
  cout << endl << "Tree Edges Size: " << G.tree_edges.size();
  cout<< endl;
  cout << "Number of Verticies: " << G.no_of_verticies << endl;
  cout << "Number of Edges: " << G.no_of_edges << endl;
  #endif
  
  
  G.resistances_accross_cycles.clear();
  
  for(unsigned int i =0; i<G.non_tree_edges.size(); i++){    
    G.resistances_accross_cycles.push_back(ST.compute_resistance_accross(G.non_tree_edges[i]));
    #ifdef VERBOSE 
    cout << i << " Resistances accross  =  " << to_double(ST.compute_resistance_accross(G.non_tree_edges[i])) << endl;
    #endif
  }
  
  
  ArcSampler sample_arc( G );
  
  
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
  
  
  #ifdef VERBOSE 
  cout << endl << "Initial state would be built " << endl;
  #endif
  
  build_initial_state(ST);
  
  #ifdef VERBOSE
  cout << endl << "Initial state built " << endl;
  #endif
  
  RationalType gap = compute_gap_by_non_tree_currents(G0, ST);
  
  
  assert( imbalances_check(ST) );
  #ifdef VERBOSE
  cout << endl << " assertion done ";
  #endif
  
  
  
  unsigned int iteration = 0;
  unsigned int nextcheck = G.no_of_edges;
  unsigned int offset = 0;
  
  #ifdef VERBOSE
  cout << endl << "while begins now" << endl;
  cout << gap << " " << delta << endl;
  #endif
  
  
  #ifdef VERBOSE
  IntegerType r_max(0);
  for(unsigned int  i = 1; i <= G.no_of_edges; i++){
    if(r_max < G.resistances[i]){
      r_max = G.resistances[i];
    }
  }
  RationalType non_tree_stretch(0);
  
  for(unsigned int i = 0; i < G.non_tree_edges.size() ; i++){
    non_tree_stretch += to_double(G.resistances[G.non_tree_edges[i]]) / to_double(G.resistances_accross_cycles[i]) ;
  }
  #endif
  
  
  
  #ifdef VERBOSE
  cout << gap  << " < - > " << delta << endl;
  #endif
  
  
  while(gap > delta){
    iteration++;
    #ifdef VERBOSE
    cout << "checking if the currents are integral" << endl;
    cout << "Edge" << "\t" << "Currents " << "\t" << "Resistances" << "\t" << "Batteries" << endl;
    for(unsigned int i = 1; i <= G.no_of_edges; i++) {
      cout << i << "\t" << to_double(G.currents[i]) << "\t" << to_double(G.resistances[i]) << "\t" << to_double(G.batteries[i]) << endl;
    }
    cout << endl;
    #endif
    
    #ifdef VERBOSE
    cout << endl << "while loop for electrical_flow_problem has begun....  gap = " << gap << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "Graph: " << endl;
    for(unsigned int a =1; a<= G.no_of_edges; a++){      
      cout << a << " =  (" << G.head(a) << " , " << G.tail(a) << ")" << endl; 
    }
    #endif
    
    #ifdef VERBOSE
    cout << endl << "Tree: " << endl;
    for(unsigned int v = 1; v <= G.no_of_verticies; v++){
      cout << v << " : " ;
      for(auto a: ST.tree_incident_edges[v]){
	cout << a << " , " << endl;
      }
    }
    cout << endl;
    #endif
    
    unsigned int non_tree_edge_index = sample_arc();
    
    
    #if !defined(NDEBUG) || defined(VERBOSE)
    const arc a = G.non_tree_edges[non_tree_edge_index];
    #endif
    assert( a > 0 );
    
    
    #ifdef VERBOSE
    cout << endl << " Current before... " << endl;
    for(unsigned int i=1; i <= G.no_of_edges; i++){
      cout << i << " : " << to_double(G.currents[i]) << endl;
    }
    #endif
    
    
    #ifdef VERBOSE 
    cout << endl << "current updation will take place now with non_tree_edge_index: " << non_tree_edge_index << endl;
    #endif
    
    update_current_DS(G0, ST, non_tree_edge_index);
    
    #ifdef VERBOSE
    cout << endl << " Current after... " << endl;
    for(unsigned int i=1; i <= G.no_of_edges; i++){
      cout << i << " : " << to_double(G.currents[i]) << endl;
    }
    #endif
    
    #ifdef VERBOSE
    cout << endl << "Current updated " << endl;
    #endif
    
    if( iteration > nextcheck + offset){
      #ifdef ElectricalFlowGap
      cout << gap << endl;
      #endif
      gap = compute_gap_by_non_tree_currents(G0, ST);
      
      
      nextcheck += G.no_of_edges;
      
      #ifdef VERBOSE
      cout << "electrical_flow_problem gap: " << gap  << " <-> " << delta << endl;  
      #endif
      
      #ifdef VERBOSE 
      bool exit = true; 
      
      long double alpha_unrounded = 0.0;
      vector<RationalType> unrounded_alphas;
      unrounded_alphas.clear();
      for( unsigned int i = 0; i < G.non_tree_edges.size(); i++){
	#ifdef VERBOSE
	cout << endl << " ------ " << i << " ---------- " << endl;
	#endif
	const IntegerType D =   ST.query(i); 
	#ifdef VERBOSE
	cout << "D = "  << to_double(D) << endl;
	#endif
	arc a = G.non_tree_edges[i]; 
	const IntegerType& r_a = G.resistances[a];
	const IntegerType& R_a = G.resistances_accross_cycles[i];
	#ifdef VERBOSE
	cout << "r_a = " << to_double(r_a) << endl;
	#endif
	const IntegerType& f_a = G.currents[a];
	#ifdef VERBOSE
	cout << endl << "f_a = " << to_double(f_a) << endl;
	#endif
	const IntegerType Delta = f_a*r_a - D;
	alpha_unrounded = to_double( Delta ) / to_double( R_a );
	unrounded_alphas.push_back(alpha_unrounded);
	#ifdef VERBOSE
	cout << "alpha (unrounded) = " << alpha_unrounded << endl;
	#endif
	
	
	#ifdef VERBOSE
	cout << endl << "currents = " << to_double(G.currents[a]) << endl;
	cout << "Delta : " << to_double(Delta) << endl;
	#endif
	
	
	
	if(alpha_unrounded >= RationalType(1)/RationalType(2) || alpha_unrounded < -RationalType(1)/RationalType(2)){
	  exit = false;
	  break;
	}
	
	#ifdef VERBOSE  
	cout << endl << "Delta = " << to_double(Delta) << endl;
	cout << endl << " alpha_unrounded = " << alpha_unrounded << endl;
	#endif    
      }
      
      if(exit == false){
	#ifdef VERBOSE
	cout << endl << "does not exit through this route. " << endl;
	#endif
      }
      
      if(exit == true){
	RationalType val(0);
	
	for(unsigned int i = 0; i < G.non_tree_edges.size(); i++){
	  arc a = G.non_tree_edges[i];
	  IntegerType R_a = G.resistances_accross_cycles[i];
	  IntegerType r_a = G.resistances[a];
	  
	  val += to_double(R_a) * to_double(R_a)/to_double(r_a);
	}
	
	
	if(RationalType(4) * gap >  val){
	  cout << gap << "  " << val << endl;
	  
	  for(auto a : unrounded_alphas){
	    cout << a << endl;
	  }
	  
	}
	
	assert (RationalType(4) * gap <=  val);
	
	
	#ifndef NDEBUG
	cout << "returned exit: " << endl;
	#endif
	
	return iteration;
      }
      #endif
      
    }
  }
  #ifdef VERBOSE
  cout << "number of electrical iterations: " << iteration << endl;
  #endif
  compute_tree_currents_by_non_tree_edge(G0, ST );
  
  assert( imbalances_check(ST) );
  
  
  return iteration;
}