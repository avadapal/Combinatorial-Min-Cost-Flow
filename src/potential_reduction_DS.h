
#include "potential_reduction.h"
#include "electrical_flow_problem_DS.h"

using namespace std;



// template<typename Network>
// void
// primal_feasibility_check(Network& N)
// {
//   typedef typename Network::GraphType Graph;
//   typedef typename Network::RationalType RationalType;
//   Graph& G0 = N.G;
//   
//   for(unsigned int v = 1; v <= G0.no_of_vertices; v++)
//   {
// #ifdef VERBOSE
//     cout << v << ": ";
// #endif
//     RationalType flow_into_v = 0;
//     for(auto a: G0.incident_edges[v])
//     {
//       
//       RationalType real_direction = a;
//       
//       if(N.arcdata[abs(a)].direction !=0)
//       {
//       real_direction = a * N.arcdata[abs(a)].direction;
//       }
//       else
//       {
// 	real_direction = a;
//       }
//       
//       if(real_direction > 0)
//       {
// 	flow_into_v -= N.arcdata[abs(a)].xlower;
//       }
//       else
//       {
// 	flow_into_v -= (N.arcdata[abs(a)].capacity- N.arcdata[abs(a)].xlower);
//       }
//       if(a > 0)
//       {
// 	flow_into_v -= N.arcdata[abs(a)].infeasibility;
//       }
//       else
//       {
// 	flow_into_v += N.arcdata[abs(a)].infeasibility;
//       }
//       
//       
//     }
//     
// #ifdef VERBOSE
//     cout << to_double(flow_into_v) << " < - > " << to_double(N.nodedata[v].demand) << endl;
// #endif
//     assert(fabs((to_double(flow_into_v) - to_double(N.nodedata[v].demand))) < 1e-4);
//   }
//   
// }

/** \brief Does a primal and dual sanity check, i.e. checks if the primal and dual constraints are satisfied
 * 
 * @param G0 The Original Graph G0
 * @param G The Auxilliary Graph G
 * @param x The primal solution x
 * @param y The Dual solution y
 * @param s The Dual slacks 
 */
template<typename IntegerType, typename RationalType>
void 
sanity_check_DS( Graph<IntegerType, RationalType>& G0, 
		      const Graph<IntegerType, RationalType> &G, 
		 std::vector<RationalType> &x,  
		 std::vector<RationalType> &y, 
		 std::vector<RationalType> &s
	       )
{
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
  //The Primal Constraint
  primal_sanity_check( G, x );
  //The Dual Constraint
  dual_sanity_check( G, y, s, x );
}

// template<typename Network, typename RationalType>
// void 
// sanity_check_DS(Network& N)
// {
//   typedef typename Network::GraphType Graph;
//   typedef typename Network::RationalType RationalType;
//   Graph& G0 =  N.G;
//   vector<RationalType> x(G0.no_of_edges * 3 + 1, 0);
//   unsigned int edge_index_of_G0 = 0;
//   for(unsigned int i = 1; i<= G0.no_of_edges; i+=3)
//   {
//     edge_index_of_G0++;
//     x[i] = N.arcdata[edge_index_of_G0].xlower;
//     x[i+1] = N.arcdata[edge_index_of_G0].capacity - N.arcdata[edge_index_of_G0].xlower;
//     x[i+2] = N.arcdata[edge_index_of_G0].infeasibility;
// //    s[i] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s1;
// //    s[i+1] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s2;
// //    s[i+2] = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s3;
//   }
// //   for(unsigned int i =1; i <= G0.no_of_edges; i++)
// //   {
// //     y[G0.tails[i]] = G0.original_auxiliary_transformed_corrospondence[i].y1;
// //     y[G0.heads[i]] = G0.original_auxiliary_transformed_corrospondence[i].y2;
// //     y[G0.no_of_verticies + i] = G0.original_auxiliary_transformed_corrospondence[i].y3;
// //   }
//   //The Primal Constraint
//   primal_sanity_check( G, x );
//   //The Dual Constraint
// //  dual_sanity_check( G, y, s, x );
// }
// 

template<typename Network>
long double 
calculate_potential_function_new( Network& N, typename Network::RationalType q)
{
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
	  
  Graph& G0 = N.G;
  
  RationalType x_t_s = 0;
  long double ln_x_s =0;
  unsigned int check_m = 0;
  q = RationalType(3) * G0.no_of_edges + ceil(sqrt(G0.no_of_edges * RationalType(3)));
  for(unsigned int i = 1; i<=G0.no_of_edges; i++)
  { 
    const auto& data = N.arcdata[i];
    //if(data.is_invalid) continue;
    if(data.xlower > 0)
    {
      const RationalType xasa1 = data.xlower * data.slower;
      assert( xasa1 > 0 );
      x_t_s+= xasa1;
      ln_x_s+= log( xasa1 );
      ++check_m;
    }
    
    if(data.capacity - data.xlower > 0)
    {
      const RationalType xasa2 = (data.capacity - data.xlower) * data.supper;
      assert( xasa2 > 0 );
      x_t_s+= xasa2;
      ln_x_s+= log( xasa2 );
      ++check_m;
    }
    
    if(data.infeasibility > 0)
    {
      const RationalType xasa3 = data.infeasibility * data.sroof;
      assert( xasa3 > 0 );
      x_t_s+= xasa3;
      ln_x_s+= log( xasa3 );
      ++check_m;
    }    
  }
  
  unsigned int m = (G0.no_of_edges * RationalType(3));
  
  #ifdef VERBOSE
  cout << m << " " << check_m << endl;
  #endif
  //m = check_m;
  assert( m == check_m );
  const long double mlogm = m*log(m);
  const long double potential = q*log(x_t_s) - ln_x_s - mlogm;
  #ifdef VERBOSE
  cout << "x_t_s: " << x_t_s << endl;
  cout<<"Potential: " << q << " * " << natural_logarithm(x_t_s) << " - " << ln_x_s << " - " << mlogm << " = " << potential << endl;
  #endif
  return potential; 
}




/** \brief calculates the potential function
 * 
 * @param G The Graph on which the min-cost-flow algorithm is being carried out
 * @param x The primal solution x
 * @param s The dual slack variables s
 * @param q The parameter q
 * 
 * @return the value of the potential
 * 
 */
template<typename IntegerType, typename RationalType>
long double 
calculate_potential_function_DS( const Graph<IntegerType, 
				 RationalType> &G0, 
				 unsigned int no_of_edges_removed,  
				 unsigned int q
			       )
{
  
  RationalType x_t_s = 0;
  long double ln_x_s =0;
  unsigned int check_m = 0;
  for(unsigned int i = 1; i<=G0.no_of_edges; i++)
  {
    const auto& data = G0.original_auxiliary_transformed_corrospondence[i];
    if(data.is_invalid) continue;
    if(data.x1 > 0)
    {
      const RationalType xasa1 = data.x1 * data.s1;
      assert( xasa1 > 0 );
      x_t_s+= xasa1;
      ln_x_s+= log( xasa1 );
      ++check_m;
    }
    
    if(data.x2 > 0)
    {
      const RationalType xasa2 = data.x2 * data.s2;
      assert( xasa2 > 0 );
      x_t_s+= xasa2;
      ln_x_s+= log( xasa2 );
      ++check_m;
    }
    
    if(data.x3 > 0)
    {
      const RationalType xasa3 = data.x3 * data.s3;
      assert( xasa3 > 0 );
      x_t_s+= xasa3;
      ln_x_s+= log( xasa3 );
      ++check_m;
    }    
  }
  
  unsigned int m = (G0.no_of_edges * 3) - no_of_edges_removed;
  
  #ifdef VERBOSE
  cout << m << " " << check_m << endl;
  #endif
  m = check_m;
  assert( m == check_m );
  const long double mlogm = m*log(m);
  const long double potential = q*log(x_t_s) - ln_x_s - mlogm;
  #ifdef VERBOSE
  cout << "x_t_s: " << x_t_s << endl;
  cout<<"Potential: " << q << " * " << natural_logarithm(x_t_s) << " - " << ln_x_s << " - " << mlogm << " = " << potential << endl;
  #endif
  return potential; 
}

/** \brief Gets the corresponding Delta-Wye resistances from the un-rounded resistances of the auxilliary graph
 * 
 * @param G0 The Original Graph
 * 
 */
template<typename IntegerType, typename RationalType>
void 
get_resistances_for_delta_wye_transformed_graph_DS(Graph<IntegerType, RationalType>& G0)
{
  for(unsigned int edge_index_of_G0 =1; edge_index_of_G0 <= G0.no_of_edges; edge_index_of_G0++){
    
    RationalType r1 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_aux; 
    RationalType r2 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_aux; 
    RationalType r3 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_aux; 
    
    RationalType x1 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1;
    RationalType x2 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2;
    RationalType x3 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3;
    
    if(x1 != 0 && x2 != 0 && x3 != 0){
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_delta_wye = (r2*r3)/(r1 + r2 + r3);
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_delta_wye = (r1*r3)/(r1 + r2 + r3);
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_delta_wye = (r2*r1)/(r1 + r2 + r3);
    }
    
    if(x1 == 0 && x2 != 0 && x3 != 0){  
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_delta_wye = 0; 
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_delta_wye = r3; 
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_delta_wye = r2; 
    }
    if(x1 != 0 && x2 == 0 && x3 != 0){     
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_delta_wye = 0; 
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_delta_wye = r3; 
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_delta_wye = r1;
    }    
    if(x1 != 0 && x2 != 0 && x3 == 0){
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_delta_wye  = r2; 
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_delta_wye = r1;      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_delta_wye = 0; 
    }   
  }
}

/** \brief Rounds the Initial flow for the electrical flow network
 * 
 * @param G0 The Original Graph G0
 * @param phi The Scaling Factor
 * 
 */
template<typename IntegerType, typename RationalType>
void 
round_f_DS(Graph<IntegerType, RationalType>&G0, 
		const RationalType& phi 
	       )
{
  
  for(unsigned int i=1; i<=G0.no_of_edges; i++){
    if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    
    if(G0.original_auxiliary_transformed_corrospondence[i].x1 != 0){
      round_flow( phi,G0.original_auxiliary_transformed_corrospondence[i].battery1_aux, G0.original_auxiliary_transformed_corrospondence[i].r1_tilde_aux, G0.original_auxiliary_transformed_corrospondence[i].f1_aux );  
    }
    else{
      G0.original_auxiliary_transformed_corrospondence[i].f1_aux = 0;
    }
    
    if(G0.original_auxiliary_transformed_corrospondence[i].x2 != 0){
      round_flow( phi,G0.original_auxiliary_transformed_corrospondence[i].battery2_aux, G0.original_auxiliary_transformed_corrospondence[i].r2_tilde_aux, G0.original_auxiliary_transformed_corrospondence[i].f2_aux);   
    }
    else{
      G0.original_auxiliary_transformed_corrospondence[i].f2_aux = 0;
    }
    
    if(G0.original_auxiliary_transformed_corrospondence[i].x3 != 0){
      round_flow( phi,G0.original_auxiliary_transformed_corrospondence[i].battery3_aux,G0.original_auxiliary_transformed_corrospondence[i].r3_tilde_aux,  G0.original_auxiliary_transformed_corrospondence[i].f3_aux);   
    }
    else{
      G0.original_auxiliary_transformed_corrospondence[i].f3_aux = 0; 
    }
    
  }
  
  #ifdef VERBOSE
  cout << "G.f_0 " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++){
    cout << i << " :  " << to_double(G.f_0[i]) << endl;
  }
  #endif
  
  for(unsigned int edge_index_of_G0 =1; edge_index_of_G0<= G0.no_of_edges; edge_index_of_G0++){     
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == true){
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux + G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_delta_wye = -G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux;
    }
    else{
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux + G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_delta_wye = -G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux;
    }
  }
  
  #ifdef VERBOSE
  cout << "G_delta_wye.f_0 " << endl;
  for(unsigned int i=1; i <= G_delta_wye.no_of_edges; i++){
    cout << i << " :  " << to_double(G_delta_wye.f_0[i]) << endl;
  }
  #endif
  
  #ifdef VERBOSE
  cout << "G.f[i] = " << to_double(G.f[i]) << endl;
  #endif
} 

template<typename Network>
void get_resistances( Network& N )
{
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  Graph& G0 = N.G;
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    N.arcdata[a].resistance_lower = RationalType(1)/(N.arcdata[a].xlower * N.arcdata[a].xlower);
    RationalType x_upper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    N.arcdata[a].resistance_upper = RationalType(1)/(x_upper * x_upper);
    N.arcdata[a].resistance_roof = RationalType(1)/(N.arcdata[a].infeasibility * N.arcdata[a].infeasibility);
  }
}


/** \brief Rounds the Delta-Wye transformed resistances for the electrical flow problem
 * 
 * @param G0 The Original Graph
 * @param rho Scaling Factor
 * 
 */
template< typename IntegerType, typename RationalType >
void 
round_r_DS( Graph<IntegerType, RationalType>& G0, 
		 const RationalType& rho 
) 
{
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
    
    #ifdef VERBOSE
    cout << "a = " << a << endl;
    #endif
    if(G0.original_auxiliary_transformed_corrospondence[a].r1_delta_wye == 0){
      G0.original_auxiliary_transformed_corrospondence[a].r1_tilde = 0;
    }
    else{
      round_resistance( rho, G0.original_auxiliary_transformed_corrospondence[a].r1_delta_wye, G0.original_auxiliary_transformed_corrospondence[a].r1_tilde );
    }
    #ifdef VERBOSE
    cout << to_double( G0.original_auxiliary_transformed_corrospondence[a].r1_tilde  ) << endl;
    #endif
    
    if(G0.original_auxiliary_transformed_corrospondence[a].r2_delta_wye == 0){
      G0.original_auxiliary_transformed_corrospondence[a].r2_tilde = 0;
    }
    else{
      round_resistance( rho, G0.original_auxiliary_transformed_corrospondence[a].r2_delta_wye, G0.original_auxiliary_transformed_corrospondence[a].r2_tilde );
    }
    #ifdef VERBOSE
    cout << to_double( G0.original_auxiliary_transformed_corrospondence[a].r2_tilde  ) << endl;
    #endif
    if(G0.original_auxiliary_transformed_corrospondence[a].r3_delta_wye == 0){
      G0.original_auxiliary_transformed_corrospondence[a].r3_tilde = 0;
    }
    else{
      round_resistance( rho, G0.original_auxiliary_transformed_corrospondence[a].r3_delta_wye, G0.original_auxiliary_transformed_corrospondence[a].r3_tilde );
    }
    #ifdef VERBOSE
    cout << to_double( G0.original_auxiliary_transformed_corrospondence[a].r3_tilde  ) << endl << endl;
    #endif
  }
}

/** \brief Gets the node imbalances for the graph
 * 
 * @param G0 The Original Graph
 * @param imbalances_original Imblances on the nodes of the original graph
 * 
 */
template<typename IntegerType, typename RationalType>
void 
get_imbalances_DS(Graph<IntegerType, RationalType>&G0, 
		       vector<IntegerType>& imbalances_original
		      )
{
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
    
    IntegerType alpha3 = G0.original_auxiliary_transformed_corrospondence[a].f3_delta_wye - G0.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde;
    G0.original_auxiliary_transformed_corrospondence[a].imb_vw = alpha3;
    G0.original_auxiliary_transformed_corrospondence[a].imb_u = -alpha3;
    
    IntegerType alpha1 = G0.original_auxiliary_transformed_corrospondence[a].f1_delta_wye - G0.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde;
    G0.original_auxiliary_transformed_corrospondence[a].imb_u -= alpha1;
    imbalances_original[G0.heads[a]] += alpha1;
    
    IntegerType alpha2 = G0.original_auxiliary_transformed_corrospondence[a].f2_delta_wye - G0.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde;
    G0.original_auxiliary_transformed_corrospondence[a].imb_u -= alpha2;
    imbalances_original[G0.tails[a]] += alpha2;
  }
}

/** \brief Warm Starts the Electrical Flow Network
 * 
 * @param G0 The Original Graph
 * @param imbalances_original Imblances of the Original Graph
 * 
 */
template<typename IntegerType, typename RationalType>
void warm_start_the_electrical_flow_network_DS(Graph<IntegerType, RationalType>& G0, 
					       vector<IntegerType>& imbalances_original
					      )
{    
  for(auto a: G0.non_tree_edges){
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid){
      G0.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde = 0;
      G0.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde = 0;
      G0.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde = 0;
      continue;
    }
    
    IntegerType delta = G0.original_auxiliary_transformed_corrospondence[a].imb_vw;
    G0.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde += delta;
    G0.original_auxiliary_transformed_corrospondence[a].imb_u += delta;
    
    delta = G0.original_auxiliary_transformed_corrospondence[a].imb_u;
    G0.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde -= delta;
    imbalances_original[G0.heads[a]] += delta;
  }
  
  for(auto a: boost::adaptors::reverse(G0.tree_edges)){
    IntegerType delta = G0.original_auxiliary_transformed_corrospondence[abs(a)].imb_vw;
    G0.original_auxiliary_transformed_corrospondence[abs(a)].electrical_flow3_tilde += delta;
    G0.original_auxiliary_transformed_corrospondence[abs(a)].imb_u += delta;
    
    if(a < 0){
      delta = imbalances_original[G0.heads[-a]];
      G0.original_auxiliary_transformed_corrospondence[-a].electrical_flow1_tilde += delta;
      G0.original_auxiliary_transformed_corrospondence[-a].imb_u += delta;
      imbalances_original[G0.heads[-a]] -= delta;
      
      delta = G0.original_auxiliary_transformed_corrospondence[-a].imb_u;
      G0.original_auxiliary_transformed_corrospondence[-a].electrical_flow2_tilde -=delta;
      G0.original_auxiliary_transformed_corrospondence[-a].imb_u -= delta;
      imbalances_original[G0.tails[-a]] += delta;
      
    }
    else{      
      delta = imbalances_original[G0.tail(a)];
      G0.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde += delta;
      G0.original_auxiliary_transformed_corrospondence[a].imb_u += delta;
      imbalances_original[G0.tails[a]] -= delta;
      
      delta = G0.original_auxiliary_transformed_corrospondence[a].imb_u;
      G0.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde -=delta;
      G0.original_auxiliary_transformed_corrospondence[a].imb_u -= delta;
      imbalances_original[G0.heads[a]] += delta; 
    }    
  }
}


template<typename Network>
void primal_step_new(Network& N, typename Network::RationalType rho, typename Network::RationalType phi){
  
  cout << "primal step" << endl;
  typedef typename Network::RationalType RationalType;

#ifdef VERBOSE
  cout << "phi = " << phi << endl;
  cout << "rho = " << rho << endl;
#endif
  
  RationalType max = rho;
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;
  
  for(unsigned int i=1; i<=G0.no_of_edges; i++){
    //if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    RationalType x_hat_lower = phi * to_double(N.arcdata[i].initial_current_lower - N.arcdata[i].current_lower);
    RationalType x_lower = N.arcdata[i].xlower;
    if(abs(x_hat_lower/x_lower) > max 
      && x_lower != 0){      
      max = abs(x_hat_lower/x_lower);
      }
      
    RationalType x_hat_upper = - to_double(x_hat_lower);
    RationalType x_upper = N.arcdata[i].capacity - N.arcdata[i].xlower;
    if(abs(x_hat_upper/x_upper) > max 
	&& x_upper != 0){
	
	max = abs(x_hat_upper/x_upper);
	}
	

   RationalType x_hat_roof = phi * to_double(N.arcdata[i].initial_current_roof - N.arcdata[i].current_roof);
   RationalType x_roof = N.arcdata[i].infeasibility;
   if(abs(x_hat_roof/x_roof) > max 
	  && x_roof != 0){
	  
	  max = abs(x_hat_roof/x_roof);
	  }
	  
  }
  
 
  
  const RationalType lambda = RationalType(1) / RationalType(4); 

#ifdef VERBOSE  
  cout << "lambda = " << lambda << endl;
#endif
  
 for(unsigned int i = 1; i <= G0.no_of_edges; i++)
  {
      RationalType x_hat_lower = phi * to_double(N.arcdata[i].initial_current_lower - N.arcdata[i].current_lower);
      RationalType x_hat_roof  = phi * to_double(N.arcdata[i].initial_current_roof - N.arcdata[i].current_roof);

      if(N.arcdata[i].direction == -1)
      {
	//x_hat_roof *= -1;
	//x_hat_lower *= -1;
      }
#ifdef VERBOSE
      cout << "x_hat_lower = " << x_hat_lower << endl;
      cout << "x_hat_roof = " << x_hat_roof << endl;
#endif
      
    if(N.arcdata[i].direction == 1 )
    {
    if(N.arcdata[i].xlower > 0) N.arcdata[i].xlower -= lambda * (x_hat_lower/max);
    if(N.arcdata[i].infeasibility > 0) N.arcdata[i].infeasibility -= lambda * (x_hat_roof/max);
    }
    else if(N.arcdata[i].direction == 0)
    {
     if(N.arcdata[i].xlower > 0) N.arcdata[i].xlower -= lambda * (x_hat_lower/max); 
    }
    else
    { 
      if(N.arcdata[i].xlower > 0) N.arcdata[i].xlower -= lambda * (x_hat_lower/max);
      if(N.arcdata[i].infeasibility > 0) N.arcdata[i].infeasibility -= lambda * (x_hat_roof/max);
    }       
  }
 
//   vector<RationalType> chi_original(G0.no_of_vertices + 1, 0);
//   vector<RationalType> chi_u(G0.no_of_edges + 1, 0);
//   
//   for(unsigned int v = 1; v <= G0.no_of_vertices; v++){
//    chi_original[v] = to_double(N.nodedata[v].demand);
//   }
//   
//    for(unsigned int a = 1; a <= G0.no_of_edges; a++){
//       chi_u[a] = to_double(N.arcdata[a].capacity);
//    }
  
//   //Non Tree Edges in the Auxilliary Graph
//   for(arc a: G0.non_tree_edges){
//     RationalType x_lower = N.arcdata[a].xlower;
//     RationalType x_roof  = N.arcdata[a].infeasibility;
//     RationalType x_hat_lower = phi * to_double(N.arcdata[a].initial_current_lower - N.arcdata[a].current_lower);
// 
// #ifdef VERBOSE
//     cout << "x_hat_lower = " << x_hat_lower << endl;
// #endif
// 
//     RationalType x_hat_roof  = phi * to_double(N.arcdata[a].initial_current_roof - N.arcdata[a].current_roof);
//     
// #ifdef VERBOSE
//     cout << "x_hat_roof = " << x_hat_roof << endl;
// #endif
//     if(x_lower != 0){
//       RationalType x_bar_a = x_lower - lambda*(x_hat_lower/max);
//       if(N.arcdata[a].direction == 1)
//       {
//       chi_original[G0.tail(a)] += x_bar_a;
//       }
//       if(N.arcdata[a].direction == -1)
//       {
//       chi_original[G0.head(a)] += x_bar_a;	
//       }
//   
//       chi_u[a] -= x_bar_a;      
//       N.arcdata[a].xlower= x_bar_a;
//       assert( x_lower > 0);
//     }
// 
//     if(N.arcdata[a].infeasibility != 0){
//       RationalType x_bar_a = x_roof - lambda*(x_hat_roof/max);
// 	{
// 	chi_original[G0.tails[abs(a)]] += x_bar_a;
//  	chi_original[G0.heads[abs(a)]] -= x_bar_a; 
// 	}
//        N.arcdata[a].infeasibility = x_bar_a;
//        assert(N.arcdata[a].infeasibility > 0);
//     }
//   }
//   
//   //Non Tree Edges in the Auxilliary Graph
//   for(arc a: G0.tree_edges){
//     RationalType x_lower = N.arcdata[abs(a)].xlower;
//     RationalType x_hat_lower = phi * to_double(N.arcdata[abs(a)].initial_current_lower - N.arcdata[abs(a)].current_lower);
//     if(x_lower != 0){
//       RationalType x_bar_a = x_lower - lambda * (x_hat_lower/max);
//       if(N.arcdata[abs(a)].direction == 1)
//       {
//       chi_original[G0.tails[abs(a)]] += x_bar_a;
//       }
//       if(N.arcdata[abs(a)].direction == -1)
//       {
//        chi_original[G0.heads[abs(a)]] += x_bar_a;
//       }
//       chi_u[abs(a)] -= x_bar_a;
//       N.arcdata[abs(a)].xlower = x_bar_a;
//      
//       assert( N.arcdata[abs(a)].xlower > 0);
//     }
//     
//   }
//   
//   //Updating the Tree Edges in the Auxilliary Graph
//   for(arc a: G0.non_tree_edges){
//      RationalType x_upper = N.arcdata[a].capacity - N.arcdata[a].xlower;
//      if(x_upper != 0){
//       RationalType x_bar_a = -chi_u[abs(a)];
//       if(N.arcdata[a].direction == 1)
//       {
//       chi_original[G0.heads[abs(a)]] -= x_bar_a;
//       }
//       else
//       {
//        chi_original[G0.tails[abs(a)]] -= x_bar_a;	
//       }
//       chi_u[abs(a)] += x_bar_a;
//       N.arcdata[abs(a)].xlower = N.arcdata[abs(a)].capacity + x_bar_a;  // N.arcdata[a].xupper = -x_bar_a;
//     }
//   }
//   
//   //Updating the Tree Edges in the Auxilliary Graph
//   for(arc a: boost::adaptors::reverse(G0.tree_edges)){
//     RationalType x_upper = N.arcdata[abs(a)].capacity - N.arcdata[abs(a)].xlower;
//     RationalType x_roof = N.arcdata[abs(a)].infeasibility;
//    
//     if(x_upper != 0){
//       RationalType x_bar_a = -chi_u[abs(a)];
//       if(N.arcdata[abs(a)].direction == 1)
//       {
//       chi_original[G0.heads[abs(a)]] -= x_bar_a;
//       }
//       if(N.arcdata[abs(a)].direction == -1)
//       {
//       chi_original[G0.tails[abs(a)]] -= x_bar_a;	    
//       }
//       chi_u[abs(a)] += x_bar_a; 
//       
//       N.arcdata[abs(a)].xlower = N.arcdata[abs(a)].capacity + x_bar_a;  // N.arcdata[a].xupper = -x_bar_a;
//     }   
//     if(x_roof != 0){
//       RationalType x_bar_a = 0;
//       x_bar_a = -chi_original[G0.tail(a)];      
//       chi_original[G0.tail(a)] += x_bar_a;
//       chi_original[G0.head(a)] -= x_bar_a; 
//   
//       { 
// 	if(N.arcdata[abs(a)].direction == 1)
// 	{
// 	N.arcdata[abs(a)].infeasibility = sign(a) * x_bar_a;
// 	}
// 	if(N.arcdata[abs(a)].direction == -1)
// 	{
// 	 N.arcdata[abs(a)].infeasibility = sign(a) * x_bar_a; 
// 	}
// 	
// 	if(N.arcdata[abs(a)].infeasibility < 0) cout << "inf = " << N.arcdata[abs(a)].infeasibility << endl;
// 	assert( N.arcdata[abs(a)].infeasibility >= 0);
//       }
//     }
//   }
//   
}

/** \brief Does a primal step
 * 
 * @param ST The Spanning Tree which was selected
 * @param x_hat_prime
 * @param g_prime
 * @param z_prime
 * @param x_hat
 * @param x
 * 
 */
template<typename IntegerType, typename RationalType>
void primal_step_DS(Graph<IntegerType,RationalType>&G0,
		    RationalType rho 
)
{
  
  #ifdef VERBOSE
  cout << endl << "primal step ()" <<endl;
  #endif
  
  RationalType max = rho;
  
  for(unsigned int i=1; i<=G0.no_of_edges; i++){
    if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    
    if(abs(G0.original_auxiliary_transformed_corrospondence[i].x_hat1/G0.original_auxiliary_transformed_corrospondence[i].x1) > max 
      && G0.original_auxiliary_transformed_corrospondence[i].x1 != 0){
      
      max = abs(G0.original_auxiliary_transformed_corrospondence[i].x_hat1/G0.original_auxiliary_transformed_corrospondence[i].x1);
      }
      
      if(abs(G0.original_auxiliary_transformed_corrospondence[i].x_hat2/G0.original_auxiliary_transformed_corrospondence[i].x2) > max 
	&& G0.original_auxiliary_transformed_corrospondence[i].x2 != 0){
	
	max = abs(G0.original_auxiliary_transformed_corrospondence[i].x_hat2/G0.original_auxiliary_transformed_corrospondence[i].x2);
	}
	
	if(abs(G0.original_auxiliary_transformed_corrospondence[i].x_hat3/G0.original_auxiliary_transformed_corrospondence[i].x3) > max 
	  && G0.original_auxiliary_transformed_corrospondence[i].x3 != 0){
	  
	  max = abs(G0.original_auxiliary_transformed_corrospondence[i].x_hat3/G0.original_auxiliary_transformed_corrospondence[i].x3);
	  }
	  
  }
  
  const RationalType lambda = RationalType(1) / RationalType(4); 
  
  vector<RationalType> chi_original(G0.no_of_verticies + 1, 0);
  for(unsigned int i = 1; i <= G0.no_of_edges; i++){
    if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    
    chi_original[G0.tails[i]] = G0.original_auxiliary_transformed_corrospondence[i].demand1;
    chi_original[G0.heads[i]] = G0.original_auxiliary_transformed_corrospondence[i].demand2;
  }
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
    G0.original_auxiliary_transformed_corrospondence[a].chi_u = G0.original_auxiliary_transformed_corrospondence[a].demand3;  
  }
  
  
  for(arc a: G0.non_tree_edges){
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
    
    if(G0.original_auxiliary_transformed_corrospondence[a].x1 != 0){
      RationalType x_bar_a = G0.original_auxiliary_transformed_corrospondence[a].x1 - lambda*(G0.original_auxiliary_transformed_corrospondence[a].x_hat1/max);
      chi_original[G0.tail(a)] += x_bar_a;
      G0.original_auxiliary_transformed_corrospondence[a].chi_u -= x_bar_a;
      
      G0.original_auxiliary_transformed_corrospondence[a].x1 = x_bar_a;
      assert( G0.original_auxiliary_transformed_corrospondence[a].x1 > 0);
    }
    
    if(G0.original_auxiliary_transformed_corrospondence[a].x3 != 0){
      RationalType x_bar_a = G0.original_auxiliary_transformed_corrospondence[a].x3 - lambda*(G0.original_auxiliary_transformed_corrospondence[a].x_hat3/max);
      if(G0.original_auxiliary_transformed_corrospondence[a].edge_reversed == true){
	chi_original[G0.tail(a)] -= x_bar_a;
	chi_original[G0.head(a)] += x_bar_a;
      }
      else{
	chi_original[G0.tail(a)] += x_bar_a;
	chi_original[G0.head(a)] -= x_bar_a;
      }
      G0.original_auxiliary_transformed_corrospondence[a].x3 = x_bar_a;
      assert( G0.original_auxiliary_transformed_corrospondence[a].x3 > 0);
    }
  }
  for(arc a: G0.tree_edges){
    if(G0.original_auxiliary_transformed_corrospondence[abs(a)].x1 != 0){
      RationalType x_bar_a = G0.original_auxiliary_transformed_corrospondence[abs(a)].x1 - lambda*(G0.original_auxiliary_transformed_corrospondence[abs(a)].x_hat1/max);
      chi_original[G0.tails[abs(a)]] += x_bar_a;
      G0.original_auxiliary_transformed_corrospondence[abs(a)].chi_u -= x_bar_a;
      G0.original_auxiliary_transformed_corrospondence[abs(a)].x1 = x_bar_a;
       assert( G0.original_auxiliary_transformed_corrospondence[abs(a)].x1 > 0);
    }
  }
  
  for(arc a: G0.non_tree_edges){
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
    
    if(G0.original_auxiliary_transformed_corrospondence[a].x2 != 0){
      RationalType x_bar_a = -G0.original_auxiliary_transformed_corrospondence[abs(a)].chi_u;
      chi_original[G0.heads[abs(a)]] -= x_bar_a;
      G0.original_auxiliary_transformed_corrospondence[abs(a)].chi_u += x_bar_a;
      G0.original_auxiliary_transformed_corrospondence[abs(a)].x2 = -x_bar_a;
      assert( G0.original_auxiliary_transformed_corrospondence[abs(a)].x2 > 0);
    }
  }
  
  for(arc a: boost::adaptors::reverse(G0.tree_edges)){
    if(G0.original_auxiliary_transformed_corrospondence[abs(a)].x2 != 0){
      RationalType x_bar_a = -G0.original_auxiliary_transformed_corrospondence[abs(a)].chi_u;
      chi_original[G0.heads[abs(a)]] -= x_bar_a;
      G0.original_auxiliary_transformed_corrospondence[abs(a)].chi_u += x_bar_a;
      G0.original_auxiliary_transformed_corrospondence[abs(a)].x2 = -x_bar_a;
      assert( G0.original_auxiliary_transformed_corrospondence[abs(a)].x2 > 0);
    }
    
    if(G0.original_auxiliary_transformed_corrospondence[abs(a)].x3 != 0){
      RationalType x_bar_a = -chi_original[G0.tail(a)];
      chi_original[G0.tail(a)] += x_bar_a;
      chi_original[G0.head(a)] -= x_bar_a;  
      
      if(G0.original_auxiliary_transformed_corrospondence[abs(a)].edge_reversed == true){
	G0.original_auxiliary_transformed_corrospondence[abs(a)].x3 = -sign(a) * x_bar_a;
	assert( G0.original_auxiliary_transformed_corrospondence[abs(a)].x3 > 0);
      }
      else{
	G0.original_auxiliary_transformed_corrospondence[abs(a)].x3 = sign(a) * x_bar_a;
	assert( G0.original_auxiliary_transformed_corrospondence[abs(a)].x3 > 0);
      }
    }
    
  }
  
}


// template<typename Network>
// void dual_feasibility_check(Network& N)
// { 
//   typedef typename Network::GraphType Graph;
//   
//   Graph& G0 = N.G;
//   
//   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
//   {
//     if(N.arcdata[a].direction == -1)
//     {
//       assert((N.nodedata[G0.heads[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].slower + N.arcdata[a].cost) < 1e-4);
//     }
//     else
//     {
//      assert((N.nodedata[G0.tails[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].slower + N.arcdata[a].cost) < 1e-4); 
//     }
//     
//     
//     if(N.arcdata[a].direction == -1)
//     {
//       assert((N.nodedata[G0.tails[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].supper) < 1e-4);
//     }
//     else
//     {
//      assert((N.nodedata[G0.heads[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].supper) < 1e-4); 
//     }
//     
//     
//     if(N.arcdata[a].infeasibility != 0)
//     {
//     assert((N.nodedata[G0.tails[a]].potential - N.nodedata[G0.heads[a]].potential - N.arcdata[a].sroof + N.arcdata[a].croof) < 1e-4);
//     }
//     
//     
//   }
//   
// }

template<typename Network>
void dual_step_new(Network& N, typename Network::RationalType duality_gap, typename Network::RationalType q, typename Network::RationalType phi){
  
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  const Graph& G0 = N.G;
#ifdef VERBOSE
  cout << "dual step" << endl;
#endif
    RationalType mu = RationalType(0);
  for(unsigned int a=1; a<=G0.no_of_edges; a++){
    mu += N.arcdata[a].slower * N.arcdata[a].xlower;
    mu += N.arcdata[a].supper * (N.arcdata[a].capacity - N.arcdata[a].xlower);
    mu += N.arcdata[a].sroof  * N.arcdata[a].infeasibility;
    //mu+=s_prime[i];
  }
  mu /= q;
  
  
  for(unsigned int i=1; i<=G0.no_of_vertices; i++){
    const RationalType difference = mu * to_double(N.nodedata[i].voltage) * phi;
    N.nodedata[i].potential -= difference;
  }
  for(unsigned int i = 1; i <= G0.no_of_edges; i++){
     const RationalType difference = mu * to_double(N.arcdata[i].voltage_vw) * phi;
     N.arcdata[i].potentialvw -= difference;
  } 
  
  for(unsigned int i=1; i <=G0.no_of_edges; i++){
    //if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].is_invalid) continue;
    
    if(N.arcdata[i].xlower != 0){
      RationalType difference = 0;
     if(N.arcdata[i].direction == 1 || N.arcdata[i].direction == 0)
      {
	difference = N.nodedata[G0.tails[i]].potential - N.arcdata[i].potentialvw;
      }
       if(N.arcdata[i].direction == -1)
       {
 	difference = N.nodedata[G0.heads[i]].potential - N.arcdata[i].potentialvw;
       }
      N.arcdata[i].slower = N.arcdata[i].cost + difference;
    }
    
    if(N.arcdata[i].capacity - N.arcdata[i].xlower != 0){
      
      RationalType difference = 0;
      
       if(N.arcdata[i].direction == 1 || N.arcdata[i].direction == 0)
       {
	difference = N.nodedata[G0.heads[i]].potential - N.arcdata[i].potentialvw;
      }     
       if(N.arcdata[i].direction == -1){
	 difference = N.nodedata[G0.tails[i]].potential - N.arcdata[i].potentialvw;
 	
       }
#ifdef VERBOSE  
      cout << i << ": " << difference << " = " << N.nodedata[G0.heads[i]].potential << " - " << N.arcdata[i].potentialvw << endl;  
#endif
      N.arcdata[i].supper = 0 + difference; 
    }
    
    if(N.arcdata[i].infeasibility != 0){
  //    if(N.arcdata[i].direction == 1)
      {
	RationalType difference = N.nodedata[G0.tails[i]].potential - N.nodedata[G0.heads[i]].potential;
	N.arcdata[i].sroof = N.arcdata[i].croof + difference;
      }
//       else{
// 	N.arcdata[i].sroof = N.arcdata[i].croof - N.nodedata[G0.tails[i]].potential + N.nodedata[G0.heads[i]].potential;
//       }
    }    
    #ifdef VERBOSE
    if(s[i] < 0){
      cout << i << ": " << s[i] << " " << s_prime[i] << " " << z_prime[i] << endl;
    }
    #endif
    
  }
}

/** \brief Does a dual step
 * 
 * @param G0 The original graph G0
 * @param q Parameter
 * @param phi Scaling Factor
 * @param duality_gap The duaity gap of the potential reduction algorithm
 *
 */
template<typename IntegerType, typename RationalType>
void dual_step_DS(Graph<IntegerType,RationalType>& G0, 
		  unsigned int q, 
		  const RationalType& phi, 
		  RationalType& duality_gap
)
{
  
  for(unsigned int i=1; i<= G0.no_of_edges; i++){
    if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    
    node v  = G0.tails[i];
    IntegerType TIV;
    
    {
      arc a = G0.incident_edges[v].front();
      if(G0.heads[abs(a)] == v)
      {
	{
	  TIV = G0.original_auxiliary_transformed_corrospondence[abs(a)].tree_induced_voltage_w;
	}
      }
      else{
	{
	  TIV = G0.original_auxiliary_transformed_corrospondence[abs(a)].tree_induced_voltage_v;
	}
      }
      
      const RationalType difference = (duality_gap/q)*to_double( TIV )*phi;
      G0.original_auxiliary_transformed_corrospondence[i].y1 -= difference;
    }
    
    v  = G0.heads[i];
    {
      arc a = G0.incident_edges[v].front();
      if(G0.heads[abs(a)] == v)
      {
	{
	  TIV = G0.original_auxiliary_transformed_corrospondence[abs(a)].tree_induced_voltage_w;
	}
      }
      else{
	{
	  TIV = G0.original_auxiliary_transformed_corrospondence[abs(a)].tree_induced_voltage_v;
	}
      }     
      const RationalType difference = (duality_gap/q)*to_double( TIV )*phi;
      G0.original_auxiliary_transformed_corrospondence[i].y2 -= difference;
    }
    
    {
      TIV = G0.original_auxiliary_transformed_corrospondence[i].tree_induced_voltage_vw;
      const RationalType difference = (duality_gap/q)*to_double( TIV )*phi;
      G0.original_auxiliary_transformed_corrospondence[i].y3 -= difference;
    }
  } 
  
  
  for(unsigned int edge_index_of_G0=1; edge_index_of_G0 <=G0.no_of_edges; edge_index_of_G0++){
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].is_invalid) continue;
    
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 != 0){
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s1 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].cost1 + 
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y1 -
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y3;
    }
    
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 != 0){
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s2 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].cost2 +
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y2 -
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y3;
    }
    
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 != 0){
      if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == false){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s3 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].cost3 +
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y1 -
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y2;
      }
      else{
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s3 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].cost3 +
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y2 -
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].y1;
      }
    }    
    #ifdef VERBOSE
    if(s[i] < 0){
      cout << i << ": " << s[i] << " " << s_prime[i] << " " << z_prime[i] << endl;
    }
    #endif
    
  }
}

/** \brief Calculates the Duality gap of the Potential Reduction Algorithm
 * 
 * @param G0 The Original Graph
 * 
 * @return Returns the Duality Gap
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType 
calculate_duality_gap_DS( const Graph<IntegerType, RationalType>& G0)
{ 
  RationalType duality_gap = 0;
  
  for(unsigned int i =1; i<=G0.no_of_edges; i++){
    if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    
    const RationalType xasa = G0.original_auxiliary_transformed_corrospondence[i].x1 * G0.original_auxiliary_transformed_corrospondence[i].s1 +
    G0.original_auxiliary_transformed_corrospondence[i].x2 * G0.original_auxiliary_transformed_corrospondence[i].s2 +
    G0.original_auxiliary_transformed_corrospondence[i].x3 * G0.original_auxiliary_transformed_corrospondence[i].s3 ;
    duality_gap += xasa;
    #ifdef VERBOSE
    cout << x[i] << " * " << s[i] << " = " << xasa << endl;
    #endif
  }  
  return duality_gap; 
}


/** \brief Computes the 2-norm of z'
 * 
 * @param G0 The Original Graph
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType 
compute_z_prime_mod_DS(Graph<IntegerType, RationalType>& G0)
{
  #ifdef VERBOSE
  cout << endl << "compute z_prime_mod" << endl;
  #endif
  RationalType z_prime_mod = 0.0;
  
  for( unsigned int a=1; a<=G0.no_of_edges; a++){
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
    
    if(G0.original_auxiliary_transformed_corrospondence[a].x1 > 0){
      z_prime_mod += G0.original_auxiliary_transformed_corrospondence[a].z_prime1 * G0.original_auxiliary_transformed_corrospondence[a].z_prime1;
    }
    if(G0.original_auxiliary_transformed_corrospondence[a].x2 > 0){
      z_prime_mod += G0.original_auxiliary_transformed_corrospondence[a].z_prime2 * G0.original_auxiliary_transformed_corrospondence[a].z_prime2;
    }
    if(G0.original_auxiliary_transformed_corrospondence[a].x3 > 0){
      z_prime_mod += G0.original_auxiliary_transformed_corrospondence[a].z_prime3 * G0.original_auxiliary_transformed_corrospondence[a].z_prime3;
    }
  }  
  
  return z_prime_mod;  
}


template<typename Network>
typename Network::RationalType compute_z_prime_mod(typename Network::RationalType& phi, Network& N)
{
  cout << "compute_z_prime_mod " << endl;
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  
  
  Graph& G0 = N.G;
  RationalType z_prime_mod(0);
  
  cout << "phi = " << phi << endl;
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    RationalType difference = 0;
    
    if(N.arcdata[a].direction == 1 || N.arcdata[a].direction == 0)
    {
    difference = phi * to_double(N.nodedata[G0.tails[a]].voltage - N.arcdata[a].voltage_vw);
    }
    else
    {
    difference = phi * to_double(N.nodedata[G0.heads[a]].voltage - N.arcdata[a].voltage_vw);   
    }
    z_prime_mod += ((N.arcdata[a].battery_lower - difference) * N.arcdata[a].xlower) *
		   ((N.arcdata[a].battery_lower - difference) * N.arcdata[a].xlower);
		   
    if(N.arcdata[a].direction == 1 || N.arcdata[a].direction == 0)
    {
    difference = phi * to_double(N.nodedata[G0.heads[a]].voltage - N.arcdata[a].voltage_vw);
    }
    else
    {
    difference = phi * to_double(N.nodedata[G0.tails[a]].voltage - N.arcdata[a].voltage_vw);  
    }
    RationalType x_upper = N.arcdata[a].capacity - N.arcdata[a].xlower;
            
    z_prime_mod += ((N.arcdata[a].battery_upper - difference) * x_upper) * 
		     ((N.arcdata[a].battery_upper - difference) * x_upper);
        
	
    difference = phi * to_double(N.nodedata[G0.tails[a]].voltage - N.nodedata[G0.heads[a]].voltage);
        
    if(N.arcdata[a].infeasibility != 0)
    {
     z_prime_mod += ((N.arcdata[a].battery_roof - difference) * N.arcdata[a].infeasibility) * 
		     ((N.arcdata[a].battery_roof - difference) * N.arcdata[a].infeasibility);
    }
    
  }
   
   return z_prime_mod;
}

template<typename Network>
typename Network::RationalType compute_z_prime_mod_new(Network& N, typename Network::RationalType phi)
{
  typedef typename Network::RationalType RationalType;
  typedef typename Network::GraphType Graph;
  const Graph& G0 = N.G;
  
  RationalType z_prime_mod = 0;
  for(unsigned int i = 1; i <= G0.no_of_edges; i++){
 
      RationalType difference = 0;
    
      {
	difference = to_double(N.nodedata[G0.tails[i]].voltage - N.arcdata[i].voltage_vw);
      }
  
  
      z_prime_mod += ((N.arcdata[i].battery_lower - phi * difference) * N.arcdata[i].xlower) * 
		     ((N.arcdata[i].battery_lower - phi * difference) * N.arcdata[i].xlower);
      
	
     
		     
      difference = to_double(N.nodedata[G0.heads[i]].voltage - N.arcdata[i].voltage_vw);
           
    
      RationalType x_upper = N.arcdata[i].capacity - N.arcdata[i].xlower;
            
      z_prime_mod += ((N.arcdata[i].battery_upper - phi * difference) * x_upper) * 
		     ((N.arcdata[i].battery_upper - phi * difference) * x_upper);
        
	
      {
	difference = to_double(N.nodedata[G0.tails[i]].voltage - N.nodedata[G0.heads[i]].voltage);
      }
  
      if(N.arcdata[i].infeasibility != 0)
      {
      z_prime_mod += ((N.arcdata[i].battery_roof - phi * difference) * N.arcdata[i].infeasibility) * 
		     ((N.arcdata[i].battery_roof - phi * difference) * N.arcdata[i].infeasibility);
      }	     
  }
  
  return z_prime_mod;
  
 }
/** \brief computes vectors g_prime, x_hat_prime, s_hat_prime, z_prime
 * 
 * @param G The graph on which the min-cost-flow algorithm is being run
 * @param g_prime
 * @param x_hat_prime
 * @param s_hat_prime
 * @param z_prime
 * @param x The primal solution x
 * 
 */
template<typename IntegerType, typename RationalType>
void 
compute_x_hat_s_hat_z_prime_DS(Graph<IntegerType, RationalType>& G0, 
				    const std::vector<RationalType>& x, 
				    const RationalType& phi
)
{
  #ifdef VERBOSE
  cout << endl << "compute_x_hat_s_hat_z_prime" << endl;  
  #endif
  
  #ifdef VERBOSE
  RationalType orthocheck = 0;
  #endif
  
  for(unsigned int edge_index_of_G0=1; edge_index_of_G0<=G0.no_of_edges; edge_index_of_G0++){
    
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].is_invalid){
      continue;
    }
    
    IntegerType D;
    RationalType cur, Diff;
    #ifdef VERBOSE
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 == 0){
      cout << "electrical_flow1 = " << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow1 << endl;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow1 = 0;
    }
    #endif
    
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 != 0){
      D = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].tree_induced_voltage_v - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].tree_induced_voltage_vw;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s_hat1 = phi * to_double(D);
      
      cur = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow1;
      Diff = to_double(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux) - cur;
      #ifdef VERBOSE
      cout << "Diff = " << (Diff) << endl << "f_0 = " << to_double(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux) << endl << "cur = " << cur << endl;
      #endif
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x_hat1 = phi*(Diff);
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].z_prime1 = (G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].battery1_aux - phi*to_double(D))*(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1);
    }
    
    
    #ifdef VERBOSE
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 == 0 ){
      cout << "electrical_flow2 = " << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow2 << endl;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow2 = 0;
    }
    #endif
    
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 != 0){
      D = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].tree_induced_voltage_w - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].tree_induced_voltage_vw;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s_hat2 = phi * to_double(D);
      
      cur = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow2;
      Diff = to_double(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux) - cur;
      #ifdef VERBOSE
      cout << "Diff = " << (Diff) << endl << "f_0 = " << to_double(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux) << endl << "cur = " << cur << endl;
      #endif
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x_hat2 = phi*(Diff);
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].z_prime2 = (G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].battery2_aux - phi*to_double(D))*(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2);
    }
    
    
    
    #ifdef VERBOSE
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 == 0){
      cout << "electrical_flow3 = " << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow3 << endl;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow3 = 0;
    }
    #endif
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 != 0){
      D = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].tree_induced_voltage_w - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].tree_induced_voltage_v;
      
      if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == false){
	D = -D;
      }
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s_hat3 = phi * to_double(D);
      
      cur = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow3;
      Diff = to_double(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux) - cur;
      #ifdef VERBOSE
      cout << "Diff = " << (Diff) << endl << "f_0 = " << to_double(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux) << endl << "cur = " << cur << endl;
      #endif
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x_hat3 = phi*(Diff);
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].z_prime3 = (G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].battery3_aux - phi*to_double(D))*(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3);
    }
    
  }
  
  #ifdef VERBOSE
  cout << "orthocheck = " << orthocheck << endl;
  #endif
  
}


/** \brief Gets the corresponding 'rounded' resistances for the auxilliary graph from the Delta-Wye resistances
 * 
 * @param G0 Original Graph
 * 
 */
template<typename IntegerType, typename RationalType>
void
get_back_transformed_resistances_for_the_auxilliary_graph_DS(Graph<IntegerType, RationalType>& G0){
  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
    
    if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
    
    const IntegerType& r3 = G0.original_auxiliary_transformed_corrospondence[a].r3_tilde;
    IntegerType r_2_times_r_3 = G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
    if(r_2_times_r_3 != 0 && r3 != 0){
      multiply_by_bit_shift(r_2_times_r_3, r3);
    }
    else{
      r_2_times_r_3 = 0;
    }
    
    IntegerType r_1_times_r_3 = G0.original_auxiliary_transformed_corrospondence[a].r1_tilde;
    if(r_1_times_r_3 != 0 && r3 != 0){
      multiply_by_bit_shift(r_1_times_r_3, r3);
    }
    else{
      r_1_times_r_3 = 0;
    }
    
    const IntegerType& r2 = G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
    IntegerType r_1_times_r_2 = G0.original_auxiliary_transformed_corrospondence[a].r1_tilde;
    if(r_1_times_r_2 != 0 && r2 != 0){
      multiply_by_bit_shift(r_1_times_r_2, r2);
    }
    else{
      r_1_times_r_2 = 0;
    }
    
    RationalType numerator = to_double( r_1_times_r_2 + r_2_times_r_3 + r_1_times_r_3);
    
    G0.original_auxiliary_transformed_corrospondence[a].r1_tilde_aux = numerator / to_double(G0.original_auxiliary_transformed_corrospondence[a].r1_tilde);
    
    G0.original_auxiliary_transformed_corrospondence[a].r2_tilde_aux  = numerator/to_double(G0.original_auxiliary_transformed_corrospondence[a].r2_tilde);  
    
    G0.original_auxiliary_transformed_corrospondence[a].r3_tilde_aux = numerator/to_double(G0.original_auxiliary_transformed_corrospondence[a].r3_tilde);
    
  }
}

/** \brief Gets the corresponding currents in the Auxilliary w.r.t the currents of the Delta-Wye graph
 * 
 * @param G0 The Original Graph
 * 
 */
template<typename IntegerType, typename RationalType>
void 
get_currents_in_the_original_network_DS(Graph<IntegerType, RationalType>& G0)
{
  
  for(unsigned int edge_index_of_G0 = 1; edge_index_of_G0 <= G0.no_of_edges; edge_index_of_G0++){
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].is_invalid) continue;
    
    const IntegerType& r_a_1 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_tilde;
    IntegerType flow_1_times_r1 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow1_tilde;
    if(flow_1_times_r1 != 0 && r_a_1 != 0){
      multiply_by_bit_shift(flow_1_times_r1, r_a_1);
    }
    else{
      flow_1_times_r1 = 0;
    }
    
    const IntegerType& r_a_2 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_tilde;
    IntegerType flow_2_times_r2 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow2_tilde;
    if(flow_2_times_r2 != 0 && r_a_2 != 0){
      multiply_by_bit_shift(flow_2_times_r2, r_a_2);
    }
    else{
      flow_2_times_r2 = 0;
    }
    
    const IntegerType& r_a_3 = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_tilde;
    IntegerType flow_3_times_r3 =  G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow3_tilde;
    if(flow_3_times_r3 != 0 && r_a_3 != 0){
      multiply_by_bit_shift(flow_3_times_r3, r_a_3);
    }
    else{
      flow_3_times_r3 = 0;
    }
    
    IntegerType x = (flow_2_times_r2) -   (flow_3_times_r3);
    
    G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow1 = to_double(x)/G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_tilde_aux;
    
    IntegerType y = flow_1_times_r1 - flow_3_times_r3;  
    
    G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow2 = to_double(y)/G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_tilde_aux;
    
    IntegerType z = flow_2_times_r2 - flow_1_times_r1; 
    
    G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow3 = to_double(z)/G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_tilde_aux;
    #ifdef VERBOSE    
    cout << "electrical_flow = " << to_double(z_) << " / " << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_tilde_aux << endl;
    #endif
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == true){
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow3 = -G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].electrical_flow3;
    }
  }
}

/** \brief Removes the Arc whose reduced cost is greater than duality gap
 * 
 * @param G0 The original graph 
 * @param edge_index_of_G0 The edge index of original graph which is being removed
 * @param a The arc number of the triangle which is being removed
 * 
 * @return bool Returns if the arc was removed
 */
template<typename IntegerType, typename RationalType>
bool remove_arc(Graph<IntegerType, RationalType>& G0,  
		unsigned int edge_index_of_G0, 
		unsigned int a
)
{
  const auto& data = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0];
  
  if(a == 1 && data.x2 != 0 && data.x3 != 0){
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == false){
      RationalType delta = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1;
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 = 0; 
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 += delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 -= delta;
	return false;
      }
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 += delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 -= delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 -= delta;
	return false;
      }
      
      return true;
      
    }
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == true){
      RationalType delta = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1;
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 = 0;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 += delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 -= delta;
	return false;
      }
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 -= delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 -= delta;
	return false;
      }
      
      return true;
      
    }  
  }
  
  if(a == 2 && data.x1 != 0 && data.x3 != 0){
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == false){
      RationalType delta = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2;
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 -= delta;
	return false;
      }
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2  = 0;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 -= delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 -= delta;
	return false;
      }
      
      return true;
      
    }
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == true){
      RationalType delta = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2;
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 < 0 ){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 -= delta;
	return false;
      }
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 = 0;
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 += delta;  
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 -= delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 -= delta;
	return false;
      }
      return true;
      
    }
  }
  
  if(a == 3 && data.x1 != 0 && data.x2 !=0){
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == false){
      RationalType delta = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3;
      assert( delta > 0 );
      if( delta >= G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 ) return false;
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 -= delta;
	return false;
      }
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 -= delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 < 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 -= delta; 
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 += delta; 
	return false;
      }
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 = 0;
      
      return true;
      cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 << " " << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 << " "
      << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 << endl;
      
    }
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == true){
      RationalType delta = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3;
      assert( delta > 0 );
      if( delta >= G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 ) return false;
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 -= delta;
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1) <= 1e-5 || G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 < 0 ){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
	return false;
      }
      
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 += delta;
      
      if(abs(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2) <= 1e-5 ||G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 < 0 ){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 += delta;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 -= delta;
	return false;
      }
      G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 = 0;
      
      return true;
      cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 << " " << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 << " "
      << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 << endl;
    }
  }
  
  return true;
}


/** \brief Finds the maximum resistance
 * 
 * @param G0 The original graph
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType find_r_max(Graph<IntegerType, RationalType>& G0)
{
  RationalType r_max(0);
  for(unsigned int i = 1; i <= G0.no_of_edges; i++){
    if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    
    if( G0.original_auxiliary_transformed_corrospondence[i].r1_delta_wye > r_max){
      r_max = G0.original_auxiliary_transformed_corrospondence[i].r1_delta_wye; 
    }
    
    if( G0.original_auxiliary_transformed_corrospondence[i].r2_delta_wye > r_max){
      r_max = G0.original_auxiliary_transformed_corrospondence[i].r2_delta_wye; 
    }
    
    if( G0.original_auxiliary_transformed_corrospondence[i].r3_delta_wye > r_max){
      r_max = G0.original_auxiliary_transformed_corrospondence[i].r3_delta_wye; 
    }      
  }
  
  return r_max;
}

/** \brief Finds the minimum resistance
 * 
 * @param G0 The Original Graph
 *
 * @return r_min Returns the minimum resistance  
 */
template<typename IntegerType, typename RationalType>
RationalType 
find_r_min(Graph<IntegerType, RationalType>& G0)
{
  RationalType r_min = numeric_limits<RationalType>::max();
  for(unsigned int i = 1; i <= G0.no_of_edges; i++){
    if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
    
    if( G0.original_auxiliary_transformed_corrospondence[i].r1_delta_wye < r_min && G0.original_auxiliary_transformed_corrospondence[i].r1_delta_wye != 0){
      r_min = G0.original_auxiliary_transformed_corrospondence[i].r1_delta_wye;
    }
    
    if( G0.original_auxiliary_transformed_corrospondence[i].r2_delta_wye < r_min && G0.original_auxiliary_transformed_corrospondence[i].r2_delta_wye != 0){
      r_min = G0.original_auxiliary_transformed_corrospondence[i].r2_delta_wye;
    }
    
    if( G0.original_auxiliary_transformed_corrospondence[i].r3_delta_wye < r_min && G0.original_auxiliary_transformed_corrospondence[i].r3_delta_wye != 0){
      r_min = G0.original_auxiliary_transformed_corrospondence[i].r3_delta_wye;
    }
    
  }
  return r_min;
}


template<typename Network>
void get_initial_currents(Network& N)
{
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  Graph& G0 = N.G; 
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    N.arcdata[a].current_lower = N.arcdata[a].battery_lower/to_double(N.arcdata[a].resistance_lower);
    N.arcdata[a].initial_current_lower = N.arcdata[a].current_lower;
    RationalType current_upper = N.arcdata[a].battery_upper/to_double(N.arcdata[a].resistance_upper);
    N.arcdata[a].cur_src_vw = N.arcdata[a].current_lower + current_upper;
    N.arcdata[a].current_roof  = N.arcdata[a].battery_roof/to_double(N.arcdata[a].resistance_roof);
    N.arcdata[a].initial_current_roof = N.arcdata[a].current_roof;
  }
}

template<typename Network>
void
round_f_DS(typename Network::RationalType phi, Network &N)
{

  cout << "round_f_DS" << endl;
  cout << "phi = " << phi << endl;
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;
  
  Graph& G0 = N.G; 
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
    
    if(N.arcdata[a].xlower != 0){
    //RationalType r_lower = RationalType(1)/(N.arcdata[a].xlower * N.arcdata[a].xlower);
    round_flow(phi, N.arcdata[a].battery_lower, N.arcdata[a].resistance_lower , N.arcdata[a].current_lower);
    N.arcdata[a].initial_current_lower = N.arcdata[a].current_lower;      
    }
    else{
      assert(false);
    }
    
    if(N.arcdata[a].capacity - N.arcdata[a].xlower != 0){
    RationalType r_upper = RationalType(1)/((N.arcdata[a].capacity - N.arcdata[a].xlower) * (N.arcdata[a].capacity - N.arcdata[a].xlower));
    IntegerType current_upper(0);
    round_flow(phi, N.arcdata[a].battery_upper, N.arcdata[a].resistance_upper, current_upper);
    N.arcdata[a].cur_src_vw = current_upper + N.arcdata[a].initial_current_lower;
    }
    else{
      assert(false);
    }
    
    if(N.arcdata[a].infeasibility != 0){
    RationalType r_roof = RationalType(1)/(N.arcdata[a].infeasibility * N.arcdata[a].infeasibility);
    round_flow(phi, N.arcdata[a].battery_roof, N.arcdata[a].resistance_roof, N.arcdata[a].current_roof);
    N.arcdata[a].initial_current_roof = N.arcdata[a].current_roof;
    }
    else{
      N.arcdata[a].current_roof = 0;
      N.arcdata[a].initial_current_roof = 0;
    }
}
  
  
//   for(unsigned int i=1; i<=G0.no_of_edges; i++){
//    
//     if(G0.original_auxiliary_transformed_corrospondence[i].x1 != 0){
//       round_flow( phi,G0.original_auxiliary_transformed_corrospondence[i].battery1_aux, G0.original_auxiliary_transformed_corrospondence[i].r1_tilde_aux, G0.original_auxiliary_transformed_corrospondence[i].f1_aux );  
//     }
//     else{
//       G0.original_auxiliary_transformed_corrospondence[i].f1_aux = 0;
//     }
//     
//     if(G0.original_auxiliary_transformed_corrospondence[i].x2 != 0){
//       round_flow( phi,G0.original_auxiliary_transformed_corrospondence[i].battery2_aux, G0.original_auxiliary_transformed_corrospondence[i].r2_tilde_aux, G0.original_auxiliary_transformed_corrospondence[i].f2_aux);   
//     }
//     else{
//       G0.original_auxiliary_transformed_corrospondence[i].f2_aux = 0;
//     }
//     
//     if(G0.original_auxiliary_transformed_corrospondence[i].x3 != 0){
//       round_flow( phi,G0.original_auxiliary_transformed_corrospondence[i].battery3_aux,G0.original_auxiliary_transformed_corrospondence[i].r3_tilde_aux,  G0.original_auxiliary_transformed_corrospondence[i].f3_aux);   
//     }
//     else{
//       G0.original_auxiliary_transformed_corrospondence[i].f3_aux = 0; 
//     }
//     
//   }
//   
//   #ifdef VERBOSE
//   cout << "G.f_0 " << endl;
//   for(unsigned int i = 1; i <= G.no_of_edges; i++){
//     cout << i << " :  " << to_double(G.f_0[i]) << endl;
//   }
//   #endif
//   
//   for(unsigned int edge_index_of_G0 =1; edge_index_of_G0<= G0.no_of_edges; edge_index_of_G0++){     
//     if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == true){
//       G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux + G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
//       G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
//       G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_delta_wye = -G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux;
//     }
//     else{
//       G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
//       G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_delta_wye = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux + G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_aux;
//       G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f3_delta_wye = -G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f1_aux - G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].f2_aux;
//     }
//   }
//   
//   #ifdef VERBOSE
//   cout << "G_delta_wye.f_0 " << endl;
//   for(unsigned int i=1; i <= G_delta_wye.no_of_edges; i++){
//     cout << i << " :  " << to_double(G_delta_wye.f_0[i]) << endl;
//   }
//   #endif
//   
//   #ifdef VERBOSE
//   cout << "G.f[i] = " << to_double(G.f[i]) << endl;
//   #endif
}

template<typename Network>
void 
round_r_DS(typename Network::RationalType rho, Network &N) 
{
  cout << "round_r_DS" << endl;
  cout << "rho = " << rho << endl;
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  const Graph& G0 = N.G;
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
  
    
    #ifdef VERBOSE
    cout << "a = " << a << endl;
    #endif
    
    RationalType r_lower = RationalType(1)/(N.arcdata[a].xlower * N.arcdata[a].xlower);
    round_resistance(rho, r_lower, N.arcdata[a].resistance_lower);

    RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    RationalType r_upper = RationalType(1)/(xupper * xupper);  
    round_resistance(rho, r_upper, N.arcdata[a].resistance_upper);
    
    RationalType r_roof = RationalType(1)/(N.arcdata[a].infeasibility * N.arcdata[a].infeasibility);
    round_resistance(rho, r_roof, N.arcdata[a].resistance_roof);
}

}


template<typename Network>
typename Network::RationalType find_r_max(Network& N){
typedef typename Network::GraphType Graph;
typedef typename Network::RationalType RationalType;
const Graph& G0 = N.G;
RationalType r_max(0);

  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
  
    RationalType r_lower = RationalType(1)/(N.arcdata[a].xlower * N.arcdata[a].xlower);
    if(r_lower > r_max){
      r_max = r_lower;
    }
    
    RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    RationalType r_upper = RationalType(1)/(xupper * xupper);
    if(r_upper > r_max){
      r_max = r_upper;
    }
    
    RationalType r_roof = RationalType(1)/(N.arcdata[a].infeasibility * N.arcdata[a].infeasibility);
    if(r_roof > r_max && N.arcdata[a].infeasibility != RationalType(0)){
      r_max = r_roof;
    }
    
  }
  
  return r_max;
}

template<typename Network>
typename Network::RationalType find_r_min(Network& N){
typedef typename Network::GraphType Graph;
typedef typename Network::RationalType RationalType;
const Graph& G0 = N.G;
RationalType r_min = numeric_limits<RationalType>::max();
  for(unsigned int a = 1; a <= G0.no_of_edges; a++){
  
    RationalType r_lower = RationalType(1)/(N.arcdata[a].xlower * N.arcdata[a].xlower);
    if(r_lower < r_min){
      r_min = r_lower;
    }
    
    RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    RationalType r_upper = RationalType(1)/(xupper * xupper);
    if(r_upper < r_min){
      r_min = r_upper;
    }
    
    RationalType r_roof = RationalType(1)/(N.arcdata[a].infeasibility * N.arcdata[a].infeasibility);
    if(r_roof < r_min){
      r_min = r_roof;
    }
    
  }
  
  return r_min;
}

template<typename Network>
typename Network::RationalType get_M(Network& N)
{
    typedef typename Network::GraphType Graph;
    typedef typename Network::RationalType RationalType;
    typedef typename Network::IntegerType IntegerType;
    Graph& G0 = N.G;
    
     RationalType M = 0;
    for(unsigned int i = 0; i < G0.non_tree_edges.size(); i++){
      arc a = G0.non_tree_edges[i];
      IntegerType R_a = G0.resistances_accross_cycles[i];
      IntegerType r_a;
      
      #ifdef VERBOSE
      cout << "a = " << a << endl;
      cout << "R_a = " << to_double(R_a) << endl;
      #endif
      
      if(N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper){
	r_a = N.arcdata[a].resistance_lower;
      }
      else{
	r_a = N.arcdata[a].resistance_upper;
      }
      
      #ifdef VERBOSE
      cout << "r_a = " << to_double(r_a) << endl;
      #endif
      RationalType num =  to_double(R_a);
      RationalType den =  to_double(r_a);;
#ifdef VERBOSE
      cout << "num = " << num << endl;
      cout << "den = " << den << endl;
#endif
      M += num/den; 
      
      IntegerType r_a_roof = N.arcdata[a].resistance_roof;
      IntegerType R_a_roof = R_a - (N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper) + N.arcdata[a].resistance_roof;
      M += to_double(R_a_roof)/to_double(r_a_roof);
    }
#ifdef VERBOSE
    cout << "M midway = " << M << endl;
#endif
    for(auto a: G0.tree_edges)
    {
      //if(N.arcdata[abs(a)].infeasibility == 0) continue;
#ifdef VERBOSE
      cout << "a = " << a << endl;
#endif
      IntegerType R_a = N.arcdata[abs(a)].resistance_roof + N.arcdata[abs(a)].resistance_lower 
			+ N.arcdata[abs(a)].resistance_upper;
#ifdef VERBOSE
      cout << "R_a = " << to_double(R_a) << endl;
#endif
      IntegerType r_a = 0;
//      if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper)
       if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
       N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper)
       {
	if(N.arcdata[abs(a)].resistance_lower > N.arcdata[abs(a)].resistance_upper)
	{
	  r_a = N.arcdata[abs(a)].resistance_lower;
	}
	else
	{
	  r_a = N.arcdata[abs(a)].resistance_upper;
	}
	 
       }
      else
      {
	r_a = N.arcdata[abs(a)].resistance_roof;
      }
      
   
#ifdef VERBOSE
      cout << "r_a = " << to_double(r_a) << endl;
#endif
      M += to_double(R_a)/to_double(r_a); 
    }
    
    cout << "M = " << M << endl;
    return M;
}

template<typename Network>
void warm_start_the_electrical_flow_network(Network& N)
{
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  Graph& G0 = N.G;
  
  vector<IntegerType> imbalances_original(G0.no_of_vertices + 1, IntegerType(0));
  vector<IntegerType> imb_vw(G0.no_of_edges + 1, IntegerType(0));
  
  for(unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
   
   IntegerType f0_lower = N.arcdata[a].initial_current_lower;

   IntegerType alpha_lower = f0_lower - N.arcdata[a].current_lower;
  
   imb_vw[a] += alpha_lower;
   if(N.arcdata[a].direction == 1)
   {
    imbalances_original[G0.tails[a]] -= alpha_lower; 
   }
   if(N.arcdata[a].direction == -1)
   {
    imbalances_original[G0.heads[a]] -= alpha_lower; 
   }
   
   IntegerType f0_upper = N.arcdata[a].cur_src_vw - N.arcdata[a].initial_current_lower;
  
   IntegerType alpha_upper = f0_upper - (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower);

   imb_vw[a] += alpha_upper;
   if(N.arcdata[a].direction == 1)
   {
    imbalances_original[G0.heads[a]] -= alpha_upper; 
   }
   if(N.arcdata[a].direction == -1)
   {
    imbalances_original[G0.tails[a]] -= alpha_upper; 
   }
   
   IntegerType f0_roof = N.arcdata[a].initial_current_roof;
   
   
   IntegerType alpha_roof = f0_roof - N.arcdata[a].current_roof;
   imbalances_original[G0.heads[a]] += alpha_roof;
   imbalances_original[G0.tails[a]] -= alpha_roof;
  }
  
//   cout << "Imbalances nodes: " << endl;
//   for(unsigned int v = 1; v <= G0.no_of_vertices; v++)
//   {
//     cout << v << " = " << imbalances_original[v] << endl;
//   }
//   
//   cout << endl << "Imbalances arcs: " << endl;
//   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
//   {
//     cout << a << " = " << imb_vw[a] << endl;
//   }
  
  
    //Updating the Tree Edges in the Auxilliary Graph
  for(arc a: G0.non_tree_edges){
      
    IntegerType delta = imb_vw[a];
    
    imb_vw[a] += delta;
    imbalances_original[G0.tails[a]] -= delta;
    if(N.arcdata[a].direction == 1)
    {
     N.arcdata[a].current_lower -= delta; 
       
    }
    if(N.arcdata[a].direction == -1)
    {
     N.arcdata[a].current_lower += delta;
    }
  }
  
  //Updating the Tree Edges in the Auxilliary Graph
  for(arc a: boost::adaptors::reverse(G0.tree_edges)){
    
    if(a > 0)
    {
      IntegerType delta = imbalances_original[G0.tails[a]];
      imb_vw[a] += delta;
      imbalances_original[G0.tails[a]] -= delta;
      
      if(N.arcdata[a].direction == 1)
      {
	N.arcdata[a].current_lower -= delta;
      }
      if(N.arcdata[a].direction == -1)
      {
	N.arcdata[a].current_lower += delta;
      }
      
      delta = imb_vw[a];
      imb_vw[a] -= delta;
      imbalances_original[G0.heads[a]] += delta;
    }
    
    if(a < 0)
    {
      IntegerType delta = imbalances_original[G0.heads[-a]];
      imb_vw[-a] += delta;
      imbalances_original[G0.heads[-a]] -= delta;
      
      if(N.arcdata[-a].direction == 1)
      {
	N.arcdata[-a].current_lower += delta;
      }
      if(N.arcdata[-a].direction == -1)
      {
	N.arcdata[-a].current_lower -= delta;
      }
      
    delta = imb_vw[-a];
    imb_vw[-a] -= delta;
    imbalances_original[G0.tails[-a]] += delta;
    }
    
  }
  
//      cout << "Imbalances nodes: " << endl;
//       for(unsigned int v = 1; v <= G0.no_of_vertices; v++)
//       {
// 	cout << v << " = " << imbalances_original[v] << endl;
//       }
//   
//       cout << endl << "Imbalances arcs: " << endl;
//       for(unsigned int a = 1; a <= G0.no_of_edges; a++)
//       {
// 	cout << a << " = " << imb_vw[a] << endl;
//       }
  
  flow_conservation_check(N);
}

// template<typename Network, typename RationalType,typename IntegerType>
// void get_spanning_tree_for_electrical_flow_problem(typename Network::GraphType& G0, 
// 						   SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original,
// 						   Network& N)
// {
//   
// 
//   
//     G0.create_low_stretch_tree_wrt_unrounded_resistances_on_original_graph(ST_original, N);
// 
// #ifdef VERBOSE    
//     cout << "root of the tree = " << ST_original.root << endl;
// #endif
//     assert( ST_original.node_tree_corrospondance[ST_original.root].back().tree_index == 0 );
// 	
//     unsigned int tree_index = 1;    
// 	
//     ST_original.get_tree_incident_edges();
//     ST_original.get_initial_state_of_the_tree_decomposition();
//     ST_original.init( ST_original.root, tree_index );	
//     ST_original.get_non_tree_edge_common_tree_index();
//     
// #ifdef VERBOSE
//     print_tree_edges(N);
//     print_non_tree_edges(N);
//     print_directions(N);
// #endif
//     
//     for(unsigned int i =0; i< G0.non_tree_edges.size(); i++){    
//      G0.resistances_accross_cycles.push_back(ST_original.compute_resistance_accross(G0.non_tree_edges[i], N));
//     
//     #ifdef VERBOSE 
//     cout << i << " Resistances accross  =  " << to_double(ST_original.compute_resistance_accross(G0.non_tree_edges[i], N)) << endl;
//     #endif
//    }
//     
//     
// }
						 


template<typename Network>
void get_batteries(Network& N, 
		   typename Network::RationalType& q, 
		   typename Network::RationalType& duality_gap)
{
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  
  Graph& G0 = N.G;  
      for(unsigned int edge = 1; edge <= G0.no_of_edges; edge++){
      
    if(N.arcdata[edge].xlower != 0)
    {
    RationalType g_prime = q*((N.arcdata[edge].slower * N.arcdata[edge].xlower)/duality_gap) - 1;
    N.arcdata[edge].battery_lower = g_prime/N.arcdata[edge].xlower;
    }
    else
    {
      assert(false);
    }
      
    if(N.arcdata[edge].capacity - N.arcdata[edge].xlower != 0)
    {
    RationalType g_prime = q*((N.arcdata[edge].supper * (N.arcdata[edge].capacity - N.arcdata[edge].xlower) )/duality_gap) - 1;
    N.arcdata[edge].battery_upper = g_prime/(N.arcdata[edge].capacity - N.arcdata[edge].xlower);
    }
    else
    {
      assert(false);
    }
    if(N.arcdata[edge].infeasibility != 0)
    {
    RationalType g_prime = q*((N.arcdata[edge].sroof * N.arcdata[edge].infeasibility )/duality_gap) - 1;
    N.arcdata[edge].battery_roof = g_prime/N.arcdata[edge].infeasibility;
    }
    else
    {
      N.arcdata[edge].battery_roof = 0;
    }
      
    }
}

template<typename Network>
void potential_reduction_algorithm_new(Network &N, typename Network::RationalType q){

  cout << "potential reduction algorithm" << endl;
  cout << "q  = " << q << endl;
  //computing the resistances
  
    typedef typename Network::GraphType Graph1; 
    typedef typename Network::IntegerType IntegerType; 
    typedef typename Network::RationalType RationalType; 
    
    Graph1& G0 = N.G;
    
    for(unsigned int  a = 1; a <= G0.no_of_edges; a++)
    {
      if(N.arcdata[a].direction == -1)
      {
      unsigned int temp = G0.heads[a];
      G0.heads[a] = G0.tails[a];
      G0.tails[a] = temp;
      }
    }
    
    for(unsigned int v = 1; v <= G0.no_of_vertices; v++)
    {
      cout << v << ": ";
      for(auto& a: G0.incident_edges[v])
      {
	cout << a << ", ";
	if(N.arcdata[abs(a)].direction == -1) a = -a;
	
      }
      cout << endl;
    }
   
   RationalType duality_gap = update_duality_gap(N);
    
    unsigned int potential_reduction_iteration = 0;
    
    while(duality_gap >= 1)
    {
   
     potential_reduction_iteration++;
     
#ifdef WithoutRounding
    get_resistances(N);  
#endif
    
    RationalType r_min = find_r_min(N);
    RationalType r_max = find_r_max(N);
#ifdef VERBOSE
    cout << "minimum resistance = " << r_min << endl;
    cout << "maximum resistance = " << r_max << endl;
#endif
    
    RationalType old_potential = calculate_potential_function_new(N, q);
    
    cout << " potential reduction duality gap = " << duality_gap << endl;
    
    duality_gap = update_duality_gap(N);
    
   
  
#ifdef WithoutRounding
    RationalType rho = RationalType(1);
#else
    RationalType epsilon(1);
    RationalType rho = epsilon*r_min;
#endif
    
    
   get_batteries(N, q, duality_gap);
    
#ifdef WithoutRounding
    get_resistances(N); 
#else
    round_r_DS(rho, N);
#endif
     
    SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType> ST_original(G0);
    
    get_spanning_tree_for_electrical_flow_problem(G0, ST_original, N);
    
      
    
#ifdef VERBOSE
    print_resistances(N);
#endif

#ifdef WithoutRounding
    RationalType phi = 1; 
#else    
    RationalType M = get_M(N);
    unsigned int m = (RationalType(3)*G0.no_of_edges);
    RationalType Mmr_max = M*m*r_max;
    RationalType theta = RationalType(1)/RationalType(8); 
    RationalType phi = exp(floor(log((theta*rho)/(sqrt(Mmr_max)))/log(10)) * log(10) );
#endif

#ifdef WithoutRounding
    if(potential_reduction_iteration <= 2)
    {
    get_initial_currents(N);
    }
    else
    {
      cout << "WARM START" << endl; 
      for(unsigned int a = 1; a <= G0.no_of_edges; a++)
      {
       N.arcdata[a].initial_current_lower = N.arcdata[a].battery_lower/to_double(N.arcdata[a].resistance_lower);
       RationalType current_upper = N.arcdata[a].battery_upper/to_double(N.arcdata[a].resistance_upper);
       N.arcdata[a].cur_src_vw = N.arcdata[a].initial_current_lower + current_upper;
       N.arcdata[a].initial_current_roof = N.arcdata[a].battery_roof/to_double(N.arcdata[a].resistance_roof);	
      }
    warm_start_the_electrical_flow_network(N);
    //get_initial_currents(N);
    flow_conservation_check(N);
    }
#else
   round_f_DS(phi, N);
#endif
  
    for(unsigned int i=1; i <= G0.no_of_edges; i++)
    {
      N.arcdata[i].current = 0;
    }
    
#ifdef WithoutRounding
   RationalType delta = RationalType(1)/RationalType(8);
#else
   RationalType delta = Mmr_max/(RationalType(2)*rho);
#endif
      electrical_flow_solver(N, ST_original, delta);
            
      RationalType z_prime_mod = compute_z_prime_mod(phi, N);

      if(z_prime_mod >= RationalType(1)/RationalType(4)){
#ifdef VERBOSE  
	  print_primal_variable(N); 
#endif
	  primal_step_new(N, rho, phi);
	
#ifdef VERBOSE
	  print_primal_variable(N);
#endif
	  primal_feasibility_check(N);
	}
	else 
	{
	dual_step_new(N, duality_gap, q, phi);
	dual_feasibility_check(N);
	}


    RationalType new_potential = calculate_potential_function_new(N, q);

    cout << "reduction in potential_reduction = " << old_potential - new_potential << endl;
    cout << "minimum expected reduction = " << RationalType(1)/RationalType(256) << endl;
    cout << "new potential = " << new_potential << endl;
    RationalType old_duality_gap = duality_gap;
    duality_gap = update_duality_gap(N);
    cout << "duality gap reduction = " << old_duality_gap - duality_gap << endl;
    assert( 256.0*(old_potential - new_potential) >= RationalType(1));
    
    
  }
}



/** \brief Runs the Potential Reduction Algorithm
 * 
 * @param G0 The Original Graph
 * @param G The Auxilliary Graph
 * @param x The Primal Solution 
 * @param y The Dual Solution
 * @param q The Parameter
 * @param minimum_dg The value to which the duality gap should reduce in the potential reduction algorithm
 * 
 */
template<typename IntegerType, typename RationalType>
void 
potential_reduction_algorithm_DS(Graph<IntegerType,RationalType>&G0, 
				      Graph<IntegerType, RationalType> &G, 
				      std::vector<RationalType> &x, 
				      std::vector<RationalType> &y, 
				      std::vector<RationalType> &s, 
				      unsigned int q, 
				      RationalType minimum_dg
)
{
  double time_spent;
  clock_t begin, end;
  begin = clock();
  
  #ifdef VERBOSE
  cout << "delta (potential_reduction_algorithm): " << minimum_dg << endl;
  #endif
  
  RationalType duality_gap = calculate_duality_gap_DS(G0);
  #ifdef VERBOSE
  cout << "potential_reduction_algorithm" << endl;
  cout << endl << "Duality Gap: " << duality_gap << endl;
  #endif
  
  
  SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType> ST_original(G0);
  
  #ifdef VERBOSE
  cout << "G0 tree edges " << endl;
  for(auto a: G0.tree_edges){
    cout << a << " , " ;
  }
  cout << endl;
  #endif
  
  int iteration = 0; 
  int primal_steps = 0;
  int dual_steps = 0;
  int electrical_flow_solver_iterations = 0;
  unsigned int no_of_edges_removed = 0;
  unsigned int no_of_invalid_arcs = 0;
  bool did_dual_step = false;
  bool arc_was_removed = false;
  #ifndef NDEBUG
  bool just_removed_arc = false;
  #endif
  #if !defined(NDEBUG) 
  long double old_potential = calculate_potential_function_DS(G0, no_of_edges_removed, q);
  #endif
  
  RationalType initial_potential = calculate_potential_function_DS(G0, no_of_edges_removed, q);
  
  
  RationalType old_tree_condition_number(0);
  RationalType new_tree_condition_number(0);
  #ifndef NDEBUG
  RationalType avg_tree_condition_number(0);
  #endif
  RationalType sum_tree_condition_number(0);
  
  
  RationalType epsilon(RationalType(1));
  
  bool first_itr = true; 
  
  unsigned int edge_index_of_G0 = 0;
  for(unsigned int i = 1; i<= G.no_of_edges; i+=3){
    edge_index_of_G0++;
    
    if(i % 3 == 1) assert(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 == x[i]);
    #ifdef VERBOSE
    cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 << " " << x[i] << endl;
    #endif
    if(i % 3 == 2) assert(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 == x[i+1]);
    #ifdef VERBOSE
    cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 << " " << x[i+1] << endl;
    #endif
    if(i % 3 == 0) assert(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 == x[i+2]);
    #ifdef VERBOSE
    cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 << " " << x[i+2] << endl;
    #endif
    
    if(i % 3 == 1) assert(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s1 == s[i]);
    #ifdef VERBOSE
    cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s1 << " " << s[i] << endl;
    #endif
    if(i % 3 == 2) assert(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s2 == s[i+1]);
    #ifdef VERBOSE
    cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s2 << " " << s[i+1] << endl;
    #endif
    if(i % 3 == 0) assert(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s3 == s[i+2]);
    #ifdef VERBOSE
    cout << G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s3 << " " << s[i+2] << endl;
    #endif
  }
  
  while(duality_gap >=minimum_dg){
    iteration++;
    
    #ifdef VERBOSE
    cout << "itertion number: " << iteration << endl;
    #endif
    no_of_invalid_arcs = 0;
    for(unsigned int i = 1; i <= G0.no_of_edges; i++){
      if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) no_of_invalid_arcs++;
    }
    q = ( ((G0.no_of_edges - no_of_invalid_arcs)*3) - no_of_edges_removed + no_of_invalid_arcs) + ceil( sqrt((((G0.no_of_edges - no_of_invalid_arcs)*3) - no_of_edges_removed + no_of_invalid_arcs)) ); 
    
    RationalType r_max(0);
    
    RationalType r_min = numeric_limits<RationalType>::max();
    for(unsigned int edge_index_of_G0=1; edge_index_of_G0<=G0.no_of_edges; edge_index_of_G0++){
      if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].is_invalid) continue;
      
      if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 != 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r1_aux =  1/(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 *G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 );
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].g_prime1 = q*(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s1 * G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1 )/duality_gap - 1;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].battery1_aux = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].g_prime1/G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x1;
      }
      
      if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 != 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r2_aux =  1/(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 *G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 );
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].g_prime2 = q*(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s2 * G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2 )/duality_gap - 1;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].battery2_aux = G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].g_prime2/G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x2;
      }
      
      if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 != 0){
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].r3_aux =  1/(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 *G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 );
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].g_prime3 = q*(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].s3 * G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3 )/duality_gap - 1;
	G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].battery3_aux =G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].g_prime3/G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].x3;
      }
    }
    
    get_resistances_for_delta_wye_transformed_graph_DS(G0);
    
    r_max = find_r_max(G0);
    
    r_min = find_r_min(G0);
        
    #ifdef VERBOSE
    for(unsigned int i = 1; i <= G_delta_wye.no_of_edges; i++){
      cout << i << " :  " << G_delta_wye.unrounded_resistances[i] << endl;
    }
    cout << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "r_min = " << r_min << endl;
    cout << "r_max = " << r_max << endl;
    #endif
    
    RationalType rho = epsilon*r_min;
    RationalType theta = RationalType(1)/RationalType(8);
    
    RationalType delta =   RationalType(G0.no_of_edges * 3)*r_max/rho;
    RationalType phi = exp(floor(log((theta*rho)/(G0.no_of_edges * 3 * sqrt(r_max)))/log(10)) * log(10) );
    
    #ifdef VERBOSE
    cout << "rho : " << rho << endl;
    cout << "phi : " << phi << endl;
    cout << "theta: " << theta << endl;
    cout << "rho: " << rho << endl;
    cout << "r_max: " << r_max << endl;
    cout << "G0.no_of_edges: " << G0.no_of_edges << endl;
    #endif
    
    round_r_DS(G0, rho);
    
    get_back_transformed_resistances_for_the_auxilliary_graph_DS(G0);
    
    
    for(unsigned int i = 1; i <= G0.no_of_edges; i++){
      if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
      G0.resistances[i] = G0.original_auxiliary_transformed_corrospondence[i].r1_tilde + G0.original_auxiliary_transformed_corrospondence[i].r2_tilde;
      G0.resistances_aux[i] = G0.original_auxiliary_transformed_corrospondence[i].r1_tilde_aux + G0.original_auxiliary_transformed_corrospondence[i].r2_tilde_aux;
    }
    
    
    if( did_dual_step == false || arc_was_removed){
      if(iteration > 5 ){
	G0.resistances_accross_cycles.clear();
	for(unsigned int i=0; i<G0.non_tree_edges.size(); i++){
	  arc a = G0.non_tree_edges[i];
	  G0.resistances_accross_cycles.push_back(ST_original.compute_resistance_accross(a));
	}
	new_tree_condition_number = ST_original.calculate_tree_condition_number();
	
	#ifdef VERBOSE
	cout << "New Stretch: " << new_tree_condition_number << " " << "old stretch: " << old_tree_condition_number << endl;
	#endif
      }
      if(new_tree_condition_number >= old_tree_condition_number || arc_was_removed){
	
	#ifdef VERBOSE 
	cout << "Stretch not good" << endl;
	#endif
	
	#ifdef VERBOSE
	cout << "did_dual_step == false" << endl;
	#endif
	
	ST_original.clear();
	for(unsigned int i = 1; i <= G0.no_of_verticies; i++){
	  assert(ST_original.node_tree_corrospondance[i].back().tree_index == 0);
	  assert(ST_original.node_tree_corrospondance[i].size() == 1);
	}
	
	#ifdef VERBOSE
	cout << "G0.no_of_edges = " << G0.no_of_edges << endl;
	#endif
	
	for(unsigned int i = 1; i <= G0.no_of_edges; i++){
	  if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
	  
	  G0.resistances[i] = G0.original_auxiliary_transformed_corrospondence[i].r1_tilde + G0.original_auxiliary_transformed_corrospondence[i].r2_tilde;
	  G0.resistances_aux[i] = G0.original_auxiliary_transformed_corrospondence[i].r1_tilde_aux + G0.original_auxiliary_transformed_corrospondence[i].r2_tilde_aux;
	}
		
	G0.create_low_stretch_tree_wrt_unrounded_resistances_on_original_graph(ST_original);
	
	#ifdef VERBOSE 
	#ifndef NDEBUG
	cout << "TCN : " << ST_original.calculate_tree_condition_number() << endl;
	cout << "TCN-original = " << ST_original.calculate_tree_condition_number() << endl;
	#endif
	#endif
		
	assert( ST_original.node_tree_corrospondance[ST_original.root].back().tree_index == 0 );
	
	unsigned int tree_index_original = 1;
	ST_original.get_tree_incident_edges();
	ST_original.get_initial_state_of_the_tree_decomposition();
	ST_original.init( ST_original.root, tree_index_original );
		
	ST_original.get_non_tree_edge_common_tree_index();
	
	G0.resistances_accross_cycles.clear();
	for(unsigned int i=0; i<G0.non_tree_edges.size(); i++){
	  arc a = G0.non_tree_edges[i];
	  if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;
	  G0.resistances_accross_cycles.push_back(ST_original.compute_resistance_accross(a)); 
	}
	
	old_tree_condition_number = ST_original.calculate_tree_condition_number();
	sum_tree_condition_number += old_tree_condition_number;
	#ifndef NDEBUG
	avg_tree_condition_number = sum_tree_condition_number/iteration;
	#endif
	
	#ifdef VERBOSE
	cout << "Old Tree Condition Number: " << old_tree_condition_number << endl;
	#endif
	
      }
      else{
	
	#ifdef VERBOSE
	cout << " stretch good enough,  " ;
	cout << "TCN-original = " << ST_original.calculate_tree_condition_number() << endl;
	#endif
	
	old_tree_condition_number = new_tree_condition_number;
	sum_tree_condition_number += old_tree_condition_number;
	#ifndef NDEBUG
	avg_tree_condition_number = sum_tree_condition_number/iteration;
	#endif
	did_dual_step = false;
	ST_original.clear_d();
	ST_original.update_sum_to_the_root(ST_original.root);
      }      
    } 
    else {
      #ifdef VERBOSE
      cout << "dual step taken" << endl;
      #endif
      sum_tree_condition_number += old_tree_condition_number;
      #ifndef NDEBUG
      avg_tree_condition_number = sum_tree_condition_number/iteration;
      #endif
      did_dual_step = false;
      ST_original.clear_d();
    }
    
    
    
    #ifdef VERBOSE
    cout << "now electrical_flow_problem would be called" << endl;
    #endif
    
    for(unsigned int i = 1; i <= G0.no_of_edges; i++){
      if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
      G0.original_auxiliary_transformed_corrospondence[i].current = 0;
    }
    
    
    RationalType M(0);
    
    no_of_invalid_arcs = 0;
    
    for(unsigned int i = 0; i < G0.non_tree_edges.size(); i++){
      arc a = G0.non_tree_edges[i];
      if(G0.original_auxiliary_transformed_corrospondence[a].is_invalid){
	no_of_invalid_arcs++;
	cout << a << " " << no_of_invalid_arcs << endl;
	
	continue;
      }
      
      IntegerType R_a = G0.resistances_accross_cycles[i];
      IntegerType r_a;
      #ifdef VERBOSE
      cout << "a = " << a << endl;
      cout << "R_a = " << to_double(R_a) << endl;
      #endif
      if(G0.original_auxiliary_transformed_corrospondence[a].r1_tilde < G0.original_auxiliary_transformed_corrospondence[a].r2_tilde){
	r_a = G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
      }
      else{
	r_a = G0.original_auxiliary_transformed_corrospondence[a].r1_tilde;
      }
      
      #ifdef VERBOSE
      cout << "r_a = " << to_double(r_a) << endl;
      #endif
      RationalType num =  to_double(R_a);
      RationalType den =  to_double(r_a);;
      M += num/den; 
    }
    r_max = find_r_max(G0);
    RationalType delta_DS = M* ( 3*(G0.no_of_edges - no_of_invalid_arcs))*r_max/(2*rho);
    
    delta = delta_DS;
    #ifdef VERBOSE
    cout << "no of invalid arcs = " << no_of_invalid_arcs << endl;
    #endif
    phi = exp(floor(log((theta*rho)/(sqrt( 3*(G0.no_of_edges - no_of_invalid_arcs))* sqrt(M) * sqrt(r_max)))/log(10)) * log(10) );
    if(G0.no_of_edges > 600) phi = phi/100;
    
    #ifdef RecordTCN
    cout << M << endl;
    #endif
    
    #ifdef RecordPhi
    cout << phi << endl;
    #endif
    #ifdef Recordr_max
    cout << r_max << endl;
    #endif
    
    #ifdef Recordr_min
    cout << r_min << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "phi = " << phi << endl;
    cout << "r_max = " << r_max << endl;
    cout << "M = " << M << endl;
    #endif
    
    round_f_DS(G0, phi);
    
    if(first_itr){
      for(unsigned int i = 1; i <= G0.no_of_edges; i++){
	if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
	
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde = G0.original_auxiliary_transformed_corrospondence[i].f1_delta_wye;
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde = G0.original_auxiliary_transformed_corrospondence[i].f2_delta_wye;
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow3_tilde = G0.original_auxiliary_transformed_corrospondence[i].f3_delta_wye;
      }
    }
    else{
      vector<IntegerType> imbalances_original(G0.no_of_verticies + 1, 0);
      get_imbalances_DS(G0, imbalances_original);
      #ifdef WarmStart
      warm_start_the_electrical_flow_network_DS(G0, imbalances_original); 
      #else
      for(unsigned int i = 1; i <= G0.no_of_edges; i++)
      {
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde = G0.original_auxiliary_transformed_corrospondence[i].f1_delta_wye;
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde = G0.original_auxiliary_transformed_corrospondence[i].f2_delta_wye;
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow3_tilde = G0.original_auxiliary_transformed_corrospondence[i].f3_delta_wye;
      } 
      #endif 
      
    }
    first_itr = false;
    
    
    int iteration_electrical = electrical_flow_problem_DS(ST_original, G0, delta, old_tree_condition_number );
    
    
    
    electrical_flow_solver_iterations += iteration_electrical;
    
    
    #ifdef RecordElectricalFlowIterations
    cout << iteration_electrical << endl;
    #endif
    
    get_currents_in_the_original_network_DS(G0);
    
    
    #ifdef VERBOSE
    cout << endl << "the electrical flow problem called" << endl;
    #endif
    
    compute_x_hat_s_hat_z_prime_DS(G0, x, phi);
    
    #ifdef VERBOSE
    const RationalType mu = duality_gap/q;
    cout << "mu = " << mu << endl;
    #endif
    
    
    RationalType z_prime_mod = compute_z_prime_mod_DS(G0);
    
    #ifdef VERBOSE
    cout<<"Z_PRIME: " << z_prime_mod << endl;
    #endif
    
    #if !defined(NDEBUG) 
    old_potential = calculate_potential_function_DS(G0, no_of_edges_removed, q);
    #endif
    
    if(4*z_prime_mod >= 1){
      #if !defined(NDEBUG) || defined(CHECKPRIMAL)
      cout << "primal-step" << endl;
      #endif
      
      
      primal_step_DS(G0, rho);
      
      
      #ifdef VERBOSE      
      RationalType new_potential = calculate_potential_function_DS(G0, no_of_edges_removed, q);
      RationalType lambda(RationalType(1)/RationalType(4));
      #endif
      
      #ifdef VERBOSE
      if(x_tilde_zero_norm <= rho){
	cout << " ||x~|| <= rho case " << endl;	
	{
	  cout << old_potential - new_potential << " <-> " << "( " << lambda << " * " <<  gTx << " / " << rho
	  << ")  - ( " << lambda << " * " << lambda << " * " << x_tilde2 << " )  / ( " << 2 << " * ( " << 1 << " - " << lambda << " ) *  " << rho << " * " << rho << " )  = " <<  
	  (lambda * gTx/rho) - (lambda * lambda*x_tilde2) / (2 * (1 -lambda) * rho *rho ) << endl;
	}
	assert(old_potential - new_potential >= (lambda * gTx/rho) - (lambda * lambda*x_tilde2) / (2 * (1 -lambda) * rho *rho ) );	
      }
      else{
	cout << " ||x~|| > rho case " << endl;
	cout << old_potential - new_potential << " >= " << " ( " << lambda << " * " << gTx << "/" << x_tilde_zero_norm <<") - (" << lambda 
	<< " * " << lambda << " * " << x_tilde2 << " ) " << " / ( " << 2 << " * ( " << 1 -lambda << " ) * " << x_tilde_zero_norm << " * " << x_tilde_zero_norm << " ) " << endl;
	assert(old_potential - new_potential >= (lambda * gTx/x_tilde_zero_norm) - (lambda * lambda*x_tilde2) / (2 * (1 -lambda) * x_tilde_zero_norm * x_tilde_zero_norm ) );
      }
      #endif
      primal_steps++;
    } else {
      
      #if !defined(NDEBUG) || defined(CHECKPRIMAL)
      cout << "dual-step" << endl;
      #endif
      
      #ifdef VERBOSE
      cout << " s- before: " ;
      for(auto a : s)
      {
	cout << a << " , " ;
      }
      cout << endl;
      #endif      
      
      dual_step_DS(G0, q, phi, duality_gap);
      
      #ifdef ArcRemoval
      
      #ifdef InvestigatePotentialAtArcRemoval
      RationalType new_potential = calculate_potential_function_DS(G0, no_of_edges_removed, q);
      #endif
      
      for(unsigned int i = 1; i <= G0.no_of_edges; i++){
	if(G0.original_auxiliary_transformed_corrospondence[i].is_invalid || G0.original_auxiliary_transformed_corrospondence[i].cant_be_removed) continue;
	RationalType reduced_cost;
	if(G0.original_auxiliary_transformed_corrospondence[i].edge_reversed == false){
	  reduced_cost = G0.original_auxiliary_transformed_corrospondence[i].cost1 - (G0.original_auxiliary_transformed_corrospondence[i].y1 - 
	  G0.original_auxiliary_transformed_corrospondence[i].y3);
	}
	else{
	  reduced_cost = G0.original_auxiliary_transformed_corrospondence[i].cost1 - (G0.original_auxiliary_transformed_corrospondence[i].y2 - 
	  G0.original_auxiliary_transformed_corrospondence[i].y3);
	}
	
	if(reduced_cost > duality_gap && G0.original_auxiliary_transformed_corrospondence[i].x1 != 0 && !G0.removed_arc[i]){
	  #ifdef VERBOSE
	  cout << "remove_arc = " << i << " , " << 1 << " " << G0.original_auxiliary_transformed_corrospondence[i].x1 << endl;
	  #endif
	  if( remove_arc(G0, i, 1) ) {
	    if(!G0.removed_arc[i]) no_of_edges_removed++;
	    
	    #ifdef InvestigatePotentialAtArcRemoval
	    cout << no_of_edges_removed << " " << initial_potential << " " << new_potential << " " << 1 << endl;
	    #endif
	    arc_was_removed = true;
	    #ifndef NDEBUG
	    just_removed_arc = true;
	    #endif
	    G0.removed_arc[i] = true;
	  } else {
	    #ifdef VERBOSE
	    cout << "no!" << endl;
	    #endif
	  }
	  #ifdef VERBOSE
	  cout << "remove_arc = " << i << " , " << 1 << " " << G0.original_auxiliary_transformed_corrospondence[i].x1 << endl;
	  #endif
	}
	
	if(G0.original_auxiliary_transformed_corrospondence[i].edge_reversed == false){
	  reduced_cost = G0.original_auxiliary_transformed_corrospondence[i].cost2 - (G0.original_auxiliary_transformed_corrospondence[i].y2 - 
	  G0.original_auxiliary_transformed_corrospondence[i].y3);
	}
	else{
	  reduced_cost = G0.original_auxiliary_transformed_corrospondence[i].cost2 - (G0.original_auxiliary_transformed_corrospondence[i].y1 - 
	  G0.original_auxiliary_transformed_corrospondence[i].y3);
	}
	
	if(reduced_cost > duality_gap && G0.original_auxiliary_transformed_corrospondence[i].x2 != 0 && !G0.removed_arc[i]){
	  #ifdef VERBOSE
	  cout << "remove_arc = " << i << " , " << 2 << " " << G0.original_auxiliary_transformed_corrospondence[i].x2 << endl;
	  #endif
	  if( remove_arc(G0, i, 2) ) {
	    
	    if(!G0.removed_arc[i]) no_of_edges_removed++;
	    
	    #ifdef InvestigatePotentialAtArcRemoval
	    cout << no_of_edges_removed << " " << initial_potential << " " << new_potential << " " << 2 << endl;
	    #endif
	    arc_was_removed = true;
	    #ifndef NDEBUG
	    just_removed_arc = true;
	    #endif
	    G0.removed_arc[i] = true;
	  } else {
	    #ifdef VERBOSE
	    cout << " no!" << endl;
	    #endif
	  }
	  #ifdef VERBOSE
	  cout << "remove_arc = " << i << " , " << 2 << " " << G0.original_auxiliary_transformed_corrospondence[i].x2 << endl;
	  #endif
	}
	
	if(G0.original_auxiliary_transformed_corrospondence[i].edge_reversed == false){
	  reduced_cost = G0.original_auxiliary_transformed_corrospondence[i].cost3 - (G0.original_auxiliary_transformed_corrospondence[i].y1 - 
	  G0.original_auxiliary_transformed_corrospondence[i].y2);
	}
	else{
	  reduced_cost = G0.original_auxiliary_transformed_corrospondence[i].cost3 - (G0.original_auxiliary_transformed_corrospondence[i].y2 - 
	  G0.original_auxiliary_transformed_corrospondence[i].y1);
	}
	
	if(reduced_cost > duality_gap && G0.original_auxiliary_transformed_corrospondence[i].x3 != 0 && !G0.removed_arc[i]){
	  #ifdef VERBOSE
	  cout << "remove_arc = " << i << " , " << 3 << " " << G0.original_auxiliary_transformed_corrospondence[i].x3 << endl;
	  #endif
	  if( remove_arc(G0, i, 3) ) {
	    if(!G0.removed_arc[i]) no_of_edges_removed++;
	    
	    #ifdef InvestigatePotentialAtArcRemoval
	    cout << no_of_edges_removed << " " << initial_potential << " " << new_potential << "	" << 3 << endl;
	    #endif
	    
	    arc_was_removed = true;
	    #ifndef NDEBUG
	    just_removed_arc = true;
	    #endif
	    G0.removed_arc[i] = true;
	    
	  } else {
	    #ifdef VERBOSE
	    cout << " no!" << endl;
	    #endif
	  }
	  #ifdef VERBOSE
	  cout << "remove_arc = " << i << " , " << 3 << " " << G0.original_auxiliary_transformed_corrospondence[i].x3 << endl;
	  #endif
	}
	
      }
      #endif
      dual_steps++;
      did_dual_step = true;
    }  
    
    duality_gap = calculate_duality_gap_DS(G0);
    
    #ifdef PontentialReductionGapVariation 
    cout << duality_gap << "	" << iteration << endl;  
    #endif
    
    
    #if !defined(NDEBUG) || defined(VERBOSE) 
    const long double new_potential = calculate_potential_function_DS(G0, no_of_edges_removed, q);
    #endif 
    
    #ifdef VERBOSE
    cout<<"updated duality gap: "<<duality_gap<< endl;
    #endif
    
    #if !defined(NDEBUG)
    sanity_check_DS(G0,  G, x, y, s );
    #endif
    //     cout << "old_potential: " << old_potential << endl;
    //     cout << "new_potential: " << new_potential << endl;
    
    
    //#ifdef VERBOSE
    #if !defined(NDEBUG) || defined(VERBOSE)
    cout << "old_potential: " << old_potential << endl;
    cout << "new_potential: " << new_potential << endl;
    cout << old_potential - new_potential << " <- > " << RationalType(1)/RationalType(256) << endl; 
    #endif
    //#endif 
    
    //     if(256.0*(old_potential - new_potential) < 1.0)
    //     {
      //       cout << "ALPHAS: " << endl;
    //       for(auto a: alphas)
    //       {
      // 	cout << a << " , ";
    //       }
    //       
    //       cout << endl;
    //       
    //       for(unsigned int i = 0 ; i < G.currents.size(); i++)
    //       {
      // 	IntegerType old_current = old_currents1[i];
    // 	IntegerType new_current = G.currents[i];
    // 	cout << to_double(old_current - new_current) << " , ";
    //       }
    //       
    //       cout << endl;
    //     }
    #ifdef VERBOSE
    cout << "no_of_edges_removed = " << no_of_edges_removed << endl;
    #endif
    //  assert(95.0*(old_potential - new_potential) >= 1.0 || did_dual_step);
    assert( 256.0*(old_potential - new_potential) >= 1.0 || just_removed_arc);
    #ifndef NDEBUG
    just_removed_arc = false;
    #endif
    //arc_was_removed = false;
    
    #if !defined(NDEBUG) 
    old_potential = new_potential;
    #endif
    
    
  }
  #ifndef NDEBUG
  long double final_potential = calculate_potential_function_DS(G0,no_of_edges_removed, q);
  #endif
  int original_no_edges = G.no_of_edges/3;
  
  #ifndef NDEBUG
  int original_no_verticies = G.no_of_verticies - original_no_edges;
  #endif
  end = clock();
  time_spent = double(end - begin)/CLOCKS_PER_SEC;
  
  #ifndef NDEBUG
  cout << "Time by Potential Reduction Algorithm = " << time_spent << endl;
  cout << "No of verticies in the original graph = " << original_no_verticies << endl;
  cout << "No. of edges in the original graph = " << original_no_edges << endl;
  cout << "needed " << iteration << " iterations" << endl;
  cout << "number of primal steps: " << primal_steps << endl;
  cout << "number of dual steps: " << dual_steps << endl;
  cout << "avg tree condition number: " << avg_tree_condition_number << " per non-tree edge: " << avg_tree_condition_number/(G.no_of_edges - G.no_of_verticies+1) << endl;
  cout << "initial potential: " << initial_potential << endl;
  cout << "final potential: " << final_potential << endl;
  cout << "electrical_flow_solver_iterations_average: " << electrical_flow_solver_iterations/iteration << endl;
  #endif
  cout << original_no_edges << "	" << initial_potential << " " << iteration << " " << primal_steps << " " << dual_steps << "	" <<
  electrical_flow_solver_iterations/iteration << "	" << time_spent << endl;
}

