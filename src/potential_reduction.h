#include<iostream>
#include<queue>
#include<vector>
#include<algorithm>
#include "graph.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/rational.hpp>
#include <boost/math/common_factor.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "random.h"
#include "electrical_flow_problem.h"


using namespace std;


template<typename IntegerType, typename T>
long double logarithm_base(IntegerType base, T number)
{
  return (log(number)/log(base));
}

/** \brief Does a primal sanity check
 * 
 * Checks if the primal constranits are satisfied
 * 
 * @param G The graph on which sanity check is being done
 * @param x The primal solution
 */ 
template<typename IntegerType, typename RationalType>
void 
primal_sanity_check( const Graph<IntegerType, RationalType> &G, 
		     const std::vector<RationalType> &x 
) 
{
  #ifdef VERBOSE
  cout << "primal_sanity_check " << endl;
  #endif
  RationalType err_primal_constraint = 0.0;
  for(unsigned int i=1; i<=G.no_of_verticies; i++) {
    RationalType flow = 0.0;
    
    for(unsigned int j=0; j<G.incident_edges[i].size(); j++ ){
      int edge = G.incident_edges[i][j];
      if(edge < 0){                 
	flow += x[-edge];
      } 
      else{     
	flow -= x[edge];
      }
    }
    
    const RationalType error = fabs( flow - G.demands[i] );
    err_primal_constraint += error;
    if(  error  > 1e-3 ){
      #ifdef VERBOSE
      cout << "primal error of node " << i << ": " << flow << " " << G.demands[i] << " " << error << endl;
      #endif
      assert( abs(error/max( static_cast<RationalType>( 1.0 ),fabs(G.demands[i]))) <= 1e-3 );
    }
    
  }
  #ifdef VERBOSE
  cout<< "Primal Error: " << err_primal_constraint << endl << endl;  
  #endif
}

/** \brief Does a dual sanity check
 * 
 * Checks if the dual constranits are satisfied
 * 
 * @param G The graph on whicht the sanity check is being done
 * @param y The dual solution
 * @param s The dual slack variables
 * 
 */
template<typename IntegerType, typename RationalType>
void dual_sanity_check( const Graph<IntegerType, RationalType> &G, 
			const std::vector<RationalType> &y, 
			const std::vector<RationalType> &s, 
			const vector<RationalType>& x
) 
{
  #ifdef VERBOSE
  cout << "dual_sanity_check" << endl;
  #endif
  RationalType err_dual_constraint = 0;
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    if(x[i] != 0 ){
      const RationalType expected_dual_soln = G.costs[i] + y[G.tails[i]] - y[G.heads[i]];
      const RationalType error = fabs( s[i] - expected_dual_soln ); 
      err_dual_constraint += error;
      #ifdef VERBOSE
      cout << s[i] << " " << expected_dual_soln << " " << err_dual_constraint << " " << error << endl;
      #endif
      
      if(1000*error > 1){
	cout << " i = " << i << " ,  " << s[i] << " " << expected_dual_soln << " " << err_dual_constraint << " " << error << endl;
      }
      assert( 1000*error < 1 );
    }
  }
  #ifdef VERBOSE
  cout << "Dual Error: " << err_dual_constraint << endl;
  #endif
}


/** \brief checks if dual and primal constraints are satisfied
 * 
 * @see primal_sanity_check( const Graph<IntegerType, RationalType> &G, const std::vector<RationalType> &x )
 * @see dual_sanity_check( const Graph<IntegerType, RationalType> &G, const std::vector<RationalType> &y, const std::vector<RationalType> &s)
 * 
 */
template<typename IntegerType, typename RationalType>
void sanity_check( const Graph<IntegerType, RationalType> &G, 
		   const std::vector<RationalType> &x, 
		   const std::vector<RationalType> &y, 
		   const std::vector<RationalType> &s
)
{
  //The Primal Constraint
  primal_sanity_check( G, x );
  //The Dual Constraint
  dual_sanity_check( G, y, s, x );
}

template< typename T >
long double natural_logarithm( const T& a ) {
  return log( a );
}

template< typename T >
long double natural_logarithm( const boost::rational<T>& r ) {
  return log( to_double( r ) );
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
long double calculate_potential_function( const Graph<IntegerType, RationalType> &G, 
					  const std::vector<RationalType>& x, 
					  const std::vector<RationalType>& s,  
					  unsigned int q
)
{
  #ifdef VERBOSE
  cout << "calculate_potential_function()" << endl;
  #endif
  RationalType x_t_s = 0;
  long double ln_x_s =0;
  int no_of_edges_removed = 0;
  for(unsigned int i = 1; i<=G.no_of_edges; i++) {
    
    if(x[i] != 0) {
      
      const RationalType xasa = x[i]*s[i];
      x_t_s+= xasa;
      
      ln_x_s+= log( xasa );
      
      #ifdef VERBOSE
      cout <<  " log( " << x[i] << " * " << s[i] << ") = " << log(x[i]*s[i]) << endl;
      #endif
    }
    else{
      no_of_edges_removed++;
    }
    
  }
  
  const unsigned int m = G.no_of_edges - no_of_edges_removed;
  
  const long double mlogm = m*log(m);
  const long double potential = q*log(x_t_s) - ln_x_s - mlogm;
  #ifdef VERBOSE
  cout << "x_t_s: " << x_t_s << endl;
  cout<<"Potential: " << q << " * " << natural_logarithm(x_t_s) << " - " << ln_x_s << " - " << mlogm << " = " << potential << endl;
  #endif
  return potential; 
}

/** \brief computes 2-norm of z-prime
 * 
 * @param G The Graph on which the min-cost-flow algorithm is being carried out
 * @param z_prime The vector z-prime
 * 
 * @return the value of 2-norm of z-prime
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType compute_z_prime_mod(Graph<IntegerType, RationalType>& G, 
				 std::vector<RationalType> z_prime
)
{
  #ifdef VERBOSE
  cout << endl << "compute z_prime_mod" << endl;
  #endif
  RationalType z_prime_mod = 0.0;
  
  for( unsigned int i=1; i<=G.no_of_edges; i++){
    z_prime_mod += z_prime[i]*z_prime[i];
  }  
  
  return z_prime_mod;  
}

/** \brief finds the maximum 
 * 
 * @param a 
 * @param x
 * 
 * @return Returns the maximum element in the vector x
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType maximum(RationalType& a, 
		     const std::vector<RationalType>& x
)
{
  RationalType max = a;
  for( auto x_i : x )  {
    if( abs(x_i) > max ){
      max = abs(x_i);
    }
  }
  return max;
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
void primal_step( const SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType>& ST, 
		  const std::vector<RationalType>&x_hat_prime,  
		  const std::vector<RationalType>& g_prime, 
		  const std::vector<RationalType>& z_prime, 
		  std::vector<RationalType>&x_hat, 
		  std::vector<RationalType>&x, 
		  vector<RationalType>&s_hat_prime, 
		  RationalType rho 
)
{
  const Graph<IntegerType, RationalType>& G = ST.G;
  #ifdef VERBOSE
  cout << endl << "primal step ()" <<endl;
  #endif
  
  RationalType xhat2 = RationalType(0);
  #ifdef VERBOSE
  RationalType gTx = 0;
  #endif
  RationalType max = rho;
  
  #ifdef VERBOSE
  RationalType x_tilde_norm_zero(0);
  #endif
  
  #ifdef VERBOSE
  
  cout << "max = " << max << endl;
  cout << "x_hat_primes" << endl;
  for(unsigned int i = 1; i<= G.no_of_edges; i++) {
    cout << i << ": " << x_hat_prime[i] << endl;
  }
  #endif
  
  for(unsigned int i=1; i<=G.no_of_edges; i++)  {
    
    #ifdef VERBOSE
    x_tilde_norm_zero += x_hat_prime[i];
    #endif
    xhat2 += x_hat_prime[i]*x_hat_prime[i];
    #ifdef VERBOSE
    gTx += x_hat_prime[i]*g_prime[i];
    #endif
    if(abs(x_hat_prime[i]) > max && x[i] != 0) {
      
      max = abs(x_hat_prime[i]);
    }
    
  }
  
  #ifdef VERBOSE
  cout << "Zero Norm x_tilde: " << x_tilde_norm_zero << endl << "rho : " << rho << endl;
  cout<<std::endl<<"2-norm of x hat prime: "<<xhat2<<std::endl;
  cout<<std::endl<<"g'^T x' = "<<gTx<<std::endl;
  #endif
  
  const RationalType lambda = RationalType(1) / RationalType(4); 
  #ifdef VERBOSE
  cout << "lambda g'^T v' - lambda^2/(2*(1-lambda))*||v'||_2^2 = " << lambda*gTx/max - lambda*lambda/2.0/(1.0-lambda)*xhat2/max/max << endl;
  #endif
  
  vector<RationalType> chi( G.demands );
  for( arc a : G.non_tree_edges ) {
    
    if(x[a] != 0){
      assert( a > 0 );
      
      const RationalType x_bar_a = x[a] - lambda*(x_hat[a]/max);
      
      #ifdef VERBOSE
      cout << " : " << lambda*((x_hat_prime[a] )/max) << " = " << lambda << " *(( " << x_hat_prime[a] << " )/ " << max << ")" << endl;
      #endif
      
      #ifdef VERBOSE
      if(x_bar_a < 0 ) {
	cout << "a = " << a << " , " << "x_bar_a = " << x_bar_a << endl;
      }
      #endif
      
      const node v = G.tail(a);
      const node w = G.head(a);
      chi[v] += x_bar_a;
      chi[w] -= x_bar_a;
      
      #ifdef VERBOSE
      cout << "x[a] before = " << x[a] << endl;
      #endif
      
      x[a] = x_bar_a;
      
      #ifdef VERBOSE
      cout << "x[a] = " << x[a] << endl;
      #endif
    }
    
  }
  for( arc a : boost::adaptors::reverse( G.tree_edges ) ) {
    
    if(x[abs(a)] != 0) {
      const node v = G.tail(a);
      const node w = G.head(a);
      #if VERBOSE
      cout << a << " = (" << v << "," << w << ") " << ST.depth[v] << " " <<  ST.depth[w] << " " << chi[v] << " " << chi[w] << endl;
      #endif
      assert( ST.depth[v] > ST.depth[w] );
      const RationalType x_bar_a = -chi[v];
      if(sign(a) * x_bar_a <= 0){
	cout << sign(a) << " * "  << x_bar_a << endl;
	
      }
      
      chi[v] += x_bar_a;
      chi[w] -= x_bar_a;
      x[abs(a)] = sign(a)*x_bar_a;
      #ifdef VERBOSE
      cout << "x[a] = " << x[abs(a)] << endl;
      #endif
    }
  }
}


/** \brief Does a Dual Step
 * 
 * @param G The Auxilliary graph
 * @param s_prime
 * @param s
 * @param y
 * @param z_prime
 * @param s_hat_prime
 * @param s_hat
 * @param x
 * @param q Parameter
 * @param phi Scaling Factor
 * @param duality_gap
 */
template<typename IntegerType, typename RationalType>
void dual_step(Graph<IntegerType, RationalType>& G, 
	       std::vector<RationalType>& s_prime, 
	       std::vector<RationalType>& s, 
	       std::vector<RationalType> &y, 
	       std::vector<RationalType>& z_prime, 
	       vector<RationalType>& s_hat_prime, 
	       std::vector<RationalType>&s_hat, 
	       std::vector<RationalType>& x, 
	       unsigned int q, 
	       const RationalType& phi, 
	       RationalType& duality_gap
)
{
  #ifdef VERBOSE
  cout << endl << "dual step" <<endl;
  #endif
  
  #ifdef VERBOSE
  RationalType mu = RationalType(0);
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    mu+=s_prime[i];
  }
  mu /= q;
  #endif
  
  #ifdef VERBOSE
  cout << "mu = " << mu << endl;
  #endif
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    #ifdef VERBOSE
    RationalType check = s_prime[i] - mu*(s_hat_prime[i]);
    #endif
    
    #ifdef VERBOSE
    cout << "check: " <<  check << "  " << s_prime[i] << endl;
    #endif
    
    s[i] -= s_hat[i] * duality_gap/q ;
    
    #ifdef VERBOSE
    cout << s[i] << " " << s_prime[i] << " " << z_prime[i] << endl;
    #endif
    assert( s[i] > 0 );
  }  
  for(unsigned int i=1; i<=G.no_of_verticies; i++){
    const RationalType difference = (duality_gap/q)*to_double( G.tree_induced_voltages[i])*phi;
    #ifdef VERBOSE
    cout << y[i] << " + " << difference << " = " << y[i] << " + " << mu << " * " <<to_double(G.tree_induced_voltages[i]) << endl;
    #endif
    y[i] -= difference;
  } 
}



/** \brief Does a dual step
 * 
 * @param G The Graph on which the min-cost-flow algorithm is being carried out
 * @param s_prime
 * @param s
 * @param y
 * @param z_prime
 * @param z
 * @param q
 * 
 */
template<typename IntegerType, typename RationalType>
void dual_step(Graph<IntegerType, RationalType>& G,
	       Graph<IntegerType, RationalType>& G_delta_wye, 
	       std::vector<RationalType>& s_prime, 
	       std::vector<RationalType>& s, 
	       std::vector<RationalType> &y, 
	       std::vector<RationalType>& z_prime, 
	       vector<RationalType>& s_hat_prime, 
	       std::vector<RationalType>&s_hat, 
	       std::vector<RationalType>& x, 
	       unsigned int q, 
	       const RationalType& phi, 
	       RationalType& duality_gap
)
{
  #ifdef VERBOSE
  cout << endl << "dual step" <<endl;
  #endif
  
  #ifdef VERBOSE
  RationalType mu = RationalType(0);
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    mu+=s_prime[i];
  }
  mu /= q;
  #endif
  
  #ifdef VERBOSE
  cout << "mu = " << mu << endl;
  #endif
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    #ifdef VERBOSE
    RationalType check = s_prime[i] - mu*(s_hat_prime[i]);
    #endif
    
    #ifdef VERBOSE
    cout << "check: " <<  check << "  " << s_prime[i] << endl;
    #endif
    
    s[i] -= s_hat[i] * duality_gap/q ;
    
    #ifdef VERBOSE
    cout << s[i] << " " << s_prime[i] << " " << z_prime[i] << endl;
    #endif
    assert( s[i] > 0 );
  }  
  for(unsigned int i=1; i<=G.no_of_verticies; i++){
    const RationalType difference = (duality_gap/q)*to_double( G_delta_wye.tree_induced_voltages[i])*phi;
    #ifdef VERBOSE
    cout << y[i] << " + " << difference << " = " << y[i] << " + " << mu << " * " <<to_double(G.tree_induced_voltages[i]) << endl;
    #endif
    y[i] -= difference;
  } 
  
}

/** \brief computes x_hat_prime, s_hat_prime, z_prime
 * 
 * @param g_prime
 * @param x_hat_prime
 * @param x_hat
 * @param s_hat_prime
 * @param s_hat
 * @param z_prime
 * @param x
 * @param phi
 * 
 */
template<typename IntegerType, typename RationalType>
void compute_x_hat_s_hat_z_prime(Graph<IntegerType, RationalType> &G , 
				 const std::vector<RationalType>& g_prime, 
				 std::vector<RationalType>& x_hat_prime, 
				 std::vector<RationalType>& x_hat,
				 std::vector<RationalType>& s_hat_prime, 
				 std::vector<RationalType>&s_hat, 
				 std::vector<RationalType>& z_prime, 
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
  
  #if !defined(NDEBUG)
  bool x_tildes_zero = true;
  #endif
  
  for(unsigned int i=1; i<=G.no_of_edges; i++) {
    const IntegerType D = G.voltage_drop( i );
    
    #ifdef VERBOSE
    cout << endl << "D = " << to_double(D) << " " << to_double(G.tree_induced_voltages[G.head(i)]) << " " << to_double(G.tree_induced_voltages[G.tail(i)]) <<  endl;
    #endif
    
    #ifdef VERBOSE
    RationalType orthocheck_s(0) , orthocheck_x(0);
    #endif
    
    s_hat[i] = phi * to_double(D);
    
    #ifdef VERBOSE
    orthocheck_s = to_double(D);
    #endif
    
    s_hat_prime[i] =  phi *to_double( D )*x[i];
    
    #ifdef VERBOSE
    cout << endl << "s_hat_prime[i] =  " << s_hat_prime[i] << endl;  
    #endif
    
    #ifdef VERBOSE
    RationalType g = g_prime[i]/x[i];
    cout << endl << "g = " << g << endl;    
    #endif
    
//     #ifdef VERBOSE
//     RationalType x_hat_val = to_double(G.f[i] - G.currents[i]);
//     cout << endl << "x_hat = " << x_hat_val << endl;
//     #endif
    
    
    #ifdef VERBOSE
    cout << "Diff = " << to_double(G.f_0[i]) << " - " << to_double(G.currents[i]) << endl;
    #endif
    
    IntegerType Diff = G.f_0[i] - G.currents[i];
    
    x_hat_prime[i] = phi*to_double(Diff)/x[i];
    
    x_hat[i] = phi * to_double(Diff);
    
    #ifdef VERBOSE
    orthocheck_x = to_double(G.f_0[i] - G.currents[i]);
    #endif
    
    #if !defined(NDEBUG)
    if(x_hat_prime[i] != 0 ) x_tildes_zero = false;
			#endif
//     #ifdef VERBOSE
//     cout << "f[i] = " << to_double(G.f[i]) << "	" << "G.currents[i] = " << to_double(G.currents[i]) << endl;
//     cout <<  endl << "x_hat_prime[i] = " << x_hat_prime[i] << endl;
//     #endif
    
    #ifdef VERBOSE
    orthocheck += orthocheck_s*orthocheck_x; 
    #endif
    
    z_prime[i] = (G.batteries[i] - phi*to_double(D))*(x[i]);
    
  }
  
  assert (!x_tildes_zero);
  #ifdef VERBOSE
  cout << "orthocheck = " << orthocheck << endl;
  #endif
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
void compute_x_hat_s_hat_z_prime( Graph<IntegerType, RationalType> &G, 
				  Graph<IntegerType, RationalType>&G_delta_wye,
				  const std::vector<RationalType>& g_prime, 
				  std::vector<RationalType>& x_hat_prime,
				  std::vector<RationalType>& x_hat,
				  std::vector<RationalType>& s_hat_prime, 
				  std::vector<RationalType>&s_hat, 
				  std::vector<RationalType>& z_prime, 
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
  
  #if !defined(NDEBUG)
  bool x_tildes_zero = true;
  #endif
  
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    const IntegerType D = G_delta_wye.tree_induced_voltages[G.tail(i)] - G_delta_wye.tree_induced_voltages[G.head(i)];                ;//G.voltage_drop( i );
      
      #ifdef VERBOSE
      cout << endl << "D = " << to_double(D) << " " << to_double(G.tree_induced_voltages[G.head(i)]) << " " << to_double(G.tree_induced_voltages[G.tail(i)]) <<  endl;
      #endif
      
      #ifdef VERBOSE
      RationalType orthocheck_s(0) , orthocheck_x(0);
      #endif
      
      s_hat[i] = phi * to_double(D);
      
      #ifdef VERBOSE
      orthocheck_s = to_double(D);
      #endif
      
      s_hat_prime[i] =  phi *to_double( D )*x[i];
      
      #ifdef VERBOSE
      cout << endl << "s_hat_prime[i] =  " << s_hat_prime[i] << endl;  
      #endif
      
      #ifdef VERBOSE
      RationalType g = g_prime[i]/x[i];
      cout << endl << "g = " << g << endl;    
      #endif
      
      #ifdef VERBOSE
      RationalType x_hat_val = to_double(G.f[i] - G.currents[i]);
      cout << endl << "x_hat = " << x_hat_val << endl;
      #endif
      
      
      #ifdef VERBOSE
      cout << "Diff = " << to_double(G.f_0[i]) << " - " << G.unrounded_currents[i] << endl;
      #endif
      
      RationalType Diff = to_double(G.f_0[i]) - G.unrounded_currents[i];
      
      x_hat_prime[i] = phi*(Diff)/x[i];
      
      x_hat[i] = phi * (Diff);
      
      #ifdef VERBOSE
      orthocheck_x = to_double(G.f_0[i] - G.currents[i]);
      #endif
      
      #if !defined(NDEBUG)
      if(x_hat_prime[i] != 0 ) x_tildes_zero = false;
      #endif
      
      #ifdef VERBOSE
      cout << "f[i] = " << to_double(G.f[i]) << "	" << "G.currents[i] = " << to_double(G.currents[i]) << endl;
      cout <<  endl << "x_hat_prime[i] = " << x_hat_prime[i] << endl;
      #endif
      
      #ifdef VERBOSE
      orthocheck += orthocheck_s*orthocheck_x; 
      #endif
      
      z_prime[i] = (G.batteries[i] - phi*to_double(D))*(x[i]);
      
  }
  
  assert (!x_tildes_zero);
  #ifdef VERBOSE
  cout << "orthocheck = " << orthocheck << endl;
  #endif
  
}


/** \brief Computes minimum value from the vector r
 * 
 * @param G The Graph
 * @param r the vector
 * 
 */
template<typename IntegerType, typename RationalType>
RationalType minimum(Graph<IntegerType, RationalType>&G,
		     vector<RationalType>& r
)
{
  RationalType min = r[1];
  
  for(unsigned int i=1; i<= G.no_of_edges; i++){
    RationalType a = r[i];
    if(a < min){
      min = a;
    }
  }
  
  return min;  
}



/** \brief Calculates the duality gap
 * 
 * @param G The Graph on which the min-cost-flow algorithm is being carried out
 * @param x The primal solution 
 * @param s The dual solution
 * 
 * @return Returns the duality gap   
 * 
 */

template<typename IntegerType, typename RationalType>
RationalType calculate_duality_gap( const Graph<IntegerType,
				    RationalType>& G, 
				    const std::vector<RationalType>& x, 
				    const std::vector<RationalType>& s
)
{ 
  RationalType duality_gap = 0;
  
  for(unsigned int i =1; i<=G.no_of_edges; i++){
    const RationalType xasa = x[i]*s[i];
    duality_gap += xasa;
    #ifdef VERBOSE
    cout << x[i] << " * " << s[i] << " = " << xasa << endl;
    #endif
  }  
  return duality_gap; 
}



/** \brief Rounds the resistances
 * 
 * @param G Auxilliary Graph
 * @param rho Parameter
 * @param resistances the vector of un-rounded resistances
 * 
 */
template< typename IntegerType, typename RationalType >
void round_r(Graph<IntegerType, RationalType>&G, 
	     const RationalType& rho, 
	     const vector<RationalType>& resistances 
) 
{
  for( unsigned a = 1; a <= G.no_of_edges; ++a ) {
    round_resistance( rho, resistances[a], G.resistances[a] );
  }
  #ifdef VERBOSE
  cout << "resistances: " << endl;
  for(unsigned int a =1; a <= G.no_of_edges; a++){
    cout << a << " " << to_double(G.resistances[a]) << endl;
  }
  cout << endl;
  #endif
}


/** \brief Finds the rounded initial flow
 * 
 * @param x The Primal Solution
 * @param G The Graph
 * @param phi The Scaling Factor
 * 
 */

template<typename IntegerType, typename RationalType>
void round_f(vector<RationalType>& x, 
	     Graph<IntegerType, RationalType>&G, 
	     const RationalType& phi 
)
{  
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    if(x[i] != 0){     
      round_flow( phi, G.batteries[i], G.resistances[i], G.f_0[i] );         
    }
  }
} 


/** \brief Finds the rounded initial flow
 * 
 * @param x The Primal Solution
 * @param G0 The Original Graph 
 * @param G The Auxilliary Graph
 * @param G_delta_wye The Delta-Wye Graph
 * @param phi The Scaling Factor
 * 
 */
template<typename IntegerType, typename RationalType>
void round_f(vector<RationalType>& x, 
	     Graph<IntegerType, RationalType>&G0, 
	     Graph<IntegerType, RationalType>&G, 
	     Graph<IntegerType, RationalType>& G_delta_wye, 
	     const RationalType& phi 
)
{
  for(unsigned int i=1; i<=G.no_of_edges; i++){
    if(x[i] != 0){      
      round_flow( phi, G.batteries[i], G.rounded_resistances_auxiallary[i], G.f_0[i] );    
    }
    if(x[i] == 0){
      G.f_0[i] = 0;
    }
    
  }
  
  #ifdef VERBOSE
  cout << "G.f_0 " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++){
    cout << i << " :  " << to_double(G.f_0[i]) << endl;
  }
  #endif
  unsigned int edge_index_of_G0 = 0;
  for(unsigned int i=1; i<= G.no_of_edges; i+=3){
    
    edge_index_of_G0++;
    
    if(G0.original_auxiliary_transformed_corrospondence[edge_index_of_G0].edge_reversed == true){
      G_delta_wye.f_0[i] = G.f_0[i+1] + G.f_0[i+2];
      G_delta_wye.f_0[i+1] = G.f_0[i] - G.f_0[i+2];
      G_delta_wye.f_0[i+2] = -G.f_0[i] - G.f_0[i+1];
    }
    else{
      G_delta_wye.f_0[i] = G.f_0[i+1] - G.f_0[i+2];
      G_delta_wye.f_0[i+1] = G.f_0[i] + G.f_0[i+2];
      G_delta_wye.f_0[i+2] = -G.f_0[i] - G.f_0[i+1];
    }
  }
  
  #ifdef VERBOSE
  cout << "G_delta_wye.f_0 " << endl;
  for(unsigned int i=1; i <= G_delta_wye.no_of_edges; i++)
  {
    cout << i << " :  " << to_double(G_delta_wye.f_0[i]) << endl;
  }
  #endif
  
//   #ifdef VERBOSE
//   cout << "G.f[i] = " << to_double(G.f[i]) << endl;
//   #endif
} 


template<typename IntegerType, typename RationalType>
void round_g(Graph<IntegerType, RationalType>& G, 
	     RationalType phi
)
{
  for(unsigned int i =1 ; i <= G.no_of_edges; i++){
    G.batteries_unrounded[i] =  G.batteries[i];
    G.g_tilde[i] = to_double(phi)*to_double(G.resistances[i]*G.f_0[i]);    
  } 
}

/** \brief Gets the node imbalances in the graph
 * 
 * @param G The Auxilliary Graph G
 * @param imbalances The Imbalances Vector
 * 
 */
template<typename IntegerType, typename RationalType>
void get_imbalances(vector<RationalType>& x, 
		    Graph<IntegerType, RationalType>& G, 
		    vector<IntegerType>& imbalances
)
{
  for(unsigned int a =1; a <= G.no_of_edges; a++){
    if(x[a] != 0){
      const IntegerType alpha = G.f_0[a] - G.currents[a];
      imbalances[G.head(a)] -= alpha;
      imbalances[G.tail(a)] += alpha;
    }
  }
}

/** \brief Balances the flow on the arcs according to the imbalances
 * 
 * @param G The Graph on which the flows have to be balanced
 * @param imbalances The node imbalances of the graph G
 * @param x The flow vector which has to be adjusted
 * 
 */
template<typename IntegerType, typename RationalType, typename T>
void balance_the_arcs(Graph<IntegerType, RationalType>& G, 
		      vector<T>& imbalances, 
		      vector<T>& x
)
{
  for(auto a:  boost::adaptors::reverse(G.tree_edges)) {
    #ifdef VERBOSE
    cout << "a: " << a << endl;
    #endif
    
    if( a < 0 ){
      T delta = imbalances[G.tail(a)];
      const node v = G.tails[-a];
      const node w = G.heads[-a];
      
      #ifdef VERBOSE
      cout << "before: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
      #endif
      //flows[-a] -= delta;
      x[-a] -= delta;
      imbalances[v] += delta;
      imbalances[w] -= delta;
      #ifdef VERBOSE
      cout << "after: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
      #endif 
    } else{
      
      T delta = imbalances[G.tail(a)];
      const node w = G.tails[a];
      const node v = G.heads[a];
      #ifdef VERBOSE
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
      
      cout << "before: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
      #endif 
      
      #ifdef VERBOSE
      cout << "current before = " << to_double(G.currents[a]) << endl;
      #endif
      x[a] += delta;
      
      imbalances[v] += delta;
      imbalances[w] -= delta;
      
      #ifdef VERBOSE
      cout << "after: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
      #endif
      
    }
  }
}


/** \brief Warm Starts the electrical flow problem
 * 
 *  @param G The Graph on which the electrical-flow problem is being run
 *  @param imbalances The node imbalances of the graph
 * 
 */
template<typename IntegerType, typename RationalType>
void warm_start_the_electrical_flow_network(Graph<IntegerType, RationalType>& G, 
					    vector<IntegerType>& imbalances
)
{
  
  balance_the_arcs(G, imbalances, G.currents);
}


/** \brief Gets the corresponding back-transformed 'rounded-resistances' for the Auxilliary Graph from the Delta-Wye graph
 * 
 * @param G The Auxilliary Graph
 * @param G_delta_wye The Delta-Wye Graph
 * 
 */
template<typename IntegerType, typename RationalType>
void get_back_transformed_resistances_for_the_auxilliary_graph(Graph<IntegerType, RationalType>&G, 
							       Graph<IntegerType, RationalType>& G_delta_wye
)
{
  unsigned int index_original_graph = 0;
  
  for(unsigned int i = 1; i <= G.no_of_edges; i+=3){
    
    RationalType x = to_double(G_delta_wye.resistances[i+1]) + to_double(G_delta_wye.resistances[i+2]) + 
    (to_double(G_delta_wye.resistances[i+1]) * to_double(G_delta_wye.resistances[i+2]))/to_double(G_delta_wye.resistances[i]);
    G.rounded_resistances_auxiallary[i] = (x);
    
    
    RationalType y = to_double(G_delta_wye.resistances[i]) + (to_double(G_delta_wye.resistances[i]) * to_double(G_delta_wye.resistances[i+2])/to_double(G_delta_wye.resistances[i+1])) + 
    to_double(G_delta_wye.resistances[i+2]);
    G.rounded_resistances_auxiallary[i+1] = (y);
    
    
    RationalType z = (to_double(G_delta_wye.resistances[i]) * to_double(G_delta_wye.resistances[i+1])/to_double(G_delta_wye.resistances[i+2])) + to_double(G_delta_wye.resistances[i])  + 
    to_double(G_delta_wye.resistances[i+1]);
    G.rounded_resistances_auxiallary[i+2] = (z);	
    
  }
}

/** \brief Gets the Delta-Wye resistances
 * 
 * @param x The primal solution
 * @param G The Auxilliary Graph
 * @param G_delta_wye The Delta-Wye Graph
 * 
 */
template<typename IntegerType, typename RationalType>
void get_resistances_for_delta_wye_transformed_graph(vector<RationalType>& x, 
						     Graph<IntegerType, RationalType>& G, 
						     Graph<IntegerType, RationalType>& G_delta_wye)
{
  for(unsigned int i =1; i <= G_delta_wye.no_of_edges; i+=3){
    if(x[i] != 0 && x[i+1] != 0 && x[i+2] != 0){
      G_delta_wye.unrounded_resistances[i] = (G.unrounded_resistances[i+1]*G.unrounded_resistances[i+2])/(G.unrounded_resistances[i] + G.unrounded_resistances[i+1] + G.unrounded_resistances[i+2]);
      G_delta_wye.unrounded_resistances[i+1] = (G.unrounded_resistances[i]*G.unrounded_resistances[i+2])/(G.unrounded_resistances[i] + G.unrounded_resistances[i+1] + G.unrounded_resistances[i+2]);
      G_delta_wye.unrounded_resistances[i+2] = (G.unrounded_resistances[i+1]*G.unrounded_resistances[i])/(G.unrounded_resistances[i] + G.unrounded_resistances[i+1] + G.unrounded_resistances[i+2]);
    }
    
    if(x[i] == 0 && x[i+1] != 0 && x[i+2] != 0){
      G_delta_wye.unrounded_resistances[i] = 0; 
      G_delta_wye.unrounded_resistances[i+1] = G.unrounded_resistances[i+2];
      G_delta_wye.unrounded_resistances[i+2] = G.unrounded_resistances[i+1]; 
    }
    if(x[i] != 0 && x[i+1] == 0 && x[i+2] != 0){
      G_delta_wye.unrounded_resistances[i] = G.unrounded_resistances[i+2]; 
      G_delta_wye.unrounded_resistances[i+1] = 0; 
      G_delta_wye.unrounded_resistances[i+2] = G.unrounded_resistances[i]; 
    }    
    if(x[i] != 0 && x[i+1] != 0 && x[i+2] == 0){
      G_delta_wye.unrounded_resistances[i] = G.unrounded_resistances[i+1]; 
      G_delta_wye.unrounded_resistances[i+1] = G.unrounded_resistances[i]; 
      G_delta_wye.unrounded_resistances[i+2] = 0; 
    }      
  }
}

/** \brief Gets the currents corresponding the to auxiliary graph from the currents in the Delta-Wye graph
 *  
 * @param G The Auxilliary Graph
 * @param G_delta_wye The Delta-Wye Graph
 * 
 */
template<typename IntegerType, typename RationalType>
void get_currents_in_the_orginal_network(Graph<IntegerType, RationalType>& G_delta_wye, 
					 Graph<IntegerType, 
					 RationalType>& G
)
{
  for(unsigned int i = 1; i <= G.no_of_edges; i+=3){
    
    IntegerType x =  ( (G_delta_wye.currents[i+1]* G_delta_wye.resistances[i+1]) - (G_delta_wye.currents[i+2] * G_delta_wye.resistances[i+2]) );
    
    
    G.unrounded_currents[i] = to_double(x)/(G.rounded_resistances_auxiallary[i]);
    
    IntegerType y = ( (G_delta_wye.currents[i] * G_delta_wye.resistances[i]) - (G_delta_wye.currents[i+2] * G_delta_wye.resistances[i+2]));
    
    G.unrounded_currents[i+1] = to_double(y)/(G.rounded_resistances_auxiallary[i+1]);
    
    IntegerType z = ( (G_delta_wye.currents[i] * G_delta_wye.resistances[i]) - (G_delta_wye.currents[i+1] * G_delta_wye.resistances[i+1]) );
    
    G.unrounded_currents[i+2] = to_double(z)/(G.rounded_resistances_auxiallary[i+2]);
    
  }
  
  #ifdef VERBOSE
  cout << "Unrounded Resistances " << endl;
  for(unsigned int i = 1; i <= G.no_of_edges; i++){
    cout << i << ": " << G.unrounded_resistances[i] << endl;
  }
  cout << endl;
  cout << "rounded resistances of delta-wye" << endl;
  for(unsigned int i=1; i <= G_delta_wye.no_of_edges; i++){
    cout << i << ": " << to_double(G_delta_wye.currents[i]) << " * " << to_double(G_delta_wye.resistances[i]) << endl;
  }
  cout << endl;
  #endif
  
}

/** \brief Runs the potential-reduction-algorithm
 * 
 * @param G The graph on which the min-cost-flow algorithm is being carried out
 * @param x the vector x to store the primal-solution
 * @param y the vector y to store the dual-solution
 * @param s the vector s to store the dual slack variables
 * @param q the parameter q
 * @param minimum_dg
 * 
 */
template<typename IntegerType, typename RationalType>
void potential_reduction_algorithm(Graph<IntegerType,RationalType>&G0, 
				   Graph<IntegerType, RationalType> &G, 
				   Graph<IntegerType, RationalType>&G_delta_wye, 
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
  
  RationalType duality_gap = calculate_duality_gap(G, x,s);
  cout << "potential_reduction_algorithm" << endl;
  cout << endl << "Duality Gap: " << duality_gap << endl;
  
  std::vector<RationalType> g_prime(G.no_of_edges+1,RationalType(0)), s_prime(G.no_of_edges+1,RationalType(0)), s_hat_prime(G.no_of_edges+1,RationalType(0)), x_hat_prime(G.no_of_edges+1,RationalType(0)), z_prime(G.no_of_edges+1,RationalType(0)), x_hat( G.no_of_edges+1, RationalType(0) );
  std::vector<RationalType> z_tilde(G.no_of_edges + 1 , RationalType(0)), x_tilde(G.no_of_edges + 1, RationalType(0)), s_hat(G.no_of_edges + 1, RationalType(0));
  
  
  SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType> ST( G_delta_wye );
  SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType> ST_auxiallary(G); 
  
  G.create_low_stretch_tree_wrt_unrounded_resistances(x, ST_auxiallary);
  
  #ifndef NDEBUG
  long double old_potential = calculate_potential_function(G,x,s,q);
  #endif
  
  int iteration = 0;
  int primal_steps = 0;
  int dual_steps = 0;
  int electrical_flow_solver_iterations = 0;
  
  bool did_dual_step = false;
  bool arc_was_removed = false;
  RationalType initial_potential = calculate_potential_function(G,x,s,q);
  
  RationalType old_tree_condition_number(0);
  RationalType new_tree_condition_number(0);
  RationalType avg_tree_condition_number(0);
  RationalType sum_tree_condition_number(0);
  
  
  RationalType epsilon(RationalType(1));
  RationalType phi_old(0);  
  vector<IntegerType>f_0_old;
  
  bool first_itr = true; 
  int count_arcs_removed = 0;
  while(duality_gap >=minimum_dg){
    iteration++;
    
    #ifdef VERBOSE
    cout << "itertion number: " << iteration << endl;
    #endif
    RationalType g_prime_sum(0);
    
    RationalType r_max(0);
    
    RationalType r_min = numeric_limits<RationalType>::max();
    for(unsigned int i=1; i<=G.no_of_edges; i++) {
      G.unrounded_resistances[i] = 1/(x[i] * x[i]);
      if(x[i] != 0){
	s_prime[i] = s[i]*x[i];
	g_prime[i] = q*s_prime[i]/duality_gap - 1;
	g_prime_sum += g_prime[i];
	
	
	#ifdef VERBOSE
	cout << "resistances updation : " << (G.unrounded_resistances[i]) << "  x[i] = " << x[i] << endl;
	#endif
	
	G.batteries[i] = g_prime[i]/x[i];
	
      }
      
    }
    
    get_resistances_for_delta_wye_transformed_graph(x, G, G_delta_wye);
    
    for(unsigned int i = 1; i <= G_delta_wye.no_of_edges; i++){
      
      if( G_delta_wye.unrounded_resistances[i] > r_max){
	r_max = G_delta_wye.unrounded_resistances[i]; 
      }
      
      if( G_delta_wye.unrounded_resistances[i] < r_min){
	r_min = G_delta_wye.unrounded_resistances[i];
      }
      
    }
    
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
    
    RationalType delta =   RationalType(G_delta_wye.no_of_edges)*r_max/rho;
    RationalType phi = exp(floor(log((theta*rho)/(G_delta_wye.no_of_edges * sqrt(r_max)))/log(10)) * log(10) );
    
    #ifdef VERBOSE
    cout << "rho : " << rho << endl;
    cout << "phi : " << phi << endl;
    #endif
    
    #ifdef VERBOSE
    cout << theta << " *  " << rho << "/ ( " << G.no_of_edges << " *  sqrt( " << r_max << " )) = " << phi << endl;
    #endif
    
    if(first_itr)
    {
      phi_old = phi;
      f_0_old = G.f_0;
    }
    
    if(!first_itr){
      RationalType val(0);
      
      for(unsigned int i = 0; i < G_delta_wye.non_tree_edges.size(); i++){
	arc a = G_delta_wye.non_tree_edges[i];
	IntegerType R_a = G_delta_wye.resistances_accross_cycles[i];
	
	#ifdef VERBOSE
	cout << "R_a = " << to_double(R_a) << endl;
	#endif
	IntegerType r_a = G_delta_wye.resistances[a];
	
	#ifdef VERBOSE
	cout << "r_a = " << to_double(r_a) << endl;
	#endif
	RationalType num =  to_double(R_a);
	RationalType den =  to_double(r_a);;
	val += num/den; 
      }
      RationalType M = val;
      
      
      
      delta = M*G_delta_wye.no_of_edges*r_max/(2*rho);
      #ifdef RecordTCN
      cout << M << endl;
      #endif
      phi = exp(floor(log((theta*rho)/(sqrt(G_delta_wye.no_of_edges)* sqrt(M) * sqrt(r_max)))/log(10)) * log(10) );
      
      
      #ifdef VERBOSE
      cout << "delta = " << delta << endl;
      #endif
    }
    
    round_r(G_delta_wye, rho, G_delta_wye.unrounded_resistances);
    
    #ifdef VERBOSE
    cout << "resistances after rounding " << endl;
    for(unsigned int i = 1; i <= G_delta_wye.no_of_edges; i++){
      cout << i << ": " << to_double(G_delta_wye.resistances[i]) << endl;
    }
    #endif
    
    get_back_transformed_resistances_for_the_auxilliary_graph(G, G_delta_wye);
    
    #ifdef VERBOSE
    cout << "rounded resistances after back-transformation" << endl;
    for(unsigned int i=1; i <= G.no_of_edges; i++){
      cout << i << ": " << G.rounded_resistances_auxiallary[i] << endl;
    }
    #endif
    round_f(x, G0, G, G_delta_wye, phi);
    #ifdef VERBOSE
    cout << "finding G.f0 and G_delta_wye.f0" << endl; 
    for(unsigned int i=1; i <= G.no_of_edges; i++){
      cout << to_double(G.f_0[i]) << "	" << to_double(G_delta_wye.f_0[i]) << endl;
    }
    cout << endl;
    #endif
    
    IntegerType phi_factor;
    get_phi_factor(phi_old, phi,phi_factor );
    
    if(first_itr){
      for(unsigned int i = 1; i <= G_delta_wye.no_of_edges; i++){
	G_delta_wye.currents[i] = G_delta_wye.f_0[i];
      }
      
      for(unsigned int i = 1; i <= G0.no_of_edges; i++){
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde = G_delta_wye.f_0[(i-1)*3 + 1];
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde = G_delta_wye.f_0[(i-1)*3 + 2];
	G0.original_auxiliary_transformed_corrospondence[i].electrical_flow3_tilde = G_delta_wye.f_0[(i-1)*3 + 3];
      }
    }
    else{
      vector<IntegerType> imbalances(G_delta_wye.no_of_verticies+1, 0);
      get_imbalances(x, G_delta_wye, imbalances);
      #ifdef VERBOSE
      cout << "imbalances got " << endl;
      #endif
      warm_start_the_electrical_flow_network( G_delta_wye, imbalances);
      #ifdef VERBOSE
      cout << "warm start done" << endl;
      #endif
    }
    
    first_itr = false;
    
    #ifdef VERBOSE
    for(unsigned int  i = 1; i <= G.no_of_verticies; i++){
      IntegerType old_currents_a(0);
      IntegerType f_0_a(0);
      IntegerType new_flow(0);
      for(auto a : G.incident_edges[i]){
	if(a > 0){
	  f_0_a += G.f_0[abs(a)];
	  new_flow += G.currents[abs(a)];
	}
	else{
	  f_0_a -= G.f_0[abs(a)];
	  new_flow -= G.currents[abs(a)];
	}	
      }
      cout << "  " << to_double(f_0_a) << " " << to_double(new_flow) << endl;
    }
    cout << endl;
    #endif
    phi_old = phi;
    f_0_old = G.f_0;
    
    #ifdef VERBOSE
    cout << "1^T g' = " << g_prime_sum << ", q - m = " << q - G.no_of_edges << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "currents" << endl;
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << to_double(G.currents[i]) << " ; ";
    }
    #endif 
    
    #ifdef VERBOSE
    cout << endl << "create_low_stretch_tree called " << endl;
    cout << "arc_was_removed = " << arc_was_removed << endl;
    #endif
    
    if( did_dual_step == false || arc_was_removed){
      if(iteration > 5 ){
	G_delta_wye.resistances_accross_cycles.clear();
	G_delta_wye.unrounded_resistances_accross_cycles.clear();
	for(unsigned int i=0; i<G_delta_wye.non_tree_edges.size(); i++){
	  arc a = G_delta_wye.non_tree_edges[i];
	  G_delta_wye.unrounded_resistances_accross_cycles.push_back(ST.compute_unrounded_resistance_accross(a));
	}
	new_tree_condition_number = ST.calculate_tree_condition_number();
	
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
	
	ST.clear(); 
	for(unsigned int i = 1; i <= G.no_of_verticies; i++){
	  assert(ST.node_tree_corrospondance[i].back().tree_index == 0);
	  assert(ST.node_tree_corrospondance[i].size() == 1);
	}
	
	G_delta_wye.create_low_stretch_tree_wrt_unrounded_resistances(x, ST ); 
	
	#ifdef VERBOSE
	cout << "tree edges " << endl;
	for(auto a: G_delta_wye.tree_edges){
	  cout << a << " , " ;
	}
	cout << endl << endl;
	
	cout << "non tree edges " << endl;
	for(auto a: G_delta_wye.non_tree_edges){
	  cout << a << " , " ;
	}
	cout << endl << endl;
	#endif
	
	#ifndef NDEBUG
	cout << "TCN : " << ST.calculate_tree_condition_number() << endl;
	#endif
	
	#ifdef VERBOSE
	cout << "LCA vector " << endl << endl;
	for(unsigned int i = 1; i <= G.no_of_edges; i++){
	  cout << i << " : " << ST.LCA[i] << endl;
	}
	#endif
	
	
	assert( ST.node_tree_corrospondance[ST.root].back().tree_index == 0 );
	
	unsigned int tree_index = 1;    
	
	ST.get_tree_incident_edges();
	ST.get_initial_state_of_the_tree_decomposition();
	ST.init( ST.root, tree_index );
	
	ST.get_non_tree_edge_common_tree_index();
	
	G_delta_wye.resistances_accross_cycles.clear();
	for(unsigned int i=0; i<G.non_tree_edges.size(); i++){
	  arc a = G.non_tree_edges[i];
	  G_delta_wye.resistances_accross_cycles.push_back(ST.compute_resistance_accross(a)); //sum_over_cycle( ST, G.non_tree_edges[i], G.resistances, 1 );
	}
	
	old_tree_condition_number = ST.calculate_tree_condition_number();
	sum_tree_condition_number += old_tree_condition_number;
	avg_tree_condition_number = sum_tree_condition_number/iteration;
	
	#ifdef VERBOSE
	cout << "Old Tree Condition Number: " << old_tree_condition_number << endl;
	#endif
	
      }
      else{
	
	#ifdef VERBOSE
	cout << " stretch good enough,  " ;
	cout << "TCN: " << ST.calculate_tree_condition_number() << endl; 
	#endif
	
	old_tree_condition_number = new_tree_condition_number;
	sum_tree_condition_number += old_tree_condition_number;
	avg_tree_condition_number = sum_tree_condition_number/iteration;
	did_dual_step = false;
	ST.clear_d();
	ST.update_sum_to_the_root(ST.root);
      }
      
    } else {
      cout << "dual step taken" << endl;
      sum_tree_condition_number += old_tree_condition_number;
      avg_tree_condition_number = sum_tree_condition_number/iteration;
      did_dual_step = false;
      ST.clear_d();
    }
    
    #ifdef VERBOSE
    cout << "now electrical_flow_problem would be called" << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "Currents before calling the electrical_flow_problem" << endl << endl;
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << i << " ; " << to_double(G.currents[i]) << endl; 
    }
    #endif
    
    
    #ifdef VERBOSE
    cout << "resistances: " << endl; 
    for(unsigned int i = 1; i <= G.no_of_edges; i++){
      cout << i << " : " << to_double(G.resistances[i]) << endl;
    }
    cout << endl << endl;
    #endif
    int iteration_electrical = electrical_flow_problem(x, G0, ST, delta, old_tree_condition_number );
    
    #ifdef VERBOSE
    cout << "currents after the electrical flow solver ... " << endl;
    for(unsigned int i = 1; i <= G_delta_wye.no_of_edges; i++){
      cout << i << ": " << to_double(G_delta_wye.currents[i]) << endl;
    }
    #endif
    electrical_flow_solver_iterations += iteration_electrical;
    
    get_currents_in_the_orginal_network(G_delta_wye, G);
    
    #ifdef VERBOSE
    cout <<  "gap after = " << compute_gap_by_non_tree_currents(ST, ST_auxiallary) << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "currents in the auxiallary network after back-transformation" << endl;
    for(unsigned int i = 1; i <= G.no_of_edges; i++){
      cout << i << " : " << G.unrounded_currents[i] << endl;
    }
    #endif
    
    ofstream file_obj;
    std::string path = "electrical_flow_log.txt";
    file_obj.open(path, std::fstream::app);
    file_obj << G.no_of_edges/3 << " " << duality_gap << " " << iteration_electrical << endl;  
    
    #ifdef VERBOSE
    cout << endl << "the electrical flow problem called" << endl;
    #endif
    
    compute_x_hat_s_hat_z_prime( G, G_delta_wye, g_prime, x_hat_prime, x_hat, s_hat_prime, s_hat, z_prime, x, phi);
    
    
    #ifdef VERBOSE
    cout << endl << "X hat prime" << endl;
    
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << x_hat_prime[i] << ";";
    }
    #endif 
    
    #ifdef VERBOSE
    cout << endl << "Z prime " << endl;
    #endif
    
    #ifdef VERBOSE
    const RationalType mu = duality_gap/q;
    cout << "mu = " << mu << endl;
    #endif
    
    #ifdef VERBOSE
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << z_prime[i] << " =  " << g_prime[i] << " - " << s_hat_prime[i] << ", " << s_prime[i] << ": " << s_prime[i] - mu*s_hat_prime[i] << " <-> " << mu*( 1 + z_prime[i] )<< endl;
      assert( RationalType(1e3) * fabs( s_prime[i] - mu*s_hat_prime[i] - mu*( 1 + z_prime[i] ) ) < RationalType(1) );
    }
    #endif
    
    RationalType z_prime_mod = compute_z_prime_mod(G, z_prime);
    
    #ifdef VERBOSE
    cout<<"Z_PRIME: " << z_prime_mod << endl;
    #endif
    
    if(4*z_prime_mod >= 1){
      #ifndef NDEBUG
      cout << "primal-step" << endl;
      #endif
      
      #ifdef VERBOSE
      RationalType old_potential = calculate_potential_function(G,x,s,q);
      #endif
      
      #ifdef VERBOSE
      RationalType x_tilde_zero_norm(0);
      RationalType gTx(0);
      RationalType x_tilde2(0);
      for(unsigned int i = 1; i <= G.no_of_edges; i++){
	x_tilde_zero_norm += x_hat_prime[i];
	gTx += G.batteries[i]*x_hat_prime[i]*x[i];
	x_tilde2 += x_hat_prime[i]*x_hat_prime[i];
      }
      #endif
      
      primal_step( ST_auxiallary, x_hat_prime, g_prime, z_prime, x_hat, x, s_hat_prime, rho);
      
      #ifdef VERBOSE      
      RationalType new_potential = calculate_potential_function(G,x,s,q);
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
      
      #ifndef NDEBUG
      cout << "dual-step" << endl;
      #endif
      
      #ifdef VERBOSE
      cout << " s- before: " ;
      for(auto a : s){
	cout << a << " , " ;
      }
      cout << endl;
      #endif      
      
      dual_step(G, G_delta_wye, s_prime, s, y, z_prime, s_hat_prime, s_hat, x, q, phi, duality_gap);
      
      
      for(unsigned int i = 1; i <= G.no_of_edges; i++){
	if(G.costs[i] - (y[G.tail(i)] - y[G.head(i)]) > duality_gap){
	  
	  #ifdef VERBOSE
	  cout << G.costs[i] << " - " << (y[G.tail(i)] - y[G.head(i)]) << " = " << duality_gap << endl;
	  cout << "removed " << i << " " << x[i] << endl;
	  #endif
	}
      }
      
      #ifdef VERBOSE
      cout << "s-after: " ;
      
      for(auto a : s){
	cout << a << " , " ;
      }
      #endif
      
      dual_steps++;
      did_dual_step = true;
    }  
    
    duality_gap = update_duality_gap(G,x,s);
    
    #ifdef PontentialReductionGapVariation
    cout << duality_gap << "	" << iteration << endl;
    #endif
    
    
    #if !defined(NDEBUG) || defined(VERBOSE)
    const long double new_potential = calculate_potential_function(G,x,s,q);
    #endif 
    
    #ifdef VERBOSE
    cout<<"updated duality gap: "<<duality_gap<< endl;
    #endif
    
    #if !defined(NDEBUG)
    sanity_check( G, x, y, s );
    #endif
    //     cout << "old_potential: " << old_potential << endl;
    //     cout << "new_potential: " << new_potential << endl;
    
    
    //#ifdef VERBOSE
    #if !defined(NDEBUG) || defined(VERBOSE)
    cout << "old_potential: " << old_potential << endl;
    cout << "new_potential: " << new_potential << endl;
    cout << old_potential - new_potential << " <- > " << RationalType(1)/RationalType(256) << endl; 
    if((old_potential - new_potential) < 0){
      
      count_arcs_removed++;
    }
    #endif
    
    //assert( 256.0*(old_potential - new_potential) >= 1.0 );
    //assert( (old_potential - new_potential) > 1e-3);
    
    
    #if !defined(NDEBUG)
    old_potential = new_potential;
    #endif
    
    
  }
  long double final_potential = calculate_potential_function(G,x,s,q);
  int original_no_edges = G.no_of_edges/3;
  int original_no_verticies = G.no_of_verticies - original_no_edges;
  end = clock();
  time_spent = double(end - begin)/CLOCKS_PER_SEC;
  cout << "Arcs Removed = "<< count_arcs_removed << endl;
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
  
}




/** \brief Runs the potential-reduction-algorithm
 * 
 * @param G The graph on which the min-cost-flow algorithm is being carried out
 * @param x the vector x to store the primal-solution
 * @param y the vector y to store the dual-solution
 * @param s the vector s to store the dual slack variables
 * @param q the parameter q
 * @param minimum_dg
 * 
 */
template<typename IntegerType, typename RationalType>
void 
potential_reduction_algorithm(Graph<IntegerType, RationalType>& G0, 
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
  
  RationalType duality_gap = calculate_duality_gap(G, x,s);
  #ifdef VERBOSE
  cout << "potential_reduction_algorithm" << endl;
  cout << endl << "Duality Gap: " << duality_gap << endl;
  #endif
  
  std::vector<RationalType> g_prime(G.no_of_edges+1,RationalType(0)), s_prime(G.no_of_edges+1,RationalType(0)), s_hat_prime(G.no_of_edges+1,RationalType(0)), x_hat_prime(G.no_of_edges+1,RationalType(0)), z_prime(G.no_of_edges+1,RationalType(0)), x_hat( G.no_of_edges+1, RationalType(0) );
  std::vector<RationalType> z_tilde(G.no_of_edges + 1 , RationalType(0)), x_tilde(G.no_of_edges + 1, RationalType(0)), s_hat(G.no_of_edges + 1, RationalType(0));
    
  SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType> ST(G); 
  
  #ifndef NDEBUG
  long double old_potential = calculate_potential_function(G,x,s,q);
  #endif
  
  int iteration = 0;
  int primal_steps = 0;
  int dual_steps = 0;
  int electrical_flow_solver_iterations = 0;
  
  bool did_dual_step = false;
  RationalType initial_potential = calculate_potential_function(G,x,s,q);
  
  
  RationalType old_tree_condition_number(0);
  RationalType new_tree_condition_number(0);
  #ifndef NDEBUG
  RationalType avg_tree_condition_number(0);
  #endif
  #ifndef NDEBUG
  RationalType sum_tree_condition_number(0);
  #endif
  
  
  RationalType epsilon(RationalType(1));
  
  vector<IntegerType>f_0_old;
  
  bool first_itr = true; 
  
  while(duality_gap >=minimum_dg){
    iteration++;
    
    #ifdef VERBOSE
    cout << "itertion number: " << iteration << endl;
    #endif
    RationalType g_prime_sum(0);
    
    RationalType r_max(0);
    RationalType r_min(RationalType(1)/RationalType(x[1]*x[1]));
    
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      s_prime[i] = s[i]*x[i];
      g_prime[i] = q*s_prime[i]/duality_gap - 1;
      g_prime_sum += g_prime[i];
      G.unrounded_resistances[i] = 1/(x[i] * x[i]);
      
      if(RationalType(1)/(x[i]*x[i]) > r_max){
	r_max = RationalType(1)/(x[i]*x[i]); 
      }
      
      if(RationalType(1)/(x[i]*x[i]) < r_min){
	r_min = RationalType(1)/(x[i]*x[i]);
      }
      
      #ifdef VERBOSE
      cout << "resistances updation : " << (G.unrounded_resistances[i]) << "  x[i] = " << x[i] << endl;
      #endif
      
      G.batteries[i] = g_prime[i]/x[i];
      
    }
    
    RationalType rho = epsilon*r_min;
    RationalType theta = RationalType(1)/RationalType(8);
    
    RationalType delta =   RationalType(G.no_of_edges)*r_max/rho;
    RationalType phi = exp(floor(log((theta*rho)/(G.no_of_edges * sqrt(r_max)))/log(10)) * log(10) );
    
    #ifdef VERBOSE
    cout << "rho : " << rho << endl;
    cout << "phi : " << phi << endl;
    #endif
    
    #ifdef VERBOSE
    cout << theta << " *  " << rho << "/ ( " << G.no_of_edges << " *  sqrt( " << r_max << " )) = " << phi << endl;
    #endif
    
    if(first_itr){
      f_0_old = G.f_0;
    }
    
    if(!first_itr){
      RationalType val(0);
      
      for(unsigned int i = 0; i < G.non_tree_edges.size(); i++){
	arc a = G.non_tree_edges[i];
	IntegerType R_a = G.resistances_accross_cycles[i];
	
	#ifdef VERBOSE
	cout << "R_a = " << to_double(R_a) << endl;
	#endif
	IntegerType r_a = G.resistances[a];
	
	#ifdef VERBOSE
	cout << "r_a = " << to_double(r_a) << endl;
	#endif
	RationalType num =  to_double(R_a);
	RationalType den =  to_double(r_a);;
	val += num/den; //to_double(R_a)/to_double(r_a);
      }
      RationalType M = val;
      #ifdef RecordTCN
      cout << M << endl;
      #endif
      
      
      delta = M*G.no_of_edges*r_max/(2*rho);
      phi = exp(floor(log((theta*rho)/(sqrt(G.no_of_edges)* sqrt(M) * sqrt(r_max)))/log(10)) * log(10) );
      
      #ifdef Recordr_max
      cout << r_max << endl;
      #endif
      
      #ifdef Recordr_min
      cout << r_min << endl;
      #endif
      
      #ifdef VERBOSE
      cout << "delta = " << delta << endl;
      #endif
    }
    
    
    round_r(G, rho, G.unrounded_resistances);
    
    round_f(x, G, phi);
    
    if(first_itr){
      for(unsigned int i = 1; i <= G.no_of_edges; i++){
	G.currents[i] = G.f_0[i];
      }
    }
    else{
      vector<IntegerType> imbalances(G.no_of_verticies+1, 0);
      get_imbalances(x, G, imbalances);
      #ifdef VERBOSE
      cout << "imbalances got " << endl;
      #endif
      
      #ifdef WarmStart    
      warm_start_the_electrical_flow_network( G, imbalances);
      #else
      for(unsigned int i = 1; i <= G.no_of_edges; i++){
	G.currents[i] = G.f_0[i];
      }
      #endif
      
      #ifdef VERBOSE
      cout << "warm start done" << endl;
      #endif
    }
    
    first_itr = false;
    
    #ifdef VERBOSE
    for(unsigned int  i = 1; i <= G.no_of_verticies; i++){
      IntegerType old_currents_a(0);
      IntegerType f_0_a(0);
      IntegerType new_flow(0);
      for(auto a : G.incident_edges[i]){
	if(a > 0){
	  f_0_a += G.f_0[abs(a)];
	  new_flow += G.currents[abs(a)];
	}
	else{
	  f_0_a -= G.f_0[abs(a)];
	  new_flow -= G.currents[abs(a)];
	}
	
      }
      cout << "  " << to_double(f_0_a) << " " << to_double(new_flow) << endl;
    }
    cout << endl;
    #endif
    
    f_0_old = G.f_0;
    
    #ifdef VERBOSE
    cout << "1^T g' = " << g_prime_sum << ", q - m = " << q - G.no_of_edges << endl;
    #endif
    
    assert( RationalType(1e3) * fabs(g_prime_sum - q + G.no_of_edges) < RationalType(1) );
    
    #ifdef VERBOSE
    cout << "currents" << endl;
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << to_double(G.currents[i]) << " ; ";
    }
    #endif 
    
    #ifdef VERBOSE
    cout << endl << "create_low_stretch_tree called " << endl;
    #endif
    
    if( did_dual_step == false){
      if(iteration > 5 ){
	G.resistances_accross_cycles.clear();
	G.unrounded_resistances_accross_cycles.clear();
	for(unsigned int i=0; i<G.non_tree_edges.size(); i++){
	  arc a = G.non_tree_edges[i];
	  
	  G.unrounded_resistances_accross_cycles.push_back(ST.compute_unrounded_resistance_accross(a));
	}
	new_tree_condition_number = ST.calculate_tree_condition_number();
	
	#ifdef VERBOSE
	cout << "New Stretch: " << new_tree_condition_number << " " << "old stretch: " << old_tree_condition_number << endl;
	#endif
      }
      if(new_tree_condition_number >= old_tree_condition_number){
	
	#ifdef VERBOSE 
	cout << "Stretch not good" << endl;
	#endif
	
	#ifdef VERBOSE
	cout << "did_dual_step == false" << endl;
	#endif
	
	ST.clear(); 
	for(unsigned int i = 1; i <= G.no_of_verticies; i++){
	  assert(ST.node_tree_corrospondance[i].back().tree_index == 0);
	  assert(ST.node_tree_corrospondance[i].size() == 1);
	}
	
	
	G.create_low_stretch_tree_wrt_unrounded_resistances(x, ST ); 
	
	#ifndef NDEBUG
	cout << "TCN : " << ST.calculate_tree_condition_number() << endl;
	#endif
	
	#ifdef VERBOSE
	cout << "LCA vector " << endl << endl;
	for(unsigned int i = 1; i <= G.no_of_edges; i++){
	  cout << i << " : " << ST.LCA[i] << endl;
	}
	#endif
	
	assert( ST.node_tree_corrospondance[ST.root].back().tree_index == 0 );
	
	unsigned int tree_index = 1;    
	
	ST.get_tree_incident_edges();
	ST.get_initial_state_of_the_tree_decomposition();
	ST.init( ST.root, tree_index );
	
	ST.get_non_tree_edge_common_tree_index();
	
	G.resistances_accross_cycles.clear();
	for(unsigned int i=0; i<G.non_tree_edges.size(); i++){
	  arc a = G.non_tree_edges[i];
	  G.resistances_accross_cycles.push_back(ST.compute_resistance_accross(a)); //sum_over_cycle( ST, G.non_tree_edges[i], G.resistances, 1 );
	}
	
	old_tree_condition_number = ST.calculate_tree_condition_number();
	#ifndef NDEBUG
	sum_tree_condition_number += old_tree_condition_number;
	#endif
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
	cout << "TCN: " << ST.calculate_tree_condition_number() << endl; 
	#endif
	
	old_tree_condition_number = new_tree_condition_number;
	#ifndef NDEBUG
	sum_tree_condition_number += old_tree_condition_number;
	avg_tree_condition_number = sum_tree_condition_number/iteration;
	#endif
	did_dual_step = false;
	ST.clear_d();
	ST.update_sum_to_the_root(ST.root);
      }
      
    } else {
      #ifndef NDEBUG
      sum_tree_condition_number += old_tree_condition_number;
      avg_tree_condition_number = sum_tree_condition_number/iteration;
      #endif
      did_dual_step = false;
      ST.clear_d();
    }
    
    #ifdef VERBOSE
    cout << "now electrical_flow_problem would be called" << endl;
    #endif
    
    #ifdef VERBOSE
    cout << "Currents before calling the electrical_flow_problem" << endl << endl;
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << i << " ; " << to_double(G.currents[i]) << endl; 
    }
    #endif
    
    
    #ifdef VERBOSE
    cout << "resistances: " << endl; 
    for(unsigned int i = 1; i <= G.no_of_edges; i++){
      cout << i << " : " << to_double(G.resistances[i]) << endl;
    }
    cout << endl << endl;
    #endif
    
    int iteration_electrical = electrical_flow_problem(x, G0, ST, delta, old_tree_condition_number );
    
//     #ifdef VERBOSE
//     cout << "currents after the electrical flow solver ... " << endl;
//     for(unsigned int i = 1; i <= G_delta_wye.no_of_edges; i++){
//       cout << i << ": " << to_double(G.currents[i]) << endl;
//     }
//     #endif
    
    #ifdef RecordElectricalFlowIterations
    cout << iteration_electrical << endl;
    #endif
    
    electrical_flow_solver_iterations += iteration_electrical;
    
    
    ofstream file_obj;
    std::string path = "electrical_flow_log.txt";
    file_obj.open(path, std::fstream::app);
    file_obj << G.no_of_edges/3 << " " << duality_gap << " " << iteration_electrical << endl;  
    
    #ifdef VERBOSE
    cout << endl << "the electrical flow problem called" << endl;
    #endif
    
    compute_x_hat_s_hat_z_prime(G, g_prime, x_hat_prime, x_hat, s_hat_prime, s_hat, z_prime, x, phi);
    
    
    #ifdef VERBOSE
    cout << endl << "X hat prime" << endl;
    
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << x_hat_prime[i] << ";";
    }
    #endif 
    
    #ifdef VERBOSE
    cout << endl << "Z prime " << endl;
    #endif
    
    #ifdef VERBOSE
    const RationalType mu = duality_gap/q;
    cout << "mu = " << mu << endl;
    #endif
    
    #ifdef VERBOSE
    for(unsigned int i=1; i<=G.no_of_edges; i++){
      cout << z_prime[i] << " =  " << g_prime[i] << " - " << s_hat_prime[i] << ", " << s_prime[i] << ": " << s_prime[i] - mu*s_hat_prime[i] << " <-> " << mu*( 1 + z_prime[i] )<< endl;
      assert( RationalType(1e3) * fabs( s_prime[i] - mu*s_hat_prime[i] - mu*( 1 + z_prime[i] ) ) < RationalType(1) );
    }
    #endif
    
    RationalType z_prime_mod = compute_z_prime_mod(G, z_prime);
    
    #ifdef VERBOSE
    cout<<"Z_PRIME: " << z_prime_mod << endl;
    #endif
    
    if(4*z_prime_mod >= 1){
      #ifdef CHECKPRIMAL
      cout << "primal-step" << endl;
      #endif
      
      #ifdef VERBOSE
      RationalType old_potential = calculate_potential_function(G,x,s,q);
      #endif
      
      #ifdef VERBOSE
      RationalType x_tilde_zero_norm(0);
      RationalType gTx(0);
      RationalType x_tilde2(0);
      for(unsigned int i = 1; i <= G.no_of_edges; i++){
	x_tilde_zero_norm += x_hat_prime[i];
	gTx += G.batteries[i]*x_hat_prime[i]*x[i];
	x_tilde2 += x_hat_prime[i]*x_hat_prime[i];
      }
      #endif
      
      primal_step( ST, x_hat_prime, g_prime, z_prime, x_hat, x, s_hat_prime, rho);
      
      #ifdef VERBOSE      
      RationalType new_potential = calculate_potential_function(G,x,s,q);
      RationalType lambda(RationalType(1)/RationalType(4));
      #endif
      
      #ifdef RecordPotential
      RationalType new_potential = calculate_potential_function(G,x,s,q);
      cout << new_potential << endl;
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
      
      #ifdef CHECKPRIMAL
      cout << "dual-step" << endl;
      #endif
      
      #ifdef VERBOSE
      cout << " s- before: " ;
      for(auto a : s){
	cout << a << " , " ;
      }
      cout << endl;
      #endif      
      
      dual_step(G, s_prime, s, y, z_prime, s_hat_prime, s_hat, x, q, phi, duality_gap);
      
      
      #ifdef VERBOSE
      cout << "s-after: " ;
      
      for(auto a : s){
	cout << a << " , " ;
      }
      #endif
      
      dual_steps++;
      did_dual_step = true;
    }  
    
    duality_gap = update_duality_gap(G,x,s);
    
    #if !defined(NDEBUG) || defined(VERBOSE)
    const long double new_potential = calculate_potential_function(G,x,s,q);
    #endif 
    
    #ifdef PontentialReductionGapVariation
    cout << duality_gap << "	" << iteration << endl;
    #endif
    
    #ifdef VERBOSE
    cout<<"updated duality gap: "<<duality_gap<< endl;
    #endif
    
    #if !defined(NDEBUG)
    sanity_check( G, x, y, s );
    #endif
    
    #if !defined(NDEBUG) || defined(VERBOSE)
    cout << "old_potential: " << old_potential << endl;
    cout << "new_potential: " << new_potential << endl;
    cout << old_potential - new_potential << " <- > " << RationalType(1)/RationalType(256) << endl; 
    #endif
    
    //assert( 256.0*(old_potential - new_potential) >= 1.0 );
    //assert( (old_potential - new_potential) > 1e-3);
    
    
    #if !defined(NDEBUG)
    old_potential = new_potential;
    #endif
    
    
  }
  
  #ifndef NDEBUG
  long double final_potential = calculate_potential_function(G,x,s,q);
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
