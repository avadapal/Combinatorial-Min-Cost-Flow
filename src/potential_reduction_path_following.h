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
#include "electrical_flow_problem_unrounded.h"
#include "electrical_flow_problem_DS.h"
#include "verbose_prints.h"


using namespace std;


// checks central path condition
template<typename Network>
bool check_central_path_condition(Network & N,
                                  typename Network::RationalType Mu,
                                  typename Network::RationalType delta
                                 ) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  typedef typename Network::IntegerType IntegerType;
  Graph& G0 = N.G;

  cout << "check_central_path_condition" << endl;
  cout << "Mu = " << to_double(Mu) << endl;
  RationalType LHS(0);
  for (unsigned int a = 1; a <= G0.no_of_edges; a++) {
    IntegerType xasa = N.arcdata[a].xlower * N.arcdata[a].slower;
    cout << "xasa = " << to_double(xasa) << endl;
    if (xasa != 0) {
      LHS += fabs( to_double(xasa) - Mu);
    }
    xasa = (N.arcdata[a].capacity - N.arcdata[a].xlower) * N.arcdata[a].supper;
    cout << "xasa = " << to_double(xasa) << endl;
    if (xasa != 0) {
      LHS += fabs( to_double(xasa) - Mu);
    }
    xasa = N.arcdata[a].infeasibility * N.arcdata[a].sroof;
    cout << "xasa = " << to_double(xasa) << endl;
    if (xasa != 0) {
      LHS += fabs( to_double(xasa) - Mu);
    }
  }
  RationalType RHS = delta * Mu;

  cout << LHS << " <-> " << RHS << endl;
  
  return LHS <= RHS;
}


template<typename Network>
void check_xasa_bounds(Network & N, typename Network::RationalType delta, typename Network::RationalType mu) {
  typedef typename Network::IntegerType IntegerType;

  for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
    IntegerType xasa = N.arcdata[a].xlower * N.arcdata[a].slower;
    assert(xasa == IntegerType(0) || 
       (to_double(xasa) >= (1 - delta) * to_double(mu) && to_double(xasa) <= (1 + delta) * to_double(mu)));

    xasa = (N.arcdata[a].capacity - N.arcdata[a].xlower) * N.arcdata[a].supper;
    assert(xasa == IntegerType(0) || 
          (to_double(xasa) >= (1 - delta) * to_double(mu) && to_double(xasa) <= (1 + delta) * to_double(mu)));

    xasa = N.arcdata[a].infeasibility * N.arcdata[a].sroof;
    assert(xasa == IntegerType(0) || 
	   (to_double(xasa) >= (1 - delta) * to_double(mu) && to_double(xasa) <= (1 + delta) * to_double(mu)));
  }
}
 
template<typename Network>
void check_lower_bound_invariant(Network & N,
                                 typename Network::RationalType delta,
                                 typename Network::IntegerType beta,
                                 typename Network:: IntegerType gamma
                                ) {
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;

  unsigned int m = RationalType(3) * N.G.no_of_edges;

  RationalType primal_lower_bound = ((1 - delta) * to_double(beta)) / ((1 + delta) * m);
  RationalType dual_lower_bound = ((1 - delta) * to_double(gamma)) / ((1 + delta) * m);

  for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
    IntegerType xlower = N.arcdata[a].xlower;
    IntegerType slower = N.arcdata[a].slower;
    assert(xlower * slower == IntegerType(0) || 
	   (to_double(xlower) >= primal_lower_bound && to_double(slower) >= dual_lower_bound));

    IntegerType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    IntegerType supper = N.arcdata[a].supper;
    assert(xupper * supper == IntegerType(0) || 
	   (to_double(xupper) >= primal_lower_bound && to_double(supper) >= dual_lower_bound));

    IntegerType xroof = N.arcdata[a].infeasibility;
    IntegerType sroof = N.arcdata[a].sroof;
    assert(xroof * sroof == IntegerType(0) ||
	  (to_double(xroof) >= primal_lower_bound && to_double(sroof) >= dual_lower_bound));
  }
}




template<typename Network>
void get_initial_currents(Network& N,
                          typename Network::RationalType mu_prime
                         ) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;
  Graph& G0 = N.G;

  for (unsigned int a = 1; a <= G0.no_of_edges; a++) {

    if (N.arcdata[a].xlower * N.arcdata[a].slower == IntegerType(0)) {
      N.arcdata[a].current_lower = IntegerType(0);
      N.arcdata[a].initial_current_lower = IntegerType(0);
    } else {
      //long long current_lower_ = to_double(N.arcdata[a].xlower) - round(mu_prime / to_double(N.arcdata[a].slower) );
      N.arcdata[a].current_lower = N.arcdata[a].xlower - convert_rounded_ratio_to_integer(mu_prime, to_double(N.arcdata[a].slower));
      N.arcdata[a].initial_current_lower = N.arcdata[a].current_lower;
    }


    IntegerType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    IntegerType current_upper;
    if (xupper * N.arcdata[a].supper == IntegerType(0)) {
      current_upper = IntegerType(0);
    } else {
      current_upper = xupper - convert_rounded_ratio_to_integer(mu_prime, to_double(N.arcdata[a].supper));
    }

    N.arcdata[a].cur_src_vw = N.arcdata[a].initial_current_lower + current_upper;

    if (N.arcdata[a].infeasibility * N.arcdata[a].sroof == IntegerType(0)) {
      N.arcdata[a].current_roof = IntegerType(0);
      N.arcdata[a].initial_current_roof = N.arcdata[a].current_roof;
    } else {
      N.arcdata[a].current_roof = N.arcdata[a].infeasibility - convert_rounded_ratio_to_integer(mu_prime, to_double(N.arcdata[a].sroof));
      N.arcdata[a].initial_current_roof = N.arcdata[a].current_roof;
    }
  }
}

template<typename Network>
unsigned int compute_no_edges_in_C(Network& N,
                                   typename Network::RationalType THRESHOLD_X,
                                   typename Network::RationalType THRESHOLD_S
                                  )
{
  typedef typename Network::RationalType RationalType;
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;
  unsigned int no_edges_in_C = 0;
  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    if (N.arcdata[a].xlower > THRESHOLD_X && N.arcdata[a].slower > THRESHOLD_S)
    {
      no_edges_in_C++;
    }

    RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    if (xupper > THRESHOLD_X && N.arcdata[a].supper > THRESHOLD_S)
    {
      no_edges_in_C++;
    }

    if (N.arcdata[a].infeasibility > THRESHOLD_X && N.arcdata[a].sroof > THRESHOLD_S)
    {
      no_edges_in_C++;
    }
  }

  return no_edges_in_C;
}


template<typename IntegerType, typename RationalType>
unsigned int compute_no_edges_in_C(
  Graph<IntegerType, RationalType> &G,
  vector<RationalType> &x,
  vector<RationalType> &s,
  RationalType THRESHOLD_X,
  RationalType THRESHOLD_S) {
  unsigned int no_edges_in_C = 0;
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    if (x[i] >= THRESHOLD_X && s[i] >= THRESHOLD_S) {
      no_edges_in_C++;
    }
  }
  return no_edges_in_C;
}


template<typename IntegerType, typename T>
long double logarithm_base(IntegerType base, T number)
{
  return (log(number) / log(base));
}

/** \brief Does a primal sanity check
 *
 * Checks if the primal constranits are satisfied
 *
 * @param G The graph on which sanity check is being done
 * @param x The primal solution
 */
template <
  typename Network,
  typename RationalType
  >
void
primal_sanity_check(
  const Network& N,
  const std::vector<RationalType> &x
) {
  // typedef typename Network::Graph Graph;
//  typedef typename Network::IntegerType IntegerType;
//  typedef typename Network::RationalType RationalType;
  const auto& G = N.G;
#ifdef VERBOSE
  cout << "primal_sanity_check " << endl;
#endif

  // this is going to be the total error
  RationalType err_primal_constraint = 0.0;

  for (unsigned int v = 1; v <= G.no_of_verticies; v++) {
    RationalType adj_flow = 0.0;

    // compute the excess at node v
    for (unsigned int j = 0; j < G.incident_edges[v].size(); j++ ) {

      int edge = G.incident_edges[v][j];

      if (edge < 0) {
        adj_flow += x[-edge];
      }
      else {
        adj_flow -= x[edge];
      }
    }

    const RationalType error = fabs( adj_flow - N.nodedata[v].demand );
    err_primal_constraint += error;
    if (  error  > 1e-3 ) {
#ifdef VERBOSE
      cout << "primal error of node " << v << ": " << flow << " " << G.demands[v] << " " << error << endl;
#endif
      //assert( abs(error/max( static_cast<RationalType>( 1.0 ),fabs(G.demands[i]))) <= 1e-3 );
    }

  }
#ifdef VERBOSE
  cout << "Primal Error: " << err_primal_constraint << endl << endl;
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
template <
  typename Network,
  typename RationalType
  >
void dual_sanity_check(
  const Network& N,
  const std::vector<RationalType> &y,
  const std::vector<RationalType> &s,
  const vector<RationalType>& x,
  RationalType THRESHOLD_X
)
{
  typedef typename Network::GraphType Graph;
//  typedef typename Network::IntegerType IntegerType;
//  typedef typename Network::RationalType RationalType;
  const Graph& G = N.G;

#ifdef VERBOSE
  cout << "dual_sanity_check" << endl;
#endif
  RationalType err_dual_constraint = 0;
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    if (x[i] >= THRESHOLD_X )
    {
      const RationalType expected_dual_soln = N.arcdata[i].cost + y[G.tails[i]] - y[G.heads[i]];
      //cout << G.costs[i] << " + " << y[G.tails[i]] << " - " << y[G.heads[i]] << " = " << s[i] << endl;
      const RationalType error = fabs( s[i] - expected_dual_soln );
      err_dual_constraint += error;
#ifdef VERBOSE
      cout << s[i] << " " << expected_dual_soln << " " << err_dual_constraint << " " << error << endl;
#endif

      if (1000 * error > 1) {
        cout << " i = " << i << " ,  " << s[i] << " " << expected_dual_soln << " " << err_dual_constraint << " " << error << endl;
      }
      assert( 1000 * error < 1 );
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
template <
  typename Network,
  typename RationalType
  >
void sanity_check(
  const Network& N,
  const std::vector<RationalType> &x,
  const std::vector<RationalType> &y,
  const std::vector<RationalType> &s,
  RationalType THRESHOLD_X
)
{
  //The Primal Constraints
  primal_sanity_check( N, x );
  //The Dual Constraints
  dual_sanity_check( N, y, s, x, THRESHOLD_X );
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
    unsigned int q,
    RationalType THRESHOLD_X
                                        )
{
#ifdef VERBOSE
  cout << "calculate_potential_function()" << endl;
#endif
  RationalType x_t_s = 0;
  long double ln_x_s = 0;
  int no_of_edges_removed = 0;
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {

    if (x[i] >= THRESHOLD_X)
    {

      const RationalType xasa = x[i] * s[i];
      x_t_s += xasa;

      ln_x_s += log( xasa );

#ifdef VERBOSE
      cout <<  " log( " << x[i] << " * " << s[i] << ") = " << log(x[i]*s[i]) << endl;
#endif
    }
    else {
      no_of_edges_removed++;
    }

  }

  const unsigned int m = G.no_of_edges - no_of_edges_removed;

  const long double mlogm = m * log(m);
  const long double potential = q * log(x_t_s) - ln_x_s - mlogm;
#ifdef VERBOSE
  cout << "x_t_s: " << x_t_s << endl;
  cout << "Potential: " << q << " * " << natural_logarithm(x_t_s) << " - " << ln_x_s << " - " << mlogm << " = " << potential << endl;
#endif
  return potential;
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
  for ( auto x_i : x )  {
    if ( abs(x_i) > max ) {
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
// template<typename IntegerType, typename RationalType>
// void primal_step( const SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType>& ST,
//      std::vector<RationalType>&h,
//      std::vector<RationalType>&x,
//       std::vector<RationalType>&s,
//      RationalType phi,
//       RationalType rho,
//      RationalType THRESHOLD_X
//    )
// {
//
//   const Graph<IntegerType, RationalType>& G = ST.G;
//
//   // print x, h
//   #ifdef VERBOSE
//     for(unsigned int i = 1; i <= G.no_of_edges; i++){
//       cout << "x["<< i << "] = " << x[i] << " , " << "h[" << i << "] = " << h[i] << endl;
//     }
//   #endif
//
//   for(unsigned int i = 1; i <= G.no_of_edges; i++){
//     if(x[i] >= THRESHOLD_X){
//       // update x by -h (Note: Reverse as in the paper, but pm(phi0-phi))
//       x[i] -= h[i];
//       // numerical issues switch
// //       if(x[i] < 0 && x[i] > -1e-5){
// //         x[i] = -x[i];
// //       }
//
//       if(x[i] < 0)
//       {
//  cout << "x[" << i << "] = " << x[i] << endl;
//       }
//       assert(x[i] >= 0);
//     }
//     else{
//       if(h[i] >= 1e-5) {
//         cout << "h[i]: " << h[i] << endl;
//       }
//       assert(h[i] < 1e-5);
//     }
//   }
// }

template<typename Network>
void arc_fixation(
  Network& N,
  typename Network::IntegerType beta,
  typename Network::IntegerType gamma,
  typename Network::RationalType delta,
  vector<typename Network::IntegerType>& demands
)
{
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  typedef typename Network::IntegerType IntegerType;
  Graph& G0 = N.G;
  unsigned int m = 3 * G0.no_of_edges;

  RationalType THRESHOLD_X = ( (1 - delta) * to_double(beta) ) / ( ( 1 + delta) * m );
  RationalType THRESHOLD_S = ( (1 - delta) * to_double(gamma) ) / ( ( 1 + delta) * m );

  cout << "THRESHOLD_X = " << THRESHOLD_X << endl;
  cout << "THRESHOLD_S = " << THRESHOLD_S << endl;
  
  for (unsigned int a = 1; a <= G0.no_of_edges; a++) {
    IntegerType xasa = N.arcdata[a].xlower * N.arcdata[a].slower;
    if (xasa != 0) {
      
      cout << "xlower = " << to_double(N.arcdata[a].xlower) << endl;
      
      if (to_double(N.arcdata[a].xlower) < THRESHOLD_X) {
	
	if(N.arcdata[a].direction == 1)
	{
	  demands[G0.tails[a]] += N.arcdata[a].xlower;
	}
	else
	{
	  demands[G0.heads[a]] += N.arcdata[a].xlower;
	}
	
        cout << "arc-fixation done, xlower set to zero" << endl;
        N.arcdata[a].xlower = 0;
      }
      if ( to_double(N.arcdata[a].slower) < THRESHOLD_S) {
        cout << "arc-fixation done, slower set to zero" << endl;
        N.arcdata[a].slower = 0;
      }
    }

    
    IntegerType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    cout << "xupper = " << to_double(xupper) << endl;
    xasa = xupper * N.arcdata[a].supper;
    if (xasa != 0) {
      if ( to_double(xupper) < THRESHOLD_X) {
        cout << "arc-fixation done, xupper set to zero" << endl;	
	if(N.arcdata[a].direction == 1)
	{
	 demands[G0.heads[a]] += xupper;
	}
	else{
	  demands[G0.tails[a]] += xupper;
	}
	
	 N.arcdata[a].xlower = N.arcdata[a].capacity;
      }
      if ( to_double(N.arcdata[a].supper) < THRESHOLD_S) {
        cout << "arc-fixation done, supper set to zero" << endl;
        N.arcdata[a].supper = 0;
      }
    }

    xasa = N.arcdata[a].infeasibility * N.arcdata[a].sroof;
    if (xasa != 0)
    {
      cout << "inf = " << to_double(N.arcdata[a].infeasibility) << endl;
      if ( to_double(N.arcdata[a].infeasibility) < THRESHOLD_X) {
        cout << "arc-fixation done, infeasibility set to zero" << endl;
       
	if(N.arcdata[a].direction == 1)
	{
	  demands[G0.tails[a]] += N.arcdata[a].infeasibility;
	  demands[G0.heads[a]] -= N.arcdata[a].infeasibility;
	}
	else
	{
	  demands[G0.heads[a]] += N.arcdata[a].infeasibility;
	  demands[G0.tails[a]] -= N.arcdata[a].infeasibility;
	}
	
	 N.arcdata[a].infeasibility = 0;
      }
      if ( to_double(N.arcdata[a].sroof) < THRESHOLD_S) {
        cout << "arc-fixation done, sroof set to zero" << endl;
        N.arcdata[a].sroof = 0;
      }
    }
  }
}


template<typename Network>
void primal_step(Network& N) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;

  Graph& G0 = N.G;


  for (unsigned int a = 1; a <= G0.no_of_edges; a++) {
    if (N.arcdata[a].xlower != 0) {
      IntegerType Diff = N.arcdata[a].initial_current_lower - N.arcdata[a].current_lower;
      IntegerType h_lower =  (Diff);
      N.arcdata[a].xlower -= h_lower;
      assert(N.arcdata[a].xlower > 0);
    }


    if (N.arcdata[a].infeasibility != 0){
      IntegerType Diff = N.arcdata[a].initial_current_roof - N.arcdata[a].current_roof;
      IntegerType h_roof = (Diff);
      N.arcdata[a].infeasibility -= h_roof;
      assert(N.arcdata[a].infeasibility > 0 || N.arcdata[a].infeasibility == 0);
    }
  }
}



template<typename Network>
void update_potentials(Network& N) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;

  Graph& G0 = N.G;

  for (unsigned int v = 1; v <= G0.no_of_vertices; v++)
  {
    IntegerType TIV = N.nodedata[v].voltage;
    N.nodedata[v].potential -= TIV;
  }
  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    IntegerType TIV_vw = N.arcdata[a].voltage_vw;
    N.arcdata[a].potentialvw -= TIV_vw;
  }
}

template<typename Network>
void update_slacks_using_potentials(Network& N) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;

  Graph& G0 = N.G;

  for (unsigned int i = 1; i <= G0.no_of_edges; i++) {
    if (N.arcdata[i].xlower != 0) {
      IntegerType difference = 0;
      if (N.arcdata[i].direction == 1 || N.arcdata[i].direction == 0) {
        difference = N.nodedata[G0.tails[i]].potential - N.arcdata[i].potentialvw;
      }
      if (N.arcdata[i].direction == -1) {
        difference = N.nodedata[G0.heads[i]].potential - N.arcdata[i].potentialvw;
      }
      N.arcdata[i].slower = N.arcdata[i].cost + difference;
    }

    if (N.arcdata[i].capacity - N.arcdata[i].xlower != 0) {
      IntegerType difference = 0;

      if (N.arcdata[i].direction == 1 || N.arcdata[i].direction == 0)
      {
        difference = N.nodedata[G0.heads[i]].potential - N.arcdata[i].potentialvw;
      }
      if (N.arcdata[i].direction == -1) {
        difference = N.nodedata[G0.tails[i]].potential - N.arcdata[i].potentialvw;
      }
      N.arcdata[i].supper = difference;
    }

    if (N.arcdata[i].infeasibility != 0) {
      IntegerType difference = N.nodedata[G0.tails[i]].potential - N.nodedata[G0.heads[i]].potential;
      N.arcdata[i].sroof = N.arcdata[i].croof + difference;
    }
  }
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

  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    RationalType a = r[i];
    if (a < min) {
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
                                    const std::vector<RationalType>& s,
                                    RationalType THRESHOLD_X,
                                    RationalType THRESHOLD_S
                                  )
{
  RationalType duality_gap = 0;

  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    if (x[i] < THRESHOLD_X || s[i] < THRESHOLD_S ) continue;
    const RationalType xasa = x[i] * s[i];
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
  for ( unsigned a = 1; a <= G.no_of_edges; ++a ) {
    round_resistance( rho, resistances[a], G.resistances[a] );
  }
#ifdef VERBOSE
  cout << "resistances: " << endl;
  for (unsigned int a = 1; a <= G.no_of_edges; a++) {
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
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    if (x[i] != 0) {
      round_flow( phi, G.batteries[i], G.unrounded_resistances[i], G.f_0[i] );
    }
  }
}




template<typename IntegerType, typename RationalType>
void round_g(Graph<IntegerType, RationalType>& G,
             RationalType phi
            )
{
  for (unsigned int i = 1 ; i <= G.no_of_edges; i++) {
    G.batteries_unrounded[i] =  G.batteries[i];
    G.g_tilde[i] = to_double(phi) * to_double(G.resistances[i] * G.f_0[i]);
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
                    vector<IntegerType>& imbalances,
                    RationalType THRESHOLD_X
                   )
{
  for (unsigned int a = 1; a <= G.no_of_edges; a++) {
    if (x[a] >= THRESHOLD_X) {
      const IntegerType alpha = G.f_0[a] - G.currents[a];
      imbalances[G.head(a)] -= alpha;
      imbalances[G.tail(a)] += alpha;
    }
  }
}


template<typename IntegerType, typename RationalType>
void update_steps()
{
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
                      vector<T>& currents
                     )
{
  for (auto a :  boost::adaptors::reverse(G.tree_edges)) {

#ifdef VERBOSE
    cout << "a: " << a << endl;
#endif

    // Depending on sign(a), x[pm a] +-=delta
    if ( a < 0 ) {
      T delta = imbalances[G.tail(a)];
      const node v = G.tails[-a];
      const node w = G.heads[-a];

#ifdef VERBOSE
      cout << "before: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
      cout << "current before = " << to_double(G.currents[-a]) << endl;
#endif

      currents[-a] -= delta;
      imbalances[v] += delta;
      imbalances[w] -= delta;
#ifdef VERBOSE
      cout << "after: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
#endif

      assert( imbalances[w] == 0 );
    } else {

      T delta = imbalances[G.tail(a)];
      const node w = G.tails[a];
      const node v = G.heads[a];
#ifdef VERBOSE
      cout << "before: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
      cout << "current before = " << to_double(G.currents[a]) << endl;
#endif

      currents[a] += delta;
      imbalances[v] += delta;
      imbalances[w] -= delta;

#ifdef VERBOSE
      cout << "after: " << endl;
      cout << v << " : " << to_double(imbalances[v]) << endl;
      cout << w << " : " << to_double(imbalances[w]) << endl << endl;
#endif

      assert( imbalances[w] == 0 );
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


template<typename Network>
typename Network::RationalType get_r_max(Network& N)
{
  typedef typename Network::RationalType RationalType;
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;

  RationalType r_max(0);

  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    if (r_max > N.arcdata[a].xlower / N.arcdata[a].slower)
    {
      r_max = N.arcdata[a].xlower / N.arcdata[a].slower;
    }

    RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    if (r_max > xupper / N.arcdata[a].supper)
    {
      r_max = xupper / N.arcdata[a].supper;
    }

    if (r_max > N.arcdata[a].infeasibility / N.arcdata[a].sroof)
    {
      r_max = N.arcdata[a].infeasibility / N.arcdata[a].sroof;
    }
  }

  return r_max;
}

template<typename Network>
void sigma_assert(Network& N,
                  typename Network::RationalType Mu,
                  typename Network::RationalType delta
                 )
{

  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;
  // Assert that sigma^2 is bounded and x_a s_a/mu constant
  RationalType sigma_2 = RationalType(0);

  for (unsigned int a = 1; a <= G0.no_of_edges; a++) {
    {

      sigma_2 +=  ( (to_double(N.arcdata[a].slower * N.arcdata[a].xlower) / to_double(Mu)) - RationalType(1) )
                  * ( ( to_double(N.arcdata[a].slower * N.arcdata[a].xlower) / to_double(Mu)) - RationalType(1) );
      assert( to_double(N.arcdata[a].xlower * N.arcdata[a].slower) / to_double(Mu) > RationalType(1) / RationalType(2));
      if ( to_double(N.arcdata[a].xlower * N.arcdata[a].slower) / to_double(Mu) >= RationalType(3) / RationalType(2)) {

        cout << "x_is_i/Mu = " << to_double(N.arcdata[a].xlower * N.arcdata[a].slower) / to_double(Mu) << endl;
      }
      cout << "xasa = " << to_double(N.arcdata[a].xlower * N.arcdata[a].slower) << endl;
      assert( to_double(N.arcdata[a].xlower * N.arcdata[a].slower) / to_double(Mu) < RationalType(3) / RationalType(2));

    }

    IntegerType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;

    {

      sigma_2 +=  ( ( to_double(N.arcdata[a].supper * xupper) / to_double(Mu)) - RationalType(1) )
                  * ( ( to_double(N.arcdata[a].supper * xupper) / to_double(Mu)) - RationalType(1) );
      assert( to_double(xupper * N.arcdata[a].supper) / to_double(Mu) > RationalType(1) / RationalType(2));
      if (to_double(xupper * N.arcdata[a].supper) / to_double(Mu) >= RationalType(3) / RationalType(2)) {

        cout << "x_is_i/Mu = " << to_double(N.arcdata[a].xlower * N.arcdata[a].slower) / to_double(Mu) << endl;
      }
      cout << "xasa = " << to_double(xupper * N.arcdata[a].supper) << endl;
      assert( to_double(xupper * N.arcdata[a].supper) / to_double(Mu) < RationalType(3) / RationalType(2));

    }


    {

      sigma_2 +=  ( ( to_double(N.arcdata[a].sroof * N.arcdata[a].infeasibility) / to_double(Mu) ) - RationalType(1) )
                  * ( ( to_double(N.arcdata[a].sroof * N.arcdata[a].infeasibility) / to_double(Mu)) - RationalType(1) );
      assert( to_double(N.arcdata[a].infeasibility * N.arcdata[a].sroof) / to_double(Mu) > RationalType(1) / RationalType(2));
      if (to_double(N.arcdata[a].infeasibility * N.arcdata[a].sroof) / to_double(Mu) >= RationalType(3) / RationalType(2)) {

        cout << "x_is_i/Mu = " << to_double(N.arcdata[a].infeasibility * N.arcdata[a].sroof) / to_double(Mu) << endl;
      }
      cout << "xasa = " << to_double(N.arcdata[a].infeasibility * N.arcdata[a].sroof) << endl;
      assert( to_double(N.arcdata[a].infeasibility * N.arcdata[a].sroof) / to_double(Mu) < RationalType(3) / RationalType(2));

    }

    if (sigma_2 > delta * delta) {

      cout << "sigma_2 = " << sigma_2 << endl;
    }
    assert(sigma_2 <= delta * delta);
  }

}

template<typename Network>
typename Network::RationalType get_r_min(Network& N)
{
  typedef typename Network::RationalType RationalType;
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;

  RationalType r_min = numeric_limits<RationalType>::max();

  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    if (r_min < N.arcdata[a].xlower / N.arcdata[a].slower)
    {
      r_min = N.arcdata[a].xlower / N.arcdata[a].slower;
    }

    RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    if (r_min < xupper / N.arcdata[a].supper)
    {
      r_min = xupper / N.arcdata[a].supper;
    }

    if (r_min < N.arcdata[a].infeasibility / N.arcdata[a].sroof)
    {
      r_min = N.arcdata[a].infeasibility / N.arcdata[a].sroof;
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
  for (unsigned int i = 0; i < G0.non_tree_edges.size(); i++) {
    arc a = G0.non_tree_edges[i];
    IntegerType R_a = G0.resistances_accross_cycles[i];
    IntegerType r_a;

#ifdef VERBOSE
    cout << "a = " << a << endl;
    cout << "R_a = " << to_double(R_a) << endl;
#endif

    if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
      r_a = N.arcdata[a].resistance_lower;
    }
    else {
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
    M += num / den;

    IntegerType r_a_roof = N.arcdata[a].resistance_roof;
    IntegerType R_a_roof = R_a - (N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper) + N.arcdata[a].resistance_roof;
    M += to_double(R_a_roof) / to_double(r_a_roof);
  }
#ifdef VERBOSE
  cout << "M midway = " << M << endl;
#endif
  for (auto a : G0.tree_edges)
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
    if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
        N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper)
    {
      if (N.arcdata[abs(a)].resistance_lower > N.arcdata[abs(a)].resistance_upper)
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
    M += to_double(R_a) / to_double(r_a);
  }

  cout << "M = " << M << endl;
  return M;
}

template<typename Network>
void get_resistances(Network& N)
{
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;

  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    //NOTE: We are ceiling here, for the case when path-following is done with scaled demands, capacity, cost.

    long long resistance_lower_ = ceil( to_double(N.arcdata[a].slower) / to_double(N.arcdata[a].xlower));
    N.arcdata[a].resistance_lower = resistance_lower_;
    long long resistance_upper_ = ceil( to_double(N.arcdata[a].supper) / to_double(N.arcdata[a].capacity - N.arcdata[a].xlower));
    N.arcdata[a].resistance_upper = resistance_upper_;
    long long resistance_roof_ = ceil( to_double(N.arcdata[a].sroof) / to_double(N.arcdata[a].infeasibility));
    if(N.arcdata[a].infeasibility == 0 ) resistance_roof_ = numeric_limits<long long>::max();
    cout << "resistance_roof_ = " << resistance_roof_ << " = " << to_double(N.arcdata[a].sroof) << "/" 
				  << to_double(N.arcdata[a].infeasibility) << endl;
				  
    N.arcdata[a].resistance_roof  = resistance_roof_;
  }
}

template<typename Network>
void get_batteries(Network& N, typename Network::RationalType Mu)
{
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;

  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    N.arcdata[a].battery_lower = N.arcdata[a].slower - (Mu / N.arcdata[a].xlower);
    N.arcdata[a].battery_upper = N.arcdata[a].supper - (Mu / (N.arcdata[a].capacity - N.arcdata[a].xlower));
    N.arcdata[a].battery_roof = N.arcdata[a].sroof - (Mu / N.arcdata[a].infeasibility);
  }
}


template<typename Network>
typename Network::RationalType find_sigma_2(Network& N,
    typename Network::RationalType Mu,
    typename Network::RationalType t
                                           )
{
  typedef typename Network::RationalType RationalType;
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;

  RationalType sigma_2(0);

  cout << "MU = " << Mu << endl;
  cout << "upper = " << RationalType(1) - RationalType(1) / (RationalType(10) * sqrt(3 * G0.no_of_edges)) << endl;

  RationalType beta = find_beta(N);
  RationalType gamma = find_gamma(N);
  RationalType delta = RationalType(1) / RationalType(10);

  unsigned int m = RationalType(3) * G0.no_of_edges;

  RationalType THRESHOLD_X = (1 - delta) * beta / ((1 + delta) * m);
  RationalType THRESHOLD_S = (1 - delta) * gamma / ((1 + delta) * m);
  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    RationalType xlower = N.arcdata[a].xlower;
    assert(xlower > THRESHOLD_X || xlower == RationalType(0));
    cout << "xlower = " << xlower << endl;
    RationalType slower = N.arcdata[a].slower;
    assert(slower > THRESHOLD_S || slower == RationalType(0));
    cout << "slower = " << slower << endl;
    RationalType xasa = xlower * slower;
    assert(xasa >= t);
    RationalType xasa_over_Mu = xasa / Mu;
    cout << "xasa_over_Mu = " << xasa_over_Mu << endl;
    assert(xasa_over_Mu > 1 - RationalType(1) / (RationalType(10) * sqrt(3 * G0.no_of_edges)));
    assert(xasa_over_Mu <= RationalType(1));
    sigma_2 += (xlower * slower / Mu - RationalType(1)) * (xlower * slower / Mu - RationalType(1));
    cout << "sigma_2 = " << sigma_2 << endl;

    RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    assert(xupper > THRESHOLD_X || xupper == RationalType(0));

    cout << "xupper = " << xupper << endl;
    RationalType supper = N.arcdata[a].supper;
    cout << "supper = " << supper << endl;
    assert(supper > THRESHOLD_S || supper == RationalType(0));
    xasa = xupper * supper;
    assert(xasa >= t);
    xasa_over_Mu = xasa / Mu;
    cout << "xasa_over_Mu = " << xasa_over_Mu << endl;
    assert(xasa_over_Mu <= RationalType(1));
    assert(xasa_over_Mu > 1 - RationalType(1) / (RationalType(10) * sqrt(3 * G0.no_of_edges)));
    sigma_2 += (xupper * supper / Mu - RationalType(1)) * (xupper * supper / Mu - RationalType(1));
    cout << "sigma_2 = " << sigma_2 << endl;

    RationalType xroof = N.arcdata[a].infeasibility;
    cout << "xroof = " << xroof << endl;
    RationalType sroof = N.arcdata[a].sroof;
    cout << "sroof = " << sroof << endl;
    assert(xroof > THRESHOLD_X || xroof == RationalType(0));
    assert(sroof > THRESHOLD_S || sroof == RationalType(0));
    cout << "sigma_2 = " << sigma_2 << endl;

    xasa_over_Mu = xroof * sroof / Mu;

    cout << "xasa_over_Mu = " << xasa_over_Mu << endl;
    assert(xasa_over_Mu <= RationalType(1));
    assert(xasa_over_Mu > 1 - RationalType(1) / (RationalType(10) * sqrt(3 * G0.no_of_edges)));
    sigma_2 += (xroof * sroof / Mu - RationalType(1)) * (xroof * sroof / Mu - RationalType(1));
  }

  return sigma_2;
}

template <
  typename Network
  > void clear_currents(Network& N)
{
  typedef typename Network::IntegerType IntegerType;
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++)
  {
    N.arcdata[a].current = IntegerType(0);
  }
}

template<typename Network> void check_asserts_before_xasa(Network& N,
    typename Network::IntegerType t,
    typename Network::IntegerType mu,
    typename Network::RationalType delta) {
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;

  for (unsigned int a = 1; a <= G0.no_of_edges; a++) {
    IntegerType xasa = N.arcdata[a].xlower * N.arcdata[a].slower;

    IntegerType gamma = find_gamma(N);
    IntegerType beta  = find_beta(N);
    IntegerType U = find_U(N, beta);
    IntegerType C = find_C(N, gamma);
    long long M_ = max(to_double(U) , to_double(C));
    IntegerType M = M_;

    assert(to_double(xasa) / to_double(mu) > to_double(t) / to_double(mu));
    assert(to_double(xasa) / to_double(mu) > 1 - to_double(beta * gamma * M) / to_double(t + beta * gamma * M) );
    assert(to_double(xasa) / to_double(mu) > 1  - to_double(M) * delta / (2 * G0.no_of_edges * 3 * to_double(U * C)));
    assert(to_double(xasa) / to_double(mu) > 1  - delta / (sqrt(3 * G0.no_of_edges)));

    xasa = (N.arcdata[a].capacity - N.arcdata[a].xlower) * N.arcdata[a].supper;
    assert(to_double(xasa) / to_double(mu) >= to_double(t) / to_double(mu));
    assert(to_double(xasa) / to_double(mu) >= 1 - to_double(beta * gamma * M) / to_double(t + beta * gamma * M) );
    assert(to_double(xasa) / to_double(mu) >= 1  - to_double(M) * delta / (2 * G0.no_of_edges * 3 * to_double(U * C)));
    assert(to_double(xasa) / to_double(mu) >= 1  - delta / (sqrt(3 * G0.no_of_edges)));

    xasa = N.arcdata[a].infeasibility * N.arcdata[a].sroof;
    assert(to_double(xasa) / to_double(mu) >= to_double(t) / to_double(mu));
    assert(to_double(xasa) / to_double(mu) >= 1 - to_double(beta * gamma * M) / to_double(t + beta * gamma * M) );
    assert(to_double(xasa) / to_double(mu) >= 1  - to_double(M) * delta / (2 * G0.no_of_edges * 3 * to_double(U * C)));
    assert(to_double(xasa) / to_double(mu) >= 1  - delta / (sqrt(3 * G0.no_of_edges)));
  }
}

/** \brief Runs the path-following-algorithm
 *
 * @param G The graph on which the min-cost-flow algorithm is being carried out
 * @param x the vector x to store the primal-solution
 * @param y the vector y to store the dual-solution
 * @param s the vector s to store the dual slack variables
 * @param q the parameter q
 * @param t The parameter t
 * @param minimum_dg
 *
 */
template <
  typename Network
  > void path_following_algorithm(
    Network& N,
    typename Network::IntegerType t,
    typename Network::IntegerType minimum_dg,
    typename Network::IntegerType beta,
    typename Network::IntegerType gamma,
    typename Network::RationalType delta,
    typename Network::IntegerType mu0,
    typename Network::RationalType tau
  )
{

  cout << "Path Following Algorithm" << endl;
  typedef typename Network::GraphType Graph1;
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;

  Graph1& G0 = N.G;

  // compute duality gap
  IntegerType duality_gap = update_duality_gap(N);
  RationalType mu         = to_double(mu0);

  // iteration counters
  int iteration = 0;




  vector<IntegerType> demands(G0.no_of_vertices + 1, 0);
  copy_demands_to_vector(N, demands);

  #ifndef NDEBUG
  primal_feasibility_check(N, demands);
  dual_feasibility_check(N);
  #endif
  
  cout << "minimum duality gap = " << to_double(minimum_dg) << endl;
  cout << "initial duality gap = " << to_double(duality_gap) << endl;
  cout << "beta = " << to_double(beta) << endl;
  cout << "gamma = " << to_double(gamma) << endl;


  // MAIN LOOP
  while (duality_gap > minimum_dg) {
    iteration++;
    cout << "iteration " << iteration << ", duality gap = " << to_double(duality_gap) 
	 << " ,  minimum_dg = " << to_double(minimum_dg) << endl;

    // Check the invariants. 
    // Note: a) The first should be guaranteed by the termination condition of the electrical flow. 
	  // b) The second should be implied by the first. 
	  // c) The third should be guaranteed by the arc-fixation.
	 
    check_central_path_condition(N, delta, mu);
    check_xasa_bounds(N, delta, mu);
    check_lower_bound_invariant(N, delta, beta, gamma);

    // Update Mu
    mu = (1 - tau) * mu;

    get_resistances(N);
    SpanningTree<Graph<IntegerType, RationalType> , IntegerType, RationalType> ST_original(G0);
    get_spanning_tree_for_electrical_flow_problem(G0, ST_original, N);
    get_initial_currents(N, mu);




    // TODO: Make sure that no removed arcs are in the tree
    // for(auto a: G0.tree_edges){
    //   assert(x[abs(a)] != IntegerType);
    // }

#ifdef VERBOSE
    print_resistances(N);
    print_currents(N);
#endif

    clear_currents(N);

    int number_of_electrical_flow_its = electrical_flow_solver(N, ST_original, mu, delta);

    cout << "number of electrical flow iterations = " << number_of_electrical_flow_its << endl;

    primal_feasibility_check(N, demands);
    dual_feasibility_check(N);

    // print_primal_variable(N);
    // print_dual_variable(N);
    arc_fixation(N, beta, gamma, delta, demands);

    // Recompute duality gap
    duality_gap = update_duality_gap(N);
  }
}