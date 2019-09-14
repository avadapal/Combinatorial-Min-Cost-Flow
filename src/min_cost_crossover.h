#include<iostream>
#include<queue>
#include<vector>
#include<algorithm>
#include "graph.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "random.h"



using namespace std;


template< typename IntegerType, typename RationalType >
void compute_s_hat_prime_g_hat_prime_s_prime_g_prime(Graph<IntegerType,RationalType> &G, vector<RationalType>& s_hat_prime, vector<RationalType>& g_hat_prime, vector<RationalType>& s_prime, vector<RationalType>& g_prime, vector<RationalType>& x, vector<RationalType>& s, unsigned int q, RationalType duality_gap)
{
  for(unsigned int i = 1; i <= G.no_of_edges; i++)
    {
      RationalType D = G.voltage_drop(i);
      s_hat_prime[i] = D*x[i];
      s_prime[i] = s[i]*x[i];
      g_prime[i] = q*s_prime[i]/duality_gap - 1;            
    }
 
     RationalType s_prime_sum = 0.0;
      for(unsigned int i = 1; i <= G.no_of_edges; i++)
      {
	s_prime_sum += s_prime[i];
      }
    
     for(unsigned int i = 1; i <= G.no_of_edges; i++)
     {
       g_hat_prime[i] = (q/s_prime_sum)*s_hat_prime[i] - 1;
     }
    
}


template< typename IntegerType, typename RationalType >
RationalType function_lamba(Graph<IntegerType,RationalType>&G, vector<RationalType>&x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v, RationalType lamba)
{
  RationalType old_potential = calculate_potential_function(G, x, s, q);
  
    vector<RationalType> y_bar1(G.no_of_verticies + 1, 0);
    vector<RationalType> s_bar1(G.no_of_edges + 1, 0);  
    vector<RationalType> y_tilde = rounding_scheme( G, v , x, y, s);
    
    for(unsigned int i = 1; i <= G.no_of_verticies; i++)
    {
      y_bar1[i] = lamba*y_tilde[i] + (1-lamba)*y[i];
      
    }
   
   s_bar1 = compute_dual_slacks( G, y_bar1);

    RationalType new_potential1 = calculate_potential_function(G, x, s_bar1, q);
    
    RationalType new_pd = new_potential1 - old_potential;
  
    return new_pd;
}


template< typename IntegerType, typename RationalType >
RationalType function_beta(Graph<IntegerType,RationalType>& G, vector<RationalType>&x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v, RationalType beta)
{
    RationalType old_potential = calculate_potential_function(G,x,s,q);
    vector<RationalType> x_bar1(G.no_of_edges + 1, 0);    
    rounding_scheme( G, v , x, y, s);
    vector<RationalType> x_tilde = find_rounding_scheme_tree_solution(G, v);
    for(unsigned int i = 1; i <= G.no_of_edges; i++)
    {
      x_bar1[i] = beta*x_tilde[i] + (1-beta)*x[i];
    }
    RationalType new_potential1 = calculate_potential_function(G, x_bar1, s, q); 
    RationalType new_pd = new_potential1 - old_potential;
    
    return new_pd;
}

template< typename IntegerType, typename RationalType >
RationalType get_derivative(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v, RationalType lamba)
{
  RationalType fun_val1 = function_lamba(G, x, s, y, q, v, lamba + RationalType(1)/RationalType(100) );
  RationalType fun_val2 = function_lamba(G, x, s, y, q, v, lamba);
  
  RationalType derivative = (fun_val1 - fun_val2)*RationalType(100);
  
  return derivative;
}

template< typename IntegerType, typename RationalType >
RationalType get_derivative_beta(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v, RationalType beta)
{
  RationalType fun_val1 = function_beta(G, x, s, y, q, v, beta + 0.01);
  cout << "fun_val1: " << fun_val1 << endl;
  
  RationalType fun_val2 = function_beta(G, x, s, y, q, v, beta);
  cout << "fun_val2: " << fun_val2 << endl;
  
  RationalType derivative = (fun_val1 - fun_val2)/0.01;
  cout << "derivative: " << derivative << endl;
  
  return derivative;  
}


template< typename IntegerType, typename RationalType >
RationalType get_double_derivative_beta(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v, RationalType beta)
{
  RationalType first_derivative1 = get_derivative_beta(G, x, s,  y, q, v, beta + 0.01);
  RationalType first_derivative2 = get_derivative_beta(G, x, s,  y, q, v, beta);
  
  RationalType second_derivative = (first_derivative1 - first_derivative2)/0.01;
  
  return second_derivative;  
}

template< typename IntegerType, typename RationalType >
RationalType get_double_derivative(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q , node v, RationalType lamba)
{
  RationalType first_derivative1 = get_derivative(G, x, s,  y, q, v, lamba + 0.01);
  RationalType first_derivative2 = get_derivative(G, x, s,  y, q, v, lamba);
  
  RationalType second_derivative = (first_derivative1 - first_derivative2)/0.01;
  
  return second_derivative;
}

template< typename IntegerType, typename RationalType >
RationalType get_optimum_lambda(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v)
{
  cout << endl <<  "get_optimum_lambda " << endl << endl;
  
  RationalType lamba = -3;
  for(unsigned int i = 0; i <= 1000; i++)
  {
    lamba += (0.01);
    cout << lamba << " : " << function_lamba(G, x, s, y, q, v, lamba) << endl;
  }
  
  
  //T lamba = 0.0;
  RationalType xn = 0.0;
  RationalType xn1 = 0.0;
 
  while(true)
  {
    xn1 = xn - (get_derivative(G,  x, s, y, q, v, xn))/(get_double_derivative(G,  x, s, y, q, v, xn));
    
    cout << xn << ":  " << get_derivative(G,  x, s, y, q, v, xn) << " / " << get_double_derivative(G,  x, s, y, q, v, xn) << endl;
    
    xn = xn1;
    
    cout << "derivative: " << get_derivative(G,  x, s, y, q, v, xn1) << endl;
    if(abs(get_derivative(G,  x, s, y, q, v, xn1)) < 0.001) break;
  }
 
 cout << "xn1: " << xn1 << endl << "xn: " << xn << endl;
 
 if(xn1 < 0)
 {
   cout << " the function value: " << function_lamba(G, x, s, y, q, v, xn1) << endl;
   RationalType val = 0.35;
   cout << " the function value at 0.35: " << function_lamba(G, x, s, y, q, v, val) << endl;
 }
 
 assert(xn1 > 0);
 
 return xn1;
}

template< typename IntegerType, typename RationalType >
RationalType get_optimum_beta(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v)
{
//  T lamba = 0.0;
  RationalType xn = 0.0;
  RationalType xn1 = 0.0;
 
  while(true)
  {
    xn1 = xn - (get_derivative_beta(G,  x, s, y, q, v, xn))/(get_double_derivative_beta(G,  x, s, y, q, v, xn));
    
    cout << get_derivative_beta(G,  x, s, y, q, v, xn) << " / " << get_double_derivative_beta(G,  x, s, y, q, v, xn) << endl;
    
    xn = xn1;
    
    if(abs(get_derivative_beta(G,  x, s, y, q, v, xn1)) < 0.001) break;
  }
 
 cout << "xn1: " << xn1 << endl << "xn: " << xn << endl;
 return xn1;
}

template< typename IntegerType, typename RationalType >
RationalType find_parameter_t(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v)
{
  RationalType old_potential = calculate_potential_function(G,x,s,q);
     
    RationalType old_pd = std::numeric_limits<RationalType>::max();
    RationalType new_pd = std::numeric_limits<RationalType>::max();

    RationalType min_lambda = 0;
    RationalType min_potential = 2*old_potential;
     
      for(unsigned int lamba = 1; lamba <= 1000; lamba++)
    {
      RationalType lamba1 = lamba * 0.001;
    
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
//     new_pd = new_potential1 - old_potential;

      new_pd = function_lamba(G, x, s, y, q, v, lamba1);
      
#ifdef VERBOSE
    cout << lamba1 << " - >  potential difference: " << new_potential1 - old_potential << endl;
#endif
     if(new_pd > old_pd || new_pd > 0)
     {
       //if( min_potential > new_potential1 ) {
      //min_potential = new_potential1;
      min_lambda = lamba1;
      break;
    }
    old_pd = new_pd;
    }
    cout << "potential difference: (lambda) " << old_pd << endl;
    return min_lambda;
}

template< typename IntegerType, typename RationalType >
RationalType find_beta(Graph<IntegerType,RationalType>& G, vector<RationalType>& x, vector<RationalType>& s, vector<RationalType>& y, unsigned int q, node v, RationalType min_potential, RationalType old_potential )
{
  RationalType old_pd = std::numeric_limits<RationalType>::max();
  RationalType new_pd = std::numeric_limits<RationalType>::max();
  RationalType min_beta = 0.0;
  for(unsigned int beta = 1; beta <= 1000; beta++)
  {
    RationalType beta1 = beta * 0.001;
    
//     vector<T> x_bar1(G.no_of_edges + 1, 0);
//     rounding_scheme( G, v , x, y, s);
//     vector<T> x_tilde = find_rounding_scheme_tree_solution(G, v);
//     for(unsigned int i = 1; i <= G.no_of_edges; i++)
//     {
//       x_bar1[i] = beta1*x_tilde[i] + (1-beta1)*x[i];
//     }
//     T new_potential1 = calculate_potential_function(G, x_bar1, s, q);
//     
//     new_pd = new_potential1 - old_potential;
    
    new_pd = function_beta(G, x, s, y, q, v, beta1);
    
#ifdef VERBOSE
    cout << beta1 << " - >  potential difference: " << new_pd << endl;
#endif
    if(new_pd > old_pd || new_pd > 0)
    {
     //if( min_potential > new_potential1 ) {
//       min_potential = new_potential1;
      min_beta = beta1;
      break;
    }
    old_pd = new_pd;
  }
  
  cout << "potential difference: (beta) " << old_pd << endl;
  return min_beta;
}