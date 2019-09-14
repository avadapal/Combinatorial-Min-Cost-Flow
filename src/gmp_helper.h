#ifndef GMP_HELPER_H
#define GMP_HELPER_H

#include <gmpxx.h>
#include <gmp.h>

#include<vector>


double to_double( const mpz_class & d)
{
  double double_value = mpz_get_d(d.get_mpz_t());   
  return double_value;
} 

void multiply_by_bit_shift(mpz_class& f_a, const mpz_class& r_a)
{
  f_a = f_a * r_a;
}
template< typename RationalType>
mpz_class round_alpha( const RationalType& alpha_unrounded, const mpz_class& Delta, const mpz_class& R_a)
{
  double alpha_intr = floor(alpha_unrounded + 0.5 );

  mpz_class alpha;
  mpz_set_d(alpha.get_mpz_t(), alpha_intr);
    
  return alpha;
}

template<typename RationalType >
void round_resistance( const RationalType& rho, const RationalType& x, mpz_class& r )
{
     mpz_class answer_upper_limit(1);

#ifdef VERBOSE
      cout << "numerator: " << numerator << endl;
      cout << "denominator: " << denominator << endl;
#endif
      while( to_double(answer_upper_limit) < x/rho)
      {
	answer_upper_limit <<= 1;
      }
      
      r = answer_upper_limit;
}

template<typename RationalType>
void round_flow( const RationalType& phi, const RationalType& battery, const mpz_class& resistance, mpz_class& flow )
{
  const RationalType flow_intr = floor( ( battery / ((phi)*to_double( resistance)) ) + RationalType(1)/RationalType(2) );
  mpz_set_d( flow.get_mpz_t(), flow_intr);
}

template<typename RationalType>
void round_flow( const RationalType& phi, const RationalType& battery, const RationalType& resistance, mpz_class& flow )
{
  const RationalType flow_intr = floor( ( battery / ((phi)* ( resistance)) ) + RationalType(1)/RationalType(2) );
  mpz_set_d( flow.get_mpz_t(), flow_intr);

}

template<typename RationalType>
void get_phi_factor( const RationalType& phi_old, const RationalType& phi, mpz_class& factor )
{
  const RationalType phi_factor_intr = floor ( (phi_old/phi) + RationalType(1)/RationalType(2) );
  mpz_init_set_d(factor.get_mpz_t() , phi_factor_intr);
}


#endif