#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<fstream>
#include<limits>
#include "random.h"
#include <boost/rational.hpp>

using namespace std;

//  long long p = 60;
//  long long base =  1LL << p; 
template<typename T, T p = 60, T base = (T(1) << p) >  
class unbounded_integer{ 
//  const static T base =  (1 << 15); //numeric_limits<T>::max() >> 2;
    
 
public:
    vector<T> digits; 
   bool sign;

  unbounded_integer() : sign(false) {}
  
  template< typename NT >
  unbounded_integer( const NT& a ) : sign(a < 0) {
    if( a == 0 ) return; 
    NT b = sign ? -a : a;
    while( b > 0 ) {
      const NT c = b % base;
      digits.push_back( c );
      b /= base;
    }
  }
  
  unbounded_integer( const unbounded_integer<T>& other ) : digits( other.digits ), sign( other.sign ) {}
  unbounded_integer<T>& operator=( const unbounded_integer<T>& other ) {
    if( this != &other ) {
      digits = other.digits;
      sign = other.sign;
    }
    return *this;
  }
  unbounded_integer( unbounded_integer<T>&& other ) : digits( other.digits ), sign( other.sign ) {}
  unbounded_integer<T>& operator=( unbounded_integer<T>&& other ) {
    if( this != &other ) {
      digits = other.digits;
      sign = other.sign;
    }
    return *this;
  }
  
  
  void swap(unbounded_integer<T>& a, unbounded_integer<T>& b)
  {
    unbounded_integer<T> temp = a;
    a = b;
    b = temp;    
  }
  
  T basis() const { return base; }
  
  void add_to_digits_at_position(T& integer_to_be_added, unsigned int position)
  {
    while(true)
    {
      T sum = digits[position] + integer_to_be_added;
      digits[position] = sum % base;
      T carry = (sum - digits[position])/base;
      if(carry !=0 )
      {
	if(digits.size() > position + 1)
	{
	  integer_to_be_added = carry;
	  position = position + 1; 
	  continue;
	}
	else{
	  digits.push_back(carry);
	  break;
	}
      }
      else break;
    }
  }
  
  
  unbounded_integer<T>& plus_assign_to(const unbounded_integer<T>& integer_to_be_added) 
  {
    #ifdef VERBOSE
    cout << "plus_assign_to() " << endl;
    cout << to_double() << "+" << integer_to_be_added.to_double() << endl;
    #endif
    T carry(0);
    auto it1 = digits.begin();
    auto end1 = digits.end();
    auto it2 = integer_to_be_added.digits.begin();
    auto end2 = integer_to_be_added.digits.end();
    for( ; it1 != end1 && it2 != end2; ++it1, ++it2 ) {
      T& a = *it1;
      #ifdef VERBOSE
      cout << "a: " << a << endl;
      #endif
      const T& b = *it2;
      #ifdef VERBOSE
      cout << "b: " << b << endl;
      #endif
      const T sum = carry + a + b;
      #ifdef VERBOSE
      cout << "sum: " << sum << endl;
      #endif
     
      #ifdef VERBOSE 
       const T d = sum % base;
      cout << "d: " << d << endl;
      #endif
      a = sum % base;
      
      carry = ( sum - a ) / base;
      #ifdef VERBOSE
      cout << "carry: " << carry << endl;
      #endif
    }
    
    if( it1 != end1 ) {
      for( ; it1 != end1; ++it1 ) {
	T& a = *it1;
	const T b = 0;
	const T sum = carry + a + b;
	
	#ifdef VERBOSE
	const T d = sum % base;
	cout << "d: " << d << endl;
	#endif
	a = sum % base; 
	carry = ( sum - a ) / base;
	if( carry == 0 ) break;
      }
    }
    
    else {
      
      for( ; it2 != end2; ++it2 ) {
	const T a = 0;
	const T& b = *it2;
	const T sum = carry + a + b;
	const T d = sum % base;
	#ifdef VERBOSE
	cout << "d: " << d << endl;
	#endif
	digits.push_back( d );
	carry = ( sum - d ) / base;
      }
    }
    if( carry > 0 ) {
      digits.push_back( carry );
      #ifdef VERBOSE
      cout << "carry: " << carry << endl;
      #endif
    }    
    
    while( !digits.empty() && digits.back() == T(0) ) {
      digits.pop_back();
    }
    
    return *this;
  }
  
  
   
  unbounded_integer<T>& minus_assign_to(const unbounded_integer<T>& integer_to_be_subtracted)
  {

    T borrow = 0;
    auto it1 = digits.begin();
    auto end1 = digits.end();
    auto it2 = integer_to_be_subtracted.digits.begin();
    auto end2 = integer_to_be_subtracted.digits.end();      
    
    
    bool we_are_subtracting_from_a_larger_number;
    if(less_than_abs(integer_to_be_subtracted))
    {
      we_are_subtracting_from_a_larger_number = false;
    }
    else{
      we_are_subtracting_from_a_larger_number = true;
    }
    for( ; it1 != end1 && it2 != end2; ++it1, ++it2 ) {
      
      T& a = *it1;
      const T& b = *it2;
      
      #ifdef VERBOSE
      cout << "  " << a << " - " << b << " - " << borrow << endl;
      #endif
      
      T difference;
      if(we_are_subtracting_from_a_larger_number)
      {
	if(a >= b + borrow){
	  difference =  a - b - borrow;
	  borrow = 0;
	} else{
	  difference = (base + a) - b - borrow;
	  borrow = 1;
	}
	
	#ifdef VERBOSE
	cout << "difference: " << difference << endl;
	#endif

	a = difference;
	
      }
      else{
	
	if(b >= a + borrow){
	  difference =  b - a - borrow;
	  borrow = 0;
	} else{
	  difference = (base + b) - a - borrow;
	  borrow = 1;
	}
	
	#ifdef VERBOSE
	cout << "difference: " << difference << endl;
	#endif
	a = difference;
	
      }     
    }
    if( it1 != end1 ) {
      assert(we_are_subtracting_from_a_larger_number);
      for( ; it1 != end1; ++it1 ) {
	T& a = *it1;
	
	T difference;
	if(a >= borrow){
	  #ifdef VERBOSE
	  cout << a << " - " << borrow << endl;
	  #endif
	  difference =  a - borrow;
	  borrow = 0;
	} else{
	  difference = (base + a) - borrow;
	  borrow = 1;
	}
	a = difference;
	
      }
    }
    if( it2 != end2 ) {
      
      for( ; it2 != end2; ++it2 ) {
	const T& b = *it2;
	T difference;
	if(b >= borrow){
	  difference =  b - borrow;
	  borrow = 0;
	} else{
	  difference = (base + b) - borrow;
	  borrow = 1;
	}
	digits.push_back(difference);
      }
    }
    
    while( !digits.empty() && digits.back() == T(0) ) {
      digits.pop_back();
    }
    if(is_zero() ) {
      sign = false;
    }
    
    return *this;
  }
  
  unbounded_integer<T>& operator+=(const unbounded_integer<T>& integer_to_be_added) 
  {
    #ifdef VERBOSE
    cout << to_double() << " + " << integer_to_be_added.to_double() << endl;
    #endif
    if( is_zero() )
    {
      
      *this = integer_to_be_added;
      return *this;
    }
    if( integer_to_be_added.is_zero() )
    {
      
      return *this;
    }
    
    else if(sign == integer_to_be_added.sign)
    {
      plus_assign_to(integer_to_be_added);
      sign = integer_to_be_added.sign;
      return *this;
    }
    
    else
    {
      
      if(greater_than_abs(integer_to_be_added))
      {
	
	#ifdef VERBOSE
	cout << to_double() << " > abs " << integer_to_be_added.to_double() << endl;
	#endif
	minus_assign_to(integer_to_be_added);
	
	#ifdef VERBOSE
	cout << "after minus_assign_to: " << to_double() << endl;
	cout << "sign: " << sign << endl;
	#endif
      }
      
      else if(less_than_abs(integer_to_be_added))
      {
	minus_assign_to(integer_to_be_added);
	sign = integer_to_be_added.sign;
      }
      
      else
      {
	
	*this = unbounded_integer<T>( T(0) );
      }
    }
    
    
    
    return *this;
  }
  
  
  
  unbounded_integer<T> operator+(const unbounded_integer<T>& integer_to_be_added) const
  {
    unbounded_integer<T> result( *this );
    result += integer_to_be_added;
    
    return result;
  }
  

  unbounded_integer<T>& operator-=(const unbounded_integer<T>& integer_to_be_subtracted)
  {
    #ifdef VERBOSE
    cout << to_double() << " - " << integer_to_be_subtracted.to_double() << endl;
    #endif
    if( integer_to_be_subtracted.is_zero() ) {
      return *this;
    }
    
    if( is_zero() )
    {
      digits = integer_to_be_subtracted.digits;
      sign = !integer_to_be_subtracted.sign;
      return  *this;
    }
    
    if(sign == integer_to_be_subtracted.sign)
    {
      if(*this == integer_to_be_subtracted)
      {
	assign_zero();
	return *this;
      }
      else
      {
	if(less_than_abs(integer_to_be_subtracted))
	{
	  minus_assign_to(integer_to_be_subtracted);
	  sign = !integer_to_be_subtracted.sign;
	  return *this;
	}
	else{
	  minus_assign_to(integer_to_be_subtracted);
	  sign = integer_to_be_subtracted.sign;
	  return *this;
	}
      }
    }
    
    if( integer_to_be_subtracted.sign )
    {
      *this = plus_assign_to(integer_to_be_subtracted);
      return *this;
    }
    
    if(sign)
    {
      *this = plus_assign_to(integer_to_be_subtracted);
      sign = true;
      return *this;
    }
    
    return *this;
  }
  
  
  
  unbounded_integer<T> operator-(const unbounded_integer<T>& integer_to_be_subtracted) const
  {
    unbounded_integer<T> result( *this );
    result -= integer_to_be_subtracted;
    return result;
  }
  
  
  
  bool is_zero() const {
    for( auto d : digits ) {
      if( d > 0 ) return false;
    }
    return true;
  }
  
  
  void assign_zero() {
    digits.clear();
    sign = false;
  }
  
  unbounded_integer<T>& operator*=(const unbounded_integer<T>& integer_to_be_multiplied)
  {
    
    
    if( is_zero() )
    {
      return *this;
    }
    if( integer_to_be_multiplied.is_zero() )
    {
      
      assign_zero();
      return *this;
    }

    if(sign == integer_to_be_multiplied.sign)
    {
      sign = false;
    }
    else{
    
      sign = true;
    }
    #ifdef VERBOSE
    cout << "digits ";
    for(unsigned int i = 0; i < digits.size(); i++)
    {
      cout << digits[i];
    }
    cout << endl;
    
    #endif
      const T& sqrt_base = 1 << (p/2); 
      const T& shift_amount = log(sqrt_base)/log(2);

    for(T i = digits.size() - 1; i >=0; --i) 
    {
      T carry(0);
      #ifdef VERBOSE
      cout << " i: " << i << endl;
      #endif
      const T a = digits[i];
      digits[i] = 0;
      
      const T&  x1 = (a >> shift_amount);
      const T& x0 = a - x1*sqrt_base;

#ifdef VERBOSE
      if(x1 > 0)
	{
	  cout << "a: " << a << endl << "x1: " << x1 << endl << "x0: " << x0 << endl << "sqrt_base: " << sqrt_base << endl << endl;
	}
#endif
      
      T shift = 0;
      #ifdef VERBOSE
      cout << "shift: " << shift << endl;
      #endif
      for(auto b: integer_to_be_multiplied.digits)
      {

#ifdef VERBOSE
	cout << "SQRT: " << sqrt(base) << endl;
#endif
	const T& shift_amount = log(sqrt_base)/log(2);
	const T& y1 = b >> shift_amount;
	const T& y0 = b - y1*sqrt_base; 

#ifdef VERBOSE
	if(y1 > 0)
	{
	  cout << "b: " << b << endl << "y1: " << y1 << endl << "y0: " << y0 << endl << "sqrt_base: " << sqrt_base << endl << endl;
	}
#endif

#ifdef VERBOSE
	cout << "a * b = " << a * b << endl;
#endif
	const T& digit0 = x0*y0;
		
	const T& digit1 = x1*y0 + x0*y1;
	
	const T& a1 = digit1 >> shift_amount;
	const T& a0 = digit1 - a1*sqrt_base; 
	
	const T& digit2 = x1*y1;
	
	const T& digit1_val = a0*sqrt_base;
	
	T d = (digit0 + digit1_val)%base;
	carry = (digit0 + digit1_val - d) / base;
	carry += digit2;
	carry += a1;
	
	T digit_size = digits.size() - 1;
	
	if( (i + shift) > digit_size)
	{
	  digits.push_back(d);
	}
	else{
	  add_to_digits_at_position(d, i + shift);  
	}
	
	
	#ifdef VERBOSE
	cout << i + shift << " < - > " << digits.size() << endl;
	#endif
	if(carry != 0)
	{
	  digit_size = digits.size() - 1;
	  if( i + shift >= digit_size)
	  {
	    digits.push_back(carry);
	  }
	  else
	  {
	    add_to_digits_at_position(carry, i + shift+1);  
	  }
	  #ifdef VERBOSE
	  cout << "digits[" << shift + 1 << "] = " << digits[shift+1] << endl;
	  #endif
	}
	
	shift++;
	
	#ifdef VERBOSE
	cout << a << " * " << b << endl;
	#endif
      }
      
      #ifdef VERBOSE
      cout << "digits: ";
      for(auto a: digits)
      {
	cout << a;
      }
      cout << endl;
      #endif
      
      
      #ifdef VERBOSE
      cout << "digits after: ";
      for(unsigned int i = 0; i < digits.size(); i++)
      {
	cout << digits[i];
      }
      cout << endl;
      #endif
    }
    
   while( !digits.empty() && digits.back() == T(0) ) {
      digits.pop_back();
    }
    
    return *this;
  }
  

  unbounded_integer<T> operator*(const unbounded_integer<T>& integer_to_be_multiplied) const
  {
    #ifdef VERBOSE
    cout << "we are multiplying: " << endl;
    cout << to_double() << " * " << integer_to_be_multiplied.to_double() << endl;
    long double checking = to_double() *  integer_to_be_multiplied.to_double() ;
    
    for(unsigned int i = 0; i < digits.size(); i++)
    {
      cout << digits[i];
    }
    cout << endl;
    
    for(unsigned int i = 0; i < integer_to_be_multiplied.digits.size(); i++)
    {
      cout << integer_to_be_multiplied.digits[i];
    }
    cout << endl;
    
    #endif
    
    
    unbounded_integer<T> final_result( integer_to_be_multiplied );
    
    final_result *= *this;
    
    
    
    if(sign == integer_to_be_multiplied.sign)
    {
      final_result.sign = false;
    }
    else{
      final_result.sign = true;
    }
    
    #ifdef VERBOSE
    cout << "Answer: " << final_result.to_double() << endl;
    
    assert(final_result.to_double() == checking);
    #endif
    return final_result;
    
  }
  
  unbounded_integer<T> operator-() const
  {
    if( is_zero() ) return unbounded_integer<T>();
    unbounded_integer<T> result( *this );
    result.sign = !sign;
    #ifdef VERBOSE
    cout << "- " << to_double() << " = " << result.to_double() << endl;
    #endif
    return result;
  }
  
  unbounded_integer<T> operator<<(unsigned int shift) const
  {
    if( is_zero() ) return *this;
    
    unbounded_integer result;
    result.sign = sign;
    result.digits.resize(shift,T(0));
    result.digits.insert(result.digits.end(), digits.begin(),digits.end() );
    
    return result;
  }
  
  
  bool less_than_abs( const unbounded_integer<T>& integer_to_be_compared ) const
  {
    if(digits.size() < integer_to_be_compared.digits.size())
    {
      return true;
    }
    else if(digits.size() > integer_to_be_compared.digits.size())
    {
      return false;
    }
    else{
      
      for(unsigned int i = digits.size(); i > 0; --i )
      {
	if(digits[i-1] < integer_to_be_compared.digits[i-1])
	{
	  return true;
	} else  if(digits[i-1] > integer_to_be_compared.digits[i-1])
	{
	  return false;
	}
      }
      
      return false;
    }
  }
  
  
  bool operator<( const unbounded_integer<T>& integer_to_be_compared ) const
  {
    if( sign && !integer_to_be_compared.sign ) {
      return true;
    }
    if( !sign && integer_to_be_compared.sign ) {
      return false;
    }
      
      
      if(sign && integer_to_be_compared.sign)
      {
	if(digits.size() < integer_to_be_compared.digits.size())
	{
	  return false;
	}
	else if(digits.size() > integer_to_be_compared.digits.size())
	{
	  return true;
	}
	else{
	  
	  for(unsigned int i = digits.size(); i > 0; --i )
	  {
	    if(digits[i-1] < integer_to_be_compared.digits[i-1])
	    {
	      return false;
	    } else  if(digits[i-1] > integer_to_be_compared.digits[i-1])
	    {
	      return true;
	    }
	  }
	  
	  return true;
	}
	
      }
      
      else{
	if(digits.size() < integer_to_be_compared.digits.size())
	{
	  return true;
	}
	else if(digits.size() > integer_to_be_compared.digits.size())
	{
	  return false;
	}
	else{
	  
	  for(unsigned int i = digits.size(); i > 0; --i )
	  {
	    if(digits[i-1] < integer_to_be_compared.digits[i-1])
	    {
	      return true;
	    } else  if(digits[i-1] > integer_to_be_compared.digits[i-1])
	    {
	      return false;
	    }
	  }
	  
	  return false;
	}
      }
  }
  
  bool operator==( const unbounded_integer<T>& integer_to_be_compared ) const
  {
    if( is_zero() && integer_to_be_compared.is_zero() ) {
      return true;
    }
    if( sign != integer_to_be_compared.sign ) {
      return false;
    }
    
    if(digits.size() != integer_to_be_compared.digits.size())
    {
      return false;
    }
    else{
      
      for(unsigned int i = 0; i < digits.size(); i++)
      {
	if(digits[i] != integer_to_be_compared.digits[i])
	{
	  return false;
	}
      }
      
      return true;
    }
  }
  
  bool operator!=( const unbounded_integer<T>& integer_to_be_compared ) const
  {
    if( is_zero() && integer_to_be_compared.is_zero() ) {
      return false;
    }
    if( sign != integer_to_be_compared.sign ) {
      return true;
    }
    
    if(digits.size() != integer_to_be_compared.digits.size())
    {
      return true;
    }
    else{
      
      for(unsigned int i = 0; i < digits.size(); i++)
      {
	if(digits[i] != integer_to_be_compared.digits[i])
	{
	  return true;
	}
      }
      
      return false;
    }
  }
  
  
  
  bool greater_than_abs( const unbounded_integer<T>& integer_to_be_compared ) const
  {
    if(digits.size() < integer_to_be_compared.digits.size())
    {
      return false;
    }
    else if(digits.size() > integer_to_be_compared.digits.size())
    {
      return true;
    }
    
    for(unsigned int i = digits.size(); i > 0; --i )
    {
      if(digits[i-1] > integer_to_be_compared.digits[i-1])
      {
	return true;
      } else if(digits[i-1] < integer_to_be_compared.digits[i-1])
      {
	return false;
      }
    }
    
    return false;
    
  }
  
  
  bool operator>( const unbounded_integer<T>& integer_to_be_compared ) const
  {
    if( sign && !integer_to_be_compared.sign ) {
      return false;
    }
    if( !sign && integer_to_be_compared.sign ) {
      return true;
    }
      
      if(sign && integer_to_be_compared.sign){
	
	if(digits.size() < integer_to_be_compared.digits.size())
	{     
	  return true;
	}
	else if(digits.size() > integer_to_be_compared.digits.size())
	{
	  return false;
	}
	
	for(unsigned int i = digits.size(); i > 0; --i )
	{
	  if(digits[i-1] > integer_to_be_compared.digits[i-1])
	  {
	    return false;
	  } else if(digits[i-1] < integer_to_be_compared.digits[i-1])
	  {
	    return true;
	  }
	}
	
	return true;
	
      }
       else{    
	if(digits.size() < integer_to_be_compared.digits.size())
	{     
	  return false;
	}
	else if(digits.size() > integer_to_be_compared.digits.size())
	{
	  return true;
	}
	
	for(unsigned int i = digits.size(); i > 0; --i )
	{
	  if(digits[i-1] > integer_to_be_compared.digits[i-1])
	  {
	    return true;
	  } else if(digits[i-1] < integer_to_be_compared.digits[i-1])
	  {
	    return false;
	  }
	}
	
	return false;
      }
  }
  
  long double to_double() const {
    long double result = 0.0;
    long double b = 1;
    for( auto d : digits ) {
      result += d*b;
      b *= base;
    }
    return sign ? -result : result;
  }
};

template< typename T, T base >
long double to_double( const unbounded_integer<T,base>& i ) {
  return i.to_double();
}


long double to_double(const long double & a){

  return a;
}

long double logarithm_base(long long base1, long double number)
{
  return (log(number)/log(base1));
}


template< typename T, T p >
void div_by_2 (unbounded_integer<T,p>& a)
{
  long long digit_size = a.digits.size();
  
  for(long long i = 0; i < digit_size-1; i++)
  {
    a.digits[i] >>= 1;
    a.digits[i] |= (a.digits[i+1] & T(1)) << (p-1);
  }
  
  a.digits[digit_size-1] >>= 1;
}

unbounded_integer<long long> convert_rounded_ratio_to_integer(long double numerator, long double denominator)
{
  assert( isfinite( numerator ) );
  assert( isfinite( denominator ) );
    long double ratio = numerator/denominator;
    unbounded_integer<long long> answer;
    if(numerator < denominator)
    {
      answer = 0;
    }
    else if(numerator == denominator)
    {
      answer = 1;
    }
    else if(numerator > denominator)
    {
      unbounded_integer<long long> answer_upper_limit = 1;
      unbounded_integer<long long> answer_lower_limit = 1;
#ifdef VERBOSE
      cout << "numerator: " << numerator << endl;
      cout << "denominator: " << denominator << endl;
#endif
      while( to_double(answer_upper_limit) < ratio)
      {
	answer_lower_limit = answer_upper_limit;
	answer_upper_limit*= 2;
      }
#ifdef VERBOSE
      cout << "answer_lower_limit: " << to_double(answer_lower_limit) << endl;
      cout << "answer_upper_limit: " << to_double(answer_upper_limit) << endl;
      cout << "num/den: " << numerator/denominator << endl;
#endif
      answer = answer_lower_limit;
#ifdef VERBOSE
      cout << to_double(answer) << " ----> " << to_double(answer_upper_limit) << endl;
#endif
      while( to_double(answer_upper_limit) - ratio >= 1 )
      {
	answer = answer_upper_limit + answer_lower_limit;
	
	div_by_2( answer );
	
#ifdef VERBOSE
	cout << "division of : " << to_double(answer_upper_limit + answer_lower_limit) << " = " << to_double(answer) << endl;
#endif
	if( to_double(answer) >= ratio )
	{
	  answer_upper_limit = answer;
	}
	if(to_double(answer) <= ratio )
	{
	  answer_lower_limit = answer;
	}
      }
      answer = answer_upper_limit;
    }
    
    return answer;
}
  
template< typename T>
double to_double( const T & d ) {

   
  return static_cast<long double>( d );
}





template< typename IntegerType, typename RationalType >
void round_flow( const RationalType& phi, const RationalType& battery, const RationalType& resistance, IntegerType& flow ) {
  
  const long long f_check = static_cast<long long>( floor( ( battery / ((phi)*( resistance )) ) + RationalType(1)/RationalType(2) ) );
  long double u_f = battery / ((phi)*( resistance ));  
 
  if(u_f + RationalType(2) < numeric_limits<long long>::max() && u_f > numeric_limits<long long>::min())
  {
    flow = f_check; 
  }
  else
  {
    unbounded_integer<long long> f_int(0);
    if(u_f > 0 )
    {
      RationalType denominator = phi*( resistance );
      f_int = convert_rounded_ratio_to_integer( battery, denominator); 
    }
    else if(u_f < 0)
    {
      RationalType denominator = phi*( resistance );
      f_int = convert_rounded_ratio_to_integer(abs( battery ), denominator);; 
      f_int *= -1;
    }
    
    if(abs(u_f - to_double(f_int)) > RationalType(1)/RationalType(2))
    {
      if(u_f > to_double(f_int))
      {
        f_int += 1;
      }
      else{
        f_int -= 1;
      }
    }
    
    if(abs(u_f - to_double(f_int)) > RationalType(1)/RationalType(2))
    {
      cout << "u_f: " << u_f << "	" << battery << " / (( " << phi << " )*( " << resistance << " )) " << endl;
      cout << "f_int: " << to_double(f_int) << endl;
      cout << "fails by: " << abs(u_f - to_double(f_int)) << endl;
    }
      
    assert(abs(u_f - to_double(f_int)) <= RationalType(1)/RationalType(2));
    
    flow = f_int;
  }
}


template< typename IntegerType, typename RationalType >
void round_flow( const RationalType& phi, const RationalType& battery, const IntegerType& resistance, IntegerType& flow ) {
  const long long f_check = static_cast<long long>( floor( ( battery / ((phi)*to_double( resistance )) ) + RationalType(1)/RationalType(2) ) );
  long double u_f = battery / ((phi)*to_double( resistance ));  
  if(u_f + RationalType(2) < numeric_limits<long long>::max() && u_f > numeric_limits<long long>::min())
  {

    flow = f_check; 
  }
  else
  {
    unbounded_integer<long long> f_int(0);
    if(u_f > 0 )
    {
      RationalType denominator = phi*to_double( resistance );
      f_int = convert_rounded_ratio_to_integer( battery, denominator); 
    }
    else if(u_f < 0)
    {
      RationalType denominator = phi*to_double( resistance );
      f_int = convert_rounded_ratio_to_integer(abs( battery ), denominator);; 
      f_int *= -1;
    }
    
    if(abs(u_f - to_double(f_int)) > RationalType(1)/RationalType(2))
    {
      if(u_f > to_double(f_int))
      {
        f_int += 1;
      }
      else{
        f_int -= 1;
      }
    }

    #if !defined(NDEBUG)
//    if(abs(u_f - to_double(f_int)) > RationalType(1)/RationalType(2))
    {
      cout << "u_f: " << u_f << endl;
      cout << "f_int: " << to_double(f_int) << endl;
      cout << "fails by: " << abs(u_f - to_double(f_int)) << endl;

    }
    #endif
      
    assert(abs(u_f - to_double(f_int)) <= RationalType(1)/RationalType(2));
    
    flow = f_int;
  }
}

void multiply_by_bit_shift(long double& f_a_dummy, const long double& r_a)
{
  f_a_dummy = f_a_dummy * r_a;
}

template<typename IntegerType>
void multiply_by_bit_shift(IntegerType& f_a_dummy, const IntegerType& r_a)
{
#ifdef VERBOSE
    cout << "computing this by bit-shift" << endl;
    cout << "f_a = " << to_double(f_a_dummy) << endl;
    
    cout << "r_a = " << to_double(r_a) << endl;
#endif
    long long shift_amount = log(to_double(r_a))/log(2) ;
#ifdef VERBOSE
    cout << "shift_amount = " << shift_amount << endl;
#endif
    long long basis = f_a_dummy.basis();
#ifdef VERBOSE
    cout << "base = " << basis << endl;
#endif
    long long expo = log(basis)/log(2);
#ifdef VERBOSE
    cout << "expo = " << expo << endl;
#endif
    if(shift_amount >= expo)
    {
      unsigned int zeros_to_be_added = shift_amount/expo;
      for(unsigned int i=1; i<= zeros_to_be_added; i++)
      {

      f_a_dummy.digits.insert(f_a_dummy.digits.begin(), 0);
      }
      shift_amount %= expo;
    }
    if(shift_amount < expo && f_a_dummy != IntegerType(0))
    {
      long long digit_back = f_a_dummy.digits.back();
      
      digit_back <<= shift_amount;
      if(digit_back >= (1 << (expo )) || digit_back <= 0)
      {
#ifdef VERBOSE
	cout << "first case" << " " << digit_back << endl;
#endif
	f_a_dummy.digits.push_back(0);
#ifdef VERBOSE
	cout << "max = " <<  numeric_limits<long long>::max() << endl;
	cout << "check = " << (numeric_limits<long long>::max() >> 3) << endl;
#endif
	long long mask =(numeric_limits<long long>::max() >> 3);
#ifdef VERBOSE
	cout << "mask = " << mask << endl;
#endif

	
	for(int i = f_a_dummy.digits.size() -1; i > 0; i--)
	{


	  f_a_dummy.digits[i] |= (( f_a_dummy.digits[i-1]) >> (expo - shift_amount));
	  f_a_dummy.digits[i-1] &= (mask >> shift_amount);
	  f_a_dummy.digits[i-1] <<= shift_amount;
	 

	}

      }
      else{
#ifdef VERBOSE
	cout << "else " << endl;
#endif
	long long mask =(numeric_limits<long long>::max() >> 3);
	for(int i = f_a_dummy.digits.size() -1; i > 0; i--)
	 {
	
	  f_a_dummy.digits[i] <<= shift_amount;
	  	  f_a_dummy.digits[i] |= (( f_a_dummy.digits[i-1]) >> (expo - shift_amount));
	  f_a_dummy.digits[i-1] &= (mask >> shift_amount);
	 }
	 f_a_dummy.digits[0] <<= shift_amount;
      }
      
      if(f_a_dummy.digits.back() == abs(0)) f_a_dummy.digits.pop_back();      
    }
}


/** \brief rounds the resistances and batteries
 * 
 * @param G The graph G on which the min-cost-flow algorithm is being applied
 * @param rho The parameter
 * @param gamma The parameter
 * 
 */
template<typename IntegerType, typename RationalType>
void round_resistance( const RationalType& rho, const RationalType& x, IntegerType& r )
{

     unbounded_integer<long long> answer_upper_limit = 1;

      while( to_double(answer_upper_limit) < x/rho)
      {
#ifdef VERBOSE
	cout << "before = " << answer_upper_limit.digits.back() << endl;
#endif
	answer_upper_limit.digits.back() <<= 1;
#ifdef VERBOSE
	cout << "after = " <<  answer_upper_limit.digits.back() << endl;
#endif
	if(answer_upper_limit.digits.back() >= answer_upper_limit.basis() )
	{
	  answer_upper_limit.digits.back() = 0;
#ifdef VERBOSE
	  cout << "here" << endl;
#endif
	  answer_upper_limit.digits.push_back(1);
#ifdef VERBOSE
	  cout << "after = " <<  answer_upper_limit.digits[answer_upper_limit.digits.size() - 1] << endl;
#endif
	}
#ifdef VERBOSE
	cout << to_double(answer_upper_limit) << " " << x/rho << endl;
#endif
      }
        
    r = answer_upper_limit; 
  
  #ifdef VERBOSE
  cout << "Resistances after rounding: " << to_double(r) << endl; 
  #endif
  
}
template< typename RationalType, typename IntegerType >
void get_phi_factor( const RationalType& phi_old, const RationalType& phi, IntegerType& factor )
{
  factor = static_cast<long long>(floor ( (phi_old/phi) + RationalType(1)/RationalType(2) ) );
}


unbounded_integer<long long> abs(unbounded_integer<long long> integer_whose_absolute_value_needed)
  {
    if(integer_whose_absolute_value_needed < 0)
    {
      integer_whose_absolute_value_needed = -integer_whose_absolute_value_needed;
    }
    
    return integer_whose_absolute_value_needed;
  }


template<typename RationalType>
unbounded_integer<long long> round_alpha( const RationalType& alpha_unrounded, const unbounded_integer<long long>& Delta, const unbounded_integer<long long>& R_a)
{
  
      unbounded_integer<long long> answer_upper_limit = 1;
      while( to_double(answer_upper_limit) < abs(to_double(Delta)/to_double(R_a)))
      {
#ifdef VERBOSE
	cout << "before = " << answer_upper_limit.digits.back() << endl;
#endif
	answer_upper_limit.digits.back() <<= 1;
#ifdef VERBOSE
	cout << "after = " <<  answer_upper_limit.digits.back() << endl;
#endif
	if(answer_upper_limit.digits.back() >= answer_upper_limit.basis() )
	{
	  answer_upper_limit.digits.back() = 0;
#ifdef VERBOSE
	  cout << "here" << endl;
#endif
	  answer_upper_limit.digits.push_back(1);
#ifdef VERBOSE
	  cout << "after = " <<  answer_upper_limit.digits[answer_upper_limit.digits.size() - 1] << endl;
#endif
	}
      }
      
  
    unbounded_integer<long long> alpha = answer_upper_limit;
 
    if(Delta > 0) alpha = -alpha;
    
  return alpha;
 
}