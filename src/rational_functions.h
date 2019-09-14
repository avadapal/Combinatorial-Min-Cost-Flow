#include<iostream>
#include<queue>
#include<vector>
#include<algorithm>
#include "graph.h"
#include<boost/rational.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "random.h"

using namespace std;


// struct LongIntegers{
// 
//   vector<long long>& operator+ (vector<long long>& num1, vector<long long>& num2)
//   {
//     std::vector<long long> answer;
//     long long limit = (std::numeric_limits<long long>::max())/2;
//     long long length = max(num1.size(), num2.size());
//     long long carry = 0;
//     for(long long i = size; i>=0; i--)     
//     {
//       long long digit = num1[i] + num2[i] + carry;
//       if(digit < limit)
//       {
// 	answer.push_back(digit);
//       }
//       else{
// 	answer.push_back(digit%limit);
// 	carry = digit / limit;
//       }
//     }
//     
//     return answer;
//   }
//   
//   vector<long long>& operator* (vector<long long>& num1, vector<long long>& num2)
//   {
//     std::vector<long long> answer;
//     long long limit = (std::numeric_limits<long long>::max())/2;
//     long long length = max(num1.size(), num2.size());
//     long long carry = 0;
//     for(long long i = size; i>=0; i--)     
//     {
//       long long digit = (num1[i] * num2[i]) + carry;
//       if(digit < limit)
//       {
// 	answer.push_back(digit);
//       }
//       else{
// 	answer.push_back(digit%limit);
// 	carry = digit / limit;
//       }
//     }
//     return answer;
//   }
//   
//   vector<long long>& operator- (vector<long long>&num1 , vector<long long>& num2)
//   {
//     std::vector<long long> answer;
//     long long limit = (std::numeric_limits<long long>::max())/2;
//     long long length = max(num1.size() , num2.size());
//     long long borrow = 0;
//     
//     for(long long i = size; i >= 0; i--)
//     {
//       long long digit = (num2[i] + borrow*limit) - num1[i];
//     }
//     
//     return answer;
//   }
// };

template<typename IntType>
boost::rational<IntType> fabs( const boost::rational<IntType>& rational_number ){
   return rational_number >= 0 ? rational_number : -rational_number;
}

template<typename IntType>
long double floor( const boost::rational<IntType>& rational_number ){
  return floor(boost::rational_cast<long double>(rational_number));
}

template<typename IntType>
long double ceil(const boost::rational<IntType>& rational_number){
  return ceil(boost::rational_cast<long double>(rational_number));
}

template<typename IntType>
long double log(const boost::rational<IntType>& rational_number){
  return log(boost::rational_cast<long double>(rational_number));
}

// template<typename IntType>
// long double max(const boost::rational<IntType>& rational_number1,const boost::rational<IntType>& rational_number2 ){
//   return max(boost::rational_cast<long double>(rational_number1), boost::rational_cast<long double>(rational_number2)  );
// }

template<typename IntType>
long double sqrt(const boost::rational<IntType>& rational_number){
  return sqrt(boost::rational_cast<long double>(rational_number));
}

