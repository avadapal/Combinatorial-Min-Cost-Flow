#ifndef RANDOM_H
#define RANDOM_H

#include <fstream>
#include <boost/random.hpp>

struct Random {
  
  std::ifstream randomSource;
  const unsigned char c1=0;
  const unsigned char c2=0;
  const unsigned char c3=0;
  const unsigned char c4=0;
  const unsigned int z1=0;
  const unsigned int z2=0;
  const unsigned long z;
  
  boost::lagged_fibonacci44497 rng;
  
  Random() : randomSource( "/dev/urandom" ),
  c1( randomSource.get() ),
  c2( randomSource.get() ),
  c3( randomSource.get() ),
  c4( randomSource.get() ),
  z1( (c1 << 8) + c2 ),
  z2( (c3 << 8) + c4 ),
  z( ( z1 << 16 ) + z2 ),
    //z(404912711),
  rng( static_cast<uint32_t>( z ) ) {
    randomSource.close();
  }
  
  Random( const unsigned long& seed ) : z( seed ), rng( z ) {}
  
  double operator()() {
    return rng();
  }
  
  unsigned long seed() const {
    return z;
  } 
};


struct ArcSampler {

  std::ifstream randomSource;
  const unsigned char c1=0;
  const unsigned char c2=0;
  const unsigned char c3=0;
  const unsigned char c4=0;
  const unsigned int z1=0;
  const unsigned int z2=0;
  const unsigned long z;
  
  
  boost::mt19937 rng;
  //boost::lagged_fibonacci44497 rng;
  boost::random::discrete_distribution<> dist;
  
  template<typename Network, typename IntegerType>
  ArcSampler( const Network& N, std::vector<IntegerType> sum_to_root) : randomSource( "/dev/urandom" ),
  c1( randomSource.get() ),
  c2( randomSource.get() ),
  c3( randomSource.get() ),
  c4( randomSource.get() ),
  z1( (c1 << 8) + c2 ),
  z2( (c3 << 8) + c4 ),
  z( ( z1 << 16 ) + z2 ),
  rng( static_cast<uint32_t>( z ) ),
  dist( N.get_distribution(sum_to_root) ) {
    randomSource.close();
    //gen.seed(static_cast<unsigned int>(z);
  }
 
  unsigned int operator()() {
    return dist(rng);
  }
  
};


#endif
