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

Random rg;

template< typename IntegerType, typename RationalType >
void update_flow(const std::vector<int> &cyc, Graph<IntegerType,RationalType> &G, const RationalType factor) 
{
  
  for( auto y : cyc )
  {
    if(y<0)
    {
      G.flows[-y] += factor;    // This is for the case when the edge is in the opposite direction - Therefore, the negative sign
    }
    else
    {	
      G.flows[y] -= factor;
    }	
  }
}


template< typename IntegerType, typename RationalType >
void update_potentials(const int ed, Graph<IntegerType,RationalType> &G, const RationalType factor)
{
  std::vector<int> cyc = G.cycle[ed];	
  for(auto x: cyc ) //temp_edges_cycle)
  {
    
    const int sign = x > 0 ? 1 : -1;
    const int e = abs(x);
    const RationalType delta = sign*G.resistances[e]*factor;
    
    #ifdef VERBOSE
    cout << "edge: " << x;
    cout << " delta: " << delta << endl;
    #endif
    for(auto z: G.cycles_corr_edges[e])
    {
      const int s = z > 0 ? 1 : -1;
      const int y = abs(z)-1;
      #ifdef VERBOSE	
      cout << "other cycles: " << z << " " << G.potential_accross_cycles[y];
      #endif
      if( y != ed ) {
	G.potential_accross_cycles[y] -= s*delta;//calculate_potential_drop(G.cycle[y], G);
      }
      #ifdef VERBOSE
      cout << " " << G.potential_accross_cycles[y] << endl;
      #endif	  
    }
    
  }
  
}


template< typename IntegerType, typename RationalType >
RationalType calculate_energy(Graph<IntegerType,RationalType> &G)
{
  RationalType energy = 0;
  for(unsigned int i=0; i<G.tree_edges.size(); i++)
  {
    int e = G.tree_edges[i];
    energy+= G.flows[e]*G.flows[e]*G.resistances[e] - 2*G.flows[e]*G.voltages[e];	
    
  }
  
  for(unsigned int i=0; i<G.non_tree_edges.size(); i++)
  {
    int e = G.non_tree_edges[i];
    energy+= G.flows[e]*G.flows[e]*G.resistances[e] - 2*G.flows[e]*G.voltages[e];	
  }  
  return energy;
}


template< typename IntegerType, typename RationalType >
void calculate_tree_induced_voltages(Graph<IntegerType,RationalType> &G)
{
  node root =0;
  for(unsigned int i =1; i<=G.no_of_verticies; i++)
  {
    if(G.parents[i] == 0)
    {        
      root = i;
      break;
    }
  }
  
  
  G.tree_induced_voltages.reserve(G.no_of_verticies+1);
  //int c =0;
  for(unsigned int i =1; i<=G.no_of_verticies; i++)
  {
    if(i!=root)
    {
      std::vector<int> path = G.get_path(i, root); 
      RationalType tree_ind_vol = 0.0;
      RationalType flow =0.0;
      RationalType voltage =0.0; 
      RationalType resistance = 0.0;
      for( auto edge : path)
      {
	if(edge<0)
	{
	  int rev_edge = -edge;
	  flow = -G.flows[rev_edge];
	  voltage = -G.voltages[rev_edge];
	  resistance = G.resistances[rev_edge];	
	}
	
	if(edge>0)
	{
	  flow = G.flows[edge];
	  voltage = G.voltages[edge];
	  resistance = G.resistances[edge];	
	}        
	tree_ind_vol += flow*resistance - voltage;        
      }
      
      G.tree_induced_voltages[i] = tree_ind_vol;
    }
    else
    {
      G.tree_induced_voltages[i]=0; 
    }  
  } 
  
  
} 

template< typename IntegerType, typename RationalType >
void calculate_tree_induced_voltages_without_batteries(Graph<IntegerType,RationalType> &G)
{
  node root =0;
  for(unsigned int i =1; i<=G.no_of_verticies; i++)
  {
    if(G.parents[i] == 0)
    {        
      root = i;
      break;
    }
  }
    
  G.tree_induced_voltages.resize(G.no_of_verticies+1);
  //int c =0;
  for(unsigned int i =1; i<=G.no_of_verticies; i++)
  {
    if(i!=root)
    {
      std::vector<int> path = G.get_path(i, root); 
      RationalType tree_ind_vol = 0.0;
      RationalType flow =0.0;
      RationalType resistance = 0.0;
      for( auto edge : path)
      {
	if(edge<0)
	{
	  int rev_edge = -edge;
	  flow = -G.flows[rev_edge];
	  resistance = G.resistances[rev_edge];	
	}
	
	if(edge>0)
	{
	  flow = G.flows[edge];
	  resistance = G.resistances[edge];	
	}        
	tree_ind_vol += flow*resistance;        
      }
      
      G.tree_induced_voltages[i] = tree_ind_vol;
    }
    else
    {
      G.tree_induced_voltages[i]=0; 
    }  
  } 
  
  
} 


template< typename IntegerType, typename RationalType >
RationalType calculate_potential_drop(const std::vector<int> &cycle, Graph<IntegerType,RationalType> &G)
{
  RationalType potential = 0.0;
  RationalType voltage = 0.0;
  RationalType resistance = 0.0;
  RationalType flow = 0;
  for( auto edge : cycle)
  {
    
    if(edge<0)
    {
      int rev_edge = -edge;
      flow = -G.flows[rev_edge];
      voltage = -G.voltages[rev_edge];
      resistance = G.resistances[rev_edge];	
    }
    
    if(edge>0)
    {
      flow = G.flows[edge];
      voltage = G.voltages[edge];
      resistance = G.resistances[edge];	
    }
    
    potential += flow*resistance - voltage;
    
  }
  
  return potential;
}


template< typename IntegerType, typename RationalType >
RationalType calculate_gap(Graph<IntegerType,RationalType> &G)
{
  RationalType gap =0;
  for( auto a : G.non_tree_edges )
  {
    assert( a > 0 );
    RationalType delta = G.tree_induced_voltages[G.tails[a]] - G.tree_induced_voltages[G.heads[a]];
    //cout<<"Potential: "<<potential<<std::endl;
    // cout<<i<<":	 "<<G.potential_accross_cycles[i]<<std::endl;
    const RationalType r_a = G.resistances[a];
    gap += delta*delta/r_a;
	cout << a << ": " << delta*delta << "^2/" << r_a << " = " << delta*delta/r_a << endl;
  }
  cout << "gap: " << gap << endl;
  return gap;
}


template< typename IntegerType, typename RationalType >
int simple_solver_AK(Graph<IntegerType,RationalType> &G, double delta )
{	

  std::vector<double> probability_distribution = G.get_distribution();
  boost::mt19937 gen; 
  gen.seed(static_cast<unsigned int>( std::time(0)) );		
  boost::random::discrete_distribution<> dist(probability_distribution);  
  
  int no_of_iterations=0;
  
  for(unsigned int i=0; i<G.non_tree_edges.size(); i++)
  {
    RationalType R_e = 0.0;
    for(auto y: G.cycle[i])
    {
      R_e += G.resistances[abs(y)];
    }
    G.resistances_accross_cycles[i] = R_e;
  }
  
  calculate_tree_induced_voltages_without_batteries(G);
  RationalType gap = calculate_gap( G );
  
  cout << "starting gap: " << gap << endl;
  
  while(gap>=delta)
    //while(energy-dual_energy>eps)
  {
    
    no_of_iterations++;
    
    if(G.non_tree_edges.size()==0)
    {
      cout<<std::endl<<"There are no cycles in the graph!"<<std::endl;
      break;
    }
    
    int ed = dist(gen);
    
    RationalType R_e=0.0;		
    std::vector<int>& cyc = G.cycle[ed];
    R_e = G.resistances_accross_cycles[ed];
    const arc a = G.non_tree_edges[ed];
    const RationalType potential = G.tree_induced_voltages[G.tails[a]] - G.tree_induced_voltages[G.heads[a]];
    
    RationalType factor = -potential/R_e; // "factor" is sum(del-C_e)/R_e - the flow that has to be subracted!
    // cout<<ed<<": "<<potential<<"	"<<R_e<<" "<<"	"<<factor<<std::endl;
    #ifdef VERBOSE
    cout << "cycle id: " << ed+1 << endl;
    cout << "potential drop before : " << fabs( calculate_potential_drop( cyc, G ) ) << " " << G.potential_accross_cycles[ed] << endl;
    #endif
    
    
    /* Here we update the flows in the cycle 
     * For every edge we store the flow in the direction of the directed edge
     * An edge having negative flow would mean that there is flow in the opposite direction */
    
    /* The function updates the flow in graph 'G', by pushing a flow of 'factor' in the cycle 'c' */
    
    
    update_flow(cyc, G, factor);
    calculate_tree_induced_voltages_without_batteries(G);
    gap = calculate_gap( G );
    
    RationalType gap_check = 0;
    for(unsigned int i=1; i<=G.count; i++)
      {
	double delta = G.tree_induced_voltages[G.tails[i]] - G.tree_induced_voltages[G.heads[i]];
	gap_check += G.flows[i]*G.resistances[i]*G.flows[i] -2*G.voltages[i]*delta/G.resistances[i] + delta*delta/G.resistances[i];
	cout << i << ": " << gap_check << endl;
      }

//     if(no_of_iterations%G.count == 0)
//     {
//       }	
//       
//       std::vector<double> chi( G.no_of_verticies+1,0 );
//       
//       for(unsigned int i =1; i<=G.no_of_verticies; i++)
//       {
// 	double chi_i = 0.0;
// 	for( auto edge : G.incident_edges[i] )
// 	{
// 	  
// 	  if(edge > 0)
// 	  {
// 	    chi_i -= g_prime[edge];
// 	  }
// 	  else{
// 	    chi_i += g_prime[-edge];
// 	  }
// 	  
// 	}
// 	cout << i << ". chi_i: " << chi_i << endl;
// 	chi[i] = chi_i;
//       }
//       
//       double piTchi = 0;
//       
//       cout << "piTchi" << endl << endl;
//       
//       for(unsigned int i=1; i<=G.no_of_verticies; i++)
//       {
// 	cout << i << endl;
// 	piTchi+= chi[i]*G.tree_induced_voltages[i];
// 	cout<< "after the assi" << endl;
// 	cout << piTchi << chi[i] << G.tree_induced_voltages[i] << endl;
//       }
//       std::vector<double> s_hat_prime(G.count+1);
//       for(unsigned int i=1; i<=G.count; i++)
//       {
// 	double delta = G.tree_induced_voltages[G.tails[i]] - G.tree_induced_voltages[G.heads[i]];
// 	s_hat_prime[i] = delta/sqrt(G.resistances[i]);
// 	gap += G.flows[i]*G.resistances[i]*G.flows[i] - 2*g_prime[i]*s_hat_prime[i] + delta*delta/G.resistances[i];
// 	cout << i << ": " << gap << endl;
//       }
//       
//       for(unsigned int i=0; i<G.tree_edges.size(); i++)
//       {
// 	int e = G.tree_edges[i];
// 	energy+= G.flows[e]*G.flows[e]*G.resistances[e] - 2*G.flows[e]*G.voltages[e];	      
//       } 
//       
//       if( energy - gap > dual_energy ) {
// 	dual_energy = energy - gap;
//       }
//       file_object<<std::scientific;
//       file_object<<gap<<" " << energy << " " << dual_energy <<std::endl;	
      
      //gap -= potential*potential/G.resistances[G.non_tree_edges[ed]]; //calculate_gap(G);
//       #ifdef VERBOSE
//       const T gap1 = calculate_gap( G );
//       cout << "gap: " << gap << " - " << gap1 << " = " << gap - gap1 << ", " << energy - dual_energy << endl;
//       #endif
//      assert( fabs( gap - calculate_gap( G ) ) < 1e-2 );
//    }       
  } 
  cout << "electrical flow gap: " << gap << endl;
  
  return no_of_iterations;
}


template< typename IntegerType, typename RationalType >
int simple_solver(Graph<IntegerType,RationalType> &G, ofstream &file_obj, int stretch, int no_of_iter, float beta, std::vector<double>& g_prime)
{	
  RationalType gap = 0;
  RationalType energy = 0;
  RationalType dual_energy = 0;
  double eps = 0.125;
  std::vector<double> probability_distribution = G.get_distribution();
  
  boost::mt19937 gen; 
  gen.seed(static_cast<unsigned int>( std::time(0)) );		
  boost::random::discrete_distribution<> dist(probability_distribution);  
  
  ofstream file_object;
  
  
  int no_of_iterations=0;
  
  G.potential_accross_cycles.reserve(G.non_tree_edges.size()); //The array potential_drop_accross is initialized!
  /* initializing the vector potential_accross_cycles */  
  
  for(unsigned int i=0; i<G.non_tree_edges.size(); i++)
  {
    RationalType R_e = 0.0;
    for(auto y: G.cycle[i])
    {
      R_e += G.resistances[abs(y)];
    }
    G.resistances_accross_cycles[i] = R_e;
  }
  
  for(unsigned int i=0; i<G.non_tree_edges.size(); i++)
  {
    assert( G.non_tree_edges[i] > 0 );
    const RationalType potential = calculate_potential_drop(G.cycle[i], G);
    //cout<<"Potential: "<<potential<<std::endl;
    G.potential_accross_cycles[i] = potential;
    // cout<<i<<":	 "<<G.potential_accross_cycles[i]<<std::endl;
    gap += potential*potential/G.resistances[G.non_tree_edges[i]];
  }
  
  
  dual_energy = energy - gap;
  
  while(gap>=eps)
    //while(energy-dual_energy>eps)
  {
    
    no_of_iterations++;
    
    if(G.non_tree_edges.size()==0)
    {
      cout<<std::endl<<"There are no cycles in the graph!"<<std::endl;
      break;
    }
    
    int ed;
    
    ed = dist(gen);
    
    RationalType R_e=0.0;		
    std::vector<int>& cyc = G.cycle[ed];
    R_e = G.resistances_accross_cycles[ed];
    const RationalType potential = G.potential_accross_cycles[ed]; //calculate_potential_drop(cyc, G);
    
    RationalType factor = potential/R_e; // "factor" is sum(del-C_e)/R_e - the flow that has to be subracted!
    // cout<<ed<<": "<<potential<<"	"<<R_e<<" "<<"	"<<factor<<std::endl;
    #ifdef VERBOSE
    cout << "cycle id: " << ed+1 << endl;
    cout << "potential drop before : " << fabs( calculate_potential_drop( cyc, G ) ) << " " << G.potential_accross_cycles[ed] << endl;
    #endif
    assert( fabs( calculate_potential_drop( cyc, G ) - G.potential_accross_cycles[ed] ) < 1e-1 ); 
    
    
    /* Here we update the flows in the cycle 
     * For every edge we store the flow in the direction of the directed edge
     * An edge having negative flow would mean that there is flow in the opposite direction */
    
    /* The function updates the flow in graph 'G', by pushing a flow of 'factor' in the cycle 'c' */
    
    
    update_flow(cyc, G, factor);
    
    /*updating the potential vector */
    
    update_potentials(ed, G, factor);
    
    
    G.potential_accross_cycles[ed] = 0;
    #ifdef VERBOSE
    cout << "potential drop after: " << calculate_potential_drop( cyc,G ) << endl;
    #endif	  
    assert( fabs( calculate_potential_drop( cyc,G ) ) < 1e-1 );
    
    
    //int div = 500;
    
    
    // if(G.no_of_verticies == 640||G.no_of_verticies == 1280)
    //  {
    //   div =1000;
    //  } 
    
    if(no_of_iterations%G.count == 0)
    {
      gap = 0;
      energy=0;
      for(unsigned int i=0; i<G.non_tree_edges.size(); i++)
      {
	const int e = G.non_tree_edges[i];
	const RationalType r_e = G.resistances[e];
	const RationalType potential = G.potential_accross_cycles[i];
	energy+= G.flows[e]*G.flows[e]*r_e - 2*G.flows[e]*G.voltages[e];
	gap += potential*potential/r_e;
	cout << e << ": " << potential << "^2/" << r_e << " = " << potential*potential/r_e << endl;
      }	
      
      std::vector<RationalType> chi( G.no_of_verticies+1,0 );
      
      for(unsigned int i =1; i<=G.no_of_verticies; i++)
      {
	RationalType chi_i = 0.0;
	for( auto edge : G.incident_edges[i] )
	{
	  
	  if(edge > 0)
	  {
	    chi_i -= g_prime[edge];
	  }
	  else{
	    chi_i += g_prime[-edge];
	  }
	  
	}
	cout << i << ". chi_i: " << chi_i << endl;
	chi[i] = chi_i;
      }
      
      RationalType piTchi = 0;
      //static double gap4 = gap;
     // cout << endl<< "checking the gap: " << gap4 << endl;
      
      cout << "piTchi" << endl << endl;
      
      for(unsigned int i=1; i<=G.no_of_verticies; i++)
      {
	cout << i << endl;
	piTchi+= chi[i]*G.tree_induced_voltages[i];
	cout<< "after the assi" << endl;
	cout << piTchi << chi[i] << G.tree_induced_voltages[i] << endl;
      }
      std::vector<RationalType> s_hat_prime(G.count+1);
      for(unsigned int i=1; i<=G.count; i++)
      {
	RationalType delta = G.tree_induced_voltages[G.tails[i]] - G.tree_induced_voltages[G.heads[i]];
	s_hat_prime[i] = delta/sqrt(G.resistances[i]);
	gap += G.flows[i]*G.resistances[i]*G.flows[i] - 2*g_prime[i]*s_hat_prime[i] + delta*delta/G.resistances[i];
	cout << i << ": " << gap << endl;
      }
      
      for(unsigned int i=0; i<G.tree_edges.size(); i++)
      {
	int e = G.tree_edges[i];
	energy+= G.flows[e]*G.flows[e]*G.resistances[e] - 2*G.flows[e]*G.voltages[e];	      
      } 
      
      if( energy - gap > dual_energy ) {
	dual_energy = energy - gap;
      }
      file_object<<std::scientific;
      file_object<<gap<<" " << energy << " " << dual_energy <<std::endl;	
      
      //gap -= potential*potential/G.resistances[G.non_tree_edges[ed]]; //calculate_gap(G);
      //#ifdef VERBOSE
      const RationalType gap1 = calculate_gap( G );
      cout << "gap: " << gap << " - " << gap1 << " = " << gap - gap1 << ", " << energy - dual_energy << endl;
      //#endif
      assert( fabs( gap - calculate_gap( G ) ) < 1e-2 );
    }       
  } 
  cout << "electrical flow gap: " << gap << endl;
  
  file_object<<std::endl<<std::endl;
  file_object<<"==============================flows  at the end of this run========================================="<<std::endl<<std::endl;
  for(unsigned int j=1; j<=G.count; j++)
  {
    file_object<<G.flows[j]<<",";
  }   
  
  file_object<<std::endl<<std::endl;
  file_object<<"============================== tree induced voltages  at the end of this run========================================="<<std::endl<<std::endl;
  calculate_tree_induced_voltages(G);
  for(unsigned int j=1; j<=G.no_of_verticies; j++)
  {
    file_object<<G.tree_induced_voltages[j]<<",";
  }   
  
  
  return no_of_iterations;
}



template< typename IntegerType, typename RationalType >
RationalType run_experiments(ofstream &file_obj, const int stretch, ofstream &file_obj1, Graph<IntegerType,RationalType> &G, float beta)
{
  
  /* running the simple solver 100 times for the same instance of the spanning tree 
   *  The first column in the file stores the stretch, the next 100 cloumns store the number of iterations in each run                         
   */	 
  int iterations[101];
  RationalType average   = 0.0;
  RationalType variance  = 0.0;
  RationalType std_error = 0.0;
  RationalType tree_condition_number = 0.0;
  
  tree_condition_number = stretch + G.count - 2*(G.no_of_verticies) +2;
  
  for(unsigned int i=1; i<=100; i++)
  {
    for(unsigned int j=1; j<=G.count; j++)
    {
      G.flows[j]=0;
    } 
    int no_of_iterations = 0; //= simple_solver(G, file_obj, stretch, i, beta);
    
    G.the_sanity_check();  
    iterations[i] = no_of_iterations;
    average+=iterations[i];		
  }
  average/=100;
  for(unsigned int i=1; i<=100; i++)
  {
    variance+= (average - iterations[i])*(average - iterations[i])/100;
  }
  std_error = sqrt(variance)/sqrt(100);
  
  iterations[0] = average;
  
  ofstream file_object;
  
  file_obj1<<stretch<<"	"<<tree_condition_number<<"	"<<average<<"	"<<std_error;
  file_obj1<<std::endl;
  
  return average;
}
