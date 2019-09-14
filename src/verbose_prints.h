

using namespace std;

template<typename Network, typename IntegerType, typename RationalType> 
 void print_queries(Network& N, SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original)
 {
    cout << "Printing Queries: " << endl;
    typedef typename Network::GraphType Graph;
    Graph& G0 = N.G;   
    for(unsigned int v = 1; v <= G0.no_of_vertices; v++)
    {
      cout << v << ": = " << to_double(ST_original.query(v,0)) << endl;
    }
 }
 
 template<typename Network>
 void print_resistances(Network& N)
 {
   cout << "Printing Resistances: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;
   
   cout << "lower , upper , roof " << endl;
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
     cout << a << ": " << to_double(N.arcdata[a].resistance_lower) << " , " << to_double(N.arcdata[a].resistance_upper) << " , "
		       << to_double(N.arcdata[a].resistance_roof) << endl;
   }
 }
 
 template<typename Network>
 void print_currents(Network& N)
 {
  cout << "Printing Currents: " << endl;
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;
  
    cout << "lower , upper , roof " << endl;
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
     cout << a << ": " << to_double(N.arcdata[a].current_lower) << " , " << to_double(N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) 
	       << " , " << to_double(N.arcdata[a].current_roof) << endl;
   }
 }
 
  template<typename Network>
 void print_initial_currents(Network& N)
 {
  cout << "Printing Initial Currents: " << endl;
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;
  
    cout << "lower , upper , roof " << endl;
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
     cout << a << ": " << to_double(N.arcdata[a].initial_current_lower) << " , " << 
			  to_double(N.arcdata[a].cur_src_vw - N.arcdata[a].initial_current_lower) 
	       << " , " << to_double(N.arcdata[a].initial_current_roof) << endl;
   }
 }
 
 template<typename Network>
 void print_original_graph_current(Network& N)
 {
   cout << "Currents: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
     cout << a << " = " << to_double(N.arcdata[a].current) << endl;
   }
   
 }
 
 template<typename Network>
 void print_tree_induced_voltages(Network& N)
 {
   cout << "Printing Tree Induced Voltages: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;
   
   cout << "head, tail, vw " << endl;
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
     cout << a << ": " << to_double(N.nodedata[G0.heads[a]].voltage) << " , " << to_double(N.nodedata[G0.tails[a]].voltage) << ", "
		       << to_double(N.arcdata[a].voltage_vw) << endl;
   }
 }
 
 template<typename Network>
 void print_tree_edges(Network& N)
 {
   cout << "Printing Tree Edges: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;
   
   for(auto a: G0.tree_edges)
   {
     cout << a << " , " ;
   }
   cout << endl;
 }
 
 template<typename Network>
 void print_non_tree_edges(Network& N)
 {
   cout << "Printing Non Tree Edges: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;
   
   for(auto a: G0.non_tree_edges)
   {
     cout << a << " , ";
   }
   cout << endl;
 }
 
 template<typename Network>
 void print_directions(Network& N)
 {
  
   cout << "Printing Directions: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;

   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
    {
      cout << a << ": = " << N.arcdata[a].direction << endl;
    }
 }
 
 template<typename Network>
 void print_primal_variable(Network& N)
 {
   cout << "Printing Primal Variable: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;
   
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
    cout << a << ": " << to_double(N.arcdata[a].xlower) << " , " << to_double(N.arcdata[a].capacity - N.arcdata[a].xlower) << " , "
	      << to_double(N.arcdata[a].infeasibility) << endl;
   }
 }


  template<typename Network>
 void print_dual_variable(Network& N)
 {
   cout << "Printing Dual Variable: " << endl;
   typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;
   cout << "Slacks: " << endl;
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
    cout << a << ": " << to_double(N.arcdata[a].slower) << " , " << to_double(N.arcdata[a].supper) << " , "
	      << to_double(N.arcdata[a].sroof) << endl;
   }
   
   cout << endl << "Potentials: " << endl;
   for(unsigned int v = 1; v <= G0.no_of_vertices; v++)
   {
     cout << "Node " << v << ":" << to_double(N.nodedata[v].potential) << endl;
   }
   for(unsigned int a = 1; a <= G0.no_of_edges; a++)
   {
     cout << "Arc " << a << "'s potential vw: " << to_double(N.arcdata[a].potentialvw) << endl;
   }
 }
 
 #ifdef VERBOSE
	  cout << endl << "######################################### " << endl;
	  
	  RationalType x_hat_norm = 0;
	  for(unsigned int i=1; i <= G0.no_of_edges; i++)
	  {
	    if(N.arcdata[i].infeasibility == RationalType(0)) continue;
	    
	    x_hat_norm += phi* to_double(N.arcdata[i].initial_current_roof - N.arcdata[i].current_roof)/
	         to_double(N.arcdata[i].infeasibility);
	  }
	  
	  cout << "x_hat_norm and rho : " << endl;
	  cout << x_hat_norm << " < - > " << rho << endl;
	  RationalType gTx(0);
	  RationalType x_tilde2(0);
	  RationalType lambda = RationalType(1)/RationalType(4);
	  
	  
	        if(x_hat_norm <= rho){    
	cout << " ||x~|| <= rho case " << endl;	
	{
	  RationalType new_potential = calculate_potential_function_new(N, q);
	  cout << old_potential - new_potential << " <-> " << "( " << lambda << " * " <<  gTx << " / " << rho
	  << ")  - ( " << lambda << " * " << lambda << " * " << x_tilde2 << " )  / ( " << 2 << " * ( " << 1 << " - " << lambda << " ) *  " << rho << " * " << rho << " )  = " <<  
	  (lambda * gTx/rho) - (lambda * lambda*x_tilde2) / (2 * (1 -lambda) * rho *rho ) << endl;
	}
	     cout << " ||x~|| <= rho case " << endl;	
	{
	  RationalType new_potential = calculate_potential_function_new(N, q);
	  cout << old_potential - new_potential << " <-> " << "( " << lambda << " * " <<  gTx << " / " << rho
	  << ")  - ( " << lambda << " * " << lambda << " * " << x_tilde2 << " )  / ( " << 2 << " * ( " << 1 << " - " << lambda << " ) *  " << rho << " * " << rho << " )  = " <<  
	  (lambda * gTx/rho) - (lambda * lambda*x_tilde2) / (2 * (1 -lambda) * rho *rho ) << endl;
	}
	  
		} 
	  cout << " ######################################### " << endl;
	  
	cout << endl << "primal variables after the update:" << endl;
	 for(unsigned int i = 1; i <= G0.no_of_edges; i++)
	 {
	  cout << i << ": " << N.arcdata[i].xlower << " , " << N.arcdata[i].capacity - N.arcdata[i].xlower << " , "
	      << N.arcdata[i].infeasibility << endl;
	 }
#endif