

using namespace std;

/** \brief Computes the tree currents using the updated non-tree edge currents
 *
 * @param ST Spanning Tree
 *
 */

template<typename Network>
void flow_conservation_check(Network& N) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  typedef typename Network::IntegerType IntegerType;
  Graph& G0 = N.G;

  for (unsigned int v = 1; v <= G0.no_of_vertices; v++) {

    IntegerType flow_into_v = 0;
    IntegerType flow_into_v_initial = 0;
    for (auto a : G0.incident_edges[v]) {

      RationalType real_direction = a;

      if (N.arcdata[abs(a)].direction != 0) {
        real_direction = a * N.arcdata[abs(a)].direction;
      } else {
        real_direction = a;
      }
      if (real_direction > 0) {
        flow_into_v -= N.arcdata[abs(a)].current_lower;
        flow_into_v_initial -= N.arcdata[abs(a)].initial_current_lower;
      } else {
        flow_into_v -= (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower);
        flow_into_v_initial -= (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].initial_current_lower);
      }
      if (a > 0) {
        flow_into_v -= N.arcdata[abs(a)].current_roof;
        flow_into_v_initial -=  N.arcdata[abs(a)].initial_current_roof;
      } else {
        flow_into_v += N.arcdata[abs(a)].current_roof;
        flow_into_v_initial +=  N.arcdata[abs(a)].initial_current_roof;
      }
    }

    RationalType difference = to_double(flow_into_v - flow_into_v_initial);
    if (fabs(difference) >= 1e-3) {
      cout << to_double(flow_into_v) << " < - > " << to_double(flow_into_v_initial) << endl;
      cout << "Assert Fails by = " << to_double(flow_into_v - flow_into_v_initial) << endl;
    }
  }
}

template<typename IntegerType, typename RationalType>
void
compute_tree_currents_by_non_tree_edge_DS(SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST )
{
  Graph<IntegerType, RationalType>& G = ST.G;
  vector<IntegerType> imbalance( G.no_of_verticies + 1, 0 );
  vector<IntegerType> carry( G.no_of_verticies + 1, 0 );
  for ( auto aa : G.non_tree_edges ) {

    if (G.original_auxiliary_transformed_corrospondence[aa].is_invalid) continue;

    assert( aa > 0 );
    const IntegerType alpha = G.original_auxiliary_transformed_corrospondence[aa].current;
    imbalance[G.head(aa)] += alpha;
    imbalance[G.tail(aa)] -= alpha;
#ifdef VERBOSE
    cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.original_auxiliary_transformed_corrospondence[aa].current) << " " << to_double(alpha) << endl;
#endif
  }

  for ( arc aa : boost::adaptors::reverse( G.tree_edges ) ) {

    if (G.original_auxiliary_transformed_corrospondence[abs(aa)].is_invalid == true) continue;

    const node v = G.tail(aa);
    const node w = G.head(aa);
    assert( ST.depth[v] > ST.depth[w] );

    if ( aa > 0 ) {
      const IntegerType& imb = imbalance[v];
      const IntegerType outflow = imb ;
      carry[w] += outflow;
      G.original_auxiliary_transformed_corrospondence[aa].current = outflow;

      imbalance[w] += imb;

#ifdef VERBOSE
      const IntegerType alpha = outflow - carry[v];
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.original_auxiliary_transformed_corrospondence[aa].current) << " " << to_double(alpha) << endl;
#endif
    } else {
      const IntegerType& imb = imbalance[v];
      const IntegerType outflow = imb ;

      carry[w] += outflow;
      G.original_auxiliary_transformed_corrospondence[-aa].current = -outflow;
      imbalance[w] += imb;

#ifdef VERBOSE
      const IntegerType alpha = outflow - carry[v];
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << -to_double(G.original_auxiliary_transformed_corrospondence[aa].current[-aa]) << " " << to_double(alpha) << endl;
#endif
    }
  }
#if VERBOSE
  cout << "imbalance of root: " << to_double( imbalance[ST.root] ) << endl;
  cout << "Tree Condition Number: " << ST.calculate_tree_condition_number();
#endif

  assert( RationalType(1e3) * fabs( to_double( imbalance[ST.root] ) )  < RationalType(1) );

}

template<typename Network, typename IntegerType, typename RationalType>
void
compute_tree_currents_by_non_tree_edge_new(Network& N,
    SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST_original
                                          )
{
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;

  //TODO: THRESHOLD_X, THRESHOLD_S needs to be incoorporated.

  vector<IntegerType> imbalance( G0.no_of_vertices + 1, 0 );
  vector<IntegerType> carry( G0.no_of_vertices + 1, 0 );
  for ( auto aa : G0.non_tree_edges ) {

    //if(G.original_auxiliary_transformed_corrospondence[aa].is_invalid) continue;

    assert( aa > 0 );
    const IntegerType alpha = N.arcdata[aa].current;
    imbalance[G0.head(aa)] += alpha;
    imbalance[G0.tail(aa)] -= alpha;
#ifdef VERBOSE
    cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.original_auxiliary_transformed_corrospondence[aa].current) << " " << to_double(alpha) << endl;
#endif
  }

  for ( arc aa : boost::adaptors::reverse( G0.tree_edges ) ) {

    //if(G.original_auxiliary_transformed_corrospondence[abs(aa)].is_invalid == true) continue;

    const node v = G0.tail(aa);
    const node w = G0.head(aa);
    assert( ST_original.depth[v] > ST_original.depth[w] );

    if ( aa > 0 ) {
      const IntegerType& imb = imbalance[v];
      const IntegerType outflow = imb ;
      carry[w] += outflow;
      N.arcdata[aa].current = outflow;

      imbalance[w] += imb;

#ifdef VERBOSE
      const IntegerType alpha = outflow - carry[v];
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << to_double(G.original_auxiliary_transformed_corrospondence[aa].current) << " " << to_double(alpha) << endl;
#endif
    } else {
      const IntegerType& imb = imbalance[v];
      const IntegerType outflow = imb ;

      carry[w] += outflow;
      N.arcdata[-aa].current = -outflow;
      imbalance[w] += imb;

#ifdef VERBOSE
      const IntegerType alpha = outflow - carry[v];
      cout << aa << " =  (" << G.tail(aa) << "," << G.head(aa) << "): " << -to_double(G.original_auxiliary_transformed_corrospondence[aa].current[-aa]) << " " << to_double(alpha) << endl;
#endif
    }
  }
#if VERBOSE
  cout << "imbalance of root: " << to_double( imbalance[ST.root] ) << endl;
  cout << "Tree Condition Number: " << ST.calculate_tree_condition_number();
#endif

  assert( RationalType(1e3) * fabs( to_double( imbalance[ST_original.root] ) )  < RationalType(1) );

}


template<typename Network>
typename Network::IntegerType compute_delta(Network& N, arc a)
{
  typedef typename Network::IntegerType IntegerType;
  
  cout << "compute delta" << endl;
  IntegerType R_a = N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper + N.arcdata[a].resistance_roof;
  IntegerType Delta(0);
  if(N.arcdata[a].resistance_roof < N.arcdata[a].resistance_lower || N.arcdata[a].resistance_roof < N.arcdata[a].resistance_upper)
  {
    if(N.arcdata[a].resistance_lower < N.arcdata[a].resistance_upper)
    {
      if(N.arcdata[a].direction == 1)
      {
	Delta = N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - 
			      N.arcdata[a].resistance_lower * N.arcdata[a].current_lower + 
			      N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;
      }
      else
      {
	Delta = N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - 
			      N.arcdata[a].resistance_lower * N.arcdata[a].current_lower - 
			      N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

      }
    }
    else
    {
      if(N.arcdata[a].direction == 1)
      {
	Delta = N.arcdata[a].resistance_lower * N.arcdata[a].current_lower - 
			    N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - 
                            N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;
      }
      else
      {
	Delta = N.arcdata[a].resistance_lower * N.arcdata[a].current_lower - 
			    N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) + 
			    N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;
      }
    }    
  }
  else
  {
  
    if (N.arcdata[abs(a)].direction == 1 || N.arcdata[abs(a)].direction == 0) {
       Delta =  N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - 
			     N.arcdata[a].resistance_lower * N.arcdata[a].current_lower + 
			     N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

      } else {
        Delta =  - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower)
                             + N.arcdata[a].resistance_lower * N.arcdata[a].current_lower
                             + N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;
      }
  }

  return Delta;
}

template<typename Network>
typename Network::IntegerType compute_alpha(Network& N, arc a)
{
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;
  
  IntegerType Delta = compute_delta(N, a);
  
  IntegerType R_a = N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper + N.arcdata[a].resistance_roof;
  
  RationalType alpha_unrounded = to_double(Delta)/to_double(R_a);
  
  long long alpha_ = round(alpha_unrounded);
  IntegerType alpha = alpha_;
  
  return alpha;
}

template<typename Network, typename IntegerType, typename RationalType>
void update_current(Network& N,
                    SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original,
                    unsigned int edge_index,
                    bool update_roof_edge
                   ){
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;
   
  cout << "Update Current, edge_index = " << edge_index << endl; 
  unsigned int edge_index_tree_edge =  edge_index - RationalType(2) * G0.non_tree_edges.size();
#ifdef VERBOSE
   print_primal_variable(N);
   print_dual_variable(N);
#endif
   arc a = 0;				    
   unsigned int edge_index_ch = 0;
   if(edge_index >= RationalType(2) * G0.non_tree_edges.size())
   {
   a = abs(G0.tree_edges[edge_index_tree_edge]);
   }
  else
  {
    edge_index_ch = edge_index;
    edge_index /= 2;
    cout << "edge_index_ch = " << edge_index_ch << endl;
    cout << "edge_index = " << edge_index << endl;
    a = G0.non_tree_edges[edge_index];
  }
    cout << "edge_index_ch = " << edge_index_ch << endl;
    cout << "edge_index = " << edge_index << endl;
    
  cout << "a = " << a << endl;
  
  cout << "xlower = " << to_double(N.arcdata[a].xlower) << endl;
  cout << "inf = " << to_double(N.arcdata[a].infeasibility) << endl;
  
  if(edge_index >= RationalType(2) * G0.non_tree_edges.size() /*&&
     N.arcdata[a].infeasibility != 0 && N.arcdata[a].sroof != 0 &&
     N.arcdata[a].xlower != 0 && N.arcdata[a].slower != 0 &&
     (N.arcdata[a].capacity - N.arcdata[a].xlower) != 0 && N.arcdata[a].supper != 0*/){
     
    
    unsigned int edge_index_tree_edge =  edge_index - RationalType(2) * G0.non_tree_edges.size();
    cout << "edge_index = " << edge_index << endl;
    cout << "comes here" << endl;
    arc a = G0.tree_edges[edge_index_tree_edge];
    if (a < 0) a = -a;

    assert(N.arcdata[a].xlower != 0 && N.arcdata[a].xlower != N.arcdata[a].capacity && N.arcdata[a].infeasibility != 0);
    
    if ((N.arcdata[a].resistance_roof < N.arcdata[a].resistance_lower || 
         N.arcdata[a].resistance_roof < N.arcdata[a].resistance_upper) ) {
      
      if ( N.arcdata[a].resistance_lower < N.arcdata[a].resistance_upper) {
        
	if (N.arcdata[a].direction == 1) {
	  IntegerType alpha = compute_alpha(N, a);
          N.arcdata[abs(a)].current_roof  -= alpha;      N.arcdata[abs(a)].current_lower += alpha;
          ST_original.update(G0.tail(a), -alpha, 0);     ST_original.update(G0.head(a), alpha, 0);
        } else {	  
	  IntegerType alpha = compute_alpha(N, a);
          N.arcdata[abs(a)].current_roof  += alpha;       N.arcdata[abs(a)].current_lower += alpha;
          ST_original.update(G0.tail(a), alpha, 0);       ST_original.update(G0.head(a), -alpha, 0);
        }
      } else {
        if (N.arcdata[a].direction == 1) {
	  IntegerType alpha = compute_alpha(N, a);
          N.arcdata[a].current_lower -= alpha;          N.arcdata[a].current_roof  += alpha;
          ST_original.update(G0.tail(a), alpha, 0);     ST_original.update(G0.head(a), -alpha, 0);
        } else {
	  IntegerType alpha = compute_alpha(N, a);
          N.arcdata[a].current_lower -= alpha;          N.arcdata[a].current_roof  -= alpha;
          ST_original.update(G0.tail(a), -alpha, 0);    ST_original.update(G0.head(a), alpha, 0);
        }
      }
    } 
    else {
      if (N.arcdata[abs(a)].direction == 1 || N.arcdata[abs(a)].direction == 0) {
 	IntegerType alpha = compute_alpha(N, a);
        N.arcdata[abs(a)].current_roof -= alpha;        N.arcdata[abs(a)].current_lower += alpha;
        ST_original.update(G0.tail(a), alpha, 0);       ST_original.update(G0.head(a), -alpha, 0);
      } else {
	IntegerType alpha = compute_alpha(N, a);
        N.arcdata[abs(a)].current_roof  -= alpha;       N.arcdata[abs(a)].current_lower -= alpha;
        ST_original.update(G0.tail(a), alpha, 0);       ST_original.update(G0.head(a), -alpha, 0);
      }

    }
    return;
  }
  else if (edge_index * RationalType(2) < RationalType(2) * G0.non_tree_edges.size())/*if( N.arcdata[a].infeasibility != 0 && N.arcdata[a].sroof != 0 &&
     N.arcdata[a].xlower != 0 && N.arcdata[a].slower != 0 &&
     (N.arcdata[a].capacity - N.arcdata[a].xlower) != 0 && N.arcdata[a].supper != 0) */{


    cout << "edge_index_ch = " << edge_index_ch << endl;
    cout << "edge_index = " << edge_index << endl;
    const arc a = G0.non_tree_edges[edge_index];
  
    const IntegerType& R_a = G0.resistances_accross_cycles[edge_index];
   
    IntegerType voltage_drop = ST_original.query(edge_index);

    cout << "Modulo = " << edge_index_ch % 2 << endl;
    if (edge_index_ch % 2 == 0) {
      cout << "Here" << endl;
      assert(N.arcdata[a].infeasibility != 0);
      
      if (N.arcdata[a].infeasibility != 0) {
        IntegerType Delta_roof =  N.arcdata[a].resistance_roof * N.arcdata[a].current_roof - voltage_drop;

        IntegerType R_a_roof = R_a - (N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper) + N.arcdata[a].resistance_roof;

#ifdef PathFollowing
        RationalType alpha_roof_ = -to_double(Delta_roof) / to_double(R_a_roof);
        long long alpha_ = round(alpha_roof_);
        IntegerType alpha_roof = alpha_;
	cout << "alpha 8 = " << to_double(alpha_roof) << endl;
#else
        RationalType alpha_unrounded =  -to_double(Delta_roof) / to_double(R_a_roof);
        IntegerType alpha_roof = round_alpha(alpha_unrounded, Delta_roof, R_a_roof);
        if (to_double(alpha_roof) * alpha_unrounded  < 0)
        {
          alpha_roof = -alpha_roof;
        }
#endif

        ST_original.update(-alpha_roof, edge_index);

        N.arcdata[a].current -= alpha_roof;
        N.arcdata[a].current_roof += alpha_roof;
      }
    } else if(N.arcdata[abs(a)].xlower != 0 && N.arcdata[abs(a)].slower != 0
              && N.arcdata[abs(a)].capacity != N.arcdata[abs(a)].xlower && N.arcdata[abs(a)].supper != 0){
      cout << "case edge_index_ch % 2 == 1" << endl;
      IntegerType r_a(0);
      IntegerType r_a_other(0);
      IntegerType f_a(0);
      IntegerType f_a_other(0);

      if((N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) && 
	  N.arcdata[a].xlower != 0 && N.arcdata[a].slower != 0){	
        r_a = N.arcdata[a].resistance_lower;
        f_a = N.arcdata[a].current_lower;
        r_a_other = N.arcdata[a].resistance_upper;
        f_a_other = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
      }else{
        r_a_other = N.arcdata[a].resistance_lower;
        f_a_other = N.arcdata[a].current_lower;
        r_a = N.arcdata[a].resistance_upper;
        f_a = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
      }

      IntegerType D = voltage_drop;
      IntegerType Delta = 0;
      if (N.arcdata[a].direction == 1) {
        if(N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper &&
	  N.arcdata[a].xlower != 0 && N.arcdata[a].slower != 0) {
          Delta = (f_a_other * r_a_other) + D;
        } else {
          Delta =  (f_a_other * r_a_other) - D;
        }
      } else {
        if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
          Delta = (f_a_other * r_a_other) - D;
        } else {
          Delta =  (f_a_other * r_a_other) + D;
        }
      }
      Delta = (f_a * r_a) - Delta;
      RationalType alpha_unrounded = to_double( Delta ) / to_double( R_a );

#ifdef WithoutRounding
      RationalType alpha = alpha_unrounded;
#endif
#ifdef PathFollowing
      long long alpha_ = round(alpha_unrounded);
      IntegerType alpha = alpha_;
      cout << "alpha 9 = " << to_double(alpha) << endl;
#else
      IntegerType alpha = round_alpha(alpha_unrounded, Delta, R_a);

      if (to_double(alpha) * alpha_unrounded  < 0) {
        alpha = -alpha;
      }

#endif

      if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper &&
	  N.arcdata[a].xlower != 0 && N.arcdata[a].slower != 0) {
        if (N.arcdata[a].direction == 1) {
          ST_original.update(alpha, edge_index);
          N.arcdata[a].current += alpha;
          N.arcdata[a].current_lower -= alpha;
        } else {
          ST_original.update(-alpha, edge_index);
          N.arcdata[a].current -= alpha;
          N.arcdata[a].current_lower -= alpha;
        }
      } else {
        if (N.arcdata[a].direction == 1) {
          ST_original.update(-alpha, edge_index);
          N.arcdata[a].current -= alpha;
          N.arcdata[a].current_lower += alpha;
        } else {
          ST_original.update(alpha, edge_index);
          N.arcdata[a].current += alpha;
          N.arcdata[a].current_lower += alpha;
        }
      }
    }
  }
  else
  {
    cout << edge_index * 2 << " >= " << RationalType(2) * G0.non_tree_edges.size() << endl;
    cout << "Do nothing" << endl;
    assert(false);
  }
  
  cout << "update current done" << endl;
}




template<typename Network, typename IntegerType, typename RationalType>
void
compute_tree_induced_voltages(Network& N,
                              SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original)
{

  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;

  node root = ST_original.root;
  N.nodedata[root].voltage = 0;

  for (auto a : G0.tree_edges) {
    assert(ST_original.depth[G0.tail(a)] > ST_original.depth[G0.head(a)]);

    if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower || N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper) {
      IntegerType offset = N.arcdata[abs(a)].current_roof * N.arcdata[abs(a)].resistance_roof;
      if (a > 0) {
        N.nodedata[G0.tail(a)].voltage    = N.nodedata[G0.head(a)].voltage + offset;
      } else {
        N.nodedata[G0.tail(a)].voltage    = N.nodedata[G0.head(a)].voltage - offset;
      }

      if (N.arcdata[abs(a)].resistance_lower < N.arcdata[abs(a)].resistance_upper) {
        IntegerType offset = N.arcdata[abs(a)].resistance_lower * N.arcdata[abs(a)].current_lower;
        if (a > 0) {
          if (N.arcdata[a].direction == 1) {
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tail(a)].voltage - offset;
          } else {
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
          }
        } else {
          if (N.arcdata[-a].direction == 1) {
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
          }
          if (N.arcdata[-a].direction == -1) {
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tail(a)].voltage - offset;
          }
        }
      } else {
        IntegerType offset = (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;
        if (a > 0) {
          if (N.arcdata[a].direction == 1) {
            N.arcdata[a].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
          } else {
            N.arcdata[a].voltage_vw = N.nodedata[G0.tail(a)].voltage - offset;
          }
        } else {
          if (N.arcdata[-a].direction == 1) {
            N.arcdata[abs(a)].voltage_vw =  N.nodedata[G0.tail(a)].voltage - offset;
          } else {
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
          }
        }
      }
    } else {
      if (a > 0) {
        if (N.arcdata[abs(a)].direction == 1) {
          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;


          N.nodedata[G0.tail(a)].voltage = N.arcdata[abs(a)].voltage_vw + N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower;
        } else {

          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - (N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_lower;


          N.nodedata[G0.tail(a)].voltage = N.arcdata[abs(a)].voltage_vw +
                                           (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) *
                                           N.arcdata[abs(a)].resistance_upper;
        }
      } else {
        if (N.arcdata[abs(a)].direction == 1) {

          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tails[abs(a)]].voltage - (N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower);


          N.nodedata[G0.heads[abs(a)]].voltage = N.arcdata[abs(a)].voltage_vw + (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;

        } else {

          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tails[abs(a)]].voltage - (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;

          N.nodedata[G0.heads[abs(a)]].voltage = N.arcdata[abs(a)].voltage_vw + ( N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower);
        }
      }
    }
  }

  for (auto a : G0.non_tree_edges) {
    if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
      if (N.arcdata[a].direction == 1)
      {
        N.arcdata[a].voltage_vw = N.nodedata[G0.head(a)].voltage - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower);
      }
      else
      {
        N.arcdata[a].voltage_vw = N.nodedata[G0.tail(a)].voltage - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower);
      }
    } else {
      if (N.arcdata[a].direction == 1) {
        N.arcdata[a].voltage_vw = N.nodedata[G0.tail(a)].voltage - N.arcdata[a].resistance_lower * N.arcdata[a].current_lower;
      } else {
        N.arcdata[a].voltage_vw = N.nodedata[G0.head(a)].voltage - N.arcdata[a].resistance_lower * N.arcdata[a].current_lower;
      }
    }
  }



#ifndef NDEBUG
  for (unsigned int v = G0.no_of_vertices; v >= 1; v--) {
    if (abs(to_double(N.nodedata[v].voltage - ST_original.query(v, 0))) >= 1e-3) {
      cout << "Node Voltage assert Fails by = " << to_double(N.nodedata[v].voltage - ST_original.query(v, 0)) << endl;
    }
    assert(abs(to_double(N.nodedata[v].voltage - ST_original.query(v, 0))) < 1e-3);
  }
#endif
}





template<typename Network, typename IntegerType, typename RationalType>
void
compute_tree_induced_voltages_and_update_potential(Network& N,
                              SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original)
{

  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;

  node root = ST_original.root;
  N.nodedata[root].voltage = 0;

  for (auto a : G0.tree_edges) {
    assert(ST_original.depth[G0.tail(a)] > ST_original.depth[G0.head(a)]);

    if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower || N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper) {
      IntegerType offset = N.arcdata[abs(a)].current_roof * N.arcdata[abs(a)].resistance_roof;
      if (a > 0) {
        IntegerType old_voltage           = N.nodedata[G0.tail(a)].voltage;
        N.nodedata[G0.tail(a)].voltage    = N.nodedata[G0.head(a)].voltage + offset;
        IntegerType change_in_voltage     = N.nodedata[G0.tail(a)].voltage - old_voltage;
        N.nodedata[G0.tail(a)].potential  -=  change_in_voltage;
      } else {
        IntegerType old_voltage           = N.nodedata[G0.tail(a)].voltage;
        N.nodedata[G0.tail(a)].voltage    = N.nodedata[G0.head(a)].voltage - offset;
        IntegerType change_in_voltage     = N.nodedata[G0.tail(a)].voltage - old_voltage;
        N.nodedata[G0.tail(a)].potential  -=  change_in_voltage;
      }

      if (N.arcdata[abs(a)].resistance_lower < N.arcdata[abs(a)].resistance_upper) {
        IntegerType offset = N.arcdata[abs(a)].resistance_lower * N.arcdata[abs(a)].current_lower;
        if (a > 0) {
          if (N.arcdata[a].direction == 1) {
            IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tail(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
            N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          } else {
            IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
            N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          }
        } else {
          if (N.arcdata[-a].direction == 1) {
            IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
            N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          }
          if (N.arcdata[-a].direction == -1) {
            IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tail(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
            N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          }
        }
      } else {
        IntegerType offset = (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;
        if (a > 0) {
          if (N.arcdata[a].direction == 1) {
            IntegerType old_voltage = N.arcdata[a].voltage_vw;
            N.arcdata[a].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[a].voltage_vw - old_voltage;
            N.arcdata[a].potentialvw -= change_in_voltage;
          } else {
            IntegerType old_voltage = N.arcdata[a].voltage_vw;
            N.arcdata[a].voltage_vw = N.nodedata[G0.tail(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[a].voltage_vw - old_voltage;
            N.arcdata[a].potentialvw -= change_in_voltage;
          }
        } else {
          if (N.arcdata[-a].direction == 1) {
            IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
            N.arcdata[abs(a)].voltage_vw =  N.nodedata[G0.tail(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
            N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          } else {
            IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
            N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - offset;
            IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
            N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          }
        }
      }
    } else {
      if (a > 0) {
        if (N.arcdata[abs(a)].direction == 1) {
          IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;

          IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
          N.arcdata[abs(a)].potentialvw -= change_in_voltage;

          old_voltage = N.nodedata[G0.tail(a)].voltage;
          N.nodedata[G0.tail(a)].voltage = N.arcdata[abs(a)].voltage_vw + N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower;
          change_in_voltage = N.nodedata[G0.tail(a)].voltage - old_voltage;
          N.nodedata[G0.tail(a)].potential -= change_in_voltage;
        } else {
          IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.head(a)].voltage - (N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_lower;

          IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
          N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          old_voltage = N.nodedata[G0.tail(a)].voltage;
          N.nodedata[G0.tail(a)].voltage = N.arcdata[abs(a)].voltage_vw +
                                           (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) *
                                           N.arcdata[abs(a)].resistance_upper;
          change_in_voltage = N.nodedata[G0.tail(a)].voltage - old_voltage;
          N.nodedata[G0.tail(a)].potential -= change_in_voltage;
        }
      } else {
        if (N.arcdata[abs(a)].direction == 1) {
          IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tails[abs(a)]].voltage - (N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower);
          IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
          N.arcdata[abs(a)].potentialvw -= change_in_voltage;

          old_voltage = N.nodedata[G0.heads[abs(a)]].voltage;
          N.nodedata[G0.heads[abs(a)]].voltage = N.arcdata[abs(a)].voltage_vw + (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;
          change_in_voltage = N.nodedata[G0.heads[abs(a)]].voltage - old_voltage;
          N.nodedata[G0.heads[abs(a)]].potential -= change_in_voltage;
        } else {
          IntegerType old_voltage = N.arcdata[abs(a)].voltage_vw;
          N.arcdata[abs(a)].voltage_vw = N.nodedata[G0.tails[abs(a)]].voltage - (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper;
          IntegerType change_in_voltage = N.arcdata[abs(a)].voltage_vw - old_voltage;
          N.arcdata[abs(a)].potentialvw -= change_in_voltage;
          old_voltage = N.nodedata[G0.heads[abs(a)]].voltage;
          N.nodedata[G0.heads[abs(a)]].voltage = N.arcdata[abs(a)].voltage_vw + ( N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower);
          change_in_voltage = N.nodedata[G0.heads[abs(a)]].voltage - old_voltage;
          N.nodedata[G0.heads[abs(a)]].potential -= change_in_voltage;
        }
      }
    }
  }

  for (auto a : G0.non_tree_edges) {
    if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
      if (N.arcdata[a].direction == 1)
      {
        IntegerType old_voltage = N.arcdata[a].voltage_vw;
        N.arcdata[a].voltage_vw = N.nodedata[G0.head(a)].voltage - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower);
        IntegerType change_in_voltage = N.arcdata[a].voltage_vw - old_voltage;
        N.arcdata[a].potentialvw -= change_in_voltage;
      }
      else
      {
        IntegerType old_voltage = N.arcdata[a].voltage_vw;
        N.arcdata[a].voltage_vw = N.nodedata[G0.tail(a)].voltage - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower);
        IntegerType change_in_voltage = N.arcdata[a].voltage_vw - old_voltage;
        N.arcdata[a].potentialvw -= change_in_voltage;
      }
    } else {
      if (N.arcdata[a].direction == 1) {
        IntegerType old_voltage = N.arcdata[a].voltage_vw;
        N.arcdata[a].voltage_vw = N.nodedata[G0.tail(a)].voltage - N.arcdata[a].resistance_lower * N.arcdata[a].current_lower;
        IntegerType change_in_voltage = N.arcdata[a].voltage_vw - old_voltage;
        N.arcdata[a].potentialvw -= change_in_voltage;
      } else {
        IntegerType old_voltage = N.arcdata[a].voltage_vw;
        N.arcdata[a].voltage_vw = N.nodedata[G0.head(a)].voltage - N.arcdata[a].resistance_lower * N.arcdata[a].current_lower;
        IntegerType change_in_voltage = N.arcdata[a].voltage_vw - old_voltage;
        N.arcdata[a].potentialvw -= change_in_voltage;
      }
    }
  }



#ifndef NDEBUG
  for (unsigned int v = G0.no_of_vertices; v >= 1; v--) {
    if (abs(to_double(N.nodedata[v].voltage - ST_original.query(v, 0))) >= 1e-3) {
      cout << "Node Voltage assert Fails by = " << to_double(N.nodedata[v].voltage - ST_original.query(v, 0)) << endl;
    }
    assert(abs(to_double(N.nodedata[v].voltage - ST_original.query(v, 0))) < 1e-3);
  }
#endif
}


template<typename Network, typename IntegerType, typename RationalType>
void
compute_tree_induced_voltages_layered_new(Network & N,
    SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original)
{
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;

  node root = ST_original.root;
  vector<IntegerType> TIV(G0.no_of_vertices + 1, 0);
  TIV.reserve(G0.no_of_vertices + 1);
  TIV[root] = 0;
  N.nodedata[root].voltage = 0;

  for ( auto a : G0.tree_edges ) {
    assert( a != 0 );
    assert( ST_original.depth[G0.tail(a)] > ST_original.depth[G0.head(a)] );

#ifdef VERBOSE
    cout << a << " = (" << G0.tail(a) << "," << G0.head(a) << ")" << endl;
#endif
    if ( a > 0 ) {

      const IntegerType& r_a_1 = N.arcdata[a].resistance_upper;
      IntegerType f_a_1 =  N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
      if (f_a_1 != 0 && r_a_1 != 0) {
        multiply_by_bit_shift(f_a_1, r_a_1);
      }
      else {
        f_a_1 = 0;
      }
      N.arcdata[a].voltage_vw = TIV[G0.heads[a]] - f_a_1;
      const IntegerType& r_a_2 = N.arcdata[a].resistance_lower;
      IntegerType f_a_2 = N.arcdata[a].current_lower;
#ifdef VERBOSE
      cout << "f_a_2 = " << to_double(f_a_2) << endl;
      cout << "r_a_2 = " << to_double(r_a_2) << endl;
#endif
      if (f_a_2 != 0 && r_a_2 != 0) {
        multiply_by_bit_shift(f_a_2, r_a_2);
      }
      else {
        f_a_2 = 0;
      }
      TIV[G0.tails[a]] = N.arcdata[a].voltage_vw + f_a_2;
#ifdef VERBOSE
      cout << "TIV[" << G0.tails[a] << "] = " << to_double(N.arcdata[a].voltage_vw) << " + " << to_double(f_a_2) << endl;
#endif

      N.nodedata[G0.heads[a]].voltage = TIV[G0.heads[a]];
#ifdef VERBOSE
      cout << "N.nodedata[" << G0.heads[a] << "].voltage = " << to_double(TIV[G0.heads[a]]) << "]" << endl;
#endif

      N.nodedata[G0.tails[a]].voltage = TIV[G0.tails[a]];
#ifdef VERBOSE
      cout << "N.nodedata[" << G0.tails[a] << "].voltage = " << to_double(TIV[G0.tails[a]]) << endl;
#endif
    } else {

      unsigned int rev = -a;
#ifdef VERBOSE
      cout << "rev = " << -a << endl;
#endif
      const IntegerType& r_a_2 = N.arcdata[rev].resistance_lower;
      IntegerType f_a_2 = N.arcdata[rev].current_lower;
      if (f_a_2 != 0 && r_a_2 != 0) {
        multiply_by_bit_shift(f_a_2, r_a_2);
      }
      else {
        f_a_2 = 0;
      }
      N.arcdata[rev].voltage_vw = TIV[G0.tails[rev]] - f_a_2;
#ifdef VERBOSE
      cout << "N.arcdata[" << rev << "].voltage_vw = " <<  to_double(TIV[G0.tails[rev]]) << " - " <<  to_double(f_a_2) << endl;
#endif
      const IntegerType& r_a_1 = N.arcdata[rev].resistance_upper;
      IntegerType f_a_1 = N.arcdata[rev].cur_src_vw - N.arcdata[rev].current_lower;
      if (f_a_1 != 0 && r_a_1 != 0) {
        multiply_by_bit_shift(f_a_1, r_a_1);
      }
      else {
        f_a_1 = 0;
      }
      TIV[G0.heads[rev]] = N.arcdata[rev].voltage_vw + f_a_2;
#ifdef VERBOSE
      cout << "TIV[" << G0.heads[rev] << "] = " << to_double(N.arcdata[rev].voltage_vw) << " + " << to_double(f_a_2) << endl;
#endif

      N.nodedata[G0.heads[rev]].voltage = TIV[G0.heads[rev]];
      N.nodedata[G0.tails[rev]].voltage = TIV[G0.tails[rev]];
    }
#ifdef VERBOSE
    cout << to_double(G.tree_induced_voltages[G.tail(a)]) << " = " << to_double(G.tree_induced_voltages[G.head(a)])
         << " + " << to_double(G.resistances[abs(a)]) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)]) << endl;
#endif
  }


  for (auto a : G0.non_tree_edges)
  {
    if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper)
    {
      const IntegerType& r_a_2 =  N.arcdata[a].resistance_upper;
      IntegerType f_a_2 = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
      if (f_a_2 != 0 && r_a_2 != 0) {
        multiply_by_bit_shift(f_a_2, r_a_2);
      }
      else {
        f_a_2 = 0;
      }

      N.arcdata[a].voltage_vw = TIV[G0.tails[a]] - f_a_2;

    }
    else {
      const IntegerType& r_a_1 = N.arcdata[a].resistance_lower;
      IntegerType f_a_1 = N.arcdata[a].current_lower;
      if (f_a_1 != 0 && r_a_1 != 0) {
        multiply_by_bit_shift(f_a_1, r_a_1);
      }
      else {
        f_a_1 = 0;
      }
      N.arcdata[a].voltage_vw = TIV[G0.heads[a]] - f_a_1;
    }
    N.nodedata[G0.heads[a]].voltage = TIV[G0.heads[a]];
    N.nodedata[G0.tails[a]].voltage = TIV[G0.tails[a]];
  }

#ifdef VERBOSE
  for (unsigned int i = 1; i <= G.no_of_verticies; i++)
  {
    cout << to_double(TIV[i]) << endl;
  }
#endif
#if !defined(NDEBUG) || defined(VERBOSE)
#ifdef VERBOSE
  for ( auto a : G.tree_edges ) {
#ifdef VERBOSE
    cout << a << " = (" << G.tail(a) << "," << G.head(a) << ")" << " " << to_double( G.tree_induced_voltages[G.tail(a)] ) << " - " <<  to_double( G.tree_induced_voltages[G.head(a)] ) << " - " << to_double( G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) << endl;
#endif

    if (fabs( to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] )  ) >= 1e-3)
    {
      cout << to_double( G.tree_induced_voltages[G.tail(a)] ) << " - " << to_double(G.tree_induced_voltages[G.head(a)]  )
           << " - " << to_double(G.resistances[abs(a)] ) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)] ) << endl;
      cout << to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) << endl;

      cout << "Assertion: " << fabs( to_double(  (G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)]) - ( G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) ) )  << " < - > " << 0.0 << endl;
      cout << to_double(G.tree_induced_voltages[G.tail(a)]) << endl;
      cout << to_double(G.tree_induced_voltages[G.head(a)]) << endl;
      cout << to_double(G.tree_induced_voltages[G.tail(a)] - G.tree_induced_voltages[G.head(a)]) << endl;
      cout << to_double(G.resistances[abs(a)]) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)]) << endl;
      cout << to_double(G.resistances[abs(a)] * sign(a)) << " * " << to_double(G.currents[abs(a)]) << endl;
      cout << to_double(G.resistances[abs(a)] * sign(a) * G.currents[abs(a)]) << endl;
    }

  }
#endif
#endif
}
/** \brief Computes the Tree Induced Voltages
 *
 * @param ST The Spanning Tree
 *
 */
template<typename IntegerType, typename RationalType>
void
compute_tree_induced_voltages_layered_DS( SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> &ST) {
  Graph<IntegerType, RationalType>& G = ST.G;

#ifdef VERBOSE
  cout << "computing tree induced voltages..." << endl;
#endif
  node root = G.root;
  vector<IntegerType> TIV(G.no_of_verticies + 1, 0);
  TIV.reserve(G.no_of_verticies + 1);
  TIV[root] = 0;
  for (auto a : G.incident_edges[root]) {
    if (G.original_auxiliary_transformed_corrospondence[abs(a)].is_invalid == true) continue;
    if (G.head(a) == root) {
      G.original_auxiliary_transformed_corrospondence[abs(a)].tree_induced_voltage_w = 0;
    }
    else if (G.tail(a) == root) {
      G.original_auxiliary_transformed_corrospondence[abs(a)].tree_induced_voltage_v = 0;
    }
    else {

      assert(false);
    }
  }

  assert( G.head( G.tree_edges.front() ) == G.root );

  for ( auto a : G.tree_edges ) {
    assert( a != 0 );
    assert( ST.depth[G.tail(a)] > ST.depth[G.head(a)] );

#ifdef VERBOSE
    cout << a << " = (" << G.tail(a) << "," << G.head(a) << ")" << endl;
#endif
    if ( a > 0 ) {

      const IntegerType& r_a_1 = G.original_auxiliary_transformed_corrospondence[a].r1_tilde;
      IntegerType f_a_1 = G.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde;
      if (f_a_1 != 0 && r_a_1 != 0) {
        multiply_by_bit_shift(f_a_1, r_a_1);
      }
      else {
        f_a_1 = 0;
      }
      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_u = TIV[G.heads[a]] - f_a_1;
      const IntegerType& r_a_2 = G.original_auxiliary_transformed_corrospondence[a].r2_tilde;
      IntegerType f_a_2 = G.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde;
      if (f_a_2 != 0 && r_a_2 != 0) {
        multiply_by_bit_shift(f_a_2, r_a_2);
      }
      else {
        f_a_2 = 0;
      }
      TIV[G.tails[a]] = G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_u + f_a_2;
      const IntegerType& r_a_3 = G.original_auxiliary_transformed_corrospondence[a].r3_tilde;
      IntegerType f_a_3 = G.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde;
      if (f_a_3 != 0 && r_a_3 != 0) {
        multiply_by_bit_shift(f_a_3, r_a_3);
      }
      else {
        f_a_3 = 0;
      }
      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_vw = G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_u + f_a_3;
      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_w = TIV[G.heads[a]];
      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_v = TIV[G.tails[a]];

    } else {

      unsigned int rev = -a;
      const IntegerType& r_a_2 = G.original_auxiliary_transformed_corrospondence[rev].r2_tilde;
      IntegerType f_a_2 =  G.original_auxiliary_transformed_corrospondence[rev].electrical_flow2_tilde;
      if (f_a_2 != 0 && r_a_2 != 0) {
        multiply_by_bit_shift(f_a_2, r_a_2);
      }
      else {
        f_a_2 = 0;
      }
      G.original_auxiliary_transformed_corrospondence[rev].tree_induced_voltage_u = TIV[G.tails[rev]] - f_a_2;
      const IntegerType& r_a_1 = G.original_auxiliary_transformed_corrospondence[rev].r1_tilde;
      IntegerType f_a_1 = G.original_auxiliary_transformed_corrospondence[rev].electrical_flow1_tilde;
      if (f_a_1 != 0 && r_a_1 != 0) {
        multiply_by_bit_shift(f_a_1, r_a_1);
      }
      else {
        f_a_1 = 0;
      }
      TIV[G.heads[rev]] = G.original_auxiliary_transformed_corrospondence[rev].tree_induced_voltage_u + f_a_1;
      const IntegerType& r_a_3 = G.original_auxiliary_transformed_corrospondence[rev].r3_tilde;
      IntegerType f_a_3 = G.original_auxiliary_transformed_corrospondence[rev].electrical_flow3_tilde;
      if (f_a_3 != 0 && r_a_3 != 0) {
        multiply_by_bit_shift(f_a_3, r_a_3);
      }
      else {
        f_a_3 = 0;
      }
      G.original_auxiliary_transformed_corrospondence[rev].tree_induced_voltage_vw = G.original_auxiliary_transformed_corrospondence[rev].tree_induced_voltage_u + f_a_3;

      G.original_auxiliary_transformed_corrospondence[rev].tree_induced_voltage_w = TIV[G.heads[rev]];
      G.original_auxiliary_transformed_corrospondence[rev].tree_induced_voltage_v = TIV[G.tails[rev]];
    }
#ifdef VERBOSE
    cout << to_double(G.tree_induced_voltages[G.tail(a)]) << " = " << to_double(G.tree_induced_voltages[G.head(a)])
         << " + " << to_double(G.resistances[abs(a)]) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)]) << endl;
#endif
  }

  for (auto a : G.non_tree_edges)
  {
    if (G.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;

    if (G.original_auxiliary_transformed_corrospondence[a].non_tree_edge == 1)
    {
      const IntegerType& r_a_2 =  G.original_auxiliary_transformed_corrospondence[a].r2_tilde;
      IntegerType f_a_2 = G.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde;
      if (f_a_2 != 0 && r_a_2 != 0) {
        multiply_by_bit_shift(f_a_2, r_a_2);
      }
      else {
        f_a_2 = 0;
      }
      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_u = TIV[G.tails[a]] - f_a_2;

      const IntegerType& r_a_3 = G.original_auxiliary_transformed_corrospondence[a].r3_tilde;
      IntegerType f_a_3 = G.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde;
      if (f_a_3 != 0 && r_a_3 != 0) {
        multiply_by_bit_shift(f_a_3, r_a_3);
      }
      else {
        f_a_3 = 0;
      }
      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_vw = G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_u + f_a_3;
    }
    else {
      const IntegerType& r_a_1 = G.original_auxiliary_transformed_corrospondence[a].r1_tilde;
      IntegerType f_a_1 = G.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde;
      if (f_a_1 != 0 && r_a_1 != 0) {
        multiply_by_bit_shift(f_a_1, r_a_1);
      }
      else {
        f_a_1 = 0;
      }
      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_u = TIV[G.heads[a]] - f_a_1;
      const IntegerType& r_a_3 = G.original_auxiliary_transformed_corrospondence[a].r3_tilde;
      IntegerType f_a_3 = G.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde;

      if (f_a_3 != 0 && r_a_3 != 0) {
        multiply_by_bit_shift(f_a_3, r_a_3);
      }
      else {

        f_a_3 = 0;
      }
      if (f_a_3 != G.original_auxiliary_transformed_corrospondence[a].r3_tilde * G.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde)
      {
        IntegerType product = G.original_auxiliary_transformed_corrospondence[a].r3_tilde * G.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde;
        cout << to_double(G.original_auxiliary_transformed_corrospondence[a].r3_tilde) << " * " << to_double(G.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde) << endl;
        cout << to_double(f_a_3) << " < - > " << to_double(product) << endl;
        cout << f_a_3.digits.size() << " < - > " << product.digits.size() << endl;
      }
      assert(f_a_3 == G.original_auxiliary_transformed_corrospondence[a].r3_tilde * G.original_auxiliary_transformed_corrospondence[a].electrical_flow3_tilde);

      G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_vw = G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_u + f_a_3;

    }

    G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_w = TIV[G.heads[a]];
    G.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_v = TIV[G.tails[a]];
  }

#ifdef VERBOSE
  for (unsigned int i = 1; i <= G.no_of_verticies; i++)
  {
    cout << to_double(TIV[i]) << endl;
  }
#endif
#if !defined(NDEBUG) || defined(VERBOSE)
#ifdef VERBOSE
  for ( auto a : G.tree_edges ) {
#ifdef VERBOSE
    cout << a << " = (" << G.tail(a) << "," << G.head(a) << ")" << " " << to_double( G.tree_induced_voltages[G.tail(a)] ) << " - " <<  to_double( G.tree_induced_voltages[G.head(a)] ) << " - " << to_double( G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) << endl;
#endif

    if (fabs( to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] )  ) >= 1e-3)
    {
      cout << to_double( G.tree_induced_voltages[G.tail(a)] ) << " - " << to_double(G.tree_induced_voltages[G.head(a)]  )
           << " - " << to_double(G.resistances[abs(a)] ) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)] ) << endl;
      cout << to_double( G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)] ) - to_double(G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) << endl;

      cout << "Assertion: " << fabs( to_double(  (G.tree_induced_voltages[G.tail(a)] -  G.tree_induced_voltages[G.head(a)]) - ( G.resistances[abs(a)]*sign(a)*G.currents[abs(a)] ) ) )  << " < - > " << 0.0 << endl;
      cout << to_double(G.tree_induced_voltages[G.tail(a)]) << endl;
      cout << to_double(G.tree_induced_voltages[G.head(a)]) << endl;
      cout << to_double(G.tree_induced_voltages[G.tail(a)] - G.tree_induced_voltages[G.head(a)]) << endl;
      cout << to_double(G.resistances[abs(a)]) << " * " << sign(a) << " * " << to_double(G.currents[abs(a)]) << endl;
      cout << to_double(G.resistances[abs(a)] * sign(a)) << " * " << to_double(G.currents[abs(a)]) << endl;
      cout << to_double(G.resistances[abs(a)] * sign(a) * G.currents[abs(a)]) << endl;
    }

  }
#endif
#endif
}

template<typename Network, typename IntegerType, typename RationalType>
RationalType compute_gap_by_non_tree_currents(Network & N,
    SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original,
    RationalType THRESHOLD_X,
    RationalType THRESHOLD_S
                                             )
{
#ifdef VERBOSE
  cout << "compute gap" << endl;
#endif
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;
  RationalType gap = 0;
  unsigned int edge_index = 0;

  for (auto a : G0.non_tree_edges)
  {
#ifdef VERBOSE
    cout << "a = " << a << endl;
    cout << "edge_index = " << edge_index << endl;
#endif

    const IntegerType& D = ST_original.query(edge_index);

    edge_index++;

    IntegerType r_a(0);
    IntegerType r_a_other(0);
    IntegerType f_a(0);
    IntegerType f_a_other(0);

    if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper)
    {
#ifdef VERBOSE
      cout << "Lower Edge is the Non Tree Edge" << endl;
#endif
      r_a = N.arcdata[a].resistance_lower;
      f_a = N.arcdata[a].current_lower;
      r_a_other = N.arcdata[a].resistance_upper;
      f_a_other = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
    }
    else
    {
#ifdef VERBOSE
      cout << "Upper Edge is the Non Tree Edge" << endl;
#endif
      r_a_other = N.arcdata[a].resistance_lower;
      f_a_other = N.arcdata[a].current_lower;
      r_a = N.arcdata[a].resistance_upper;
      f_a = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
    }


    if (N.arcdata[a].infeasibility != 0)
    {

#ifdef VERBOSE
      cout << "Computation of gap, infeasibility not zero.. direction is " << N.arcdata[a].direction << endl;
#endif
      IntegerType Delta_roof = N.arcdata[a].current_roof * N.arcdata[a].resistance_roof - D;
      RationalType contribution = to_double(Delta_roof * Delta_roof) / to_double(N.arcdata[a].resistance_roof);

#ifdef VERBOSE
      cout << "contribution = " << contribution << endl;
#endif

#ifdef VERBOSE
      cout << "gap = " << gap << endl;
      cout << gap << " + = " << contribution << endl;
      cout << " * " << a << endl;
#endif
      gap += contribution;
    }

    IntegerType Delta = 0;
    if (N.arcdata[a].direction == 1)
    {
      if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
        Delta = (f_a_other * r_a_other) + D;
      }
      else {
        Delta =  (f_a_other * r_a_other) - D;
      }
    }
    else
    {
      if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
        Delta = (f_a_other * r_a_other) - D;
      }
      else {
        Delta =  (f_a_other * r_a_other) + D;
      }
    }


    Delta = (f_a * r_a) - Delta;
#ifdef VERBOSE
    cout << endl << "Delta = " << to_double(Delta) << endl << "r_a = " << to_double(r_a) << endl;
#endif

    const IntegerType DeltaSquare = Delta * Delta;

#ifdef VERBOSE
    cout << "DeltaSquare = " << to_double(DeltaSquare) << " / " << "r_a = " << to_double(r_a) << endl;
#endif
    const RationalType contribution = to_double( DeltaSquare ) / to_double( r_a );

#ifdef VERBOSE
    const RationalType& g_a = G.batteries[i];
    const IntegerType KVL_violation = D - f_a * r_a;
    cout << i << fixed << ":\t\t" << to_double(f_a) << "\t\t" << to_double(r_a) << "\t\t" << g_a << "\t\t" << to_double(D) << "\t\t" << contribution << "\t\t" << to_double(KVL_violation) << endl;
#endif

#ifdef VERBOSE
    cout << gap << " += " << contribution << endl;
    cout << " * " << a << endl;
#endif

    gap += contribution;
#ifdef VERBOSE
    cout << "gap so far = " << gap << endl;
#endif
  }

  for (auto a : G0.tree_edges)
  {
    if (a < 0) a = -a;
    //if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper)
    if (N.arcdata[a].resistance_roof < N.arcdata[a].resistance_lower ||
        N.arcdata[a].resistance_roof < N.arcdata[a].resistance_upper)
    {
      if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper)
      {
#ifdef VERBOSE
        cout << "Lower Edge is the non tree edge" << endl;
#endif
        if (N.arcdata[a].xlower < THRESHOLD_X || N.arcdata[a].slower < THRESHOLD_S)
        {
          continue;
        }

        if (N.arcdata[a].direction == 1)
        {

          assert(N.arcdata[a].direction != 0);

          IntegerType Delta = N.arcdata[a].resistance_lower * N.arcdata[a].current_lower -
                              N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) -
                              N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          IntegerType DeltaSquare = Delta * Delta;
          IntegerType r_a = N.arcdata[a].resistance_lower;
          RationalType contribution = to_double(DeltaSquare) / to_double(r_a);
          gap += contribution;
        }
        else
        {
          IntegerType Delta = N.arcdata[a].resistance_lower * N.arcdata[a].current_lower -
                              N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) +
                              N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          IntegerType DeltaSquare = Delta * Delta;
          IntegerType r_a = N.arcdata[a].resistance_lower;
          RationalType contribution = to_double(DeltaSquare) / to_double(r_a);
          gap += contribution;
        }

      }
      else
      {
#ifdef VERBOSE
        cout << "Upper Edge is the non tree edge" << endl;
#endif
        RationalType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
        if (xupper < THRESHOLD_X || N.arcdata[a].supper < THRESHOLD_S)
        {
          continue;
        }
        if (N.arcdata[a].direction == 1)
        {
          IntegerType Delta = N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) -
                              N.arcdata[a].resistance_lower * N.arcdata[a].current_lower
                              + N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          IntegerType DeltaSquare = Delta * Delta;
          IntegerType r_a = N.arcdata[a].resistance_upper;
          RationalType contribution = to_double(DeltaSquare) / to_double(r_a);
          gap += contribution;
        }

        if (N.arcdata[a].direction == -1)
        {
          IntegerType Delta = N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) -
                              N.arcdata[a].resistance_lower * N.arcdata[a].current_lower
                              - N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          IntegerType DeltaSquare = Delta * Delta;
          IntegerType r_a = N.arcdata[a].resistance_upper;
          RationalType contribution = to_double(DeltaSquare) / to_double(r_a);
          gap += contribution;
        }
      }
    }
    else
    {
#ifdef VERBOSE
      cout << "Roof Edge is the tree edge" << endl;
#endif
      //assert(N.arcdata[a].direction != 0);
      if (N.arcdata[a].infeasibility < THRESHOLD_X || N.arcdata[a].sroof < THRESHOLD_S)
      {
        continue;
      }

      if (N.arcdata[a].direction == 1)
      {
        IntegerType Delta = N.arcdata[a].resistance_roof * N.arcdata[a].current_roof -
                            N.arcdata[a].resistance_lower * N.arcdata[a].current_lower +
                            N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower);


        IntegerType DeltaSquare = Delta * Delta;
        IntegerType r_a = N.arcdata[a].resistance_roof;
        RationalType contribution = to_double(DeltaSquare) / to_double(r_a);
        gap += contribution;
      }

      if (N.arcdata[a].direction == -1)
      {
        IntegerType Delta = N.arcdata[a].resistance_roof * N.arcdata[a].current_roof +
                            N.arcdata[a].resistance_lower * N.arcdata[a].current_lower -
                            N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower);

        IntegerType DeltaSquare = Delta * Delta;
        IntegerType r_a = N.arcdata[a].resistance_roof;
        RationalType contribution = to_double(DeltaSquare) / to_double(r_a);
        gap += contribution;
      }
    }


  }

#ifdef VERBOSE
  cout << "gap = " << gap << endl;
#endif

  return gap;
}

template<typename Network, typename IntegerType, typename RationalType>
RationalType compute_gap_by_non_tree_currents_new(Network & N,
    SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType> & ST_original)
{

//#ifdef VERBOSE
  cout << endl << "compute gap by non tree currents" << endl;
//#endif

  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;
  RationalType gap = 0;
  unsigned int edge_index = 0;
#ifndef NDEBUG
  unsigned int edge_index_with_max_query = 0;
#endif
  IntegerType max_query(0);
  for ( auto a : G0.non_tree_edges ) {
    //if (G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;

    cout << "a = " << a << endl;
    cout << "edge_index = " << edge_index << endl;

    const IntegerType& D = ST_original.query(edge_index);

    if (D > max_query) {
      max_query = D;
#ifndef NDEBUG
      edge_index_with_max_query = edge_index;
#endif
    }

#ifdef VERBOSE
    cout << "D = " << to_double(D) << " " << "ST  = " << to_double(ST_original.query(edge_index)) << endl;
#endif
    edge_index++;

    IntegerType r_a(0);
    IntegerType r_a_other(0);
    IntegerType f_a(0);
    IntegerType f_a_other(0);
    if (N.arcdata[a].resistance_lower >   N.arcdata[a].resistance_upper) {
      r_a = N.arcdata[a].resistance_lower;
      f_a = N.arcdata[a].current_lower;
      r_a_other = N.arcdata[a].resistance_upper;
      f_a_other = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
    }
    else {
      r_a_other = N.arcdata[a].resistance_lower;
      f_a_other = N.arcdata[a].current_lower;
      r_a = N.arcdata[a].resistance_upper;
      f_a = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
    }

#ifdef VERBOSE
    cout << "Delta-offset = " << to_double(r_a_other * f_a_other) << " or " << to_double(f_a * r_a) << endl;
#endif
#ifdef VERBOSE
    cout << "f_a = " << to_double(f_a) << endl;
#endif

    IntegerType Delta = 0;
    if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
      if (f_a_other != 0 && r_a_other != 0) {
        multiply_by_bit_shift(f_a_other, r_a_other);
      }
      else {

        f_a_other = 0;
      }
      Delta = f_a_other - D;
    }
    else {
      if (f_a_other != 0 && r_a_other != 0) {
        multiply_by_bit_shift(f_a_other, r_a_other);
      }
      else {

        f_a_other = 0;
      }
      Delta =  f_a_other + D;
    }
    if (f_a != 0 && r_a != 0) {

      multiply_by_bit_shift(f_a, r_a);
    }
    else {
      f_a = 0;
    }
    Delta = f_a - Delta;

    IntegerType f_a_dummy = f_a;

    //#ifdef VERBOSE
    cout << endl << "Delta = " << to_double(Delta) << endl << "r_a = " << to_double(r_a) << endl;
    //#endif

    const IntegerType DeltaSquare = Delta * Delta;

//#ifdef VERBOSE
    cout << "DeltaSquare = " << to_double(DeltaSquare) << " / " << "r_a = " << to_double(r_a) << endl;
//#endif
    const RationalType contribution = to_double( DeltaSquare ) / to_double( r_a );

#ifdef VERBOSE
    const RationalType& g_a = G.batteries[i];
    const IntegerType KVL_violation = D - f_a * r_a;
    cout << i << fixed << ":\t\t" << to_double(f_a) << "\t\t" << to_double(r_a) << "\t\t" << g_a << "\t\t" << to_double(D) << "\t\t" << contribution << "\t\t" << to_double(KVL_violation) << endl;
#endif
    gap += contribution;
//#ifdef VERBOSE
    cout << "gap so far = " << gap << endl;
//#endif


  }

  for (auto a : G0.tree_edges)
  {

    if (a < 0) a = -a;
    //if(N.arcdata[abs(a)].infeasibility == 0) continue;
    cout <<  "a = " << a << endl;
    IntegerType Delta;

    if (N.arcdata[a].direction == 1)
    {
      Delta =  N.arcdata[a].resistance_roof * N.arcdata[a].current_roof +
               N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) -
               N.arcdata[a].resistance_lower * N.arcdata[a].current_lower;
    }
    else {
      Delta =  N.arcdata[a].resistance_roof * N.arcdata[a].current_roof -
               N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) +
               N.arcdata[a].resistance_lower * N.arcdata[a].current_lower;
    }

    IntegerType DeltaSquare = Delta * Delta;
#ifdef VERBOSE
    cout << "DeltaSquare = " << to_double(DeltaSquare) << endl;
#endif
    RationalType contribution = to_double(DeltaSquare) / to_double(N.arcdata[a].resistance_roof);
#ifdef VERBOSE
    cout << "contribution = " << contribution << endl;
#endif
    gap += contribution;
  }
#ifdef VERBOSE
  cout << "edge_index_with_max_query = " << edge_index_with_max_query << endl;
#endif

  assert(edge_index_with_max_query >= 0);

#ifdef VERBOSE
  cout << "gapcheck = " << gap << endl;
#endif
  return gap;
}

/** \brief Computes the Electrical Flow gap
 *
 * @param G0 The Original Graph
 * @param ST Spanning Tree
 *
 */
template<typename IntegerType, typename RationalType>
RationalType
compute_gap_by_non_tree_currents_DS( Graph<IntegerType, RationalType>& G0,
                                     SpanningTree<Graph<IntegerType, RationalType>,
                                     IntegerType, RationalType> &ST
                                   )
{

#ifdef VERBOSE
  cout << "compute gap by non-tree currents" << endl;
#endif

#ifdef VERBOSE
  cout << "currents in non_tree_edges: " <<  endl;
  for (auto a : G.non_tree_edges) {
    cout << to_double(G.currents[a]) << " , " ;
  }
  cout << endl;
#endif


#ifdef VERBOSE
  cout << "arc\t\tflow\t\t\tresistance\t\tbatteries\t\tb-flow\t\tvoltage\t\tgap-contrib.\t\tKVL viol." << endl;
#endif

  RationalType gap = 0;
  unsigned int edge_index = 0;
#ifndef NDEBUG
  unsigned int edge_index_with_max_query = 0;
#endif
  IntegerType max_query(0);
  for ( auto i : G0.non_tree_edges ) {
    if (G0.original_auxiliary_transformed_corrospondence[i].is_invalid) continue;
#ifdef VERBOSE
    cout << i << " :  " << G0.original_auxiliary_transformed_corrospondence[i].non_tree_edge << endl;
#endif

    const IntegerType& D = ST.query(edge_index);

    if (D > max_query) {
      max_query = D;
#ifndef NDEBUG
      edge_index_with_max_query = edge_index;
#endif
    }

#ifdef VERBOSE
    cout << "D = " << to_double(D) << " " << "ST  = " << to_double(ST.query(edge_index)) << endl;
#endif
    edge_index++;

#ifdef VERBOSE
    cout << "non tree edge = " << G0.original_auxiliary_transformed_corrospondence[i].non_tree_edge << endl;
#endif
    const IntegerType& r_a = G0.original_auxiliary_transformed_corrospondence[i].resistance_non_tree_edge;

    IntegerType f_a = G0.original_auxiliary_transformed_corrospondence[i].current_non_tree_edge;

    const IntegerType& r_a_other = G0.original_auxiliary_transformed_corrospondence[i].resistance_other_edge;
    IntegerType f_a_other = G0.original_auxiliary_transformed_corrospondence[i].current_other_edge;
#ifdef VERBOSE
    cout << "Delta-offset = " << to_double(r_a_other * f_a_other) << " or " << to_double(f_a * r_a) << endl;
#endif
#ifdef VERBOSE
    cout << "f_a = " << to_double(f_a) << endl;
#endif

    IntegerType Delta = 0;
    if (G0.original_auxiliary_transformed_corrospondence[i].non_tree_edge == 1) {
      if (f_a_other != 0 && r_a_other != 0) {
        multiply_by_bit_shift(f_a_other, r_a_other);
      }
      else {

        f_a_other = 0;
      }
      Delta = f_a_other - D;
    }
    else {
      if (f_a_other != 0 && r_a_other != 0) {
        multiply_by_bit_shift(f_a_other, r_a_other);
      }
      else {

        f_a_other = 0;
      }
      Delta =  f_a_other + D;
    }
    if (f_a != 0 && r_a != 0) {

      multiply_by_bit_shift(f_a, r_a);
    }
    else {
      f_a = 0;
    }
    Delta = f_a - Delta;

    IntegerType f_a_dummy = f_a;

#ifdef VERBOSE
    cout << endl << "Delta = " << to_double(Delta) << endl << "r_a = " << to_double(r_a) << endl;
#endif

    const IntegerType DeltaSquare = Delta * Delta;
    const RationalType contribution = to_double( DeltaSquare ) / to_double( r_a );

#ifdef VERBOSE
    const RationalType& g_a = G.batteries[i];
    const IntegerType KVL_violation = D - f_a * r_a;
    cout << i << fixed << ":\t\t" << to_double(f_a) << "\t\t" << to_double(r_a) << "\t\t" << g_a << "\t\t" << to_double(D) << "\t\t" << contribution << "\t\t" << to_double(KVL_violation) << endl;
#endif
    gap += contribution;

  }

#ifdef VERBOSE
  cout << "edge_index_with_max_query = " << edge_index_with_max_query << endl;
#endif

  assert(edge_index_with_max_query >= 0);

#ifdef VERBOSE
  cout << "gapcheck = " << gap << endl;
#endif
  return gap;
}



template<typename Network, typename IntegerType, typename RationalType>
void build_initial_state(Network & N,
                         SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original)
{
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;

  ST_original.clear_d();
  vector<IntegerType> carry(G0.no_of_vertices + 1, 0);

  for (auto a : boost::adaptors::reverse( G0.tree_edges)){
    node v = G0.tail(a);
    node w = G0.head(a);

    if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower || N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper) {
      IntegerType outflow = N.arcdata[abs(a)].current_roof;
      if (a < 0){
        outflow = -outflow;
      }
      IntegerType alpha = outflow - carry[v];

#ifdef VERBOSE
      cout << "a = " << a << endl;
      cout << "outflow = " << to_double(outflow) << endl;
      cout << "carry[" << v << "] = " << to_double(carry[v]) << endl;
      cout << "alpha = " << to_double(alpha) << endl;
      cout << "v = " << v << endl;
#endif

      IntegerType height_offset = 0;
      ST_original.update_with_DS(v, alpha, 0, height_offset);
      carry[w] += outflow;


#ifdef VERBOSE
      print_queries(N, ST_original);
#endif

    }
    else
    {
#ifdef VERBOSE
      cout << "Lower Edge and Upper Edge are tree edges" << endl;
      cout << "updating the current of the lower edge" << endl;
#endif
      IntegerType height_offset = 0;
      IntegerType outflow = 0;
      if (a > 0)
      {
        if (N.arcdata[abs(a)].direction == 1)
        {
          outflow = N.arcdata[abs(a)].current_lower;
        }
        else
        {
          outflow = N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower;
        }
      }
      else
      {
        if (N.arcdata[abs(a)].direction == -1)
        {
          outflow = N.arcdata[abs(a)].current_lower;
        }
        else
        {
          outflow = N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower;
        }
      }
#ifdef VERBOSE
      cout << "outflow = " << to_double(outflow) << endl;
#endif
      IntegerType alpha = outflow - carry[v];
      carry[w] += outflow;
#ifdef VERBOSE
      cout << "alpha = " << to_double(alpha) << endl;
#endif
      ST_original.update_with_DS(v, alpha, 0, height_offset);

#ifdef VERBOSE
      print_queries(N, ST_original);
#endif

#ifdef VERBOSE
      cout << "updating the current of the upper edge" << endl;
#endif
      if (a > 0 )
      {
        if (N.arcdata[abs(a)].direction == 1)
        {
          height_offset = N.arcdata[abs(a)].resistance_lower;
          outflow = - N.arcdata[abs(a)].cur_src_vw;   // ) - lower
        }
        else
        {
          height_offset = N.arcdata[abs(a)].resistance_upper;
          outflow = - N.arcdata[abs(a)].cur_src_vw;
        }
      }
      else
      {

        if (N.arcdata[abs(a)].direction == -1)
        {
          height_offset = N.arcdata[abs(a)].resistance_lower;
          outflow = - N.arcdata[abs(a)].cur_src_vw;   // ) - lower
        }
        else
        {
          height_offset = N.arcdata[abs(a)].resistance_upper;
          outflow = - N.arcdata[abs(a)].cur_src_vw;
        }
      }
#ifdef VERBOSE
      cout << "outflow = " << to_double(outflow) << endl;
#endif
      alpha = outflow;// - carry[v];
#ifdef VERBOSE
      cout << "alpha = " << to_double(alpha) << endl;
#endif
      ST_original.update_with_DS(v, alpha, 0, height_offset);
      carry[w] += outflow;
#ifdef VERBOSE
      print_queries(N, ST_original);
#endif
    }
  }


#ifdef VERBOSE
  cout << "Build Initial State Complete" << endl;
  print_queries(N, ST_original);
#endif

}
/** \brief the initial state to begin the update using datastructure is built
 *
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 *
 */
template<typename IntegerType, typename RationalType>
void build_initial_state_DS(SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original)
{
#ifdef VERBOSE

  cout << "non tree edges of G0 . " << endl;

  for (auto a : ST_original.G.non_tree_edges) {
    cout << a << " , " ;
  }
  cout << endl << endl;
#endif

  Graph<IntegerType, RationalType>& G0 = ST_original.G;
  ST_original.clear_d();

  vector<IntegerType> carry_DS(G0.no_of_verticies + 1, 0);

  vector<IntegerType> carry(G0.no_of_verticies + 1, 0);
#ifdef VERBOSE
  cout << "resistances of G0 " << endl;
  for (unsigned int i = 1; i <= G0.no_of_edges; i++) {
    cout << i << " = " << to_double(G0.resistances[i]) << " " << to_double(G0.original_auxiliary_transformed_corrospondence[i].r1_tilde) << " " <<
         to_double(G0.original_auxiliary_transformed_corrospondence[i].r2_tilde) << " " << endl;
  }
  cout << endl;

  cout << "electrical flows of G0 " << endl;
  for (unsigned int i = 1; i <= G0.no_of_edges; i++) {
    cout << i << " = " << " " << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde) << " " <<
         to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde) << " " <<
         to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow3_tilde ) << endl;
  }
  cout << endl;
#endif
  for (auto a : (G0.non_tree_edges)) {
    if (G0.original_auxiliary_transformed_corrospondence[a].is_invalid) continue;

    IntegerType alpha = 0;
    IntegerType alpha1 = 1;
    IntegerType flow1_tilde = G0.original_auxiliary_transformed_corrospondence[abs(a)].electrical_flow1_tilde;
    IntegerType flow2_tilde = G0.original_auxiliary_transformed_corrospondence[abs(a)].electrical_flow2_tilde;
    IntegerType flow3_tilde = G0.original_auxiliary_transformed_corrospondence[abs(a)].electrical_flow3_tilde;

    if (G0.original_auxiliary_transformed_corrospondence[a].r1_tilde < G0.original_auxiliary_transformed_corrospondence[a].r2_tilde) {
      alpha1 = flow3_tilde;
#ifdef VERBOSE
      cout << a << " : = " << to_double(alpha1) << " " << to_double(G.currents[3 * a]);
#endif
      alpha = - flow2_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].current_non_tree_edge = flow2_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].current_other_edge = flow1_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].resistance_non_tree_edge = G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].resistance_other_edge = G0.original_auxiliary_transformed_corrospondence[a].r1_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].non_tree_edge = 2;
#ifdef VERBOSE
      cout << "we choose r2_tilde" << endl;
#endif

      IntegerType height_offset = -G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
      IntegerType height_offset1 = -G0.original_auxiliary_transformed_corrospondence[a].r2_tilde - G0.original_auxiliary_transformed_corrospondence[a].r3_tilde;

      {
#ifdef VERBOSE
        cout << "edge_index_transformed_graph = " << edge_index_transformed_graph << endl;
        cout << "v = " << G.tail(edge_index_transformed_graph) << endl;
#endif
        node v = G0.tail(a);

        ST_original.update_with_DS(v, alpha, 0,  height_offset);

        carry[v] -= G0.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde;
      }

    }
    else {
#ifdef VERBOSE
      cout << "we choose r1_tilde" << endl;
#endif
      alpha =  - flow1_tilde;
      alpha1 = flow3_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].current_non_tree_edge = flow1_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].resistance_non_tree_edge = G0.original_auxiliary_transformed_corrospondence[a].r1_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].current_other_edge = flow2_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].resistance_other_edge = G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
      G0.original_auxiliary_transformed_corrospondence[a].non_tree_edge = 1;
#ifdef VERBOSE
      cout << a << " : = " << to_double(alpha1) << " " << to_double(G.currents[3 * a]);
#endif

      IntegerType height_offset =  -G0.original_auxiliary_transformed_corrospondence[a].r1_tilde;
      IntegerType height_offset1 =  - (G0.original_auxiliary_transformed_corrospondence[a].r1_tilde + G0.original_auxiliary_transformed_corrospondence[a].r3_tilde);

#ifdef VERBOSE
      cout << "edge_index_transformed_graph = " << edge_index_transformed_graph << endl;
      cout << "v = " << G.tail(edge_index_transformed_graph) << endl;
#endif
      node v = G0.head(a);
      ST_original.update_with_DS(v, alpha, 0,  height_offset);

      carry[v] -= G0.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde;
    }

#ifdef VERBOSE
    cout << "Query (G0) after the first round of update: " << a << endl;

    for (unsigned int i = 1; i <= G0.no_of_verticies; i++) {
      cout << i << " :  " << to_double(ST_original.query(i, 0)) << endl;
    }
    cout << "done" << endl;
#endif
  }

  for (auto a : boost::adaptors::reverse( G0.tree_edges)) {

    node v = G0.tail(a);
    node w = G0.head(a);
#ifdef VERBOSE
    cout << a << " = ( " << v << " , " << w << " ) " << endl;
    cout << ST_original.depth[v] << " " << ST_original.depth[w] << endl;
#endif
    assert(ST_original.depth[v] > ST_original.depth[w]);

    IntegerType flow1_tilde = G0.original_auxiliary_transformed_corrospondence[abs(a)].electrical_flow1_tilde;
    IntegerType flow2_tilde = G0.original_auxiliary_transformed_corrospondence[abs(a)].electrical_flow2_tilde;
    if (a <  0) {
      flow1_tilde *= -1;
      flow2_tilde *= -1;
    }
    IntegerType outflow = IntegerType(0);
    if (G0.head(a) == v) {
      outflow = flow1_tilde;
    }
    else {
      outflow = flow2_tilde;
    }
#ifdef VERBOSE
    cout << a << " outflow = " << to_double(outflow) << endl;
    cout << "carry[ " << v << to_double(carry[v]) << endl;
#endif
    IntegerType alpha = outflow - carry[v];

    IntegerType height_offset =  0;

#ifdef VERBOSE
    cout << a << "-: " << v << " " << to_double(alpha) << " " << to_double(height_offset) << " " << to_double(G0.resistances[abs(a)]) << " = "
         << to_double(G0.original_auxiliary_transformed_corrospondence[abs(a)].r1_tilde) << " + " << to_double(G0.original_auxiliary_transformed_corrospondence[abs(a)].r2_tilde) << endl;
#endif
    ST_original.update_with_DS(v, alpha, 0,  height_offset);
    carry[w] += outflow;
#ifdef VERBOSE
    cout << "carry[" << w << "] = " << to_double(carry[w]) << endl;
#endif
#ifdef VERBOSE
    cout << to_double(outflow) << " =  - " << to_double( flow2_tilde ) << " - " << to_double(alpha) << endl;
#endif
    if (G0.head(a) == v) {
      outflow = -flow2_tilde - alpha - carry[v];
    }
    else {
      outflow = -flow1_tilde - alpha - carry[v];
    }
#ifdef VERBOSE
    cout << to_double(outflow) << " - " << to_double(carry[v]) << endl;
#endif
    alpha = outflow;
#ifdef VERBOSE
    cout << "outflow = " << to_double(outflow) << endl;
#endif

    if (G0.head(a) == v) {
      height_offset = G0.original_auxiliary_transformed_corrospondence[abs(a)].r1_tilde;
    }
    else {
      height_offset = G0.original_auxiliary_transformed_corrospondence[abs(a)].r2_tilde;
    }

#ifdef VERBOSE
    cout << a << " " << v << " " << to_double(alpha) << " " << to_double(height_offset) << " " << to_double(G0.resistances[abs(a)]) << endl;
#endif
    ST_original.update_with_DS(v, alpha, 0 , height_offset);

    carry[w] += outflow;
#ifdef VERBOSE
    cout << "carry[" << w << "] = " << to_double(carry[w]) << endl;
#endif
  }

#ifdef VERBOSE
  cout << "Query (G0) after the first round of update: " << a << endl;

  for (unsigned int i = 1; i <= G0.no_of_verticies; i++) {
    cout << i << " :  " << to_double(ST_original.query(i, 0)) << endl;
  }
  cout << "done" << endl;
#endif




#ifdef VERBOSE
  cout << "build_initial_state() " << endl;
  cout << "currents ... " << endl;
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    cout << i << " : " << to_double(G.currents[i]) << endl;
  }
#endif

#ifdef VERBOSE
  cout << "init datastructure with initial currents" << endl;
  cout << "tree edges first" << endl;
#endif


  compute_tree_induced_voltages_layered_DS( ST_original);

#ifdef VERBOSE
  cout << "non-tree edges" << endl;
#endif

#ifdef VERBOSE
  for ( arc a : G.non_tree_edges ) {
    assert( a > 0 );
    const IntegerType alpha = G.currents[a];

    const IntegerType pi_head = ST.query( G.head(a) , 0 );
    const IntegerType pi_tail = ST.query( G.tail(a) , 0 );
    const IntegerType voltage_drop =  pi_tail - pi_head;

    const IntegerType ohm_violation = alpha * G.resistances[a] - voltage_drop;
    cout << a << " = (" << G.tail(a) << "," << G.head(a) << "): " << to_double(alpha) << " * " << to_double(G.resistances[a])
         << " - (" << to_double(pi_head) << " - " << to_double(pi_tail) << ") = " << to_double(ohm_violation) << endl;
  }
#endif

}


template<typename Network>
void adjust_tree_currents(Network & N)
{
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;
  for (auto a : G0.tree_edges) {
    if (N.arcdata[abs(a)].direction == 1)
    {
      //if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper)
      if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
          N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper)
      {
        N.arcdata[abs(a)].current_roof -= N.arcdata[abs(a)].current;
      }
      else
      {
        N.arcdata[abs(a)].current_lower -= N.arcdata[abs(a)].current;

      }
    }
    else
    {
      //if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper)
      if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
          N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper)
      {
        N.arcdata[abs(a)].current_roof -= N.arcdata[abs(a)].current;
      }
      else
      {
        N.arcdata[abs(a)].current_lower += N.arcdata[abs(a)].current;

      }
    }
  }

}

// checks if |r_a phi_a - (pi_w-pi_v)| \le r(C_a) for all a notin T
template<typename Network>
bool electrical_flow_termination_condition(Network & N)
{
  cout << "Checking the Electrical Flow Termination Condition" << endl;

  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;
  typedef typename Network::GraphType Graph;

  Graph& G0 = N.G;

  bool exit = true;

  RationalType LHS = 0;

  for (auto a : G0.tree_edges)
  {
    if (a < 0) a = -a;
    RationalType RHS = to_double(N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper + N.arcdata[a].resistance_roof) / 2;

    if (N.arcdata[a].resistance_roof < N.arcdata[a].resistance_lower ||
        N.arcdata[a].resistance_roof < N.arcdata[a].resistance_upper)
    {
      if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper)
      {
        cout << "Lower is the non-tree edge" << endl;

        IntegerType pi_head = N.arcdata[a].voltage_vw;
        IntegerType pi_tail(0);

        if (N.arcdata[a].direction == 1)
        {
          pi_tail = N.nodedata[G0.tails[a]].voltage;
        }
        else
        {
          pi_tail = N.nodedata[G0.heads[a]].voltage;
        }

        LHS = fabs( to_double(N.arcdata[a].resistance_lower * N.arcdata[a].current_lower + (pi_head - pi_tail)));
      }
      else
      {
        cout << "Upper is the non-tree edge" << endl;
        IntegerType pi_head = N.arcdata[a].voltage_vw;
        IntegerType pi_tail(0);

        if (N.arcdata[a].direction == 1)
        {
          pi_tail = N.nodedata[G0.heads[a]].voltage;
        }
        else
        {
          pi_tail = N.nodedata[G0.tails[a]].voltage;
        }

        IntegerType current_upper = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
        LHS = fabs(to_double(N.arcdata[a].resistance_upper * current_upper +  (pi_head - pi_tail)));
      }
    }
    else
    {
      cout << "Roof is the non-tree edge" << endl;
      IntegerType pi_head = N.nodedata[G0.heads[a]].voltage;
      IntegerType pi_tail = N.nodedata[G0.tails[a]].voltage;
      cout << "pi_head = " << to_double(pi_head) << endl;
      cout << "pi_tail = " << to_double(pi_tail) << endl;
      LHS = fabs(to_double(N.arcdata[a].resistance_roof * N.arcdata[a].current_roof + (pi_head - pi_tail)));
    }

    cout << "elec termination: " <<  LHS << " < - > " << RHS << endl;
    if (LHS > RHS)
    {
      exit = false;
    }
  }
  cout << "Non Tree Edges in the Original Graph " << endl;
  for (unsigned int i = 0; i < G0.non_tree_edges.size(); i++)
    //for(unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    cout << "i = " << i << endl;
    RationalType LHS(0);
    RationalType RHS = to_double(G0.resistances_accross_cycles[i]) / RationalType(2);
    arc a = G0.non_tree_edges[i];
    IntegerType pi_head = N.arcdata[a].voltage_vw;
    IntegerType pi_tail(0);
    if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper)
    {

      if (N.arcdata[a].direction == 1)
      {
        pi_tail = N.nodedata[G0.tails[a]].voltage;
      }
      else
      {
        pi_tail = N.nodedata[G0.heads[a]].voltage;
      }
      LHS = fabs( to_double(N.arcdata[a].resistance_lower * N.arcdata[a].current_lower + (pi_head - pi_tail)));

      cout << "elec termination: " << to_double(LHS) << " < - > " << RHS << endl;
      if (to_double(LHS) > RHS)
      {
        exit = false;
      }
    }

    else
    {
      if (N.arcdata[a].direction == 1)
      {
        pi_tail = N.nodedata[G0.heads[a]].voltage;
      }
      else
      {
        pi_tail = N.nodedata[G0.tails[a]].voltage;
      }
      IntegerType current_upper = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
      LHS = fabs(to_double(N.arcdata[a].resistance_upper * current_upper + (pi_head - pi_tail)));
      cout << "elec termination: " << to_double(LHS) << " < - > " << RHS << endl;
      if (to_double(LHS) > RHS)
      {
        exit = false;
      }
    }

    pi_head = N.nodedata[G0.heads[a]].voltage;
    pi_tail = N.nodedata[G0.tails[a]].voltage;
    LHS = fabs(to_double(N.arcdata[a].resistance_roof * N.arcdata[a].current_roof + (pi_head - pi_tail)));



    cout << "elec termination: " << to_double(LHS) << " < - > " << RHS << endl;
    if (to_double(LHS) > RHS)
    {
      exit = false;
    }

  }

  return exit;

}

template <typename Network>
void copy_currents_to_initial_currents(Network & N)
{
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++)
  {
    N.arcdata[a].initial_current_lower = N.arcdata[a].current_lower;
    N.arcdata[a].initial_current_roof = N.arcdata[a].current_roof;
  }
}

template<typename Network>
void asser_all_alphas_not_zero(Network& N)
{
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;
  
  bool all_alphas_zero = true;
  
  for(auto a: N.G.tree_edges)
  {
    if(a < 0) a = -a;
    
    if(N.arcdata[a].xlower == 0 || N.arcdata[a].infeasibility == 0 || N.arcdata[a].capacity - N.arcdata[a].xlower == 0)
    {
      continue;
    }
    
    IntegerType R_a = N.arcdata[a].resistance_lower + N.arcdata[a].resistance_roof + N.arcdata[a].resistance_upper;
    cout << "R_a = " << to_double(N.arcdata[a].resistance_lower) << " + "  << to_double(N.arcdata[a].resistance_roof) << " + " <<
			to_double(N.arcdata[a].resistance_upper) << endl;
			
    if (N.arcdata[a].resistance_roof < N.arcdata[a].resistance_lower || N.arcdata[a].resistance_roof < N.arcdata[a].resistance_upper
       ) {
      if (N.arcdata[a].resistance_lower < N.arcdata[a].resistance_upper) {
        if (N.arcdata[a].direction == 1) {
          IntegerType Delta = N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - N.arcdata[a].resistance_lower * N.arcdata[a].current_lower + N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          RationalType alpha_unrounded = to_double(Delta) / to_double(R_a);
          long long alpha_ = round(alpha_unrounded);
          IntegerType alpha = alpha_;
	  if(alpha != 0)
	  {
	    cout << "ALPHA 1 = " << to_double(alpha) << endl;
	    all_alphas_zero = false;
	  }
	  
        } else {
          IntegerType Delta = N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - N.arcdata[a].resistance_lower * N.arcdata[a].current_lower - N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          RationalType alpha_unrounded = to_double(Delta) / to_double(R_a);
          long long alpha_ = round(alpha_unrounded);
          IntegerType alpha = alpha_;
	  if(alpha != 0)
	  {
	     cout << "ALPHA 2 = " << to_double(alpha) << endl;
	    all_alphas_zero = false;
	  }
        }
      } else {
        if (N.arcdata[a].direction == 1) {
          IntegerType Delta = N.arcdata[a].resistance_lower * N.arcdata[a].current_lower - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          RationalType alpha_unrounded = to_double(Delta) / to_double(R_a);

          long long alpha_ = round(alpha_unrounded);
          IntegerType alpha = alpha_;
	  if(alpha != 0)
	  {
	     cout << "ALPHA 3 = " << to_double(alpha) << endl;
	    all_alphas_zero = false;
	  }

        } else {
          IntegerType Delta = N.arcdata[a].resistance_lower * N.arcdata[a].current_lower - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) + N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

          RationalType alpha_unrounded = to_double(Delta) / to_double(R_a);
          long long alpha_ = round(alpha_unrounded);
          IntegerType alpha = alpha_;
	  
	  if(alpha != 0)
	  {
	    cout << "ALPHA 4 = " << to_double(alpha) << endl;
	    all_alphas_zero = false;
	  }
        }
      }
    } else {
            

      if (N.arcdata[abs(a)].direction == 1 || N.arcdata[abs(a)].direction == 0) {
        IntegerType Delta =  N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) - N.arcdata[a].resistance_lower * N.arcdata[a].current_lower + N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

        RationalType alpha_unrounded = to_double(Delta) / to_double(R_a);

        long long alpha_ = round(alpha_unrounded);
        IntegerType alpha = alpha_;
	
	if(alpha != 0)
	{
	  cout << "ALPHA 6 = " << to_double(alpha) << endl;
	  all_alphas_zero = false;
	}
      } else {
        IntegerType Delta =  - N.arcdata[a].resistance_upper * (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower)
                             + N.arcdata[a].resistance_lower * N.arcdata[a].current_lower
                             + N.arcdata[a].resistance_roof * N.arcdata[a].current_roof;

        RationalType alpha_unrounded = to_double(Delta) / to_double(R_a);

        long long alpha_ = round(alpha_unrounded);
        IntegerType alpha = alpha_;
	if(alpha != 0)
	{
	  cout << "ALPHA 7 = " << to_double(alpha) << endl;
	  all_alphas_zero = false;
	}
      }
    }
  }
  
  for(unsigned int edge_index = 0; edge_index < N.G.non_tree_edges.size(); edge_index++)
  {
    const arc a = N.G.non_tree_edges[edge_index];
    const IntegerType& R_a = N.G.resistances_accross_cycles[edge_index];
    IntegerType voltage_drop = N.nodedata[N.G.tails[a]].voltage - N.nodedata[N.G.heads[a]].voltage;
    {
      if (N.arcdata[a].infeasibility != 0) {
        IntegerType Delta_roof =  N.arcdata[a].resistance_roof * N.arcdata[a].current_roof - voltage_drop;

        IntegerType R_a_roof = R_a - (N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper) + N.arcdata[a].resistance_roof;

        RationalType alpha_roof_ = -to_double(Delta_roof) / to_double(R_a_roof);
        long long alpha_ = round(alpha_roof_);
        IntegerType alpha_roof = alpha_;
	
	if(alpha_roof != 0)
	{
	   cout << "ALPHA 8 = " << to_double(alpha_roof) << endl;
	 all_alphas_zero = false; 
	}
      }
    } 
     
    {
      IntegerType r_a(0);
      IntegerType r_a_other(0);
      IntegerType f_a(0);
      IntegerType f_a_other(0);

      if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
        r_a = N.arcdata[a].resistance_lower;
        f_a = N.arcdata[a].current_lower;
        r_a_other = N.arcdata[a].resistance_upper;
        f_a_other = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
      } else {
        r_a_other = N.arcdata[a].resistance_lower;
        f_a_other = N.arcdata[a].current_lower;
        r_a = N.arcdata[a].resistance_upper;
        f_a = N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower;
      }

      IntegerType D = voltage_drop;
      IntegerType Delta = 0;
      if (N.arcdata[a].direction == 1) {
        if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
          Delta = (f_a_other * r_a_other) + D;
        } else {
          Delta =  (f_a_other * r_a_other) - D;
        }
      } else {
        if (N.arcdata[a].resistance_lower > N.arcdata[a].resistance_upper) {
          Delta = (f_a_other * r_a_other) - D;
        } else {
          Delta =  (f_a_other * r_a_other) + D;
        }
      }
      Delta = (f_a * r_a) - Delta;
      RationalType alpha_unrounded = to_double( Delta ) / to_double( R_a );

      long long alpha_ = round(alpha_unrounded);
      IntegerType alpha = alpha_;
      
      if(alpha != 0)
      {
	 cout << "ALPHA 9 = " << to_double(alpha) << endl;
	all_alphas_zero = false;
      }
    }
    
  }
  
  assert(!all_alphas_zero);
  cout << "Asserted that all alphas are not zero " << endl;
  
}




template<typename Network>
void compute_epsilon_2_norm(Network & N,
                            typename Network::IntegerType initial_cur_lower,
                            typename Network::IntegerType initial_cur_roof,
                            typename Network::RationalType mu,
                            typename Network::RationalType delta
                           )
{
  typedef typename Network::GraphType Graph;
  typedef typename Network::RationalType RationalType;
  typedef typename Network::IntegerType IntegerType;

  Graph& G0 = N.G;

  RationalType epsion_2(0);
  for (unsigned int a = 1; a <= G0.no_of_edges; a++)
  {
    IntegerType delta_x_lower = N.arcdata[a].current_lower - initial_cur_lower;
    IntegerType delta_x_upper = -delta_x_lower;
    IntegerType delta_x_roof  = N.arcdata[a].current_roof - initial_cur_roof;

    IntegerType delta_s_lower(0);
    IntegerType delta_s_upper(0);
    if (N.arcdata[a].direction == 1)
    {
      delta_s_lower = N.arcdata[a].voltage_vw - N.nodedata[G0.tails[a]].voltage;
      delta_s_upper = N.arcdata[a].voltage_vw - N.nodedata[G0.heads[a]].voltage;
    }
    else
    {
      delta_s_lower = N.arcdata[a].voltage_vw - N.nodedata[G0.heads[a]].voltage;
      delta_s_upper = N.arcdata[a].voltage_vw - N.nodedata[G0.tails[a]].voltage;
    }

    IntegerType delta_s_roof = N.nodedata[G0.heads[a]].voltage - N.nodedata[G0.tails[a]].voltage;

    RationalType epsilon_a_lower = ( to_double(N.arcdata[a].slower * delta_x_lower + N.arcdata[a].xlower * delta_s_lower
                                     + N.arcdata[a].xlower * N.arcdata[a].slower) - 1) / mu;

    IntegerType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;
    RationalType epsilon_a_upper = ( to_double(N.arcdata[a].supper * delta_x_upper + xupper * delta_s_upper
                                     + xupper * N.arcdata[a].supper) - 1) / mu;

    IntegerType xroof = N.arcdata[a].infeasibility;
    RationalType epsilon_a_roof = ( to_double(N.arcdata[a].sroof * delta_x_roof + xroof * delta_s_roof
                                    + xroof * N.arcdata[a].sroof) - 1) / mu;

    epsion_2 +=  epsilon_a_lower * epsilon_a_lower  +
                 epsilon_a_upper * epsilon_a_upper  +
                 epsilon_a_roof  * epsilon_a_roof;

  }

  // RationalType epsion_2_norm = sqrt(epsion_2);
  // cout << "epsion_2_norm = " << epsilon_2_norm << endl;
  // cout << "delta_2 / 4 = " << delta * delta / RationalType(4) << endl;
}



template<typename Network, typename IntegerType, typename RationalType>
int electrical_flow_solver(Network & N,
                           SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original,
                           RationalType Mu,
                           RationalType delta
                          )
{
  typedef typename Network::GraphType Graph;
  Graph& G0 = N.G;


  build_initial_state(N, ST_original);

  ArcSampler sample_arc_original(N, ST_original.sum_of_the_resistances_from_root_to);

  // recompute resistances accross cylces
  G0.resistances_accross_cycles.clear();
  for (unsigned int i = 0; i < G0.non_tree_edges.size(); i++) {
    G0.resistances_accross_cycles.push_back(ST_original.compute_resistance_accross(G0.non_tree_edges[i], N));
  }

  compute_tree_induced_voltages(N, ST_original);
  update_potentials(N);
  // do one initial dual step
  update_slacks_using_potentials(N);

  cout << "electrical flow solver" << endl;

  compute_tree_induced_voltages_and_update_potential(N, ST_original);
  bool exit = check_central_path_condition(N, Mu, delta);

 
  cout << "check_central_path_condition before the start of the loop = " << exit << endl;


//#ifdef VERBOSE
  print_primal_variable(N);
  print_resistances(N);
  print_currents(N);
  print_initial_currents(N);
  print_tree_induced_voltages(N);
//#endif

  if (exit) return 0;

  flow_conservation_check(N);
  unsigned int iteration_number = 0;
  // unsigned int offset = 0;
  // unsigned int nextcheck_of_gap = 0;

  IntegerType initial_cur_lower = N.arcdata[1].current_lower;
  IntegerType initial_cur_roof  = N.arcdata[1].current_roof;

  compute_epsilon_2_norm(N, initial_cur_lower, initial_cur_roof, Mu, delta);

  while (true) {

    cout << "while loop" << endl;
    
    print_primal_variable(N);
    print_dual_variable(N);
    print_tree_edges(N);
    print_non_tree_edges(N);
    
    unsigned int non_tree_edge_index =  sample_arc_original();
    
    cout << "non_tree_edge_index  = " << non_tree_edge_index << endl;
    
    bool update_roof_edge = true;

    flow_conservation_check(N);
    asser_all_alphas_not_zero(N);
    update_current(N, ST_original, non_tree_edge_index, update_roof_edge);
    cout << "Update Current Done" << endl;

    //if (iteration_number >= offset + nextcheck_of_gap)
    {
      compute_tree_currents_by_non_tree_edge_new(N, ST_original);
      adjust_tree_currents(N);

      flow_conservation_check(N);
      compute_tree_induced_voltages_and_update_potential(N, ST_original);
      primal_step(N);
      update_slacks_using_potentials(N);
      for (unsigned int a = 1; a <= G0.no_of_edges; a++) {
        N.arcdata[a].current = 0;
      }
      copy_currents_to_initial_currents(N);
      bool exit = check_central_path_condition(N, Mu, delta);

      // cout << "exit = " << exit << endl;
      // cout << "other exit = " << electrical_flow_termination_condition(N) << endl;
      //exit = electrical_flow_termination_condition(N);
      if (exit) {
        break;
      }
#ifdef VERBOSE
      cout << "Tree Induced Voltages: " << endl;
      for (unsigned int v = 1; v <= N.G.no_of_vertices; v++) {
        cout << v << ": = " << to_double(N.nodedata[v].voltage) << endl;
      }

      for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
        cout << a << ": = " << to_double(N.arcdata[a].voltage_vw) << endl;
      }

      cout << endl << "R * C" << endl;
      for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
        cout << a << ": " << to_double(N.arcdata[a].current_lower * N.arcdata[a].resistance_lower) << endl;
        cout << a << ": " << to_double(N.arcdata[a].current_roof * N.arcdata[a].resistance_roof) << " = "
             << to_double(N.arcdata[a].current_roof) << " * " << to_double(N.arcdata[a].resistance_roof) << endl;
      }
#endif
    }
    iteration_number++;
  }


  compute_tree_currents_by_non_tree_edge_new(N, ST_original);
  adjust_tree_currents(N);


  flow_conservation_check(N);
  compute_tree_induced_voltages(N, ST_original);

  return iteration_number;
}

/*template<typename Network, typename IntegerType, typename RationalType>
int electrical_flow_problem_new(Network& N,
        SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original,
        RationalType delta)
{

  cout << "electrical_flow_problem" << endl;

    typedef typename Network::GraphType Graph;
   Graph& G0 = N.G;

#ifdef VERBOSE

   cout << "Non Tree Edges:  ";
   for(auto a: G0.non_tree_edges)
   {
     cout << a << " , ";
   }
   cout << endl;

   cout << "currents before the electrical flow problem" << endl;
  for(unsigned int i =1; i <= G0.no_of_edges; i++){

    cout << i << ": " << to_double(N.arcdata[i].initial_current_lower) << " , " <<
             to_double(N.arcdata[i].cur_src_vw - N.arcdata[i].initial_current_lower) << " , " <<
       to_double(N.arcdata[i].initial_current_roof) << endl;
  }

#endif

   ArcSampler sample_arc_original(N);

   cout << "sample_arc_original(N)" << endl;

   G0.resistances_accross_cycles.clear();

   cout << "Non tree edges size: " << G0.non_tree_edges.size() << endl;

   for(unsigned int i =0; i<G0.non_tree_edges.size(); i++){
     G0.resistances_accross_cycles.push_back(ST_original.compute_resistance_accross(G0.non_tree_edges[i], N));

    #ifdef VERBOSE
    cout << i << " Resistances accross  =  " << to_double(ST_original.compute_resistance_accross(G0.non_tree_edges[i], N)) << endl;
    #endif
   }

   cout << "compute_tree_induced_voltages_layered_new" << endl;

     compute_tree_induced_voltages_layered_new(N, ST_original);
#ifdef VERBOSE
     cout << "Tree Induced done" << endl;
     for(unsigned int i = 1; i <= G0.no_of_edges; i++){

    cout << to_double(N.nodedata[G0.heads[i]].voltage) << " , " << to_double(N.nodedata[G0.tails[i]].voltage) << endl;
  }
#endif


   build_initial_state_new(N, ST_original);

   cout << endl << "Querys of Vertices: " << endl;

   for(unsigned int i = 1; i <= G0.no_of_vertices; i++)
   {
     cout << i << ": " << to_double(ST_original.query(i, 0)) << endl;
   }

   for(auto a: G0.tree_edges){
    if(a < 0) a *= -1;

    IntegerType diff_by_currents = N.arcdata[a].current_lower * N.arcdata[a].resistance_lower -
      (N.arcdata[a].cur_src_vw - N.arcdata[a].current_lower) * N.arcdata[a].resistance_upper;
    IntegerType diff_by_query = ST_original.query(G0.head(a), 0) - ST_original.query(G0.tail(a), 0);
    assert(diff_by_query - diff_by_currents == 0);
  }




   RationalType gap = compute_gap_by_non_tree_currents_new(N, ST_original);
   cout << "gap original = " << gap << endl;


     cout << gap << " < - > " << delta << endl;
     int no_of_elec_flow_itrs = 0;

     unsigned int iteration = 0;
     unsigned int nextcheck = G0.no_of_edges;
        while(true)
   //while(gap >= delta)
   {
   unsigned int non_tree_edge_index =  sample_arc_original();
   cout << gap << " < - > " << delta << endl;
   cout << "Non Tree Edge index = " << non_tree_edge_index << endl;
#ifdef VERBOSE
  cout << "printing currents before" << endl;
    for(unsigned int i = 1; i <= N.G.no_of_edges; i++){

      cout << i << ": " << to_double(N.arcdata[i].current_lower) << " , " << to_double(N.arcdata[i].cur_src_vw - N.arcdata[i].current_lower)
                     << endl;
    }
#endif

  update_current_DS_new(N, ST_original, non_tree_edge_index);

  #ifdef VERBOSE
  cout << "printing currents new" << endl;
    for(unsigned int i = 1; i <= N.G.no_of_edges; i++){

      cout << i << ": " << to_double(N.arcdata[i].current_lower) << " , " << to_double(N.arcdata[i].cur_src_vw - N.arcdata[i].current_lower)
                     << endl;
    }
#endif


  compute_tree_induced_voltages_layered_new(N, ST_original);

  gap = compute_gap_by_non_tree_currents_new(N, ST_original);

   bool should_break = true;
   for(auto a: G0.tree_edges)
   {
     cout << a << ": ";
     IntegerType Delta(0);
     if(N.arcdata[abs(a)].direction == 1)
     {
     Delta = N.arcdata[abs(a)].current_roof * N.arcdata[abs(a)].resistance_roof +
       (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper -
       N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower;


     }
     else{
       Delta = N.arcdata[abs(a)].current_roof * N.arcdata[abs(a)].resistance_roof -
       (N.arcdata[abs(a)].cur_src_vw - N.arcdata[abs(a)].current_lower) * N.arcdata[abs(a)].resistance_upper +
       N.arcdata[abs(a)].current_lower * N.arcdata[abs(a)].resistance_lower;
     }

     RationalType alpha = to_double(Delta)/
     to_double(N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper + N.arcdata[abs(a)].resistance_roof);

     cout << to_double(Delta)/
     to_double(N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper + N.arcdata[abs(a)].resistance_roof) << endl;

     if(alpha > RationalType(1) || alpha < -RationalType(1))
     {
      should_break == false;
     }

  }

  if(should_break == true && gap <= delta)
  {
    break;
  }


 //   if(iteration > nextcheck)
//    {
//   gap = compute_gap_by_non_tree_currents_new(N, ST_original);
//
//     nextcheck += G0.no_of_edges;
//   }
  iteration++;
  cout << " updated electrical flow gap = " << gap << endl;

   }



   cout << "currents before we do this compute_tree_currents_by_non_tree_edge_new" << endl;
   for(unsigned int i =1; i <= G0.no_of_edges; i++){

    cout << i << ": " << to_double(N.arcdata[i].current_lower) << " , " <<
             to_double(N.arcdata[i].cur_src_vw - N.arcdata[i].current_lower) << " , " <<
       to_double(N.arcdata[i].current_roof) << endl;
  }

   cout << endl << "Initial Currents" << endl;
   for(unsigned int i =1; i <= G0.no_of_edges; i++){

    cout << i << ": " << to_double(N.arcdata[i].initial_current_lower) << " , " <<
             to_double(N.arcdata[i].cur_src_vw - N.arcdata[i].initial_current_lower) << " , " <<
       to_double(N.arcdata[i].initial_current_roof) << endl;
  }

   cout << endl << "currents[a] " << endl;
   for(unsigned int i = 1; i <= G0.no_of_edges; i++)
   {
     cout << i << ": " << to_double(N.arcdata[i].current) << endl;
   }

   compute_tree_currents_by_non_tree_edge_new(N, ST_original);

      cout << "currents[a] later " << endl;
   for(unsigned int i = 1; i <= G0.no_of_edges; i++)
   {

     if(N.arcdata[i].direction == -1) N.arcdata[i].current *= -1;

     cout << i << ": " << to_double(N.arcdata[i].current) << endl;
   }

    for(auto a: G0.tree_edges){

      if(N.arcdata[abs(a)].direction == 1)
      {
          N.arcdata[abs(a)].current_lower += N.arcdata[abs(a)].current;
      }
      else
      {
      N.arcdata[abs(a)].current_lower -= N.arcdata[abs(a)].current;
       }
      }



//   cout << "A few asserts; " << endl;
//
//   cout << to_double(N.arcdata[3].current_roof - N.arcdata[2].current_roof - N.arcdata[2].current_lower -
//             (N.arcdata[3].cur_src_vw - N.arcdata[3].current_lower)) << endl;
//   cout << to_double(N.arcdata[3].initial_current_roof - N.arcdata[2].initial_current_roof - N.arcdata[2].initial_current_lower -
//             (N.arcdata[3].cur_src_vw - N.arcdata[3].initial_current_lower)) << endl;
//
//   cout << "3 : "; assert(N.arcdata[3].current_roof - N.arcdata[2].current_roof - N.arcdata[2].current_lower -
//             (N.arcdata[3].cur_src_vw - N.arcdata[3].current_lower) ==
//       N.arcdata[3].initial_current_roof - N.arcdata[2].initial_current_roof - N.arcdata[2].initial_current_lower -
//             (N.arcdata[3].cur_src_vw - N.arcdata[3].initial_current_lower) );
//
//
//    cout << endl;
//    cout << to_double(N.arcdata[1].current_roof - (N.arcdata[1].cur_src_vw - N.arcdata[1].current_lower)
//                       +  N.arcdata[2].current_roof - (N.arcdata[2].cur_src_vw - N.arcdata[2].current_lower)
//     ) << endl;
//   cout << to_double(N.arcdata[1].initial_current_roof - (N.arcdata[1].cur_src_vw - N.arcdata[1].initial_current_lower)
//                  + N.arcdata[2].initial_current_roof - (N.arcdata[2].cur_src_vw - N.arcdata[2].initial_current_lower)
//   ) << endl;
//
//    assert(
//  N.arcdata[1].current_roof - (N.arcdata[1].cur_src_vw - N.arcdata[1].current_lower)
//                       -  N.arcdata[2].current_roof - (N.arcdata[2].cur_src_vw - N.arcdata[2].current_lower)
//          ==
//           N.arcdata[1].initial_current_roof - (N.arcdata[1].cur_src_vw - N.arcdata[1].initial_current_lower)
//                  - N.arcdata[2].initial_current_roof - (N.arcdata[2].cur_src_vw - N.arcdata[2].initial_current_lower)
//        );
  cout << "currents after the electrical flow problem" << endl;
  for(unsigned int i =1; i <= G0.no_of_edges; i++){


    cout << i << ": " << to_double(N.arcdata[i].current_lower) << " , " <<
             to_double(N.arcdata[i].cur_src_vw - N.arcdata[i].current_lower) << " , " <<
       to_double(N.arcdata[i].current_roof) << endl;
  }

  cout << endl << "Resistances: " << endl;
  for(unsigned int i = 1; i <= G0.no_of_edges; i++)
  {
    cout << i << ": " << to_double(N.arcdata[i].resistance_lower) << " , " << to_double(N.arcdata[i].resistance_upper) << " , "
   << to_double(N.arcdata[i].resistance_roof) << endl;
  }
  compute_tree_induced_voltages_layered_new(N, ST_original);


  cout << endl << "Tree Induced voltage: " << endl;
  for(unsigned int i = 1; i <= G0.no_of_edges; i++)
  {
    cout << i << " " << to_double(N.nodedata[G0.heads[i]].voltage) << " , " << to_double(N.nodedata[G0.tails[i]].voltage) << " , " <<
          to_double(N.arcdata[i].voltage_vw) << endl;
   }

   cout << "Query nodes " << endl;
   for(unsigned int v =1; v <= G0.no_of_vertices; v++)
   {
     cout << v << ": " << to_double(ST_original.query(v, 0)) << endl;
   }

   cout << endl << "Query of non tree edges: " << endl;
   for(unsigned int i = 0; i < G0.non_tree_edges.size(); i++)
   {
     cout << i << ": " << to_double(ST_original.query(i)) << endl;
   }


   cout << "End of Electrical Flow problem" << endl;




    return no_of_elec_flow_itrs;
}*/

/** \brief Runs the electrical flow problem
 *
 * @param ST The spanning tree which was selected for the electrical-flow-problem
 * @param delta
 *
 * @see build_initial_state(ST)
 * @see compute_gap(G)
 *
 */
// template<typename IntegerType, typename RationalType>
// int electrical_flow_problem_DS(SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original,
//             Graph<IntegerType, RationalType>& G0,
//             const RationalType& delta, RationalType& TCN)
// {
//
//   #ifdef VERBOSE
//   cout << "electrical_flow_problem   " << "delta: " << delta <<endl;
//   #endif
//
//
//   #ifdef VERBOSE
//   cout << "G.no_of_edges: " <<  G.no_of_edges << endl << "Number of Tree Edges: " << G.tree_edges.size() << endl << "Non Tree Edges: " << G.non_tree_edges.size() << endl;
//   #endif
//
//   #ifdef VERBOSE
//   cout << endl << "D_EXT and D_DROP " << endl;
//   for(unsigned int v = 0; v < ST.d_ext.size(); v++){
//     cout << to_double(ST.d_ext[v]) << " = " << to_double(ST.d_drop[v]) << endl;
//   }
//
//   cout << endl << "electrical_flow_problem " << endl;
//   cout<< "Non tree edges size: " << G.non_tree_edges.size();
//   cout<< endl;
//   for(auto a: G.non_tree_edges){
//     cout<< a << ", ";
//   }
//   cout << endl;
//   cout << endl << "Tree Edges Size: " << G.tree_edges.size();
//   cout<< endl;
//   cout << "Number of Verticies: " << G.no_of_verticies << endl;
//   cout << "Number of Edges: " << G.no_of_edges << endl;
//   #endif
//
//
//   G0.resistances_accross_cycles.clear();
//   for(unsigned int i =0; i<G0.non_tree_edges.size(); i++){
//     G0.resistances_accross_cycles.push_back(ST_original.compute_resistance_accross(G0.non_tree_edges[i]));
//
//     #ifdef VERBOSE
//     cout << i << " Resistances accross  =  " << to_double(ST.compute_resistance_accross(G.non_tree_edges[i])) << endl;
//     #endif
//   }
//
//   ArcSampler sample_arc_original(G0);
//
//   #ifdef VERBOSE
//   cout << "arc sampler called.." << endl;
//   cout << "currents: " << endl;
//
//   for(unsigned int i=1; i<=G.no_of_edges; i++){
//     cout << i << " : " << to_double(G.currents[i]) << "  " << G.head(i) << " , " << G.tail(i) << endl;
//   }
//
//   cout << "demands: " << endl;
//   for(unsigned int i=1; i<=G.no_of_verticies; i++){
//     cout << i << " : " << G.demands[i] << endl;
//   }
//
//   cout << endl << "Node Tree Corrospondance" << endl << endl;
//
//   for(unsigned int i=1; i<=G.no_of_verticies; i++){
//     cout << i << " : ";
//     for(auto tree_height: ST.node_tree_corrospondance[i]){
//       cout << " [ " << tree_height.tree_index <<" , " << tree_height.above_d << " , " << " , " << tree_height.lca << " ] " << " , " ;
//     }
//     cout << endl;
//   }
//
//   cout << endl << "Tree Decomposition Information" << endl << endl;
//
//   unsigned int count = 0;
//   for(auto TDI: ST.tree_decomposition){
//     cout << count << " : " ;
//     cout << " [ " << TDI.size <<" , " << TDI.s << " , " << TDI.d << " ] " << " , " ;
//     count++;
//     cout << endl;
//   }
//   #endif
//
//
//   #ifdef VERBOSE
//   cout << endl << "Initial state would be built " << endl;
//   #endif
//   #ifdef VERBOSE
//   cout << "depths G0" << endl;
//   for(unsigned int i = 1; i <= G0.no_of_verticies; i++){
//     cout << i << " = " << ST_original.depth[i] << endl;
//   }
//
//   cout << "depths " << endl;
//   for(unsigned int i = 1; i <= G.no_of_verticies; i++){
//     cout << i << " = " << ST.depth[i] << endl;
//   }
//   #endif
//
//   #ifdef VERBOSE
//   cout << "current comparison before build_initial_state" << endl;
//
//   for(unsigned int i = 1; i <= G0.no_of_edges; i++){
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde) << " ";
//     cout << to_double(G.currents[3*(i-1) + 1]) << endl;
//
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde) << " ";
//     cout << to_double(G.currents[3*(i-1) + 2]) << endl;
//
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow3_tilde) << " ";
//     cout << to_double(G.currents[3*(i-1) + 3]) << endl;
//   }
//
//   cout << endl << endl;
//   #endif
//
//   build_initial_state_DS(ST_original);
//
//   #ifdef VERBOSE
//   cout << "build_initial_state_DS done" << endl;
//   #endif
//
//   cout << "Query: " << endl;
//
//   #ifdef VERBOSE
//   cout << "QUERY G0" << endl;
//   for(unsigned int v = 1; v <= G0.no_of_verticies; v++){
//     cout << v << " " << to_double(ST_original.query(v, 0)) << " <-> " << to_double(G0.tree_induced_voltages[v]) <<  endl;
//   }
//   #endif
//
//   #ifdef VERBOSE
//   cout << endl << endl;
//
//   cout << "QUERY G" << endl;
//   for(unsigned int v = 1; v <= G.no_of_verticies; v++){
//     cout << v << " " << to_double(ST.query(v, 0)) << " <-> " << to_double(G.tree_induced_voltages[v]) << endl;
//   }
//
//   cout << "QUERY comparison " << endl;
//
//   for(unsigned int v = 1; v < G0.no_of_verticies; v++){
//     cout << to_double(ST_original.query(v, 0) - ST_original.query(v + 1, 0)) << " < - > " << to_double(ST.query(v, 0) - ST.query(v + 1, 0) ) << endl;
//   }
//
//   cout << "Query with edge index" << endl << endl;
//
//   for(unsigned int i = 0; i < G0.non_tree_edges.size(); i++){
//     arc a = G0.non_tree_edges[i];
//     node v = G0.head(a);
//     node w = G0.tail(a);
//     cout << to_double(ST_original.query(i)) << " =  " << to_double(ST_original.query(w,0)) << " - " << to_double(ST_original.query(v,0)) << " = " << to_double(ST_original.query(w,0) - ST_original.query(v,0)) << endl;
//   }
//   cout << endl << endl;
//
//   cout << "Query with edge index of ST and ST_original" << endl << endl;
//
//   for(unsigned int i = 0; i < G0.non_tree_edges.size(); i++){
//     arc a = G0.non_tree_edges[i];
//     node v = G0.head(a);
//     node w = G0.tail(a);
//     cout << "("<< v << " , " << w << ") " << endl;
//     cout << "non_tree_edge_common: " << ST_original.non_tree_edge_common_tree_index[i] << " " << ST.non_tree_edge_common_tree_index[i] << endl;
//     cout << to_double(ST_original.query(i)) << " =  " << to_double(ST.query(i)) << " = " << to_double(ST.query(v,0) - ST.query(w,0)) << endl;
//   }
//   cout << endl << endl;
//
//   cout << "Ohm's Law Check" << endl;
//   #endif
//   for(auto a: G0.tree_edges){
//     if(a < 0) a *= -1;
//     IntegerType diff_by_currents = G0.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde * G0.original_auxiliary_transformed_corrospondence[a].r1_tilde -
//     G0.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde * G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
//     IntegerType diff_by_query = ST_original.query(G0.head(a), 0) - ST_original.query(G0.tail(a), 0);
//     assert(diff_by_query - diff_by_currents == 0);
//   }
//
//
//   #ifdef VERBOSE
//   cout << "resistances_accross_cycles" << endl;
//   for(unsigned int i = 0; i < G0.non_tree_edges.size(); i++){
//     cout << to_double(ST.G.resistances_accross_cycles[i]) << " <-> " << to_double(ST_original.G.resistances_accross_cycles[i]) << endl;
//   }
//   #endif
//
//   #ifdef VERBOSE
//   cout << endl << "Initial state built " << endl;
//   #endif
//
//
//   RationalType gapcheck = compute_gap_by_non_tree_currents_DS( G0, ST_original);
//
//   #ifdef VERBOSE
//   RationalType gap = gapcheck;
//   cout << "gap and gapcheck ="  << gap << "<->" << gapcheck << endl;
//   #endif
//
//   #ifdef VERBOSE
//   cout << endl << " assertion done ";
//   #endif
//
//   unsigned int iteration = 0;
//   unsigned int nextcheck = G0.no_of_edges;
//   unsigned int offset = 0;
//
//   #ifdef VERBOSE
//   cout << endl << "while begins now" << endl;
//   cout << gap << " " << delta << endl;
//   #endif
//
//
//   #ifdef VERBOSE
//   IntegerType r_max(0);
//   for(unsigned int  i = 1; i <= G.no_of_edges; i++){
//     if(r_max < G.resistances[i]){
//       r_max = G.resistances[i];
//     }
//   }
//   RationalType non_tree_stretch(0);
//
//   for(unsigned int i = 0; i < G.non_tree_edges.size() ; i++){
//     non_tree_stretch += to_double(G.resistances[G.non_tree_edges[i]]) / to_double(G.resistances_accross_cycles[i]) ;
//   }
//   #endif
//
//
//   #if  defined(VERBOSE) || defined(ElectricalFlowGap)
//   cout << gapcheck  << " < - > " << delta << endl;
//   #endif
//
//   while(gapcheck > delta){
//     iteration++;
//     #ifdef VERBOSE
//     cout << "checking if the currents are integral" << endl;
//     cout << "Edge" << "\t" << "Currents " << "\t" << "Resistances" << "\t" << "Batteries" << endl;
//     for(unsigned int i = 1; i <= G.no_of_edges; i++){
//       cout << i << "\t" << to_double(G.currents[i]) << "\t" << to_double(G.resistances[i]) << "\t" << to_double(G.batteries[i]) << endl;
//     }
//     cout << endl;
//     #endif
//
//     #ifdef VERBOSE
//     cout << endl << "while loop for electrical_flow_problem has begun....  gap = " << gap << endl;
//     #endif
//
//     #ifdef VERBOSE
//     cout << "Graph: " << endl;
//     for(unsigned int a =1; a<= G.no_of_edges; a++){
//       cout << a << " =  (" << G.head(a) << " , " << G.tail(a) << ")" << endl;
//     }
//     #endif
//
//     #ifdef VERBOSE
//     cout << endl << "Tree: " << endl;
//     for(unsigned int v = 1; v <= G.no_of_verticies; v++){
//       cout << v << " : " ;
//       for(auto a: ST.tree_incident_edges[v]){
//  cout << a << " , " << endl;
//       }
//     }
//     cout << endl;
//     #endif
//
//
//     unsigned int non_tree_edge_index_original = sample_arc_original();
//
//
//     #ifdef VERBOSE
//     cout << endl << " Current before... " << endl;
//     for(unsigned int i=1; i <= G.no_of_edges; i++){
//       cout << i << " : " << to_double(G.currents[i]) << endl;
//     }
//     #endif
//
//
//     #ifdef VERBOSE
//     cout << endl << "current updation will take place now with non_tree_edge_index: " << non_tree_edge_index_original << "  edge is  " << ST_original.G.non_tree_edges[non_tree_edge_index_original] <<
//     endl << "resistance_non_tree_edge = " << to_double(G0.original_auxiliary_transformed_corrospondence[ST_original.G.non_tree_edges[non_tree_edge_index_original]].resistance_non_tree_edge) << endl;
//     #endif
//
//
//     update_current_DS_DS(ST_original, non_tree_edge_index_original);
//     #ifdef VERBOSE
//     cout << endl << " Current after... " << endl;
//     for(unsigned int i=1; i <= G.no_of_edges; i++){
//       cout << i << " : " << to_double(G.currents[i]) << endl;
//     }
//     #endif
//
//
//
//     #ifdef VERBOSE
//     cout << endl << "Current updated " << endl;
//     #endif
//
//     if( iteration > nextcheck + offset ){
//
//       gapcheck = compute_gap_by_non_tree_currents_DS( G0, ST_original);
//
//       #ifdef ElectricalFlowGap
//       cout << gapcheck << endl;
//       #endif
//
//       #ifdef VERBOSE
//
//       cout << "gapcheck = " << gapcheck << endl;
//       #endif
//
//       nextcheck += G0.no_of_edges;
//       #ifdef VERBOSE
//       cout << "electrical_flow_problem gap: " << gapcheck << " <-> " << delta << endl;
//       #endif
//
//       #ifdef VERBOSE
//
//
//       Graph<IntegerType, RationalType>& G0 = ST_original.G;
//
//       vector<RationalType> unrounded_alphas;
//       bool exit = true;
//       #ifdef VERBOSE
//       cout << "non tree edge size = "<< G0.non_tree_edges.size() << endl;
//       #endif
//       for(unsigned int edge_index1 = 0; edge_index1 < G0.non_tree_edges.size(); edge_index1++){
//  #ifdef VERBOSE
//  cout << edge_index1 << " " << G0.non_tree_edges.size() << endl;
//  #endif
//  const arc a = G0.non_tree_edges[edge_index1];
//
//
//
//
//  #ifdef VERBOSE
//  cout << endl << "Edge Index Sampled: (with DS) " << edge_index1 << endl;
//  cout << "the non tree edge: (with DS) " << G0.non_tree_edges[edge_index1] << endl;
//  #endif
//
//  const IntegerType& R_a = G0.resistances_accross_cycles[edge_index1];
//  IntegerType voltage_drop = ST_original.query(edge_index1); // v_u - v_w;
//
//  if(G0.original_auxiliary_transformed_corrospondence[a].non_tree_edge == 1){
//    IntegerType f_a = G0.original_auxiliary_transformed_corrospondence[a].current_other_edge;
//    const IntegerType& r_a =  G0.original_auxiliary_transformed_corrospondence[a].r2_tilde;
//    if(f_a != 0 && r_a != 0){
//      multiply_by_bit_shift(f_a, r_a);
//    }
//    else{
//      f_a = 0;
//    }
//    voltage_drop = f_a - voltage_drop;
//  }
//  if(G0.original_auxiliary_transformed_corrospondence[a].non_tree_edge == 2){
//    IntegerType f_a = G0.original_auxiliary_transformed_corrospondence[a].current_other_edge;
//    const IntegerType& r_a =  G0.original_auxiliary_transformed_corrospondence[a].r1_tilde;
//    if(f_a != 0 && r_a != 0){
//      multiply_by_bit_shift(f_a, r_a);
//    }
//    else{
//      f_a = 0;
//    }
//    voltage_drop = f_a + voltage_drop;
//  }
//  #ifdef VERBOSE
//  cout << "Voltage Drop: " << to_double(voltage_drop) << endl;
//  #endif
//
//
//  IntegerType f_a = G0.original_auxiliary_transformed_corrospondence[a].current_non_tree_edge;
//  const IntegerType& r_a = G0.original_auxiliary_transformed_corrospondence[a].resistance_non_tree_edge;
//  if(f_a != 0 && r_a != 0){
//    multiply_by_bit_shift(f_a, r_a);
//  }
//  else{
//    f_a = 0;
//  }
//  const IntegerType& Delta = voltage_drop - f_a;
//
//  #ifdef VERBOSE
//  cout << "Delta : " << to_double(Delta) << endl;
//  cout << "R_a: " << to_double(R_a) << endl;
//  cout << to_double(G0.original_auxiliary_transformed_corrospondence[a].current_non_tree_edge) << " * " <<  to_double(G0.original_auxiliary_transformed_corrospondence[a].resistance_non_tree_edge) << endl;
//  #endif
//  RationalType alpha_unrounded = -to_double( Delta )/ to_double( R_a );
//  unrounded_alphas.push_back(alpha_unrounded);
//
//  if(alpha_unrounded >= RationalType(1)/RationalType(2) || alpha_unrounded < -RationalType(1)/RationalType(2)){
//    exit = false;
//    break;
//  }
//
//  #ifdef VERBOSE
//  cout << endl << "Delta = " << to_double(Delta) << endl;
//  cout << endl << " alpha_unrounded = " << alpha_unrounded << endl;
//  #endif
//       }
//
//       if(exit == false){
//  #ifdef VERBOSE
//  cout << endl << "does not exit through this route. " << endl;
//  #endif
//       }
//
//       if(exit == true){
//  RationalType val(0);
//
//  for(unsigned int i = 0; i < G0.non_tree_edges.size(); i++){
//    arc a = G0.non_tree_edges[i];
//    IntegerType R_a = G0.resistances_accross_cycles[i];
//    IntegerType r_a = G0.resistances[a];
//
//    val += to_double(R_a) * to_double(R_a)/to_double(r_a);
//  }
//
//  cout << "unrounded alpha " << endl;
//  for(auto a : unrounded_alphas){
//    cout << a << endl;
//  }
//
//  if(RationalType(4) * gapcheck >  val){
//    cout << gapcheck << "  " << val << endl;
//
//    for(auto a : unrounded_alphas){
//      cout << a << endl;
//    }
//  }
//
//  assert (RationalType(4) * gapcheck <=  val);
//
//  #ifndef NDEBUG
//  cout << "returned exit: " << endl;
//  #endif
//  return iteration;
//       }
//       #endif
//
//     }
//   }
//   #ifdef VERBOSE
//   cout << "number of electrical iterations: " << iteration << endl;
//   #endif
//   #ifdef VERBOSE
//   cout << "gap = " << gap << endl << "gapcheck = " << gapcheck << endl;
//   #endif
//
//
//   #ifdef VERBOSE
//   cout << "sum of alphas = " << endl;
//   for(auto a: G0.non_tree_edges){
//     cout << a << " = " << to_double(G0.original_auxiliary_transformed_corrospondence[a].current) << endl;
//   }
//   cout << endl ;
//   #endif
//   compute_tree_currents_by_non_tree_edge_DS(ST_original );
//
//   #ifdef VERBOSE
//   for(unsigned int a = 1; a <= G0.no_of_edges; a++){
//     cout <<a << " = " << to_double(G0.original_auxiliary_transformed_corrospondence[a].current) << endl;
//   }
//   #endif
//
//   for(unsigned int a = 1; a <= G0.no_of_edges; a++){
//     G0.original_auxiliary_transformed_corrospondence[a].electrical_flow2_tilde -= G0.original_auxiliary_transformed_corrospondence[a].current;
//     G0.original_auxiliary_transformed_corrospondence[a].electrical_flow1_tilde += G0.original_auxiliary_transformed_corrospondence[a].current;
//   }
//   compute_tree_induced_voltages_layered_DS(ST_original);
//
//   #ifdef VERBOSE
//   cout << "check TVL: " << endl;
//   for(unsigned int i =1; i <= G0.no_of_edges; i++){
//     cout << i << " = " ;
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].tree_induced_voltage_w) << " , ";
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].tree_induced_voltage_u) << " , ";
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].tree_induced_voltage_vw) << " , ";
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].tree_induced_voltage_v) << " , ";
//     cout << endl;
//   }
//
//   cout << endl;
//   for(auto a: G0.tree_edges){
//     cout << a << ", " ;
//   }
//   cout << endl;
//   for(unsigned int a = 1; a <= G0.no_of_edges; a++){
//     cout << a << " = " << "( " << to_double(G.tree_induced_voltages[G0.head(a)]) - to_double(G.tree_induced_voltages[G0.tail(a)]) << ") " << endl;
//     cout << a << " = " << "( " << to_double(G0.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_w) -to_double(G0.original_auxiliary_transformed_corrospondence[a].tree_induced_voltage_v)
//     << ") " << endl << endl;
//   }
//   #endif
//
//   #ifdef VERBOSE
//   cout << "current comparison" << endl;
//
//   cout << "G's non_tree_edges" << endl;
//   for(auto a: G.non_tree_edges){
//     cout << a << " , ";
//   }
//   cout << endl << endl;
//
//   for(unsigned int i = 1; i <= G0.no_of_edges; i++){
//
//     cout << "the non tree edge is " <<  G0.original_auxiliary_transformed_corrospondence[i].non_tree_edge << endl;
//     cout << i << endl << "------------" << endl;
//
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde) << " ";
//     cout << to_double(G.currents[3*(i-1) + 1]) << " = " << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow1_tilde) - to_double(G.currents[3*(i-1) + 1]) << endl;
//
//     cout << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde) << " ";
//     cout << to_double(G.currents[3*(i-1) + 2]) << " = " << to_double(G0.original_auxiliary_transformed_corrospondence[i].electrical_flow2_tilde) - to_double(G.currents[3*(i-1) + 2]) << endl;
//   }
//
//   cout << endl << endl;
//
//   assert( imbalances_check(ST) );
//   #endif
//
//
//   return iteration;
// }