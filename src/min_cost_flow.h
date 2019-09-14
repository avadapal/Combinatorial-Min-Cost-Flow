#ifndef MIN_COST_FLOW_H
#define MIN_COST_FLOW_H

using namespace std;


template <
  typename Network,
  typename IntegerType
  >
void compute_BFS_solution(
  Network& N,
  node s,
  vector<IntegerType>&z
) {
  const auto& G = N.G;

  typedef typename Network:: RationalType RationalType;
  vector<IntegerType>(G.no_of_edges + 1, 0).swap( z );
  vector<RationalType> b( G.no_of_vertices + 1 );

  for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
    b[i] = to_double(N.nodedata[i].demand);
    cout << "b[i] = " << to_double(b[i]) << " " << to_double(N.nodedata[i].demand) << endl;
  }

  vector<int> visited( G.no_of_vertices + 1, 0 );
  vector<arc> bfstree( G.no_of_vertices + 1, 0 );
  vector<arc> non_tree;//(G.no_of_edges+1,0);
  vector<bool> is_tree_edge(G.no_of_edges + 1, 0);
  vector<bool> is_non_tree(G.no_of_edges + 1, 0);
  deque<node> order;
  deque<node> Q;
  Q.push_back( s );
  visited[s] = 1;


  while ( !Q.empty() ) {
    const node v = Q.front();
    for ( auto a : G.incident_edges[v] ) {
      node w;
      if ( a > 0 ) {
        assert( G.tails[a] == v );
        w = G.heads[a];
      } else {
        assert( G.heads[-a] == v );
        w = G.tails[-a];
      }
      if ( !visited[w] ) {
        Q.push_back( w );
        visited[w] = 1;
        bfstree[w] = a;
        is_tree_edge[abs(a)] = true;
      }
      if (visited[w] && !is_tree_edge[abs(a)] && !is_non_tree[abs(a)]) {
        is_non_tree[abs(a)] = true;
        non_tree.push_back(abs(a));
      }
    }
    Q.pop_front();
    order.push_front( v );
  }

  // print tree
  cout << "BFS tree arcs:" << endl;
  for (auto a : bfstree) {
    cout << a << " , " ;
  }
  cout << endl;
  cout << "Non tree arcs: " << endl;
  for (auto a : non_tree) {
    cout << a << " , " ;
  }
  cout << endl;





#ifndef NDEBUG
  const node r = order.back();
//#ifdef VERBOSE
  cout << "root r: " << r << endl;
  cout << "b[r] = " << b[r] << endl;
//#endif
#endif

  assert( bfstree[r] == 0 );
  order.pop_back();
  for ( auto v : order ) {

    cout << "v = " << v << endl;
    const arc a = bfstree[v];

    const RationalType delta = b[v];
    cout << "b[v] before = " << v << ", " << to_double(b[v]) << endl;
    cout << "b[r] before =" << r << " , " << to_double(b[r]) << endl;
    N.pull_flow( -a, delta , b );
    cout << "b[v] after = " << v << ", " << to_double(b[v]) << endl;
    cout << "b[r] after = " << r << " , " << to_double(b[r]) << endl;
    assert( RationalType(1e10)*fabs(to_double(b[v])) < RationalType(1) );

    const arc aa = abs( a );
    long long z_ = N.arcdata[aa].flow;
    z[aa] = z_;
  }

  if(b[r] != 0) cout << "b[r] = " << b[r] << endl;
  assert( b[r] == 0 );

#ifdef VERBOSE
  cout << endl << "The z vector" << endl;
  G.print_arc_data(z);
#endif

#ifndef NDEBUG
//    primal_sanity_check( N, z );
#endif
}






/** \brief Computes the BFS solution such that the non-tree edges have flow za/2
 *
 * @param G Graph
 * @param s The Dual Slacks
 * @param z The BFS solution
 *
 */
template <
  typename Network,
  typename RationalType
  >
void compute_BFS_solution_new_initialization(
  Network& N,
  node s,
  vector<RationalType>&z
) {
  const auto& G = N.G;

  vector<RationalType>(G.no_of_edges + 1, 0).swap( z );
  vector<RationalType> b( G.no_of_vertices + 1 );

  for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
    b[i] = to_double(N.nodedata[i].demand);
    cout << "b[i] = " << b[i] << endl;
  }

  vector<int> visited( G.no_of_vertices + 1, 0 );
  vector<arc> bfstree( G.no_of_vertices + 1, 0 );
  vector<arc> non_tree;//(G.no_of_edges+1,0);
  vector<bool> is_tree_edge(G.no_of_edges + 1, 0);
  vector<bool> is_non_tree(G.no_of_edges + 1, 0);
  deque<node> order;
  deque<node> Q;
  Q.push_back( s );
  visited[s] = 1;


  while ( !Q.empty() ) {
    const node v = Q.front();
    for ( auto a : G.incident_edges[v] ) {
      node w;
      if ( a > 0 ) {
        assert( G.tails[a] == v );
        w = G.heads[a];
      } else {
        assert( G.heads[-a] == v );
        w = G.tails[-a];
      }
      if ( !visited[w] ) {
        Q.push_back( w );
        visited[w] = 1;
        bfstree[w] = a;
        is_tree_edge[abs(a)] = true;
      }
      if (visited[w] && !is_tree_edge[abs(a)] && !is_non_tree[abs(a)]) {
        is_non_tree[abs(a)] = true;
        non_tree.push_back(abs(a));
      }
    }
    Q.pop_front();
    order.push_front( v );
  }

  // print tree
  cout << "BFS tree arcs:" << endl;
  for (auto a : bfstree) {
    cout << a << " , " ;
  }
  cout << endl;
  cout << "Non tree arcs: " << endl;
  for (auto a : non_tree) {
    cout << a << " , " ;
  }
  cout << endl;


  // set imbalances due to non tree flow
  for (auto a : non_tree)
  {
    z[a] = N.arcdata[a].capacity / 2;
    b[G.heads[a]] -= N.arcdata[a].capacity / 2;
    b[G.tails[a]] += N.arcdata[a].capacity / 2;
  }


#ifndef NDEBUG
  const node r = order.back();
#ifdef VERBOSE
  cout << "root r: " << r << endl;
#endif
#endif

  assert( bfstree[r] == 0 );
  order.pop_back();
  for ( auto v : order ) {

    const arc a = bfstree[v];

    const RationalType delta = b[v];
    cout << "b[v] before = " << b[v] << endl;
    N.pull_flow( -a, delta , b );
    cout << "b[v] after = " << b[v] << endl;
    assert( RationalType(1e10)*fabs(b[v]) < RationalType(1) );

    const arc aa = abs( a );
    z[aa] = N.arcdata[aa].flow;
  }

  assert( b[r] == 0 );

#ifdef VERBOSE
  cout << endl << "The z vector" << endl;
  G.print_arc_data(z);
#endif

#ifndef NDEBUG
//    primal_sanity_check( N, z );
#endif
}

/** \brief Computes a Initial Feasible solution on a BFS Tree
 *
 * @param G The Graph
 * @param s The Dual Slacks
 * @param z Feasible Initial Solution
 *
 */
template<typename IntegerType, typename RationalType>
void compute_BFS_solution(Graph<IntegerType, RationalType> &G,
                          node s,
                          vector<RationalType>&z
                         )
{
  vector<RationalType>(G.no_of_edges + 1, 0).swap( z );

  vector<RationalType> b( G.demands );


#ifdef VERBOSE
  cout << endl << endl;
  cout << "-----------------DEMANDS--------------------" << std::endl;
  for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
    cout << "b[" << i << "]: ";
    cout << b[i] << " , " << endl;
  }
  cout << endl;
#endif

  vector<int> visited( G.no_of_vertices + 1, 0 );
  vector<arc> bfstree( G.no_of_vertices + 1, 0 );
  vector<arc> non_tree;//(G.no_of_edges+1,0);
  vector<bool> is_tree_edge(G.no_of_edges + 1, 0);
  vector<bool> is_non_tree(G.no_of_edges + 1, 0);
  deque<node> order;
  deque<node> Q;
  Q.push_back( s );
  visited[s] = 1;

  while ( !Q.empty() ) {
    const node v = Q.front();
    for ( auto a : G.incident_edges[v] ) {
      node w;
      if ( a > 0 ) {
        assert( G.tails[a] == v );
        w = G.heads[a];
      } else {
        assert( G.heads[-a] == v );
        w = G.tails[-a];
      }
      if ( !visited[w] ) {
        Q.push_back( w );
        visited[w] = 1;
        bfstree[w] = a;
        is_tree_edge[abs(a)] = true;
      }
      if (visited[w] && !is_tree_edge[abs(a)] && !is_non_tree[abs(a)]) {
        is_non_tree[abs(a)] = true;
        non_tree.push_back(abs(a));
      }
    }
    Q.pop_front();
    order.push_front( v );
  }

  cout << "BFS tree" << endl;
  for (auto a : bfstree) {
    cout << a << " , " ;
  }
  cout << endl;
  cout << "Non Tree Arcs " << endl;
  for (auto a : non_tree)  {
    cout << a << " , " ;
  }

  cout << endl;

#ifndef NDEBUG
  const node r = order.back();
#ifdef VERBOSE
  cout << "r: " << r << endl;
#endif
#endif
  assert( bfstree[r] == 0 );
  order.pop_back();
  for ( auto v : order ) {

#ifdef VERBOSE
    cout << "v: " << v << endl;
#endif
    const arc a = bfstree[v];
#ifdef VERBOSE
    cout << "demand before: b[" << v << "] = " << b[v] << endl;
#endif

#ifdef VERBOSE
    cout << "flow being pulled: " << b[v] << endl;
#endif
    const RationalType delta = b[v];
    G.pull_flow( -a, delta , b );
#ifdef VERBOSE
    cout << "demand after: b[" << v << "] = " << b[v] << endl;
#endif
    assert( RationalType(1e10)*fabs(b[v]) < RationalType(1) );
    const arc aa = abs( a );
    z[aa] = G.flows[aa];
  }

#ifdef VERBOSE
  cout << "b[" << r << "] = " << b[r] << endl;
#endif
  assert( b[r] == 0 );

#ifdef VERBOSE
  cout << endl << "The z vector" << endl;
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    cout << "z[" << i << "]: " << z[i] << endl;
  }
#endif

#ifndef NDEBUG
  primal_sanity_check( G, z );
#endif
}

/** \brief Finds a tree solution formed by the rounding scheme
 *
 * @param G The Original Graph
 * @param s The Dual Slacks
 *
 */
template <
  typename Network
  >
vector<typename Network::RationalType> find_rounding_scheme_tree_solution(
  Network& N,
  node s
) {
  typedef typename Network::GraphType Graph;
//  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;
  const Graph& G = N.G;

  vector<RationalType> z (G.no_of_edges, 0);
  vector<RationalType>(G.no_of_edges + 1, 0).swap( z );
  vector<RationalType> b( G.demands );
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    G.flows[i] = 0;
  }

#ifdef VERBOSE
  cout << "-----------------DEMANDS--------------------" << std::endl;
  for (unsigned int i = 1; i <= G.no_of_vertices; i++)
  {
    cout << "b[" << i << "]: ";
    cout << b[i] << " , " << endl;
  }
  cout << endl;
#endif

  vector<int> visited( G.no_of_vertices + 1, 0 );
  vector<arc> bfstree( G.no_of_vertices + 1, 0 );
  deque<node> order;
  deque<node> Q;
  Q.push_back( s );
  visited[s] = 1;


  while ( !Q.empty() ) {
    const node v = Q.front();
#ifdef VERBOSE
    cout << "v : " << v << endl;
#endif
    for ( auto a : G.rounding_scheme_incident_edges[v] ) {
      node w;
#ifdef VERBOSE
      cout << "a: " << a << endl;
#endif
      if ( a > 0 ) {
        assert( G.tails[a] == v );
        w = G.heads[a];
      } else {
        assert( G.heads[-a] == v );
        w = G.tails[-a];
      }
      if ( !visited[w] ) {
        Q.push_back( w );
        visited[w] = 1;
        bfstree[w] = a;
      }
    }


    Q.pop_front();
#ifdef VERBOSE
    cout << "pushed into order : " << v << endl;
#endif
    order.push_front( v );
  }


#if !defined(NDEBUG) || defined(VERBOSE)
  const node r = order.back();
#ifdef VERBOSE
  cout << "r: " << r << endl;
#endif
#endif
  assert( bfstree[r] == 0 );
  order.pop_back();
  for ( auto v : order ) {

#ifdef VERBOSE
    cout << "v: " << v << endl;
#endif
    const arc a = bfstree[v];
#ifdef VERBOSE
    cout << "demand before: b[" << v << "] = " << b[v] << endl;
#endif

#ifdef VERBOSE
    cout << "flow being pulled: " << b[v] << endl;
#endif
    const RationalType delta = b[v];
    G.pull_flow( -a, delta , b );
#ifdef VERBOSE
    cout << "demand after: b[" << v << "] = " << b[v] << endl;
#endif
    assert( RationalType(1e10)*fabs(b[v]) < RationalType(1) );
    const arc aa = abs( a );
    z[aa] = G.flows[aa];
  }

#ifdef VERBOSE
  cout << "b[" << r << "] = " << b[r] << endl;
#endif
  assert( b[r] == 0 );

#ifdef VERBOSE
  cout << endl << "The z vector" << endl;
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    cout << "z[" << i << "]: " << z[i] << endl;
  }
#endif
  primal_sanity_check( G, z );

  return z;
}

/** \brief Writes the Random Graph used in the dimacs format
 *
 * @param G The Graph
 *
 */
template< typename Network >
void write_the_graph_in_dimacs_format( const Network& N )
{
  ofstream graph_file;
  graph_file.open("graph.txt");
  graph_file << "p" << " " << "min" << " " << N.G.no_of_vertices << " " << N.G.no_of_edges << endl;
  for (unsigned int v = 1; v <= N.G.no_of_vertices; v++) {
    graph_file << "n" << " " << v << " " << -N.nodedata[v].demand << endl;
  }
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
    graph_file << "a" << " " << N.G.tail(a) << " " << N.G.head(a) << " " << 0 << " " << N.arcdata[a].capacity << " " << N.arcdata[a].cost << endl;
  }
}

template<typename Network>
typename Network::IntegerType find_beta(Network& N)
{
  typedef typename Network::IntegerType IntegerType;
  long long gcd = to_double(N.arcdata[1].capacity);

  for (unsigned int a = 1; a <= N.G.no_of_edges; a++)
  {
    long long capacity = to_double(N.arcdata[a].capacity);
    gcd = boost::math::gcd(gcd, capacity);
  }
  for (unsigned int v = 1; v <= N.G.no_of_vertices; v++)
  {
    long long demand = to_double(N.nodedata[v].demand);
    gcd = boost::math::gcd(gcd, demand);
  }

  IntegerType beta = gcd;
  return beta;
}

template<typename Network>
typename Network::IntegerType find_gamma(Network& N)
{
  typedef typename Network::IntegerType IntegerType;
  long long gcd = to_double(N.arcdata[1].cost);
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++)
  {
    long long cost = to_double(N.arcdata[a].cost);
    gcd = boost::math::gcd(gcd, cost);
  }
  IntegerType gamma = gcd;

  return gamma;
}

template<typename Network>
typename Network::IntegerType find_U(Network& N,
                                     typename Network::IntegerType beta)
{
  typedef typename Network::IntegerType IntegerType;

  IntegerType u_inf(0);
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++)
  {
    if (u_inf < N.arcdata[a].capacity)
    {
      u_inf = N.arcdata[a].capacity;
    }
  }

  IntegerType b_1(0);
  for (unsigned int v = 1; v <= N.G.no_of_vertices; v++)
  {
    b_1 += abs(N.nodedata[v].demand);
  }

  long long U = max(to_double(u_inf) / to_double(beta), to_double(b_1) / ((2) * to_double(beta)));

  return U;
}

template<typename Network>
typename Network::IntegerType find_C(Network& N,
                                     typename Network::IntegerType gamma)
{
  typedef typename Network::IntegerType IntegerType;
  IntegerType c_inf(0);
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++)
  {
    if (c_inf < N.arcdata[a].cost)
    {
      c_inf = N.arcdata[a].cost;
    }
  }

  long long C_ = ceil(to_double(c_inf) / to_double(gamma));
  IntegerType C =  C_;

  return C;
}

template <
  typename Network
  >
typename Network::RationalType
find_Gamma(
  const Network& N
) {
  typename Network::RationalType Gamma = 0.0;
  // Gamma is set to the max of c_i, u_i
  for (unsigned int i = 1; i <= N.G.no_of_edges; i++) {
    const auto cap = N.arcdata[i].capacity;
    if ( cap > Gamma) {
      Gamma = cap;
    }
    const auto cost = N.arcdata[i].cost;
    if ( cost > Gamma) {
      Gamma = cost;
    }
  }

  // Compute b_1/2
  typename Network::RationalType sum_abs_demands = 0.0;
  for (unsigned int i = 1; i <= N.G.no_of_vertices; i++) {
    sum_abs_demands += fabs( N.nodedata[i].demand );
  }

  // Gamma = max{c,u,b/2}
  if (2 * Gamma < sum_abs_demands) {
    Gamma = sum_abs_demands / 2;
  }

  return Gamma;
}



template<typename Network>
void scale_b_u_c(Network& N,
                 typename Network::IntegerType beta_scaling_factor,
                 typename Network::IntegerType gamma_scaling_factor,
                 typename Network::RationalType delta
                )
{

  unsigned int n = N.G.no_of_vertices;
  unsigned int m = N.G.no_of_edges;

  // demand scaling
  for (unsigned int v = 1; v <= n; v++) {
    N.nodedata[v].demand *= beta_scaling_factor;
  }
  // capacity scaling
  for (unsigned int a = 1; a <= m; a++) {
    N.arcdata[a].capacity *= beta_scaling_factor;
  }
  // cost scaling
  for (unsigned int a = 1; a <= m; a++) {
    N.arcdata[a].cost *= gamma_scaling_factor;
  }
}


template<typename Network>
typename Network::IntegerType compute_t(Network& N,
                                        typename Network::IntegerType beta,
                                        typename Network::IntegerType gamma,
                                        typename Network::RationalType delta)
{
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;
  typedef typename Network::RationalType RationalType;

  Graph& G0 = N.G;
  unsigned int m = 3 * G0.no_of_edges;

  IntegerType U = find_U(N, beta);
  IntegerType C = find_C(N, gamma);

  IntegerType t_path_following_ =  convert_rounded_ratio_to_integer(RationalType(2) * to_double(beta * gamma * m * U * C) , delta);
  cout << "compute_t" << endl;
  cout << "t_path_following = " << to_double(t_path_following_) << " = " << RationalType(2) * to_double(beta * gamma * m * U * C)
						     << "/" << delta << endl;
  IntegerType t_path_following = t_path_following_;
  t_path_following = t_path_following_;
  return t_path_following;
}


/** \brief Find an initial solution
 *
 * Takes in the original graph G0, feasible integral solution z and constructs a graph G1 and finds initial solution to G1
 *
 * @param Graph orginial graph G0
 * @param Graph new graph G1
 * @param vector initial primal solution x
 * @param vector initial dual solution y
 * @param vector initial slack variables s
 * @param vector integral feasible solution to G0
 * @param template parameter q
 *
 */
template <
  typename Network,
  typename RationalType,
  typename IntegerType
  >
IntegerType find_initial_solution(
  Network& N,
  vector<IntegerType> &z,
  RationalType &q,
  typename Network::IntegerType t_path_following
) {
  typedef typename Network::GraphType Graph;
  const Graph& G0 = N.G;

  const unsigned int nedges_G0  = G0.no_of_edges;
//  const unsigned int n          = G0.no_of_vertices;

#ifdef VERBOSE
  cout << "find_initial_solution()" << endl;
#endif

#ifndef NDEBUG
  cout << "G0 looks like: " << endl;
  unsigned int i, j;
  for (i = 1; i <= G0.no_of_vertices; i++) {
    cout << "Vertex " << i; cout << " :";
    for (j = 0; j < G0.incident_edges[i].size(); j++ ) {
      cout << G0.incident_edges[i][j];
      cout << ",";
    }
    cout << std::endl;
  }

  cout << "Demands" << endl;
  for (unsigned int i = 1; i <= G0.no_of_vertices; i++) {
    cout << i << ": " << to_double(N.nodedata[i].demand) << endl;
  }

  cout << endl << endl;
#endif


#ifdef VERBOSE
  cout << "z and capacities of G0: " << endl;
  for (unsigned int i = 1; i <= nedges_G0; i++) {
    cout << "z[" << i << "]: " << z[i] << ", "
         << "cap[" << i << "]: " << N.arcdata[i].capacity << ", ";
  }
  cout << endl;
#endif


//   for(unsigned int i=1; i<=n; i++){
//     G1.demands[i] = G0.demands[i];
//   }


  // determine t, new m, n, Gamma...

  //#ifdef New_Initialization
  int no_edges_aux = 3 * G0.no_of_edges; // RationalType(2) * nedges_G0 +  n - 1;

  cout << "no_edges_aux = " << no_edges_aux << endl;

  q = no_edges_aux + ceil( sqrt(no_edges_aux) );

  // Gamma = max{c_i,u_i,b/2}
  RationalType Gamma = 0.0;
  RationalType sum_costs = 0.0;

  for (unsigned int i = 1; i <= G0.no_of_edges; i++) {

    const auto& data = N.arcdata[i];
    sum_costs += to_double(data.cost);

    if ( to_double(data.capacity) > Gamma) {
      Gamma = to_double(data.capacity);
    }
    if ( to_double(data.cost) > Gamma ) {
      Gamma = to_double(data.cost);
    }
  }

  RationalType Sum_abs_demands = 0.0;
  for (unsigned int i = 1; i <= G0.no_of_vertices; i++) {
    Sum_abs_demands += fabs( to_double(N.nodedata[i].demand) );
  }
  RationalType Sum_abs_demands_half = Sum_abs_demands / 2;

  if (Gamma < Sum_abs_demands_half) {
    Gamma = Sum_abs_demands_half;
  }
  cout << "Gamma = " << Gamma << endl;

  IntegerType t = t_path_following; // max(sum_c_times_sum_b_2, 4*ceil(sqrt(no_edges_aux) - 1)* Gamma * Gamma);
//
//     cout << "Parameter t: " << t << endl;
//
//     // determine_t(G0,z,q - no_edges_aux );
//     // cout << "sum_costs = " << sum_costs << endl;
//     // cout << "other = " << ceil(sqrt(G0.no_of_edges * 3))* Gamma * Gamma << endl;
//     // cout << "det t = " << determine_t(G0,z,q - 3*nedges_G0 ) << endl;
//   #else

  q = (RationalType(3) * G0.no_of_edges) + ceil( sqrt( RationalType(3) * G0.no_of_edges) );
  //RationalType t =  determine_t(N,z,q - 3*nedges_G0 );

  //#endif

#ifdef VERBOSE
  cout << "t =  " << t << endl;
#endif


  // resize x,s,y
#ifdef New_Initialization
//     x.resize(no_edges_aux  + 1, 0);
//     s.resize(no_edges_aux  + 1, 0);
//     y.resize(n + nedges_G0 + 1, 0);
//
//   #else
  x.resize(3 * nedges_G0 + 1, 0 );
  s.resize(3 * nedges_G0 + 1, 0 );
  y.resize(n + nedges_G0 + 1, 0 );
#endif

  // find new u,c,b
  // node cur_node_index = n + 1;
  // cout << "n = " << n << endl;

//   node cur_node_index_delta_wye = n + 1;
//   node cur_node_index_delta_wye_extra = G1.no_of_vertices + 1;
  for (unsigned int i = 1; i <= nedges_G0; i++) {

    // create new node
    // node vw = cur_node_index;
    // cur_node_index ++;
    node v  = G0.tails[i];
    node w  = G0.heads[i];

    // // for delta-wye construction
    // node vw1 = cur_node_index_delta_wye_extra++;
    // node vw2 = cur_node_index_delta_wye++;
    // node v1 = G0.tails[i];
    // node w1 = G0.heads[i];
    // G_delta_wye.new_edge(w1, vw1);
    // G_delta_wye.new_edge(v1, vw1);


    IntegerType cap  = N.arcdata[i].capacity;
    // RationalType cost = N.arcdata[i].cost;
    assert(cap != 0);

    // Setting demands of G1
    N.nodedata[w].demand  +=  - cap;
    // G1.demands[vw] = cap;

    // // This is the data structure for representing G1 in G0
    // G0.original_auxiliary_transformed_corrospondence[i].demand1 = G1.demands[v];
    // G0.original_auxiliary_transformed_corrospondence[i].demand2 = G1.demands[w];
    // G0.original_auxiliary_transformed_corrospondence[i].demand3 = G1.demands[vw];

    // G0.original_auxiliary_transformed_corrospondence[i].capacity = cap;

    // // insert new arcs
    // arc vvw = G1.new_edge(v, vw);
    // arc wvw = G1.new_edge(w, vw);

    // Setting costs of G1
    // G1.costs[wvw] = 0;
    // G1.costs[vvw] = cost;

    // // For G0-G1 data structure
    // G0.original_auxiliary_transformed_corrospondence[i].cost1 = G1.costs[vvw];


    // // #ifdef VERBOSE
    //   cout << "costs of " << vvw << " = " << G1.costs[vvw] << endl;
    //   cout << "costs of " << wvw << " = " << G1.costs[wvw] << endl;
    // // #endif

    // Setting x of vvw and wvw
    long long xlower_ = to_double(cap) / 2.0;
    N.arcdata[i].xlower = xlower_;
    // x[wvw] = cap/2.0;

    // // G0-G1 data structure
    // G0.original_auxiliary_transformed_corrospondence[i].x1 = x[vvw];
    // G0.original_auxiliary_transformed_corrospondence[i].x2 = x[wvw];
    // G0.original_auxiliary_transformed_corrospondence[i].cost2 = G1.costs[wvw];


    // Setting y
    N.nodedata[v].potential  = 0;
    N.nodedata[w].potential  = 0;
    long long potentialvw_ = -2.0 * to_double(t) / to_double(cap);
    N.arcdata[i].potentialvw = potentialvw_;


    // // G0-G1 data structure
    // G0.original_auxiliary_transformed_corrospondence[i].y1 = y[v];
    // G0.original_auxiliary_transformed_corrospondence[i].y2 = y[w];
    // G0.original_auxiliary_transformed_corrospondence[i].y3 = y[vw];


    // Setting s of vvw and wvw
    N.arcdata[i].slower = N.arcdata[i].cost - N.arcdata[i].potentialvw;
    N.arcdata[i].supper =                   - N.arcdata[i].potentialvw;


    // // G0-G1 data structure
    // G0.original_auxiliary_transformed_corrospondence[i].s1 = s[vvw];
    // G0.original_auxiliary_transformed_corrospondence[i].s2 = s[wvw];
    // // fill arc_map for reconstruction of original instance
    // G1.arc_map[i].orig = i;
    // G1.arc_map[i].vvw  = vvw;
    // G1.arc_map[i].wvw  = wvw;

    // a hat and x and s of a hat
    if (to_double(z[i]) == to_double(cap) / 2.0) {
      N.arcdata[i].direction = 0;
      // // G0-G1 data structure
      // G0.original_auxiliary_transformed_corrospondence[i].x3 = 0;
    } else {
      if (to_double(z[i]) > to_double(cap) / 2.0) {
        N.arcdata[i].direction = 1;
      } else {
        N.arcdata[i].direction = -1;
      }
      // arc vw = G1.new_edge(v,w);
      // G_delta_wye.new_edge(vw2, vw1);
      // Setting x, s, cost of a hat
      if (N.arcdata[i].direction == 1)
      {
        long long infeasibility_ = fabs(to_double(z[i]) - to_double(cap) / 2.0 );
        N.arcdata[i].infeasibility = infeasibility_;
      }
      else
      {
        long long infeasibility_ = fabs(to_double(z[i]) - to_double(cap) / 2.0 );
        N.arcdata[i].infeasibility = infeasibility_;
      }

      long long croof_ = N.arcdata[i].direction == 0 ? 1 : ceil(to_double(t) / fabs( to_double(z[i]) - to_double(cap) / 2.0 ) );
      N.arcdata[i].croof = croof_;
      N.arcdata[i].sroof = N.arcdata[i].croof;
      cout << "x^[" << i << "] = " << to_double(N.arcdata[i].infeasibility) << " c^[" << i << "] = " << to_double(N.arcdata[i].croof) << endl;
      assert( N.arcdata[i].croof > 0 );

      // // G0-G1 data structure
      // G0.original_auxiliary_transformed_corrospondence[i].edge_reversed = false;
      // G0.original_auxiliary_transformed_corrospondence[i].x3 = x[vw];
      // G0.original_auxiliary_transformed_corrospondence[i].cost3 = G1.costs[vw];
      // G0.original_auxiliary_transformed_corrospondence[i].s3 = s[vw];

      // // data structure for getting original solution in the end
      // G1.arc_map[i].ahat = vw;
    }
  }

  // #ifndef NDEBUG
  //   cout << "A look at the Graph later: " << endl;
  //   for(unsigned int i = 1; i <= G1.no_of_vertices; i++){
  //     cout << "Vertex: " << i << " Arcs: " ;
  //     for(auto a: G1.incident_edges[i]){
  //       cout << a << " , " ;
  //     }
  //     cout << endl;
  //   }

  //   cout << " Arc  " << " Head   " << " Tail" << endl;
  //   for(unsigned int i = 1; i <= G1.no_of_edges; i++){

  //     cout << i <<"  " << G1.heads[i] << " " << G1.tails[i] << endl;
  //   }
  // #endif

  // #ifndef NDEBUG
  //   cout << "A look at the Graph_delta_wye later: " << endl;
  //   for(unsigned int i = 1; i <= G_delta_wye.no_of_vertices; i++){
  //     cout << "Vertex: " << i << " Arcs: " ;
  //     for(auto a: G_delta_wye.incident_edges[i]){
  //       cout << a << " , " ;
  //     }
  //     cout << endl;
  //   }
  // #endif

  // #ifdef VERBOSE
  //   cout << "Arc map looks like: " << endl;
  //   for(auto am : G1.arc_map){
  //     cout << "orig: " << am.orig << ", ";
  //     cout << "vvw: "  << am.vvw  << ", ";
  //     cout << "wvw: "  << am.wvw  << ", ";
  //     cout << "ahat: " << am.ahat << endl;
  //   }
  //   cout << endl;

  //   cout << "G1 no of verticies: " << G1.no_of_vertices << endl;

  //   cout<<"-----------------DEMANDS--------------------" << endl;
  //   for(unsigned int i=1; i<=G1.no_of_vertices; i++){
  //     cout << "b["<< i << "]: ";
  //     cout<<G1.demands[i]<<" , " << endl;
  //   }
  //   cout << endl;
  // #endif


  // checks and asserts

  // This is commented out, because THRESHOLD is currently not known here...
  // sanity_check( G1, x, y, s );

  // #ifdef VERBOSE
  //   cout << endl;
  //   cout << "vector X " <<  endl;
  //   for(unsigned i = 1; i <=G1.no_of_edges; i++){

  //     cout << x[i] << " , " ;

  //   }
  //   cout << endl;
  //   cout << "vector S " <<  endl;
  //   for(unsigned i = 1; i <=G1.no_of_edges; i++){
  //     cout << s[i] << " , " ;
  //   }
  // #endif

  // compute duality gap and log of it
  IntegerType duality_gap = 0;
  RationalType log_x_s     = 0.0;

  for (unsigned int i = 1; i <= G0.no_of_edges; i++) {
    IntegerType xisi = N.arcdata[i].xlower * N.arcdata[i].slower;
    duality_gap += xisi;
    log_x_s += log(to_double(xisi));

    xisi =  (N.arcdata[i].capacity - N.arcdata[i].xlower) * N.arcdata[i].supper;
    duality_gap += xisi;
    log_x_s += log(to_double(xisi));

    if (N.arcdata[i].direction != 0) {
      xisi =  N.arcdata[i].infeasibility * N.arcdata[i].sroof;
      duality_gap += xisi;
      log_x_s += log(to_double(xisi));
    }
  }

  // #ifdef VERBOSE
  //   for(unsigned int j=1; j<=G1.no_of_vertices; j++){
  //     cout << "b[" << j << "]=" << G1.demands[j] << ", ";
  //   }
  //   const unsigned int m = G1.no_of_edges;
  //   const unsigned int q1 = m + static_cast<unsigned int>( ceil( sqrt(m) ) );
  //   const RationalType mlogm = m*log(m);
  //   cout<<"Duality Gap: "<<duality_gap << std::endl;
  //   cout<<"Potential: " << q1 << " * " << log(duality_gap) << " - " << log_x_s << " - " << mlogm << " = " << q1*log(duality_gap)-log_x_s - mlogm << endl;

  //   cout << "G1 looks like: " << endl;
  //   cout << "x s of G1 look like: " << endl;

  //   for(unsigned int i=1; i<=G1.no_of_edges;i++){
  //     cout << x[i] << ", ";
  //   }

  //   cout << endl;
  // #endif

  return t;
}

template<typename Network>
typename Network::IntegerType update_duality_gap(Network& N) {

  typedef typename Network::IntegerType IntegerType;

  IntegerType duality_gap = 0;
  for (unsigned int i = 1; i <= N.G.no_of_edges; i++) {
#ifdef VERBOSE
    cout << "duality_gap = " << duality_gap << endl;
#endif
    IntegerType xisi = N.arcdata[i].xlower * N.arcdata[i].slower;
#ifdef VERBOSE
    cout << to_double(N.arcdata[i].xlower) << " * " << to_double(N.arcdata[i].slower) << endl;
#endif
    duality_gap += xisi;
#ifdef VERBOSE
    cout << "duality_gap = " << duality_gap << endl;
#endif
    xisi =  (N.arcdata[i].capacity - N.arcdata[i].xlower) * N.arcdata[i].supper;
#ifdef VERBOSE
    cout << (N.arcdata[i].capacity - N.arcdata[i].xlower) << " * " << N.arcdata[i].supper << endl;
#endif
    duality_gap += xisi;
#ifdef VERBOSE
    cout << "duality_gap = " << duality_gap << endl;
#endif
    xisi =  N.arcdata[i].infeasibility * N.arcdata[i].sroof;
#ifdef VERBOSE
    cout << N.arcdata[i].infeasibility << " * " << N.arcdata[i].sroof << endl;
#endif

    duality_gap += xisi;

  }

  return duality_gap;
}

template<typename Network, typename RationalType, typename IntegerType>
void get_spanning_tree_for_electrical_flow_problem(typename Network::GraphType& G0,
    SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST_original,
    Network& N)
{



  G0.create_low_stretch_tree_wrt_unrounded_resistances_on_original_graph(ST_original, N);

#ifdef VERBOSE
  cout << "root of the tree = " << ST_original.root << endl;
#endif
  assert( ST_original.node_tree_corrospondance[ST_original.root].back().tree_index == 0 );

  unsigned int tree_index = 1;

  ST_original.get_tree_incident_edges();
  ST_original.get_initial_state_of_the_tree_decomposition();
  ST_original.init( ST_original.root, tree_index );
  ST_original.get_non_tree_edge_common_tree_index();

#ifdef VERBOSE
  print_tree_edges(N);
  print_non_tree_edges(N);
  print_directions(N);
#endif

  for (unsigned int i = 0; i < G0.non_tree_edges.size(); i++) {
    G0.resistances_accross_cycles.push_back(ST_original.compute_resistance_accross(G0.non_tree_edges[i], N));

#ifdef VERBOSE
    cout << i << " Resistances accross  =  " << to_double(ST_original.compute_resistance_accross(G0.non_tree_edges[i], N)) << endl;
#endif
  }


}



/** \brief Calculate the duality gap
 *
 * Returns the duality gap for the current primal and dual solutions
 *
 * @param Graph G
 * @param vector the current primal solution x
 * @param vector the current slack variables s
 *
 * @return the duality gap
 *
 */
template<typename IntegerType, typename RationalType>
RationalType update_duality_gap( const Graph<IntegerType, RationalType>&G,
                                 const std::vector<RationalType> &x,
                                 const std::vector<RationalType> &s
                               )
{
  RationalType duality_gap = 0;

  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    const RationalType xs = x[i] * s[i];
#ifdef VERBOSE
    cout << i << ": " << x[i] << " * " << s[i] << " = " << xs << endl;
#endif
    duality_gap += xs;
  }

  return duality_gap;
}
/** \brief Returns the sum of the demands of S_k
 *
 * @param G The Graph
 * @param S_k The Set S_k
 *
 * @return Sum of the demands
 *
 */
template<typename IntegerType, typename RationalType>
RationalType b( const Graph<IntegerType, RationalType> &G,
                const std::vector<node>& S_k
              )
{
  RationalType demand = 0.0;
  for (unsigned int i = 0; i < S_k.size(); i++) {
    demand += G.demands[S_k[i]];
  }
  return demand;
}

/** \brief Returns the incoming arcs for the set of nodes S_k
 *
 * @param G The Graph
 * @param S_k The set of nodes
 *
 * @return A vector of incoming arcs
 *
 */
template<typename IntegerType, typename RationalType>
std::vector<arc> del_in( const Graph<IntegerType, RationalType> &G,
                         const std::vector<node> S_k
                       )
{
  std::vector<arc> incoming_arcs;
  for (unsigned int i = 0; i < S_k.size(); i++) {
    for (unsigned int j = 0; j < G.incident_edges[i].size(); j++) {
      if (G.incident_edges[i][j] < 0) {
        incoming_arcs.push_back(G.incident_edges[i][j]);
      }
    }
  }
  return incoming_arcs;
}

/** \brief Returns the outgoing arcs for the set of nodes S_k
 *
 * @param G The Graph
 * @param S_k The set of nodes
 *
 * @return A vector of outgoing arcs
 *
 */
template<typename IntegerType, typename RationalType>
std::vector<arc>
del_out( const Graph<IntegerType, RationalType>&G,
         const std::vector<node> S_k
       )
{
  std::vector<arc> outgoing_arcs;
  for (unsigned int i = 0; i < S_k.size(); i++) {
    for (unsigned int j = 0; j < G.incident_edges[i].size(); j++) {
      if (G.incident_edges[i][j] > 0) {
        outgoing_arcs.push_back(G.incident_edges[i][j]);
      }
    }
  }
  return outgoing_arcs;
}

template<typename IntegerType, typename RationalType>
RationalType
minimal( const Graph<IntegerType, RationalType>& G,
         const std::vector<arc>& arcs,
         const std::vector<RationalType>& y,
         const std::vector<RationalType>& y_0,
         arc& a_k
       )
{
  RationalType min = std::numeric_limits<RationalType>::max();;
  a_k = 0;
  for (auto a : arcs)
  {
    int w = G.head(a);
    int v = G.tail(a);
    if (min > G.costs[a] + y[v] - y_0[w]) {
      min = G.costs[a] + y[v] - y_0[w];
      a_k = a;
    }
  }
  return min;
}

template< typename IntegerType, typename RationalType >
struct node_data {
  arc a;
  node v;
  RationalType data;
  RationalType x;
};

template< typename IntegerType, typename RationalType >
class compare_data {
public:
  bool operator()(const node_data<IntegerType, RationalType>& d1, const node_data<IntegerType, RationalType>& d2) const
  {
    if (d1.data > d2.data)
    {
      return true;
    }
    else
    {
      if (d1.data == d2.data) {

        if (d1.x > d2.x) return true;
      }

      return false;
    }
  }

};

/** \brief Runs the Rounding Scheme Algorithm
 *
 * @param G The Graph
 * @param s The Dual Slacks
 * @param x The Primal Solution
 * @param y_0 Integral dual solution
 * @param s_0 Integral Dual Slacks
 */
template <
  typename Network,
  typename RationalType
  >
void
rounding_scheme(
  Network& N,
  const node& s,
  const std::vector<RationalType> &x,
  const std::vector<RationalType> &y_0,
  const std::vector<RationalType> &s_0,
  const RationalType THRESHOLD_X,
  vector<int>& visited,
  vector<bool>& inside_set_S,
  vector<RationalType>& y_tilde
) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;
//  typedef typename Network::RationalType RationalType;
  const Graph& G = N.G;

  priority_queue<node_data<IntegerType, RationalType>, vector<node_data<IntegerType, RationalType> >, compare_data<IntegerType, RationalType> > incoming_arcs_pq;
  priority_queue<node_data<IntegerType, RationalType>, vector<node_data<IntegerType, RationalType> >, compare_data<IntegerType, RationalType> > outgoing_arcs_pq;

  std::vector<RationalType> in(G.no_of_vertices + 1, numeric_limits<RationalType>::max());
  std::vector<RationalType> out(G.no_of_vertices + 1, numeric_limits<RationalType>::max());
  //std::vector<int> visited(G.no_of_vertices+1 , 0);
#ifdef VERBOSE
  cout << "rounding_scheme begins..." << endl;
  cout << "number of verticies: " << G.no_of_vertices << endl;
#endif
  std::vector<RationalType> y(G.no_of_vertices + 1, 0);

//   for(unsigned int i = 1; i <= G.no_of_vertices; i++){
//
//     y[i] = y_0[i] - y_0[s];
//     cout << y[i] << " = " << y_0[i] << " - " << y_0[s] << endl;
//   }

  std::vector<node> S;
  y[s] = 0;
  RationalType delta_k = -y_0[s];
  S.push_back(s);
  inside_set_S[s] = true;
  visited[s] = 1;

  // push arcs in queues
  for (auto a : G.incident_edges[s]) {

    cout << "THRESHOLD_X = " << THRESHOLD_X << endl;
    if ( x[abs(a)] < THRESHOLD_X ) continue;

    if (a < 0) {
      node w = G.head(a);
      RationalType data = G.costs[-a] + y_0[w] - y[s];

      if (data < in[w]) {
        in[w] = data;
        node_data<IntegerType, RationalType> w_v_data = {a, w, data, x[abs(a)]};
        incoming_arcs_pq.push(w_v_data);
      }
    }
    if (a > 0) {
      node w = G.head(a);
      RationalType data = G.costs[a] + y[s] - y_0[w];

      if (data < out[w]) {
        out[w] = data;
        node_data<IntegerType, RationalType> v_w_data = {a, w, data, x[abs(a)]};
        outgoing_arcs_pq.push(v_w_data);
      }
    }
  }


  while (!(outgoing_arcs_pq.empty() && incoming_arcs_pq.empty()))
    //for(unsigned int k =1; k< G.no_of_vertices; k++)
  {

    arc a_k = 0;
    node v_k = 0;
    node w_k = 0;

    RationalType b_S_k = b(G, S);
    if (b_S_k < 0 || incoming_arcs_pq.empty()) {
      assert(!outgoing_arcs_pq.empty());
      node_data<IntegerType, RationalType> minimal_data = outgoing_arcs_pq.top();
      delta_k = minimal_data.data;

      if (delta_k != out[minimal_data.v] || visited[minimal_data.v] == 1) {
        outgoing_arcs_pq.pop();
        continue;
      }

      w_k = minimal_data.v;
      visited[w_k] = 1;
      a_k = minimal_data.a;
      if (a_k > 0 ) {
        G.rounding_scheme_incident_edges[G.tails[a_k]].push_back(a_k);
        G.rounding_scheme_incident_edges[G.heads[a_k]].push_back(-a_k);
      }

      if (a_k < 0 ) {
        cout << G.tails[-a_k] << " -> " << -a_k << endl;
        G.rounding_scheme_incident_edges[G.tails[-a_k]].push_back(-a_k);
        cout << G.heads[-a_k] << " -> " << a_k << endl;
        G.rounding_scheme_incident_edges[G.heads[-a_k]].push_back(a_k);
      }

#ifdef VERBOSE
      cout << "Rounding: " << k << " - > " << abs(a_k) << endl;
#endif
      v_k = G.tail(a_k);
      assert( G.head(a_k) == w_k );
      cout << "delta_k: " << delta_k << endl;
      cout << "v_k: " << v_k << endl;
      cout << "w_k: " << w_k << endl;
      cout << "a_k: " << a_k << endl;
      //#endif
    }

    else {
      assert(!incoming_arcs_pq.empty());
      node_data<IntegerType, RationalType> minimal_data = incoming_arcs_pq.top();
      delta_k = -minimal_data.data;
      if (-delta_k != in[minimal_data.v] || visited[minimal_data.v] == 1) {
        incoming_arcs_pq.pop();
        continue;
      }

      w_k = minimal_data.v;
      visited[w_k] = 1;
      a_k = minimal_data.a;
      if (a_k > 0 ) {
        G.rounding_scheme_incident_edges[G.tails[a_k]].push_back(a_k);
        G.rounding_scheme_incident_edges[G.heads[a_k]].push_back(-a_k);
      }

      if (a_k < 0 ) {
        cout << G.tails[-a_k] << " - > " << -a_k << endl;
        G.rounding_scheme_incident_edges[G.tails[-a_k]].push_back(-a_k);
        cout << G.heads[-a_k] << " -> " << a_k << endl;
        G.rounding_scheme_incident_edges[G.heads[-a_k]].push_back(a_k);
      }
#ifdef VERBOSE
      cout << "Rounding: " << k << " - > " << abs(a_k) << endl;
#endif
      v_k = G.tail(a_k);
      assert(G.head(a_k) == w_k);
      cout << "delta_k: " << delta_k << endl;
      cout << "v_k: " << v_k << endl;
      cout << "w_k: " << w_k << endl;
      cout << "a_k: " << a_k << endl;
    }

    assert( w_k != 0 );
    y[w_k] = y_0[w_k] + delta_k;

    cout << "updation of y: , w_k = " << w_k << " :  ";
    cout << y[w_k] << " = " << y_0[w_k] << " + " << delta_k << endl;

    RationalType b_T_y0 = 0.0;
    RationalType b_T_y = 0.0;

#ifdef VERBOSE
    cout << w_k << ": " << y[w_k] << " = " << y_0[w_k] << " + " << delta_k << endl;
    cout << endl << "assertion val: " << G.costs[abs(a_k)] + y[v_k] - y[w_k] << " = " << G.costs[abs(a_k)] << " + " << y[v_k] << " - " << y[w_k] << endl;
#endif

    assert( RationalType(1e6) * (G.costs[abs(a_k)] + y[v_k] - y[w_k]) > - RationalType(1));
    S.push_back(w_k);
    inside_set_S[w_k] = true;

#ifdef VERBOSE
    cout << "w_k = " << w_k << endl;
#endif

    for (unsigned int i = 1; i <= G.no_of_vertices; i++) {

#ifdef VERBOSE
      cout << "node in S = " << i << endl;
#endif

      b_T_y0 += G.demands[i] * y_0[i];
      b_T_y += G.demands[i] * y[i];

#ifdef VERBOSE
      cout << G.demands[i] << " * " << y[i] << "      " << G.demands[i] << " * " << y_0[i] << endl;
#endif
    }

    if (b_T_y < b_T_y0) {

#ifdef VERBOSE
      cout << b_T_y << " " << b_T_y0 << endl;
      cout << "difference = " << b_T_y - b_T_y0 << endl;
#endif

      RationalType sum_dem = 0.0;
      for (auto i : S) {
        sum_dem += G.demands[i];
      }
#ifdef VERBOSE
      cout << "sum_dem = " << sum_dem << endl;
#endif
    }

    for (auto a : G.incident_edges[w_k]) {
      if (x[abs(a)] < THRESHOLD_X) {
        continue;
      }

      if (a < 0) {
        node w = G.head(a);
        RationalType data = G.costs[-a] + y_0[w] - y[w_k];
        if (data < in[w] && visited[w] == 0) {
          node_data<IntegerType, RationalType> w_v_data = {a, w, data, x[abs(a)]};
          incoming_arcs_pq.push(w_v_data);
          in[w] = data;
        }
#ifdef VERBOSE
        cout << "arc: " << a << " tail: " << w_k << " " << G.head(a) << " " << w << " " << in[w] << " " << data << endl;
#endif
      }

      if (a > 0) {
        node w = G.head(a);
        RationalType data = G.costs[a] + y[w_k] - y_0[w];
        if (data < out[w] && visited[w] == 0) {
          node_data<IntegerType, RationalType> v_w_data = {a, w, data, x[abs(a)]};
          outgoing_arcs_pq.push(v_w_data);
          out[w] = data;
#ifdef VERBOSE
          cout << "updated: " << out[w] << endl;
#endif
        }
#ifdef VERBOSE
        cout << "arc: " << a << " tail: " << w_k  << " " << G.tail(a) << " " << w << " " << out[w] << " " << data << endl;
#endif
      }

    }
#ifdef VERBOSE
    cout << "s_k+1 later " << S.size() << endl;
#endif
  }

  RationalType b_T_y0 = 0.0;
  RationalType b_T_y = 0.0;
  for (auto i : S)
    //for(unsigned int i=1; i<=G.no_of_vertices; i++)
  {
    //if(inside_set_S[i])
    {
      b_T_y0 += G.demands[i] * y_0[i];
      b_T_y += G.demands[i] * y[i];
      cout << i << ": " << b_T_y << " = " << G.demands[i] << " * " << y[i] << endl;
    }

  }


  cout << "b_T_y0 = " << b_T_y0 << endl;
  cout << "b_T_y = " << b_T_y << endl;
  //if(check_b_T_y_assert)
  {
    assert(b_T_y >= b_T_y0 );
  }

}


/** \brief Compute the Dual Slacks
 *
 * Computes and returns the dual slack variables
 *
 * @param Graph G
 * @param vector y Dual Variables
 * @return Dual Slack Variables
 */
template <
  typename Network,
  typename RationalType
  >
vector<RationalType> compute_dual_slacks(
  const Network& N,
  const vector<RationalType>& y
) {
  const auto& G = N.G;
  vector<RationalType> s( G.no_of_edges + 1 );

  for ( unsigned int i = 1; i <= G.no_of_edges; ++i ) {
    const auto& cost = N.arcdata[i].cost;
    s[i] = cost + y[G.tails[i]] - y[G.heads[i]];
    if (s[i] < 1e-6) {
      s[i] = 0;
    }

#ifdef VERBOSE
    cout << "s[" << i << "] = " << s[i] << endl;
    cout << costs << " + " << y[G.tails[i]] << " - " << y[G.heads[i]] << endl;
#endif

    assert( RationalType(1e6) *s[i] > -RationalType(1) );
  }
  return s;
}

/** \brief Compute Primal Objective Value
 *
 * Computes and returns the Primal Objective Value
 *
 * @param Graph G
 * @param vector x Primal Solution
 * @return Primal Objective Value
 */
template <
  typename IntegerType,
  typename RationalType
  >
RationalType primal_obj_value(
  const Graph<IntegerType,
  RationalType>& G,
  const vector<RationalType>& x
) {
  RationalType obj = 0;
  for ( unsigned i = 1; i <= G.no_of_edges; ++i ) {
    const RationalType caxa = G.costs[i] * x[i];
    obj += caxa;

#ifdef VERBOSE
    cout << G.costs[i] << " * " << x[i] << " = " << caxa << endl;
#endif
  }
  return obj;
}


/** \brief Compute Dual Objective Value
 *
 * Computes and returns the dual objective value
 *
 * @param Graph G
 * @param vector y Dual Solution
 * @return Dual objective value
 */
template< typename IntegerType, typename RationalType >
RationalType
dual_obj_value( const Graph<IntegerType, RationalType>& G,
                const vector<RationalType>& y
              )
{
  RationalType obj = 0;
  for ( unsigned i = 1; i <= G.no_of_vertices; ++i ) {
    const RationalType bvyv = G.demands[i] * y[i];
    obj += bvyv;
#ifdef VERBOSE
    cout << G.demands[i] << " * " << y[i] << " = " << bvyv << endl;
#endif
  }
  return obj;
}

template<typename Network, typename RationalType>
RationalType determine_t(Network& N, const vector<RationalType>& z, unsigned int p)
{
  RationalType sum = 0;
  RationalType cTu = 0;
  RationalType Tu = 0;
  RationalType max1 = 1;
  for ( unsigned int i = 1; i <= N.G.no_of_edges; ++i ) {
    Tu += N.arcdata[i].cost;
    cTu += N.arcdata[i].cost * N.arcdata[i].capacity;
    const RationalType d = fabs(z[i] - (N.arcdata[i].capacity) / 2.0 );
    if ( d > max1 ) {
      max1 = d;
#ifdef VERBOSE
      cout << "fabs(" << z[i] << " - " << (G.capacities[i]) << " / " << 2.0  << ") " << endl;
#endif
    }
    sum += d;
  }
#ifdef VERBOSE
  cout << "max1 = " << max1 << endl;
  cout << endl << "cTu = " << cTu << endl;
#endif
  return max( Tu * max1, (cTu + sum) / p );
}


/** \brief Find the parameter 't'
 *
 * Calculates and returns the parameter 't'
 *
 * @param Graph G
 * @param vector z
 * @param template p
 * @return parameter 't'
 */
template< typename IntegerType, typename RationalType >
RationalType
determine_t( const Graph<IntegerType, RationalType>& G,
             const vector<RationalType>& z,
             unsigned int p
           )
{
  RationalType sum = 0;
  RationalType cTu = 0;
  RationalType Tu = 0;
  RationalType max1 = 1;
  for ( unsigned i = 1; i <= G.no_of_edges; ++i ) {
    Tu += G.costs[i];
    cTu += G.costs[i] * G.capacities[i];
    const RationalType d = fabs(z[i] - (G.capacities[i]) / 2.0 );
    if ( d > max1 ) {
      max1 = d;
#ifdef VERBOSE
      cout << "fabs(" << z[i] << " - " << (G.capacities[i]) << " / " << 2.0  << ") " << endl;
#endif
    }
    sum += d;
  }
#ifdef VERBOSE
  cout << "max1 = " << max1 << endl;
  cout << endl << "cTu = " << cTu << endl;
#endif
  return max( Tu * max1, (cTu + sum) / p );
}

int read_n(string filepath) {
  int n = 0;
  int nodeid;
  ptree tree;
  read_xml(filepath, tree);
  const ptree & graphml = tree.get_child("graphml", empty_ptree());
  const ptree & graph   = graphml.get_child("graph", empty_ptree());

  BOOST_FOREACH(const ptree::value_type & nore, graph) {
    const ptree & nore_attrs = nore.second.get_child("<xmlattr>", empty_ptree());
    BOOST_FOREACH(const ptree::value_type & nore_attr, nore_attrs) {
      if (strncmp(nore_attr.first.data(), "id", 2) == 0) {
        nodeid = stoi(nore_attr.second.data());
        n = max(n, nodeid);
      }
    }
  }
  return n;
}

int read_m(string filepath) {
  int m = 0;
  ptree tree;
  read_xml(filepath, tree);
  const ptree & graphml = tree.get_child("graphml", empty_ptree());
  const ptree & graph   = graphml.get_child("graph", empty_ptree());

  BOOST_FOREACH(const ptree::value_type & nore, graph) {
    const ptree & nore_attrs = nore.second.get_child("<xmlattr>", empty_ptree());
    BOOST_FOREACH(const ptree::value_type & nore_attr, nore_attrs) {
      if (strncmp(nore_attr.first.data(), "source", 6) == 0) {
        m++;
      }
    }
  }
  cout << "file contains " << m << " edges " << endl;
  return m;
}


// for ssp variant
template <
  typename Network
  >
void set_initial_x_s_y( Network &N ) {

  typedef typename Network::RationalType RationalType;

  for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {

    N.arcdata[a].xlower = RationalType(0);

    // assuming cost is nonzero
    //N.arcdata[a].slower = N.arcdata[a].cost;
  }

  for (unsigned int v = 1; v <= N.G.no_of_vertices; v++) {
    N.nodedata[v].potential = RationalType(0);
  }

}

template <
  typename Network,
  typename RationalType
  >
void update_x_by_s_cuts(Network &N, arc a_hat, RationalType deficit_s, bool forward) {
  a_hat = std::abs(a_hat);
  auto x_old = N.arcdata[a_hat].xlower;

  if (deficit_s < 0)
    // a_hat is an outgoing edge
    N.arcdata[a_hat].xlower = forward ? std::min(N.arcdata[a_hat].capacity, x_old - deficit_s) : std::max(static_cast<decltype(x_old)>(0), x_old + deficit_s);

  if (deficit_s > 0)
    // a_hat is an incoming edge
    N.arcdata[a_hat].xlower = forward ? std::min(N.arcdata[a_hat].capacity, x_old + deficit_s) : std::max(static_cast<decltype(x_old)>(0), x_old - deficit_s);
}


template<typename Network>
void update_x_by_tree(Network &N, std::vector<int> depth_count) {
  auto depth_sums = std::vector<int>(N.G.no_of_vertices, 0);
  depth_sums[0] = 0;
  for (auto i = 1u; i < depth_sums.size(); ++i)
    depth_sums[i] = depth_sums[i - 1] + depth_count[i - 1];
  auto tree = std::vector<arc>(N.G.no_of_vertices, 0);
  auto num_nodes = 0;
  for (auto v = 1u; v <= N.G.no_of_vertices; ++v) {
    if (N.nodedata[v].depth == -1)
      continue;
    auto depth_v = N.nodedata[v].depth;
    tree[depth_sums[depth_v] + --depth_count[depth_v]] = v;
    ++num_nodes;
  }
  tree.resize(num_nodes);

#ifndef NDEBUG
  for (auto i = 1u; i < tree.size(); ++i)
    assert(N.nodedata[tree[i - 1]].depth <= N.nodedata[tree[i]].depth && "tree should be ordered by depth");
#endif

#ifdef USE_IMBALANCE
  for (auto &nodedata : N.nodedata)
    nodedata.imbalance = nodedata.deficit;
#endif

  for (auto v : boost::adaptors::reverse(tree)) {
    auto a = abs(N.nodedata[v].parent);
    auto i = N.G.tails[a];
    auto j = N.G.heads[a];

    auto x_old = N.arcdata[a].xlower;

#ifdef USE_IMBALANCE
    if (N.nodedata[i].depth > N.nodedata[j].depth) {
      N.nodedata[j].imbalance += N.nodedata[i].imbalance;
      N.arcdata[a].xlower = min(max(static_cast<decltype(x_old)>(0), x_old - N.nodedata[i].imbalance), N.arcdata[a].capacity);
    } else {
      N.nodedata[i].imbalance += N.nodedata[j].imbalance;
      N.arcdata[a].xlower = min(max(static_cast<decltype(x_old)>(0), x_old + N.nodedata[j].imbalance), N.arcdata[a].capacity);
    }
#else
    if (N.nodedata[i].depth > N.nodedata[j].depth) {
      N.arcdata[a].xlower = min(max(static_cast<decltype(x_old)>(0), x_old - N.nodedata[i].deficit), N.arcdata[a].capacity);
    } else {
      N.arcdata[a].xlower = min(max(static_cast<decltype(x_old)>(0), x_old + N.nodedata[j].deficit), N.arcdata[a].capacity);
    }
#endif

    auto flow_incr = N.arcdata[a].xlower - x_old;
    N.nodedata[i].deficit += flow_incr;
    N.nodedata[j].deficit -= flow_incr;
  }
}



template <
  typename Network,
  typename QueueType
  >
void update_s_in_s_out(const Network &N, node v, QueueType &s_in, QueueType &s_out) {
  for (const auto edge : N.G.incident_edges[v]) {
    const auto w = edge > 0 ? N.G.heads[edge] : N.G.tails[-edge];
    if (N.nodedata[w].visited)
      continue;
    const auto &arcdata = N.arcdata[abs(edge)];
    if (arcdata.xlower < arcdata.capacity) {
      if (edge > 0) {
        s_out.push({arcdata.cost + N.nodedata[v].potential - N.nodedata[w].potential, edge});
      } else {
        s_in.push({arcdata.cost + N.nodedata[w].potential - N.nodedata[v].potential, edge});
      }
    }
    if (arcdata.xlower > 0) {
      if (edge > 0) {
        s_in.push({ -arcdata.cost - N.nodedata[v].potential + N.nodedata[w].potential, edge});
      } else {
        s_out.push({ -arcdata.cost - N.nodedata[w].potential + N.nodedata[v].potential, edge});
      }
    }
  }
}



template<typename Network>
void successive_shortest_path(Network &N) {

  auto a_n = std::vector<node>();
  a_n.reserve(N.G.no_of_vertices);
  for (auto i = 1u; i <= N.G.no_of_vertices; ++i) {
    auto &nodedata = N.nodedata[i];
    nodedata.deficit = nodedata.demand;
    if (nodedata.deficit < 0)
      a_n.push_back(i);
  }

  using PotentialType = decltype(N.nodedata.front().potential);
  using DeficitType = decltype(N.nodedata.front().deficit);

  auto step_count = 0u;
  while (!a_n.empty()) {
    auto deficit_1 = std::accumulate(N.nodedata.begin() + 1, N.nodedata.end(), static_cast<DeficitType>(0), [](DeficitType deficit_1, const typename decltype(N.nodedata)::value_type & nodedata) {
      return deficit_1 + std::abs(nodedata.deficit);
    });
    std::cout << "starting step " << ++step_count << ", deficit_1 = " << deficit_1 << std::endl;

    auto v_start = *std::min_element(std::begin(a_n), std::end(a_n), [&N](node lhs, node rhs) {
      return N.nodedata[lhs].deficit < N.nodedata[rhs].deficit || (N.nodedata[lhs].deficit == N.nodedata[rhs].deficit && lhs < rhs);
    });
    //std::cout << "  selected node: " << v_start << std::endl;

    auto delta = -N.nodedata[v_start].potential;
    auto deficit_s = N.nodedata[v_start].deficit;

    using PriorityQueueElementType = std::pair<PotentialType, arc>;
    //auto compare_first = [](const PriorityQueueElementType &lhs, const PriorityQueueElementType &rhs) {
    //  return lhs.first > rhs.first;
    //};
    auto compare_first = [&N](const PriorityQueueElementType & lhs, const PriorityQueueElementType & rhs) {
      if (lhs.first != rhs.first)
        return lhs.first > rhs.first;
      auto lhs_tail = N.G.tails[std::abs(lhs.second)];
      auto lhs_head = N.G.heads[std::abs(lhs.second)];
      auto rhs_tail = N.G.tails[std::abs(rhs.second)];
      auto rhs_head = N.G.heads[std::abs(rhs.second)];
      if (lhs_tail != rhs_tail)
        return lhs_tail > rhs_tail;
      return lhs_head > rhs_head;
    };
    auto s_in = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);
    auto s_out = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);

    for (auto i = 1u; i <= N.G.no_of_vertices; ++i) {
      N.nodedata[i].visited = false;
      N.nodedata[i].depth = -1;
      N.nodedata[i].parent = 0;
    }

    N.nodedata[v_start].visited = true;
    N.nodedata[v_start].depth = 0;

    auto num_visited = 1u;

    auto depth_count = std::vector<int>(N.G.no_of_vertices, 0);
    ++depth_count[0];

    update_s_in_s_out(N, v_start, s_in, s_out);

    while (num_visited < N.G.no_of_vertices && !(s_in.empty() && s_out.empty())) {
      PotentialType s_a_hat;
      arc a_hat;
      node v, w;
#ifdef UPDATE_S_CUTS
      bool out_it_goes;
#endif

      if (deficit_s < 0 || s_in.empty()) {
        // flow is supposed to go out of S
        assert(!s_out.empty());
        std::tie(s_a_hat, a_hat) = s_out.top();
        //std::cout << "    from s_out: edge " << a_hat << " (from " << N.G.tails[std::abs(a_hat)] << " to " << N.G.heads[std::abs(a_hat)] << ") with value " << s_a_hat << std::endl;
        s_out.pop();
        delta = s_a_hat;
#ifdef UPDATE_S_CUTS
        out_it_goes = true;
#endif
      } else {
        // flow is supposed to go into S
        assert(!s_in.empty());
        std::tie(s_a_hat, a_hat) = s_in.top();
        //std::cout << "    from s_in:  edge " << a_hat << " (from " << N.G.tails[std::abs(a_hat)] << " to " << N.G.heads[std::abs(a_hat)] << ") with value " << s_a_hat << std::endl;
        s_in.pop();
        delta = -s_a_hat;
#ifdef UPDATE_S_CUTS
        out_it_goes = false;
#endif
      }

      if (a_hat > 0) {
        v = N.G.tails[a_hat];
        w = N.G.heads[a_hat];
      } else {
        w = N.G.tails[-a_hat];
        v = N.G.heads[-a_hat];
      }

      if (N.nodedata[w].visited)
        // the edge is not a cut edge (anymore)
        continue;

      N.nodedata[w].depth = N.nodedata[v].depth + 1;
      ++depth_count[N.nodedata[w].depth];
      N.nodedata[w].parent = a_hat;

      N.nodedata[w].potential += delta;


#ifdef UPDATE_S_CUTS
      // update x[a_hat]
      auto x_old = N.arcdata[std::abs(a_hat)].xlower;

      update_x_by_s_cuts(N, a_hat, deficit_s, forward);

      auto flow_incr = N.arcdata[std::abs(a_hat)].xlower - x_old;
#endif


      update_s_in_s_out(N, w, s_in, s_out);

      // put w in S
      deficit_s += N.nodedata[w].deficit;
      N.nodedata[w].visited = true;
      ++num_visited;

#ifdef UPDATE_S_CUTS
      // update deficit vector
      if (out_it_goes == forward) {
        deficit[v] += flow_incr;
        deficit[w] -= flow_incr;
      } else {
        deficit[v] -= flow_incr;
        deficit[w] += flow_incr;
      }
#endif
    }

#ifndef UPDATE_S_CUTS
    update_x_by_tree(N, depth_count);
#endif

    // update A_n
    a_n.clear();
    for (auto v = 1u; v <= N.G.no_of_vertices; ++v)
      if (N.nodedata[v].deficit < 0)
        a_n.push_back(v);
  }

  //std::cout << "  primal function value = " << N.calculate_primal_objective_value() << std::endl;

  //for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
  //  cout << "  x[" << a << "] = " << N.arcdata[a].xlower << endl;
  //}
}




template <
  typename Network,
  typename QueueType
  >
void update_s_in_s_out_node(const Network &N, node v, QueueType &s_in_n, QueueType &s_out_n, QueueType &s_in_a, QueueType &s_out_a) {
  const auto &nodedata_v = N.nodedata[v];
  for (const auto edge : N.G.incident_edges[v]) {
    const auto &arcdata = N.arcdata[std::abs(edge)];
    if (arcdata.transformed) {
      if (arcdata.visited)
        continue;
      if (edge > 0) {
        s_out_n.push({arcdata.cost + nodedata_v.potential - arcdata.potential, {edge, false}});
        if (arcdata.xlower > 0)
          s_in_n.push({ -arcdata.cost - nodedata_v.potential + arcdata.potential, { -edge, false}});
      } else {
        s_out_n.push({nodedata_v.potential - arcdata.potential, { -edge, true}});
        if (arcdata.xupper > 0)
          s_in_n.push({ -nodedata_v.potential + arcdata.potential, {edge, true}});
      }
    } else {
      const auto w = edge > 0 ? N.G.heads[edge] : N.G.tails[-edge];
      const auto &nodedata_w = N.nodedata[w];
      if (nodedata_w.visited)
        continue;
      if (edge > 0) {
        s_out_a.push({arcdata.cost + nodedata_v.potential - nodedata_w.potential, {edge, true}});
      } else {
        s_in_a.push({arcdata.cost + nodedata_w.potential - nodedata_v.potential, {edge, false}});
      }
      if (arcdata.xlower > 0) {
        if (edge > 0) {
          s_in_a.push({ -arcdata.cost - nodedata_v.potential + nodedata_w.potential, {edge, true}});
        } else {
          s_out_a.push({ -arcdata.cost - nodedata_w.potential + nodedata_v.potential, {edge, false}});
        }
      }
    }
  }
}


template <
  typename Network,
  typename QueueType
  >
void update_s_in_s_out_arc(const Network &N, arc a, QueueType &s_in_a, QueueType &s_out_a) {
  a = std::abs(a);
  auto i = N.G.tails[a];
  auto j = N.G.heads[a];
  const auto &nodedata_i = N.nodedata[i];
  const auto &nodedata_j = N.nodedata[j];
  const auto &arcdata = N.arcdata[a];
  assert(arcdata.transformed);
  if (!nodedata_i.visited) {
    s_in_a.push({arcdata.cost + nodedata_i.potential - arcdata.potential, {a, false}});
    if (arcdata.xlower > 0)
      s_out_a.push({ -arcdata.cost - nodedata_i.potential + arcdata.potential, { -a, false}});
  }
  if (!nodedata_j.visited) {
    s_in_a.push({nodedata_j.potential - arcdata.potential, {a, true}});
    if (arcdata.xupper > 0)
      s_out_a.push({ -nodedata_j.potential + arcdata.potential, { -a, true}});
  }
}


template <
  typename Network,
  typename QueueType,
  typename TreeType
  >
auto explore_node(Network &N, typename QueueType::value_type top_element, QueueType &s_in_n, QueueType &s_out_n, QueueType &s_in_a, QueueType &s_out_a, TreeType &tree, bool out)->std::pair<bool, decltype(N.nodedata.front().deficit)> {
  auto delta = out ? top_element.first : -top_element.first;
  auto a_hat = top_element.second.first;
  auto upper = top_element.second.second;

  const auto &arcdata = N.arcdata[std::abs(a_hat)];

  auto v = upper ? N.G.heads[std::abs(a_hat)] : N.G.tails[std::abs(a_hat)];
  assert((!N.arcdata[std::abs(a_hat)].transformed || N.arcdata[std::abs(a_hat)].visited) && "v is a node that was added after exploring the edge a_hat");
  assert((N.arcdata[std::abs(a_hat)].transformed || N.nodedata[upper ? N.G.tails[std::abs(a_hat)] : N.G.heads[std::abs(a_hat)]].visited) && "v is a node that was added after exploring the other node attached to the edge a_hat");

  auto &nodedata = N.nodedata[v];

  if (nodedata.visited)
    return {false, 0};

  auto depth = N.arcdata[std::abs(a_hat)].transformed ?
               N.arcdata[std::abs(a_hat)].depth + 1 :
               N.nodedata[upper ? N.G.tails[std::abs(a_hat)] : N.G.heads[std::abs(a_hat)]].depth + 1;
  nodedata.depth = depth;
  nodedata.visited = true;
  if (tree.size() < static_cast<decltype(tree.size())>(depth))
    tree.resize(tree.size() + 1);
  // TODO: verify that the sign of a_hat is correct
  tree[depth - 1].push_back({v, a_hat});
  nodedata.parent = a_hat;

  nodedata.potential += delta;

#ifdef UPDATE_S_CUTS
  // TODO: ...
  // update x[a_hat]
  //auto x_old = N.arcdata[std::abs(a_hat)].xlower;

  //update_x_by_s_cuts(N, a_hat, deficit_s, forward);

  //auto flow_incr = N.arcdata[std::abs(a_hat)].xlower - x_old;
#endif

  update_s_in_s_out_node(N, v, s_in_n, s_out_n, s_in_a, s_out_a);

  return {true, nodedata.deficit};
}


template <
  typename Network,
  typename QueueType,
  typename TreeType
  >
auto explore_arc(Network &N, typename QueueType::value_type top_element, QueueType &s_in_a, QueueType &s_out_a, TreeType &tree, bool out) -> std::pair<bool, decltype(N.arcdata.front().deficit)> {
  auto delta = out ? top_element.first : -top_element.first;
  auto a_hat = top_element.second.first;
  auto upper = top_element.second.second;

  auto parent = upper ? N.G.heads[std::abs(a_hat)] : N.G.tails[std::abs(a_hat)];
  assert(N.nodedata[parent].visited && "a_hat is an arc that was added after exploring the parent node");

  auto &arcdata = N.arcdata[std::abs(a_hat)];

  assert(arcdata.transformed);

  if (arcdata.visited)
    return {false, 0};

  auto depth = N.nodedata[parent].depth + 1;
  arcdata.depth = depth;
  arcdata.visited = true;
  if (tree.size() < static_cast<decltype(tree.size())>(depth))
    tree.resize(tree.size() + 1);
  // TODO: verify that the sign of a_hat is correct
  tree[depth - 1].push_back({parent, a_hat});
  arcdata.parent = parent;

  arcdata.potential += delta;

#ifdef UPDATE_S_CUTS
  // TODO: ...
  // update x[a_hat]
  //auto x_old = N.arcdata[std::abs(a_hat)].xlower;

  //update_x_by_s_cuts(N, a_hat, deficit_s, forward);

  //auto flow_incr = N.arcdata[std::abs(a_hat)].xlower - x_old;
#endif

  update_s_in_s_out_arc(N, a_hat, s_in_a, s_out_a);

  return {true, arcdata.deficit};
}


template <
  typename Network,
  typename TreeType
  >
void update_x_by_tree_resp_cap(Network &N, const TreeType &tree) {
  for (const auto &depth : boost::adaptors::reverse(tree)) {
    for (const auto &edge : depth) {
      node v;
      arc a;
      std::tie(v, a) = edge;

      auto &nodedata = N.nodedata[v];
      auto &arcdata = N.arcdata[std::abs(a)];

      assert(((nodedata.parent == a) != (arcdata.parent == v)) && "either a is v's parent in the tree or the other way around");
      assert(((N.G.heads[std::abs(a)] == v) != (N.G.tails[std::abs(a)] == v)) && "v must be either a's head or tail");

      if (arcdata.transformed) {
        auto upper = N.G.heads[std::abs(a)] == v;
        auto forward = a > 0;
        auto node_is_parent = nodedata.depth < arcdata.depth;

        decltype(arcdata.xlower) incr;

        if (upper) {
          auto x_old = arcdata.xupper;
          assert(nodedata.potential - arcdata.potential == 0);
          if (forward) {
            if (node_is_parent) {
              arcdata.xupper = std::max(static_cast<decltype(x_old)>(0), x_old + arcdata.deficit);
            } else {
              arcdata.xupper = std::max(static_cast<decltype(x_old)>(0), x_old - nodedata.deficit);
            }
          } else {
            if (node_is_parent) {
              arcdata.xupper = std::max(static_cast<decltype(x_old)>(0), x_old + arcdata.deficit);
            } else {
              arcdata.xupper = std::max(static_cast<decltype(x_old)>(0), x_old - nodedata.deficit);
            }
          }
          incr = arcdata.xupper - x_old;
        } else {
          auto x_old = arcdata.xlower;
          assert(arcdata.cost + nodedata.potential - arcdata.potential == 0);
          if (forward) {
            if (node_is_parent) {
              arcdata.xlower = std::max(static_cast<decltype(x_old)>(0), x_old + arcdata.deficit);
            } else {
              arcdata.xlower = std::max(static_cast<decltype(x_old)>(0), x_old - nodedata.deficit);
            }
          } else {
            if (node_is_parent) {
              arcdata.xlower = std::max(static_cast<decltype(x_old)>(0), x_old + arcdata.deficit);
            } else {
              arcdata.xlower = std::max(static_cast<decltype(x_old)>(0), x_old - nodedata.deficit);
            }
          }
          incr = arcdata.xlower - x_old;
        }
        nodedata.deficit += incr;
        arcdata.deficit -= incr;
      } else {
        auto i = N.G.tails[std::abs(a)];
        auto j = N.G.heads[std::abs(a)];

        assert(arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential == 0);

        auto x_old = arcdata.xlower;

        if (N.nodedata[i].depth > N.nodedata[j].depth) {
          arcdata.xlower = max(static_cast<decltype(x_old)>(0), x_old - N.nodedata[i].deficit);
        } else {
          arcdata.xlower = max(static_cast<decltype(x_old)>(0), x_old + N.nodedata[j].deficit);
        }

        auto flow_incr = arcdata.xlower - x_old;
        N.nodedata[i].deficit += flow_incr;
        N.nodedata[j].deficit -= flow_incr;
      }
    }
  }
}


template<typename Network>
void successive_shortest_path_resp_cap(Network &N) {

  using CapacityType = decltype(N.arcdata.front().capacity);
  using PotentialType = decltype(N.nodedata.front().potential);
  using DeficitType = decltype(N.nodedata.front().deficit);

  auto num_transformed_edges = 0u;

#ifndef INFINITE_CAPACITY
#define INFINITE_CAPACITY std::numeric_limits<CapacityType>::max()
#endif

  for (auto a = 1u; a <= N.G.no_of_edges; ++a) {
    auto &arcdata = N.arcdata[a];
    arcdata.transformed = arcdata.capacity != INFINITE_CAPACITY;
    if (arcdata.transformed) {
      arcdata.deficit = arcdata.capacity;
      ++num_transformed_edges;
    }
  }

  for (auto v = 1u; v <= N.G.no_of_vertices; ++v) {
    auto &nodedata = N.nodedata[v];
    // initialize in_capacity as the sum of the capacities of all ingoing edges
    nodedata.in_capacity = std::accumulate(std::begin(N.G.incident_edges[v]), std::end(N.G.incident_edges[v]), static_cast<CapacityType>(0), [&N](CapacityType in_capacity, arc a) { return in_capacity + (a < 0 && N.arcdata[-a].transformed ? N.arcdata[-a].capacity : 0); });
    nodedata.deficit = nodedata.demand - nodedata.in_capacity;
  }

  do {
    auto v_start = 0u;
    auto v_start_deficit = static_cast<DeficitType>(0);
    auto a_start = 0u;
    auto a_start_deficit = static_cast<DeficitType>(0);

    for (auto v = 1u; v <= N.G.no_of_vertices; ++v) {
      if (N.nodedata[v].deficit < v_start_deficit) {
        v_start = v;
        v_start_deficit = N.nodedata[v].deficit;
      }
    }

    for (auto a = 1u; a <= N.G.no_of_edges; ++a) {
      if (N.arcdata[a].transformed && N.arcdata[a].deficit < a_start_deficit) {
        a_start = a;
        a_start_deficit = N.arcdata[a].deficit;
      }
    }

#ifndef NEDBUG
    auto deficit_sum = std::accumulate(std::begin(N.nodedata) + 1, std::end(N.nodedata), static_cast<DeficitType>(0), [](DeficitType deficit_sum, const typename decltype(N.nodedata)::value_type & nodedata) {
      return deficit_sum + nodedata.deficit;
    });
    deficit_sum = std::accumulate(std::begin(N.arcdata) + 1, std::end(N.arcdata), deficit_sum, [](DeficitType deficit_sum, const typename decltype(N.arcdata)::value_type & arcdata) {
      return deficit_sum + arcdata.deficit;
    });
    assert(deficit_sum == 0);
#endif

    assert(a_start == 0u || N.arcdata[a_start].transformed);

    if (v_start == 0u && a_start == 0u)
      break;

    auto step_count = 0u;
    auto deficit_1 = std::accumulate(std::begin(N.nodedata) + 1, std::end(N.nodedata), static_cast<DeficitType>(0), [](DeficitType deficit_1, const typename decltype(N.nodedata)::value_type & nodedata) {
      return deficit_1 + std::abs(nodedata.deficit);
    });
    deficit_1 = std::accumulate(std::begin(N.arcdata) + 1, std::end(N.arcdata), deficit_1, [](DeficitType deficit_1, const typename decltype(N.arcdata)::value_type & arcdata) {
      return deficit_1 + std::abs(arcdata.deficit);
    });
    std::cout << "starting step " << ++step_count << ", deficit_1 = " << deficit_1 << std::endl;

    using PriorityQueueElementType = std::pair<PotentialType, std::pair<arc, bool>>;
    //auto compare_first = [](const PriorityQueueElementType &lhs, const PriorityQueueElementType &rhs) {
    //  return lhs.first > rhs.first;
    //};
    auto compare_first = [&N](const PriorityQueueElementType & lhs, const PriorityQueueElementType & rhs) {
      if (lhs.first != rhs.first)
        return lhs.first > rhs.first;
      auto lhs_tail = N.G.tails[std::abs(lhs.second.first)];
      auto lhs_head = N.G.heads[std::abs(lhs.second.first)];
      auto rhs_tail = N.G.tails[std::abs(rhs.second.first)];
      auto rhs_head = N.G.heads[std::abs(rhs.second.first)];
      if (lhs_tail != rhs_tail)
        return lhs_tail > rhs_tail;
      return lhs_head > rhs_head;
    };
    auto s_in_n = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);
    auto s_out_n = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);
    auto s_in_a = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);
    auto s_out_a = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);

    for (auto i = 1u; i <= N.G.no_of_vertices; ++i) {
      N.nodedata[i].visited = false;
      N.nodedata[i].depth = -1;
      N.nodedata[i].parent = 0;
    }

    for (auto i = 1u; i <= N.G.no_of_edges; ++i) {
      if (N.arcdata[i].transformed) {
        N.arcdata[i].visited = false;
        N.arcdata[i].depth = -1;
        N.arcdata[i].parent = 0;
      } else {
        assert(N.arcdata[i].deficit == 0);
      }
    }

    bool use_v_start = v_start != 0u && (a_start == 0u || v_start_deficit < a_start_deficit);
    if (use_v_start) {
      std::cout << "  selected node: " << v_start << std::endl;
    } else {
      std::cout << "  selected arc: " << a_start << std::endl;
    }

    auto deficit_s = use_v_start ? v_start_deficit : a_start_deficit;

    if (use_v_start) {
      N.nodedata[v_start].visited = true;
      N.nodedata[v_start].depth = 0;
      update_s_in_s_out_node(N, v_start, s_in_n, s_out_n, s_in_a, s_out_a);
    } else {
      N.arcdata[a_start].visited = true;
      N.arcdata[a_start].depth = 0;
      update_s_in_s_out_arc(N, a_start, s_in_a, s_out_a);
    }

    auto num_visited = 1u;

    // 4 cases
    // node, arc; node is arc's head ===> forward; upper arc ===> arc positive; node is arc's head
    // node, arc; node is arc's tail ===> forward; lower arc ===> arc positive; node is arc's tail
    // arc, node; node is arc's head ===> backward; upper arc ==> arc negative; node is arc's head
    // arc, node; node is arc's tail ===> backward; lower arc ==> arc negative; node is arc's tail
    std::vector<std::vector<std::pair<node, arc>>> tree;

    //auto inner_step_count = 0u;
    while (num_visited < N.G.no_of_vertices + num_transformed_edges && !(s_in_n.empty() && s_in_a.empty() && s_out_n.empty() && s_out_a.empty())) {
      //std::cout << "  starting inner step " << ++inner_step_count << ", deficit_s = " << deficit_s << std::endl;

      auto pop_top = [](decltype(s_in_n) &queue) -> PriorityQueueElementType {
        auto top = queue.top();
        queue.pop();
        return top;
      };

      bool visited;
      DeficitType new_deficit;

      if (deficit_s < 0 || (s_in_n.empty() && s_in_a.empty())) {
        assert(!(s_out_n.empty() && s_out_a.empty()));
        if (!s_out_n.empty() && (s_out_a.empty() || s_out_n.top().first < s_out_a.top().first)) {
          std::tie(visited, new_deficit) = explore_arc(N, pop_top(s_out_n), s_in_a, s_out_a, tree, true);
        } else {
          std::tie(visited, new_deficit) = explore_node(N, pop_top(s_out_a), s_in_n, s_out_n, s_in_a, s_out_a, tree, true);
        }
      } else {
        assert(!(s_in_n.empty() && s_in_a.empty()));
        if (!s_in_n.empty() && (s_in_a.empty() || s_in_n.top().first < s_in_a.top().first)) {
          std::tie(visited, new_deficit) = explore_arc(N, pop_top(s_in_n), s_in_a, s_out_a, tree, false);
        } else {
          std::tie(visited, new_deficit) = explore_node(N, pop_top(s_in_a), s_in_n, s_out_n, s_in_a, s_out_a, tree, false);
        }
      }

      if (visited) {
        deficit_s += new_deficit;
        //if (deficit_s == 0)
        //  break;
        ++num_visited;
      }

#ifdef UPDATE_S_CUTS
      // update deficit vector
      if (out_it_goes == forward) {
        deficit[v] += flow_incr;
        deficit[w] -= flow_incr;
      } else {
        deficit[v] -= flow_incr;
        deficit[w] += flow_incr;
      }
#endif
    }

#ifndef UPDATE_S_CUTS
    update_x_by_tree_resp_cap(N, tree);
#endif
  } while (true);

  //std::cout << "  primal function value = " << N.calculate_primal_objective_value() << std::endl;

  //for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
  //  cout << "  x[" << a << "] = " << N.arcdata[a].xlower << endl;
  //}
}


/** \brief Assigns Costs and Resistances to the Arcs of the Graph
 *
 * @param G The Graph
 *
 */
template <
  typename Network
  >
void assign_costs_resistances_capacities( Network& N ) {
  typedef typename Network::RationalType RationalType;
  Random rg;
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
    RationalType random_cost      = static_cast<RationalType>( rg.rng() );
    RationalType random_capacity  = static_cast<RationalType>( rg.rng() );

    RationalType cap = static_cast<RationalType>( 2 * ceil(random_capacity * 200) + 1 );
    RationalType cost = static_cast<RationalType>( ceil( random_cost * 20 ) );

    N.arcdata[a].capacity = cap;
    N.arcdata[a].cost = cost;
  }
}


template <
  typename Network
  >
void assign_demands( Network& N ) {
  typedef typename Network::RationalType RationalType;
  Random rg;
  RationalType sum_demands(0);
  for (unsigned int v = 1; v < N.G.no_of_vertices; v++) {
    double d = rg.rng();
    N.nodedata[v].demand = ceil(d * 30);
    sum_demands += N.nodedata[v].demand;
  }
  N.nodedata[N.G.no_of_vertices].demand = -sum_demands;
}


/** \brief Computes the maximum capacity and the maximum cost
 *
 * @param max_cap Maximum Capacity
 * @param max_cost Maximum Cost
 * @param capacities Capacities Vector
 * @param costs Costs Vector
 *
 */
template< typename IntegerType, typename RationalType >
void comp_max_cap_cost(RationalType& max_cap,
                       RationalType& max_cost,
                       const vector<RationalType>& capacities,
                       const vector<RationalType>& costs) {
  // assumes costs and capacities have the same length
  assert(capacities.size() == costs.size());

  for (auto cap : capacities) {
    if (cap > max_cap) {
      max_cap = cap;
    }
  }
  for (auto cost : costs) {
    if (cost > max_cost) {
      max_cost = cost;
    }
  }
}

/** \brief Prints the Primal and Dual Solutions
 *
 * @param x Primal Solution
 * @param y Dual Solution
 * @param G The Graph
 *
 */
template< typename IntegerType, typename RationalType >
void print_x_y(const vector<RationalType>& x,
               const vector<RationalType>& y,
               Graph<IntegerType, RationalType> &G
              )
{
  cout << "Vector X: " << endl;

  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    cout << i << ": " << x[i] << endl;
  }

  cout << "Vector Y: " << endl;

  for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
    cout << i << ": " << y[i] << endl;
  }
  cout << endl;
}

/** \brief Prints the Primal and Dual Objective Values
 *
 * @param G The Graph
 * @param x The primal solution
 * @param y The dual solution
 * @param s The Dual Slacks
 *
 */
template<typename IntegerType, typename RationalType>
void print_obj_values_and_assert(Graph<IntegerType, RationalType> &G,
                                 vector<RationalType> &x,
                                 vector<RationalType> &y,
                                 vector<RationalType> &s)
{
  const RationalType cTx = primal_obj_value( G, x );
  const RationalType bTy = dual_obj_value( G, y );
  cout << "primal objective value: ";
  cout << cTx << endl;
  cout << "dual objective value: ";
  cout << bTy << endl;
  cout << "duality gap: " << cTx - bTy << endl;
//   assert( RationalType(1e3) * fabs(cTx - bTy - calculate_duality_gap( G, x, s ) ) < RationalType(1) );
//   assert( cTx - bTy < 1 );

  cout << endl << "optimum = " << floor( cTx )  << " " << ceil( bTy ) << endl;
  //  assert( floor(cTx)  == ceil(bTy)  );

  cout << "assertion done " << endl;
}

template<typename IntegerType, typename RationalType>
void print_obj_values_and_assert(Graph<IntegerType, RationalType> &G0,
                                 Graph<IntegerType, RationalType> &G,
                                 vector<RationalType> &x,
                                 vector<RationalType> &y,
                                 vector<RationalType> &s)
{
  const RationalType cTx = primal_obj_value( G0, x );
  const RationalType bTy = dual_obj_value( G, y );
  cout << "primal objective value: ";
  cout << cTx << endl;
  cout << "dual objective value: ";
  cout << bTy << endl;
  cout << "duality gap: " << cTx - bTy << endl;
//   assert( RationalType(1e3) * fabs(cTx - bTy - calculate_duality_gap( G, x, s ) ) < RationalType(1) );
//   assert( cTx - bTy < 1 );

  cout << endl << "optimum = " << floor( cTx )  << " " << ceil( bTy ) << endl;
  //  assert( floor(cTx)  == ceil(bTy)  );

  cout << "assertion done " << endl;
}


template<typename Network>
void
copy_demands_to_vector(Network& N, vector<typename Network::IntegerType>& demands)
{
  for(unsigned int v = 1; v <= N.G.no_of_vertices; v++)
  {
    demands[v] = N.nodedata[v].demand;
  }
}

template<typename Network>
void
primal_feasibility_check(Network& N, vector<typename Network::IntegerType>& demands) {

  cout << "Primal Feasiblility Check" << endl;
  
  typedef typename Network::IntegerType IntegerType;

  for (unsigned int v = 1; v <= N.G.no_of_vertices; v++) {
    IntegerType flow_into_v = 0;
    for (auto a : N.G.incident_edges[v]) {
      IntegerType real_direction = a;
      if (N.arcdata[abs(a)].direction != 0) {
        real_direction = a * N.arcdata[abs(a)].direction;
      } else {
        real_direction = a;
      }
      if (real_direction > 0) {
        flow_into_v -= N.arcdata[abs(a)].xlower;
      } else {
        flow_into_v -= (N.arcdata[abs(a)].capacity - N.arcdata[abs(a)].xlower);
      }
      if (a > 0) {
        flow_into_v -= N.arcdata[abs(a)].infeasibility;
      } else {
        flow_into_v += N.arcdata[abs(a)].infeasibility;
      }
    }
   
   cout << to_double(flow_into_v - demands[v]) << "< - >" << 0 << endl;
   
   //assert(flow_into_v - demands[v] == IntegerType(0));
  }
}

template<typename Network>
void dual_feasibility_check(Network& N) {
  
  for (unsigned int a = 1; a <= N.G.no_of_edges; a++) {
    if (N.arcdata[a].direction == -1) {
      assert(to_double(N.nodedata[N.G.heads[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].slower +
                       N.arcdata[a].cost) < 1e-3);
    } else {
      assert(to_double(N.nodedata[N.G.tails[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].slower +
                       N.arcdata[a].cost) < 1e-3);
    }

    if (N.arcdata[a].direction == -1) {
      assert( to_double(N.nodedata[N.G.tails[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].supper) < 1e-3);
    } else {
      assert( to_double(N.nodedata[N.G.heads[a]].potential - N.arcdata[a].potentialvw - N.arcdata[a].supper) < 1e-3);
    }

    if (N.arcdata[a].infeasibility != 0) {
      assert(to_double(N.nodedata[N.G.tails[a]].potential - N.nodedata[N.G.heads[a]].potential - N.arcdata[a].sroof +
                       N.arcdata[a].croof ) < 1e-3);
    }
  }
}


template <
  typename Network
  >
void flip_arcs(Network& N) {
  for (unsigned int  a = 1; a <= N.G.no_of_edges; a++)
  {
    if (N.arcdata[a].direction == -1)
    {
      unsigned int temp = N.G.heads[a];
      N.G.heads[a] = N.G.tails[a];
      N.G.tails[a] = temp;
    }
  }

  for (unsigned int v = 1; v <= N.G.no_of_vertices; v++)
  {
    for (auto& a : N.G.incident_edges[v])
    {
      if (N.arcdata[abs(a)].direction == -1) a = -a;
    }
  }
}

/** \brief Computes the Initial Solution after converting the original graph into an auxilliary graph
 *
 * @param G0 The original graph
 * @param G The Auxilliary Graph
 * @param G_delta_wye The Delta-Wye Graph
 * @param x The Primal Solution
 * @param y The Dual Solution
 * @param s The Dual Slacks
 * @param q Parameter
 *
 */
template <
  typename Network,
  typename RationalType
  >
typename Network::IntegerType preprocessing(
  Network& N,
  RationalType &q,
  typename Network::IntegerType t_path_following
) {
  typedef typename Network::GraphType Graph;
  typedef typename Network::IntegerType IntegerType;
  // typedef typename Network::RationalType RationalType;
  const Graph& G0 = N.G;

#ifdef VERBOSE
  cout << "pre_processing() " << endl;
#endif

  SpanningTree< Graph, IntegerType, RationalType> ST(G0);


#ifdef VERBOSE
  cout << "Capacities: " << std::endl;
  G0.print_arc_data(G0.capacities);
  cout << std::endl << "Costs: " << std::endl;
  G0.print_arc_data(G0.costs);
#endif

  vector<IntegerType> z(G0.no_of_edges + 1, IntegerType(0));


  compute_BFS_solution(N, 1, z);

#ifdef VERBOSE
  for (unsigned int i = 1; i <= G0.no_of_edges; i++) {
    cout << "z[" << i << "] = " << z[i] << endl;
  }
#endif


  IntegerType t = find_initial_solution( N, z, q, t_path_following);

  flip_arcs(N);

  return t;
}



/** \brief Runs the ford-fulkerson algorithm to do a max-flow computation
 *
 * @param G The Original Graph
 * @param G1 The Auxilliary Graph
 * @param x The primal solution
 * @param x_integral The primal integral solution
 * @param s_tilde Integral dual slacks
 *
 */
template<typename IntegerType, typename RationalType>
void
ford_fulkerson(
  Graph<IntegerType, RationalType>& G,
  Graph<IntegerType, RationalType>& G1,
  vector<RationalType>& x,
  vector<RationalType>& x_integral,
  vector<RationalType>& s_tilde
)
{
  for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
    G1.demands[i] = G.demands[i];
    G1.incident_edges[i] = G.incident_edges[i];
  }
  for (unsigned int i = 1; i <= G.no_of_edges; i++) {
    G1.capacities[i] = G.capacities[i];
    G1.heads[i] = G.heads[i];
    G1.tails[i] = G.tails[i];
  }

  unsigned int n = G.no_of_vertices;
  node source = G.no_of_vertices + 1;
  node target = G.no_of_vertices + 2;

  G1.demands[source] = 0;
  G1.demands[target] = 0;
#ifdef VERBOSE
  cout << "source: " << source << endl << "target: " << target << endl;
#endif
  G1.no_of_edges = G.no_of_edges;
  for (unsigned int i = 1; i <= n; i++) {
    if (G1.demands[i] < 0) {
#ifdef VERBOSE
      cout << source << " , " << i << endl;
#endif
      arc a = G1.new_edge(source, i);
#ifdef VERBOSE
      cout << "a = " << a << endl;
#endif
      x_integral[a] = 0;
      G1.capacities[a] =  abs(G.demands[i]);
      G1.demands[i] = 0;
    }

    if (G1.demands[i] > 0) {
      arc a = G1.new_edge(i, target);
#ifdef VERBOSE
      cout << "a = " << a << endl;
#endif
      x_integral[a] = 0;
      G1.capacities[a] =  abs(G.demands[i]);
      G1.demands[i] = 0;
    }
  }

#ifdef VERBOSE
  cout << "source " << G1.demands[source] << endl;
  cout << "target " << G1.demands[target] << endl;
  cout << "dem " << endl;
#endif

#ifdef VERBOSE

  for (unsigned int i = 1; i <= G1.no_of_vertices; i++) {
    cout << G1.demands[i] <<  endl;
  }

#endif

#ifdef VERBOSE
  cout << "cap " << endl;
#endif
#ifdef VERBOSE
  for (unsigned int i = 1; i <= G1.no_of_edges; i++) {
    cout << G1.capacities[i] << endl;
  }
#endif
  for (unsigned int i = 1 ; i <= G.no_of_edges; i++) {
    x_integral[i] = 0;
  }

  while (true) {
    vector<int> visited( G1.no_of_vertices + 1, 0 );
    vector<arc> bfstree( G1.no_of_vertices + 1, 0 );
    deque<node> order;
    deque<node> Q;
    Q.push_back( source );
    visited[source] = 1;

    while ( !Q.empty() ) {
      const node v = Q.front();
      if ( v == target) {
        break;
      }
#ifdef VERBOSE
      cout << "v: " << v << endl;
#endif
      for ( auto a : G1.incident_edges[v] ) {
#ifdef VERBOSE
        cout << "a: " << a << endl;
#endif
        if (s_tilde[abs(a)] != 0 && abs(a) <= G.no_of_edges) continue;
        if (a > 0) {
          if (G1.capacities[abs(a)] - x_integral[abs(a)] <= 0 ) continue;
        }
        if (a < 0) {
          if (x_integral[abs(a)] <= 0) continue;
        }
        node w;
        if ( a > 0 ) {
          assert( G1.tails[a] == v );
          w = G1.heads[a];
        } else {
          assert( G1.heads[-a] == v );
          w = G1.tails[-a];
        }

        if ( !visited[w] ) {
          Q.push_back( w );
          visited[w] = 1;
          bfstree[w] = a;

#ifdef VERBOSE
          cout << "w = "  << w << endl;
          cout << "bfstree[w] = " << a << " = " << bfstree[w] << endl;
#endif
        }
      }
      Q.pop_front();
      order.push_front( v );
    }

#ifdef VERBOSE
    cout << "Q.front() = " << target << endl;
#endif
    if (Q.front() != target) {
      break;
    }

    node u = target;
    vector<arc> path;
    while (true) {
      arc a = bfstree[u];
      path.push_back(a);
      node w;
      if (G1.head(a) == u) {
        w = G1.tail(a);
      }
      else {
        w = G1.head(a);
      }
      if (w == source) {
        break;
      }
      u = w;
    }
    unsigned int del_P = numeric_limits<unsigned int>::max();
    for (auto a1 : path) {
      if (a1 > 0) {
        if (G1.capacities[a1] - x_integral[a1] < del_P) {
          del_P = G1.capacities[a1] - x_integral[a1];
        }
      }
      if (a1 < 0) {
        if (x_integral[-a1] < del_P) {
          del_P = x_integral[-a1];
        }
      }
    }
    for (auto a : path) {
#ifdef VERBOSE
      cout << " path " << a << endl;
#endif

#ifdef VERBOSE
      cout << "del_P = " << del_P << endl;
#endif
      if (a > 0) {
        x_integral[a] += del_P;
      }
      else {
        x_integral[-a] -= del_P;
      }
    }
  }
#ifdef VERBOSE
  cout << " x_integral : " << endl;
  for (auto a : x_integral) {
    cout << a << " , " ;
  }
#endif
}

/** \brief Does a post processing
 *
 * @param G The Graph
 * @param x The Primal Solution
 * @param y The Dual Solution
 * @param s The Dual Slacks
 * @param y_tilde Integral Dual Solution
 * @param s_tilde Integral Dual Slacks
 *
 */
template <
  typename Network,
  typename RationalType
  >
void postprocessing(
  Network& N,
  vector<RationalType> &x,
  vector<RationalType> &y,
  vector<RationalType> &s,
  vector<RationalType> &y_tilde,
  vector<RationalType> &s_tilde,
  const RationalType& THRESHOLD_X
) {
  typedef typename Network::GraphType Graph;
//  typedef typename Network::IntegerType IntegerType;
//  typedef typename Network::RationalType RationalType;
  const Graph& G = N.G;

  vector<bool> inside_set_S (G.no_of_vertices + 1, false);
  std::vector<int> visited(G.no_of_vertices + 1 , 0);
  for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
    for (unsigned int i = 1; i <= G.no_of_vertices; i++) {

      cout << "inside_set_S[" << i << "] = " << inside_set_S[i] << endl;
    }
    if (!inside_set_S[i]) {
      cout << "i = " << i << endl;

      rounding_scheme( N, i, x, y, s, THRESHOLD_X, visited, inside_set_S, y_tilde);
    }
  }
  s_tilde = compute_dual_slacks( G, y_tilde );

  for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
#ifdef VERBOSE
    cout << y_tilde[i] << " " << y_tilde[i] - floor( y_tilde[i] + 1e-6 ) << endl;
#endif
    assert( RationalType(1e6) * ( y_tilde[i] - floor( y_tilde[i] + 1e-6 )) < RationalType(1) );
  }
}

template<typename IntegerType, typename RationalType>
void check_slack_null(Graph<IntegerType, RationalType> &G,
                      vector<RationalType> &y
                     )
{
  for (auto am : G.arc_map ) {

    RationalType s_ahat = G.costs[am.ahat] - y[G.tails[am.ahat]] + y[G.heads[am.ahat]];
    assert(s_ahat == 0);
  }
}

/** \brief Reconstructs the primal solution for the Original Graph
 *
 * @param G0 The Original Graph
 * @param G1 The Auxilliary Graph
 * @param x0 The Re-constructed solution
 * @param x1 Solution on the Auxilliary graph
 *
 */
template<typename IntegerType, typename RationalType>
void reconstruct_orig_solution(Graph<IntegerType, RationalType> &G0,
                               Graph<IntegerType, RationalType> &G1,
                               vector<RationalType> &x0,
                               vector<RationalType> &x1
                              )
{
  for (unsigned int i = 0; i <= G0.no_of_edges; ++i) {
#ifdef VERBOSE
    cout << G1.arc_map[i].orig << " ";
#endif
    if (G1.arc_map[i].orig != 0) {
      x0[G1.arc_map[i].orig] = x1[G1.arc_map[i].vvw];
    }
  }
}

/** \brief Reconstructs the Dual Slacks for the Original Graph
 *
 * @param G0 Original Graph
 * @param G1 The Auxilliary Graph
 * @param s_tilde0 The Reconstructed Dual Slacks
 * @param s_tilde Dual Slacks on the Auxilliary Graph
 *
 */
template<typename IntegerType, typename RationalType>
void reconstruct_orig_solution_s(Graph<IntegerType, RationalType> &G0,
                                 Graph<IntegerType, RationalType> &G1,
                                 vector<RationalType> &s_tilde0,
                                 vector<RationalType> &s_tilde
                                )
{
  for (unsigned int i = 0; i <= G0.no_of_edges; ++i) {
#ifdef VERBOSE
    cout << G1.arc_map[i].orig << " ";
#endif
    if (G1.arc_map[i].orig != 0) {
      s_tilde0[G1.arc_map[i].orig] = s_tilde[G1.arc_map[i].vvw];
    }
  }
}

#endif
