#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include "graph.h"
#include "simple_solver.h"
#include "random.h"
#include "potential_reduction.h"

using namespace std;

template<class T>
void compute_BFS_solution(Graph<T> &G, node s, std::vector<T>&z)
{
  vector<T>(G.no_of_edges+1,0).swap( z );
  
  vector<T> b( G.demands );
  cout<<std::endl<<std::endl;
  
  #ifdef VERBOSE
  cout<<"-----------------DEMANDS--------------------"<<std::endl; 
  for(unsigned int i=1; i<=G.no_of_verticies; i++)
  {
    cout<<b[i]<<" , ";
  }
  cout << endl << endl;
  #endif
  
  vector<int> visited( G.no_of_verticies+1,0 );
  vector<arc> bfstree( G.no_of_verticies+1, 0 );
  deque<node> order;
  deque<node> Q;
  Q.push_back( s );
  visited[s] = 1;
  while( !Q.empty() )
  {
    const node v = Q.front();
    
    for( auto a : G.incident_edges[v] ) {
      node w;
      if( a > 0 ) {
        assert( G.tails[a] == v );
        w = G.heads[a];
      } else {
        assert( G.heads[-a] == v );
        w = G.tails[-a];
      }
      if( !visited[w] ) {
        Q.push_back( w );
        visited[w] = 1;
        bfstree[w] = a;
      }
    }
    Q.pop_front();
    order.push_front( v );
  }
  #ifndef NDEBUG
  const node r = order.back();
  #endif
  assert( bfstree[r] == 0 );
  order.pop_back();
  for( auto v : order ) {
    const arc a = bfstree[v];
    #ifdef VERBOSE
    cout << "demand before: " << b[v] << endl;
    #endif
    G.pull_flow( -a, b[v], b );
    #ifdef VERBOSE
    cout << "demand after: " << b[v] << endl;
    #endif
    assert( fabs( b[v] ) < 1e-10 );
    const arc aa = abs( a );
    z[aa] = G.flows[aa];
  }
  assert( b[r] == 0 );
  #ifdef VERBOSE
  cout << endl << "The Z vector" << endl;
  for(unsigned int i=1; i<=G.no_of_edges; i++)
  {
    cout << "z " << i << ": " << z[i] << endl;
  }
  #endif
  primal_sanity_check( G, z );
}


template<class T>
void find_initial_solution_new( const Graph<T> &G0, Graph<T> &G1, unsigned int n, std::vector<double>& x, std::vector<double> &y, std::vector<double> &s, std::vector<double> &z, T q)
{ 
  const unsigned int nedges_G0 = G0.no_of_edges;
  const double t = determine_t(G0,z,sqrt(nedges_G0));

  #ifdef VERBOSE
  cout << "G0 looks like: " << endl;
  unsigned int i,j;
  for(i=1; i<=G0.no_of_verticies; i++)
  {
    cout<<"Vertex "<<i; cout<<" :";
    for(j=0; j<G0.incident_edges[i].size(); j++ )
    {
      cout<<G0.incident_edges[i][j];
      cout<<",";
    }
    cout<<std::endl;  
  }  
  cout << "z s of G0 look like: " << endl;
  for(unsigned int i=1; i<=nedges_G0;i++){
    cout << "z[" << i << "]: " << z[i] << ", " << "cap[" << i << "]: " << G0.capacities[i] << ",";
  }
  cout << endl;
  #endif
  
  for(unsigned int i=1; i<=n; i++){
    G1.demands[i] = G0.demands[i];
  }
  x.resize(3 * nedges_G0 + 1);
  s.resize(3 * nedges_G0 + 1);
  y.resize(n + nedges_G0 + 1);
  
  for(unsigned int i=1; i<=nedges_G0; i++){
    
    node vw = G1.create_node();
    node v  = G0.tails[i];
    node w  = G0.heads[i];
    
    T cap  = G0.capacities[i];
    T cost = G0.costs[i];
    
    assert(cap!=0);
    
    // G1.demands[v]  += G0.demands[v];
    G1.demands[w]  +=  - cap;
    G1.demands[vw] += cap;   
    
    arc vvw = G1.new_edge(v, vw);
    arc wvw = G1.new_edge(w, vw);
    
    G1.costs[vvw] = cost;
    G1.costs[wvw] = 0;
    
    x[vvw] = cap/2.0;
    x[wvw] = cap/2.0;
    
    y[v]  = 0.0;
    y[w]  = 0.0;
    y[vw] = -2.0*t/cap;
    
    s[vvw] = G1.costs[vvw] + y[v] - y[vw];
    s[wvw] = G1.costs[wvw] + y[w] - y[vw];
    
    if (z[i] > cap/2.0){
      arc vw = G1.new_edge(v,w);
      G1.costs[vw] = ceil( t / fabs( z[i] - cap/2.0 ) ); 
      assert(G1.costs[vw]!=0);
      x[vw] = fabs( z[i] - cap/2.0 );
      s[vw] = G1.costs[vw];
    }
    if (z[i] < cap/2.0){
      arc wv = G1.new_edge(w,v);
      G1.costs[wv] = ceil( t / fabs( z[i] - cap/2.0 ) ); 
      assert(G1.costs[wv]!=0);
      x[wv] = fabs( z[i] - cap/2.0 );
      s[wv] = G1.costs[wv];
    }
    if (z[i] == cap/2.0){
      cout << "I do nothing z=u/2..." << endl;
    }
    
  }
  
  // checks and asserts
  sanity_check( G1, x, y, s );
  #ifdef VERBOSE  
  T duality_gap = 0.0;
  T log_x_s     = 0.0;
  for(unsigned int i=1; i<=G1.no_of_edges; i++)
  {
    cout << "s[" << i << "]: " << s[i] << ", ";
    assert(x[i]>0);
    assert(s[i]>0);
    const T xs = x[i]*s[i];
    duality_gap += xs;
    log_x_s += log(xs);
  }
  for(unsigned int j=1; j<=G1.no_of_verticies; j++){
    cout << "b[" << j << "]=" << G1.demands[j] << ", ";
  }
  
  const unsigned int m = G1.no_of_edges;
  const unsigned int q = m + static_cast<unsigned int>( ceil( sqrt(m) ) );
  const double mlogm = m*log(m);
  cout<<"Duality Gap: "<<duality_gap << std::endl;
  cout<<"Potential: " << q << " * " << log(duality_gap) << " - " << log_x_s << " - " << mlogm << " = " << q*log(duality_gap)-log_x_s - mlogm << endl;
  
  cout << "G1 looks like: " << endl;
  // G1.print();
  for(i=1; i<=G1.no_of_verticies; i++)
  {
    cout<<"Vertex "<<i; cout<<" :";
    for(j=0; j<G1.incident_edges[i].size(); j++ )
    {
      cout<<G1.incident_edges[i][j];
      cout<<",";
    }
    cout<<std::endl;  
  }
  
  cout << "x s of G1 look like: " << endl;
  for(unsigned int i=1; i<=G1.no_of_edges;i++){
    cout << x[i] << ", ";
  }
  cout << endl;
  #endif
  
}

template<class T>
Graph<T> find_initial_solution(Graph<T> &G, int n, double t, std::vector<double>& x, std::vector<double> &y, std::vector<double> &s, std::vector<double> &z)
{
  unsigned int new_count = G.no_of_edges*3; //new_count represents the number of edges in the graph G_1
  unsigned int no_of_edges_g0 = G.no_of_edges; //no_of_edges represents the number of edges in the graphc G_0
  x.resize(new_count+1);
  y.resize(n + no_of_edges_g0 + 1); 
  
  unsigned int no_of_verticies = n;
  
  for(unsigned int i=1; i<=no_of_edges_g0; i++)
  {   
    unsigned int new_node = G.create_node();
    #ifdef VERBOSE
    cout << "new_node " << new_node << endl;
    #endif
    assert( new_node == no_of_verticies + i );
    G.demands[new_node]= G.capacities[i];   
    node head = G.heads[i];
    G.demands[head] = G.demands[head] - G.capacities[i];
    //  cout<<head<< "Demand Changed"<<endl;
    node tail = G.tails[i]; 
    G.new_edge(head, new_node);
    #ifdef VERBOSE
    cout << no_of_edges_g0+2*i-1 << " " << G.tails[no_of_edges_g0+2*i-1] << " " << G.heads[no_of_edges_g0+2*i-1] << endl;
    #endif
    assert( G.heads[no_of_edges_g0+2*i-1] == new_node );
    
    G.new_edge(tail, new_node);
    
    #ifdef VERBOSE
    cout << no_of_edges_g0+2*i << " " << G.tails[no_of_edges_g0+2*i] << " " << G.heads[no_of_edges_g0+2*i] << endl;
    #endif
    
    assert( G.heads[no_of_edges_g0+2*i] == new_node );
    
    x[no_of_edges_g0+2*i-1]= G.capacities[i]/2.0;
    x[no_of_edges_g0+2*i] = G.capacities[i]/2.0;
    G.costs[no_of_edges_g0+ 2*i -1] = 0;
    G.costs[no_of_edges_g0 +2*i] = G.costs[i]; 
    if(z[i] > G.capacities[i]/2.0)
    {
      G.costs[i] = ceil(t/fabs(z[i] - (G.capacities[i]/2.0))); 
      x[i] = fabs((z[i] - (G.capacities[i]/2.0)));
      y[G.heads[i]] = 0;
      y[G.tails[i]] = 0;
      y[no_of_verticies + i] = -2*t/G.capacities[i];  
      //cout<<y[no_of_verticies+i]<<"	"<<G.capacities[i]<<"	"<<std::endl;     
    }
    else
    {
      G.costs[i] = ceil(t/fabs(z[i] - (G.capacities[i]/2.0)));     
      G.flip_arc(i);
      int cap = G.capacities[i];
      if (cap % 2==0) {
        #ifdef VERBOSE
        cout << "capacity is " << cap << endl;
        #endif
      }
      assert(cap % 2!=0);
      x[i] = fabs((z[i] - (G.capacities[i]/2.0)));
      y[G.heads[i]] = 0;
      y[G.tails[i]] = 0;
      y[no_of_verticies + i] = -2*t/G.capacities[i];     
      //cout<<y[no_of_verticies+i]<<"	"<<G.capacities[i]<<"	"<<std::endl;
    }  
  } 
  assert( G.no_of_edges == 3*no_of_edges_g0 );
  s.resize(G.no_of_edges+1);
  cout<<std::endl<<std::endl;
  for(unsigned int i=1; i<= G.no_of_edges; i++)
  {
    s[i] = G.costs[i] + y[G.tails[i]] - y[G.heads[i]];
    #ifdef VERBOSE
    cout << i << " = (" << G.tails[i] << "," << G.heads[i] << "): " << x[i]*s[i] << "; " << x[i] << ", " <<s[i]<< " = " <<G.costs[i]<< " + " << y[G.tails[i]] << " - " << y[G.heads[i]] << std::endl;
    #endif
    assert( s[i] > 0 );
    assert( x[i] > 0 );
  } 
  
  sanity_check( G, x, y, s );
  
  double duality_gap = 0.0;
  double log_x_s = 0.0;
  for(unsigned int i=1; i<=G.no_of_edges; i++)
  {
    const double xs = x[i]*s[i];
    duality_gap += xs;
    log_x_s += log(xs);
  }
  
  
  
  #ifdef VERBOSE
  const unsigned int m = G.no_of_edges;
  const unsigned int q = G.no_of_edges + static_cast<unsigned int>( ceil( sqrt(m) ) );
  const double mlogm = m*log(m);
  cout<<"Duality Gap: "<<duality_gap << std::endl;
  cout<<"Potential: " << q << " * " << log(duality_gap) << " - " << log_x_s << " - " << mlogm << " = " << q*log(duality_gap)-log_x_s - mlogm << endl;
  #endif
  return G;
}



template<typename T>
T update_duality_gap(Graph<T>&G, std::vector<double> &x, std::vector<double> &s)
{
  T duality_gap=0;
  
  for(unsigned int i=1; i<=G.no_of_edges; i++)
  {
    const T xs = x[i]*s[i];
    #ifdef VERBOSE
    cout<< i << ": " << x[i]<<" * " << s[i] << " = " << xs << endl;
    #endif
    duality_gap += xs;  
  }
  
  return duality_gap;
}

template<typename T>
T b( const Graph<T> &G, const std::vector<node>& S_k)
{
  T demand = 0.0;
  for(unsigned int i=0; i<S_k.size(); i++)
  {
    demand+=G.demands[S_k[i]];
  }
  return demand;
}

template<typename T>
std::vector<arc> del_in( const Graph<T> &G, const std::vector<node> S_k)
{
  std::vector<arc> incoming_arcs;
  for(unsigned int i=0; i<S_k.size(); i++)
  {
    for(unsigned int j=0; j< G.incident_edges[i].size(); j++)
    {
      if(G.incident_edges[i][j] < 0)
      {
        incoming_arcs.push_back(G.incident_edges[i][j]);
      }
    }
  }
  return incoming_arcs;
}


template<typename T>
std::vector<arc> del_out( const Graph<T>&G, const std::vector<node> S_k)
{
  std::vector<arc> outgoing_arcs;
  for(unsigned int i=0; i<S_k.size(); i++)
  {
    for(unsigned int j=0; j< G.incident_edges[i].size(); j++)
    {
      if(G.incident_edges[i][j] > 0)
      {
        outgoing_arcs.push_back(G.incident_edges[i][j]);
      }
    }
  }
  return outgoing_arcs;
}

template<typename T>
T minimal( const Graph<T>& G, const std::vector<arc>& arcs, const std::vector<T>& y, const std::vector<T>& y_0, arc& a_k)
{
  T min = std::numeric_limits<T>::max();;
  a_k = 0;
  for(auto a: arcs)
  {
    int w = G.head(a);
    int v = G.tail(a);
    if(min > G.costs[a] + y[v] - y_0[w])
    {
      min = G.costs[a] + y[v] - y_0[w];
      a_k = a;
      
    }
  }
  return min;
}





struct node_data{
  arc a;  
  node v;
  double data;
};

class compare_data {
public:
  bool operator()(const node_data& d1, const node_data& d2) const
  {
    if(d1.data > d2.data)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  
};


template<typename T>
std::vector<T> rounding_scheme( const Graph<T>&G, const node& s, const std::vector<T> &x, const std::vector<T> &y_0, const std::vector<T> &s_0)
{
  priority_queue<node_data, vector<node_data>, compare_data> incoming_arcs_pq;
  priority_queue<node_data, vector<node_data>, compare_data> outgoing_arcs_pq;
  
  std::vector<T> in(G.no_of_verticies + 1, numeric_limits<T>::max());
  std::vector<T> out(G.no_of_verticies + 1, numeric_limits<T>::max());
  std::vector<int> visited(G.no_of_verticies+1 , 0);
  #ifdef VERBOSE
  cout << "rounding_scheme begins..." << endl;
  cout << "number of verticies: " << G.no_of_verticies << endl;
  #endif
  std::vector<T> y(G.no_of_verticies + 1);
  std::vector<node> S;
  y[s] =0;
  T delta_k = -y_0[s];
  S.push_back(s);
  visited[s]=1;
  
  for(auto a: G.incident_edges[s])
  {
    if(a < 0)
    {
      node w = G.head(a);
      double data = G.costs[-a] + y_0[w] - y[s];
      
      if(data < in[w])
      {
        in[w] = data;    
        node_data w_v_data = {a, w, data};
        incoming_arcs_pq.push(w_v_data);
        #ifdef VERBOSE
        cout << "arc: " << a << " tail: " << s << " " << G.head(a) << " " << w << " " <<in[w] << " "<< data << endl; 
        #endif
        
      }
    }
    if(a > 0)
    {
      node w = G.head(a);
      double data = G.costs[a] + y[s] - y_0[w];
      
      if(data < out[w])
      {
        out[w] = data; 
        node_data v_w_data = {a, w, data};
        outgoing_arcs_pq.push(v_w_data);
        #ifdef VERBOSE
        cout << "arc: " << a << " tail: " << s << " " << G.head(a) << " " << w << " " <<out[w] << " "<< data << endl; 
        #endif
        
      }
    }
  }
  
  
  for(unsigned int k =1; k< G.no_of_verticies; k++)
  {
    #if !defined(NDEBUG) || defined(VERBOSE)
    arc a_k = 0;
    node v_k= 0;
    #endif
    node w_k =0;
    T b_S_k = b(G,S);
    if(b_S_k < 0||incoming_arcs_pq.empty())
    { 
      assert(!outgoing_arcs_pq.empty());
      
      node_data minimal_data = outgoing_arcs_pq.top();
      #ifdef VERBOSE
      cout << "-ve branch:  " << minimal_data.a << " " <<minimal_data.v << " " << minimal_data.data << " " << visited[minimal_data.v] << " " << out[minimal_data.v] << endl;
      #endif
      delta_k = minimal_data.data;
      if(delta_k!=out[minimal_data.v]||visited[minimal_data.v]==1)
      {
        outgoing_arcs_pq.pop();
        --k;
        continue;
      }
      
      w_k = minimal_data.v;
      visited[w_k]=1;
      #if !defined(NDEBUG) || defined(VERBOSE)
      a_k = minimal_data.a;
      v_k = G.tail(a_k);
      #endif
      assert( G.head(a_k)==w_k ); 
      #ifdef VERBOSE
      cout << "delta_k: " << delta_k <<endl;
      cout << "v_k: " << v_k << endl;
      cout << "w_k: " << w_k << endl;
      cout << "a_k: " << a_k << endl;
      #endif
    }
    else
    { 
      assert(!incoming_arcs_pq.empty());
      node_data minimal_data = incoming_arcs_pq.top();
      #ifdef VERBOSE
      cout << "+ve branch:  " << minimal_data.a << " " <<minimal_data.v << " " << minimal_data.data << " " << visited[minimal_data.v] << " " << in[minimal_data.v] << endl;
      #endif
      delta_k = -minimal_data.data;
      if(-delta_k!=in[minimal_data.v]||visited[minimal_data.v]==1)
      {
        incoming_arcs_pq.pop();
        --k;
        continue;
      }
      w_k = minimal_data.v;
      visited[w_k] =1;
      #if !defined(NDEBUG) || defined(VERBOSE)
      a_k = minimal_data.a;
      v_k = G.tail(a_k);
      #endif
      assert(G.head(a_k) == w_k);
      #ifdef VERBOSE
      cout << "delta_k: " << delta_k <<endl;
      cout << "v_k: " << v_k << endl;
      cout << "w_k: " << w_k << endl;
      cout << "a_k: " << a_k << endl;
      #endif
    }
    assert( w_k != 0 );
    y[w_k] = y_0[w_k] + delta_k;
    #ifdef VERBOSE
    cout << w_k << ": " << y[w_k] << " = "<< y_0[w_k] << " + " << delta_k << endl;   
    cout << endl << "assertion val: " <<G.costs[abs(a_k)] + y[v_k] - y[w_k] << " = " << G.costs[abs(a_k)] << " + " << y[v_k] << " - " << y[w_k]<< endl;
    #endif
    assert( G.costs[abs(a_k)] + y[v_k] - y[w_k] > -1e-6 );
    S.push_back(w_k);
    
    for(auto a: G.incident_edges[w_k])
    {
      if(a < 0)
      {
        
        node w = G.head(a);
        
        double data = G.costs[-a] + y_0[w] - y[w_k];
        if(data < in[w]&&visited[w]==0)
        {
          node_data w_v_data = {a, w, data};
          incoming_arcs_pq.push(w_v_data);
          in[w] = data;
          #ifdef VERBOSE
          cout << "updated: " << in[w] << endl;
          #endif
        }
        #ifdef VERBOSE
        cout << "arc: " << a << " tail: " << w_k << " " << G.head(a) << " " << w << " " <<in[w] << " "<< data << endl; 
        #endif
      }
      if(a > 0)
      {
        node w = G.head(a);
        double data = G.costs[a] + y[w_k] - y_0[w];
        if(data < out[w]&&visited[w]==0)
        { 
          node_data v_w_data = {a, w, data};
          outgoing_arcs_pq.push(v_w_data);
          out[w] = data;    
          #ifdef VERBOSE
          cout << "updated: " << out[w] << endl;
          #endif   
        }
        #ifdef VERBOSE 
        cout << "arc: " << a << " tail: " << w_k  << " " << G.tail(a) << " " << w << " " <<out[w] << " "<< data << endl; 	   
        #endif	  
      }
    }
    #ifdef VERBOSE
    cout << "s_k+1 later " << S.size() << endl;
    #endif
  }
  
  T b_T_y0 = 0.0;
  T b_T_y = 0.0;
  for(unsigned int i=1; i<=G.no_of_verticies; i++)
  {
    b_T_y0 += G.demands[i]*y_0[i];
    b_T_y += G.demands[i]*y[i]; 
  }
  assert(b_T_y >= b_T_y0 );
  return y;
}

template< typename T > 
vector<T> compute_dual_slacks( const Graph<T>& G, const vector<T>& y ) {
  vector<T> s( G.no_of_edges+1 );
  for( unsigned i = 1; i <= G.no_of_edges; ++i ) {
    s[i] = G.costs[i] + y[G.tails[i]] - y[G.heads[i]];
    #ifdef VERBOSE
    cout << s[i] << endl;
    #endif
    assert( s[i] > -1e-6 );
  }
  return s;
}

template< typename T > 
T primal_obj_value( const Graph<T>& G, const vector<T>& x ) {
  T obj = 0;
  for( unsigned i = 1; i <= G.no_of_edges; ++i ) {
    const T caxa = G.costs[i]*x[i];
    obj += caxa;
    #ifdef VERBOSE
    cout << G.costs[i] << " * " << x[i] << " = " << caxa << endl;
    #endif
  }
  return obj;
}

template< typename T > 
T dual_obj_value( const Graph<T>& G, const vector<T>& y ) {
  T obj = 0;
  for( unsigned i = 1; i <= G.no_of_verticies; ++i ) {
    const T bvyv = G.demands[i]*y[i];
    obj += bvyv;
    #ifdef VERBOSE
    cout << G.demands[i] << " * " << y[i] << " = " << bvyv << endl;
    #endif
  }
  return obj;
}

template< typename T >
T determine_t( const Graph<T>& G, const vector<T>& z, const T& p ) {
  T sum = 0;
  for( unsigned i = 1; i <= G.no_of_edges; ++i ) {
    sum += G.capacities[i]*G.costs[i] + fabs( z[i] - G.capacities[i]/2.0 );
  }
  return sum/p;
}


int read_n(string filepath)
{
  int n = 0;
  int nodeid;
  ptree tree;
  read_xml(filepath, tree);  
  const ptree & graphml = tree.get_child("graphml", empty_ptree());
  const ptree & graph   = graphml.get_child("graph", empty_ptree());
  
  BOOST_FOREACH(const ptree::value_type & nore, graph){
    const ptree & nore_attrs = nore.second.get_child("<xmlattr>", empty_ptree());
    BOOST_FOREACH(const ptree::value_type & nore_attr, nore_attrs){
      if (strncmp(nore_attr.first.data(),"id",2) == 0){
        nodeid = stoi(nore_attr.second.data());  
        n = max(n,nodeid);
      }
    }
  }
return n;
}


int main()
{
  // Comment that in for creating random graphs
  #ifdef READGRAPH
    cout << "Enter filepath to input graph.. " << endl;
    string filepath;
    cin >> filepath; 
    int n = read_n(filepath);
    Graph<double> G0(n);
    G0.read_graph(filepath);
  #else
    cout << "Enter the number of nodes in the graph.. " << endl;
    int n = 0;
    cin>>n;
    Graph<double> G0(n);
    G0.create_graph();
  #endif
  
  
  SpanningTree<Graph<double>,double> ST(G0);
  
  G0.create_low_stretch_tree( ST );
  //G0.generate_cycles();
  //G0.print();

  #ifdef VERBOSE
    cout<<"Capacities: "<<std::endl;
    for(unsigned int i=1; i<=G0.no_of_edges; i++)
    {
      cout<<G0.capacities[i]<<" , ";
    } 
    
    cout<<std::endl<<"Costs: "<<std::endl;
    for(unsigned int i=1; i<=G0.no_of_edges; i++)
    {
      cout<<G0.costs[i]<<" , ";
    }
  #endif

  double max_cap = 0.0;
  double max_cost = 0.0;   
  for(unsigned int i=1; i<=G0.no_of_edges; i++)
  {
    if(G0.capacities[i]>max_cap)
    {
      max_cap = G0.capacities[i];
    }
    
    if(G0.costs[i]>max_cost)
    {
      max_cost = G0.costs[i];
    }
  }

  
  cout<<std::endl<<std::endl;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> s;
  std::vector<double> s_hat_prime;
  std::vector<double> s_prime;
  std::vector<double> x_prime;
  std::vector<double> z;
  
  compute_BFS_solution(G0, 1, z);
  #ifdef VERBOSE
    for(unsigned int i=1; i<=G0.no_of_edges; i++)
    {
      cout<<i<<" : "<< z[i]<< endl;     
    }
  #endif
  // vector<double> y0( G.no_of_verticies+1, 0 );
  vector<double> y_start( G0.no_of_verticies+1, 0 ); 
  //    for( node s = 1; s <= G0.no_of_verticies; ++s ) {
  //      y_start = rounding_scheme(G0, s, x, y0, G0.costs );
  //      cout << "initial rounding from " << s << " yields: " << dual_obj_value( G, y_start ) << endl;
  //    }
   
  cout << "initial rounding yields: " << dual_obj_value( G0, y_start ) << endl;
  
  for(unsigned int i=1; i<=G0.no_of_edges; i++)
  {
    G0.costs[i] = G0.costs[i] + y_start[G0.tail(i)] - y_start[G0.head(i)];
  }
   
   
  const double q = 3*G0.no_of_edges + ceil( sqrt(3*G0.no_of_edges) ); 
  
  Graph<double> G(n);
  find_initial_solution_new(G0, G, n , x, y, s, z, q); 
  
  #ifdef VERBOSE
    cout << "Vector X: " << endl;
  
    for(unsigned int i=1; i<=G.no_of_edges; i++)
    {
      cout << i << ": " << x[i] << endl;
    }
  
    cout << "Vector Y: " << endl;
  
    for(unsigned int i=1; i<=G.no_of_verticies; i++)
    {
      cout<<y[i]<<" , ";
    }
  #endif

  #ifdef VERBOSE
     cout << endl;
  #endif
  sanity_check(G,x,y,s);
  
  #ifdef VERBOSE
     cout<<endl; 
  #endif
   
  potential_reduction_algorithm(G, x, y, s, q);  
  sanity_check(G,x,y,s);
  const double cTx = primal_obj_value( G, x );
  const double bTy = dual_obj_value( G, y );
  cout << "primal objective value: " << cTx << endl;
  cout << "dual objective value: " << bTy << endl;
  cout << "duality gap: " << cTx - bTy << endl;
  assert( fabs( cTx - bTy - calculate_duality_gap( G, x, s ) ) < 1e-3 );
  assert( cTx - bTy < 1 );
  
  cout << endl << "optimum = " << floor( cTx ) << endl;
  assert( floor( cTx ) == ceil( bTy ) );
  
 
  vector<double> y_tilde = rounding_scheme(G, 1, x, y, s); 
  vector<double> s_tilde = compute_dual_slacks( G, y_tilde );
  for(unsigned int i=1; i<= G.no_of_verticies; i++)
  {
    #ifdef VERBOSE
      cout << y_tilde[i] << " " << y_tilde[i] - floor( y_tilde[i]+1e-6 ) << endl;
    #endif
    assert( y_tilde[i] - floor( y_tilde[i]+1e-6 ) < 1e-6 );
  }
 
  sanity_check(G,x,y_tilde,s_tilde);
  const double bTy_tilde = dual_obj_value( G, y_tilde );
  cout << "primal objective value: " << cTx << endl;
  cout << "dual objective value: " << bTy_tilde << endl;
  cout << "duality gap: " << cTx - bTy_tilde << endl;
  assert( fabs( cTx - bTy_tilde - calculate_duality_gap( G, x, s_tilde ) ) < 1e-3 );
  assert( cTx - bTy_tilde < 1 );
 
  return 0;
}

