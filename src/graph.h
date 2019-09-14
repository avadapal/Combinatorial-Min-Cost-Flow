#ifndef GRAPH_H
#define GRAPH_H

#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<fstream>
#include<limits>
#include <boost/range/adaptor/reversed.hpp>
#include "random.h"
#include "unbounded_integers.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/rational.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

using namespace std;
using namespace boost::property_tree;

const ptree& empty_ptree() {
	static ptree t;
	return t;
}

typedef unsigned int node;
typedef int arc;

struct TreeDecompositionInfo {
	int size;
	node s;
	node d;
	TreeDecompositionInfo() : size(0), s(0), d(0) {}
};

struct ArcConstruction {
	arc orig;
	arc vvw;
	arc wvw;
	arc ahat;
	ArcConstruction() : orig(0), vvw(0), wvw(0), ahat(0) {}
};

struct NodeInformation {
	unsigned int tree_index;
	int above_d;
	node lca;
	NodeInformation() : tree_index(0), above_d(-2), lca(0) {}
};


template <
  typename IntegerType,
  typename RationalType
  >
struct BasicNodeData {
	BasicNodeData() :
		demand(0),
		potential(0) {}

	// input parameter
	IntegerType demand;

	IntegerType potential;
};

template <
  typename IntegerType,
  typename RationalType
  >
struct InteriorPointMethodNodeData : public BasicNodeData<IntegerType, RationalType> {
};


template <
  typename IntegerType,
  typename RationalType
  >
struct SSPVariantNodeData : public BasicNodeData<IntegerType, RationalType> {
	SSPVariantNodeData() :
		BasicNodeData<IntegerType, RationalType>(),
		deficit(0),
#ifdef USE_IMBALANCE
		imbalance(0),
#endif
#ifdef RESTORE_BALANCED_NODES
		deficit_delta(0),
#endif
		depth(-1),
		parent(0),
		visited(false) {}

	IntegerType deficit;
#ifdef USE_IMBALANCE
	IntegerType imbalance;
#endif
#ifdef RESTORE_BALANCED_NODES
	IntegerType deficit_delta;
#endif
	int depth;
	arc parent;
	bool visited;
};


template <
  typename IntegerType,
  typename RationalType
  >
struct SSPVariantRCNodeData : public SSPVariantNodeData<IntegerType, RationalType> {
	SSPVariantRCNodeData() :
		SSPVariantNodeData<IntegerType, RationalType>(),
		in_capacity(0) {}

	RationalType in_capacity;
};


template <
  typename IntegerType,
  typename RationalType
  >
struct FullNodeData : public BasicNodeData<IntegerType, RationalType> {
	FullNodeData() :
		BasicNodeData<IntegerType, RationalType>(),
		voltage(0) {}

	// electrical flow variables
	IntegerType voltage;
};



template <
  typename IntegerType,
  typename RationalType
  >
struct BasicArcData {
	BasicArcData() :
		capacity(0),
		cost(0),
		xlower(0) {}

	// input parameter
	IntegerType capacity;
	IntegerType cost;

	// primal variables - for x
	IntegerType xlower;
};


template <
  typename IntegerType,
  typename RationalType
  >
struct SSPVariantRBNArcData : public BasicArcData<IntegerType, RationalType> {
	SSPVariantRBNArcData() :
		BasicArcData<IntegerType, RationalType>(),
		flow_delta(0) {}

	RationalType flow_delta;
};


template <
  typename IntegerType,
  typename RationalType
  >
struct SSPVariantRCArcData : public BasicArcData<IntegerType, RationalType> {
	SSPVariantRCArcData() :
		BasicArcData<IntegerType, RationalType>(),
		xupper(0),
		potential(0),
		deficit(0),
		parent(0),
		depth(-1),
		visited(false),
#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
		transformed(false),
		lower_in_tree(false),
		upper_in_tree(false) {}
#else
		transformed(false) {}
#endif

	RationalType xupper;
	RationalType potential;
	IntegerType deficit;
	node parent;
	int depth;
	bool visited;
	bool transformed;
#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
	bool lower_in_tree;
	bool upper_in_tree;
#endif
};


template <
  typename IntegerType,
  typename RationalType
  >
struct InteriorPointMethodArcData : public BasicArcData<IntegerType, RationalType> {

	IntegerType resistance;
	IntegerType battery;
	IntegerType current;
	IntegerType f_0;

	RationalType flows;
	//dual varible
	RationalType s;

	vector<IntegerType> resistances_accross_cycle;
};


template <
  typename IntegerType,
  typename RationalType
  >
struct FullArcData : public BasicArcData<IntegerType, RationalType> {
	FullArcData() :
		BasicArcData<IntegerType, RationalType>(),
		f_0(0),
		flow(0),
		lower(0),
		infeasibility(0),
		slower(0),
		supper(0),
		sroof(0),
		potentialvw(0),
		croof(0),
		direction(0),
		voltage_vw(0),
		cur_src_vw(0),
		current(0),
		current_lower(0),
		current_roof(0),
		initial_current_lower(0),
		initial_current_roof(0),
		battery_lower(0),
		battery_upper(0),
		battery_roof(0),
		resistance_lower(0),
		resistance_upper(0),
		resistance_roof(0) {}

	// used for initialization
	RationalType f_0;
	RationalType flow;

	// primal variables - for x
	RationalType lower;
	IntegerType infeasibility;

	// dual variables
	IntegerType slower;
	IntegerType supper;
	IntegerType sroof;
	IntegerType potentialvw;

	// cost and direction of ahat
	IntegerType croof;
	int direction;

	// electrical flow variables
	IntegerType voltage_vw;
	IntegerType cur_src_vw;

	IntegerType current;

	IntegerType current_lower;
	IntegerType current_roof;

	IntegerType initial_current_lower;
	IntegerType initial_current_roof;

	RationalType battery_lower;
	RationalType battery_upper;
	RationalType battery_roof;

	IntegerType resistance_lower;
	IntegerType resistance_upper;
	IntegerType resistance_roof;
};

struct Original_Auxialiary_Transformed {

public:

	long double x1, x2, x3;
	long double s1, s2, s3;
	long double y1, y2, y3;
	long double cost1, cost2, cost3;
	long double demand1, demand2, demand3;
	long double capacity;
	long double g_prime1, g_prime2, g_prime3;
	long double x_hat1, x_hat2, x_hat3;
	long double s_hat1, s_hat2, s_hat3;
	long double z_prime1, z_prime2, z_prime3;
	long double r1_tilde_aux, r2_tilde_aux, r3_tilde_aux;
	long double r1_delta_wye, r2_delta_wye, r3_delta_wye;
	long double r1_aux, r2_aux, r3_aux;
	long double battery1_aux, battery2_aux, battery3_aux;
	long double chi_u;
	long double resistance_effective;
	unbounded_integer<long long> r1_tilde, r2_tilde, r3_tilde;
	unbounded_integer<long long> imb_vw;
	unbounded_integer<long long> imb_u;
	unbounded_integer<long long>  electrical_flow1_tilde, electrical_flow2_tilde,
	                  electrical_flow3_tilde;
	unbounded_integer<long long> current_non_tree_edge;
	unbounded_integer<long long> current_other_edge;
	unbounded_integer<long long> current_third_edge;
	unbounded_integer<long long> resistance_non_tree_edge;
	unbounded_integer<long long> tree_ind_vol_1;
	unbounded_integer<long long> tree_ind_vol_2;
	unbounded_integer<long long> current;
	unbounded_integer<long long> resistance_other_edge;
	unbounded_integer<long long> tree_induced_voltage_v;
	unbounded_integer<long long> tree_induced_voltage_w;
	unbounded_integer<long long> tree_induced_voltage_vw;
	unbounded_integer<long long> tree_induced_voltage_u;
	unbounded_integer<long long> f1_aux, f2_aux, f3_aux;
	unbounded_integer<long long> f1_delta_wye, f2_delta_wye, f3_delta_wye;
	unsigned int non_tree_edge;
	bool edge_reversed;
	long double electrical_flow1, electrical_flow2, electrical_flow3;
	bool arc_removed;
	bool is_invalid;
	bool cant_be_removed;
	Original_Auxialiary_Transformed() :
		x1(0), x2(0), x3(0), current_non_tree_edge(0),
		current_other_edge(0), current(0), edge_reversed(false),
		arc_removed(false), is_invalid(false), cant_be_removed(false) {}
};


// Network
template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData = FullNodeData<IType, RType>,
  typename ArcData = FullArcData<IType, RType>
  >
class Network {
public:
	typedef Graph GraphType;
	typedef IType IntegerType;
	typedef RType RationalType;

	Graph &G;
	vector<NodeData> nodedata;
	vector<ArcData> arcdata;

	Network(Graph& G) : G( G ), nodedata( G.no_of_vertices + 1 ), arcdata( G.no_of_edges + 1 ) {};
	Network(Graph& G, const string &filename);
	vector<double> get_distribution(vector<IntegerType>&) const;
	auto calculate_primal_objective_value() -> RationalType;

	void pull_flow(
	  arc a,
	  const RationalType& delta,
	  vector<RationalType>& b
	) ;

	void write_dimacs(ostream &os);

private:
	void read_graphml(const string &filename);
	void read_dimacs(const string &filename);
};


template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData,
  typename ArcData
  >
Network<Graph, IType, RType, NodeData, ArcData>::Network(Graph &G, const string &filename) :
	G(G),
	nodedata(G.no_of_vertices + 1),
	arcdata(G.no_of_edges + 1) {

	if (boost::ends_with(filename, ".graphml")) {
		// if the file name ends with .graphml, read it as graphml
		read_graphml(filename);
	} else {
		// otherwise assume the graph is in dimacs format
		read_dimacs(filename);
	}

#ifndef NDEBUG
	// sanity checks
	auto sum_demands = static_cast<decltype(nodedata.front().demand)>(0);
	for (const auto &data : nodedata)
		sum_demands += data.demand;
	assert(sum_demands == 0 && "the sum of all demands must be zero.");
#endif
};


template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData,
  typename ArcData
  >
void Network<Graph, IType, RType, NodeData, ArcData>::read_graphml(const string &filename) {
	ptree tree;
	read_xml(filename, tree);

	string costs;
	string capacities;

	// find out keys representing costs/capacities
	for (const auto &element : tree.get_child("graphml")) {
		if (element.first != "key")
			continue;
		auto attribute_name = element.second.get<string>(ptree::path_type("<xmlattr>/attr.name", '/'));
		if (attribute_name == "costs" || attribute_name == "cost")
			costs = element.second.get<string>("<xmlattr>.id");
		if (attribute_name == "capacities" || attribute_name == "capacity")
			capacities = element.second.get<string>("<xmlattr>.id");
	}

	if (costs.empty()) {
		std::cerr << "cost attribute not found!" << std::endl;
		exit(1);
	}

	if (capacities.empty()) {
		std::cerr << "capacity attribute not found!" << std::endl;
		exit(1);
	}

	const auto &graph = tree.get_child("graphml.graph");

	// set costs/capacities/demands
	int edgeid = 0;
	for (const auto &element : graph) {
		if (element.first == "node") {
			auto nodeid = element.second.get<int>("<xmlattr>.id");
			nodedata[nodeid].demand = element.second.get<int>("data");
		}
		if (element.first == "edge") {
			++edgeid;
			for (const auto &child : element.second) {
				if (child.first != "data")
					continue;
				if (child.second.get<string>("<xmlattr>.key") == costs)
					arcdata[edgeid].cost = child.second.get_value<int>();
				if (child.second.get<string>("<xmlattr>.key") == capacities)
					arcdata[edgeid].capacity = child.second.get_value<int>();
			}
		}
	}
}


template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData,
  typename ArcData
  >
void Network<Graph, IType, RType, NodeData, ArcData>::read_dimacs(const string &filename) {
	ifstream dimacs_file(filename);
	string line;

	// parse node descriptors
	while (getline(dimacs_file, line)) {
		if (boost::starts_with(line, "c") || boost::starts_with(line, "p"))
			continue;
		if (boost::starts_with(line, "a"))
			break;
		assert(boost::starts_with(line, "n") && "here should be a node descriptor");

		auto node_tokens = vector<string>();
		boost::split(node_tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);

		assert(stoul(node_tokens[1]) > 0 && stoul(node_tokens[1]) <= G.no_of_vertices && "node ids must be in the range [1, n]");

		nodedata[stoul(node_tokens[1])].demand = -stoi(node_tokens[2]);
		//cout << "d = " << to_double(nodedata[stoul(node_tokens[1])].demand) << endl;
	}
  
	int edgeid = 0;

	// parse arc descriptors
	do {
		if (boost::starts_with(line, "c") || line.empty())
			continue;
		assert(boost::starts_with(line, "a") && "here should be an arc descriptor");
		++edgeid;

		auto arc_tokens = vector<string>();
		boost::split(arc_tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);

		if (stoi(arc_tokens[3]) != 0) {
			cerr << "currently only networks with minimum flow of 0 are supported." << endl;
			exit(1);
		}

		arcdata[edgeid].capacity = stoi(arc_tokens[4]);
		arcdata[edgeid].cost = stoi(arc_tokens[5]);
	} while (getline(dimacs_file, line));
}


template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData,
  typename ArcData
  >
void Network<Graph, IType, RType, NodeData, ArcData>::write_dimacs(ostream &os) {
	// problem line
	os << "p min " << G.no_of_vertices << " " << G.no_of_edges << '\n';

	// node descriptors
	for (auto n = 1u; n <= G.no_of_vertices; ++n) {
		const auto &data = nodedata[n];
		if (data.demand != 0)
			os << "n " << n << " " << -data.demand << '\n';
	}

	// arc descriptors
	for (auto a = 1u; a <= G.no_of_edges; ++a) {
		const auto &data = arcdata[a];
		os << "a " << G.tails[a] << " " << G.heads[a] << " 0 " << data.capacity << " " << data.cost << '\n';
	}
	os.flush();
}


template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData,
  typename ArcData
  >
void Network<Graph, IType, RType, NodeData, ArcData>::pull_flow (
  arc a,
  const RationalType& delta,
  vector<RationalType>& b
) {

	if ( a > 0 ) {

		const node v = G.tails[a];
		const node w = G.heads[a];
		arcdata[a].flow -= delta;
		b[v] -= delta;
		b[w] += delta;

	} else {

		const node w = G.tails[-a];
		const node v = G.heads[-a];

		arcdata[-a].flow += delta;
		b[v] -= delta;
		b[w] += delta;

	}
}

template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData,
  typename ArcData
  >
auto Network<Graph, IType, RType, NodeData, ArcData>::calculate_primal_objective_value() -> RType {
	return accumulate(begin(arcdata) + 1, end(arcdata), static_cast<RType>(0), [](RType value, const ArcData & arcdata) { return value + arcdata.xlower * arcdata.cost; });
}

// Spanning Tree
template <
  typename Graph,
  typename IntegerType,
  typename RationalType
  >
class SpanningTree {
public:
	const Graph &G;

	//! A constructor
	/**
	 * A constructor
	 */
	SpanningTree( const Graph& G );

	std::vector<int> depth;
	node root;
	std::vector<unsigned int> size_sub_tree_rooted_at;
	std::vector<std::vector<NodeInformation>> node_tree_corrospondance;
	std::vector<std::vector<int> > tree_incident_edges;
	std::vector<std::vector<node> > children;
	std::vector<unsigned int> discovery;
	std::vector<unsigned int> finished;
	std::vector<arc> dfs_euler_tour;
	std::vector<IntegerType> d_drop;
	std::vector<IntegerType> d_ext;
	std::vector<IntegerType> sum_of_the_resistances_from_root_to;
	std::vector<RationalType> sum_of_the_unrounded_resistances_from_root_to;
	std::vector<node> LCA;
	node tree_decompose(node, unsigned int&);
	std::vector<arc> arc_to_parent;
	std::vector<node> parent;
	std::vector<unsigned int> non_tree_edge_common_tree_index;
	std::vector<node> find_adjacent_nodes (node);
	void find_LCA_and_above_d_flags(node d, unsigned int);
	long double calculate_tree_condition_number(Network<Graph, IntegerType, RationalType>& N);
	long double calculate_tree_condition_number_wrt_unrounded_resistances();

	/** \brief finds the vertex separator
	 *
	 */
	node find_vertex_separator(node, unsigned int, unsigned int&);

	/**
	 * A member function taking two arguments and returning void
	 * @param a node arguement
	 * @param an unsigned integer arguement
	 */
	void init(node, unsigned int&);


	/**
	 * A constant member taking two arguments and returning template
	 * @param a node arguement
	 * @param an unsigned int arguement
	 */
	IntegerType query(node, unsigned int) const;

	/**
	 * A constant member function taking one argument and returning template
	 * @param an unsigned int argument
	 */
	IntegerType query(unsigned int) const;

	/**
	 * A member function taking two arguements and returning void
	 * @param a node arguement
	 * @param an unsigned int arguement
	 */
	void update( node, const IntegerType& );

	/**
	 * A member function taking three arguements and returning void
	 * @param a node arguement
	 * @param a template arguement
	 * @param an unsigned int arguement
	 */
	void update(node, const IntegerType&, unsigned int);

	void update_with_DS(node,  IntegerType&, unsigned int, IntegerType);
	/**
	 * A member function taking two arguements and returning void
	 * @param a template arguement
	 * @param an unsigned int arguement
	 */
	void update( const IntegerType&, unsigned int );

	/** \brief finds the least common ancestor
	 *
	 * @param a the arc
	 *
	 * @return Returns the least common ancestor of the two end points of the arc a
	 *
	 */
	node find_LCA(arc a);

	/** \brief finds the least common ancestor
	 *
	 * @param v node
	 * @param w node
	 *
	 * @return Returns the least common ancestor of v and w
	 */
	node find_LCA(node v, node w );

	/** \brief store LCA for each edge in the graph
	 *
	 */
	void fill_LCA_vector(vector<RationalType>&);
	void fill_LCA_vector();

	/** \brief stores the edges incident to each node in the spanning tree
	 *
	 */
	void get_tree_incident_edges();

	void print_children();

	void print_sizes();

	/** \brief Stores the sum of the resistances from the node to the root
	 *
	 */
	void update_sum_to_the_root(node, Network<Graph, IntegerType, RationalType>& N);
	void update_sum_to_the_root_wrt_unrounded_resistances(node, Network<Graph, IntegerType, RationalType>& N);

	void get_initial_state_of_the_tree_decomposition();

	void get_non_tree_edge_common_tree_index();

	/** \brief labels
	 */
	void label_the_tree_rooted_at(node, node, unsigned int, unsigned int &);

	/**
	 * A member function taking one arguement and returning template
	 * @param an arc arguement
	 */
	IntegerType compute_resistance_accross(arc, Network<Graph, IntegerType, RationalType>& N );

	RationalType compute_unrounded_resistance_accross(arc, Network<Graph, IntegerType, RationalType>& N );



	IntegerType voltage_drop( arc a ) const  {
#ifdef VERBOSE
		cout << "voltage_drop function " << a << " = (" << G.tail(a) << " , " << G.head(a) << ")" << endl;
#endif
		return query(G.tail(a)) - query(G.head(a));
	}

	/** \brief calculates the voltage drop
	 *
	 * @param v
	 * @param w
	 *
	 * @return Returns the voltage drop beween nodes v and w
	 *
	 */
	IntegerType voltage_drop(node v , node w) const {
		return query(v) - query(w);
	}

	struct Stack_Data {
		node d;
		unsigned int tree_index;
		node s;
	};


	//std::vector<bool> root_is_split_node;
	//std::vector<node> root_of_split_node;


	std::vector<TreeDecompositionInfo> tree_decomposition;
	stack<Stack_Data> S;

	/** \brief Runs the init function for the part of the tree below 'd'
	 *
	 * @param Stack_Data
	 * @param node
	 * @param node
	 * @param int
	 */
	void init_below( node, node, unsigned int&);

	void clear();

	/** \brief Clears the vectors d_ext and d_drop
	 *
	 */
	void clear_d();

};


template <
  typename Graph,
  typename IntegerType,
  typename RationalType
  >
SpanningTree<Graph, IntegerType, RationalType>::SpanningTree( const Graph& g ) : G(g),
	depth(G.no_of_vertices + 1, -1),
	root(0),
	size_sub_tree_rooted_at(G.no_of_vertices + 1, 0),
	node_tree_corrospondance(G.no_of_vertices + 1),
	tree_incident_edges(G.no_of_vertices + 1),
	children(G.no_of_vertices + 1),
	discovery(G.no_of_vertices + 1, 0),
	finished(G.no_of_vertices + 1, 0),
	d_drop(3 * G.no_of_vertices + 1, 0),
	d_ext(3 * G.no_of_vertices + 1, 0),
	sum_of_the_resistances_from_root_to(G.no_of_vertices + 1, 0),
	sum_of_the_unrounded_resistances_from_root_to(G.no_of_vertices + 1, 0),
	LCA(G.no_of_edges + 1, 0),
	arc_to_parent(G.no_of_vertices + 1, 0),
	parent(G.no_of_vertices + 1, 0),
//root_is_split_node(G.no_of_vertices + 1,0),
//root_of_split_node(G.no_of_vertices + 1,0),
	tree_decomposition(3 * G.no_of_vertices + 1)
{
	for (unsigned int i = 1; i <= G.no_of_vertices; i++)
	{
		node_tree_corrospondance[i].push_back( NodeInformation() );
	}
	//    for(unsigned int i= 0; i<=G.no_of_edges; i++)
	//    {
	//      G.rounding_scheme_tree_edges[i] = i;
	//    }
}


/** \brief clear
 *
 * clears the information in the object
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::clear() {

	assert( node_tree_corrospondance.size() == G.no_of_vertices + 1 );
	vector<NodeInformation>().swap( node_tree_corrospondance[0] );
	vector<node>( G.no_of_edges + 1, 0).swap( LCA );
	for (unsigned int i = 1; i <= G.no_of_vertices; i++)
	{

		vector<arc>().swap( tree_incident_edges[i] );
		vector<NodeInformation>().swap( node_tree_corrospondance[i] );
		node_tree_corrospondance[i].push_back( NodeInformation() );
		root = 0;

		depth[i] = -1;
		children[i].clear();
		size_sub_tree_rooted_at[i] = 0;
		discovery[i] = 0;
		finished[i] = 0;
		sum_of_the_resistances_from_root_to[i] = 0;
		arc_to_parent[i] = 0;
		parent[i] = 0;
		//root_is_split_node[i] = 0;
		//root_of_split_node[i] = 0;
	}

	clear_d();
}

/** \brief
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::get_non_tree_edge_common_tree_index()
{
	non_tree_edge_common_tree_index.clear();
	non_tree_edge_common_tree_index.reserve(G.non_tree_edges.size() );
	for (auto a : G.non_tree_edges)
	{
		node u = G.tail(a);
		node v = G.head(a);

		unsigned int index = 0;
		for ( auto iter1 = node_tree_corrospondance[u].begin(), iter2 = node_tree_corrospondance[v].begin(),
		      end1 = node_tree_corrospondance[u].end(), end2 = node_tree_corrospondance[v].end();
		      iter1 != end1 && iter2 != end2; ++iter1, ++iter2 )
		{
			if ( iter1->tree_index == iter2->tree_index && iter1->above_d == iter2->above_d && iter1->lca == iter2->lca ) {
				++index;
			} else {
				non_tree_edge_common_tree_index.push_back( index );
				break;
			}
		}
	}

}



/** \brief Clear d_ext and d_drop
 *
 * Clears the vectors d_ext and d_drop
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::clear_d() {

	assert(d_ext.size() == d_drop.size());

	for (unsigned int i = 0; i < d_ext.size(); i++) {
		d_ext[i] = 0;
		d_drop[i] = 0;
	}

}

/** \brief create a minimum spanning tree
 *
 * Finds the Minimum Spanning Tree of the given Graph
 *
 * @param SpanningTree ST
 */
template< typename SpanningTree, typename NumberType >
void
create_MST_wrt_given_vector(
  vector<NumberType>&x,
  SpanningTree& ST)
{

//  tree_edges.clear();
//  non_tree_edges.clear();

	std::vector<unsigned int> permutation( x.size() );
	for (unsigned int i = 0; i < permutation.size(); i++) {
		permutation[i] = i;
	}

	std::sort(permutation.begin() + 1, permutation.end(), [x] (unsigned int a, unsigned int b) -> bool { return x[a] < x[b]; } );


	make_tree(ST, permutation);


}


// Graph
template<typename IntegerType, typename RationalType>
class Graph {

	Random rg;

public:

	void swap( RationalType&, RationalType& );
	void swap (int &, int &);
	unsigned int no_of_vertices;
	unsigned int no_of_edges;
	unsigned int count_tree_edges;
	unsigned int no_of_removed_edges;
	//unsigned int count;
	unsigned int i, j;
	vector<node> heads;
	vector<node> tails;
	//std::vector<node> parents;
	vector<node> parent;
	vector<vector<int> > incident_edges;
	//vector<bool> removed_arc;
	//vector<int> visited;
	vector<int> rank;
	vector<arc> arcs_to_root;
	//node root;
	vector<arc> tree_edges;
	vector<int> non_tree_edges;
//   std::vector<int> probability_distribution;
//   std::vector<RationalType> demands;
//   std::vector<Original_Auxialiary_Transformed> original_auxiliary_transformed_corrospondence;
//   std::vector<IntegerType> resistances;
//   std::vector<RationalType>resistances_aux;
//   std::vector<RationalType> unrounded_resistances;
	std::vector<RationalType> max_val;
//   std::vector<RationalType> rounded_resistances_auxiallary;
//   std::vector<RationalType> weights;
	//std::vector<RationalType> voltages;
//   std::vector<RationalType> capacities;
//   std::vector<RationalType> costs;
//   std::vector<RationalType> flows;
	//std::vector<IntegerType> f;
//   std::vector<RationalType> unrounded_currents;
//   std::vector<IntegerType> currents;
//   std::vector<RationalType> batteries;
	//std::vector<RationalType> batteries_unrounded;
	//std::vector<RationalType> g_tilde;
	//std::vector<IntegerType> r_tilde;
//  std::vector<IntegerType> f_0;
	//std::vector<RationalType> shortest_distances;
//   std::vector<IntegerType> tree_induced_voltages;
//   std::vector<RationalType> probability;
	//std::vector<RationalType> tree_induced_voltages_with_batteries;
	vector<IntegerType> resistances_accross_cycles;
	vector<RationalType> resistances_accross_cycles_aux;
	vector<RationalType> unrounded_resistances_accross_cycles;
	//std::vector<int> already_pushed;
//   std::vector<unsigned int> rounding_scheme_tree_edges;
//   std::vector<vector<arc>> rounding_scheme_incident_edges;
//   std::vector<unsigned int> cluster;
//   std::vector< vector<node> > nodes_in_cluster;
//   vector<ArcConstruction> arc_map;
	int create_edge(node , node);
	arc new_edge( node, node );
	arc new_edge( arc, node, node );
	void remove_arc(arc);

	Graph(const unsigned int n, const unsigned int m = 0) :

		heads(m + 1, 0),
		tails(m + 1, 0),
		//parents( n+1, 0 ),
		//parent( n+1, 0 ),
		incident_edges( n + 1 ),
		//removed_arc(m+1, false),
		//visited( n+1,0 ),
		rank(n + 1, 0),
		arcs_to_root( n + 1, 0 ),
		//demands(n+1,RationalType(0) ),
		//original_auxiliary_transformed_corrospondence(m +1),
		//resistances(m+1,IntegerType(0) ),
		//resistances_aux(m+1,RationalType(0) ),
		//unrounded_resistances(m+1, RationalType(0) ),
		max_val(m + 1, RationalType(0) )
		//rounded_resistances_auxiallary(m+1, RationalType(0) ),
		//weights(m+1, RationalType(0) ),
		//voltages(m+1,RationalType(0)), // Voltages refer to the battaries that are attached in the circuit
		//capacities(m+1,RationalType(0)),
		//costs(m+1,RationalType(0)),
		//flows(m+1,RationalType(0)),
		//f(m+1, IntegerType(0)),
		//unrounded_currents(m+1, RationalType(0)),
		//currents(m+1,IntegerType(0)),
		//batteries(m+1,RationalType(0) ),
		//batteries_unrounded(m+1, RationalType(0)),
		//g_tilde(m+1, RationalType(0)),
		//r_tilde(m+1, IntegerType(0) ),
		//f_0(m + 1, IntegerType(0) ),
		//tree_induced_voltages( n+1, IntegerType(0) ),
		//probability(m , RationalType(1)),
		//rounding_scheme_tree_edges(m + 1),
		//rounding_scheme_incident_edges(m + 1),
		//cluster(n+1),
		//nodes_in_cluster(n+1),
		//arc_map(m/3+1,ArcConstruction())
	{
#ifndef NDEBUG
		cout << "The seed which generated the random numbers for the graph: " << rg.z << endl;
#endif
		no_of_vertices = n;
		no_of_removed_edges = 0;
		no_of_edges = 0;
		count_tree_edges = 0; // The number of tree edges is initially set to 0
		resistances_accross_cycles.reserve( m + 1 );
		resistances_accross_cycles_aux.reserve( m + 1 );
		unrounded_resistances_accross_cycles.reserve(m + 1);
		//visited.reserve(n+1);
		for (unsigned int i = 0; i <= n; i++)
		{
			//parents[i] = 0;
			parent[i] = i;
		}

//    demands.resize(n+1);

		//     RationalType demand_n(0);
		//
		//     for( unsigned int i=1; i<n; i++){
		//       double d = rg.rng(); //random number to generate random demands
		//       demands[i] = ceil(d*30);
		//       demand_n += demands[i];
		//     }
		//     demands[n] = -demand_n;

	}

	Graph(const string &filename) :
		no_of_vertices(0),
		no_of_edges(0),
		count_tree_edges(0),
		no_of_removed_edges(0),
		i(), //?
		j(), //?
		heads(),
		tails(),
		parent(),
		incident_edges(),
		rank(),
		arcs_to_root(),
		tree_edges(),
		non_tree_edges(),
		resistances_accross_cycles(),
		resistances_accross_cycles_aux(),
		unrounded_resistances_accross_cycles() {

		if (boost::ends_with(filename, ".graphml")) {
			// if the file name ends with .graphml, read it as graphml
			read_graphml(filename);
		} else {
			// otherwise assume the graph is in dimacs format
			read_dimacs(filename);
		}
	}

	void print();
	void print_lp();
	void print_node_data();
	void print_arc_data();
	void create_graph();
	void read_graph(string filepath, int m);
	node find(node);
	bool unite( arc );


	void form_petal_decomposition_spanning_tree(SpanningTree<Graph, IntegerType, RationalType>& );
	RationalType size_del(vector<node>&);
	vector<arc> del(vector<node>&, vector<node>&);
	RationalType size_E(vector<node>&);
	RationalType vol_x(vector<node>&, vector<node>&, node, RationalType);
	RationalType size_E_plus(vector<node>&);
	RationalType dijkstra(node, node, vector<node>&, vector<node>&);
	RationalType dijkstra(node, node);
	vector<RationalType> dijkstra(node);
	vector<RationalType> dijkstra(node, vector<node>&);
	void Ball(node, node, RationalType, vector<node>&, vector<node>&);
	void Ball_x(node, RationalType, vector<node>& , vector<node>&);
	vector<node> subtract_vectors(vector<node>&, vector<node>&);
	unsigned int subtract_vectors_size(vector<node>&, vector<node>&);
	vector<node> get_W_r(node, node, RationalType, vector<node>& , vector<node>&);
	RationalType cost(vector<arc>);
	RationalType rad(node, vector<node>& );
	void add_path(node, unsigned int);
	void create_petal(vector<node>& , vector<node>& , node , node, RationalType , RationalType, vector<node>& , node& );
	unsigned int create_remaining_petals(
	  unsigned int ,
	  node ,
	  node ,
	  RationalType,
	  vector<node>&,
	  vector<vector<node>>& ,
	  vector<vector<node>>&,
	  vector<node>&,
	  vector<node>&,
	  vector<node>&
	);
	void create_stigma(
	  vector<vector<node>>&,
	  vector<vector<node>>&,
	  node,
	  vector<node>&,
	  vector<node>&
	);
	vector<node> get_shortest_path(node, node, vector<node>&);
	vector<unsigned int> create_BFS_tree(vector<node>& , node );
	unsigned int petal_decomposition(
	  vector<node>& ,
	  node ,
	  unsigned int,
	  vector<vector<node>>& ,
	  vector<node>& ,
	  vector<node>&,
	  vector<node>&
	);
	vector<unsigned int> hierarchial_petal_decomposition(
	  SpanningTree<Graph,
	  IntegerType,
	  RationalType>&ST,
	  vector<node>&,
	  node,
	  node
	);

	vector<RationalType> dijkstra_clusters(
	  node,
	  unsigned int,
	  unsigned int,
	  vector<vector<arc>>&,
	  vector<unsigned int>&,
	  vector<unsigned int>&
	);
	vector< vector<arc> > contract_clusters(
	  unsigned int,
	  vector< vector<arc>> &,
	  vector<unsigned int>& ,
	  vector<unsigned int>&
	);
	void form_spanning_tree_alon(SpanningTree<Graph, IntegerType, RationalType>& ST, RationalType);
	vector<RationalType> scale_resistances();
	int get_clusters(
	  unsigned int,
	  vector<vector<arc> > &,
	  vector<unsigned int>& ,
	  vector<unsigned int>&,
	  unsigned int,
	  RationalType
	);
	void get_shortest_path_spanning_tree(
	  unsigned int ,
	  node,
	  vector<unsigned int>& ,
	  vector<vector<arc>>& ,
	  unsigned int,
	  vector<unsigned int>&,
	  vector<unsigned int>&
	);
	void clear_cluster_details();
	/** \brief flips the direction of the arc
	 *
	 */
	void flip_arc(arc);

	/** \brief
	 *
	 */
	vector<arc> get_path(const node, const node, SpanningTree<Graph, IntegerType, RationalType>& ST);
	vector<arc> get_cycle(const arc, SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>&);
	// int create_node();
	/** \brief gets the edge between given two nodes
	 */
	int get_edge(const node , const node );
	node get_the_other_end_of_arc(node , arc);
	void simple_solver();
	//T calculate_gap(std::vector<T> &);
	//vector<double> get_distribution() const;
	void the_sanity_check();

	void change_spanning_tree();

	void get_random_resistances();

	/** \brief creates a low stretch spanning tree for the given Graph
	 *
	 * @param ST The spanning tree to be created
	 *
	 */
	void create_low_stretch_tree(
	  SpanningTree<Graph,
	  IntegerType,
	  RationalType>& ST,
	  Network<Graph<IntegerType, RationalType>,
	  IntegerType,
	  RationalType>& N
	);
	void  create_low_stretch_tree_wrt_unrounded_resistances(
	  vector<RationalType>&x,
	  SpanningTree<Graph,
	  IntegerType,
	  RationalType>& ST,
	  Network<Graph,
	  IntegerType,
	  RationalType>& N
	);

	node get_l_star( vector<vector<int>>&, unsigned int , unsigned int , RationalType);

	/** \brief creates a Mininum Spanning Tree
	 *
	 */
	void create_MST(SpanningTree<Graph, IntegerType, RationalType>& ST);

	void create_MST_wrt_unrounded_resistances(
	  vector<RationalType>&,
	  SpanningTree<Graph,
	  IntegerType,
	  RationalType>& ST
	);
	void create_low_stretch_tree_wrt_unrounded_resistances_on_original_graph(
	  SpanningTree<Graph,
	  IntegerType,
	  RationalType>& ST,
	  Network<Graph<IntegerType, RationalType>,
	  IntegerType,
	  RationalType>& N
	);
	void create_MST_wrt_unrounded_resistances_on_original_graph(
	  SpanningTree<Graph,
	  IntegerType,
	  RationalType>& ST,
	  Network<Graph, IntegerType, RationalType>& N
	);
	void create_SPT(SpanningTree<Graph, IntegerType, RationalType>& ST);
	void create_BFS_tree_with_random_root(
	  SpanningTree< Graph,
	  IntegerType,
	  RationalType>& ST
	);
	/** \brief creates a tree corresponding to the one in rounding scheme algorithm
	 *
	 */
	void create_rounding_scheme_tree(SpanningTree<Graph, IntegerType, RationalType>&);
	void make_tree(SpanningTree<Graph, IntegerType, RationalType>&, std::vector<unsigned int> &);
	void euler_tour_and_get_subtree_sizes(SpanningTree<Graph, IntegerType, RationalType>& ST);
	void get_depth_and_children(SpanningTree<Graph, IntegerType, RationalType>& ST);

	/**
	 *
	 */
	void get_layers_in_spanning_tree(
	  node ,
	  int,
	  vector<vector<arc>>&,
	  unsigned int,
	  vector<unsigned int>&,
	  vector<unsigned int>&,
	  unsigned int,
	  RationalType
	);



	void print_arc_data(
	  const vector<RationalType>& array
	) {
		for (unsigned int i = 1; i <= no_of_edges; i++) {
			cout << i << ": " << array[i] << "; ";
		}
		cout << endl;
	}

	void print_node_data(
	  const vector<RationalType>& array
	) {
		for (unsigned int i = 1; i <= no_of_vertices; i++) {
			cout << i << ": " << array[i] << "; ";
		}
		cout << endl;
	}

	void initialize_DS()
	{
//    original_auxiliary_transformed_corrospondence.reserve( no_of_edges+1 );
	}


//   template< typename IT> void push_flow( arc a, const IT& delta, vector<RationalType>& b ) {
//     if( a > 0 ) {
//       const node v = tails[a];
//       const node w = heads[a];
//       flows[a] += delta;
//       b[v] -= delta;
//       b[w] += delta;
//     } else {
//       const node w = tails[-a];
//       const node v = heads[-a];
//       flows[-a] -= delta;
//       b[v] -= delta;
//       b[w] += delta;
//     }
//   }
//   template< typename IT> void push_flow( arc a, const IT& delta ) {
//     push_flow( a, delta, demands );
//   }

	node tail( arc a ) const {
		return a > 0 ? tails[a] : heads[-a];
	}
	node head( arc a ) const {
		return a > 0 ? heads[a] : tails[-a];
	}
//   IntegerType voltage_drop( arc a ) const {
//     return tree_induced_voltages[tail(a)] - tree_induced_voltages[head(a)];
//   }

//   IntegerType voltage_drop(node v , node w) const {
//     return tree_induced_voltages[v] - tree_induced_voltages[w];
//   }

private:
	void read_graphml(const string &filename);
	void read_dimacs(const string &filename);

	void initialize_internal_data_structures(int n, int m);
};


template <
  typename IntegerType,
  typename RationalType
  >
void Graph<IntegerType, RationalType>::read_graphml(const string &filename) {
	ptree tree;
	read_xml(filename, tree);

	const ptree &graphml = tree.get_child("graphml", empty_ptree());
	const ptree &graph = graphml.get_child("graph", empty_ptree());

	int n = 0;
	int m = 0;
	int nodeid;

	// count nodes and edges
	BOOST_FOREACH(const ptree::value_type & nore, graph) {
		const ptree & nore_attrs = nore.second.get_child("<xmlattr>", empty_ptree());
		BOOST_FOREACH(const ptree::value_type & nore_attr, nore_attrs) {
			if (strncmp(nore_attr.first.data(), "id", 2) == 0) {
				nodeid = stoi(nore_attr.second.data());
				n = max(n, nodeid);
			}
			if (strncmp(nore_attr.first.data(), "source", 6) == 0)
				++m;
		}
	}

	initialize_internal_data_structures(n, m);

	// parse edges
	BOOST_FOREACH(const ptree::value_type & nore, graph) {

		const ptree & nore_attrs = nore.second.get_child("<xmlattr>", empty_ptree());
		bool edge = false;
		int source = 0;
		int target = 0;
		BOOST_FOREACH(const ptree::value_type & nore_attr, nore_attrs) {
			if (strncmp(nore_attr.first.data(), "id", 2) != 0) {
				if (strncmp(nore_attr.first.data(), "source", 6) == 0) {
					edge = true;
					source = stoi(nore_attr.second.data());
				}
				if (strncmp(nore_attr.first.data(), "target", 6) == 0) {
					assert(edge);
					target = stoi(nore_attr.second.data());
				}
			}
		}

		if (edge) {
			++no_of_edges;
			new_edge(no_of_edges, source, target);
		}
	}
}


template <
  typename IntegerType,
  typename RationalType
  >
void Graph<IntegerType, RationalType>::read_dimacs(const string &filename) {
	ifstream dimacs_file(filename);
	string line;

	// skip to problem description line
	while (getline(dimacs_file, line))
		if (boost::starts_with(line, "p"))
			break;

	auto problem_tokens = vector<string>();
	boost::split(problem_tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);
	assert(problem_tokens[0] == "p");
	assert(problem_tokens[1] == "min");
	auto n = stoul(problem_tokens[2]);
	auto m = stoul(problem_tokens[3]);

	initialize_internal_data_structures(n, m);

	// skip to arc descriptors
	while (getline(dimacs_file, line))
		if (boost::starts_with(line, "a"))
			break;

	// parse arc descriptors
	do {
		if (boost::starts_with(line, "c") || line.empty())
			continue;
		assert(boost::starts_with(line, "a") && "here should be an arc descriptor");

		auto arc_tokens = vector<string>();
		boost::split(arc_tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);

		assert(arc_tokens.size() == 6);

		++no_of_edges;
		new_edge(no_of_edges, stoul(arc_tokens[1]), stoul(arc_tokens[2]));
	} while (getline(dimacs_file, line));
}


template <
  typename IntegerType,
  typename RationalType
  >
void Graph<IntegerType, RationalType>::initialize_internal_data_structures(int n, int m) {
	// workaround for ids starting at 1
	no_of_vertices = n;

	heads.resize(m + 1, 0);
	tails.resize(m + 1, 0);

	incident_edges.resize(n + 1);
	arcs_to_root.resize(n + 1, 0);
	rank.resize(n + 1, 0);
	resistances_accross_cycles.reserve(m + 1);
	resistances_accross_cycles_aux.reserve(m + 1);
	unrounded_resistances_accross_cycles.reserve(m + 1);

	parent.resize(n + 1, 0);
	for (auto i = 0; i <= n; ++i)
		parent[i] = i;
}


/*
template<typename Graph, typename IntegerType, typename RationalType>
RationalType
SpanningTree<Graph, IntegerType, RationalType>::compute_unrounded_resistance_accross(arc a)
{
  const node v = G.head(a);
  const node w = G.tail(a);
  const node lca =  LCA[a];
  RationalType resistances_accross_cycle =
    sum_of_the_unrounded_resistances_from_root_to[v]
    + sum_of_the_unrounded_resistances_from_root_to[w]
    - RationalType(2) * sum_of_the_unrounded_resistances_from_root_to[lca]
    + G.unrounded_resistances[a];

  return resistances_accross_cycle;
}*/


template<typename Graph, typename IntegerType, typename RationalType>
RationalType
SpanningTree<Graph, IntegerType, RationalType>::compute_unrounded_resistance_accross(arc a, Network<Graph, IntegerType, RationalType>& N )
{
	const node v = G.head(a);
	const node w = G.tail(a);
	const node lca =  LCA[a];
	IntegerType xupper = N.arcdata[a].capacity - N.arcdata[a].xlower;

	RationalType resistance_a(0);
//  if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper)
	if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
	    N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper)
	{
		resistance_a =  to_double(N.arcdata[abs(a)].resistance_roof);
	}

	resistance_a =  RationalType(1) / to_double(N.arcdata[a].xlower * N.arcdata[a].xlower) +
	                RationalType(1) / to_double(xupper * xupper);
	RationalType resistances_accross_cycle =
	  sum_of_the_unrounded_resistances_from_root_to[v]
	  + sum_of_the_unrounded_resistances_from_root_to[w]
	  - RationalType(2) * sum_of_the_unrounded_resistances_from_root_to[lca]
	  + resistance_a;

	return resistances_accross_cycle;
}


/** \brief compute the resistances across cycle
 *
 * @param[in]
 * @return Returns the sum of resistances across the cycle
 */
template<typename Graph, typename IntegerType, typename RationalType>
IntegerType
SpanningTree<Graph, IntegerType, RationalType>::compute_resistance_accross(arc a, Network<Graph, IntegerType, RationalType>& N )
{

	const node v = G.head(a);
	const node w = G.tail(a);
	const node lca =  LCA[a];

#ifdef VERBOSE
	cout << "v = " << v << " , " << "w = " << w << ", " << "lca = " << lca << endl;

	cout << "resistance accross cycle = "
	     << to_double(sum_of_the_resistances_from_root_to[v])
	     << " + " << to_double(sum_of_the_resistances_from_root_to[w])
	     << " - " << 2 << " * " << to_double(sum_of_the_resistances_from_root_to[lca])
	     << " + " << to_double(G.resistances[abs(a)]) << endl;
#endif


	IntegerType resistances_accross_cycle(0);

//  if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper)
	if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
	    N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper)
	{
		resistances_accross_cycle = sum_of_the_resistances_from_root_to[v]
		                            + sum_of_the_resistances_from_root_to[w]
		                            - IntegerType(2) * sum_of_the_resistances_from_root_to[lca]
		                            //+ N.arcdata[abs(a)].resistances;
		                            + N.arcdata[abs(a)].resistance_upper + N.arcdata[abs(a)].resistance_lower;
	}
	else
	{
		resistances_accross_cycle = sum_of_the_resistances_from_root_to[v]
		                            + sum_of_the_resistances_from_root_to[w]
		                            - IntegerType(2) * sum_of_the_resistances_from_root_to[lca]
		                            //+ N.arcdata[abs(a)].resistances;
		                            + N.arcdata[abs(a)].resistance_upper + N.arcdata[abs(a)].resistance_lower;
	}


	return resistances_accross_cycle;
}


/** \brief calculates the tree condition number
 *
 * Returns the tree condition number
 *
 */
template <
  typename Graph,
  typename IntegerType,
  typename RationalType
  >
long double SpanningTree< Graph, IntegerType, RationalType>::calculate_tree_condition_number(
  Network<Graph,
  IntegerType,
  RationalType>& N
) {
#ifdef VERBOSE
	cout << "calculate_tree_condition_number() " << endl << endl;
#endif

	long double tree_condition_number = 0.0;
	for (unsigned int i = 0; i < G.non_tree_edges.size(); i++)
	{

		arc edge = G.non_tree_edges[i];

		const long double r_e = to_double( N.arcdata[edge].resistance );
		const long double R_e = to_double( G.resistances_accross_cycles[i] );


		const RationalType stretch_e = R_e / r_e;
#ifdef VERBOSE
		cout << R_e << "  " << r_e << " " << R_e / r_e << " " << stretch_e << endl;
#endif
		tree_condition_number += stretch_e;
#ifdef VERBOSE
		cout << "TCN (in) : " << tree_condition_number << endl;
#endif
	}
	return tree_condition_number;
}

/** \brief Calculates the Tree Condition Number
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
long double
SpanningTree<Graph, IntegerType, RationalType>::calculate_tree_condition_number_wrt_unrounded_resistances()
{

#ifdef VERBOSE
	cout << "calculate_tree_condition_number() " << endl << endl;
#endif

	long double tree_condition_number = 0.0;
	for (unsigned int i = 0; i < G.non_tree_edges.size(); i++)
	{

		arc edge = G.non_tree_edges[i];

		const long double r_e =  G.unrounded_resistances[edge] ;
		const long double R_e =  G.unrounded_resistances_accross_cycles[i] ;


		const RationalType stretch_e = R_e / r_e;
#ifdef VERBOSE
		cout << R_e << "  " << r_e << " " << R_e / r_e << " " << stretch_e << endl;
#endif
		tree_condition_number += stretch_e;
#ifdef VERBOSE
		cout << "TCN (in) : " << tree_condition_number << endl;
#endif
	}
	return tree_condition_number;
}

/** \brief finds the Least Common Ancestor
 *
 * Returns the Least Common Ancestor
 *
 * @param node v
 * @param node w
 * @return returns the LCA with node v and node w
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
node
SpanningTree<Graph, IntegerType, RationalType>::find_LCA(
  node v,
  node w)
{
#ifdef VERBOSE
	cout << "Assertion v: " << v << endl;
#endif
	assert( v > 0 );
	assert( w > 0 );

	while ( v != w ) {
		assert( v > 0 );
		assert( w > 0 );

		if ( depth[v] > depth[w] ) {
			const arc a = G.arcs_to_root[v];
			assert( a != 0 );
			assert( v == G.tail(a) );
			v = G.head( a );
		} else {
			const arc a = G.arcs_to_root[w];
			assert( a != 0 );
			assert( w == G.tail(a) );
			w = G.head( a );
		}
	}
	return v;
}


/** \brief finds the Least Common Ancestor
 *
 * Returns the Least Common Ancestor
 *
 * @param arc a
 * @return returns the Least Common Ancestor of the arc a
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
node
SpanningTree<Graph, IntegerType, RationalType>::find_LCA(
  arc a)
{
	node v = G.tail(a);
	node w = G.head(a);

	assert( v > 0 );
	assert( w > 0 );

	while ( v != w ) {

		if ( depth[v] > depth[w] ) {
			const arc a = G.arcs_to_root[v];
			assert( a != 0 );
			assert( v == G.tail(a) );
			v = G.head( a );
		} else {
			const arc a = G.arcs_to_root[w];
			assert( a != 0 );
			assert( w == G.tail(a) );
			w = G.head( a );
		}
	}
	return v;

}

/** \brief Fills in the LCA vector
 *
 * @param x The Primal Solution
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::fill_LCA_vector(vector<RationalType>& x)
{
	for (unsigned int i = 1; i <= G.no_of_edges; i++)
	{
		if (x[i] != 0)
		{
			LCA[i] = find_LCA(i);
		}
	}
}

/** \brief Fills in the LCA vector
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::fill_LCA_vector()
{

	for (unsigned int i = 1; i <= G.no_of_edges; i++)
	{
		LCA[i] = find_LCA(i);
	}
}

/** \brief Gets the incident edges in the tree for each node
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>:: get_tree_incident_edges()
{
	for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
		arc a = G.arcs_to_root[i];
		node u = G.head(a);
		node v = G.tail(a);
		tree_incident_edges[u].push_back(-a);
		tree_incident_edges[v].push_back(a);
	}

#ifdef VERBOSE
	cout << endl << "Tree incident edges..." << endl;
	for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
		cout << i << " : ";
		for (unsigned int j = 0; j < tree_incident_edges[i].size(); j++) {
			cout << tree_incident_edges[i][j] << " , ";
		}
		cout << endl;
	}
#endif

}

/** \brief Set up the tree for Tree Decomposition
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::get_initial_state_of_the_tree_decomposition()
{
	assert(root != 0 );
	const node s = root;
	unsigned int initial_tree_index = 1;
	const node d = find_vertex_separator(s, size_sub_tree_rooted_at[s], initial_tree_index);
	TreeDecompositionInfo TDI;
	TDI.s = s;
	TDI.d = d;
	TDI.size = size_sub_tree_rooted_at[s] + 1;
#ifdef VERBOSE
	cout << "TDI: " << TDI.s << " " << " " << TDI.d << " " << TDI.size << endl;
#endif

	tree_decomposition[0] = TDI;
	assert( tree_decomposition[0].s != 0 );

	unsigned int tree_index = 0;

#ifdef VERBOSE
	const node vertex_separator = d;
	cout << "vertex separator (before anything):" << vertex_separator << endl;
#endif

	find_LCA_and_above_d_flags(TDI.d, tree_index);
}


/** \brief Find the vertex separator
 *
 * Finds and returns the vertex separator of the given Graph
 *
 * @param node s
 * @param int size_at_s
 * @param int tree_index
 * @return vertex separator
 */
template<typename Graph, typename IntegerType, typename RationalType>
node
SpanningTree<Graph, IntegerType, RationalType>::find_vertex_separator(
  node s,
  unsigned int size_at_s,
  unsigned int& tree_index)
{
#if VERBOSE
	cout << s << " " << node_tree_corrospondance[s].back().tree_index << " " << tree_index << endl;
#endif
	assert( node_tree_corrospondance[s].back().tree_index == tree_index - 1 );

	while (true) {
#ifdef VERBOSE
		cout << "find vertex separator called.." << endl;
		cout << " s: " << s << "   " << "size_at_s: " << size_at_s << endl;
#endif
		if (size_at_s == 1) {
			return children[s].front();

		}
#ifdef VERBOSE
		cout << "s: (vertex separator) " << s << " " << "size: " << size_at_s << endl;
#endif
		bool flag = false;
		unsigned int max = 0;
		node v = 0;

		for (auto w : children[s]) {

			NodeInformation current_tree_index = node_tree_corrospondance[w].back();
#ifdef VERBOSE
			cout << "tree index comparison: " << current_tree_index.tree_index << " = " << tree_index - 1 << " , " << w << endl;
#endif
			if (current_tree_index.tree_index != (tree_index - 1)) continue;

			if (2 * (size_sub_tree_rooted_at[w] + 1) > size_at_s) {

				assert( !flag );
				flag = true;
				if (size_sub_tree_rooted_at[w] + 1 > max) {

					max = size_sub_tree_rooted_at[w];
					v = w;
				}
			}
		}

		if (flag == false) {
			assert( s != 0 );
			return s;
		} else {
			s = v;
		}

	}
}

/** \brief Given a node v, Finds the nodes adjacent to it
 *
 * @param v Node v
 *
 * @return Returns a vector of nodes adjacent to v
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
std::vector<node>
SpanningTree<Graph , IntegerType, RationalType>::find_adjacent_nodes(
  node v)
{
	assert( !tree_incident_edges[v].empty() );
	vector<node> adjacent_to_v;
	for ( arc a : tree_incident_edges[v] ) {
		node u = 0;
#ifdef VERBOSE
		cout << "a: " << a << endl;
#endif
		if (G.head(a) == v) {
			u = G.tail(a);
		}
		else {
			u = G.head(a);
		}
		adjacent_to_v.push_back(u);
	}
	return adjacent_to_v;
}


/** \brief Updates the sume of the un-rounded resistances from each node to the root of the tree
 *
 * @param root The Root of the Tree
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void SpanningTree<Graph, IntegerType, RationalType>::update_sum_to_the_root_wrt_unrounded_resistances(node root, Network<Graph, IntegerType, RationalType>& N)
{
	for (unsigned int i = 0; i <= G.no_of_vertices; i++) {
		sum_of_the_unrounded_resistances_from_root_to[i] = 0;
	}
	for ( auto iter = dfs_euler_tour.begin() + discovery[root], end = dfs_euler_tour.begin() + finished[root]; iter != end; ++iter ) {
		const arc a = *iter;
		node v = G.tail(a);
		node w = G.head(a);
//#ifdef VERBOSE
		cout << a << " = ( " << v << " , " << w << " ) " << endl;
//#endif
		if (sum_of_the_unrounded_resistances_from_root_to[w] == RationalType( 0 ) && w != root ) {
//#ifdef VERBOSE
			cout << w << " (before) : " << sum_of_the_unrounded_resistances_from_root_to[w] << endl;
//#endif

			IntegerType xupper = N.arcdata[abs(a)].capacity - N.arcdata[abs(a)].xlower;

			RationalType resistance_a(0);
			if ((N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
			    N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper) && 
			    N.arcdata[abs(a)].infeasibility != 0 && N.arcdata[abs(a)].sroof != 0 )
			{
				resistance_a =   to_double(N.arcdata[abs(a)].resistance_roof);
				
				cout << "resistance_a () = " << to_double(resistance_a) << endl;
			}
			else
			{
#ifdef PathFollowing
			  if(N.arcdata[abs(a)].xlower  != 0 && N.arcdata[abs(a)].slower != 0 && 
			     (N.arcdata[abs(a)].capacity - N.arcdata[abs(a)].xlower != 0) && N.arcdata[abs(a)].supper != 0)
			  {
			resistance_a = to_double(N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper);
			cout << "resistance_a = " << to_double(N.arcdata[abs(a)].resistance_lower) << " + "
						  << to_double(N.arcdata[abs(a)].resistance_upper) << endl;
			  }
			  else
			  {
			    resistance_a =  to_double(N.arcdata[abs(a)].resistance_roof);
			  }
#else
				resistance_a = RationalType(1) / to_double(N.arcdata[abs(a)].xlower * N.arcdata[abs(a)].xlower) +
				               RationalType(1) / to_double(xupper * xupper);
#endif
					  
					       
				cout << "resistance_a = " << resistance_a << endl;
			}
			sum_of_the_unrounded_resistances_from_root_to[w] +=  sum_of_the_unrounded_resistances_from_root_to[v] +
			    resistance_a;

//#ifdef VERBOSE

		      cout << to_double(sum_of_the_resistances_from_root_to[w]) << " += " << 
			      to_double(sum_of_the_resistances_from_root_to[v]) << " + " <<
			      to_double(resistance_a) << endl;
			      
//#endif
		}
	}
//#ifdef VERBOSE
	print_resistances(N);
	
	cout << "Sum of UnRounded Resistances from Root: " << endl;
	
	for (unsigned int i = 1; i <= G.no_of_vertices; i++) {
		cout << i << " : " << to_double(sum_of_the_unrounded_resistances_from_root_to[i]) << endl;
	}
//#endif
}





/** \brief Finds the sum of the resistances from the nodes to the root
 *
 * @param root The root of the tree
 *
 */

template<typename Graph, typename IntegerType, typename RationalType>
void SpanningTree<Graph, IntegerType, RationalType>::update_sum_to_the_root(node root, Network<Graph, IntegerType, RationalType>& N)
{


	vector<unsigned int>visited(N.G.no_of_vertices + 1, 0);
	for (unsigned int i = 0; i <= N.G.no_of_vertices; i++) {
		sum_of_the_resistances_from_root_to[i] = 0;
	}


	for (auto iter = dfs_euler_tour.begin() + discovery[root], end = dfs_euler_tour.begin() + finished[root]; iter != end; ++iter) {
		const arc a = *iter;
		node v = G.tail(a);
		node w = G.head(a);
#ifdef VERBOSE
		cout << a << " = ( " << v << " , " << w << " ) " << endl;
#endif
		if (visited[w] == 0 && w != root) {
#ifdef VERBOSE
			cout << w << " (before) : " << to_double(sum_of_the_resistances_from_root_to[w]) << endl;
#endif
//      if(N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower + N.arcdata[abs(a)].resistance_upper)
			if (N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_lower ||
			    N.arcdata[abs(a)].resistance_roof < N.arcdata[abs(a)].resistance_upper)
			{
				sum_of_the_resistances_from_root_to[w] +=  sum_of_the_resistances_from_root_to[v] +
				    N.arcdata[abs(a)].resistance_roof;
			}
			else
			{
				sum_of_the_resistances_from_root_to[w] +=  sum_of_the_resistances_from_root_to[v]
				    + (N.arcdata[abs(a)].resistance_upper + N.arcdata[abs(a)].resistance_lower);
			}
#ifdef VERBOSE
			cout << w << " : " << to_double(sum_of_the_resistances_from_root_to[w]);
			cout << " = " << to_double(sum_of_the_resistances_from_root_to[v]);
			cout << " +  " << to_double(N.arcdata[abs(a)].resistance_upper + N.arcdata[abs(a)].resistance_lower) <<  endl;
#endif
			visited[w] = 1;
		}
	}
#ifdef VERBOSE
	cout << endl << "number of verticies = " << N.G.no_of_vertices << "  Root = " << root << endl;
	for (unsigned int i = 1; i <= N.G.no_of_vertices; i++) {
		cout << "sum to root: " << i << " : " << to_double(sum_of_the_resistances_from_root_to[i]) << endl;
	}
#endif
}

/** \brief Find the LCA and above_d flags
 *
 * Finds the LCA and the above_d flags for nodes in the tree with the given tree index
 *
 * @param node d
 * @param int tree_index
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::find_LCA_and_above_d_flags(
  node d,
  unsigned int tree_index)
{
	const node s = tree_decomposition[tree_index].s;
	node current_node = d;
	node previous_node = 0;
	node next_node = parent[current_node];
#ifdef VERBOSE
	cout << "current node: " << current_node << " next node: " << next_node << endl;
	cout << "s: " << s << " @ " << finished[s] << ", d: " << d << " @ " << discovery[d] << endl;
#endif
	assert( s != 0 );
	assert( d != 0 );

	while ( current_node != 0 ) {
		for ( auto iter = dfs_euler_tour.begin() + discovery[current_node], end = dfs_euler_tour.begin() + finished[current_node]; iter != end; ++iter ) {
			const arc a = *iter;

#ifdef VERBOSE
			cout << a << " = (" << G.tail(a) << "," << G.head(a) << ") " << endl;
#endif
			assert( a != 0 );
			const node v = G.tail(a);
			const node w = G.head(a);
			assert( v != 0 );
			assert( w != 0 );

			node_tree_corrospondance[v].back().lca = current_node;
			if (current_node == d) {
				node_tree_corrospondance[v].back().above_d = -1;
			} else {
				node_tree_corrospondance[v].back().above_d = 0;
			}
			if ( w == previous_node || node_tree_corrospondance[w].back().tree_index != tree_index) {

#ifdef VERBOSE
				cout << "jump " << finished[w] - discovery[w] + 1 << endl;
#endif

				iter += finished[w] - discovery[w] + 1;
				continue;
			}
			node_tree_corrospondance[w].back().lca = current_node;
			if (current_node == d) {
				node_tree_corrospondance[w].back().above_d = -1;
			} else {
				node_tree_corrospondance[w].back().above_d = 0;
			}
		}
		previous_node = current_node;
		current_node = next_node;
		next_node = parent[current_node];

#ifdef VERBOSE
		cout << "current node: " << current_node << " next node: " << next_node << endl;
#endif
	}
	node_tree_corrospondance[s].back().above_d = 0;
	node_tree_corrospondance[s].back().lca = s;
	node_tree_corrospondance[d].back().above_d = 1;
	node_tree_corrospondance[d].back().lca = d;
}

/** \brief Do the Tree Decomposition
 *
 * Does a tree decomposition of the tree rooted at the given node and labeled with the given tree index and returns the vertex separator
 *
 * @param node s Tree rooted at
 * @param int tree_index
 * @return the vertex separator
 */
template<typename Graph, typename IntegerType, typename RationalType>
node
SpanningTree<Graph, IntegerType, RationalType>::tree_decompose(
  node s,
  unsigned int& tree_index)
{
	node d = tree_decomposition[tree_index - 1].d;

#ifdef VERBOSE
	cout << "d : " << d << endl;
#endif

	std::vector<node> adjacent_to_d;

	if (d == s)
	{
		return d;
	}
	int req_tree_index = tree_index + 1;

	label_the_tree_rooted_at(s, d, req_tree_index, tree_index);

#ifdef VERBOSE
	print_sizes();
	print_children();
#endif

	tree_index++;
	assert( tree_index < 3 * G.no_of_vertices );

	TreeDecompositionInfo TDI;
	TDI.size = size_sub_tree_rooted_at[s] + 1;
	TDI.s = s;
	TDI.d = find_vertex_separator(TDI.s, TDI.size - 1, tree_index);
#ifdef VERBOSE
	cout << "TDI , " << tree_index - 1 << " : " << " ( " << TDI.size - 1 << " , " << TDI.s << " , " << TDI.d << " ) " ;
#endif
	tree_decomposition[tree_index - 1] = TDI;
	assert( tree_decomposition[tree_index - 1].s != 0 );
	find_LCA_and_above_d_flags(TDI.d, tree_index - 1);

	return d;
}

/** \brief Prints the children of each node
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::print_children()
{
	cout << "children " << endl;
	for (unsigned int i = 1; i <= G.no_of_vertices; i++)
	{	cout << i << " : " ;
		for (unsigned int j = 0; j < children[i].size(); j++)
		{
			cout << children[i][j] << ", ";
		}
		cout << endl;
	}

}

/** \brief Prints the sizes of the sub-tree rooted at each node
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::print_sizes()
{
	cout << endl << "Sizes of the sub trees" << endl;
	for (unsigned int i = 1; i <= G.no_of_vertices; i++)
	{
		cout << i << " : " << size_sub_tree_rooted_at[i] << endl;
	}
}


/** \brief Label the part of the tree which is rooted at the given node
 *
 * Labels the part of the tree rooted at the given node with the appropriate tree index
 *
 * @param node y
 * @param node d1
 * @param int
 * @param int
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::label_the_tree_rooted_at(
  node y,
  node d1,
  unsigned int req_tree_index,
  unsigned int &tree_index
) {

#ifdef VERBOSE
	print_children();
	cout << "The Euler Tour: " << y << " " << discovery[y] << "  " << finished[y] << endl;

	for ( auto iter = dfs_euler_tour.begin() + discovery[y], end = dfs_euler_tour.begin() + finished[y]; iter != end; ++iter ) {
		const arc a = *iter;
		cout << a << " = (" << G.tail(a) << "," << G.head(a) << ") " << endl;
	}
	cout << endl << endl;
#endif


	children[y].clear();
	size_sub_tree_rooted_at[y] = 0;
	children[d1].clear();
	size_sub_tree_rooted_at[d1] = 0;

#ifdef VERBOSE
	cout << y << " " << d1 << " " << discovery[y] << " " << finished[y] << endl;
#endif

	if (discovery[y] == finished[y]) {
#ifdef VERBOSE
		cout << endl << "y : " << y << endl;
#endif
		NodeInformation node_info;
		node_info.tree_index = tree_index;
		node_info.above_d = -2;
		node_tree_corrospondance[y].push_back(node_info);
#ifdef VERBOSE
		cout << "Update Node Correspondance: " << y << " -> " << tree_index << endl;
#endif
	}

	for ( auto iter = dfs_euler_tour.begin() + discovery[y], end = dfs_euler_tour.begin() + finished[y]; iter != end; ++iter ) {
		const arc a = *iter;

#ifdef VERBOSE
		cout << a << " = (" << G.tail(a) << "," << G.head(a) << ") " << endl;
#endif
		assert( a != 0 );
		const node v = G.tail(a);
		const node w = G.head(a);
		assert( v != 0 );
		assert( w != 0 );

		if (node_tree_corrospondance[v].back().tree_index != tree_index) {
#ifdef VERBOSE
			cout << "Update Node Correspondance: " << v << " -> " << tree_index << endl;
#endif
			NodeInformation node_info;
			node_info.tree_index = tree_index;
			node_info.above_d = -2;
			assert(v < node_tree_corrospondance.size());
			node_tree_corrospondance[v].push_back(node_info);
			children[v].clear();
			size_sub_tree_rooted_at[v] = 0;
		}

#ifdef VERBOSE
		cout << w << " == " << d1 << " || ";
		cout <<  " ( " << node_tree_corrospondance[w].back().tree_index;
		cout << "!=" << tree_index << " && " << node_tree_corrospondance[w].back().tree_index;
		cout << " != " << req_tree_index - 2 << " ) " << endl;
#endif
		if (w == d1 || ( node_tree_corrospondance[w].back().tree_index != tree_index &&  node_tree_corrospondance[w].back().tree_index != req_tree_index - 2) ) {

#ifdef VERBOSE
			cout << "skip " << w << " by " << finished[w] << " - " << discovery[w] << " + 1 = " << finished[w] - discovery[w] + 1 << endl;
#endif
			iter += finished[w] - discovery[w] + 1;

			continue;
		}
#ifdef VERBOSE
		cout << endl << "w : " << w << endl;
#endif

		if (node_tree_corrospondance[w].back().tree_index != tree_index) {
			NodeInformation node_info;
			node_info.tree_index = tree_index;
			node_info.above_d = -2;
			node_tree_corrospondance[w].push_back(node_info);
#ifdef VERBOSE
			cout << "Update Node Correspondance: " << w << " -> " << tree_index << endl;
#endif
			children[w].clear();
			size_sub_tree_rooted_at[w] = 0;
		} else {
			size_sub_tree_rooted_at[w] += size_sub_tree_rooted_at[v] + 1;
			children[w].push_back(v);
		}
	}
}

/** \brief Runs the init function on the graph below the vertex separator d
 *
 * @param d The vertex separator
 * @param s The root
 * @param tree_index
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::init_below(
  node d,
  node s,
  unsigned int &tree_index)
{
#ifdef VERBOSE
	cout << "processing the other (not the T_0) part of the tree  " << S.size() << endl;
#endif

	Stack_Data d_s = S.top();
	node d1 = d_s.d;

	unsigned int req_tree_index = d_s.tree_index;

#ifdef VERBOSE
	cout << "d1: " << d1 << endl;
#endif
	std::vector<node> incident_to_d1;

	for (auto a : tree_incident_edges[d1]) {
		node v = (G.head(a) == d1) ? G.tail(a) : G.head(a);
		incident_to_d1.push_back(v);
	}

#ifdef VERBOSE
	cout << "d1: " << d1 << " ,  " << req_tree_index << endl ;
#endif

	if (d == s) req_tree_index++;
	assert( req_tree_index < 3 * G.no_of_vertices );

	for (auto y : incident_to_d1) {
#ifdef VERBOSE
		cout << (node_tree_corrospondance[y].back()).tree_index << " = " << req_tree_index - 2 << " ,  " << y << endl;
#endif
		if ((node_tree_corrospondance[y].back()).tree_index != req_tree_index - 2) continue;

#ifdef VERBOSE
		cout << "Y: " << y << endl;
#endif

		label_the_tree_rooted_at(y, d1, req_tree_index, tree_index);

		size_sub_tree_rooted_at[d1] = 1;

		for (auto v : incident_to_d1) {
#ifdef VERBOSE
			cout << (node_tree_corrospondance[v].back()).tree_index << "!= " << tree_index << endl;
#endif
			if ((node_tree_corrospondance[v].back()).tree_index != tree_index) continue;
			size_sub_tree_rooted_at[d1] += size_sub_tree_rooted_at[v];
#ifdef VERBOSE
			cout << "intermediate size: " << size_sub_tree_rooted_at[d1] << endl;
#endif
			children[d1].push_back(v);
		}

		NodeInformation node_info;
		node_info.tree_index = tree_index;
		node_info.above_d = -2;

		node_tree_corrospondance[d1].push_back(node_info);

#ifdef VERBOSE
		cout << "Update Node Correspondance: " << d1 << "-> " << tree_index << endl;
#endif

		tree_index++;
		assert( tree_index < 3 * G.no_of_vertices );

		assert( d1 != 0 );
		TreeDecompositionInfo TDI;
		TDI.size = size_sub_tree_rooted_at[y] + 2;
		TDI.s = d1;
		TDI.d = find_vertex_separator(d1, TDI.size - 1, tree_index);
#ifdef VERBOSE
		cout << "TDI , " << tree_index - 1 << " : " << " ( " << TDI.size - 1 << " , " << TDI.s << " , " << TDI.d << " ) " ;
#endif
		tree_decomposition[tree_index - 1] = TDI;
		assert( tree_decomposition[tree_index - 1].s != 0 );


		find_LCA_and_above_d_flags(TDI.d, tree_index - 1);


#ifdef VERBOSE
		print_sizes();
#endif
		if (size_sub_tree_rooted_at[d1] >= 2) {
			init(d1, tree_index);
		}
	}

	S.pop();
#ifdef VERBOSE
	cout << "(size: ) " << S.size();
	print_sizes();
#endif
}

/** \brief Runs the init function for the tree-datastructe
 *
 * @param s Root
 * @param tree_index
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::init(node s,
    unsigned int& tree_index
                                                    )
{

#ifdef VERBOSE
	print_children();
	cout << "tree decompose called: " << s << " , "  << tree_index << endl;
#endif

	node d = tree_decompose(s, tree_index);

	Stack_Data stack_data;
	stack_data.d = d;
	stack_data.s = s;
	stack_data.tree_index = tree_index;

#ifdef VERBOSE
	cout << "pushed into the stack S: " << "( " << stack_data.d << " , " << stack_data.s << " , " << stack_data.tree_index << " ) " ;
#endif

	S.push(stack_data);

	if (size_sub_tree_rooted_at[s] >= 2 && d != s) {
		init(s, tree_index);
	}
#ifdef VERBOSE
	cout << "size: " << S.size() << endl;
#endif

	if (!S.empty()) {
		init_below( d, s, tree_index);
	}
}


/** \brief Returns the voltage drop accross the the arc
 *
 * @param edge_index The edge_index of the arc
 *
 * @return Returns the voltage drop
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
IntegerType
SpanningTree<Graph, IntegerType, RationalType>::query (
  unsigned int edge_index) const
{
	const arc a = G.non_tree_edges[edge_index];
	const unsigned int index = non_tree_edge_common_tree_index[edge_index];
	const node u = G.tail(a);
	const node v = G.head(a);

	const IntegerType v_u = query(u, index);
	const IntegerType v_w = query(v, index);

	IntegerType voltage_drop = v_u - v_w;

	return voltage_drop;
}

template<typename Graph, typename IntegerType, typename RationalType>
IntegerType
SpanningTree<Graph, IntegerType, RationalType>::query (
  node v,
  unsigned int index) const
{
#ifdef VERBOSE
	cout << "query called,  v = " << v << endl;
#endif
	IntegerType drop = 0;
#ifdef VERBOSE
	const TreeDecompositionInfo& TDI0 = tree_decomposition[0];
	const int& size0 = TDI0.size;
	cout << "size: " << size0 << endl;
	cout << "d : " << TDI0.d << endl;
	cout << "s : " << TDI0.s << endl;
#endif
	for (auto iter = node_tree_corrospondance[v].begin() + index, end = node_tree_corrospondance[v].end(); iter != end; ++iter) {

		const NodeInformation& current_tree = *iter;
		unsigned int current_tree_index = current_tree.tree_index;

		const TreeDecompositionInfo& TDI  = tree_decomposition[current_tree_index];
		const int& size = TDI.size;
		if (size <= 1) break;
#ifdef VERBOSE
		cout << "current_tree: " << current_tree_index << endl;
		cout << "size: " << size << endl;
#endif
		node d = TDI.d;
#ifdef VERBOSE
		cout << "d: " << d << endl;
		cout << "drop: " << to_double(drop) << endl;
#endif
		if (d == v) {
			drop += d_drop[current_tree_index];
#ifdef VERBOSE
			cout << "case d == v" << endl;
			cout << endl << "d_drop: [" << current_tree_index << " ] =  " << to_double(drop) << endl;
#endif
			break;
		}

		if (size == 2) {
			break;
		}

		const node lca = current_tree.lca;
		const node s = TDI.s;
#ifdef VERBOSE
		cout << "lca: " << lca << "  " << "s: " << s << endl;
#endif
		assert( lca != 0 );
		assert( s != 0 );

		if (current_tree.above_d == 0) {
			IntegerType height = sum_of_the_resistances_from_root_to[lca] - sum_of_the_resistances_from_root_to[s];
#ifdef VERBOSE
			cout << "case T_0: " << endl;
			IntegerType drop_before = drop;
			cout << to_double( d_ext[current_tree_index] ) << " * " << to_double( height ) << endl;
#endif

			height *= d_ext[current_tree_index];

			//cout << "height later: " << to_double(height) << endl;

			drop += height;

#ifdef VERBOSE
			cout << "drop: " << to_double( drop ) << " = " << to_double( drop_before ) << " + " << to_double( height ) << endl;
#endif
		} else {
#ifdef VERBOSE
			cout << "else case: " << endl;
#endif
			drop += d_drop[current_tree_index];
#ifdef VERBOSE
			cout << "drop: " << to_double(drop) << endl;
#endif
		}
	}
#ifdef VERBOSE
	cout << "value returned from query: " << to_double(drop) << endl;
#endif
	return drop;
}

/** \brief Runs the update function
 *
 * @param v Node on which update being done
 * @param alpha Amount by which update being done
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::update(
  node v,
  const IntegerType& alpha)
{
	update( v, alpha, 0 );
}


/** \brief Runs the update function
 *
 * @param v The node on which the update function is run
 * @param alpha The amount by which update has to be done
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::update(
  node v,
  const IntegerType& alpha,
  unsigned int index
)
{

#ifdef VERBOSE
	cout << "update called, " << v << endl;
	cout << "alpha: " << to_double(alpha) << endl;
#endif

	for (auto iter = node_tree_corrospondance[v].begin() + index, end = node_tree_corrospondance[v].end(); iter != end; ++iter) {
		const NodeInformation& current_tree = *iter;
		unsigned current_tree_index = current_tree.tree_index;

#ifdef VERBOSE
		cout << "current tree (in update): " << current_tree_index << endl;
#endif

		assert( current_tree_index < tree_decomposition.size() );
		const node lca = current_tree.lca;
		const TreeDecompositionInfo& TDI = tree_decomposition[current_tree_index];
		const node s = TDI.s;
#ifdef VERBOSE
		cout << "lca: " << lca << "  " << "s: " << s << endl;
#endif
		IntegerType height = sum_of_the_resistances_from_root_to[lca];
		height -= sum_of_the_resistances_from_root_to[s];

#ifdef VERBOSE
		cout << "Height picked: " << to_double(height) << endl;
#endif


		height *= alpha;

#ifdef VERBOSE
		cout << "height = " << to_double(height) << endl;
		cout << "d_drop[" << current_tree_index << "] = " << d_drop[current_tree_index] << endl;
#endif
		d_drop[current_tree_index] += height;

#ifdef VERBOSE
		cout << "new d_drop: " << to_double(d_drop[current_tree_index]);
#endif
		const unsigned size = TDI.size;

		if (size <= 2) {
			break;
		}

		const node d = TDI.d;

		if ( current_tree.above_d != 0) {
#ifdef VERBOSE
			cout << "case not in T_0: " << endl;
#endif
			d_ext[current_tree_index] += alpha;
#ifdef VERBOSE
			cout << "updated d_ext: " << to_double(d_ext[current_tree_index]) << endl;
#endif
		}
#ifdef VERBOSE
		cout << endl << v << " != " << d << endl;
#endif

		if (v != d) {
			continue;
		}
		else {
			break;
		}
	}
}

/** \brief Runs the update function when the Data Structure to store the three graphs is used
 *
 * @param v node
 * @param index Index of the edge being updated
 * @param height_offset
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::update_with_DS(
  node v,
  IntegerType& alpha,
  unsigned int index,
  IntegerType height_offset
)

{
#ifdef VERBOSE
	cout << "update called, " << v << endl;
	cout << "alpha: " << to_double(alpha) << endl;
#endif

	for (auto iter = node_tree_corrospondance[v].begin() + index, end = node_tree_corrospondance[v].end(); iter != end; ++iter ) {
		const NodeInformation& current_tree = *iter;
		unsigned current_tree_index = current_tree.tree_index;
#ifdef VERBOSE
		cout << "current tree (in update): " << current_tree_index << endl;
#endif
		assert( current_tree_index < tree_decomposition.size() );
		const node lca = current_tree.lca;
		const TreeDecompositionInfo& TDI = tree_decomposition[current_tree_index];
		const node s = TDI.s;
#ifdef VERBOSE
		cout << "lca: " << lca << "  " << "s: " << s << endl;
#endif
		IntegerType height = sum_of_the_resistances_from_root_to[lca];
#ifdef VERBOSE
		cout << "sum_of_the_resistances_from_root_to[lca] = " << to_double(sum_of_the_resistances_from_root_to[lca]) << endl;
#endif
		height -= sum_of_the_resistances_from_root_to[s];
#ifdef VERBOSE
		cout << " sum_of_the_resistances_from_root_to[s] = " << to_double(sum_of_the_resistances_from_root_to[s]) << endl;
#endif
#ifdef VERBOSE
		cout << "size = " << TDI.size;
#endif

#ifdef VERBOSE
		cout << "height before = " << to_double(height) << endl;
		cout << "d = " << TDI.d << endl;
		cout << "s = " << TDI.s << endl;
#endif

		if (height_offset > 0) {
			node lca1 = find_LCA(v, TDI.d);
			if ( (height != 0 && (current_tree.above_d == 0 || v == TDI.d)) && v == lca1) {
				height -= height_offset;
			}
		}

#ifdef VERBOSE
		cout << "Height picked: " << to_double(height) << endl << "height_offset = " << to_double(height_offset) << endl;
#endif

		height *= alpha;
#ifdef VERBOSE
		cout << "height = " << to_double(height) << endl;
		cout << "d_drop[" << current_tree_index << "] = " << to_double(d_drop[current_tree_index]) << endl;
#endif
		d_drop[current_tree_index] += height;

#ifdef VERBOSE
		cout << "new d_drop: " << to_double(d_drop[current_tree_index]) << endl;
#endif
		const unsigned size = TDI.size;

		if (size <= 2) {
			break;
		}

		const node d = TDI.d;

		if ( current_tree.above_d != 0) {
#ifdef VERBOSE
			cout << "case not in T_0: " << endl;
#endif
			d_ext[current_tree_index] += alpha;
#ifdef VERBOSE
			cout << "updated d_ext: " << to_double(d_ext[current_tree_index]) << endl;
#endif
		}
#ifdef VERBOSE
		cout << endl << v << " != " << d << endl;
#endif

		if (v != d) {
			continue;
		}
		else {
			break;
		}
	}
}

/** \brief Runs the update function over the arc
 *
 * @param alpha Amount to be updated
 * @param edge_index Index of the edge being updated
 *
 */
template<typename Graph, typename IntegerType, typename RationalType>
void
SpanningTree<Graph, IntegerType, RationalType>::update(
  const IntegerType& alpha,
  unsigned int edge_index
) {

	const arc a = G.non_tree_edges[edge_index];

	const node v = G.tail( a );
	const node w = G.head( a );

	const unsigned int index = non_tree_edge_common_tree_index[edge_index];

	update( v,  alpha, index );
	update( w, -alpha, index );
}

/** \brief Given the node v and arc (v,u) returns u
 *
 * @return u
 *
 */
template<typename IntegerType, typename RationalType>
node Graph<IntegerType, RationalType>:: get_the_other_end_of_arc(node v, arc a)
{
	if (head(a) == v) {
		return tail(a);
	}
	else {
		return head(a);
	}
}


template<typename IntegerType, typename RationalType>
struct node_distance {
	node u;
	RationalType shortest_distance;
};


template<typename IntegerType, typename RationalType>
class compare_distance {


public:

	bool operator() (const node_distance<IntegerType, RationalType>& d1, const node_distance<IntegerType, RationalType>& d2) const
	{
		if (d1.shortest_distance > d2.shortest_distance)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};


//     /** \brief Runs the Dijkstra
//      *
//      * @param source The Source Node
//      * @param target The Target Node
//      *
//      * @return Returns the shortest distance between source and target
//      *
//      */
//     template<typename IntegerType, typename RationalType>
//     RationalType
//     Graph<IntegerType, RationalType>::dijkstra(node source,
// 					       node target
//     )
//     {
//       #ifdef VERBOSE
//       cout << "dijkstra()" << endl;
//       #endif
//
//       #ifdef VERBOSE
//       cout << "source :  " << source << "   target: " << target << endl;
//       #endif
//
//       vector<unsigned int> mark_nodes(2*no_of_vertices + 1 , 0);
//
//       vector<RationalType> dist(2*no_of_vertices + 1, numeric_limits<RationalType>::max());
//       vector<bool> visited(2*no_of_vertices + 1, false);
//       dist[source] = 0;
//
//       priority_queue<node_distance<IntegerType, RationalType>, vector<node_distance<IntegerType, RationalType> >, compare_distance<IntegerType, RationalType> > Q;
//
//       node_distance<IntegerType, RationalType> first = {source , 0};
//       Q.push(first);
//
//       #ifdef VERBOSE
//       cout << "while will begin" << endl;
//       #endif
//       while (!Q.empty()){
//
// 	node_distance<IntegerType, RationalType> temp = Q.top();
// 	Q.pop();
// 	node u = temp.u;
//
// 	for(arc a : incident_edges[u]){
// 	  #ifdef VERBOSE
// 	  cout << "u: " << u << endl;
// 	  cout << "a: " << a << endl;
// 	  #endif
// 	  node v = get_the_other_end_of_arc(u,a);
// 	  #ifdef VERBOSE
// 	  cout << "v: " << v << endl;
// 	  #endif
//
//
// 	  #ifdef VERBOSE
// 	  cout << (dist[u]) << " + " << to_double(resistances[abs(a)]) << " < " << (dist[v]) << endl;
// 	  #endif
// 	  if( dist[u] + to_double(resistances[abs(a)]) < dist[v])
// 	  {
//
// 	    dist[v] =  dist[u] + to_double(resistances[abs(a)]);
// 	    #ifdef VERBOSE
// 	    cout << "dist[" << v << "] = " << (dist[v]) << endl;
// 	    #endif
// 	    node_distance<IntegerType, RationalType> new_node;
// 	    new_node.u = v;
// 	    new_node.shortest_distance = dist[v];
// 	    Q.push(new_node);
// 	  }
// 	}
//       }
//
//
//       RationalType d = dist[target];
//
//       #ifdef VERBOSE
//       cout << "before returning let us check the distances ... " << endl;
//       for(unsigned int i = 1; i <= no_of_vertices; i++)
//       {
// 	cout << i << ": " << to_double(dist[i]) << endl;
//       }
//       #endif
//
//       // cout << "returns .. " << d << endl;
//
//       return d;
//
//     }
//
//
//     template<typename IntegerType, typename RationalType>
//     vector<RationalType> Graph<IntegerType, RationalType>::dijkstra(node source)
//     {
//       #ifdef VERBOSE
//       cout << "dijkstra()" << endl;
//       #endif
//
//       #ifdef VERBOSE
//       cout << "source :  " << source  << endl;
//       #endif
//
//       vector<unsigned int> mark_nodes(2*no_of_vertices + 1 , 0);
//
//       vector<RationalType> dist(2*no_of_vertices + 1, numeric_limits<RationalType>::max());
//       vector<bool> visited(2*no_of_vertices + 1, false);
//       dist[source] = 0;
//
//       priority_queue<node_distance<IntegerType, RationalType>, vector<node_distance<IntegerType, RationalType> >,
// compare_distance<IntegerType, RationalType> > Q;
//
//       node_distance<IntegerType, RationalType> first = {source , 0};
//       Q.push(first);
//
//       #ifdef VERBOSE
//       cout << "while will begin" << endl;
//       #endif
//       while (!Q.empty())
//       {
//
// 	node_distance<IntegerType, RationalType> temp = Q.top();
// 	Q.pop();
// 	node u = temp.u;
//
// 	for(arc a : incident_edges[u])
// 	{
// 	  #ifdef VERBOSE
// 	  cout << "u: " << u << endl;
// 	  cout << "a: " << a << endl;
// 	  #endif
// 	  node v = get_the_other_end_of_arc(u,a);
// 	  #ifdef VERBOSE
// 	  cout << "v: " << v << endl;
// 	  #endif
//
// 	  #ifdef VERBOSE
// 	  cout << (dist[u]) << " + " << to_double(resistances[abs(a)]) << " < " << (dist[v]) << endl;
// 	  #endif
// 	  if( dist[u] + to_double(resistances[abs(a)]) < dist[v])
// 	  {
//
// 	    dist[v] =  dist[u] + to_double(resistances[abs(a)]);
// 	    #ifdef VERBOSE
// 	    cout << "dist[" << v << "] = " << (dist[v]) << endl;
// 	    #endif
//
// 	    node_distance<IntegerType, RationalType> new_node;
// 	    new_node.u = v;
// 	    new_node.shortest_distance = dist[v];
// 	    Q.push(new_node);
// 	  }
// 	}
//       }
//
//
//       #ifdef VERBOSE
//       cout << "before returning let us check the distances ... " << endl;
//       for(unsigned int i = 1; i <= no_of_vertices; i++)
//       {
// 	cout << i << ": " << to_double(dist[i]) << endl;
//       }
//       #endif
//
//
//       return dist;
//
//     }
//
//
//
//     template<typename IntegerType, typename RationalType>
//     vector<RationalType> Graph<IntegerType, RationalType>::dijkstra(node source,
// 								    vector<node>& X
//     )
//     {
//       #ifdef VERBOSE
//       cout << "dijkstra()" << endl;
//       #endif
//
//       #ifdef VERBOSE
//       cout << "source :  " << source << endl;
//       #endif
//
//       vector<unsigned int> mark_nodes(2*no_of_vertices + 1 , 0);
//       #ifdef VERBOSE
//       cout << "checking the resistances: " << endl;
//       for(unsigned int i=1; i <= no_of_edges; i++)
//       {
// 	cout << i << " : "  << weights[i] << " " << to_double(resistances[i]) << endl;
//       }
//       cout << endl;
//       #endif
//
//       cout << "mark nodes ... " << source << endl << endl;
//       for(auto x: X)
//       {
// 	cout << x << " , ";
//       }
//       cout << endl;
//       for(node v : X)
//       {
// 	mark_nodes[v] = 1;
//       }
//       cout << " node are marked" << endl;
//       vector<RationalType> dist(2*no_of_vertices + 1, numeric_limits<RationalType>::max());
//       vector<bool> visited(2*no_of_vertices + 1, false);
//       dist[source] = 0;
//
//       priority_queue<node_distance<IntegerType, RationalType>, vector<node_distance<IntegerType, RationalType> >, compare_distance<IntegerType, RationalType> > Q;
//
//       node_distance<IntegerType, RationalType> first = {source , 0};
//       Q.push(first);
//
//       #ifdef VERBOSE
//       cout << "while will begin" << endl;
//       #endif
//       while (!Q.empty())
//       {
//
// 	node_distance<IntegerType, RationalType> temp = Q.top();
// 	Q.pop();
// 	node u = temp.u;
//
// 	for(arc a : incident_edges[u])
// 	{
// 	  #ifdef VERBOSE
// 	  cout << "u: " << u << endl;
// 	  cout << "a: " << a << endl;
// 	  #endif
// 	  node v = get_the_other_end_of_arc(u,a);
// 	  #ifdef VERBOSE
// 	  cout << "v: " << v << endl;
// 	  #endif
//
// 	  if(mark_nodes[v] != 1) continue;
//
// 	  #ifdef VERBOSE
// 	  cout << (dist[u]) << " + " << weights[abs(a)] << " < " << (dist[v]) << endl;
// 	  #endif
// 	  if( dist[u] + weights[abs(a)] < dist[v])
// 	  {
//
// 	    dist[v] =  dist[u] + weights[abs(a)];
// 	    #ifdef VERBOSE
// 	    cout << "dist[" << v << "] = " << (dist[v]) << endl;
// 	    #endif
// 	    node_distance<IntegerType, RationalType> new_node;
// 	    new_node.u = v;
// 	    new_node.shortest_distance = dist[v];
// 	    Q.push(new_node);
// 	  }
// 	}
//       }
//
//
//       #ifdef VERBOSE
//       cout << "source: " << source << endl;
//       cout << "before returning let us check the distances ... " << endl;
//       for(unsigned int i = 1; i <= no_of_vertices; i++)
//       {
// 	cout << i << ": " <<  dist[i] << endl;
//       }
//       #endif
//
//
//       return dist;
//
//     }
//
//
//
//
//     template<typename IntegerType, typename RationalType>
//     RationalType Graph<IntegerType, RationalType>::dijkstra(node source,
// 							    node target,
// 							    vector<node>& X,
// 							    vector<node>& previous
//     )
//     {
//       #ifdef VERBOSE
//       cout << "dijkstra()" << endl;
//       #endif
//
//       #ifdef VERBOSE
//       cout << "source :  " << source << "   target: " << target << endl;
//       #endif
//
//       vector<unsigned int> mark_nodes(2*no_of_vertices + 1 , 0);
//       #ifdef VERBOSE
//       cout << "checking the resistances: " << endl;
//       for(unsigned int i=1; i <= no_of_edges; i++)
//       {
// 	cout << i << " : "  << weights[i] << " " << to_double(resistances[i]) << endl;
//       }
//       cout << endl;
//       #endif
//
//       for(node v : X)
//       {
// 	mark_nodes[v] = 1;
//       }
//       vector<RationalType> dist(2*no_of_vertices + 1, numeric_limits<RationalType>::max());
//       vector<bool> visited(2*no_of_vertices + 1, false);
//       dist[source] = 0;
//
//       priority_queue<node_distance<IntegerType, RationalType>, vector<node_distance<IntegerType, RationalType> >, compare_distance<IntegerType, RationalType> > Q;
//
//       node_distance<IntegerType, RationalType> first = {source , 0};
//       Q.push(first);
//
//       #ifdef VERBOSE
//       cout << "while will begin" << endl;
//       #endif
//       while (!Q.empty())
//       {
//
// 	node_distance<IntegerType, RationalType> temp = Q.top();
// 	Q.pop();
// 	node u = temp.u;
//
// 	for(arc a : incident_edges[u])
// 	{
// 	  #ifdef VERBOSE
// 	  cout << "u: " << u << endl;
// 	  cout << "a: " << a << endl;
// 	  #endif
// 	  node v = get_the_other_end_of_arc(u,a);
// 	  #ifdef VERBOSE
// 	  cout << "v: " << v << endl;
// 	  #endif
//
// 	  if(mark_nodes[v] != 1) continue;
//
// 	  #ifdef VERBOSE
// 	  cout << (dist[u]) << " + " << weights[abs(a)] << " < " << (dist[v]) << endl;
// 	  #endif
// 	  if( dist[u] + weights[abs(a)] < dist[v])
// 	  {
//
// 	    dist[v] =  dist[u] + weights[abs(a)];
// 	    #ifdef VERBOSE
// 	    cout << "dist[" << v << "] = " << (dist[v]) << endl;
// 	    #endif
//
// 	    previous[v] = u;
// 	    #ifdef VERBOSE
// 	    cout << "previous[" << v << "] = " << previous[v] << endl;
// 	    #endif
//
// 	    node_distance<IntegerType, RationalType> new_node;
// 	    new_node.u = v;
// 	    new_node.shortest_distance = dist[v];
// 	    Q.push(new_node);
//
// 	  }
// 	}
//       }
//
//
//       RationalType d = dist[target];
//
//       #ifdef VERBOSE
//       cout << "vector X: ";
//       for(auto a : X)
//       {
// 	cout << a << " , " ;
//       }
//       cout << endl;
//       cout << "before returning let us check the distances ... " << endl;
//       for(unsigned int i = 1; i <= no_of_vertices; i++)
//       {
// 	cout << i << ": " << to_double(dist[i]) << endl;
//       }
//       #endif
//
//       // cout << "returns .. " << d << endl;
//
//       return d;
//
//     }

//     /** \brief Creates a BFS tree with a random root
//      *
//      * @param ST The Spanning Tree
//      *
//      */
//     template<typename IntegerType, typename RationalType>
//     void
//     Graph<IntegerType, RationalType>::create_BFS_tree_with_random_root(SpanningTree<Graph, IntegerType, RationalType> &ST)
//     {
//       tree_edges.clear();
//       non_tree_edges.clear();
//       Random rg;
//       RationalType node_rand = rg.rng();
//       node x0 = node_rand * no_of_vertices;
//
//       vector<int> visited( no_of_vertices+1,1 );
//       for(unsigned int v = 1; v <= no_of_vertices; v++)
//       {
// 	visited[v] = 0;
//       }
//       vector<unsigned int> bfstree;
//       bfstree.push_back(0);
//       deque<node> order;
//       deque<node> Q;
//
//       Q.push_back( x0 );
//       visited[x0] = 1;
//
//       while( !Q.empty() )
//       {
// 	const node v = Q.front();
//
// 	for( auto a : incident_edges[v] ) {
// 	  node w;
// 	  if( a > 0 ) {
// 	    assert( tails[a] == v );
// 	    w = heads[a];
// 	  } else {
// 	    assert( heads[-a] == v );
// 	    w = tails[-a];
// 	  }
// 	  if( !visited[w] ) {
// 	    Q.push_back( w );
// 	    visited[w] = 1;
// 	    bfstree.push_back(abs(a));
// 	  }
// 	}
// 	Q.pop_front();
// 	order.push_front( v );
//       }
//
//       make_tree(ST , bfstree);
//       for(unsigned int i=1; i <=no_of_edges; i++){
// 	bool flag = false;
//
// 	for(auto a: bfstree){
// 	  if(a == i) {
// 	    flag = true;
// 	  }
// 	}
//
// 	if(flag == false){
// 	  non_tree_edges.push_back(i);
// 	}
//
//       }
//     }

//     /** \brief Creates a shortest path tree
//      *
//      * @param ST The Spanning Tree
//      *
//      */
//     template<typename IntegerType, typename RationalType>
//     void
//     Graph<IntegerType, RationalType>::create_SPT(SpanningTree<Graph, IntegerType, RationalType> &ST)
//     {
//       tree_edges.clear();
//       non_tree_edges.clear();
//       #ifdef VERBOSE
//       cout << "get_shortest_path_spanning_tree() " << endl;
//       #endif
//
//       Random rg;
//       RationalType node_rand = rg.rng();
//       node x0 = node_rand * no_of_vertices;
//
//       vector<RationalType> dist = dijkstra(x0);
//
//       vector<unsigned int> spt_tree;
//       spt_tree.push_back(0);
//
//       for(node v = 1; v <= no_of_vertices; v++)
//       {
// 	{
// 	  #ifdef VERBOSE
// 	  cout << "v : " << v << endl;
// 	  #endif
//
// 	  for(arc a: incident_edges[v])
// 	  {
// 	    node u = 0;
// 	    if(a > 0)
// 	    {
// 	      u = heads[a];
// 	    }
// 	    else{
// 	      u = tails[-a];
// 	    }
// 	    #ifdef VERBOSE
// 	    cout << "u: " << u  << endl  ;
//
// 	    assert(u != 0);
// 	    assert(v != u);
// 	    #endif
//
// 	    #ifdef VERBOSE
// 	    cout << x0 << " " << u << " " << v << endl;
// 	    cout << a << " : " << dist[u] << " + " << to_double(resistances[abs(a)]) << " = " << dist[v] << endl;
// 	    #endif
// 	    if( dist[u] + to_double(resistances[abs(a)]) == dist[v])
// 	    {
// 	      spt_tree.push_back(abs(a));
// 	      break;
// 	    }
// 	  }
//
// 	}
//       }
//
//       make_tree(ST , spt_tree);
//       for(unsigned int i=1; i <=no_of_edges; i++)
//       {
// 	bool flag = false;
//
// 	for(auto a: spt_tree)
// 	{
// 	  if(a == i)
// 	  {
// 	    flag = true;
// 	  }
// 	}
//
// 	if(flag == false)
// 	{
// 	  non_tree_edges.push_back(i);
// 	}
//
//       }
//
//     }
//
//
//
//
//     template<typename IntegerType, typename RationalType>
//     vector<unsigned int>
//     Graph<IntegerType, RationalType>::create_BFS_tree(vector<node>& X, node x0)
//     {
//       #ifdef VERBOSE
//       cout << "get_shortest_path_spanning_tree()   x0 = " << x0 << endl;
//       #endif
//
//       vector<RationalType> dist = dijkstra(x0, X);
//
//       #ifdef VERBOSE
//       cout << "dijkstra done " << endl;
//       #endif
//       vector<node> parents_current(no_of_vertices + 1, 0);
//
//       #ifdef VERBOSE
//       cout << "distances ... " << endl << endl;
//       for(node v = 1; v <= no_of_vertices; v++)
//       {
// 	cout << v << " : " << dist[v] << endl;
//       }
//
//       cout << endl << endl;
//       #endif
//
//       vector<unsigned int> spt_tree;
//
//       for(auto v: X){
// 	//if(cluster[v] == cluster_number)
// 	{
// 	  #ifdef VERBOSE
// 	  cout << "v : " << v << endl;
// 	  #endif
// 	  for(arc a: incident_edges[v]){
// 	    //node u = get_the_other_end_of_arc(v, a);
// 	    node u = 0;
// 	    if(a > 0){
// 	      u = heads[a];
// 	    }
// 	    else{
// 	      u = tails[-a];
// 	    }
// 	    #ifdef VERBOSE
// 	    cout << "u: " << u  << endl  ;
// 	    assert(u != 0);
// 	    assert(v != u);
// 	    #endif
// 	    //if(cluster[u] != cluster_number) continue;
//
// 	    bool is_there_in_X = false;
// 	    for(auto v: X){
// 	      if(u == v) is_there_in_X = true;
// 	    }
// 	    if(is_there_in_X == false) continue;
//
// 	    #ifdef VERBOSE
// 	    cout << a << " : " << dist[u] << " + " << resistances[abs(a)] << " = " << dist[v] << endl;
// 	    #endif
// 	    if(dist[u] + weights[abs(a)] == dist[v])/* || (dist[v] + resistances[abs(a)] == dist[u])*/
// 	    {
// 	      #ifdef VERBOSE
// 	      cout << a << " : " << dist[u] << " + " << weights[abs(a)] << " = " << dist[v] << endl;
// 	      cout << "tree_edges push - " << a << endl;
// 	      #endif
// 	      spt_tree.push_back(abs(a));
// 	      parents_current[v] = u;
// 	      #ifdef VERBOSE
// 	      cout << "parents_current[" << v << "] = " << u << endl;
// 	      #endif
// 	      // cout << v << " : : tree edge added... " << abs(a) << endl;
// 	      break;
// 	    }
// 	  }
//
//
// 	}
//       }
//
//       #ifdef VERBOSE
//       cout << "The SPT edges : " << endl;
//
//       for(auto a: spt_tree)
//       {
// 	cout << a << " , " ;
//       }
//
//       cout << endl;
//       #endif
//       return spt_tree;
//
//     }
//
//
//     template<typename IntegerType, typename RationalType>
//     vector<arc>
//     Graph<IntegerType, RationalType> :: del(vector<node>& X, vector<node>& B)
//     {
//       #ifdef VERBOSE
//       cout << "del( - ) " << endl;
//       #endif
//       vector<arc> del_X;
//
//
//
//
//
//       for(node b: B)
//       {
// 	for(auto a: incident_edges[b])
// 	{
// 	  node u = 0;
// 	  if(head(a) == b)
// 	  {
// 	    u = tail(a);
// 	  }
// 	  else{
// 	    u = head(a);
// 	  }
//
// 	  //      bool b_in_X = false;
// 	  bool u_in_X = false;
// 	  for(auto x: B)
// 	  {
// 	    if(x == b)
// 	    {
// 	      //	  b_in_X = true;
// 	    }
// 	    if(x == u)
// 	    {
// 	      u_in_X = true;
// 	    }
// 	  }
//
// 	  if(u_in_X == false)
// 	    //if(b_in_X == false && u_in_X == true)
// 	  {
// 	    del_X.push_back(a);
// 	  }
//
// 	}
//       }
//
//       return del_X;
//     }
//
//     template<typename IntegerType, typename RationalType>
//     RationalType
//     Graph<IntegerType, RationalType> :: cost(vector<arc> X)
//     {
//       #ifdef VERBOSE
//       cout << "cost() " << endl;
//       #endif
//       RationalType total_cost(0);
//       for(auto a: X)
//       {
// 	#ifdef VERBOSE
// 	cout << a << " - " << RationalType(1)/weights[abs(a)] << endl;
// 	#endif
// 	total_cost += RationalType(1)/weights[abs(a)];
//       }
//       #ifdef VERBOSE
//       cout << "total_cost: " << total_cost << endl;
//       #endif
//       return total_cost;
//     }


template<typename IntegerType, typename RationalType>
RationalType
Graph<IntegerType, RationalType>:: size_del(vector<node>& X)
{
#ifdef VERBOSE
	cout << "size_del() " << endl;
#endif

	unsigned int size_of_del = 0;

	for (unsigned int a = 1; a <= no_of_edges; a++) {
		node v = head(a);
		node u = tail(a);

		bool v_in_X = false;
		bool u_in_X = false;

		for (node x : X) {
			if (v == x) {
				v_in_X = true;
			}

			if (u == x) {
				u_in_X = true;
			}
		}

		if ( (u_in_X == true && v_in_X == false) || (u_in_X == false && v_in_X == true)) {
			size_of_del++;
		}

	}

#ifdef VERBOSE
	cout << "size_of_del: " << size_of_del << endl;
#endif

	return size_of_del;
}


template<typename IntegerType, typename RationalType>
RationalType
Graph<IntegerType, RationalType>:: size_E_plus( vector<node>& X)
{
	//   cout << "size_E_plus() " << endl;

	unsigned int size_of_E_plus = 0;

	for (unsigned int a = 1; a <= no_of_edges; a++) {
		node v = head(a);
		node u = tail(a);

		bool v_in_X = false;
		bool u_in_X = false;

		for (node x : X)
		{
			if (v == x)
			{
				v_in_X = true;
			}

			if (u == x)
			{
				u_in_X = true;
			}
		}

		if ( u_in_X == true || v_in_X == true )
		{
			size_of_E_plus ++;
		}

	}
	return size_of_E_plus;

}


//     template<typename IntegerType, typename RationalType>
//     RationalType
//     Graph<IntegerType, RationalType>::vol_x(vector<node> &X,
// 					    vector<node>& Ball,
// 					    node x0,
// 					    RationalType r
//     )
//     {
//       #ifdef VERBOSE
//       cout << "vol_x() " << endl;
//       #endif
//       RationalType alpha_sum(0);
//
//       #ifdef VERBOSE
//       cout << "r: " << r << endl;
//       #endif
//
//       //   }
//
//       vector<arc> del_X = del(X, Ball);
//
//       //cout << " ,,,, " << endl;
//
//       vector<RationalType> dist = dijkstra(x0, X);
//
//       #ifdef VERBOSE
//       cout << "blubb .. " << endl;
//
//       cout << "del_X: " << endl;
//       for(auto a : del_X)
//       {
// 	cout << a << " , " ;
//       }
//       cout << endl;
//       #endif
//
//
//       for(auto a: del_X)
//       {
// 	#ifdef VERBOSE
// 	cout << a << "  weight = " << weights[abs(a)] << endl;
// 	#endif
// 	node v = head(a);
// 	node u = tail(a);
//
// 	#ifdef VERBOSE
// 	cout << "v = " << v << endl << "u = " << u << endl;
// 	#endif
//
// 	if(r < dist[v])
// 	{
// 	  #ifdef VERBOSE
// 	  cout << endl <<  r << "   " << dist[u] << "  " << weights[abs(a)] << endl;
// 	  cout << "dist[v] = " << dist[v] << endl;
// 	  #endif
//
// 	  alpha_sum += (r - dist[u] )/weights[abs(a)];
// 	}
// 	else
// 	{
// 	  #ifdef VERBOSE
// 	  cout << endl <<  r << "   " << dist[v] << "  " << weights[abs(a)] << endl;
// 	  cout << "dist[u] = " << dist[u] << endl;
// 	  #endif
//
// 	  alpha_sum += (r - dist[v])/weights[abs(a)];
// 	}
// 	#ifdef VERBOSE
// 	cout << "jjjj" << endl;
// 	#endif
//
//       }
//
//       alpha_sum = 0;
//       alpha_sum += size_E(Ball);
//       alpha_sum += 1;
//       return alpha_sum;
//       }

template<typename IntegerType, typename RationalType>
RationalType
Graph<IntegerType, RationalType>:: size_E( vector<node> &X)
{

	unsigned int size_of_E = 0;

	for (unsigned int a = 1; a <= no_of_edges; a++)
	{

		node v = head(a);

		node u = tail(a);

		bool v_in_X = false;
		bool u_in_X = false;

		for (node x : X)
		{

			if (v == x)
			{
				v_in_X = true;
			}

			if (u == x)
			{
				u_in_X = true;
			}
		}

		if ( u_in_X == true && v_in_X == true )
		{
			size_of_E ++;
		}

	}
	return size_of_E;

}

//       template<typename IntegerType, typename RationalType>
//       vector<RationalType>
//       Graph<IntegerType, RationalType>::dijkstra_clusters(node source,
// 							  unsigned int cluster_number,
// 							  unsigned int current_no_of_vertices,
// 							  vector<vector<arc>>& incident_edges_current,
// 							  vector<unsigned int>& head_current,
// 							  vector<unsigned int>& tail_current
//       )
//       {
// 	#ifdef VERBOSE
// 	cout << "Resistance initially: " << endl ;
// 	for(unsigned int i = 1; i <= no_of_edges; i++)
// 	{
// 	  cout << i << " : "<< resistances[i] << endl;
// 	}
// 	#endif
//
// 	#ifdef VERBOSE
// 	cout << cluster[source] << " == " << cluster_number << endl;
// 	#endif
// 	assert(cluster[source] == cluster_number);
//
// 	vector<RationalType> dist(current_no_of_vertices + 1, numeric_limits<long long>::max());
// 	vector<bool> visited(current_no_of_vertices + 1, false);
//
// 	dist[source] = 0;
//
// 	priority_queue<node_distance<IntegerType, RationalType>, vector<node_distance<IntegerType, RationalType> >, compare_distance<IntegerType, RationalType> > Q;
//
// 	node_distance<IntegerType, RationalType> first = {source , 0};
// 	Q.push(first);
//
//
// 	while (!Q.empty())
// 	{
// 	  node_distance<IntegerType, RationalType> temp = Q.top();
// 	  Q.pop();
// 	  node u = temp.u;
//
// 	  for(arc a : incident_edges_current[u])
// 	  {
// 	    node v = 0;
// 	    if(a > 0 )
// 	    {
// 	      v = head_current[a];
// 	    }
// 	    else{
// 	      v = tail_current[-a];
// 	    }
//
// 	    if(cluster[v] != cluster_number) continue;
//
// 	    if( dist[u] + to_double(resistances[abs(a)]) < dist[v])
// 	    {
// 	      dist[v] =  dist[u] + to_double(resistances[abs(a)]);
// 	      node_distance<IntegerType, RationalType> new_node;
// 	      new_node.u = v;
// 	      new_node.shortest_distance = dist[v];
// 	      Q.push(new_node);
//
// 	    }
// 	  }
// 	}
//
// 	#ifdef VERBOSE
// 	cout << endl << "checking... " << endl << endl;
// 	for(unsigned int i =1; i<=current_no_of_vertices; i++)
// 	{
// 	  cout << dist[i] << endl;
// 	}
// 	#endif
//
// 	return dist;
//       }
//
//       template<typename IntegerType, typename RationalType>
//       void Graph<IntegerType, RationalType> :: clear_cluster_details()
//       {
// 	for(node v = 1; v <= no_of_vertices; v++)
// 	{
// 	  cluster[v] = 0;
// 	  nodes_in_cluster[v].clear();
// 	}
//       }



//       template<typename IntegerType, typename RationalType>
//       void Graph<IntegerType, RationalType> :: get_shortest_path_spanning_tree(unsigned int cluster_number, node root, vector<unsigned int>& tree_edges_alon,
// 									       vector<vector<arc>>& incident_edges_current, unsigned int current_no_of_vertices,
// 									       vector<unsigned int>& head_current, vector<unsigned int>& tail_current)
//       {
//
// 	#ifdef VERBOSE
// 	cout << "printing the graph ... " << endl;
//
// 	for(node v=1; v <= current_no_of_vertices; v++)
// 	{
// 	  cout << v << " : " ;
// 	  for(arc a : incident_edges_current[v])
// 	  {
// 	    cout << a << " , ";
// 	  }
//
// 	  cout << endl;
// 	}
//
// 	cout << endl << endl;
//
// 	cout << "the resistances" << endl;
//
// 	for(unsigned int i =1; i <=no_of_edges; i++)
// 	{
// 	  cout << i << " : " << resistances[i] << endl;
// 	}
// 	#endif
//
// 	#ifdef VERBOSE
// 	cout << "get_shortest_path_spanning_tree() " << endl;
// 	#endif
// 	vector<RationalType> dist = dijkstra_clusters(root, cluster_number, current_no_of_vertices, incident_edges_current, head_current, tail_current);
//
// 	#ifdef VERBOSE
// 	cout << "dijkstra done " << endl;
// 	#endif
// 	vector<node> parents_current(current_no_of_vertices + 1, 0);
//
// 	#ifdef VERBOSE
// 	cout << "distances ... " << endl << endl;
// 	for(node v = 1; v <= no_of_vertices; v++)
// 	{
// 	  cout << v << " : " << dist[v] << endl;
// 	}
//
// 	cout << endl << endl;
// 	#endif
// 	assert(cluster[root] == cluster_number);
// 	// cout << "current_no_of_vertices: " << current_no_of_vertices << endl << endl;
//
// 	#ifdef VERBOSE
// 	cout << "current_no_of_vertices = " << current_no_of_vertices << endl;
//
// 	cout << "cluster_number = " << cluster_number << endl;
// 	#endif
//
// 	for(node v = 1; v <= current_no_of_vertices; v++)
// 	{
// 	  if(cluster[v] == cluster_number)
// 	  {
// 	    #ifdef VERBOSE
// 	    cout << "v : " << v << endl;
// 	    #endif
// 	    for(arc a: incident_edges_current[v])
// 	    {
// 	      node u = 0;
// 	      if(a > 0)
// 	      {
// 		u = head_current[a];
// 	      }
// 	      else{
// 		u = tail_current[-a];
// 	      }
// 	      #ifdef VERBOSE
// 	      cout << "u: " << u  << endl ;
//
// 	      assert(u != 0);
// 	      assert(v != u);
// 	      #endif
// 	      if(cluster[u] != cluster_number) continue;
//
// 	      #ifdef VERBOSE
// 	      cout << a << " : " << dist[u] << " + " << resistances[abs(a)] << " = " << dist[v] << endl;
// 	      #endif
// 	      if((dist[u] + to_double(resistances[abs(a)]) == dist[v])/* || (dist[v] + resistances[abs(a)] == dist[u])*/)
// 	      {
// 		#ifdef VERBOSE
// 		cout << "tree_edges push - " << a << endl;
// 		#endif
// 		tree_edges_alon.push_back(abs(a));
// 		parents_current[v] = u;
// 		#ifdef VERBOSE
// 		cout << "parents_current[" << v << "] = " << u << endl;
// 		#endif
// 		break;
// 	      }
// 	    }
//
//
// 	  }
// 	}
//
// 	#ifdef VERBOSE
// 	cout << "Parents ... " << endl << endl;
// 	for(node v=1; v <= no_of_vertices; v++)
// 	{
// 	  cout << v << " : " << parents_current[v] << endl;
// 	}
// 	#endif
//       }
//
//       template<typename IntegerType, typename RationalType>
//       void Graph<IntegerType, RationalType>::form_spanning_tree_alon(SpanningTree<Graph, IntegerType, RationalType>& ST, RationalType x)
//       {
// 	#ifdef VERBOSE
// 	cout << "form_spanning_tree_alon() " << endl;
// 	#endif
// 	vector<vector<arc> > incident_edges_current(no_of_vertices + 1);
// 	vector<unsigned int> tree_edges_alon (no_of_vertices);
// 	vector<unsigned int> head_current(no_of_edges + 1);
// 	vector<unsigned int> tail_current(no_of_edges + 1);
// 	tree_edges_alon.clear();
// 	tree_edges.clear();
// 	non_tree_edges.clear();
// 	clear_cluster_details();
// 	tree_edges_alon.push_back(0);
//
// 	for(unsigned int v=1; v<= no_of_vertices; v++)
// 	{
// 	  incident_edges_current[v] = incident_edges[v];
// 	}
//
// 	for(unsigned int a = 1; a <= no_of_edges; a++)
// 	{
// 	  head_current[a] = heads[a];
// 	  tail_current[a] = tails[a];
// 	}
//
// 	unsigned int current_no_of_vertices = no_of_vertices;
//
// 	unsigned int iteration_number = 1;
//
// 	while(tree_edges_alon.size() != no_of_vertices)
// 	{
//
// 	  #ifdef VERBOSE
// 	  cout << endl <<  endl << "==========================================================================" << endl << endl;
//
// 	  for(unsigned int a = 1; a <= no_of_edges; a++)
// 	  {
// 	    cout << a << " = ( " << heads[a] << " , " << tails[a] << " )  [ "<< cluster[heads[a]] << " , " << cluster[tails[a]] << " ] " << endl;
// 	  }
// 	  #endif
//
// 	  #ifdef VERBOSE
// 	  cout << "tree-edges size: " << tree_edges_alon.size() << " " << no_of_vertices << endl << endl;
// 	  #endif
// 	  unsigned int no_of_clusters = get_clusters(current_no_of_vertices, incident_edges_current, head_current, tail_current, iteration_number, x);
// 	  #ifdef VERBOSE
// 	  cout << "Number of Clusters.. " << no_of_clusters << endl;
// 	  #endif
//
// 	  #ifdef VERBOSE
// 	  cout << endl << endl << "Current Clusters ... " << endl << endl;
// 	  #endif
//
// 	  #ifdef VERBOSE
// 	  for(node v = 1; v <= current_no_of_vertices; v++)
// 	  {
// 	    cout << v << "  : " << cluster[v] << endl;
// 	  }
//
//
// 	  cout << endl << "checking nodes in cluster.. " << endl << endl;
// 	  for(unsigned int i = 1; i <= no_of_clusters; i++)
// 	  {
// 	    cout << i << " :  " ;
// 	    for(node v : nodes_in_cluster[i])
// 	    {
// 	      cout << v << " , " ;
// 	    }
// 	    cout << endl;
// 	  }
// 	  #endif
//
//
// 	  for(unsigned int i=1; i<= no_of_clusters; i++)
// 	  {
// 	    node root = nodes_in_cluster[i].front();
// 	    assert(cluster[root] == i);
// 	    get_shortest_path_spanning_tree(i, root,tree_edges_alon, incident_edges_current, current_no_of_vertices, head_current, tail_current);
// 	  }
//
// 	  incident_edges_current = contract_clusters(no_of_clusters, incident_edges_current, head_current, tail_current);
//
// 	  iteration_number++;
//
// 	  current_no_of_vertices = no_of_clusters;
//
// 	  #ifdef VERBOSE
// 	  cout << " Incident Edges after the verticies are contracted ... " << current_no_of_vertices << "  " << endl << endl;
//
// 	  for(node v = 1; v <= current_no_of_vertices; v++)
// 	  {
// 	    cout << v << " : " ;
// 	    for(auto a : incident_edges_current[v])
// 	    {
// 	      cout << a << " , " ;
// 	    }
// 	    cout << endl;
// 	  }
// 	  #endif
// 	}
//
// 	make_tree(ST, tree_edges_alon);
//
// 	for(unsigned int i=1; i <=no_of_edges; i++)
// 	{
// 	  bool flag = false;
//
// 	  for(auto a: tree_edges_alon)
// 	  {
// 	    if(a == i)
// 	    {
// 	      flag = true;
// 	    }
// 	  }
//
// 	  if(flag == false)
// 	  {
// 	    non_tree_edges.push_back(i);
// 	  }
//
// 	}
//
// 	#ifdef VERBOSE
//
// 	cout << "Non Tree Edges: " << non_tree_edges.size() << endl;
// 	cout << "non tree_edges" << endl;
// 	for(unsigned int i = 0; i < non_tree_edges.size() ; i++)
// 	{
// 	  cout << non_tree_edges[i] << endl;
// 	}
//
// 	cout << endl << endl;
//
// 	cout << "Tree Edges: " << tree_edges.size() << endl;
// 	cout << "Alon Edges: " << tree_edges_alon.size() << endl;
//
//
// 	cout << "The Tree Edges ... " << endl << endl;
//
// 	for(unsigned int i = 0; i < tree_edges_alon.size(); i++)
// 	{
// 	  cout << tree_edges[i] << " " << tree_edges_alon[i+1] << endl;
// 	}
// 	#endif
//
//       }

//       template<typename IntegerType, typename RationalType>
//       int Graph<IntegerType, RationalType>::get_clusters(unsigned int current_no_of_vertices, vector< vector<arc>>& incident_edges_current, vector<unsigned int>& head_current,
// 							 vector<unsigned int>& tail_current, unsigned int iteration_number, RationalType x)
//
//       {
// 	clear_cluster_details();
//
// 	for(node v=1; v< no_of_vertices; v++)
// 	{
// 	  assert(cluster[v] == 0);
// 	  assert(nodes_in_cluster[v].size() == 0) ;
// 	}
// 	int cluster_number = 0 ;
//
// 	bool clustering_not_complete = true;
//
// 	while(clustering_not_complete)
// 	{
// 	  cout << "while clustering not complete ... " << endl << endl;
//
// 	  node root  = 0;
// 	  clustering_not_complete = false;
// 	  for(node v = 1; v <=current_no_of_vertices; v++)
// 	  {
// 	    if(cluster[v] == 0 )
// 	    {
// 	      clustering_not_complete = true;
// 	      root = v;
// 	      break;
// 	    }
// 	  }
//
// 	  if(clustering_not_complete == false) break;
//
// 	  cluster_number++;
// 	  #ifdef VERBOSE
// 	  cout << "Cluster Numbers: " << cluster_number << "   clustering_not_complete:" << clustering_not_complete << endl;
// 	  cout << "number of verticies... " << no_of_vertices << endl;
// 	  cout << "get_layers_in_spanning_tree() would be called,  root: " << root << " , " << cluster[root] << endl << endl;
// 	  #endif
// 	  assert(cluster[root] == 0);
//
// get_layers_in_spanning_tree(root,
//    cluster_number, incident_edges_current,
//    current_no_of_vertices, head_current, tail_current, iteration_number, x);
// 	  #ifdef VERBOSE
// 	  cout << "get_layers_in_spanning_tree() has been called" << endl << endl;
// 	  #endif
// 	}
//
// 	return cluster_number;
//       }
//
//
//       template<typename IntegerType, typename RationalType>
// vector< vector<arc> > Graph<IntegerType,
// RationalType>::contract_clusters(unsigned int no_of_clusters,
//   vector<vector<arc>>& incident_edges_spanning_tree,
//   vector<unsigned int>& head_current, vector<unsigned int>& tail_current)
//       {
//
// 	cout << endl << "contract_clusters() " << endl;
//
// 	vector<vector<arc> > incident_edges_spanning_tree_temp(no_of_vertices + 1);
// 	unsigned int contracted_verticies = 0;
// 	for(unsigned int i =1 ; i <= no_of_clusters; i ++)
// 	{
// 	  contracted_verticies++;
// 	  for(node v = 1; v <= no_of_vertices; v++)
// 	  {
// 	    if(cluster[v] == i )
// 	    {
// 	      for(arc a: incident_edges_spanning_tree[v])
// 	      {
// 		if(cluster[head(a)] != cluster[tail(a)])
// 		{
// 		  incident_edges_spanning_tree_temp[contracted_verticies].push_back(a);
//
// 		  if(a > 0)
// 		  {
// 		    tail_current[a] = contracted_verticies;
// 		  }
// 		  if(a < 0)
// 		  {
// 		    head_current[-a] = contracted_verticies;
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//
// 	return incident_edges_spanning_tree_temp;
//       }
//
//       template<typename IntegerType, typename RationalType>
//       void
//       Graph<IntegerType, RationalType>::flip_arc(
// 	arc arc)
//       {
// 	int head = heads[arc];
// 	heads[arc] = tails[arc];
// 	tails[arc] = head;
//
// 	head = heads[arc];
// 	int tail = tails[arc];
// 	flows[arc]= -flows[arc];
// 	for(unsigned int i=0; i<incident_edges[head].size(); i++)
// 	{
// 	  if(arc==abs(incident_edges[head][i]))
// 	  {
// 	    incident_edges[head][i] = incident_edges[head][i]*(-1);
// 	  }
// 	}
//
// 	for(unsigned int i=0; i<incident_edges[tail].size(); i++)
// 	{
// 	  if(arc==abs(incident_edges[tail][i]))
// 	  {
// 	    incident_edges[tail][i] = incident_edges[tail][i]*(-1);
// 	  }
// 	}
//
//       }



template<typename IntegerType, typename RationalType>

node Graph<IntegerType, RationalType> :: get_l_star( vector <vector<int>> & countE , unsigned int j, unsigned int no_of_levels, RationalType x)
{

	cout << "get_l_star() " << endl;

	cout << "printing counE ... " << endl;

	for (unsigned int i = 1; i <= j; i++)
	{
		cout << i << ": ";
		for (unsigned int l = 1; l <= no_of_levels; l++)
		{
			cout << countE[i][l] << " , ";
		}

		cout << endl;
	}

	cout << endl;

	node l_star = 0;

	for (unsigned int i = 1; i <= j ; i++)
	{
		int sumE_i = 0;
		node l_star_i = 0;
		for (unsigned int l = 1; l <= no_of_levels; l++)
		{

			sumE_i += countE[i][l];

			cout << "x = " << x << endl;
			cout << sumE_i << " <= " << x*countE[i][l + 1] << endl;

			if (sumE_i >= x * countE[i][l + 1])
			{
				l_star_i = l;
				break;
			}
		}

		if (l_star_i > l_star)
		{
			l_star = l_star_i;
		}
	}

	return l_star;
}



//       template<typename IntegerType, typename RationalType>
//       void Graph<IntegerType, RationalType>::get_layers_in_spanning_tree(node root, int cluster_number,
// 									 vector< vector<arc>>& incident_edges_current,
// 									 unsigned int current_no_of_vertices,
// 									 vector<unsigned int>& head_current,
// 									 vector<unsigned int>& tail_current,
// 									 unsigned int iteration_number,
// 									 RationalType x
//       )
//       {
//
// 	#ifdef VERBOSE
// 	cout << "printing the graph ... " << endl << endl;
//
// 	for(node v=1; v <= current_no_of_vertices; v++)
// 	{
// 	  cout << v << " : " ;
// 	  for(auto a: incident_edges_current[v])
// 	  {
// 	    cout << a << " , " ;
// 	  }
// 	  cout << endl;
// 	}
// 	#endif
//
// 	#ifdef VERBOSE
// 	cout << endl << "get_layers_in_spanning_tree()" << endl;
// 	#endif
//
// 	#ifdef VERBOSE
// 	cout << "Root : " << root << endl;
// 	#endif
//
// 	cluster[root] = cluster_number;
// 	nodes_in_cluster[cluster_number].push_back(root);
//
// 	#ifdef VERBOSE
// 	cout << "nodes_in_cluster[" << cluster_number << "]  pushed:  " << root << "size:  " << nodes_in_cluster[cluster_number].size() << endl;
// 	#endif
//
// 	for(unsigned int i =1; i < nodes_in_cluster[cluster_number].size(); i++){
// 	  cout << nodes_in_cluster[cluster_number][i] << " , " ;
// 	}
// 	cout << endl << endl;
//
// 	vector<int> visited(current_no_of_vertices + 1, 0);
// 	vector<arc> bfstree( current_no_of_vertices + 1 , 0);
// 	vector<int> level (current_no_of_vertices +1,  numeric_limits<int>::max());
// 	vector<vector<arc> > E(current_no_of_vertices + 1);
//
//
// 	vector<vector<int> > countE( no_of_vertices + 1, vector<int>( no_of_vertices + 1) );
//
// 	for(unsigned int i = 0; i <=no_of_vertices; i++){
// 	  for(unsigned int j = 0; j <= no_of_vertices; j++){
// 	    countE[i][j] = 0;
// 	  }
// 	}
//
// 	vector<vector<node> > V(current_no_of_vertices + 1);
// 	vector<int> mark_edges(no_of_edges + 1, 0);
//
//
// 	int rho = ceil(3*log(current_no_of_vertices)/log(x));
// 	RationalType mu = 9*rho*log(current_no_of_vertices);
// 	RationalType y = x*mu;
//
// 	cout << "x: " << x << endl;
// 	cout << "rho: " << rho << endl;
// 	cout << "mu: " << mu << endl;
// 	cout << "y: " << y << endl;
//
// 	#ifdef VERBOSE
//
// 	cout << "lognloglogn : " << log(current_no_of_vertices)*log(log(current_no_of_vertices)) << endl;
// 	cout << "current_no_of_vertices: " << current_no_of_vertices << endl;
//
// 	#endif
//
// 	level[root] = 0;
// 	deque<node> order;
// 	deque<node> Q;
// 	Q.push_back(root);
// 	visited[root] = 1;
// 	unsigned int no_of_levels = 0;
// 	#ifdef VERBOSE
// 	cout << "while loop will begin" << endl;
//
// 	cout << "root: " << root << endl;
// 	#endif
// 	while(!Q.empty()){
// 	  const node v = Q.front();
//
// 	  #ifdef VERBOSE
// 	  cout << "vertex v: " << v << endl;
// 	  #endif
//
// 	  for( auto a : incident_edges_current[v] ) {
//
// 	    #ifdef VERBOSE
// 	    cout << "vertex v : " << v << endl;
// 	    cout << "arc a:  " << a << endl;
// 	    #endif
// 	    node w;
// 	    if( a > 0 ) {
// 	      assert( tail_current[a] == v );
// 	      w = head_current[a];
// 	    } else {
// 	      assert( head_current[-a] == v );
// 	      w = tail_current[-a];
// 	    }
//
// 	    #ifdef VERBOSE
// 	    cout << "vertex w : " << w << "  " << visited[w] << " " << cluster[w] << endl;
// 	    #endif
//
// 	    if(cluster[w] != 0) continue;
// 	    if( !visited[w] ) {
// 	      Q.push_back( w );
// 	      visited[w] = 1;
// 	      bfstree[w] = a;
// 	      level[w] = level[v] + 1;
// 	      V[level[w]].push_back(w);
// 	    }
//
// 	    if(visited[w]){
// 	      if(mark_edges[abs(a)] == 0){
// 		E[level[w]].push_back(a);
// 		#ifdef VERBOSE
// 		cout << "Pushing into E[" << level[w] << "]:  " << a << endl;
// 		#endif
//
//
// 		for(unsigned int i = 1; i <= iteration_number; i++){
//
// 		  #ifdef VERBOSE
// 		  cout << "checking the resistance range ... " << endl;
// 		  cout << "y: " << y << endl;
// 		  cout << i << ":   arc: " << a << "  " << pow(y ,i-1) << " < " << resistances[abs(a)] << " < " << pow(y, i) << endl;
// 		  #endif
// 		  long long pow_y_i_minus_1A = pow(y,i-1);
// 		  long long pow_y_iA = pow(y, i);
//
// 		  unbounded_integer<long long> pow_y_i_minus_1= pow_y_i_minus_1A;
// 		  unbounded_integer<long long> pow_y_i = pow_y_iA;
//
// 		  if( (( pow_y_i_minus_1 < resistances[abs(a)] )
// || ( pow_y_i_minus_1 == resistances[abs(a)] )) && ( ( resistances[abs(a)] < pow_y_i )
// || ( resistances[abs(a)] ==  pow_y_i ) ) )
// 		  {
// 		    countE[i][level[w]]++;
// 		    #ifdef VERBOSE
// 		    cout << " ( "  << i << ", " << level[w] << " ) :  " << countE[i][level[w]] << endl;
// 		    #endif
// 		  }
// 		}
//
// 		mark_edges[abs(a)] = 1;
// 	      }
// 	      no_of_levels = level[w];
// 	    }
//
// 	  }
//
// 	  Q.pop_front();
//
// 	  order.push_front( v );
// 	}
//
//
// 	cout << "printing countE: " <<  endl;
//
// 	for(unsigned int i=1; i <= iteration_number; i++)
// 	{
// 	  cout << i << ": ";
// 	  for(unsigned int l = 1; l <= no_of_levels; l++)
// 	  {
// 	    cout << countE[i][l] << " , ";
// 	  }
//
// 	  cout << endl;
// 	}
//
//
//
// 	no_of_levels ++;
// 	#ifdef VERBOSE
// 	cout << "Number of Levels: " << no_of_levels << endl << endl;
// 	cout << "while loop ended ... " << endl;
// 	#endif
// 	node l_star = 0;
//
//
// 	RationalType sum_E = 0.0;
// 	bool no_l_star = true;
//
// 	#ifdef VERBOSE
// 	cout << endl << "checking the 'E' ... " << endl;
// 	int sum = 0;
// 	for(unsigned int i = 1; i <= current_no_of_vertices; i++)
// 	{
// 	  sum += E[i].size();
// 	  cout << i << " : (  size = " << E[i].size() << " ) " ;
// 	  for(auto a : E[i])
// 	  {
// 	    cout << a << " , " ;
// 	  }
// 	  cout << endl;
// 	}
//
// 	cout << "checking the 'V' ... " << endl;
// 	for(unsigned int i=1; i<= no_of_levels; i++)
// 	{
// 	  cout << i << " : " ;
// 	  for(auto a: V[i])
// 	  {
// 	    cout << a << " , ";
// 	  }
// 	  cout << endl;
// 	}
// 	#endif
//
//
// 	cout << "number of levels: " << no_of_levels << endl;
// 	for(unsigned int i = 1; i <= no_of_levels; i++)
// 	{
// 	  sum_E += E[i].size();
// 	  #ifdef VERBOSE
// 	  cout << sum_E << " >= " << x*E[i+1].size() << "  " << E[i+1].size() << endl;
// 	  #endif
// 	  if(sum_E >= x*E[i+1].size())
// 	  {
// 	    l_star = i;
// 	    no_l_star = false;
// 	    break;
// 	  }
// 	}
//
// 	cout << "current number of verticies before get_l_star() is called ... " <<  current_no_of_vertices << endl << endl;
// 	l_star = get_l_star(countE , iteration_number, no_of_levels, x);
//
// 	cout << "l_star : " << l_star << endl;
//
// 	if(l_star == 0)
// 	{
// 	  no_l_star = true;
// 	}
//
//
// 	if(no_l_star)
// 	{
// 	  for(node v =1; v <= no_of_vertices; v++)
// 	  {
// 	    if(cluster[v] == 0 )
// 	    {
// 	      cluster[v] = cluster_number;
// 	      nodes_in_cluster[cluster_number].push_back(v);
// 	    }
// 	  }
// 	}
// 	else
// 	{
// 	  for(unsigned int i = 1; i <= l_star; i++)
// 	  {
// 	    for(node v: V[i])
// 	    {
// 	      cluster[v] = cluster_number;
// 	      nodes_in_cluster[cluster_number].push_back(v);
// 	    }
// 	  }
// 	}
//
//
//
// 	#ifdef VERBOSE
//
// 	cout << endl << "l_star :  " << l_star << endl << endl;
// 	cout << endl << "checking the levels..." << endl;
// 	for(unsigned int i = 1; i <= current_no_of_vertices; i++)
// 	{
// 	  cout << i << " : " << level[i] << endl;
// 	}
//
//
// 	#endif
//       }


template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>::swap(
  RationalType &x,
  RationalType &y)
{
	RationalType temp;
	temp = x;
	x = y;
	y = temp;

}

template<typename IntegerType, typename RationalType>
void Graph<IntegerType, RationalType>::swap(int &x, int &y) {
	int temp;
	temp = x;
	x = y;
	y = temp;
}


template<typename IntegerType, typename RationalType>
node
Graph<IntegerType, RationalType>::find(
  node v)
{
	node w = v;
	while (v != parent[v])
	{
		v = parent[v];
	}

	while (w != v)
	{
		node p = parent[w];
		parent[w] = v;
		w = p;
	}

	return v;

}

/** \brief unite function
 *
 * @param a Arc a
 *
 */
template<typename IntegerType, typename RationalType>
bool
Graph<IntegerType, RationalType>::unite(
  arc a)
{
	assert( a > 0 );
	const node v = tails[a];
#ifdef VERBOSE
	cout << " v = " << v << endl;
#endif
	const node w = heads[a];
#ifdef VERBOSE
	cout << "w = " << w << endl;
#endif
	const int r = find(v);
#ifdef VERBOSE
	cout << "r = " << r << endl;
#endif
	const int s = find(w);

#ifdef VERBOSE
	cout << "s = " << s << endl;
	cout << "rank[s] = " << rank[s] << endl;
	cout << "rank[r] = " << rank[r] << endl;
#endif

	if (r == s) {
		{
			non_tree_edges.push_back( a );

#ifdef VERBOSE
			cout << endl << a << ": set already_pushed to 1 (case 1)" << endl;
#endif
		}
		return false;
	}
	else {
		if (rank[r] > rank[s]) {
			parent[s] = r;

#if VERBOSE
			cout << "arcs_to_root of " << s << " is " << arcs_to_root[s] << endl;
			cout << -a << " is now a tree edge " << endl;
#endif

			assert( arcs_to_root[s] == 0 );
			arcs_to_root[s] = -a;
			tree_edges.push_back( -a );
		}
		else {
			parent[r] = s;

#if VERBOSE
			cout << "arcs_to_root of " << r << " is " << arcs_to_root[r] << endl;
			cout << endl << a << " is now a tree edge " << endl;
#endif

			assert( arcs_to_root[r] == 0 );
			arcs_to_root[r] = a;
			tree_edges.push_back( a );

			if (rank[r] == rank[s]) {
				rank[s] += 1;
			}
		}
		count_tree_edges++;

#ifdef VERBOSE
		cout << count_tree_edges << ": " << a << " = (" << tail(a) << "," << head(a) << ") " << r << " " << s << endl;
#endif

		return true;
	}
}

/** \brief Creates the tree
 *
 * @param ST The Spanning Tree
 * @param permutation
 *
 */
template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>::make_tree(
  SpanningTree<Graph, IntegerType, RationalType>& ST, vector<unsigned int>& permutation)
{

	count_tree_edges = 0;
	tree_edges.clear();
	for (unsigned int i = 0; i <= no_of_vertices; i++)
	{
		arcs_to_root[i] = 0;
		ST.depth[i] = -1;
		//parents[i] = 0;
		parent[i] = i;
		ST.size_sub_tree_rooted_at[i] = 0;
		ST.children[i].clear();
	}


	for (unsigned int i = 1; i < permutation.size(); i++) {

		const arc a = permutation[i];

		unite( a );


	}

	for ( auto& a : boost::adaptors::reverse( tree_edges ) ) {
		assert( a != 0 );
		const node s = tail(a);
		const node t = head(a);
		const arc as = arcs_to_root[s];
		const arc at = arcs_to_root[t];
#ifdef VERBOSE
		cout << "tree arc " << a << " = (" << s << "," << t << "), "
		     << as << " = (" << tail(as) << "," << head(as) << "), "
		     << at << " = (" << tail(at) << "," << head(at) << ") " << endl;
#endif
		if ( as == a ) continue;
		if ( at == -a ) {
			a = -a;
			continue;
		}
		node v = head(as);
		arc uv = as;
		while ( uv != a ) {
			const arc vw = arcs_to_root[v];
#ifdef VERBOSE
			cout << v << "--> " << uv << " = (" << tail(uv) << "," << head(uv) << ") -> " << vw << " = (" << tail(vw) << "," << head(vw) << ")" << endl;
#endif
			assert( vw != 0 );
			arcs_to_root[v] = -uv;
			uv = vw;
			v = head(uv);
		}
		arcs_to_root[s] = a;
	}
	for ( auto &a : tree_edges ) {
		if ( a == -arcs_to_root[head(a)] ) {
#ifdef VERBOSE
			cout << "flip " << a << endl;
#endif
			a = -a;
		}
	}

}


//       /** \brief Creates a minimum spanning tree
//        *
//        * @param ST The Spanning Tree
//        *
//        */
//       template<typename IntegerType, typename RationalType>
//       void
//       Graph<IntegerType, RationalType>::create_MST(
// 	SpanningTree<Graph, IntegerType, RationalType>& ST)
//       {
// 	tree_edges.clear();
// 	non_tree_edges.clear();
//
// 	std::vector<unsigned int> permutation( no_of_edges+1 );
// 	for(unsigned int i=0; i<=no_of_edges; i++){
// 	  permutation[i] = i;
// 	}
// 	std::sort(permutation.begin()+1, permutation.end(),[this] (unsigned int a, unsigned int b) -> bool { return resistances[a] < resistances[b]; } );
//
// 	make_tree(ST, permutation);
//
//       }

/** \brief Run an Euler tour through the Spanning Tree
 *
 * Runs and euler tour and stores in the arcs in a DFS order. Finds and stores the sizes of the sub tree rooted at each vertex
 *
 * @param SpanningTree ST
 */
template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>::euler_tour_and_get_subtree_sizes(
  SpanningTree<Graph, IntegerType, RationalType>& ST)
{
	ST.dfs_euler_tour.clear();
	ST.dfs_euler_tour.reserve( 2 * no_of_vertices );
	vector<int> visited(no_of_vertices + 1, 0 );
	stack <node> order;
	stack<node> L;
	L.push( ST.root );
	visited[ST.root] = 1;
	ST.discovery[ST.root] = 0;
	unsigned int time = 0;

#ifdef VERBOSE
	cout << "order: " << endl << endl;
#endif

	while ( !L.empty() ) {
		const node v = L.top();

#ifdef VERBOSE
		cout << time << ": " << v << " " << visited[v] << " " << ST.arc_to_parent[v] << " = (" << tail(ST.arc_to_parent[v]) << "," << head(ST.arc_to_parent[v]) << ")" << endl;
#endif

		if ( visited[v] == 1 ) {
			order.push( v );
			visited[v] = -1;
			if ( ST.arc_to_parent[v] != 0 ) {
				ST.discovery[v] = time;
				ST.dfs_euler_tour.push_back( -ST.arc_to_parent[v] );
			}
			for ( auto w : ST.children[v]) {
				L.push( w );
				visited[w] = 1;
			}
		} else {
			ST.finished[v] = time - 1;
			if ( ST.arc_to_parent[v] != 0 ) {
				ST.dfs_euler_tour.push_back( ST.arc_to_parent[v] );
			}
			L.pop();
		}

#ifdef VERBOSE
		cout << endl;
		cout << v << "; " ;
#endif

		++time;
	}

	while (!order.empty())
	{
		node w = order.top();
		if (ST.children[w].size() == 0) {
			ST.size_sub_tree_rooted_at[w] = 0;
#ifdef VERBOSE
			cout << "Size at " << w << " : " << ST.size_sub_tree_rooted_at[w] << endl;
#endif
			order.pop();
		}
		else {
			int size = 0;
			std::vector<int> children_rooted_at;
			children_rooted_at.clear();
			for (auto c : ST.children[w]) {
				size += ST.size_sub_tree_rooted_at[c];
				children_rooted_at.push_back(c);
			}

			ST.size_sub_tree_rooted_at[w] = size + ST.children[w].size();

#ifdef VERBOSE
			cout << endl << "Size at " << w << " : " << ST.size_sub_tree_rooted_at[w] << endl;
#endif
			order.pop();
		}
	}
}

/** \brief Find depth and children of each node
 *
 * Finds and stores the depth and children of each node in the Spanning Tree
 *
 * @param SpanningTree ST
 */
template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>:: get_depth_and_children(
  SpanningTree<Graph, IntegerType, RationalType>& ST)
{
	stack<arc> S;
	for ( auto a : tree_edges ) {
		S.push( a );
	}
	while ( !S.empty() ) {
		const arc a = S.top();
		const node v = tail(a);
		const node w = head(a);
		if ( arcs_to_root[w] == 0 ) {
			//root = w;
			ST.root = w;
			ST.depth[w] = 0;
			ST.parent[w] = 0;
			ST.arc_to_parent[w] = 0;
		}
		if ( ST.depth[w] >= 0 ) {
			if (ST.depth[v] == -1)
			{
				ST.children[w].push_back(v);
				ST.parent[v] = w;
				ST.arc_to_parent[v] = a;
			}
			ST.depth[v] = ST.depth[w] + 1;
			S.pop();
		} else {
			S.push(arcs_to_root[w]);
		}
	}
}



//       template<typename IntegerType, typename RationalType>
//       void
//       Graph<IntegerType, RationalType>::form_petal_decomposition_spanning_tree(SpanningTree<Graph, IntegerType, RationalType>& ST)
//       {
// 	#ifdef VERBOSE
// 	cout << "form_petal_decomposition_spanning_tree() " << endl;
//
// 	cout << "The Graph: " << endl;
// 	for(unsigned int i = 1; i <= no_of_vertices; i++){
// 	  cout << i << " : ";
// 	  for(auto a : incident_edges[i]){
// 	    cout << a << " , ";
// 	  }
// 	  cout << endl;
// 	}
// 	cout << endl;
// 	#endif
//
// 	#ifdef VERBOSE
// 	cout << "Resistances: " << endl;
// 	for(unsigned int i = 1; i <= no_of_edges; i++){
// 	  cout << i << ": " << to_double(resistances[i]) << endl;
// 	}
// 	cout << endl;
// 	#endif
//
// 	vector<node> X (no_of_vertices, 0);
// 	for(unsigned int i = 0; i < no_of_vertices; i++){
// 	  X[i] = i+1;
// 	}
//
// 	for(unsigned int i = 1; i <= no_of_edges; i++)
// 	{
// 	  weights[i] = to_double(resistances[i]);
// 	}
//
// 	vector<unsigned int> tree_edges_petal;
// 	vector<unsigned int> tree_edges_petal1;
// 	tree_edges_petal1.push_back(0);
//
// 	Random rg;
// 	RationalType node_rand = rg.rng();
// 	node x0 = node_rand * no_of_vertices;
// 	if( x0 == 0)
// 	{
// 	  x0 ++;
// 	}
// 	#ifdef VERBOSE
// 	cout << "x0: " << x0 << endl;
// 	#endif
//
// 	tree_edges_petal = hierarchial_petal_decomposition(ST, X, x0, x0);
//
// 	for(auto a: tree_edges_petal)
// 	{
// 	  if(a != 0)
// 	  {
// 	    tree_edges_petal1.push_back(a);
// 	  }
// 	}
// 	cout << endl;
//
// 	make_tree(ST, tree_edges_petal1);
//
// 	#ifdef VERBOSE
// 	cout << "after make tree" << endl;
// 	cout << "tree_edges_petal: ";
// 	for(auto a: tree_edges_petal)
// 	{
// 	  cout << a << " , " ;
// 	}
// 	cout << endl;
// 	cout << "tree_edges: " ;
//
// 	for(auto a: tree_edges)
// 	{
// 	  cout << a << " , " ;
// 	}
// 	cout << endl;
// 	#endif
// 	non_tree_edges.clear();
// 	for(unsigned int i=1; i <=no_of_edges; i++)
// 	{
// 	  bool flag = false;
//
// 	  for(auto a: tree_edges_petal)
// 	  {
//
// 	    if(a == i)
// 	    {
// 	      flag = true;
// 	    }
// 	  }
//
// 	  if(flag == false)
// 	  {
// 	    non_tree_edges.push_back(i);
// 	  }
//
// 	}
//
// 	#ifdef VERBOSE
// 	cout << "tree edges: " << tree_edges.size() << endl;
// 	for(auto a : tree_edges)
// 	{
// 	  cout << a << " , " ;
// 	}
// 	cout << endl;
// 	cout << "non tree edges: " ;
// 	for(auto a : non_tree_edges)
// 	{
// 	  cout << a << " , ";
// 	}
// 	cout << endl;
// 	#endif
//       }


template<typename IntegerType, typename RationalType>
RationalType
Graph<IntegerType, RationalType>:: rad(node x0, vector<node>& X)
{
#ifdef VERBOSE
	cout << "rad( " << x0 << " , " << " X  ) " << endl;
#endif
	RationalType radius = 0;

	vector<RationalType> dist = dijkstra(x0, X);
	for (node v : X)
	{
#ifdef VERBOSE
		cout << "x0: " << x0 << "    v: " << v << endl;
#endif
		vector<node> previous(no_of_vertices + 1, 0);

#ifdef VERBOSE
		cout << dist[v] << " >= " << (radius) << endl;
#endif

		previous.clear();
		if ( dist[v] >= radius)
		{
			radius = dist[v];
		}
	}

	return radius;
}


template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>::add_path(node t, unsigned int l)
{

	for (unsigned int i = 1; i <= l; i++)
	{

		no_of_vertices++;

#ifdef VERBOSE
		cout << "Number of verticies: " << no_of_vertices << endl;
#endif

#ifdef VERBOSE
		cout << "new_edge between ... " << t << " and " << no_of_vertices << endl;
#endif
		incident_edges.resize(no_of_vertices + 1);
		incident_edges[no_of_vertices].clear();

#ifdef VERBOSE
		cout << "size: " << incident_edges[no_of_vertices].size() << endl;
#endif

		new_edge(t, no_of_vertices);
		t = no_of_vertices;
	}
}


template<typename IntegerType, typename RationalType>
vector<unsigned int>
Graph<IntegerType, RationalType>::hierarchial_petal_decomposition(SpanningTree<Graph, IntegerType, RationalType>&ST, vector<node>& X, node x0, node t)
{

#ifdef VERBOSE
	cout << "hierarchial_petal_decomposition called.. : " << endl;
#endif
	vector<unsigned int> tree_edges_petal;
	vector<unsigned int> tree_edges_BFS;
	RationalType radius = rad(x0, X);
	RationalType lognloglogn = log(no_of_vertices) * log(log(no_of_vertices));

#ifdef VERBOSE
	cout << "x0 : " << x0 << endl;
	cout << "radius: " << radius << endl << "10 * lognloglogn: " << 10 * lognloglogn << endl;
#endif
	if (radius <= 10 * lognloglogn)
	{
#ifdef VERBOSE
		cout << "create_BFS_tree() " << endl;
#endif
		tree_edges_BFS = create_BFS_tree(X, x0);

		return tree_edges_BFS;
	}
	else
	{
#ifdef VERBOSE
		cout << "else case: " << endl;
#endif
		vector<vector<node>> petals;
		petals.push_back(X);
		vector<node> petal_centers;
		petal_centers.push_back(x0);
		vector<node> y(no_of_vertices + 1 , 0);
		vector<node> t_(no_of_vertices + 1, 0);
		unsigned int s = petal_decomposition(X, x0, t , petals, petal_centers, y, t_);
#ifdef VERBOSE
		cout << "petal decomposition done: " << endl;
		cout << "number of petals: " << s << endl;

		cout << "petal centers: " ;
		for (auto v : petal_centers)
		{
			cout << v << " , " ;
		}
		cout << endl;

		cout << "y: " ;
		for (auto a : y)
		{
			cout << a << " , ";
		}
		cout << endl;

		cout << "t_: " ;
		for (auto a : t_)
		{
			cout << a << " , " ;
		}
		cout << endl;

		cout << "s: " << s << endl;
#endif
		vector<vector<unsigned int> > tree_edges_petal_loop(s + 1);

		for (unsigned int j = 0; j <= s; j++)
		{
#ifdef VERBOSE
			cout << "hierarchial_petal_decomposition() " << j << endl;
			cout << "X: " ;
			for (auto a : petals[j])
			{
				cout << a << " , ";
			}
			cout << endl;
#endif
			tree_edges_petal_loop[j] = hierarchial_petal_decomposition(ST, petals[j], petal_centers[j], t_[j]);

			for (auto a : tree_edges_petal_loop[j])
			{
#ifdef VERBOSE
				cout << a << " , " ;
#endif
				tree_edges_petal.push_back(a);
			}

			if (j != 0)
			{
#ifdef VERBOSE
				cout << get_edge(y[j], petal_centers[j]) << " =  edge b/w " << y[j] << " and " << petal_centers[j] << endl;
#endif
				tree_edges_petal.push_back(abs(get_edge(y[j], petal_centers[j])));
			}
			cout << endl;
		}
	}

	return tree_edges_petal;

}


template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>::create_stigma(vector<vector<node>>& petals, vector<vector<node>>& Y, node s, vector<node>& y, vector<node>& t_ )
{
#ifdef VERBOSE
	cout << "create_stigma() called ... " << endl;
#endif
	petals[0] = Y[s];
#ifdef VERBOSE
	cout << "petals[0] = " << endl;
	for (auto a : petals[0])
	{
		cout << a << " , " ;
	}
	cout << endl;
#endif
}

//       template<typename IntegerType, typename RationalType>
//       unsigned int
//       Graph<IntegerType, RationalType>::create_remaining_petals(unsigned int j, node r0, node x0, RationalType delta, vector<node>& X, vector<vector<node>>& Y,
// 								vector<vector<node>>& petals, vector<node>& petal_centers, vector<node>& y, vector<node>& t_)
//       {
// 	#ifdef VERBOSE
// 	cout << "create_remaining_petals() " << endl;
// 	#endif
// 	vector<node> nodes_in_ball;
//
// 	cout << endl;
//
// 	Ball_x(x0, 3*delta/4, X, nodes_in_ball);
//
// 	#ifdef VERBOSE
// 	cout << "Ball: " ;
// 	for(auto a : nodes_in_ball)
// 	{
// 	  cout << a << " , " ;
// 	}
// 	#endif
//
// 	while( subtract_vectors_size(Y[j-1], nodes_in_ball) != 0 )
// 	{
// 	  t_[j] = 0;
// 	  vector<node> previous(2*no_of_vertices  + 1, 0);
// 	  vector<RationalType> dist = dijkstra(x0, X);
// 	  for(auto a: Y[j-1])
// 	  {
// 	    #ifdef VERBOSE
// 	    cout << "a out of Y1 : " << a << endl;
// 	    cout << "subtract_vectors: --- " << ceil(r0) << " , " << floor(r0) << " = " << dist[a]<< endl;
// 	    cout << "r0: " << r0 << endl << "dijkstra: " << dist[a] << endl;
// 	    cout << "rad(x0, X); " << rad(x0, X) << endl;
// 	    #endif
// 	    if( dist[a] > 3*delta/4  )
// 	    {
// 	      t_[j] = a;
// 	      #ifdef VERBOSE
// 	      cout << "X: ";
// 	      for(auto a: X)
// 	      {
// 		cout << a << " , " ;
// 	      }
// 	      cout << endl;
// 	      cout << "t_[" << j << "] = "  << t_[j] << endl;
// 	      #endif
// 	      break;
// 	    }
// 	  }
//
// 	  #ifdef VERBOSE
// 	  cout << "t_[j] : " << t_[j] << endl;
// 	  #endif
//
// 	  assert(t_[j] != 0);
//
// 	  vector<node> W_r;
// 	  node p_r;
//
// 	  #ifdef VERBOSE
// 	  RationalType r = delta/8;
// 	  cout << r << " " << x0 << endl;
// 	  #endif
//
// 	  RationalType l0 = 0;
// 	  cout << " OK here" << endl;
// 	  create_petal(X, Y[j-1], t_[j], x0, l0 , delta/8, W_r, p_r);
// 	  cout << "Reaches " << endl;
// 	  #ifdef VERBOSE
// 	  cout << "petal pushed: " << endl;
// 	  for(auto a: W_r)
// 	  {
// 	    cout << a << " , ";
// 	  }
// 	  cout << endl;
// 	  #endif
//
// 	  petals.push_back(W_r);
//
// 	  #ifdef VERBOSE
// 	  cout << "p_r: " << p_r << endl;
// 	  #endif
//
// 	  petal_centers.push_back(p_r);
// 	  Y[j] = subtract_vectors(Y[j-1], W_r);
//
// 	  cout << " here " << endl;
// 	  vector<node> path_xj_tj = get_shortest_path(p_r, t_[j], X);
// 	  cout << " now here " << endl;
// 	  for(unsigned int i = 0; i < path_xj_tj.size()-1; i++)
// 	  {
// 	    arc a = get_edge(path_xj_tj[i], path_xj_tj[i+1]);
// 	    weights[abs(a)] /= 2;
// 	  }
//
// 	  vector<node> path_x0_tj = get_shortest_path(x0, t_[j], X);
//
// 	  #ifdef VERBOSE
// 	  cout << "path_x0_tj-: x0: " << x0 << "	p_r: " << p_r << endl;
// 	  for(node v: path_x0_tj)
// 	  {
// 	    cout << v << " , " ;
// 	  }
// 	  cout << endl;
// 	  #endif
// 	  for(unsigned int i=0; i < path_x0_tj.size(); i++)
// 	  {
// 	    if(path_x0_tj[i] == p_r)
// 	    {
// 	      y[j] = path_x0_tj[i+1];
// 	      break;
// 	    }
// 	  }
//
// 	  #ifdef VERBOSE
// 	  cout << "y[j] = " << y[j] << endl;
// 	  #endif
// 	  j++;
// 	}
//
// 	return j;
//       }
//
//       template<typename IntegerType, typename RationalType>
//       unsigned int
//       Graph<IntegerType, RationalType>::petal_decomposition(vector<node>& X,
// 							    node x0,
// 							    unsigned int t,
// 							    vector<vector<node> >& petals,
// 							    vector<node>& petal_centers,
// 							    vector<node>& y,
// 							    vector<node>& t_
//       )
//       {
//
// 	cout << endl;
//
// 	RationalType delta = rad(x0, X);
//
// 	#ifdef VERBOSE
// 	cout << "delta: " << delta << endl;
// 	#endif
// 	RationalType r0 = delta/2;
// 	vector< vector<node>> Y (no_of_vertices + 1);
//
// 	Y[0] = X;
// 	#ifdef VERBOSE
// 	cout << "here: " << endl;
// 	#endif
// 	unsigned int j = 1;
//
// 	/* creating the first petal */
// 	vector<node> previous(no_of_vertices + 1, 0);
//
// 	#ifdef VERBOSE
// 	cout << "x0: " << x0 << "   " << "t: " << t << endl;
// 	cout << "X: " << endl;
// 	for(auto a: X)
// 	{
// 	  cout << a << " , " ;
// 	}
// 	cout << endl;
// 	cout << "Graph: " << endl;
// 	for(unsigned int i = 1; i <= no_of_vertices; i++)
// 	{
// 	  cout << i << ": ";
// 	  for(auto a: incident_edges[i])
// 	  {
// 	    cout << a << " , ";
// 	  }
// 	  cout << endl;
// 	}
// 	#endif
//
// 	RationalType d_x0_t = dijkstra(x0, t, X, previous);
// 	#ifdef VERBOSE
// 	cout << "d_x0_t: " << d_x0_t << endl;
// 	//node t1_prime = 0;
// 	cout << "Previous: " << endl;
//
// 	for(auto a: previous)
// 	{
// 	  cout << a << " , " ;
// 	}
// 	cout << endl;
// 	#endif
//
// 	#ifdef VERBOSE
// 	cout << d_x0_t << " >= " <<  5 << " * " << delta << " / " << 8 << endl;
// 	#endif
// 	if(d_x0_t >= 5*delta/8)
// 	{
// 	  vector<node> X1;
// 	  node x1;
// 	  create_petal(X, X, t, x0, d_x0_t - 5*delta/8, d_x0_t - delta/2, X1, x1);
// 	  vector<node> path_x0_t = get_shortest_path( x0, t, X );
//
//
//
// 	  for(unsigned int i = 0; i < path_x0_t.size(); i++)
// 	  {
// 	    if(path_x0_t[i] == x1)
// 	    {
// 	      y[1] = path_x0_t[i+1];
// 	      break;
// 	    }
// 	  }
//
// 	  t_[0] = y[1];
// 	  t_[1] = t;
// 	  j = 2;
//
// 	}
//
// 	else
// 	{
//
// 	  t_[0] = t;
// 	}
//
//
// 	j = create_remaining_petals(j, r0, x0, delta, X, Y, petals, petal_centers, y, t_);
//
//
// 	unsigned int s = j - 1;
//
// 	#ifdef VERBOSE
// 	cout << endl << "s: " << s << endl;
// 	#endif
// 	create_stigma(petals, Y, s, y, t_);
//
// 	return s;
//       }
//
//       template<typename IntegerType, typename RationalType>
//       vector<node>
//       Graph<IntegerType, RationalType>::get_W_r(node x0,
// 						node t,
// 						RationalType r,
// 						vector<node>& Y,
// 						vector<node>& path_x0_t
//       )
//       {
//
// 	#ifdef VERBOSE
// 	cout << "get_W_r():    " << "x0: " << x0 << " , t: " << t << " , r: " << r << endl;
// 	#endif
//
// 	vector<node> W_r;
// 	W_r.clear();
//
// 	#ifdef VERBOSE
// 	cout << "printing the path_x0_t: " ;
// 	for(node p: path_x0_t)
// 	{
// 	  cout << p << " , ";
// 	}
// 	cout << endl;
// 	#endif
//
// 	vector<RationalType> dist_t = dijkstra(t, Y);
// 	for(node p: path_x0_t)
// 	{
// 	  vector<node> previous(2*no_of_vertices+1, 0);
//
// 	  #ifdef VERBOSE
// 	  cout << "p: " << p << endl;
// 	  cout << "t: " << t << endl;
// 	  #endif
//
// 	  RationalType dypt = dist_t[p];
//
// 	  #ifdef VERBOSE
// 	  cout << "dypt: " << (dypt) << endl;
// 	  #endif
//
// 	  if(dypt <= r)
// 	  {
// 	    #ifdef VERBOSE
// 	    cout << "r: " << (r) << endl << "dypt: " << (dypt) << endl << "Ball will be called.. " << endl;
// 	    #endif
// 	    RationalType r_minus_dypt_over_2 = ((r - dypt) ) / 2;
// 	    #ifdef VERBOSE
// 	    cout << "r_minus_dypt_over_2 = " << r_minus_dypt_over_2 << endl;
// 	    #endif
// 	    #ifdef VERBOSE
// 	    cout << "Y: " ;
// 	    for(auto a: Y)
// 	    {
// 	      cout << a << " , " ;
// 	    }
// 	    cout << endl;
// 	    #endif
// 	    Ball(x0, p, r_minus_dypt_over_2, Y, W_r);
// 	  }
//
// 	  #ifdef VERBOSE
// 	  cout << "check of ball : " << endl;
// 	  cout << "r : " << r << endl << "dypt: " << dypt << endl;
// 	  cout << "x0: " << x0 << endl << "p: " << p << endl << "(r - dypt)/2 : " << (r - dypt)/2 << endl;
//
// 	  cout << " Y: " ;
// 	  for(auto a: Y)
// 	  {
// 	    cout << a << " , " ;
// 	  }
// 	  cout << endl;
// 	  cout << "W_r: " ;
// 	  for(auto a: W_r)
// 	  {
// 	    cout << a << " , " ;
// 	  }
// 	  cout << endl;
// 	  #endif
// 	}
//
// 	cout << "get_W_r " << endl;
//
// 	return W_r;
//       }
//
//
//       template<typename IntegerType, typename RationalType>
//       vector<node>
//       Graph<IntegerType, RationalType>::get_shortest_path(node source,
// 							  node target,
// 							  vector<node>& X
//       )
//       {
//
// 	#ifdef VERBOSE
// 	cout << "get_shortest_path() " << endl;
// 	cout << "source: " << target << endl << "target: " << target << endl;
// 	#endif
// 	vector<node> previous(2*no_of_vertices + 1, 0);
//
//
// 	vector<node> path_x0_t;
//
//
// 	dijkstra(source, target, X, previous);
// 	path_x0_t.push_back(target);
// 	if(source == target)
// 	{
// 	  return path_x0_t;
// 	}
//
// 	node v = target;
//
// 	if(previous[v] == 0 )
// 	{
// 	  v = source;
// 	}
//
// 	while(true)
// 	{
//
// 	  v = previous[v];
// 	  path_x0_t.push_back(v);
//
// 	  if(v == source) break;
// 	}
//
// 	return path_x0_t;
//
//       }
//
//
//       template<typename IntegerType, typename RationalType>
//       void
//       Graph<IntegerType, RationalType>::create_petal(vector<node>& X,
// 						     vector<node>& Y,
// 						     node t,
// 						     node x0,
// 						     RationalType l0,
// 						     RationalType hi,
// 						     vector<node>& W_r,
// 						     node& p_r
//       )
//       {
// 	RationalType R = hi - l0;
//
// 	#ifdef VERBOSE
// 	cout << "create_petal() " << endl;
// 	#endif
//
// 	#ifdef VERBOSE
// 	cout << "R: " << R << endl;
// 	#endif
//
// 	vector<node> all_verticies(no_of_vertices,0);
//
// 	for(unsigned int i = 1; i <= no_of_vertices; i++)
// 	{
// 	  all_verticies[i-1] = i;
// 	}
//
// 	vector<node> path_x0_t = get_shortest_path(x0, t, all_verticies);
// 	vector<node> previous(2*no_of_vertices + 1, 0);
//
// 	#ifdef VERBOSE
// 	cout << "no_of_vertices: " << no_of_vertices << endl;
// 	#endif
//
// 	long long L = ceil(log(log(no_of_vertices)));
// 	long long q = 1;
//
// 	long long r =0;
//
// 	while(true)
// 	{
//
// 	  #ifdef VERBOSE
// 	  cout << " while loop ... " << endl;
// 	  #endif
//
// 	  RationalType a = l0 + (q -1)*R/L;
// 	  #ifdef VERBOSE
// 	  cout << "l0: " << l0 << endl;
// 	  cout << "q: " << q << endl;
// 	  cout << "R: " << R << endl;
// 	  cout << "a: " << a << endl;
// 	  #endif
// 	  vector<node> W_a = get_W_r(x0, t, a, Y, path_x0_t );
// 	  #ifdef VERBOSE
// 	  cout << "W_a: ";
//
// 	  for(auto a: W_a)
// 	  {
// 	    cout << a << " , ";
// 	  }
//
// 	  cout << endl;
// 	  #endif
//
// 	  #ifdef VERBOSE
// 	  cout << "t = " << t << endl;
// 	  #endif
//
// 	  RationalType rad_W_a = rad(t, W_a);
//
// 	  #ifdef VERBOSE
// 	  cout << "rad_W_a: " << rad_W_a << endl;
// 	  #endif
//
// 	  RationalType chi = (size_E(X) + 1)/vol_x(Y, W_a, t, rad_W_a );
// 	  #ifdef VERBOSE
// 	  cout << "chi = " << (size_E(X) + 1) << " / " << vol_x(Y, W_a, t, rad_W_a ) << " = "<< chi << endl;
// 	  #endif
// 	  assert(chi > 0);
// 	  for(long long r1 = a; r1 <= a + R/L; r1++)
// 	  {
//
// 	    #ifdef VERBOSE
// 	    cout << "for loop .. " << a << " to " << a + R/L << endl;
// 	    cout << "r1 = " << r1 << endl;
// 	    #endif
//
// 	    vector<node> W_r = get_W_r(x0, t, r1, Y, path_x0_t);
//
// 	    RationalType rad_W_r = rad(t, W_r);
// 	    #ifdef VERBOSE
// 	    cout << "del(Y, W_r): " << endl;
// 	    #endif
// 	    vector<arc> check = del(Y,W_r);
//
// 	    #ifdef VERBOSE
// 	    for( auto a: check)
// 	    {
// 	      cout << a << " - > " << weights[abs(a)] << endl;
// 	    }
//
// 	    cout <<  endl;
//
// 	    cout << "r1 : " << r1 << endl;
// 	    cout << "rad_W_r " << rad_W_r << endl;
//
// 	    #endif
//
// 	    #ifdef VERBOSE
// 	    cout << R  << " * " << cost( check ) << " <= " << vol_x(Y, W_r,t , rad_W_r) << " * " << 8 << " * " << L << " * " << log(chi) << endl;
// 	    #endif
// 	    if( R * cost( check ) <= vol_x(Y, W_r,t , rad_W_r) * 8 * L * log(chi) )
// 	    {
// 	      r = r1;
// 	      break;
// 	    }
// 	  }
//
// 	  vector<node> W_r = get_W_r(x0, t, r, Y, path_x0_t);
//
// 	  RationalType rad_W_r = rad(t, W_r);
// 	  #ifdef VERBOSE
// 	  cout << "W_r: ";
// 	  for(auto a: W_r)
// 	  {
// 	    cout << a << " , ";
// 	  }
// 	  cout << endl;
// 	  cout <<  (vol_x(Y, W_r, t, rad_W_r) - 1)  << " * " << pow(2, pow(log(no_of_edges) , 1 - q/L) ) << " <= " << 2*size_E(X) << endl;
// 	  #endif
// 	  if( (vol_x(Y, W_r, t, rad_W_r) - 1) * pow(2, pow(log(no_of_edges) , 1 - q/L) ) <= 2*size_E(X))
// 	  {
// 	    break;
// 	  }
// 	  else{
//
// 	    q++;
// 	  }
// 	}
//
// 	#ifdef VERBOSE
// 	cout << "r = " << r << endl;
// 	cout << "finding pr': " << endl;
// 	cout << "t[j] = " << t << endl;
// 	#endif
//
// 	vector<RationalType> dist_t = dijkstra(t, X);
// 	for(auto aa: boost::adaptors::reverse(path_x0_t))
// 	{
// 	  #ifdef VERBOSE
// 	  cout << "aa: " << aa << endl;
// 	  cout << dist_t[aa] << " <= " << r << endl;
// 	  #endif
// 	  if(dist_t[aa] <= r)
// 	  {
// 	    p_r = aa;
// 	    break;
// 	  }
// 	}
//
// 	#ifdef VERBOSE
// 	cout << "p_r = " << p_r << endl;
// 	#endif
//
// 	W_r = get_W_r(x0, t, r, Y, path_x0_t);
// 	//p_r = aa;
//       }
//
//
//       template<typename IntegerType, typename RationalType>
//       void
//       Graph<IntegerType, RationalType>::Ball_x(node x0,
// 					       RationalType r0,
// 					       vector<node>& X,
// 					       vector<node>& nodes_in_ball
//       )
//       {
// 	vector<RationalType> dist = dijkstra(x0, X);
// 	for(node x: X)
// 	{
// 	  if(dist[x] <= r0)
// 	  {
// 	    nodes_in_ball.push_back(x);
// 	  }
// 	}
//       }

template<typename IntegerType, typename RationalType>
vector<node>
Graph<IntegerType, RationalType>::subtract_vectors (vector<node>& A,
    vector<node>& B
                                                   )
{
	vector<node> subtracted_result;
	for (auto a : A) {
		bool not_there = false;
		for (auto b : B) {
			if (a == b ) not_there = true;
		}

		if (not_there == false) subtracted_result.push_back(a);
	}

	return subtracted_result;
}


template<typename IntegerType, typename RationalType>
unsigned int
Graph<IntegerType, RationalType>::subtract_vectors_size (vector<node>& A,
    vector<node>& B
                                                        )
{
	unsigned int number = 0;
	for (auto a : A) {
		bool not_there = false;
		for (auto b : B) {
			if (a == b ) not_there = true;
		}
		if (not_there == false) number++;
	}

	return number;
}

template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>::Ball(node x,
                                       node y,
                                       RationalType r,
                                       vector<node>& X,
                                       vector<node>& nodes_in_ball
                                      )
{
	assert( x != 0);
	assert( y != 0);
	assert( log(10) != 1 );

	vector<RationalType> dist_x = dijkstra(x, X);
	vector<RationalType> dist_y = dijkstra(y, X);

#ifdef VERBOSE
	cout << "Ball() " << x << " , " << y << endl;
#endif

	for (node z : X) {

		RationalType dXxy = dist_x[y];
#ifdef VERBOSE
		cout << "dXxy: " << dXxy << endl;
#endif
		RationalType dXyz = dist_y[z];
#ifdef VERBOSE
		cout << "dXyz: " << dXyz << endl;
#endif
		RationalType dXxz = dist_x[z];
#ifdef VERBOSE
		cout << "dXxz: " << dXxz << endl;
#endif

#ifdef VERBOSE
		cout << dXxy + dXyz - dXxz << " <= " << r << endl;
#endif
		if (dXxy + dXyz - dXxz <= r)
		{
			nodes_in_ball.push_back(z);
#ifdef VERBOSE
			cout << "pushed into the ball :  " << z << endl;
#endif
		}

	}
}

/** \brief Creates a minimum spanning tree w.r.t to the unrounded_resistances
 *
 * @param ST The Spanning Tree
 *
 */
template<typename IntegerType, typename RationalType>
void
Graph<IntegerType, RationalType>::create_MST_wrt_unrounded_resistances_on_original_graph
(SpanningTree<Graph, IntegerType, RationalType>& ST, Network<Graph, IntegerType, RationalType>& N)
{
	tree_edges.clear();
	non_tree_edges.clear();

	std::vector<unsigned int> permutation( no_of_edges + 1 );

	for (unsigned int i = 0; i <= no_of_edges; i++)
	{

		permutation[i] = i;

	}

	max_val.resize( no_of_edges + 1 );

	for (unsigned int a = 1; a <= no_of_edges; a++) {
		RationalType resistance_effective =
		  to_double(N.arcdata[a].resistance_roof * (N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper))
		  / to_double(N.arcdata[a].resistance_roof + N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper);
		if (N.arcdata[a].resistance_roof < N.arcdata[a].resistance_lower + N.arcdata[a].resistance_upper)
		{
			max_val[a] = resistance_effective;
		}
		else
		{
			max_val[a] = resistance_effective;
		}
	}

	std::sort(permutation.begin() + 1, permutation.end(), [this] (unsigned int a, unsigned int b) -> bool { return max_val[a] <
	          max_val[b]; } );


	make_tree(ST, permutation);
}



template<typename IntegerType, typename RationalType>
void Graph<IntegerType, RationalType>::create_low_stretch_tree_wrt_unrounded_resistances_on_original_graph(
  SpanningTree<Graph, IntegerType, RationalType>& ST,
  Network<Graph<IntegerType, RationalType>,
  IntegerType, RationalType>& N)
{

	create_MST_wrt_unrounded_resistances_on_original_graph( ST, N );

	get_depth_and_children(ST);

	euler_tour_and_get_subtree_sizes(ST);


// 	for(auto a: tree_edges){
// 	  if(original_auxiliary_transformed_corrospondence[abs(a)].is_invalid) original_auxiliary_transformed_corrospondence[abs(a)].cant_be_removed = true;
// 			       original_auxiliary_transformed_corrospondence[abs(a)].is_invalid = false;
//
// 	}

	sort( tree_edges.begin(), tree_edges.end(), [this, &ST]( arc a, arc b ) -> bool { return ST.depth[this->head(a)] < ST.depth[this->head(b)]; } );

	ST.fill_LCA_vector();



#ifdef VERBOSE
	cout << ST.root << endl;
#endif

	ST.update_sum_to_the_root_wrt_unrounded_resistances(ST.root, N);


	unrounded_resistances_accross_cycles.clear();

	for (unsigned int i = 0; i < non_tree_edges.size(); i++) {
		arc a = non_tree_edges[i];
		unrounded_resistances_accross_cycles.push_back(ST.compute_unrounded_resistance_accross(a, N));
	}

	ST.update_sum_to_the_root(ST.root, N);

	resistances_accross_cycles.clear();

#ifdef VERBOSE
	cout << non_tree_edges.size() << endl;
#endif
	for (unsigned int i = 0; i < non_tree_edges.size(); i++) {
#ifdef VERBOSE
		cout << "[" << i << "]" << endl;
#endif


#ifdef VERBOSE
		arc a = non_tree_edges[i];
		cout << "a = " << a << endl;
		cout << "lca = " << ST.LCA[a] << endl;
#endif
		//resistances_accross_cycles.push_back(ST.compute_resistance_accross(a, N));

	}



}

//       /** \brief Creates a Low stretch tree with respect to unrounded-resistances
//        *
//        * @param ST The Spanning Tree
//        *
//        */
//       template<typename IntegerType, typename RationalType>
//       void
// Graph<IntegerType,
// RationalType>::create_low_stretch_tree_wrt_unrounded_resistances_on_original_graph(SpanningTree<Graph,
//   IntegerType, RationalType>& ST, Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N)
//       {
// 	create_MST_wrt_unrounded_resistances_on_original_graph( ST );
//
// 	get_depth_and_children(ST);
//
// 	euler_tour_and_get_subtree_sizes(ST);
//
// 	for(auto a: tree_edges){
// 	  if(original_auxiliary_transformed_corrospondence[abs(a)].is_invalid) original_auxiliary_transformed_corrospondence[abs(a)].cant_be_removed = true;
// 			       original_auxiliary_transformed_corrospondence[abs(a)].is_invalid = false;
//
// 	}
//
// 	sort( tree_edges.begin(), tree_edges.end(), [this,&ST]( arc a, arc b ) -> bool { return ST.depth[this->head(a)] < ST.depth[this->head(b)]; } );
//
// 	ST.fill_LCA_vector();
//
//
//
// 	#ifdef VERBOSE
// 	cout << ST.root << endl;
// 	#endif
//
// 	ST.update_sum_to_the_root_wrt_unrounded_resistances(ST.root);
//
//
// 	unrounded_resistances_accross_cycles.clear();
//
// 	for(unsigned int i=0; i< non_tree_edges.size(); i++){
// 	  arc a = non_tree_edges[i];
// 	  unrounded_resistances_accross_cycles.push_back(ST.compute_unrounded_resistance_accross(a));
// 	}
//
// 	ST.update_sum_to_the_root(ST.root, N);
//
// 	resistances_accross_cycles.clear();
//
// 	#ifdef VERBOSE
// 	cout << non_tree_edges.size() << endl;
// 	#endif
// 	for(unsigned int i=0; i< non_tree_edges.size(); i++){
// 	  #ifdef VERBOSE
// 	  cout <<"[" << i << "]" << endl;
// 	  #endif
//
// 	  arc a = non_tree_edges[i];
// 	  #ifdef VERBOSE
// 	  cout << "a = " << a << endl;
// 	  cout << "lca = " << ST.LCA[a] << endl;
// 	  #endif
// 	  resistances_accross_cycles.push_back(ST.compute_resistance_accross(a, N));
//
// 	}
//
//       }
//
//       /** \brief Creates a Low Stretch tree with respect to unrounded resistances
//        *
//        * @param x The Primal Solution
//        * @param ST The Spanning Tree
//        *
//        */
//       template<typename IntegerType, typename RationalType>
//       void
//       Graph<IntegerType, RationalType>::create_low_stretch_tree_wrt_unrounded_resistances(vector<RationalType>& x,
// 											  SpanningTree<Graph, IntegerType, RationalType>& ST,
// 											  Network<Graph, IntegerType, RationalType>& N
//       )
//       {
// 	//create_SPT(ST);
// 	cout << "1" << endl;
// 	create_MST_wrt_unrounded_resistances(x, ST);
// 	cout << "2" << endl;
// 	//create_BFS_tree_with_random_root(ST);
// 	//create_MST(ST);
// 	//form_petal_decomposition_spanning_tree(ST);
// 	//form_spanning_tree_alon(ST, 20.0);
//
// 	get_depth_and_children(ST);
// 	cout << "3" << endl;
// 	#ifdef VERBOSE
// 	ST.print_children();
// 	#endif
//
// 	euler_tour_and_get_subtree_sizes(ST);
// 	cout << "4" << endl;
// 	sort( tree_edges.begin(), tree_edges.end(), [this,&ST]( arc a, arc b ) -> bool { return ST.depth[this->head(a)] < ST.depth[this->head(b)]; } );
// 	cout << "5" << endl;
// 	#ifdef VERBOSE
// 	assert( ST.root != 0 );
// 	int size_at_s = ST.size_sub_tree_rooted_at[ST.root];
// 	unsigned int initial_tree_index = 1;
// 	const node d = ST.find_vertex_separator(ST.root, size_at_s, initial_tree_index);
// 	cout << " vertex separator: " << d << endl;
// 	cout << "root: " << ST.root << endl;
//
//
// 	for( auto a : tree_edges ) {
// 	  cout << a << " = (" << tail(a) << "," << head(a) << ") " << ST.depth[tail(a)] << " " << ST.depth[head(a)] << endl;
// 	}
// 	cout<< endl << "tree edges are sorted..." << endl << endl;
// 	#endif
//
// 	ST.fill_LCA_vector();
// 	cout << "6" << endl;
// 	ST.update_sum_to_the_root_wrt_unrounded_resistances(ST.root);
//
// 	cout << "7" << endl;
//
// 	unrounded_resistances_accross_cycles.clear();
// 	cout << "8" << endl;
// 	for(unsigned int i=0; i< non_tree_edges.size(); i++){
// 	  arc a = non_tree_edges[i];
// 	  unrounded_resistances_accross_cycles.push_back(ST.compute_unrounded_resistance_accross(a));
// 	}
// 	cout << "9" << endl;
//
// 	ST.update_sum_to_the_root(ST.root, N);
// 	cout << "10" << endl;
// 	resistances_accross_cycles.clear();
// 	cout << "11" << endl;
// 	for(unsigned int i=0; i< non_tree_edges.size(); i++){
// 	  arc a = non_tree_edges[i];
// 	  resistances_accross_cycles.push_back(ST.compute_resistance_accross(a, N));
// 	}
// 	cout << "12" << endl;
//       }
//
//
//
//
//
//
//       /** \brief Create a low stretch tree
//        *
//        * Creates a low stretch spanning tree. Finds the depth, children of ech node
//        *
//        * @param SpanningTree ST
//        *
//        */
//       template<typename IntegerType, typename RationalType>
//       void
//       Graph<IntegerType, RationalType>:: create_low_stretch_tree(
// 	SpanningTree<Graph, IntegerType, RationalType>& ST,
// 	Network<Graph<IntegerType,RationalType>, IntegerType, RationalType>& N
//       )
//       {
//
// 	#ifdef VERBOSE
// 	cout << "creat LST!! " << endl;
// 	#endif
//
// 	form_petal_decomposition_spanning_tree(ST);
// 	get_depth_and_children(ST);
//
// 	#ifdef VERBOSE
// 	ST.print_children();
// 	#endif
//
// 	euler_tour_and_get_subtree_sizes(ST);
//
// 	sort( tree_edges.begin(), tree_edges.end(), [this,&ST]( arc a, arc b ) -> bool { return ST.depth[this->head(a)] < ST.depth[this->head(b)]; } );
//
// 	#ifdef VERBOSE
// 	assert( ST.root != 0 );
// 	int size_at_s = ST.size_sub_tree_rooted_at[ST.root];
// 	unsigned int initial_tree_index = 1;
// 	const node d = ST.find_vertex_separator(ST.root, size_at_s, initial_tree_index);
// 	cout << " vertex separator: " << d << endl;
// 	cout << "root: " << ST.root << endl;
//
//
// 	for( auto a : tree_edges ) {
// 	  cout << a << " = (" << tail(a) << "," << head(a) << ") " << ST.depth[tail(a)] << " " << ST.depth[head(a)] << endl;
// 	}
// 	cout<< endl << "tree edges are sorted..." << endl << endl;
// 	#endif
//
// 	ST.fill_LCA_vector();
//
// 	ST.update_sum_to_the_root(ST.root, N);
//
//
//
// 	resistances_accross_cycles.clear();
//
// 	for(unsigned int i=0; i< non_tree_edges.size(); i++){
// 	  arc a = non_tree_edges[i];
// 	  resistances_accross_cycles.push_back(ST.compute_resistance_accross(a, N));
// 	}
//
//       }


/** \brief Create a Random Graph
 *
 * Creates a Random instance of a graph
 *
 */
template <
  typename IntegerType,
  typename RationalType
  >
void  Graph<IntegerType, RationalType>::create_graph() {
	while (true) {

		// sample two nodes
		node node_1, node_2;

		double p = rg.rng();
		node_1 = static_cast<int>( ceil( p * no_of_vertices ) );

		p = rg.rng();
		node_2 = static_cast<int>( ceil( p * no_of_vertices ) );

		// create arc between them
		create_edge(node_1, node_2);

		// print inserted arc
		cout  << "Arc " << count_tree_edges << ":"
		      << " (" << node_1 << "," << node_2 << ")" << endl;

		// break out if graph got connected
		if (count_tree_edges == no_of_vertices - 1) {
			break;
		}
	}
}



/** \brief Reads the graph from a file
 *
 * @param filepath The path of the graph file
 * @param m The number of edges
 *
 */
//      template<typename IntegerType, typename RationalType>
//      void Graph<IntegerType, RationalType>::read_graph(string filepath, int m)
//      {
// costs.reserve(m+1);
// capacities.reserve(m+1);

// ptree tree;
// read_xml(filepath, tree);

// const ptree & graphml = tree.get_child("graphml", empty_ptree());
// const ptree & graph   = graphml.get_child("graph", empty_ptree());

// BOOST_FOREACH(const ptree::value_type & nore, graph){

//   const ptree & nore_attrs = nore.second.get_child("<xmlattr>", empty_ptree());
//   bool node = false;
//   bool edge = false;
//   int nodeid = 0;
//   int source = 0;
//   int target = 0;
//   BOOST_FOREACH(const ptree::value_type & nore_attr, nore_attrs){
//     if (strncmp(nore_attr.first.data(),"id",2) == 0){
//       node   = true;
//       nodeid = stoi(nore_attr.second.data());
//     }
//     else {
//       if (strncmp(nore_attr.first.data(),"source",6) == 0){
// 	edge = true;
// 	source = stoi(nore_attr.second.data());
//       }
//       if (strncmp(nore_attr.first.data(),"target",6) == 0){
// 	assert(edge);
// 	target = stoi(nore_attr.second.data());
//       }
//     }
//   }

//   if (edge) {
//     ++no_of_edges;
//     new_edge(no_of_edges,source,target);
//   }


//   const ptree & nore_tree = nore.second.get_child("", empty_ptree());
//   int id = 0;
//   BOOST_FOREACH(const ptree::value_type & nore_data, nore_tree){
//     if (strncmp(nore_data.first.data(),"data",4) == 0){
//       if (node) {
// 	#ifdef VERBOSE
// 	cout << "b[" << nodeid << "]: " << nore_data.second.data() << endl;
// 	#endif
// 	RationalType demand = stoi(nore_data.second.data());
// 	demands[nodeid] = demand;
//       }
//       if (edge) {
// 	if (id==0){
// 	  id = 1;
// 	  RationalType capacity = stoi(nore_data.second.data());
// 	  // assert(stoi(nore_data.second.data()) % 2 != 0);
// 	  #ifdef VERBOSE
// 	  cout << "(" << source << "," << target << "): cap: " << capacity << endl;
// 	  #endif
// 	  capacities[no_of_edges] = capacity;
// 	}
// 	else {
// 	  if (id==1){
// 	    id = 2;
// 	    RationalType cost = stoi(nore_data.second.data());
// 	    #ifdef VERBOSE
// 	    cout << "(" << source << "," << target << "): cost: " << cost << endl;
// 	    #endif
// 	    costs[no_of_edges] = cost;

// 	  }
// 	}
//       }
//     }
//   }
// }
//      }


//       /** \brief Assign Random Resistances
//        *
//        * Assigns Random Resistances to the edges in the graph
//        *
//        */
//       template<typename IntegerType, typename RationalType>
//       void Graph<IntegerType, RationalType>::get_random_resistances()
//       {
// 	resistances.resize( no_of_edges+1 );
// 	for( auto& r : resistances )
// 	{
// 	  RationalType random_resistance = rg.rng();
// 	  r = random_resistance*400 + 1;
// 	}
//       }


/** \brief Create a new Edge
 *
 * Adds an arc between the given nodes and returns the arc_map
 *
 * @param node u
 * @param node v
 * @return arc
 */
template<typename IntegerType, typename RationalType>
arc Graph<IntegerType, RationalType>::new_edge( node u, node v) {
	no_of_edges++;
	heads.resize( no_of_edges + 1 );
	tails.resize( no_of_edges + 1 );
	return new_edge( no_of_edges, u, v );
}

template<typename IntegerType, typename RationalType>
void Graph<IntegerType, RationalType>::remove_arc(arc a)
{

	node u = tail(a);
	node v = head(a);

	for ( std::vector<int>::iterator iter = incident_edges[v].begin(); iter != incident_edges[v].end(); ++iter) {
		if ( abs(*iter) == abs(a)) {
			incident_edges[v].erase( iter );
			break;
		}
	}

	for ( std::vector<int>::iterator iter = incident_edges[u].begin(); iter != incident_edges[u].end(); ++iter ) {
		if ( abs(*iter) == abs(a)) {
			incident_edges[u].erase( iter );
			break;
		}
	}
}

/** \brief creates a new edge
 *
 * @param a Arc
 * @param u node u
 * @param v node v
 *
 * @return Returns the arc created
 *
 */
template<typename IntegerType, typename RationalType>
arc
Graph<IntegerType, RationalType>::new_edge( arc a,
    node u,
    node v
                                          )
{
	assert( a > 0 );
	assert( static_cast<unsigned int>( a ) <= no_of_edges );

#ifdef VERBOSE
	cout << "new edge: no_of_edges " << no_of_edges << " = (" << u << "," << v << ")" << endl;
#endif

	tails[a] = u;
	heads[a] = v;

	incident_edges[u].push_back( a );
	incident_edges[v].push_back( -a);

#ifdef VERBOSE
	cout << "tails[" << a << "]: " << u << endl;
	cout << "heads[" << a << "]: " << v << endl;
	cout << "incident edges[" << u << "] pushed back: " << a << " length of incidentedges[" << u << "]: " << incident_edges[u].size() << endl;
	cout << "incident edges[" << v << "] pushed back: " << -a << " length of incidentedges[" << v << "]: " << incident_edges[v].size() << endl;
#endif

	return a;
}


/** \brief creates arc between the given the given nodes
 *
 * @param u The head
 * @param v The tail
 *
 * @return Returns the arc created
 *
 */
template<typename IntegerType, typename RationalType>
arc
Graph<IntegerType, RationalType>::create_edge(node u,
    node v
                                             )
{
	int flag = 0;
	for (j = 1; j <= no_of_edges; j++) {

		if (heads[j] == u && tails[j] == v) {
			flag = 1;
		}
		if (heads[j] == v && tails[j] == u) {
			flag = 1;
		}
	}

	if (u != v && flag == 0) {
		const double p = rg.rng();
		const arc a = p > 0.5 ? new_edge( v, u ) : new_edge( u, v );
		unite( a );
		return a;
	} else {
		return 0;
	}
}


/** \brief Given the nodes, returns the edge joining them
 *
 * @param a Node
 * @param b Node
 *
 * @return Returns the edge joining a and b
 *
 */
template<typename IntegerType, typename RationalType>
int
Graph<IntegerType, RationalType>::get_edge(const node a,
    const node b
                                          )
{
	int edge = 0;
	for (unsigned int i = 1; i <= no_of_edges; i++) {
		if (heads[i] == a && tails[i] == b) {
			edge = -i;
			break;
		}

		if (heads[i] == b && tails[i] == a) {
			edge = i;
			break;
		}
	}
	return edge;
}


/* get_path takes the arguement as two verticies and returns a path between the two verticies in the spanning tree
 * The implementation of this take linear time in the number of edges but the way it is implemented is quite cumbursome in my opinion */
template<typename IntegerType, typename RationalType>
vector<arc> Graph<IntegerType, RationalType>::get_path( node s, node t, SpanningTree<Graph, IntegerType, RationalType> &ST )
{
	node v = s;
	node w = t;
	assert( v > 0 );
	assert( w > 0 );
	deque<arc> uppath, downpath;

	while ( v != w ) {
		assert( v > 0 );
		assert( w > 0 );
		if ( ST.depth[v] > ST.depth[w] ) {
			const arc a = arcs_to_root[v];
			assert( a != 0 );
			assert( v == tail(a) );
			uppath.push_back( a );
			v = head( a );
		} else {
			const arc a = arcs_to_root[w];
			assert( a != 0 );
			assert( w == tail(a) );
			downpath.push_front( -a );
			w = head( a );
		}
	}
	vector<arc> path( uppath.begin(), uppath.end() );
	path.insert( path.end(), downpath.begin(), downpath.end() );
	return path;
}


/* get_cycle takes a non-tree edge as an arguement and returns the cycle that comes when it is addded to the spanning tree */
template< typename IntegerType, typename RationalType>
std::vector<arc> Graph<IntegerType, RationalType>::get_cycle( arc a, SpanningTree<Graph<IntegerType, RationalType>, IntegerType, RationalType>& ST )
{
	const node v = tail(a);
	const node w = head(a);
	vector<arc> circuit = get_path( w, v, ST  );
	assert( head(circuit.back()) == tail(a) );
	assert( head(a) == tail(circuit.front()) );
	circuit.push_back( a );

	return circuit;
}


template <
  typename Graph,
  typename IType,
  typename RType,
  typename NodeData,
  typename ArcData
  >

vector<double> Network<Graph, IType, RType, NodeData, ArcData>:: get_distribution(vector<IntegerType>& sum_to_root) const {
	std::vector<double> distribution( RationalType(2) * G.non_tree_edges.size() + G.tree_edges.size());
	int counter = 0;
	for (unsigned int i = 0; i < G.non_tree_edges.size(); i++) {

		const int edge = G.non_tree_edges[i];
		IntegerType R_e = G.resistances_accross_cycles[i];
		IntegerType r_e_roof = arcdata[edge].resistance_roof;
#ifdef VERBOSE
		cout << "r_e_roof = " << to_double(r_e_roof) << endl;
#endif
		IntegerType R_e_roof = R_e - arcdata[edge].resistance_lower - arcdata[edge].resistance_upper + arcdata[edge].resistance_roof;
		RationalType prob = to_double(R_e_roof) / to_double(r_e_roof);
		distribution[counter] = prob;
		if(arcdata[edge].infeasibility == 0){
		  distribution[counter] = 0;
		}

#ifdef VERBOSE
		cout << "distribution[" << counter << "] = " << distribution[counter] << " = " <<
		     to_double(R_e_roof) << "/" << to_double(r_e_roof) << endl;
#endif
		counter++;

		IntegerType r_e(0); //= arcdata[edge].resistance_lower + arcdata[edge].resistance_upper;

		if (arcdata[edge].resistance_lower > arcdata[edge].resistance_upper){
#ifdef VERBOSE
			cout << "Lower Edge is the Tree Edge:" << endl;
#endif
			r_e = arcdata[edge].resistance_lower;
		}
		else{
#ifdef VERBOSE
			cout << "Upper Edge is the Tree Edge:" << endl;
#endif
			r_e = arcdata[edge].resistance_upper;
		}


		R_e = G.resistances_accross_cycles[i];
		prob = to_double( R_e ) / to_double( r_e );
#ifdef VERBOSE
		cout << "prob = " << to_double(R_e) << "/" << to_double(r_e) << endl;
#endif
		distribution[counter] = prob;
		if(arcdata[edge].xlower == 0 || arcdata[edge].capacity - arcdata[edge].xlower == 0){
		  distribution[counter] = 0;
		}
#ifdef VERBOSE
		cout << "distribution[" << counter << "] = " << distribution[counter] << endl;
#endif
		counter++;


	}
	for (auto a : (G.tree_edges))
	{
//	  if(arcdata[abs(a)].resistance_roof < arcdata[abs(a)].resistance_lower + arcdata[abs(a)].resistance_upper)
		if (arcdata[abs(a)].resistance_roof < arcdata[abs(a)].resistance_lower ||
		    arcdata[abs(a)].resistance_roof < arcdata[abs(a)].resistance_upper)
		{
			if (arcdata[abs(a)].resistance_lower < arcdata[abs(a)].resistance_upper)
			{
#ifdef VERBOSE
				cout << "The Upper Edge is the Non Tree Edge" << endl;
#endif
				IntegerType sum_to_root_tail = sum_to_root[G.tails[abs(a)]];;
				IntegerType sum_to_root_head = sum_to_root[G.heads[abs(a)]];
				IntegerType R_a = arcdata[abs(a)].resistance_roof + arcdata[abs(a)].resistance_lower +
				                  arcdata[abs(a)].resistance_upper;
				IntegerType r_a = arcdata[abs(a)].resistance_upper;
				RationalType prob = to_double(R_a) / to_double(r_a);
				distribution[counter] = prob;

#ifdef VERBOSE
				cout << "distribution[" << counter << "] = " << distribution[counter] << endl;
#endif
			}
			else
			{
#ifdef VERBOSE
				cout << "The Lower Edge is the Non Tree Edge" << endl;
#endif
				IntegerType sum_to_root_tail = sum_to_root[G.tails[abs(a)]];;
				IntegerType sum_to_root_head = sum_to_root[G.heads[abs(a)]];
				IntegerType R_a = arcdata[abs(a)].resistance_roof + arcdata[abs(a)].resistance_lower +
				                  arcdata[abs(a)].resistance_upper;
				IntegerType r_a = arcdata[abs(a)].resistance_lower;
				RationalType prob = to_double(R_a) / to_double(r_a);
				distribution[counter] = prob;

#ifdef VERBOSE
				cout << "distribution[" << counter << "] = " << distribution[counter] << endl;
#endif
			}
		}
		else
		{
#ifdef VERBOSE
			cout << "The Roof is the Non Tree Edge" << endl;
#endif
			IntegerType sum_to_root_tail = sum_to_root[G.tails[abs(a)]];;
			IntegerType sum_to_root_head = sum_to_root[G.heads[abs(a)]];
			IntegerType R_a = arcdata[abs(a)].resistance_lower + arcdata[abs(a)].resistance_upper
			                  + arcdata[abs(a)].resistance_roof;

			IntegerType r_a = arcdata[abs(a)].resistance_roof;
			RationalType prob = to_double(R_a) / to_double(r_a);
			distribution[counter] = prob;

#ifdef VERBOSE
			cout << "distribution[" << counter << "] = " << distribution[counter] << endl;
#endif
		}

// 	  if (edge < 0) edge *= -1;
// 	  double prob = to_double(arcdata[edge].resistance_lower + arcdata[edge].resistance_upper + arcdata[edge].resistance_roof)/
// 			to_double(arcdata[edge].resistance_roof);
//
//
// 	  if(arcdata[edge].infeasibility != 0)
// 	  {
// 	  distribution[counter] = 1.0;
// 	  }
// 	  else
// 	  {
// 	    distribution[counter] = 0.0;
// 	  }
#ifdef VERBOSE
		cout << "distribution[" << counter << "] = " << distribution[counter] << endl;
#endif
		if(arcdata[abs(a)].supper == 0 || arcdata[abs(a)].capacity - arcdata[abs(a)].xlower == 0){
		 distribution[counter] = 0;  
		}
		
		if(arcdata[abs(a)].sroof == 0 || arcdata[abs(a)].infeasibility == 0){
		  distribution[counter] = 0; 
		}
		
		if(arcdata[abs(a)].xlower == 0 || arcdata[abs(a)].slower == 0){
         	  distribution[counter] = 0;
		}
		
		counter++;

	}

//#ifdef VERBOSE
	cout << "printing the distribution" << endl;
	for (auto a : distribution)
	{
		cout << a << " , " ;
	}
//#endif

	return distribution;
}


//       /** \brief Print the Linear Program
//        *
//        * Prints the linear program corrosponding to the problem
//        *
//        */
//       template<typename IntegerType, typename RationalType>
//       void Graph<IntegerType, RationalType>::print_lp(){
// 	cout << " LP " << endl;
// 	ofstream LP_file;
// 	LP_file.open("LP.lp");
// 	cout << "Minimize" << endl;
// 	LP_file << "Minimize" << endl;
// 	cout << " obj: ";
// 	LP_file << " obj: " ;
// 	for(unsigned int a=1; a<=no_of_edges; a++){
// 	  cout << " + " << costs[a] << " x" << a << endl;
//
// 	  LP_file << " + " << costs[a] << " x" << a << endl;
// 	}
//
// 	cout << "Subject To" << endl;
// 	LP_file << "Subject To" << endl;
// 	for(unsigned int v=1; v<=no_of_vertices; v++){
// 	  cout << " c" << v << ": ";
// 	  LP_file << " c" << v << ": ";
// 	  unsigned int inc_arc_counter = 0;
// 	  for(auto a: incident_edges[v]){
// 	    if (inc_arc_counter == 0){
// 	      if (a<0){
// 		cout << "x" << abs(a) << " ";
// 		LP_file << "x" << abs(a) << " ";
// 	      } else {
// 		cout << "- x" << abs(a) << " ";
// 		LP_file << "- x" << abs(a) << " ";
// 	      }
// 	    } else {
// 	      if (a<0){
// 		cout << "+ x" << abs(a) << " ";
// 		LP_file << "+ x" << abs(a) << " ";
// 	      } else {
// 		cout << "- x" << abs(a) << " ";
// 		LP_file << "- x" << abs(a) << " ";
// 	      }
// 	    }
// 	    inc_arc_counter++;
// 	  }
// 	  cout << " = " << demands[v] << endl;
// 	  LP_file << " = " << demands[v] << endl;
// 	}
//
// 	cout << "Bounds" << endl;
// 	LP_file << "Bounds" << endl;
// 	for(unsigned int a=1; a<=no_of_edges; a++){
// 	  cout << " 0 <= x" << a << " <= " << capacities[a] << endl;
// 	  LP_file << " 0 <= x" << a << " <= " << capacities[a] << endl;
// 	}
// 	LP_file << "Generals" << endl;
// 	for(unsigned int a = 1; a <= no_of_edges; a++)
// 	{
// 	  LP_file << "x" << a << " ";
// 	}
// 	cout << "END" << endl;
// 	LP_file << endl << "END" << endl;
//       }

#endif

