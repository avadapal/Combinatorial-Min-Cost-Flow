#pragma once

#include "min_cost_flow_sspvariant.h"
#include "graph.h"


using RationalType = long long int;
using IntegerType = long long int;
using GraphType = Graph<IntegerType, RationalType>;
using NodeDataType = SSPVariantNodeData<IntegerType, RationalType>;

#ifdef RESP_CAP
using ArcDataType = SSPVariantRCArcData<IntegerType, RationalType>;
#else
#ifdef RESTORE_BALANCED_NODES
using ArcDataType = SSPVariantRBNArcData<RationalType>;
#else
using ArcDataType = BasicArcData<RationalType>;
#endif
#endif

using NetworkType = Network<GraphType, IntegerType, RationalType, NodeDataType, ArcDataType>;


class DefaultNodeIterator;
class DefaultEdgeIterator;
class DefaultNodeAccessor;
class DefaultEdgeAccessor;


class DefaultNodeIterator {

public:
	DefaultNodeIterator() = delete;
	DefaultNodeIterator(const DefaultNodeIterator &) = default;

	DefaultNodeIterator(node v) : current_node(v) {}

	auto operator++() -> DefaultNodeIterator & {
		++current_node;
		return *this;
	}

	auto operator*() const -> node {
		return current_node;
	}

	auto operator==(const DefaultNodeIterator &other) const -> bool {
		return current_node == other.current_node;
	}

	auto operator!=(const DefaultNodeIterator &other) const -> bool {
		return !operator==(other);
	}

private:
	node current_node;

};


class DefaultEdgeIterator {

public:
	DefaultEdgeIterator() = delete;
	DefaultEdgeIterator(const DefaultEdgeIterator &) = default;

	DefaultEdgeIterator(arc a) : current_edge(a) {}

	auto operator++() -> DefaultEdgeIterator & {
		++current_edge;
		return *this;
	}

	auto operator*() const -> arc {
		return current_edge;
	}

	auto operator==(const DefaultEdgeIterator &other) const -> bool {
		return current_edge == other.current_edge;
	}

	auto operator!=(const DefaultEdgeIterator &other) const -> bool {
		return !operator==(other);
	}

private:
	arc current_edge;

};


class DefaultNodeAccessor {

public:
	DefaultNodeAccessor(NetworkType &N, const DefaultNodeIterator &node_it) : nodedata(N.nodedata[*node_it]) {}

	auto get_potential() const -> RationalType {
		return nodedata.potential;
	}

	void set_potential(RationalType potential) {
		nodedata.potential = potential;
	}

	auto get_deficit() const -> IntegerType {
		return nodedata.deficit;
	}

	void set_deficit(IntegerType deficit) {
		nodedata.deficit = deficit;
	}

#ifdef NEW_OPTIMIZATION
	auto get_demand() const -> IntegerType {
		return nodedata.demand;
	}
#endif

	auto get_depth() const -> int {
		return nodedata.depth;
	}

	void set_depth(int depth) {
		nodedata.depth = depth;
	}

	auto is_visited() const -> bool {
		return nodedata.visited;
	}

	void set_visited(bool visited) {
		nodedata.visited = visited;
	}

#ifdef RESTORE_BALANCED_NODES
	auto get_deficit_delta() const -> IntegerType {
		return nodedata.deficit_delta;
	}

	void set_deficit_delta(IntegerType deficit_delta) {
		nodedata.deficit_delta = deficit_delta;
	}
#endif

private:
	NodeDataType &nodedata;

};


class DefaultEdgeAccessor {

public:
	DefaultEdgeAccessor(NetworkType &N, const DefaultEdgeIterator &edge_it) : N(N), head(N.G.heads[std::abs(*edge_it)]), tail(N.G.tails[std::abs(*edge_it)]), arcdata(N.arcdata[std::abs(*edge_it)]) {}

	auto get_head() const -> DefaultNodeIterator {
		return DefaultNodeIterator(head);
	}

	auto get_tail() const -> DefaultNodeIterator {
		return DefaultNodeIterator(tail);
	}

	auto increase_flow(RationalType amount) -> RationalType {
		auto delta = amount > 0 ?
			std::min(arcdata.capacity - arcdata.xlower, amount) :
			-std::min(arcdata.xlower, -amount);
		arcdata.xlower += delta;
		return delta;
	}

	auto increase_flow_to_capacity() -> RationalType {
		assert(arcdata.xlower <= arcdata.capacity);
		auto flow_incr = arcdata.capacity - arcdata.xlower;
		arcdata.xlower = arcdata.capacity;
		return flow_incr;
	}

	auto can_increase_flow() const -> bool {
#ifdef RESP_CAP
		return true;
#else
		return arcdata.xlower < arcdata.capacity;
#endif
	}

#ifdef NEW_OPTIMIZATION
	auto get_capacity() const -> IntegerType {
		return arcdata.capacity;
	}

	auto get_residual_capacity() const -> IntegerType {
		return arcdata.capacity - arcdata.xlower;
	}

	auto get_residual_backward_capacity() const -> IntegerType {
		return arcdata.xlower;
	}

	auto decrease_flow_to_zero() -> RationalType {
		assert(arcdata.xlower <= arcdata.capacity);
		auto flow_incr = -arcdata.xlower;
		arcdata.xlower = 0;
		return flow_incr;
	}
#endif

	auto can_decrease_flow() const -> bool {
		return arcdata.xlower > 0;
	}

	auto get_forward_cost() const -> RationalType {
		return arcdata.cost;
	}

	auto get_backward_cost() const -> RationalType {
		return -get_forward_cost();
	}

#ifdef RESTORE_BALANCED_NODES
	// TODO: why is the flow a rational type but deficit is integer type?
	auto get_flow_delta() const -> RationalType {
		return arcdata.flow_delta;
	}

	void set_flow_delta(RationalType flow_delta) {
		arcdata.flow_delta = flow_delta;
	}
#endif

#ifdef RESP_CAP
	auto is_transformed() const -> bool {
		return arcdata.transformed;
	}

	auto can_decrease_upper_flow() const -> bool {
		return arcdata.xupper > 0;
	}

	auto increase_upper_flow(RationalType amount) -> RationalType {
		auto delta = amount > 0 ? amount : -std::min(arcdata.xupper, -amount);
		arcdata.xupper += delta;
		return delta;
	}

	auto get_potential() const -> RationalType {
		return arcdata.potential;
	}

	void set_potential(RationalType potential) {
		arcdata.potential = potential;
	}

	auto get_deficit() const -> IntegerType {
		return arcdata.deficit;
	}

	void set_deficit(IntegerType deficit) {
		arcdata.deficit = deficit;
	}

	auto get_depth() const -> int {
		return arcdata.depth;
	}

	void set_depth(int depth) {
		arcdata.depth = depth;
	}

	auto is_visited() const -> bool {
		return arcdata.visited;
	}

	void set_visited(bool visited) {
		arcdata.visited = visited;
	}

#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
	auto is_lower_in_tree() const -> bool {
		return arcdata.lower_in_tree;
	}

	void set_lower_in_tree(bool lower_in_tree) {
		arcdata.lower_in_tree = lower_in_tree;
	}

	auto is_upper_in_tree() const -> bool {
		return arcdata.upper_in_tree;
	}

	void set_upper_in_tree(bool upper_in_tree) {
		arcdata.upper_in_tree = upper_in_tree;
	}

	auto decrease_flow_to_zero() -> RationalType {
		assert(arcdata.xlower >= 0);
		auto flow_incr = -arcdata.xlower;
		arcdata.xlower = 0;
		return flow_incr;
	}

	auto increase_upper_flow_to_capacity() -> RationalType {
		assert(arcdata.xupper <= arcdata.capacity);
		auto flow_incr = arcdata.capacity - arcdata.xupper;
		arcdata.xupper = arcdata.capacity;
		return flow_incr;
	}
#endif
#endif

protected:
	NetworkType &N;
	node head;
	node tail;
	ArcDataType &arcdata;

};


class DefaultIncidentEdgesIterator {

public:
	DefaultIncidentEdgesIterator() = delete;
	DefaultIncidentEdgesIterator(const DefaultIncidentEdgesIterator &) = default;

	DefaultIncidentEdgesIterator(const NetworkType &N, DefaultNodeIterator node_it) : N(N), v(*node_it), pos(0) {}
	DefaultIncidentEdgesIterator(const NetworkType &N, DefaultNodeIterator node_it, size_t pos) : N(N), v(*node_it), pos(pos) {}

	auto operator++() -> DefaultIncidentEdgesIterator & {
		++pos;
		return *this;
	};

	auto operator*() const -> DefaultEdgeIterator {
		return DefaultEdgeIterator(std::abs(N.G.incident_edges[v][pos]));
	};

	auto operator==(const DefaultIncidentEdgesIterator &other) const -> bool {
		return v == other.v && pos == other.pos;
	}

	auto operator!=(const DefaultIncidentEdgesIterator &other) const -> bool {
		return !operator==(other);
	}

private:
	const NetworkType &N;
	const node v;
	size_t pos;

};


template<>
void initialize_internal_data_structures(NetworkType &N) {
	for (auto a = 1u; a <= N.G.no_of_edges; ++a) {
		// saturate negative cost edges
		if (N.arcdata[a].cost < 0) {
			N.arcdata[a].xlower = N.arcdata[a].capacity;
			N.nodedata[N.G.tails[a]].deficit += N.arcdata[a].capacity;
			N.nodedata[N.G.heads[a]].deficit -= N.arcdata[a].capacity;
		}
	}
};

template<>
void initialize_deficit(NetworkType &N) {
	for (node v = 1u; v <= N.G.no_of_vertices; ++v)
		N.nodedata[v].deficit += N.nodedata[v].demand;
};

template<>
auto get_number_of_nodes(const NetworkType &N) -> unsigned int {
	return N.G.no_of_vertices;
};


template<>
auto begin_incident_edges(const NetworkType &N, const DefaultNodeIterator &node) -> DefaultIncidentEdgesIterator {
	return DefaultIncidentEdgesIterator(N, node);
};

template<>
auto end_incident_edges(const NetworkType &N, const DefaultNodeIterator &node) -> DefaultIncidentEdgesIterator {
	return DefaultIncidentEdgesIterator(N, node, N.G.incident_edges[*node].size());
};

template<>
#ifndef NDEBUG
void run_debug_checks(const NetworkType &N) {

#ifndef INFINITE_CAPACITY
#define INFINITE_CAPACITY std::numeric_limits<decltype(N.arcdata.front().capacity)>::max()
#endif

	for (auto a = 1u; a < N.G.no_of_edges; ++a) {
		const auto &arcdata = N.arcdata[a];
		const auto &i = N.G.tails[a];
		const auto &j = N.G.heads[a];

		// check primal feasibility
		assert(arcdata.xlower <= arcdata.capacity && "capacity constraints check failed");
#ifdef RESP_CAP
		if (arcdata.transformed)
			assert(arcdata.xlower + arcdata.xupper == arcdata.capacity && "primal feasibility check failed");
#endif

		// check dual feasibility
		assert((arcdata.capacity != INFINITE_CAPACITY || arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential >= 0) && "dual feasibility check failed");

		// check complementary slackness
		assert((arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential <= 0 || arcdata.xlower == 0) && "complementary slackness check 1 failed");
		assert((arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential >= 0 || arcdata.xlower == arcdata.capacity) && "complementary slackness check 2 failed");
		assert((!(arcdata.xlower > 0 && arcdata.xlower < arcdata.capacity) || arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential == 0) && "complementary slackness check 3 failed");
	}
}
#else
void run_debug_checks(const NetworkType &) {}
#endif



template<>
auto compute_dual_objective_value(const NetworkType &N) -> IntegerType {
	auto dual_objective_value = static_cast<IntegerType>(0);

	for (auto n = 1u; n <= N.G.no_of_vertices; ++n)
		dual_objective_value += N.nodedata[n].demand * N.nodedata[n].potential;

	for (auto a = 1u; a <= N.G.no_of_edges; ++a)
		dual_objective_value += N.arcdata[a].capacity * std::min(static_cast<IntegerType>(0), N.arcdata[a].cost + N.nodedata[N.G.tails[a]].potential - N.nodedata[N.G.heads[a]].potential);

	return dual_objective_value;
}

template<>
auto compute_dual_objective_value_change(const NetworkType &N, const IntegerType& Delta) -> IntegerType {
	auto dual_objective_value = static_cast<IntegerType>(0);

//	cout << "compute dual objective value change" << endl;
	
	for (auto n = 1u; n <= N.G.no_of_vertices; ++n) {
		const auto &nodedata = N.nodedata[n];
		if( !nodedata.visited ) {
//			cout << n << ": " << nodedata.demand << " * " << Delta << " = " << nodedata.demand * Delta << endl;
			dual_objective_value += nodedata.demand * Delta;
		}
	}
	for (auto a = 1u; a <= N.G.no_of_edges; ++a) {
		const node& v = N.G.tails[a];
		const node& w = N.G.heads[a];
//		cout << N.G.tails[a] << " -> " << N.G.heads[a] << ": ";
		if( N.nodedata[v].visited ) {
			if( !N.nodedata[w].visited ) {
				const auto diff = N.arcdata[a].capacity * std::min(static_cast<IntegerType>(0), N.arcdata[a].cost + N.nodedata[v].potential - N.nodedata[w].potential - Delta);
// 				cout << N.arcdata[a].capacity << " * min{ 0, " << N.arcdata[a].cost << " + " << N.nodedata[v].potential << " - " << N.nodedata[w].potential << " - " << Delta << " } = " << diff << endl;
				dual_objective_value += diff;
			} else {
//				cout << endl;
			}
		} else {
			if( N.nodedata[w].visited ) {
				const auto diff = N.arcdata[a].capacity * std::min(static_cast<IntegerType>(0), N.arcdata[a].cost + N.nodedata[v].potential - N.nodedata[w].potential + Delta);
// 				cout << N.arcdata[a].capacity << " * min{ 0, " << N.arcdata[a].cost << " + " << N.nodedata[v].potential << " - " << N.nodedata[w].potential << " + " << Delta << " } = " << diff << endl;
				dual_objective_value += diff;
			} else {
//				cout << endl;
			}
		}
	}
	
	return dual_objective_value;
}


template<>
auto compute_demand_s(const NetworkType &N) -> IntegerType {
	auto demand_s = static_cast<IntegerType>(0);

	for (auto n = 1u; n <= N.G.no_of_vertices; ++n) {
		const auto &nodedata = N.nodedata[n];
		if (nodedata.visited)
			demand_s += nodedata.demand;
	}

	return demand_s;
}


#ifdef RESP_CAP

template<>
void initialize_internal_data_structures_rc(NetworkType &N) {
	initialize_internal_data_structures(N);
};

template<>
void initialize_deficit_rc(NetworkType &N) {
	initialize_deficit(N);

#ifndef INFINITE_CAPACITY
#define INFINITE_CAPACITY std::numeric_limits<decltype(N.arcdata.front().capacity)>::max()
#endif

	for (auto a = 1u; a <= N.G.no_of_edges; ++a) {
		auto &arcdata = N.arcdata[a];
		arcdata.transformed = arcdata.capacity != INFINITE_CAPACITY;
		if (arcdata.transformed)
			arcdata.xupper = arcdata.capacity - arcdata.xlower;
	}
};

#endif

