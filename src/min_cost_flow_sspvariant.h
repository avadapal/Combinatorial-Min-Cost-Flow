#pragma once

#include "min_cost_flow_sspvariant_queue.h"
#include "hybrid_queue.h"

#include <utility>
#include <vector>
#include <queue>

#include <boost/range/adaptor/reversed.hpp>


/* illegal combinations of preprocessor directives */

#if START_NODE_SELECTION_METHOD == 2 && !defined(RESTORE_BALANCED_NODES)
static_assert(false, "If START_NODE_SELECTION_METHOD is set to 2, RESTORE_BALANCED_NODES must be used.");
#endif

#if START_NODE_SELECTION_METHOD == 1 && defined(STOP_IF_DEFICIT_IS_ZERO)
static_assert(false, "If START_NODE_SELECTION_METHOD is set to 1, STOP_IF_DEFICIT_IS_ZERO must not be used.");
#endif

#if START_NODE_SELECTION_METHOD == 1 && defined(STOP_IF_DEFICIT_SIGN_CHANGES)
static_assert(false, "If START_NODE_SELECTION_METHOD is set to 1, STOP_IF_DEFICIT_SIGN_CHANGES must not be used.");
#endif

#if START_NODE_SELECTION_METHOD == 1 && defined(STOP_IF_NODE_DEFICIT_SIGN_CHANGES)
static_assert(false, "If START_NODE_SELECTION_METHOD is set to 1, STOP_IF_NODE_DEFICIT_SIGN_CHANGES must not be used.");
#endif

#if defined(ONLY_IF_PREVIOUSLY_BALANCED) && !defined(RESTORE_BALANCED_NODES)
static_assert(false, "ONLY_IF_PREVIOUSLY_BALANCED can only be used in conjunction with RESTORE_BALANCED_NODES.");
#endif

#if defined(NEW_OPTIMIZATION) && defined(START_FROM_TREE)
static_assert(false, "NEW_OPTIMIZATION is (currently) incompatible with START_FROM_TREE");
#endif

#if defined(NEW_OPTIMIZATION) && defined(HYBRID_QUEUE)
static_assert(false, "NEW_OPTIMIZATION is incompatible with HYBRID_QUEUE");
#endif


/* functions that must be implemented depending on the graph representation */

template<typename Network>
void initialize_internal_data_structures(Network &N);


template<typename Network>
void initialize_deficit(Network &N);


template<typename Network>
void initialize_internal_data_structures_rc(Network &N);


template<typename Network>
void initialize_deficit_rc(Network &N);


template<typename Network>
auto get_number_of_nodes(const Network &N) -> unsigned int;


template<typename Network>
void run_debug_checks(const Network &N);


template<typename Network, typename NumberType>
auto compute_dual_objective_value(const Network &N) -> NumberType;

template<typename Network, typename NumberType>
auto compute_dual_objective_value_change(const Network &N, const NumberType& Delta) -> NumberType;

/*
template<typename Network, typename NumberType>
auto compute_delta_dual_objective_value(const Network &N, NumberType delta) -> NumberType;
*/

template<typename Network, typename NumberType>
auto compute_demand_s(const Network &N)->NumberType;


template<
	typename Network,
	typename NodeIterator,
	typename IncidentEdgesIterator
>
auto begin_incident_edges(const Network &N, const NodeIterator &node) -> IncidentEdgesIterator;


template<
	typename Network,
	typename NodeIterator,
	typename IncidentEdgesIterator
>
auto end_incident_edges(const Network &N, const NodeIterator &node) -> IncidentEdgesIterator;


/* generic functions */

template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor
>
auto get_node_with_most_absolute_deficit(Network &N, const NodeIterator &nbegin, const NodeIterator &nend) -> NodeIterator {
	auto n = nend;
	auto deficit_n = static_cast<decltype(NodeAccessor(N, nbegin).get_deficit())>(0);

	for (auto node_it = nbegin; node_it != nend; ++node_it) {
		const auto node_accessor = NodeAccessor(N, node_it);
		const auto abs_deficit = std::abs(node_accessor.get_deficit());
		if (abs_deficit > deficit_n) {
			n = node_it;
			deficit_n = abs_deficit;
		}
	}

	return n;
}

template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor
>
#if !defined(START_NODE_SELECTION_METHOD) || START_NODE_SELECTION_METHOD == 0
auto get_start_node(Network &N, const NodeIterator &nbegin, const NodeIterator &nend, const NodeIterator &) -> NodeIterator {
	// choose the one with the most absolute deficit
	return get_node_with_most_absolute_deficit<Network, NodeIterator, NodeAccessor>(N, nbegin, nend);
#elif START_NODE_SELECTION_METHOD == 1
auto get_start_node(Network &, const NodeIterator &nbegin, const NodeIterator &, const NodeIterator &) -> NodeIterator {
	// always choose the first node
	return nbegin;
#elif START_NODE_SELECTION_METHOD == 2
auto get_start_node(Network &N, const NodeIterator &nbegin, const NodeIterator &nend, const NodeIterator &previous) -> NodeIterator {
	// choose the one from last iteration, or the next imbalanced node if the previous one is now balanced
#ifndef NDEBUG
	if (previous != nend)
		for (auto node_it = nbegin; node_it != previous; ++node_it)
			assert(NodeAccessor(N, node_it).get_deficit() == 0);
#endif
	if (previous != nend && NodeAccessor(N, previous).get_deficit() != 0)
		return previous;
	auto next = previous;
	while (next != nend && NodeAccessor(N, next).get_deficit() == 0)
		++next;
	return next;
#else
auto get_start_node(Network &N, const NodeIterator &nbegin, const NodeIterator &nend, const NodeIterator &previous) -> NodeIterator {
	// choose the one from last iteration, or the one with the most absolute deficit if the previous one is now balanced
	static_assert(START_NODE_SELECTION_METHOD == 3, "START_NODE_SELECTION_METHOD must be one of:\n    0 (choose the node with the most absolute deficit)\n    1 (always choose the first node)\n    2 (choose the one from last iteration, or the next imbalanced node if the previous one is now balanced)/n    3 (choose the one from last iteration, or the one with the most absolute deficit if the previous one is now balanced)");
	if (previous != nend && NodeAccessor(N, previous).get_deficit() != 0)
		return previous;
	return get_node_with_most_absolute_deficit<Network, NodeIterator, NodeAccessor>(N, nbegin, nend);
#endif
}

template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename QueueType,
#if defined(STOP_IF_DEFICIT_CHANGES_SIGN) || defined(STOP_IF_NODE_DEFICIT_CHANGES_SIGN)
	typename PotentialType,
	typename DeficitType
>
void update_s_in_s_out(Network &N, const NodeIterator &v, QueueType &queue, PotentialType last_delta, DeficitType deficit_s) {
#else
	typename PotentialType
>
void update_s_in_s_out(Network &N, const NodeIterator &v, QueueType &queue, PotentialType last_delta) {
#endif
	// TODO: last_delta is only used for debugging and may be removed
	for (auto incident_edge_it = begin_incident_edges<Network, NodeIterator, IncidentEdgesIterator>(N, v); incident_edge_it != end_incident_edges<Network, NodeIterator, IncidentEdgesIterator>(N, v); ++incident_edge_it) {
		auto edge = *incident_edge_it;
		auto edge_accessor = EdgeAccessor(N, edge);
		bool is_outgoing = edge_accessor.get_tail() == v;
		auto w = is_outgoing ? edge_accessor.get_head() : edge_accessor.get_tail();

		auto v_accessor = NodeAccessor(N, v);
		auto w_accessor = NodeAccessor(N, w);

		auto potential_v = v_accessor.get_potential();
		auto potential_w = w_accessor.get_potential();

		if (w_accessor.is_visited())
			continue;

		if (is_outgoing) {
#if NEW_OPTIMIZATION
			auto key = edge_accessor.get_forward_cost() + potential_v - potential_w;
			queue.push_back({key, edge});
#else
#if defined(STOP_IF_DEFICIT_CHANGES_SIGN) || defined(STOP_IF_NODE_DEFICIT_CHANGES_SIGN)
			if (deficit_s < 0 && edge_accessor.can_increase_flow()) {
#else
			if (edge_accessor.can_increase_flow()) {
#endif
				auto key = edge_accessor.get_forward_cost() + potential_v - potential_w;
				assert(key - last_delta >= 0);
				queue.push({key, edge}, true);
			}
#if defined(STOP_IF_DEFICIT_CHANGES_SIGN) || defined(STOP_IF_NODE_DEFICIT_CHANGES_SIGN)
			if (deficit_s > 0 && edge_accessor.can_decrease_flow()) {
#else
			if (edge_accessor.can_decrease_flow()) {
#endif
				auto key = edge_accessor.get_backward_cost() + potential_w - potential_v;
				assert(key + last_delta >= 0);
				queue.push({key, edge}, false);
			}
#endif
		} else {
#if NEW_OPTIMIZATION
			auto key = -(edge_accessor.get_forward_cost() + potential_w - potential_v);
			queue.push_back({key, edge});
#else
#if defined(STOP_IF_DEFICIT_CHANGES_SIGN) || defined(STOP_IF_NODE_DEFICIT_CHANGES_SIGN)
			if (deficit_s > 0 && edge_accessor.can_increase_flow()) {
#else
			if (edge_accessor.can_increase_flow()) {
#endif
				auto key = edge_accessor.get_forward_cost() + potential_w - potential_v;
				assert(key + last_delta >= 0);
				queue.push({key, edge}, false);
			}
#if defined(STOP_IF_DEFICIT_CHANGES_SIGN) || defined(STOP_IF_NODE_DEFICIT_CHANGES_SIGN)
			if (deficit_s < 0 && edge_accessor.can_decrease_flow()) {
#else
			if (edge_accessor.can_decrease_flow()) {
#endif
				auto key = edge_accessor.get_backward_cost() + potential_v - potential_w;
				assert(key - last_delta >= 0);
				queue.push({key, edge}, true);
			}
#endif
		}
	}
}


template<
	typename Network,
	typename NodeIterator,
#ifdef START_FROM_TREE
	typename EdgeIterator,
#endif
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename DeficitType,
	typename QueueType,
	typename TreeType
>
#ifdef START_FROM_TREE
auto explore_node(Network &N, const NodeIterator &nend, typename QueueType::value_type top_element, QueueType &queue, TreeType &tree, std::vector<std::tuple<NodeIterator, bool, EdgeIterator, int>> &nodes_in_visited_order) -> NodeIterator {
#else
auto explore_node(Network &N, const NodeIterator &nend, typename QueueType::value_type top_element, QueueType &queue, TreeType &tree) -> NodeIterator {
#endif
//  cout << "out = " << out << endl;
	auto delta = top_element.first;
	auto a_hat = top_element.second;

	auto edge_accessor = EdgeAccessor(N, a_hat);

	auto v = edge_accessor.get_tail();
	auto w = edge_accessor.get_head();

	auto v_accessor = NodeAccessor(N, v);
	auto w_accessor = NodeAccessor(N, w);

	if (v_accessor.is_visited() && w_accessor.is_visited())
		return nend;

	auto &new_node = v_accessor.is_visited() ? w : v;

	auto &old_accessor = v_accessor.is_visited() ? v_accessor : w_accessor;
	auto &new_accessor = v_accessor.is_visited() ? w_accessor : v_accessor;

	auto depth = old_accessor.get_depth() + 1;
	new_accessor.set_depth(depth);
	new_accessor.set_visited(true);
	new_accessor.set_potential(new_accessor.get_potential() + delta);


#ifdef START_FROM_TREE
	nodes_in_visited_order.push_back(std::make_tuple(new_node, out, a_hat, depth - 1));
#endif

	if (tree.size() < static_cast<decltype(tree.size())>(depth))
		tree.resize(tree.size() + 1);
	tree[depth - 1].push_back({new_node, a_hat});

	assert(edge_accessor.get_forward_cost() + v_accessor.get_potential() - w_accessor.get_potential() == 0
		|| edge_accessor.get_backward_cost() + w_accessor.get_potential() - v_accessor.get_potential() == 0);

	update_s_in_s_out<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, QueueType, decltype(v_accessor.get_potential())>(N, new_node, queue, delta);

	return new_node;
}



template<
	typename Network,
	typename EdgeIterator,
	typename NodeAccessor,
	typename EdgeAccessor
>
void update_x_by_slack(Network &N, const EdgeIterator &ebegin, const EdgeIterator &eend) {
	for (auto edge_it = ebegin; edge_it != eend; ++edge_it) {
		auto edge_accessor = EdgeAccessor(N, edge_it);
		auto v_accessor = NodeAccessor(N, edge_accessor.get_tail());
		auto w_accessor = NodeAccessor(N, edge_accessor.get_head());
		const auto reduced_cost = edge_accessor.get_forward_cost() + v_accessor.get_potential() - w_accessor.get_potential();

		if (reduced_cost < 0) {
			const auto flow_incr = edge_accessor.increase_flow_to_capacity();
			v_accessor.set_deficit(v_accessor.get_deficit() + flow_incr);
			w_accessor.set_deficit(w_accessor.get_deficit() - flow_incr);
		} else if (reduced_cost > 0) {
			const auto flow_incr = edge_accessor.decrease_flow_to_zero();
			v_accessor.set_deficit(v_accessor.get_deficit() + flow_incr);
			w_accessor.set_deficit(w_accessor.get_deficit() - flow_incr);
		}
	}
}


#ifdef RESTORE_BALANCED_NODES
template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename TreeType
>
void update_x_by_tree(Network &N, const TreeType &tree, const NodeIterator &nbegin, const NodeIterator &nend) {
	for (auto node_it = nbegin; node_it != nend; ++node_it) {
		auto node_accessor = NodeAccessor(N, node_it);
		node_accessor.set_deficit_delta(0);
	}
#else
template<
	typename Network,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename TreeType
>
void update_x_by_tree(Network &N, const TreeType &tree) {
#endif
	for (const auto &depth : boost::adaptors::reverse(tree)) {
		for (const auto &tree_node : depth) {
			// the node v has been explored over the edge a
#ifndef NDEBUG
			// TODO: maybe don't store the node in the tree as it is only used in debugging
			auto v = tree_node.first;
#endif
			auto a = tree_node.second;

			auto edge_accessor = EdgeAccessor(N, a);

			assert(((edge_accessor.get_head() == v) != (edge_accessor.get_tail() == v)) && "the node must be either the edge's head or its tail");
			
			auto i = edge_accessor.get_tail();
			auto j = edge_accessor.get_head();

			auto i_accessor = NodeAccessor(N, i);
			auto j_accessor = NodeAccessor(N, j);

			assert(edge_accessor.get_forward_cost() + i_accessor.get_potential() - j_accessor.get_potential() == 0);

			auto flow_incr = i_accessor.get_depth() > j_accessor.get_depth() ? -i_accessor.get_deficit() : j_accessor.get_deficit();

			// only update edges with slack 0
			if (!((flow_incr < 0 && edge_accessor.get_backward_cost() + j_accessor.get_potential() - i_accessor.get_potential() == 0)
				|| (flow_incr > 0 && edge_accessor.get_forward_cost() + i_accessor.get_potential() - j_accessor.get_potential() == 0))) {
#ifdef RESTORE_BALANCED_NODES
				edge_accessor.set_flow_delta(0);
#endif
				continue;
			}

			flow_incr = edge_accessor.increase_flow(flow_incr);

			i_accessor.set_deficit(i_accessor.get_deficit() + flow_incr);
			j_accessor.set_deficit(j_accessor.get_deficit() - flow_incr);

#ifdef RESTORE_BALANCED_NODES
			edge_accessor.set_flow_delta(flow_incr);
			i_accessor.set_deficit_delta(i_accessor.get_deficit_delta() + flow_incr);
			j_accessor.set_deficit_delta(j_accessor.get_deficit_delta() - flow_incr);
#endif
		}
	}

#ifdef RESTORE_BALANCED_NODES

	for (const auto &depth : tree) {
		for (const auto &tree_node : depth) {
			auto a = tree_node.second;
			auto edge_accessor = EdgeAccessor(N, a);

			auto i = edge_accessor.get_tail();
			auto j = edge_accessor.get_head();

			auto i_accessor = NodeAccessor(N, i);
			auto j_accessor = NodeAccessor(N, j);

			using DeficitType = decltype(i_accessor.get_deficit());

			if (i_accessor.get_depth() > j_accessor.get_depth()) {
				// i is the child node and j the parent node

				auto current_deficit = j_accessor.get_deficit();
				auto previous_deficit = current_deficit - j_accessor.get_deficit_delta();

#ifdef ONLY_IF_PREVIOUSLY_BALANCED
				if (!(previous_deficit == 0 && current_deficit != 0))
					continue;

				auto target_deficit = 0;
#else
				auto target_deficit = previous_deficit < 0 ?
					std::min(static_cast<DeficitType>(0), std::max(previous_deficit, current_deficit)) :
					std::max(static_cast<DeficitType>(0), std::min(previous_deficit, current_deficit));
#endif

				assert(std::abs(target_deficit) <= std::abs(previous_deficit));
				assert(std::abs(target_deficit) <= std::abs(current_deficit));

				if (current_deficit == target_deficit)
					// nothing changed
					continue;

				auto flow_incr = -(target_deficit - current_deficit);

				if (!((flow_incr > 0 && edge_accessor.get_flow_delta() < 0) || (flow_incr < 0 && edge_accessor.get_flow_delta() > 0)))
					// flow was updated in the other direction
					continue;

				flow_incr = flow_incr > 0 ?
					std::min(flow_incr, -edge_accessor.get_flow_delta()) :
					std::max(flow_incr, -edge_accessor.get_flow_delta());

#ifndef NDEBUG
				auto actual_flow_incr = edge_accessor.increase_flow(flow_incr);
				assert(actual_flow_incr == flow_incr);
#else
				edge_accessor.increase_flow(flow_incr);
#endif

				i_accessor.set_deficit(i_accessor.get_deficit() + flow_incr);
				i_accessor.set_deficit_delta(i_accessor.get_deficit_delta() + flow_incr);
				j_accessor.set_deficit(j_accessor.get_deficit() - flow_incr);
				j_accessor.set_deficit_delta(j_accessor.get_deficit_delta() - flow_incr);

				assert(std::abs(j_accessor.get_deficit()) <= std::abs(current_deficit));
			} else {
				// j is the child node and i the parent node

				auto current_deficit = i_accessor.get_deficit();
				auto previous_deficit = current_deficit - i_accessor.get_deficit_delta();

#ifdef ONLY_IF_PREVIOUSLY_BALANCED
				if (!(previous_deficit == 0 && current_deficit != 0))
					continue;

				auto target_deficit = 0;
#else
				auto target_deficit = previous_deficit < 0 ?
					std::min(static_cast<DeficitType>(0), std::max(previous_deficit, current_deficit)) :
					std::max(static_cast<DeficitType>(0), std::min(previous_deficit, current_deficit));
#endif

				assert(std::abs(target_deficit) <= std::abs(previous_deficit));
				assert(std::abs(target_deficit) <= std::abs(current_deficit));

				if (current_deficit == target_deficit)
					// nothing changed
					continue;

				auto flow_incr = target_deficit - current_deficit;

				if (!((flow_incr > 0 && edge_accessor.get_flow_delta() < 0) || (flow_incr < 0 && edge_accessor.get_flow_delta() > 0)))
					// flow was updated in the other direction
					continue;

				flow_incr = flow_incr > 0 ?
					std::min(flow_incr, -edge_accessor.get_flow_delta()) :
					std::max(flow_incr, -edge_accessor.get_flow_delta());

#ifndef NDEBUG
				auto actual_flow_incr = edge_accessor.increase_flow(flow_incr);
				assert(actual_flow_incr == flow_incr);
#else
				edge_accessor.increase_flow(flow_incr);
#endif

				i_accessor.set_deficit(i_accessor.get_deficit() + flow_incr);
				i_accessor.set_deficit_delta(i_accessor.get_deficit_delta() + flow_incr);
				j_accessor.set_deficit(j_accessor.get_deficit() - flow_incr);
				j_accessor.set_deficit_delta(j_accessor.get_deficit_delta() - flow_incr);

				assert(std::abs(i_accessor.get_deficit()) <= std::abs(current_deficit));
			}
		}
	}

#ifndef NDEBUG
	for (const auto &depth : tree) {
		for (const auto &tree_node : depth) {
			auto a = tree_node.second;
			auto edge_accessor = EdgeAccessor(N, a);

			auto i = edge_accessor.get_tail();
			auto j = edge_accessor.get_head();

			auto i_accessor = NodeAccessor(N, i);
			auto j_accessor = NodeAccessor(N, j);

			auto old_deficit = i_accessor.get_depth() > j_accessor.get_depth() ?
				j_accessor.get_deficit() - j_accessor.get_deficit_delta() :
				i_accessor.get_deficit() - i_accessor.get_deficit_delta();
			auto new_deficit = i_accessor.get_depth() > j_accessor.get_depth() ?
				j_accessor.get_deficit() :
				i_accessor.get_deficit();

			assert(old_deficit > 0 || (new_deficit <= 0 && new_deficit >= old_deficit));
			assert(old_deficit < 0 || (new_deficit >= 0 && new_deficit <= old_deficit));
		}
	}
#endif
#endif

}


template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor,
	typename DeficitType
> auto compute_deficit_1(Network &N, const NodeIterator &nbegin, const NodeIterator &nend) -> DeficitType {
	auto deficit_1 = static_cast<DeficitType>(0);

	for (auto node_it = nbegin; node_it != nend; ++node_it) {
		auto node_accessor = NodeAccessor(N, node_it);
		deficit_1 += std::abs(node_accessor.get_deficit());
	}

	return deficit_1;
}


template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename DeficitType,
	typename QueueType
> auto build_initial_s(Network &N, const NodeIterator &nbegin, const NodeIterator &nend, NodeIterator &start_node, QueueType &queue, bool is_first_iteration) -> NodeIterator {
	start_node = get_start_node<Network, NodeIterator, NodeAccessor>(N, nbegin, nend, start_node);

#if START_NODE_SELECTION_METHOD == 1
	if (compute_deficit_1<Network, NodeIterator, NodeAccessor, DeficitType>(N, nbegin, nend) == 0 && !is_first_iteration)
		return nend;
#else
	if (start_node == nend) {
		if (is_first_iteration) {
			start_node = nbegin;
		} else {
			return nend;
		}
	}
#endif

	auto start_node_accessor = NodeAccessor(N, start_node);
	start_node_accessor.set_visited(true);
	start_node_accessor.set_depth(0);

	update_s_in_s_out<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, QueueType, decltype(start_node_accessor.get_potential())>(N, start_node, queue, 0);

	return start_node;
}


template<
	typename Network,
	typename NodeIterator,
	typename EdgeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename DeficitType,
	typename QueueType,
	typename TreeType
> auto build_initial_s_as_tree(Network &N, std::vector<std::tuple<NodeIterator, bool, EdgeIterator, int>> &nodes_in_visited_order, const NodeIterator &nbegin, const NodeIterator &nend, QueueType &queue, TreeType &tree) -> std::pair<unsigned int, DeficitType> {
	// TODO: cleanup
	
	// first iteration
	if (nodes_in_visited_order.empty()) {
		auto v_start = nbegin;
		nodes_in_visited_order.push_back(std::make_tuple(v_start, true, *begin_incident_edges<Network, NodeIterator, IncidentEdgesIterator>(N, v_start), 0));

		auto v_start_accessor = NodeAccessor(N, v_start);
		v_start_accessor.set_visited(true);
		v_start_accessor.set_depth(0);

		update_s_in_s_out<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, QueueType>(N, v_start, queue);

		return {1, v_start_accessor.get_deficit()};
	}

	assert(!nodes_in_visited_order.empty());
	auto v_start_accessor = NodeAccessor(N, std::get<0>(nodes_in_visited_order.front()));
	v_start_accessor.set_visited(true);
	v_start_accessor.set_depth(0);
	auto deficit_s = v_start_accessor.get_deficit();
	auto num_visited = 1u;

	for (auto v = std::begin(nodes_in_visited_order) + 1; v != std::end(nodes_in_visited_order); ++v) {
		if ((std::get<bool>(*v) != (deficit_s < 0))) {
			nodes_in_visited_order.erase(v, std::end(nodes_in_visited_order));
			break;
		}

		auto v_accessor = NodeAccessor(N, std::get<NodeIterator>(*v));
		auto edge_accessor = EdgeAccessor(N, std::get<EdgeIterator>(*v));

		auto outgoing = edge_accessor.get_head() == std::get<NodeIterator>(*v);
		auto increase_flow = std::get<bool>(*v) == outgoing;
		if ((increase_flow && !edge_accessor.can_increase_flow()) || (!increase_flow && !edge_accessor.can_decrease_flow())) {
			nodes_in_visited_order.erase(v, std::end(nodes_in_visited_order));
			break;
		}

		++num_visited;
		v_accessor.set_visited(true);
		deficit_s += v_accessor.get_deficit();

		auto depth = std::get<int>(*v);
		v_accessor.set_depth(depth);
		if (depth != 0) {
			if (tree.size() < static_cast<decltype(tree.size())>(depth))
				tree.resize(tree.size() + 1);
			tree[depth - 1].push_back({std::get<NodeIterator>(*v), std::get<EdgeIterator>(*v)});
		}
	}

	//std::cout << "kept " << num_visited << " nodes from previous iteration" << std::endl;

	for (const auto &v : nodes_in_visited_order)
		update_s_in_s_out<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, QueueType>(N, std::get<NodeIterator>(v), queue);
	
	return {num_visited, deficit_s};
}


template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor
> void reset_data(Network &N, const NodeIterator &nbegin, const NodeIterator &nend) {
	for (auto node_it = nbegin; node_it != nend; ++node_it) {
		auto node_accessor = NodeAccessor(N, node_it);
		node_accessor.set_visited(false);
		node_accessor.set_depth(-1);
	}
}


template<
	typename Network,
	typename NodeIterator,
	typename EdgeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator
>
#ifndef NEW_OPTIMIZATION
void successive_shortest_path(Network &N, const NodeIterator &nbegin, const NodeIterator &nend) {
#else
void successive_shortest_path(Network &N, const NodeIterator &nbegin, const NodeIterator &nend, const EdgeIterator &ebegin, const EdgeIterator &eend) {
#endif

#if (defined(STOP_IF_DEFICIT_IS_ZERO) || defined(STOP_IF_DEFICIT_SIGN_CHANGES) || defined(STOP_IF_NODE_DEFICIT_SIGN_CHANGES)) && !defined(EXPLORE_ENTIRE_GRAPH_EVERY_X_ITERATIONS)
#define EXPLORE_ENTIRE_GRAPH_EVERY_X_ITERATIONS 0
#endif

	using PotentialType = decltype(NodeAccessor(N, nbegin).get_potential());
	using DeficitType = decltype(NodeAccessor(N, nbegin).get_deficit());

	initialize_internal_data_structures(N);
	initialize_deficit(N);

	auto step_count = 0u;

#ifdef START_FROM_TREE
	auto nodes_in_visited_order = std::vector<std::tuple<NodeIterator, bool, EdgeIterator, int>>();
	nodes_in_visited_order.reserve(get_number_of_nodes(N));
#else
	auto start_node = nend;
#endif

#ifndef NDEBUG
	auto last_dual_objective_value = compute_dual_objective_value<Network, DeficitType>(N);
#endif

	do {
		++step_count;

#ifndef NDEBUG
		run_debug_checks(N);
#endif

		using PriorityQueueElementType = std::pair<PotentialType, EdgeIterator>;

#ifdef NEW_OPTIMIZATION
		auto queue = std::vector<PriorityQueueElementType>();
#else
#ifdef USE_HYBRID_QUEUE
		using InternalQueueType = hybrid_queue<PotentialType, EdgeIterator>;

		auto queue = SSPVariantQueue<InternalQueueType>(InternalQueueType(), InternalQueueType());
#else
		auto compare_first = [](const PriorityQueueElementType &lhs, const PriorityQueueElementType &rhs) {
			return lhs.first > rhs.first;
		};

		using InternalQueueType = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>;
		auto queue = SSPVariantQueue<InternalQueueType>(InternalQueueType(compare_first), InternalQueueType(compare_first));
#endif
#endif

		auto tree = std::vector<std::vector<std::pair<NodeIterator, EdgeIterator>>>();

		reset_data<Network, NodeIterator, NodeAccessor>(N, nbegin, nend);

		unsigned int num_visited;
		DeficitType deficit_s = 0;

#ifndef NDEBUG
		std::cout << "starting step " << step_count;
		std::cout << ", deficit_1 = " << compute_deficit_1<Network, NodeIterator, NodeAccessor, DeficitType>(N, nbegin, nend);
		auto new_dual_objective_value = compute_dual_objective_value<Network, DeficitType>(N);
		std::cout << ", dual objective value = " << new_dual_objective_value << std::endl;
		assert(last_dual_objective_value <= new_dual_objective_value);
		last_dual_objective_value = new_dual_objective_value;
#endif

#ifdef START_FROM_TREE
		// TODO: cleanup
		std::tie(num_visited, deficit_s) = build_initial_s_as_tree<Network, NodeIterator, EdgeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, QueueType>(N, nodes_in_visited_order, nbegin, nend, queue, tree);
		if (num_visited == get_number_of_nodes(N)) {
			auto deficit_1 = static_cast<DeficitType>(0);
			for (auto node_it = nbegin; node_it != nend; ++node_it)
				deficit_1 += std::abs(NodeAccessor(N, node_it).get_deficit());
			// TODO: fix this
			assert(deficit_1 == 0);
			if (deficit_1 == 0)
				break;

			nodes_in_visited_order.erase(std::begin(nodes_in_visited_order) + 1, std::end(nodes_in_visited_order));
			tree.clear();
		}
		if (num_visited == 0)
			break;
#else
		auto v_start = build_initial_s<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, decltype(queue)>(N, nbegin, nend, start_node, queue, step_count == 1u);
		if (v_start == nend)
			break;
		auto v_start_accessor = NodeAccessor(N, v_start);
		num_visited = 1;
		deficit_s = v_start_accessor.get_deficit();
#ifdef NEW_OPTIMIZATION
		auto demand_s = v_start_accessor.get_demand();

		auto queue_compare = [&N, &deficit_s](const PriorityQueueElementType& lhs, const PriorityQueueElementType& rhs) -> bool {
			if (lhs.first != rhs.first)
				return lhs.first < rhs.first;

			const auto lhs_accessor = EdgeAccessor(N, lhs.second);
			const auto rhs_accessor = EdgeAccessor(N, rhs.second);

			const auto lhs_v_accessor = NodeAccessor(N, lhs_accessor.get_tail());
			const auto lhs_w_accessor = NodeAccessor(N, lhs_accessor.get_head());

			const auto rhs_v_accessor = NodeAccessor(N, rhs_accessor.get_tail());
			const auto rhs_w_accessor = NodeAccessor(N, rhs_accessor.get_head());

			assert(lhs_v_accessor.is_visited() != lhs_w_accessor.is_visited());
			assert(rhs_v_accessor.is_visited() != rhs_w_accessor.is_visited());

			const auto lhs_is_outgoing = lhs_v_accessor.is_visited();
			const auto rhs_is_outgoing = rhs_v_accessor.is_visited();

			// tie breaker: residual capacity
			if (deficit_s < 0) {
				// prefer edge over which we can send more flow out of the network
				const auto lhs_residual_capacity = lhs_is_outgoing ?
					lhs_accessor.get_residual_capacity() :
					lhs_accessor.get_residual_backward_capacity();
				const auto rhs_residual_capacity = rhs_is_outgoing ?
					rhs_accessor.get_residual_capacity() :
					rhs_accessor.get_residual_backward_capacity();
				return lhs_residual_capacity < rhs_residual_capacity;
			} else {
				// prefer edge over which we can send more flow into the network
				const auto lhs_residual_capacity = lhs_is_outgoing ?
					lhs_accessor.get_residual_backward_capacity() :
					lhs_accessor.get_residual_capacity();
				const auto rhs_residual_capacity = rhs_is_outgoing ?
					rhs_accessor.get_residual_backward_capacity() :
					rhs_accessor.get_residual_capacity();
				return lhs_residual_capacity < rhs_residual_capacity;
			}
		};

		std::sort(std::begin(queue), std::end(queue), queue_compare);
#endif
#endif
#ifdef STOP_IF_NODE_DEFICIT_SIGN_CHANGES
		const auto v_start_deficit = deficit_s;
		assert(v_start_deficit != 0);
#endif

		while (num_visited < get_number_of_nodes(N) && !queue.empty()) {
	
#if defined(STOP_IF_DEFICIT_IS_ZERO) || defined(STOP_IF_DEFICIT_SIGN_CHANGES) || defined(STOP_IF_NODE_DEFICIT_SIGN_CHANGES)
			auto delta = queue.top(deficit_s).first;
#endif

#ifdef NEW_OPTIMIZATION
			//assert(std::is_sorted(std::begin(queue), std::end(queue), queue_compare));
			// TODO: remove sort stuff after explore_node since it has to be sorted after changing deficit_s anyway
			std::sort(std::begin(queue), std::end(queue), queue_compare);

			if (queue.empty())
				break;

			// build pre- and postfix sums
			auto prefix_sums = std::vector<DeficitType>(queue.size(), 0);
			auto postfix_sums = std::vector<DeficitType>(queue.size(), 0);

			auto prefix_summand = [&N, &queue](std::size_t i) -> DeficitType {
				auto edge_accessor = EdgeAccessor(N, queue[i].second);
				return NodeAccessor(N, edge_accessor.get_head()).is_visited() ? 0 : -edge_accessor.get_capacity();
			};

			auto postfix_summand = [&N, &queue](std::size_t i) -> DeficitType {
				auto edge_accessor = EdgeAccessor(N, queue[i].second);
				return NodeAccessor(N, edge_accessor.get_tail()).is_visited() ? 0 : edge_accessor.get_capacity();
			};

			prefix_sums.front() = prefix_summand(0);
			postfix_sums.back() = postfix_summand(queue.size() - 1);

			for (auto i = 1u; i < queue.size(); ++i)
				prefix_sums[i] = prefix_sums[i - 1] + prefix_summand(i);

			for (int i = queue.size() - 2; i >= 0; --i)
				postfix_sums[i] = postfix_sums[i + 1] + postfix_summand(i);

			// compute diff
			auto max_point = std::end(queue);
			//std::cout << "values:" << std::endl;
			for (auto i = 0u; i < queue.size() - 1; ++i) {
				if (queue[i + 1].first == queue[i].first)
					continue;
				//std::cout << "  " << i << ": " << (-demand_s + prefix_sums[i] + postfix_sums[i + 1]) << std::endl;
				if (-demand_s + prefix_sums[i] + postfix_sums[i + 1] <= 0) {
					max_point = std::begin(queue) + i + 1;
					break;
				}
			}

			auto queue_size = queue.size();

			auto new_node = explore_node<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, decltype(queue), decltype(tree)>(N, nend, *(max_point - 1), queue, tree);
			
			std::sort(std::begin(queue) + queue_size, std::end(queue), queue_compare);

			// remove arcs where both nodes are visited
			auto both_visited = [&N](const PriorityQueueElementType &queue_element) -> bool {
				auto edge_accessor = EdgeAccessor(N, queue_element.second);
				return NodeAccessor(N, edge_accessor.get_tail()).is_visited() && NodeAccessor(N, edge_accessor.get_head()).is_visited();
			};

			auto mid = std::remove_if(std::begin(queue), std::begin(queue) + queue_size, both_visited);
			mid = queue.erase(mid, std::begin(queue) + queue_size);

			assert(std::is_sorted(std::begin(queue), mid, queue_compare));

			std::inplace_merge(std::begin(queue), mid, std::end(queue), queue_compare);
#else
			auto pop_top = [](decltype(queue) &queue, const DeficitType &deficit_s) -> PriorityQueueElementType {
				auto top = queue.top(deficit_s);
				queue.pop(deficit_s);
				return top;
			};

#ifdef START_FROM_TREE
			auto new_node = explore_node<Network, NodeIterator, EdgeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, decltype(queue), decltype(tree)>(N, nend, pop_top(queue, deficit_s), queue, tree, nodes_in_visited_order);
#else
			auto new_node = explore_node<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, decltype(queue), decltype(tree)>(N, nend, pop_top(queue, deficit_s), queue, tree);
#endif

#endif

			if (new_node != nend) {
				auto node_accessor = NodeAccessor(N, new_node);
#ifdef STOP_IF_DEFICIT_SIGN_CHANGES
				auto previous_deficit_s = deficit_s;
#endif
				deficit_s += node_accessor.get_deficit();
#ifdef NEW_OPTIMIZATION
				demand_s += node_accessor.get_demand();
#endif
				++num_visited;
#if defined(STOP_IF_DEFICIT_IS_ZERO) || defined(STOP_IF_DEFICIT_SIGN_CHANGES) || defined(STOP_IF_NODE_DEFICIT_SIGN_CHANGES)
				if (EXPLORE_ENTIRE_GRAPH_EVERY_X_ITERATIONS == 0 || step_count % EXPLORE_ENTIRE_GRAPH_EVERY_X_ITERATIONS != 0) {
#ifdef STOP_IF_DEFICIT_IS_ZERO
					if (step_count > 1 && deficit_s == static_cast<DeficitType>(0)) {
#elif defined(STOP_IF_DEFICIT_SIGN_CHANGES)
					if (step_count > 1 && ((previous_deficit_s > 0 && !(deficit_s > 0)) || (previous_deficit_s < 0 && !(deficit_s < 0)))) {
#else
					if (step_count > 1 && ((node_accessor.get_deficit() > 0 && v_start_deficit < 0) || (node_accessor.get_deficit() < 0 && v_start_deficit > 0))) {

						// remove everything but the path from the start node to the final node from the tree
						tree.erase(std::begin(tree) + node_accessor.get_depth(), std::end(tree));
						auto node_that_remains_in_the_tree = new_node;

						for (auto &depth : boost::adaptors::reverse(tree)) {
							depth.erase(std::remove_if(std::begin(depth), std::end(depth), [&node_that_remains_in_the_tree](const std::pair<NodeIterator, EdgeIterator> &tree_node) {
								return tree_node.first != node_that_remains_in_the_tree;
							}), std::end(depth));

							assert(depth.size() == 1);

							const auto edge_accessor = EdgeAccessor(N, depth.front().second);

							const auto &tail = edge_accessor.get_tail();
							const auto &head = edge_accessor.get_head();

							node_that_remains_in_the_tree = tail == node_that_remains_in_the_tree ? head : tail;
						}
#endif
						for (auto node_it = nbegin; node_it != nend; ++node_it) {
							auto node_accessor = NodeAccessor(N, node_it);
							if (!node_accessor.is_visited())
								node_accessor.set_potential(node_accessor.get_potential() + delta);
						}
						// std::cout << "size of S: " << num_visited << std::endl;
						break;
					}
				}
#endif
			}

		}

#if START_NODE_SELECTION_METHOD == 1 && !defined(STOP_IF_DEFICIT_IS_ZERO) && !defined(STOP_IF_DEFICIT_SIGN_CHANGES) && !defined(STOP_IF_NODE_DEFICIT_SIGN_CHANGES)
		if (num_visited != get_number_of_nodes(N)) {
			std::cerr << "Disconnected graphs are not supported with this configuration." << std::endl;
			exit(1);
		}
#endif

#ifdef NEW_OPTIMIZATION
		update_x_by_slack<Network, EdgeIterator, NodeAccessor, EdgeAccessor>(N, ebegin, eend);
#endif

		
#ifdef RESTORE_BALANCED_NODES
		update_x_by_tree<Network, NodeIterator, NodeAccessor, EdgeAccessor, decltype(tree)>(N, tree, nbegin, nend);
#else
		update_x_by_tree<Network, NodeAccessor, EdgeAccessor, decltype(tree)>(N, tree);
#endif

	} while (true);

#ifndef NDEBUG
	auto deficit_1 = compute_deficit_1<Network, NodeIterator, NodeAccessor, DeficitType>(N, nbegin, nend);
	assert(deficit_1 == 0);
#endif
	std::cout << "found solution after " << (step_count - 1) << ((step_count - 1) == 1 ? " step" : " steps") << ", dual objective value = " << compute_dual_objective_value<Network, DeficitType>(N) << std::endl;

}



template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename QueueType
>
void update_s_in_s_out_node_rc(Network &N, const NodeIterator &v, QueueType &s_in_n, QueueType &s_out_n, QueueType &s_in_a, QueueType &s_out_a) {
	for (auto incident_edge_it = begin_incident_edges<Network, NodeIterator, IncidentEdgesIterator>(N, v); incident_edge_it != end_incident_edges<Network, NodeIterator, IncidentEdgesIterator>(N, v); ++incident_edge_it) {
		auto edge = *incident_edge_it;
		auto edge_accessor = EdgeAccessor(N, edge);
		bool is_outgoing = edge_accessor.get_tail() == v;

		if (edge_accessor.is_transformed()) {
			if (edge_accessor.is_visited())
				continue;

			auto potential_v = NodeAccessor(N, v).get_potential();
			auto potential_a = edge_accessor.get_potential();

			if (is_outgoing) {
				assert(edge_accessor.can_increase_flow());
				s_out_n.push({edge_accessor.get_forward_cost() + potential_v - potential_a, {edge, false}});
				if (edge_accessor.can_decrease_flow())
					s_in_n.push({edge_accessor.get_backward_cost() + potential_a - potential_v, {edge, false}});
			} else {
				s_out_n.push({potential_v - potential_a, {edge, true}});
				if (edge_accessor.can_decrease_upper_flow())
					s_in_n.push({potential_a - potential_v, {edge, true}});
			}
		} else {
			auto w = is_outgoing ? edge_accessor.get_head() : edge_accessor.get_tail();
			auto w_accessor = NodeAccessor(N, w);

			if (w_accessor.is_visited())
				continue;

			auto potential_v = NodeAccessor(N, v).get_potential();
			auto potential_w = w_accessor.get_potential();

			if (is_outgoing) {
				assert(edge_accessor.can_increase_flow());
				s_out_a.push({edge_accessor.get_forward_cost() + potential_v - potential_w, {edge, false}});
				if (edge_accessor.can_decrease_flow())
					s_in_a.push({edge_accessor.get_backward_cost() + potential_w - potential_v, {edge, false}});
			} else {
				assert(edge_accessor.can_increase_flow());
				s_in_a.push({edge_accessor.get_forward_cost() + potential_w - potential_v, {edge, false}});
				if (edge_accessor.can_decrease_flow())
					s_out_a.push({edge_accessor.get_backward_cost() + potential_v - potential_w, {edge, false}});
			}
		}
	}
}


template<
	typename Network,
	typename EdgeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename QueueType
>
void update_s_in_s_out_edge_rc(Network &N, const EdgeIterator &a, QueueType &s_in_a, QueueType &s_out_a) {
	auto edge_accessor = EdgeAccessor(N, a);
	assert(edge_accessor.is_transformed());
	auto i = edge_accessor.get_tail();
	auto j = edge_accessor.get_head();
	auto i_accessor = NodeAccessor(N, i);
	auto j_accessor = NodeAccessor(N, j);
	auto potential_a = edge_accessor.get_potential();

	if (!i_accessor.is_visited()) {
		auto potential_i = i_accessor.get_potential();
		assert(edge_accessor.can_increase_flow());
		s_in_a.push({edge_accessor.get_forward_cost() + potential_i - potential_a, {a, false}});
		if (edge_accessor.can_decrease_flow())
			s_out_a.push({edge_accessor.get_backward_cost() + potential_a - potential_i, {a, false}});
	}
	if (!j_accessor.is_visited()) {
		auto potential_j = j_accessor.get_potential();
		s_in_a.push({potential_j - potential_a, {a, true}});
		if (edge_accessor.can_decrease_upper_flow())
			s_out_a.push({potential_a - potential_j, {a, true}});
	}
}


template<
	typename Network,
	typename NodeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename DeficitType,
	typename QueueType,
	typename TreeType
>
auto explore_node_rc(Network &N, typename QueueType::value_type top_element, QueueType &s_in_n, QueueType &s_out_n, QueueType &s_in_a, QueueType &s_out_a, TreeType &tree, bool out) -> std::pair<bool, DeficitType> {
	auto delta = out ? top_element.first : -top_element.first;
	auto a_hat = top_element.second.first;
	auto upper = top_element.second.second;

	auto edge_accessor = EdgeAccessor(N, a_hat);

	auto v = edge_accessor.get_tail();
	auto w = edge_accessor.get_head();

	auto v_accessor = NodeAccessor(N, v);
	auto w_accessor = NodeAccessor(N, w);

	decltype(v_accessor.get_depth()) depth;

	auto &new_node = edge_accessor.is_transformed() ? (upper ? w : v) : (v_accessor.is_visited() ? w : v);
	auto &new_accessor = edge_accessor.is_transformed() ? (upper ? w_accessor : v_accessor) : (v_accessor.is_visited() ? w_accessor : v_accessor);

	if (new_accessor.is_visited())
		return {false, static_cast<DeficitType>(0)};

	if (edge_accessor.is_transformed()) {
		assert((upper ? v_accessor : w_accessor).is_visited());
#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
		auto &old_accessor = upper ? v_accessor : w_accessor;
		assert(!upper || edge_accessor.is_lower_in_tree());
		assert(upper || edge_accessor.is_upper_in_tree());
		depth = old_accessor.get_depth() + 1;
#else
		depth = edge_accessor.get_depth() + 1;
#endif
	} else {
		auto &old_accessor = v_accessor.is_visited() ? v_accessor : w_accessor;
		assert(old_accessor.is_visited());
		depth = old_accessor.get_depth() + 1;
	}

	if (tree.size() < static_cast<decltype(tree.size())>(depth))
		tree.resize(tree.size() + 1);

#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
	if (upper) {
		edge_accessor.set_upper_in_tree(true);
	} else {
		edge_accessor.set_lower_in_tree(true);
	}
	tree[depth - 1].push_back({new_node, a_hat});
#else
	tree[depth - 1].push_back({new_node, {a_hat, upper}});
#endif

	new_accessor.set_depth(depth);
	new_accessor.set_visited(true);
	new_accessor.set_potential(new_accessor.get_potential() + delta);

#ifndef NDEBUG
	if (edge_accessor.is_transformed()) {
		if (upper) {
			assert(w_accessor.get_potential() - edge_accessor.get_potential() == 0);
		} else {
			assert(edge_accessor.get_forward_cost() + v_accessor.get_potential() - edge_accessor.get_potential() == 0
				|| edge_accessor.get_backward_cost() + edge_accessor.get_potential() - v_accessor.get_potential() == 0);
		}
	} else {
		assert(edge_accessor.get_forward_cost() + v_accessor.get_potential() - w_accessor.get_potential() == 0
			|| edge_accessor.get_backward_cost() + w_accessor.get_potential() - v_accessor.get_potential() == 0);
	}
#endif

	update_s_in_s_out_node_rc<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, QueueType>(N, new_node, s_in_n, s_out_n, s_in_a, s_out_a);

	return {true, new_accessor.get_deficit()};
}


template<
	typename Network,
	typename EdgeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator,
	typename DeficitType,
	typename QueueType,
	typename TreeType
>
auto explore_arc_rc(Network &N, typename QueueType::value_type top_element, QueueType &s_in_a, QueueType &s_out_a, TreeType &tree, bool out) -> std::pair<bool, DeficitType> {
	auto delta = out ? top_element.first : -top_element.first;
	auto a_hat = top_element.second.first;
	auto upper = top_element.second.second;

	auto edge_accessor = EdgeAccessor(N, a_hat);
	assert(edge_accessor.is_transformed());

	if (edge_accessor.is_visited())
		return {false, static_cast<DeficitType>(0)};

	auto parent = upper ? edge_accessor.get_head() : edge_accessor.get_tail();
	auto parent_accessor = NodeAccessor(N, parent);
	assert(parent_accessor.is_visited() && "a_hat is an arc that was added after exploring the parent node");

	edge_accessor.set_visited(true);
	edge_accessor.set_potential(edge_accessor.get_potential() + delta);

#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
	if (upper) {
		edge_accessor.set_upper_in_tree(true);
	} else {
		edge_accessor.set_lower_in_tree(true);
	}
#else
	auto depth = parent_accessor.get_depth() + 1;
	edge_accessor.set_depth(depth);
	edge_accessor.set_visited(true);
	if (tree.size() < static_cast<decltype(tree.size())>(depth))
		tree.resize(tree.size() + 1);
	tree[depth - 1].push_back({parent, {a_hat, upper}});
#endif

#ifndef NDEBUG
	if (upper) {
		assert(edge_accessor.get_potential() - parent_accessor.get_potential() == 0);
	} else {
		assert(edge_accessor.get_forward_cost() + parent_accessor.get_potential() - edge_accessor.get_potential() == 0
			|| edge_accessor.get_backward_cost() + edge_accessor.get_potential() - parent_accessor.get_potential() == 0);
	}
#endif

	update_s_in_s_out_edge_rc<Network, EdgeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, QueueType>(N, a_hat, s_in_a, s_out_a);

	return {true, edge_accessor.get_deficit()};
}


template<
	typename Network,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename TreeType
>
void update_x_by_tree_rc(Network &N, const TreeType &tree) {
	for (const auto &depth : boost::adaptors::reverse(tree)) {
		for (const auto &tree_node : depth) {
			auto v = tree_node.first;
			auto a = tree_node.second.first;

			auto edge_accessor = EdgeAccessor(N, a);

			assert(((edge_accessor.get_head() == v) != (edge_accessor.get_tail() == v)) && "the node must be either the edge's head or its tail");

			if (edge_accessor.is_transformed()) {
				auto upper = tree_node.second.second;
				auto node_accessor = NodeAccessor(N, v);

				auto flow_incr = edge_accessor.get_depth() > node_accessor.get_depth() ? edge_accessor.get_deficit() : -node_accessor.get_deficit();

				if (upper) {
					assert(node_accessor.get_potential() - edge_accessor.get_potential() == 0);
					flow_incr = edge_accessor.increase_upper_flow(flow_incr);
				} else {
					// only update edges with slack 0
					if (!((flow_incr < 0 && edge_accessor.get_backward_cost() + edge_accessor.get_potential() - node_accessor.get_potential() == 0)
						|| (flow_incr > 0 && edge_accessor.get_forward_cost() + node_accessor.get_potential() - edge_accessor.get_potential() == 0)))
						continue;
					flow_incr = edge_accessor.increase_flow(flow_incr);
				}

				node_accessor.set_deficit(node_accessor.get_deficit() + flow_incr);
				edge_accessor.set_deficit(edge_accessor.get_deficit() - flow_incr);
			} else {
				auto i = edge_accessor.get_tail();
				auto j = edge_accessor.get_head();

				auto i_accessor = NodeAccessor(N, i);
				auto j_accessor = NodeAccessor(N, j);

				auto flow_incr = i_accessor.get_depth() > j_accessor.get_depth() ? -i_accessor.get_deficit() : j_accessor.get_deficit();

				// only update edges with slack 0
				if (!((flow_incr < 0 && edge_accessor.get_backward_cost() + j_accessor.get_potential() - i_accessor.get_potential() == 0)
					|| (flow_incr > 0 && edge_accessor.get_forward_cost() + i_accessor.get_potential() - j_accessor.get_potential() == 0)))
					continue;

				flow_incr = edge_accessor.increase_flow(flow_incr);

				i_accessor.set_deficit(i_accessor.get_deficit() + flow_incr);
				j_accessor.set_deficit(j_accessor.get_deficit() - flow_incr);
			}
		}
	}
}


template<
	typename Network,
	typename EdgeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename TreeType
>
void update_x_by_tree_rc_in_original_graph(Network &N, const TreeType &tree, const EdgeIterator &ebegin, const EdgeIterator &eend) {
	for (auto edge_it = ebegin; edge_it != eend; ++edge_it) {
		auto edge_accessor = EdgeAccessor(N, edge_it);
		if (!edge_accessor.is_transformed())
			continue;

		if (edge_accessor.is_lower_in_tree() && !edge_accessor.is_upper_in_tree()) {
			auto v_accessor = NodeAccessor(N, edge_accessor.get_tail());
			auto w_accessor = NodeAccessor(N, edge_accessor.get_head());
			assert(edge_accessor.get_forward_cost() + v_accessor.get_potential() - edge_accessor.get_potential() == 0);
			auto flow_incr = edge_accessor.increase_flow_to_capacity();
			edge_accessor.increase_upper_flow(-flow_incr);
			v_accessor.set_deficit(v_accessor.get_deficit() + flow_incr);
			w_accessor.set_deficit(w_accessor.get_deficit() - flow_incr);
			assert(edge_accessor.get_deficit() == 0);
		}

		if (!edge_accessor.is_lower_in_tree() && edge_accessor.is_upper_in_tree()) {
			auto v_accessor = NodeAccessor(N, edge_accessor.get_tail());
			auto w_accessor = NodeAccessor(N, edge_accessor.get_head());
			assert(w_accessor.get_potential() - edge_accessor.get_potential() == 0);
			auto flow_incr = edge_accessor.increase_upper_flow_to_capacity();
			edge_accessor.increase_flow(-flow_incr);
			v_accessor.set_deficit(v_accessor.get_deficit() - flow_incr);
			w_accessor.set_deficit(w_accessor.get_deficit() + flow_incr);
			assert(edge_accessor.get_deficit() == 0);
		}
	}

	for (const auto &depth : boost::adaptors::reverse(tree)) {
		for (const auto &tree_node : depth) {
			auto v = tree_node.first;
			auto a = tree_node.second;

			auto edge_accessor = EdgeAccessor(N, a);

			assert(((edge_accessor.get_head() == v) != (edge_accessor.get_tail() == v)) && "the node must be either the edge's head or its tail");

			auto i = edge_accessor.get_tail();
			auto j = edge_accessor.get_head();

			auto i_accessor = NodeAccessor(N, i);
			auto j_accessor = NodeAccessor(N, j);

			auto flow_incr = i_accessor.get_depth() > j_accessor.get_depth() ? -i_accessor.get_deficit() : j_accessor.get_deficit();

			// only update edges with slack 0
			if (!((flow_incr < 0 && edge_accessor.get_backward_cost() + j_accessor.get_potential() - i_accessor.get_potential() == 0)
				|| (flow_incr > 0 && edge_accessor.get_forward_cost() + i_accessor.get_potential() - j_accessor.get_potential() == 0)))
				continue;

			flow_incr = edge_accessor.increase_flow(flow_incr);
			if (edge_accessor.is_transformed())
				edge_accessor.increase_upper_flow(-flow_incr);

			i_accessor.set_deficit(i_accessor.get_deficit() + flow_incr);
			j_accessor.set_deficit(j_accessor.get_deficit() - flow_incr);
		}
	}
}


template<
	typename Network,
	typename NodeIterator,
	typename EdgeIterator,
	typename NodeAccessor,
	typename EdgeAccessor,
	typename IncidentEdgesIterator
>
void successive_shortest_path_rc(Network &N, const NodeIterator &nbegin, const NodeIterator &nend, const EdgeIterator &ebegin, const EdgeIterator &eend) {

#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
#undef ALLOW_START_FROM_EDGES
#endif

	using PotentialType = decltype(NodeAccessor(N, nbegin).get_potential());
	using DeficitType = decltype(NodeAccessor(N, nbegin).get_deficit());

	initialize_internal_data_structures_rc(N);
	initialize_deficit_rc(N);

	auto num_nodes = get_number_of_nodes(N);
	auto num_transformed_edges = 0u;

	for (auto edge_it = ebegin; edge_it != eend; ++edge_it)
		if (EdgeAccessor(N, edge_it).is_transformed())
			++num_transformed_edges;

	auto step_count = 0u;

	do {
#ifndef NDEBUG
		auto sum_deficit = static_cast<DeficitType>(0);

		for (auto node_it = nbegin; node_it != nend; ++node_it) {
			auto node_accessor = NodeAccessor(N, node_it);
			sum_deficit += node_accessor.get_deficit();
		}

		for (auto edge_it = ebegin; edge_it != eend; ++edge_it) {
			auto edge_accessor = EdgeAccessor(N, edge_it);
			if (edge_accessor.is_transformed())
				sum_deficit += edge_accessor.get_deficit();
		}

		assert(sum_deficit == 0 && "the sum of all deficits must be 0");

		run_debug_checks(N);
#endif

		auto v_start = nend;
		auto v_start_deficit = static_cast<DeficitType>(0);

		for (auto node_it = nbegin; node_it != nend; ++node_it) {
			auto node_accessor = NodeAccessor(N, node_it);
			if (node_accessor.get_deficit() < v_start_deficit) {
				v_start = node_it;
				v_start_deficit = node_accessor.get_deficit();
			}
		}

#ifdef ALLOW_START_FROM_EDGES
		auto a_start = eend;
		auto a_start_deficit = static_cast<DeficitType>(0);

		for (auto edge_it = ebegin; edge_it != eend; ++edge_it) {
			auto edge_accessor = EdgeAccessor(N, edge_it);
			if (edge_accessor.is_transformed() && edge_accessor.get_deficit() < a_start_deficit) {
				a_start = edge_it;
				a_start_deficit = edge_accessor.get_deficit();
			}
		}
#endif

#ifdef ALLOW_START_FROM_EDGES
		if (v_start == nend && a_start == eend) {
#else
		if (v_start == nend) {
#endif
			if (step_count == 0u) {
				v_start = nbegin;
			} else {
				break;
			}
		}

		auto deficit_1 = static_cast<DeficitType>(0);

		for (auto node_it = nbegin; node_it != nend; ++node_it) {
			auto node_accessor = NodeAccessor(N, node_it);
			deficit_1 += std::abs(node_accessor.get_deficit());
		}

#ifdef ALLOW_START_FROM_EDGES
		for (auto edge_it = ebegin; edge_it != eend; ++edge_it) {
			auto edge_accessor = EdgeAccessor(N, edge_it);
			if (edge_accessor.is_transformed())
				deficit_1 += std::abs(edge_accessor.get_deficit());
		}
#endif

		std::cout << "starting step " << ++step_count << ", deficit_1 = " << deficit_1 << std::endl;

		using PriorityQueueElementType = std::pair<PotentialType, std::pair<EdgeIterator, bool>>;

		auto compare_first = [](const PriorityQueueElementType &lhs, const PriorityQueueElementType &rhs) {
			return lhs.first > rhs.first;
		};

		using PriorityQueueType = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>;

		auto s_in_n = PriorityQueueType(compare_first);
		auto s_out_n = PriorityQueueType(compare_first);
		auto s_in_a = PriorityQueueType(compare_first);
		auto s_out_a = PriorityQueueType(compare_first);

		for (auto node_it = nbegin; node_it != nend; ++node_it) {
			auto node_accessor = NodeAccessor(N, node_it);
			node_accessor.set_visited(false);
			node_accessor.set_depth(-1);
		}

		for (auto edge_it = ebegin; edge_it != eend; ++edge_it) {
			auto edge_accessor = EdgeAccessor(N, edge_it);
			if (edge_accessor.is_transformed()) {
				edge_accessor.set_visited(false);
#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
				edge_accessor.set_lower_in_tree(false);
				edge_accessor.set_upper_in_tree(false);
#else
				edge_accessor.set_depth(-1);
#endif
			}
		}

#ifdef ALLOW_START_FROM_EDGES
		bool use_v_start = v_start != nend && (a_start == eend || v_start_deficit < a_start_deficit);
		auto deficit_s = use_v_start ? v_start_deficit : a_start_deficit;

		if (use_v_start) {
			auto v_start_accessor = NodeAccessor(N, v_start);
			v_start_accessor.set_visited(true);
			v_start_accessor.set_depth(0);
			update_s_in_s_out_node_rc<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, PriorityQueueType>(N, v_start, s_in_n, s_out_n, s_in_a, s_out_a);
		} else {
			auto a_start_accessor = EdgeAccessor(N, a_start);
			a_start_accessor.set_visited(true);
			a_start_accessor.set_depth(0);
			update_s_in_s_out_edge_rc<Network, EdgeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, PriorityQueueType>(N, a_start, s_in_a, s_out_a);
		}
#else
		auto deficit_s = v_start_deficit;
		auto v_start_accessor = NodeAccessor(N, v_start);
		v_start_accessor.set_visited(true);
		v_start_accessor.set_depth(0);
		update_s_in_s_out_node_rc<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, PriorityQueueType>(N, v_start, s_in_n, s_out_n, s_in_a, s_out_a);
#endif

		auto num_visited = 1u;

#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
		auto tree = std::vector<std::vector<std::pair<NodeIterator, EdgeIterator>>>();
#else
		auto tree = std::vector<std::vector<std::pair<NodeIterator, std::pair<EdgeIterator, bool>>>>();
#endif

		while (num_visited < num_nodes + num_transformed_edges && !(s_in_n.empty() && s_in_a.empty() && s_out_n.empty() && s_out_a.empty())) {

			bool visited;
			DeficitType new_deficit;

			auto pop_top = [](PriorityQueueType &queue) -> PriorityQueueElementType {
				auto top = queue.top();
				queue.pop();
				return top;
			};

			if (deficit_s < 0 || (s_in_n.empty() && s_in_a.empty())) {
				assert(!(s_out_n.empty() && s_out_a.empty()));
				if (!s_out_n.empty() && (s_out_a.empty() || s_out_n.top().first < s_out_a.top().first)) {
					std::tie(visited, new_deficit) = explore_arc_rc<Network, EdgeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, PriorityQueueType, decltype(tree)>(N, pop_top(s_out_n), s_in_a, s_out_a, tree, true);
				} else {
					std::tie(visited, new_deficit) = explore_node_rc<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, PriorityQueueType, decltype(tree)>(N, pop_top(s_out_a), s_in_n, s_out_n, s_in_a, s_out_a, tree, true);
				}
			} else {
				assert(!(s_in_n.empty() && s_in_a.empty()));
				if (!s_in_n.empty() && (s_in_a.empty() || s_in_n.top().first < s_in_a.top().first)) {
					std::tie(visited, new_deficit) = explore_arc_rc<Network, EdgeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, PriorityQueueType, decltype(tree)>(N, pop_top(s_in_n), s_in_a, s_out_a, tree, false);
				} else {
					std::tie(visited, new_deficit) = explore_node_rc<Network, NodeIterator, NodeAccessor, EdgeAccessor, IncidentEdgesIterator, DeficitType, PriorityQueueType, decltype(tree)>(N, pop_top(s_in_a), s_in_n, s_out_n, s_in_a, s_out_a, tree, false);
				}
			}

			if (visited) {
				deficit_s += new_deficit;
				++num_visited;
			}

		}

#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
		update_x_by_tree_rc_in_original_graph<Network, EdgeIterator, NodeAccessor, EdgeAccessor, decltype(tree)>(N, tree, ebegin, eend);
#else
		update_x_by_tree_rc<Network, NodeAccessor, EdgeAccessor, decltype(tree)>(N, tree);
#endif

	} while (true);

}
