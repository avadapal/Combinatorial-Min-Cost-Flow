#pragma once

#include "min_cost_flow_sspvariant.h"

#define SSP_VARIANT
#include "gridflow.h"


// node as i and j coordinates; the apex is -1, -1
// cannont use the Node struct from apex_grid.h since the operator== only compares keys
using ApexGridNodeType = std::pair<int, int>;

using ApexGridType = ApexGrid<8>;
using ApexGridNodeDataType = NodeData<8>;


class ApexGridNodeIterator;
class ApexGridEdgeIterator;
class ApexGridNodeAccessor;
class ApexGridEdgeAccessor;



class ApexGridNodeIterator {

public:
	ApexGridNodeIterator() = delete;
	ApexGridNodeIterator(const ApexGridNodeIterator &) = default;

	ApexGridNodeIterator(ApexGridType &apex_grid, int i, int j) : width(apex_grid.width), current_node(i, j) {}

	auto operator++() -> ApexGridNodeIterator & {
		if (current_node.first == -1 && current_node.second == -1) {
			current_node.first = 0;
			current_node.second = 0;
			return *this;
		}

		if (current_node.second == width - 1) {
			++current_node.first;
			current_node.second = 0;
			return *this;
		}

		++current_node.second;
		return *this;
	}

	auto operator*() const -> std::pair<int, int> {
		return current_node;
	}

	auto operator==(const ApexGridNodeIterator &other) const -> bool {
		return current_node == other.current_node;
	}

	auto operator!=(const ApexGridNodeIterator &other) const -> bool {
		return !operator==(other);
	}

private:
	int width;
	ApexGridNodeType current_node;

};


class ApexGridEdgeIterator {

public:
	ApexGridEdgeIterator() = delete;
	ApexGridEdgeIterator(const ApexGridEdgeIterator &) = default;

	ApexGridEdgeIterator(const ApexGridType &apex_grid, int i, int j, int edge_id = 0) : width(apex_grid.width), current_node(i, j), edge_id(edge_id) {}

	auto operator++() -> ApexGridEdgeIterator & {
		if (edge_id == 2) {
			// reset edge_id and go to next node
			edge_id = 0;

			if (current_node.second == width - 1) {
				++current_node.first;
				current_node.second = 0;
				return *this;
			}

			++current_node.second;
			return *this;
		}

		++edge_id;
		return *this;
	}

	auto operator*() const -> std::pair<ApexGridNodeType, int> {
		return {current_node, edge_id};
	}

	auto operator==(const ApexGridEdgeIterator &other) const -> bool {
		return current_node == other.current_node && edge_id == other.edge_id;
	}

	auto operator!=(const ApexGridEdgeIterator &other) const -> bool {
		return !operator==(other);
	}

	friend class ApexGridEdgeAccessor;

private:
	int width;
	ApexGridNodeType current_node;

	/* edge_id ranges from 0 to 2 and each value can represent multiple edges:
	 * 0 are the edges to the right
	 * 1 are the edges downwards
	 * 2 are the edges connecting this node to the apex */
	int edge_id;

};


class ApexGridNodeAccessor {

public:
	ApexGridNodeAccessor(ApexGridType &apex_grid, const ApexGridNodeIterator &node_it) : nodedata((*node_it).first == -1 ? apex_grid.apex : apex_grid.grid[(*node_it).first][(*node_it).second]) {}

	auto get_potential() const -> int {
		return nodedata.potential;
	}

	void set_potential(int potential) {
		nodedata.potential = potential;
	}

	auto get_deficit() const -> int {
		return nodedata.deficiency;
	}

	void set_deficit(int deficit) {
		nodedata.deficiency = deficit;
	}

#ifdef RESTORE_BALANCED_NODES
	auto get_deficit_delta() const -> int {
		assert(false && "not implemented");
		return -1;
	}

	void set_deficit_delta(int) {
		assert(false && "not implemented");
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

private:
	ApexGridNodeDataType &nodedata;

};


class ApexGridEdgeAccessor {

public:
	ApexGridEdgeAccessor(ApexGridType &apex_grid, const ApexGridEdgeIterator &edge_it) :
		apex_grid(apex_grid),
		edge_id(edge_it.edge_id),
		tail(edge_it.current_node),
		head() {
		switch (edge_id) {
		case 0:
			head = {tail.first, tail.second + 1};
			break;
		case 1:
			head = {tail.first + 1, tail.second};
			break;
		case 2:
			head = {-1, -1};
			break;
		default:
			assert(false);
			break;
		}
	}

	auto get_head() const -> ApexGridNodeIterator {
		return ApexGridNodeIterator(apex_grid, head.first, head.second);
	}

	auto get_tail() const -> ApexGridNodeIterator {
		return ApexGridNodeIterator(apex_grid, tail.first, tail.second);
	}

	auto increase_flow(int amount) -> int {
		if (edge_id == 2) {
			auto &flow = apex_grid.grid[tail.first][tail.second].flowtoapex;
			if (amount > 0) {
				if (flow < 0) {
					auto delta = std::min(amount, -flow);
					flow += delta;
					return delta;
				} else {
					flow += amount;
					return amount;
				}
			} else {
				if (flow > 0) {
					auto delta = -std::min(-amount, flow);
					flow += delta;
					return delta;
				} else {
					flow += amount;
					return amount;
				}
			}
		}
		assert(edge_id == 0 || edge_id == 1);
		auto &flow = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right.flow : apex_grid.grid[tail.first][tail.second].down.flow;
		if (amount > 0 && flow < apex_grid.n) {
			++flow;
			return 1;
		}
		if (amount < 0 && flow > -apex_grid.n) {
			--flow;
			return -1;
		}
		return 0;
	}

	auto can_increase_flow() const -> bool {
#ifdef RESP_CAP
		return true;
#else
		switch (edge_id) {
		case 0:
			return apex_grid.grid[tail.first][tail.second].right.flow < apex_grid.n;
		case 1:
			return apex_grid.grid[tail.first][tail.second].down.flow < apex_grid.n;
		case 2:
			return true;
		default:
			assert(false);
			return false;
		}
#endif
	}

	auto can_decrease_flow() const -> bool {
		switch (edge_id) {
		case 0:
			return apex_grid.grid[tail.first][tail.second].right.flow > -apex_grid.n;
		case 1:
			return apex_grid.grid[tail.first][tail.second].down.flow > -apex_grid.n;
		case 2:
			return true;
		default:
			assert(false);
			return false;
		}
	}

#ifdef NEW_OPTIMIZATION
	auto get_capacity() const -> int {
		assert(false && "not implemented");
		return -1;
	}

	auto get_residual_capacity() const -> int {
		assert(false && "not implemented");
		return -1;
	}

	auto get_residual_backward_capacity() const -> int {
		assert(false && "not implemented");
		return -1;
	}
#endif

	auto get_forward_cost() const -> int {
		int flow;
		switch (edge_id) {
		case 0:
			flow = apex_grid.grid[tail.first][tail.second].right.flow;
			break;
		case 1:
			flow = apex_grid.grid[tail.first][tail.second].down.flow;
			break;
		case 2:
			flow = apex_grid.grid[tail.first][tail.second].flowtoapex;
			return flow >= 0 ? apex_grid.grid[tail.first][tail.second].apexcost : 0;
		default:
			assert(false);
			return -1;
		}
		auto r = std::min(apex_grid.n - 1, apex_grid.ReturnR(flow));
		//auto r = apex_grid.ReturnR(flow);
		//assert(r < apex_grid.n);
		return apex_grid.grid[tail.first][tail.second].right.difference[r];
	}

	auto get_backward_cost() const -> int {
		int flow;
		switch (edge_id) {
		case 0:
			flow = apex_grid.grid[tail.first][tail.second].right.flow;
			break;
		case 1:
			flow = apex_grid.grid[tail.first][tail.second].down.flow;
			break;
		case 2:
			flow = apex_grid.grid[tail.first][tail.second].flowtoapex;
			return flow <= 0 ? 0 : -apex_grid.grid[tail.first][tail.second].apexcost;
		default:
			assert(false);
			return -1;
		}
		auto l = std::max(0, apex_grid.ReturnL(flow));
		//auto l = apex_grid.ReturnL(flow);
		//assert(l >= 0);
		return -apex_grid.grid[tail.first][tail.second].right.difference[l];
	}

#ifdef RESTORE_BALANCED_NODES
	auto get_flow_delta() const -> int {
		assert(false && "not implemented");
		return -1;
	}

	void set_flow_delta(int) {
		assert(false && "not implemented");
	}
#endif

#ifdef RESP_CAP
	auto is_transformed() const -> bool {
		return tail.first != -1 && edge_id < 2;
	}

	auto can_decrease_upper_flow() const -> bool {
		assert(is_transformed());
		const auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		return arcdata.xupper > 0;
	}

	auto increase_upper_flow(int amount) -> int {
		assert(is_transformed());
		auto &flow = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right.xupper : apex_grid.grid[tail.first][tail.second].down.xupper;
		if (amount < 0 && flow > 0) {
			--flow;
			return -1;
		}
		if (amount > 0) {
			++flow;
			return 1;
		}
		return 0;
	}

	auto get_potential() const -> int {
		assert(is_transformed());
		const auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		return arcdata.potential;
	}

	void set_potential(int potential) {
		assert(is_transformed());
		auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		arcdata.potential = potential;
	}

	auto get_deficit() const -> int {
		assert(is_transformed());
		const auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		return arcdata.deficit;
	}

	void set_deficit(int deficit) {
		assert(is_transformed());
		auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		arcdata.deficit = deficit;
	}

	auto get_depth() const -> int {
		assert(is_transformed());
		const auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		return arcdata.depth;
	}

	void set_depth(int depth) {
		assert(is_transformed());
		auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		arcdata.depth = depth;
	}

	auto is_visited() const -> bool {
		assert(is_transformed());
		const auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		return arcdata.visited;
	}

	void set_visited(bool visited) {
		assert(is_transformed());
		auto &arcdata = edge_id == 0 ? apex_grid.grid[tail.first][tail.second].right : apex_grid.grid[tail.first][tail.second].down;
		arcdata.visited = visited;
	}

#ifdef BUILD_TREE_IN_ORIGINAL_GRAPH
	auto is_lower_in_tree() const -> bool {
		// TODO
		assert(false);
		return false;
	}

	void set_lower_in_tree(bool lower_in_tree) {
		// TODO
		assert(false);
	}

	auto is_upper_in_tree() const -> bool {
		// TODO
		assert(false);
		return false;
	}

	void set_upper_in_tree(bool upper_in_tree) {
		// TODO
		assert(false);
	}

	auto increase_flow_to_capacity() -> int {
		// TODO
		assert(false);
		return 0;
	}

	auto increase_upper_flow_to_capacity() -> int {
		// TODO
		assert(false);
		return 0;
	}
#endif
#endif

private:
	ApexGridType &apex_grid;

	/* edge_id ranges from 0 to 2 and each value can represent multiple edges:
	 * 0 are the edges to the right
	 * 1 are the edges downwards
	 * 2 are the edges connecting this node to the apex */
	int edge_id;
	ApexGridNodeType tail;
	ApexGridNodeType head;

};


/* iteration order: 
 * 1. right
 * 2. down
 * 3. left
 * 4. up
 * 5. to apex
 * 6. from apex */
class ApexGridIncidentEdgesIterator {

public:
	ApexGridIncidentEdgesIterator() = delete;
	ApexGridIncidentEdgesIterator(const ApexGridIncidentEdgesIterator &) = default;

	ApexGridIncidentEdgesIterator(const ApexGridType &apex_grid, ApexGridNodeIterator node_it, int edge_id = 0) : apex_grid(apex_grid), v(*node_it), edge_id(edge_id) {
		advance_if_at_border();
	}

	auto operator++() -> ApexGridIncidentEdgesIterator & {
		++edge_id;
		advance_if_at_border();
		return *this;
	};

	auto operator*() const -> ApexGridEdgeIterator {
		if (v.first != -1) {
			switch (edge_id) {
			case 0:
				return ApexGridEdgeIterator(apex_grid, v.first, v.second, 0);
			case 1:
				return ApexGridEdgeIterator(apex_grid, v.first, v.second, 1);
			case 2:
				return ApexGridEdgeIterator(apex_grid, v.first, v.second - 1, 0);
			case 3:
				return ApexGridEdgeIterator(apex_grid, v.first - 1, v.second, 1);
			case 4:
				return ApexGridEdgeIterator(apex_grid, v.first, v.second, 2);
			default:
				assert(false);
				return ApexGridEdgeIterator(apex_grid, -1, -1, 0);
			}
		} else {
			return ApexGridEdgeIterator(apex_grid, edge_id / apex_grid.width, edge_id % apex_grid.width, 2);
		}
	};

	auto operator==(const ApexGridIncidentEdgesIterator &other) const -> bool {
		return v == other.v && edge_id == other.edge_id;
	}

	auto operator!=(const ApexGridIncidentEdgesIterator &other) const -> bool {
		return !operator==(other);
	}

private:
	void advance_if_at_border() {
		if (v.first != -1) {
			// right border
			if (edge_id == 0 && v.second == apex_grid.width - 1)
				++edge_id;

			// lower border
			if (edge_id == 1 && v.first == apex_grid.height - 1)
				++edge_id;

			// left border
			if (edge_id == 2 && v.second == 0)
				++edge_id;

			// upper border
			if (edge_id == 3 && v.first == 0)
				++edge_id;
		}
	}

	const ApexGridType &apex_grid;
	const ApexGridNodeType v;

	/* for normal nodes, edge_id ranges from 0 to 4 and each value can represent multiple edges:
	 *   0 are the edges to the right
	 *   1 are the edges downwards
	 *   2 are the edges to the left
	 *   3 are the edges upwards
	 *   4 are the edges connecting this node to the apex
	 * for the apex node, edge_id ranges from 0 to width * height - 1
	 */
	int edge_id;

};


template<>
void initialize_internal_data_structures(ApexGridType &apex_grid) {
	apex_grid.initialize();
};


template<>
void initialize_deficit(ApexGridType &apex_grid) {
	std::vector<Node> Sources;
	assert(!Sources.size());
	apex_grid.FlowInitialization(0, apex_grid.height, 0, apex_grid.width);
	apex_grid.DeficiencyInitialization(Sources, 0, apex_grid.height, 0, apex_grid.width);

	apex_grid.for_each_nodedata([](decltype(apex_grid.apex) &nodedata) {nodedata.deficiency = -nodedata.deficiency; },
		0, apex_grid.height, 0, apex_grid.width);
};


template<>
auto get_number_of_nodes(const ApexGridType &apex_grid) -> unsigned int {
	return apex_grid.width * apex_grid.height + 1;
};


template<>
auto begin_incident_edges(const ApexGridType &apex_grid, const ApexGridNodeIterator &node) -> ApexGridIncidentEdgesIterator {
	return ApexGridIncidentEdgesIterator(apex_grid, node);
};


template<>
auto end_incident_edges(const ApexGridType &apex_grid, const ApexGridNodeIterator &node) -> ApexGridIncidentEdgesIterator {
	if ((*node).first == -1) {
		return ApexGridIncidentEdgesIterator(apex_grid, node, apex_grid.width * apex_grid.height);
	} else {
		return ApexGridIncidentEdgesIterator(apex_grid, node, 5);
	}
};


template<>
void run_debug_checks(const ApexGridType &) {
	// TODO
}


template<>
auto compute_dual_objective_value(const ApexGridType &) -> int {
	assert(false && "not implemented");
	// TODO: implement
	return -1;
}


template<>
auto compute_demand_s(const ApexGridType &N) -> int {
	assert(false && "not implemented");
	// TODO: implement
	return -1;
}


#ifdef RESP_CAP

template<>
void initialize_internal_data_structures_rc(ApexGridType &apex_grid) {
	initialize_internal_data_structures(apex_grid);
};

template<>
void initialize_deficit_rc(ApexGridType &apex_grid) {
	initialize_deficit(apex_grid);

	// set initial xupper
	for (auto i = 0; i < apex_grid.height; ++i) {
		for (auto j = 0; j < apex_grid.width; ++j) {
			if (i < apex_grid.height)
				apex_grid.grid[i][j].down.xupper = 2 * apex_grid.n;
			if (j < apex_grid.width)
				apex_grid.grid[i][j].right.xupper = 2 * apex_grid.n;
		}
	}
};

#endif



