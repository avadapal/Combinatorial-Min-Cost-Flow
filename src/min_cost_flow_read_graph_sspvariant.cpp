#include "graph.h"
#include "random.h"
#include "min_cost_flow.h"

#include <iostream>
#include <queue>
#include <stack>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <boost/rational.hpp>
#include <tuple>
#include <chrono>

using namespace std;

int main(){

	cout << "Enter filepath to input graph.. " << endl;
	string filename;
	cin >> filename;

	using RationalType = long long int;
	using IntegerType = long long int;

	Graph<IntegerType, RationalType> G(filename);
#ifdef RESP_CAP
	Network<decltype(G), IntegerType, RationalType, SSPVariantRCNodeData<IntegerType, RationalType>, SSPVariantRCArcData<IntegerType, RationalType>> N(G, filename);
#else
	Network<decltype(G), IntegerType, RationalType, SSPVariantNodeData<IntegerType, RationalType>, BasicArcData<RationalType>> N(G, filename);
#endif

  // does not assign resistances, name is stupid
  //assign_costs_resistances_capacities(N);
  //assign_demands(N);

	cout << "Number of Edges: " << N.G.no_of_edges << endl;
	cout << "Number of Nodes: " << N.G.no_of_vertices << endl;

//  unsigned int q;

//  G0.print_lp();
//  write_the_graph_in_dimacs_format( N );  

	auto start = std::chrono::steady_clock::now();

	//set_initial_x_s_y(N);

#ifdef RESP_CAP
	successive_shortest_path_resp_cap(N);
#else
	successive_shortest_path(N);
#endif

	auto elapsed_time = std::chrono::duration<double>(std::chrono::steady_clock::now() - start);
	std::cout << "total time: " << elapsed_time.count() << "s" << std::endl;

	
	long long total_cost = 0;
	for(unsigned int i = 1; i <= N.G.no_of_edges; i++)
	{
	  total_cost += N.arcdata[i].xlower * N.arcdata[i].cost;
	  if(N.arcdata[i].xlower != 0)
	  {
	    cout << "x[" << i << "] = " << N.arcdata[i].xlower << endl;
	  }
	}
	cout << "Total Cost = " << total_cost << endl;
#ifndef NDEBUG
#ifdef RESP_CAP
	// check flow conservation
	for (auto v = 1u; v < N.G.no_of_vertices; ++v) {
		auto sum_v = static_cast<decltype(N.arcdata.front().xlower)>(0);
		for (auto edge : N.G.incident_edges[v]) {
			if (edge > 0) {
				sum_v -= N.arcdata[edge].xlower;
			} else {
				sum_v += N.arcdata[-edge].xlower;
			}
		}
		assert(sum_v == N.nodedata[v].demand && "flow conservation check failed");
	}

	// check capacity constraints
	for (auto a = 1u; a < N.G.no_of_edges; ++a) {
		if (N.arcdata[a].transformed) {
			assert(N.arcdata[a].xlower + N.arcdata[a].xupper == N.arcdata[a].capacity && "primal feasibility check failed");
		} else {
			assert(N.arcdata[a].xlower <= N.arcdata[a].capacity && "capacity constraints check failed");
		}
	}

	// check complementary slackness
	for (auto a = 1u; a < N.G.no_of_edges; ++a) {
		const auto &arcdata = N.arcdata[a];
		if (arcdata.transformed) {
			const auto s_lower = arcdata.cost + N.nodedata[N.G.tails[a]].potential - arcdata.potential;
			const auto s_upper = N.nodedata[N.G.heads[a]].potential - arcdata.potential;
			assert(s_lower * arcdata.xlower == 0 && "complementary slackness check failed");
			assert(s_upper * arcdata.xupper == 0 && "complementary slackness check failed");
		} else {
			const auto &i = N.G.tails[a];
			const auto &j = N.G.heads[a];
			assert(arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential <= 0 || arcdata.xlower == 0);
			assert(arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential >= 0 || arcdata.xlower == arcdata.capacity);
			assert(!(arcdata.xlower > 0 && arcdata.xlower < arcdata.capacity) || arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential == 0);
		}
	}
#else
	// check flow conservation
	for (auto v = 1u; v < N.G.no_of_vertices; ++v) {
		auto sum_v = static_cast<decltype(N.arcdata.front().xlower)>(0);
		for (auto edge : N.G.incident_edges[v]) {
			if (edge > 0) {
				sum_v -= N.arcdata[edge].xlower;
			} else {
				sum_v += N.arcdata[-edge].xlower;
			}
		}
		assert(sum_v == N.nodedata[v].demand && "flow conservation check failed");
	}

	// check capacity constraints
	for (auto a = 1u; a < N.G.no_of_edges; ++a)
		assert(N.arcdata[a].xlower <= N.arcdata[a].capacity && "capacity constraints check failed");

	// check complementary slackness
	for (auto a = 1u; a < N.G.no_of_edges; ++a) {
		const auto &arcdata = N.arcdata[a];
		const auto &i = N.G.tails[a];
		const auto &j = N.G.heads[a];
		assert(arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential <= 0 || arcdata.xlower == 0);
		assert(arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential >= 0 || arcdata.xlower == arcdata.capacity);
		assert(!(arcdata.xlower > 0 && arcdata.xlower < arcdata.capacity) || arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential == 0);
	}
#endif
#endif

	/*
  cout << "Primal and Dual Solution: " << endl;
  for(unsigned int a = 1; a <= N.G.no_of_edges; a++){
	  cout << "x[" << a << "] = " << N.arcdata[a].xlower << endl;
  }
  for(unsigned int v = 1; v<= N.G.no_of_vertices; v++){
    cout << "y[" << v << "] = " << N.nodedata[v].potential << endl;
  }
  cout << endl;
  
  cout << "Graph: " << endl;
  cout << "Vertex: Incident Edges" << endl;
  for(unsigned int v = 1; v <= N.G.no_of_vertices; v++){
    cout << v << " :  ";
    for(auto a: N.G.incident_edges[v]){
      cout << a << ", ";
    }
    cout << endl;
  }
  */
  return 0;
}

