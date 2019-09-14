#include "graph.h"
#include "random.h"
#include "min_cost_flow_sspvariant.h"
#include "min_cost_flow_sspvariant_default.h"
#include "min_cost_flow_sspvariant_apex_grid.h"

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

	// check if file exists
	if (!std::ifstream(filename)) {
		std::cerr << "File not found: " << filename << std::endl;
		exit(1);
	}

	if (boost::ends_with(filename, ".pgm") || boost::ends_with(filename, ".ppm")) {
#ifndef NEW_OPTIMIZATION
		ifstream file(filename);
		PPMImage image(file);
		file.close();

		InstanceConversion<8> instance(image);

		auto start = std::chrono::steady_clock::now();

#ifdef RESP_CAP
		successive_shortest_path_rc<decltype(instance.apexgrid), ApexGridNodeIterator, ApexGridEdgeIterator, ApexGridNodeAccessor, ApexGridEdgeAccessor, ApexGridIncidentEdgesIterator>(
			instance.apexgrid,
			ApexGridNodeIterator(instance.apexgrid, -1, -1),
			ApexGridNodeIterator(instance.apexgrid, instance.apexgrid.height, 0),
			ApexGridEdgeIterator(instance.apexgrid, -1, -1),
			ApexGridEdgeIterator(instance.apexgrid, instance.apexgrid.height, 0));
#else
		successive_shortest_path<decltype(instance.apexgrid), ApexGridNodeIterator, ApexGridEdgeIterator, ApexGridNodeAccessor, ApexGridEdgeAccessor, ApexGridIncidentEdgesIterator>(
			instance.apexgrid,
			ApexGridNodeIterator(instance.apexgrid, -1, -1),
			ApexGridNodeIterator(instance.apexgrid, instance.apexgrid.height, 0));
#endif

		auto elapsed_time = std::chrono::duration<double>(std::chrono::steady_clock::now() - start);
		std::cout << "total time: " << elapsed_time.count() << "s" << std::endl;

		assert(instance.apexgrid.check(0, instance.apexgrid.height, 0, instance.apexgrid.width));

#endif
		
		return 0;
	}


	using RationalType = long long int;
	using IntegerType = long long int;

	Graph<IntegerType, RationalType> G(filename);
#ifdef RESP_CAP
	Network<decltype(G), IntegerType, RationalType, SSPVariantNodeData<IntegerType, RationalType>, SSPVariantRCArcData<IntegerType, RationalType>> N(G, filename);
#else
#ifdef RESTORE_BALANCED_NODES
	Network<decltype(G), IntegerType, RationalType, SSPVariantNodeData<IntegerType, RationalType>, SSPVariantRBNArcData<RationalType>> N(G, filename);
#else
	Network<decltype(G), IntegerType, RationalType, SSPVariantNodeData<IntegerType, RationalType>, BasicArcData<RationalType>> N(G, filename);
#endif
#endif

  // does not assign resistances, name is stupid
  //assign_costs_resistances_capacities(N);
  //assign_demands(N);

	cout << "Number of Edges: " << N.G.no_of_edges << endl;
	cout << "Number of Nodes: " << N.G.no_of_vertices << endl;

//	unsigned int runs = 0;
//	cout << "How many runs? "; cin >> runs;
	
//  unsigned int q;

//  G0.print_lp();
//  write_the_graph_in_dimacs_format( N );  

	auto start_chrono = std::chrono::steady_clock::now();
	auto start_clock = clock();

//	for( unsigned int r = 0; r < runs; ++r ) {
	  
#ifdef RESP_CAP
	successive_shortest_path_rc<decltype(N), DefaultNodeIterator, DefaultEdgeIterator, DefaultNodeAccessor, DefaultEdgeAccessor, DefaultIncidentEdgesIterator>(N, DefaultNodeIterator(1), DefaultNodeIterator(N.G.no_of_vertices + 1), DefaultEdgeIterator(1), DefaultEdgeIterator(N.G.no_of_edges + 1));
#else
#ifndef NEW_OPTIMIZATION
	successive_shortest_path<decltype(N), DefaultNodeIterator, DefaultEdgeIterator, DefaultNodeAccessor, DefaultEdgeAccessor, DefaultIncidentEdgesIterator>(N, DefaultNodeIterator(1), DefaultNodeIterator(N.G.no_of_vertices + 1));
#else
	successive_shortest_path<decltype(N), DefaultNodeIterator, DefaultEdgeIterator, DefaultNodeAccessor, DefaultEdgeAccessor, DefaultIncidentEdgesIterator>(N, DefaultNodeIterator(1), DefaultNodeIterator(N.G.no_of_vertices + 1), DefaultEdgeIterator(1), DefaultEdgeIterator(N.G.no_of_edges + 1));
#endif
#endif

//	}

	const auto time_chrono = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_chrono).count();
	const auto time_clock = ((float) (clock() - start_clock)) / CLOCKS_PER_SEC;

	std::cout << "total time with chrono: " << time_chrono << "s" << std::endl;
	std::cout << "total time with clock: " << time_clock << "s" << std::endl;

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
		if (!(sum_v == N.nodedata[v].demand)) {
			std::cerr << "flow conservation check failed" << std::endl;
			exit(1);
		}
	}

	// check complementary slackness et al
	for (auto a = 1u; a < N.G.no_of_edges; ++a) {
		const auto &arcdata = N.arcdata[a];
		const auto &i = N.G.tails[a];
		const auto &j = N.G.heads[a];

		// check primal feasibility
		if (!(arcdata.xlower <= arcdata.capacity)) {
			std::cerr << "capacity constraints check failed" << std::endl;
			exit(1);
		}

		// check dual feasibility
		if (!(arcdata.capacity != INFINITE_CAPACITY || arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential >= 0)) {
			std::cerr << "dual feasibility check failed" << std::endl;
			exit(1);
		}

		// check complementary slackness
		if (!(arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential <= 0 || arcdata.xlower == 0)) {
			std::cerr << "complementary slackness check 1 failed" << std::endl;
			exit(1);
		}
		if (!(arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential >= 0 || arcdata.xlower == arcdata.capacity)) {
			std::cerr << "complementary slackness check 2 failed" << std::endl;
			exit(1);
		}
		if (!(!(arcdata.xlower > 0 && arcdata.xlower < arcdata.capacity) || arcdata.cost + N.nodedata[i].potential - N.nodedata[j].potential == 0)) {
			std::cerr << "complementary slackness check 3 failed" << std::endl;
			exit(1);
		}
	}

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

