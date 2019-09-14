#include "graph.h"
#include "random.h"
#include "min_cost_flow_sspvariant.h"
#include "min_cost_flow_sspvariant_default.h"
//#include "min_cost_flow_sspvariant_apex_grid.h"

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
	/*
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
	*/

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

	ofstream dimacs_out(filename + ".min");

	N.write_dimacs(dimacs_out);

	dimacs_out.close();
	exit(0);

	//cout << "Number of Edges: " << N.G.no_of_edges << endl;
	//cout << "Number of Nodes: " << N.G.no_of_vertices << endl;

//  unsigned int q;

  return 0;
}

