#include <iostream>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <signal.h>

#include <queue>
#include <stack>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <boost/rational.hpp>
#include <tuple>
#include <chrono>
#include <dirent.h>


#include "graph.h"
#include "random.h"
#include "min_cost_flow_sspvariant.h"
#include "min_cost_flow_sspvariant_default.h"
#include "min_cost_flow_sspvariant_apex_grid.h"

#ifndef OUR_IMPL

#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/cycle_canceling.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/capacity_scaling.h>
#include <lemon/dimacs.h>

using namespace lemon;

#endif

using namespace std;

/*
int getdir (string dir, vector<string> &files){
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

bool starts_with(const string& s1, const string& s2) {
    return s2.size() <= s1.size() && s1.compare(0, s2.size(), s2) == 0;
}

bool ends_with(const string& full_string, const string& ending) {
    if (full_string.length() >= ending.length()) {
        return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

*/


template <class Network>
struct Worker {

	Network &N;

	Worker(Network &N) : N(N) {}


  void operator() () {
#ifdef OUR_IMPL
		auto start_chrono = std::chrono::steady_clock::now();
		auto start_clock = clock();
#ifdef RESP_CAP
		successive_shortest_path_rc<typename std::remove_reference<decltype(N)>::type, DefaultNodeIterator, DefaultEdgeIterator, DefaultNodeAccessor, DefaultEdgeAccessor, DefaultIncidentEdgesIterator>(N, DefaultNodeIterator(1), DefaultNodeIterator(N.G.no_of_vertices + 1), DefaultEdgeIterator(1), DefaultEdgeIterator(N.G.no_of_edges + 1));
#else
		// auto x = begin_incident_edges<decltype(N), DefaultNodeIterator, DefaultIncidentEdgesIterator>(N,DefaultNodeIterator(1));
		successive_shortest_path<typename std::remove_reference<decltype(N)>::type, DefaultNodeIterator, DefaultEdgeIterator, DefaultNodeAccessor, DefaultEdgeAccessor, DefaultIncidentEdgesIterator>(N, DefaultNodeIterator(1), DefaultNodeIterator(N.G.no_of_vertices + 1));
#endif
    const auto time_chrono = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_chrono).count();
    const auto time_clock = ((float) (clock() - start_clock)) / CLOCKS_PER_SEC;

    std::cout << "total time with chrono: " << time_chrono << "s" << std::endl;
    std::cout << "total time with clock: " << time_clock << "s" << std::endl;
#else

#ifdef COST_SCALING
	  auto start = chrono::steady_clock::now();
	  N.run();
	  auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
	  cout << elapsed_time.count() << " s" << endl;
	  cout << endl;
#endif
#ifdef NETWORK_SIMPLEX
	  auto start = chrono::steady_clock::now();
	  N.run();
	  auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
	  cout << elapsed_time.count() << " s" << endl;
	  cout << endl;
#endif
#ifdef SUCC_SHORTEST_PATH
	  auto start = chrono::steady_clock::now();
	  N.run(false);
	  auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
	  cout << elapsed_time.count() << " s" << endl;
	  cout << endl;
#endif
#ifdef CAP_SCALING
	  auto start = chrono::steady_clock::now();
	  N.run();
	  auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
	  cout << elapsed_time.count() << " s" << endl;
	  cout << endl;
#endif
#ifdef CYCLE_CANCELING
	  auto start = chrono::steady_clock::now();
	  N.run();
	  auto elapsed_time = chrono::duration<double>(chrono::steady_clock::now() - start);
	  cout << elapsed_time.count() << " s" << endl;
	  cout << endl;
#endif

#endif

    //std::cerr << "Work done. Bye!" << std::endl;
    exit (0);
  }

  void initialize () {}
};

struct Kill {
  template< typename T>
  void operator() (T t) {
    //std::cerr << "Kill" << std::endl;
    exit (1);
  }
};


template <class Network>
class timed_job
{
private:
  boost::asio::io_service io_service_;
  boost::asio::deadline_timer timer_;
  Worker<Network> worker;
  Kill kill;
  boost::thread t;

public:
  timed_job( Network &N, int timeout ) :
    timer_( io_service_, boost::posix_time::seconds( timeout ) ),
    worker(N)// Deadline timer
  {
  }

  void start()
  {
    worker.initialize();
    // timer_.async_wait (kill);
    timer_.async_wait(boost::bind(&timed_job::stop, this));
  
    // Post your work
    io_service_.post(boost::bind(&timed_job::do_work, this));

    //std::cout << "Not run yet." << std::endl;
    io_service_.run();
    //std::cout << "stopped." << std::endl;
  }

private:
  void stop()
  {
    std::cout << "Timeout!" << std::endl;
    //io_service_.stop();
    exit(1);
  }

  void do_work ()
  {

   // std::cerr << "constructing..." << std::endl;

    t = boost::thread (worker);
    //std::cerr << "constructed..." << std::endl;
    t.detach();
    //std::cerr << "detached..." << std::endl;

    // // Keep posting the work.
    // io_service_.post
    //   (
    //    boost::bind
    //    (
    //     &timed_job::do_work, this
    //     )
    //    );
  }
};

int main(int argc, char *argv[]) {

	if (argc != 3) {
		std::cout << "usage: " << argv[0] << " <filename> <timeout>" << std::endl;
		exit(1);
	}

	const auto filename = argv[1];
	const auto timeout = atoi(argv[2]);

	if (!ifstream(filename)) {
		std::cerr << "File not found: " << filename << std::endl;
		exit(1);
	}

#ifdef OUR_IMPL

	using RationalType = long long int;
	using IntegerType = long long int;
	Graph<IntegerType, RationalType> G(filename);


#ifdef RESP_CAP
	using ArcDataType = SSPVariantRCArcData<IntegerType, RationalType>;
#else
#ifdef RESTORE_BALANCED_NODES
	using ArcDataType = SSPVariantRBNArcData<RationalType>;
#else
	using ArcDataType = BasicArcData<RationalType>;
#endif
#endif


	Network<decltype(G), IntegerType, RationalType, SSPVariantNodeData<IntegerType, RationalType>, ArcDataType> N(G, filename);

	cout << "Number of Edges: " << N.G.no_of_edges << endl;
	cout << "Number of Nodes: " << N.G.no_of_vertices << endl;

#else

	DIGRAPH_TYPEDEFS(SmartDigraph);
	SmartDigraph graph;
	IntArcMap lower(graph), capacity(graph), cost(graph);
	IntNodeMap supply(graph);

	// Read DIMACS input file
	ifstream input(filename);
	readDimacsMin(input, graph, lower, capacity, cost, supply);
	input.close();

	// compute m and n for lemon implementations
	auto arcs = 0;
	for (SmartDigraph::ArcIt a(graph); a != INVALID; ++a) {
		arcs++;
	}
	auto nodes = 0;
	for (SmartDigraph::NodeIt n(graph); n != INVALID; ++n) {
		nodes++;
	}
	cout << "Number of Edges: " << arcs << endl;
	cout << "Number of Nodes: " << nodes << endl;
#ifdef COST_SCALING
	CostScaling<SmartDigraph> N(graph);
#endif
#ifdef NETWORK_SIMPLEX
	NetworkSimplex<SmartDigraph> N(graph);
#endif
#ifdef SUCC_SHORTEST_PATH
	CapacityScaling<SmartDigraph> N(graph);
#endif
#ifdef CAP_SCALING
	CapacityScaling<SmartDigraph> N(graph);
#endif
#ifdef CYCLE_CANCELING
	CycleCanceling<SmartDigraph> N(graph);
#endif
	N.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);
#endif

    timed_job<decltype(N)> job( N, timeout );
    job.start();
    //cout << "Normal exit" << endl;
  return 0;
}
