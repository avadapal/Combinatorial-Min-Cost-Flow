#include <iostream>
#include <fstream>
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/cycle_canceling.h>
#include <lemon/dimacs.h>

using namespace std;
using namespace lemon;

int main() { 

  // Create a directed graph
  //DIGRAPH_TYPEDEFS(SmartDigraph);
  //SmartDigraph g;
//   ListGraph g;
//   // Create data structures (i.e. maps) associating values to nodes and arcs of the graph
//   ListGraph::EdgeMap<int> lower(g), capacity(g), cost(g);
//   ListGraph::NodeMap<int> supply(g);


  DIGRAPH_TYPEDEFS(SmartDigraph);
  SmartDigraph g;
  
  // Create data structures (i.e. maps) associating values to nodes and arcs of the graph
  IntArcMap lower(g), capacity(g), cost(g);
  IntNodeMap supply(g);

  // Read DIMACS input file
  ifstream input("graph.txt");
  readDimacsMin(input, g, lower, capacity, cost, supply);
  input.close();

  // Initialize NetworkSimplex algorithm object
   CycleCanceling<SmartDigraph> ns(g);
  ns.lowerMap(lower).upperMap(capacity).costMap(cost).supplyMap(supply);

  // Run NetworkSimplex
  ns.run();
  

  
std::ifstream file("graph.txt");
std::string   line;
int count = 0;
vector<int> heads(1200000, 0);
vector<int> tails(1200000, 0);
while(std::getline(file, line))
{
    std::stringstream   linestream(line);
    char                data;
    int                 val1;
    int                 val2;

    // If you have truly tab delimited data use getline() with third parameter.
    // If your data is just white space separated data
    // then the operator >> will do (it reads a space separated word into a string).
    //std::getline(linestream, data, '\t');  // read up-to the first tab (discard tab).

    // Read the integers using the operator >>
    linestream >> data >> val1 >> val2;
    if(data == 'a')
    {
      count++;
      heads[count] = val1;
      tails[count] = val2;
    }
}
  
  // Print total flow cost
  printf("Total flow cost: %d\n\n", ns.totalCost());
  
 // Print flow values on the arcs
  printf("Flow values on arcs:\n");
  vector<long double> x (1200000, 0);
  for (ArcIt a(g); a != INVALID; ++a) {
    printf("Arc %d: %d/%d\n", g.id(a), ns.flow(a), capacity[a]);
     x[g.id(a)+1] = ns.flow(a);
  }
  
  vector<int> imb(1200000, 0);
  for(int i = 1; i <= 2000; i++)
  {
    imb[heads[i]] += x[i];
    imb[tails[i]] -= x[i];
  }
  
  for(int i = 1; i <= 65; i++)
  {
    cout << i << " = " <<  imb[i] << endl;
  }
  
  return 0;
}
