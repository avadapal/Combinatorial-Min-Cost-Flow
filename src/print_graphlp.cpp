#include<iostream>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>
#include "graph.h"
#include "random.h"
#include "potential_reduction.h"
#include "min_cost_flow.h"

using namespace std;


int main(){
  cout << "Enter filepath to input graph.. " << endl;
  string filepath;
  cin >> filepath; 
  int n = read_n(filepath);
  int m = read_m(filepath);
  //typedef long double T;
  
  typedef long double RationalType;
  typedef long long IntegerType;
  
  Graph<IntegerType, RationalType> G0(n);
  G0.read_graph(filepath,m);

  G0.print_lp();

  return 0;
}