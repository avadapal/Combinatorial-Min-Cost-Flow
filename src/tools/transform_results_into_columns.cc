#include <iostream>
#include <sstream>

using namespace std;

int main() {

  string capture;
  unsigned long n = 0;
  unsigned long m = 0;
  unsigned long C = 0;
  unsigned long U = 0;
  unsigned long b1 = 0;
  double time = 0;

  cout << "n m C U b1 time" << endl;
  
  while( !cin.eof() ) {
    cin >> capture; cerr << "!" << capture << endl;
    if( capture == "graph:" ) {
      getline( cin, capture, '_' ); cerr << "!" << capture << endl; 
      getline( cin, capture, '_' ); cerr << "#" << capture << endl;
      stringstream( capture ) >> C;
      getline( cin, capture, '_' ); cerr << "#" << capture << endl;
      stringstream( capture ) >> U;
      getline( cin, capture, '_' ); cerr << "!" << capture << endl;
      getline( cin, capture, '.' ); cerr << "#" << capture << endl;
      stringstream( capture ) >> b1;
    } else if( capture == "Edges:" ) {
      cin >> m;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> n;
      cin >> time;
      cin >> capture; cerr << "!" << capture << endl;
      cout << n << " " << m << " " << C << " " << U << " " << b1 << " " << time << endl;
    } 
  }
  
  return 0;
}