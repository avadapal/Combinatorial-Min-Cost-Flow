#include <iostream>
#include <sstream>

using namespace std;

int main() {

  string capture;
  string graph;
  string method;
  unsigned long n = 0;
  unsigned long m = 0;
  unsigned long steps = 0;
  unsigned long obj = 0;
  double time = 0;

  cout << "graph n m obj method steps time" << endl;
  
  while( !cin.eof() ) {
    cin >> capture; cerr << "!" << capture << endl;
    if( capture == "graph:" ) {
      cin >> graph;
    } else if ( capture[0] == '.' ) {
      stringstream capturestream( capture );
      getline( capturestream, method, '/' );
      getline( capturestream, method );
    }else if( capture == "Edges:" ) {
      cin >> m;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> n;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> steps;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> obj;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> capture; cerr << "!" << capture << endl;
      cin >> time;
      cin >> capture; cerr << "!" << capture << endl;
      cout << graph << " " << n << " " << m << " " << obj << " " << method << " " << steps << " " << time << endl;
    } 
  }
  
  return 0;
}