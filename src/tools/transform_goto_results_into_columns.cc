#include <iostream>
#include <sstream>
#include <assert.h>

using namespace std;

int main() {

  string capture;
  string graph;
  unsigned long n = 0;
  unsigned long m = 0;
  double time = 0;

  cout << "graph n m time" << endl;
  
  while( !cin.eof() ) {
    cin >> capture; cerr << "!" << capture << endl;
    if( capture == "graph:" ) {
      cin >> graph;
    } else if( capture == "Edges:" ) {
      cin >> m;
      cin >> capture; cerr << "!" << capture << endl; assert( capture == "Number" );
      cin >> capture; cerr << "!" << capture << endl; assert( capture == "of" );
      cin >> capture; cerr << "!" << capture << endl; assert( capture == "Nodes:" );
      cin >> n;
      getline( cin, capture ); 
      getline( cin, capture ); 
      if( capture == "Timeout!" ) {
	cout << graph << " " << n << " " << m << " Timeout" << endl;
      } else {
	if( capture[0] == 'f' ) {
	  getline( cin, capture );
	}
	stringstream capturestream( capture );
	string subcapture;
	if( capture[0] == 'o' ) {
	  capturestream >> subcapture; cerr << "!" << subcapture << endl; assert( subcapture == "our" );
	  capturestream >> subcapture; cerr << "!" << subcapture << endl; assert( subcapture == "implementation:" );
	}
	capturestream >> time; cerr << "#" << capture << endl;
	cin >> capture; cerr << "!" << capture << endl;
	cout << graph << " " << n << " " << m << " " << time << endl;
      }
    } 
  }
  
  return 0;
}