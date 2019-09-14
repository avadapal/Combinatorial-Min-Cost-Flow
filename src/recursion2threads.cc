#include <pthread.h>
#include <getopt.h>

#include "gridflow.h"
#include "ppmimage.h"

using namespace std;

void printUsage( char* name ) {
  cerr << "usage: " << name << " [-l] [-r] <input> <output>" << endl;
  cerr << "\t\t-l number of levels" << endl;
  cerr << "\t\t-r number of runs" << endl;
}

struct ThreadData {
  ApexGrid<8> apexgrid;
  int levels,top,bottom,left,right;
  ThreadData( const ApexGrid<8>& g, int lev, int t, int b, int l, int r ) : apexgrid(g),levels(lev),top(t),bottom(b),left(l),right(r) {}
};

void *multithread( void *arg ) {
    ThreadData* data = (ThreadData*) arg;
    data->apexgrid.min_cost_circulation(data->levels,data->top,data->bottom,data->left,data->right);

  return arg;
}

int main( int argc, char* argv[] ) {
	
	if( argc < 3 ) {
	  printUsage( argv[0] );
		return -1;
	}

	int l = 1;
	int r = 1;
	while (true)
	{
	  int c = getopt(argc, argv, "l:r:");
	  if (c == -1) break;
			 switch (c)
			 {
			   
			   case 'l': l = atoi(optarg); break;
			   case 'r': r = atoi(optarg); break;
			   default : printUsage( argv[0] ); return -1;
			 }
	}
	ifstream file( argv[argc-2] );
	PPMImage image( file );
	file.close();
	
	ApexGrid<8> apexgrid( image );
	bool success = true;
	clock_t start = clock();
	for( int rr = 0; rr < r; ++rr ) {
	  apexgrid.initialize();
	  ThreadData data1( apexgrid, l-1, 0, apexgrid.height/2, 0, apexgrid.width );
	  ThreadData data2( apexgrid, l-1, apexgrid.height/2, apexgrid.height, 0, apexgrid.width );
	  data1.apexgrid.owner = false;
	  data2.apexgrid.owner = false;
	  pthread_t threads[2];
	  pthread_create(&threads[0], NULL, multithread, &data1 );
	  pthread_create(&threads[1], NULL, multithread, &data2 );
	  void *status;
	  pthread_join(threads[0], &status);
	  pthread_join(threads[1], &status);

	  apexgrid.apex.potential=0;
	  apexgrid.UpdatePotentialsWithMax(data1.apexgrid.apex.potential,0,apexgrid.height/2,0,apexgrid.width);
	  apexgrid.UpdatePotentialsWithMax(data2.apexgrid.apex.potential,apexgrid.height/2,apexgrid.height,0,apexgrid.width);
	  
	  
	  
	  vector <Node> Sources;
	  assert(!Sources.size());
	  
	  apexgrid.HorizontalLengthBalance(apexgrid.height/2-1,0,apexgrid.height,0,apexgrid.width);
	  apexgrid.DeficiencyInitialization(Sources,0,apexgrid.height,0,apexgrid.width);
	  
	  apexgrid.successive_shortest_path(Sources,0,apexgrid.height,0,apexgrid.width);
	  assert(!apexgrid.apex.deficiency);
// 	  cout << data1.apexgrid.Q.pushes << endl;
// 	  cout << data2.apexgrid.Q.pushes << endl;
// 	  cout << apexgrid.Q.pushes << endl;
	  
	}
	clock_t end = clock();
	if( success ) {
	  const double diff = (end - start) / static_cast<double>( CLOCKS_PER_SEC );
 	  cout << diff <<" Succeeded"<<endl;
		if( !apexgrid.toImage() ) {
			cerr << "Failure in constructing image from instance!" << endl;
			return -1;
		}

		ofstream output( argv[argc-1] );
		output << apexgrid.image;
		output.close();
	}
	else
	{	
		cout<<"Failed"<<endl;
	}
	
	return 0;
}