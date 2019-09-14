#include <getopt.h>

#include "gridflow.h"
#include "ppmimage.h"

using namespace std;

void printUsage( char* name ) {
  cerr << "usage: " << name << " [-l] [-r] <input> <output>" << endl;
  cerr << "\t\t-l number of levels" << endl;
  cerr << "\t\t-r number of runs" << endl;
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
#ifdef WIN32
    int c = getopt(argc, (TCHAR * const *)argv, (const TCHAR *)"l:r:");
    if (c == -1) break;
    switch (c)
    {

      case 'l': l = atoi((const char *)optarg); break;
      case 'r': r = atoi((const char *)optarg); break;
      default : printUsage( argv[0] ); return -1;
    }
#else
    int c = getopt(argc, argv, "l:r:");
    if (c == -1) break;
    switch (c)
    {

      case 'l': l = atoi(optarg); break;
      case 'r': r = atoi(optarg); break;
      default : printUsage( argv[0] ); return -1;
    }
#endif
  }
  ifstream file( argv[argc-2] );
  PPMImage image( file );
  file.close();

  InstanceConversion<8> instance( image );
  bool success = false;
  clock_t start = clock();
  for( int rr = 0; rr < r; ++rr ) {
    instance.apexgrid.initialize();
    success = instance.apexgrid.solve_min_cost_flow(l);
  }
  clock_t end = clock();
  
  cout << "r = " << r << endl;
  cout << "n = " << image.n << endl;
  cout << "m = " << image.m << endl;
  cout << "height = " << instance.apexgrid.height << endl;
  cout << "width = " << instance.apexgrid.width << endl;
  cout << "blocksize = " << instance.apexgrid.blocksize << endl;
  
  for(int i=0; i<instance.apexgrid.height; i++)
  {
    for(int j = 0; j < instance.apexgrid.width; j++)
    {
      cout << instance.apexgrid.grid[i][j].down.backward_length << " , "  << instance.apexgrid.grid[i][j].down.forward_length << " | " ;
    }
    cout << endl;
  }
  
  cout << endl;
    for(int i=0; i<instance.apexgrid.height; i++)
  {
    for(int j = 0; j < instance.apexgrid.width; j++)
    {
      cout << instance.apexgrid.grid[i][j].right.backward_length << " , "  << instance.apexgrid.grid[i][j].right.forward_length << " | " ;
    }
    cout << endl;
  }
  cout << endl;
//   for(int i = 0; i < instance.apexgrid.height; i++)
//   {
//     for(int j = 0; j < instance.apexgrid.width; j++)
//     {
//        for(auto a: instance.apexgrid.grid[i][j].down.difference)
//        {
// 	 cout << a << " , ";
//        }
//        cout << endl;
//     }
//   }
  int N = image.n;
  ofstream graph_file;
  graph_file.open("graph.txt");
  graph_file << "p" << " " << "min" << " " << N << " " << N * 10 << endl;
  for(int v=1; v<= N + 1; v++)
  {
    graph_file << "n" << " " << v << " " << 0 << endl;
  }

  for(int v =1; v <= N; v++)
  {
    int i = (v-1)/8;
    int j = v % 8 - 1;
    if(v%8 == 0) j=7;
    
    graph_file << "a" << " " << N + 1 << " " << v << " " << 0 << " " << 8000000 << " " << 0 << endl;
    graph_file << "a" << " " << v << " " << N + 1 << " " << 0 << " " << 8000000 << " " << instance.apexgrid.grid[i][j].apexcost << endl;
    if(v % 8 != 0)
    {
      for(int k = 0; k < 8; k++)
      {
      graph_file << "a" << " " << v << " " << v + 1 << " " << 0 << " " << 1 << " " << instance.apexgrid.grid[i][j].right.difference[k] << endl;
      graph_file << "a" << " " << v + 1 << " " << v << " " << 0 << " " << 1 << " " << -instance.apexgrid.grid[i][j].right.difference[k] << endl;
      }
    
    }
    if(v+8 <= N)
    {
      for(int k=0; k < 8; k++)
      {
      graph_file << "a" << " " << v + 8 << " " << v  << " " << 0 << " " << 1 << " " << -instance.apexgrid.grid[i][j].down.difference[k] << endl;
      graph_file << "a" << " " << v << " " << v + 8  << " " << 0 << " " << 1 << " " << instance.apexgrid.grid[i][j].down.difference[k] << endl;
      }
    }
    
  }
  
  for(int i = 0; i < r; i++){
    
    for(int Bi = 0 ;  Bi < instance.apexgrid.height; Bi++)
    {
      for(int Bj = 0; Bj < instance.apexgrid.width; Bj++)
      {
    	      cout << instance.apexgrid.grid[Bi][Bj].deficiency << " " << " " << instance.apexgrid.grid[Bi][Bj].apexcost << " " <<
    	      instance.apexgrid.grid[Bi][Bj].potential << " " << instance.apexgrid.grid[Bi][Bj].distance << endl << endl;
      }
  
    }
    cout << instance.apexgrid.apex.deficiency << " " << instance.apexgrid.apex.apexcost << " " << instance.apexgrid.apex.potential << " " <<
    instance.apexgrid.apex.distance << " " << instance.apexgrid.Apex.x << " " << instance.apexgrid.Apex.y << " " << endl;
  }
  
  if( success ) {
    const double diff = (end - start) / static_cast<double>( CLOCKS_PER_SEC );

    cout << image.n*image.m/64 << " " << l << " " << diff << " " << instance.apexgrid.Q.pushes << endl;

    if( !instance.apexgrid.toImage() ) {
      cerr << "Failure in constructing image from instance!" << endl;
      return -1;
    }
    ofstream output( argv[argc-1] );
    output << instance.apexgrid.image;
    output.close();
  }
  else
  {
    cout<<"Failed"<<endl;
  }

  return 0;
}