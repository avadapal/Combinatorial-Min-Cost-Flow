#include <getopt.h>
#ifdef WIN32
#include <tchar.h>
#endif

#include "gridflow.h"
#include "ppmimage.h"

using namespace std;

#ifndef INSTANCE_CONVERSION_FACTOR
#define INSTANCE_CONVERSION_FACTOR 8
#endif

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

  InstanceConversion<INSTANCE_CONVERSION_FACTOR> instance( image );
  bool success = false;
  clock_t start = clock();
  for( int rr = 0; rr < r; ++rr ) {
    instance.apexgrid.initialize();
    success = instance.apexgrid.solve_min_cost_flow(l);
  }
  clock_t end = clock();
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