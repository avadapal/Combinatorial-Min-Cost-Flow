/******************************************************************************* 
 * Copyright (c) 2011, George Nomikos (nomikos88@gmail.com)
 * Copyright (c) 2011, Andreas Karrenbauer (andreas.karrenbauer@uni-konstanz.de)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The name of the author may not be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *******************************************************************************/

#ifndef GRIDFLOW_H
#define GRIDFLOW_H

#include <vector>
#include <queue>
#include <algorithm>
#include <tuple>
#include <boost/range/adaptor/reversed.hpp>

#include "ppmimage.h"
#include "bucketqueue.h"

#ifdef WIN32
#include <ctime>
#endif

#define	INFINITE 8000000

using namespace std;

/** \struct ArcData
 *  \brief Contains information concerning the arcs of a node
 */
template< int N >
struct ArcData {
	int flow = 0;
	int difference[N];
	int forward_length, backward_length;
#ifdef RESP_CAP
	int xupper = 0;
	int potential = 0;
	int deficit = 0;
	int depth = 0;
	bool visited = false;
#endif
};

/** \struct Node
 * \brief Contains information about the position of a node and his distance from the source
 */
struct Node {
	int y;
	int x;
	int key;

	Node( const int& i = -1, const int& j = -1, const int& k = 0 ) : y(i), x(j), key(k) {}

	bool operator<( const Node& v ) const {
		return key > v.key;
	}

	bool operator==(const Node& other) const {
		return key == other.key;
	}
};

/** \struct NodeData
 * \brief Contains crucial information for each node in order to solve the min cost flow problem
 */
template< int N >
struct NodeData {
	int deficiency = 0;
	int apexcost = 0;
	int potential = 0;
#ifdef SSP_VARIANT
	bool visited = false;
	int depth = -1;
#ifdef USE_IMBALANCE
	int imbalance;
#endif
#else
	int distance;
#endif
	int flowtoapex;
	Node parent;
	ArcData<N> down;
	ArcData<N> right;
};

/** \struct ApexGrid
 * \brief Represent the grid flow
 */
template< int N >
struct ApexGrid {
	static const int n = N;

	PPMImage image;
  
	const int height;
	const int width;
	const int blocksize;

	int** blockmin;
	int** blockmax;

	NodeData<N>** grid;
	
	Node Apex;
	NodeData<N> apex;
	
	//priority_queue< Node > Q;
	BucketQueue< Node > Q;
	vector< Node > visited;

	int framenumber;
	int framerate;
	bool owner;
	
	ApexGrid( const int& h, const int& w ) : width(w), height(h), blocksize(N), Q(2*image.maxval+1), framenumber(100000), framerate(100000), owner(true) {
		grid = new NodeData<N>*[height];
		for( int i = 0; i < height; ++i ) {
			grid[i] = new NodeData<N>[width];
		}
	}

	ApexGrid( const PPMImage& img ) : image(img), height((image.n+N-1)/N), width((image.m+N-1)/N), blocksize(N), Q(2*image.maxval+1), framenumber(100000), framerate(100000), owner(true) {
	  blockmin = new int*[height];
	  blockmax = new int*[height];

//	  visited.reserve( height*width );

	  grid = new NodeData<N>*[height];
	  for( int i = 0; i < height; ++i ) {
	    grid[i] = new NodeData<N>[width];
	  }
	  
	  apex.deficiency = 0;
	  apex.apexcost = 0;
	  apex.potential = 0;
#ifdef SSP_VARIANT
	  apex.visited = false;
	  apex.depth = -1;
#ifdef USE_IMBALANCE
	  apex.imbalance = -1;
#endif
#else
	  apex.distance = INFINITE;
#endif
	  Apex.x = -1;
	  Apex.y = -1;
	  
	  const int maxval = image.maxval;
	  
	  for( int Bi = 0; Bi < height; ++Bi )
	  {
	    blockmin[Bi] = new int[width];
	    blockmax[Bi] = new int[width];
	  }
	  for( int Bi = 0; Bi < height; ++Bi )
	  {
	    for( int Bj = 0; Bj < width; ++Bj )
	    {
	      blockmax[Bi][Bj] = 0;
	      blockmin[Bi][Bj] = maxval;
	      for( int u = 0; u < N; ++u )
	      {
		for( int v = 0; v < N; ++v )
		{
		  const int i = Bi*N + u;
		  const int j = Bj*N + v;
		  const int r = i < image.n && j < image.m ? image.r(i,j) : 0;
		  if( r < blockmin[Bi][Bj] )
		  {
		    blockmin[Bi][Bj] = r;
		  }
		  if( r > blockmax[Bi][Bj] )
		  {
		    blockmax[Bi][Bj] = r;
		  }
		}
	      }

	      for( int u = 0; u < N; ++u )
	      {
		for( int v = 0; v < N; ++v )
		{
		  const int i = Bi*N + u;
		  const int j = Bj*N + v;
		  if( i < image.n && j < image.m ) {
		    image.r(i,j) -= blockmin[Bi][Bj];
		  }
		}
	      }
	      
	      grid[Bi][Bj].deficiency = 0;
	      grid[Bi][Bj].apexcost = maxval - blockmax[Bi][Bj] + blockmin[Bi][Bj];
	      grid[Bi][Bj].potential = 0;
#ifdef SSP_VARIANT
		  grid[Bi][Bj].visited = false;
		  grid[Bi][Bj].depth = -1;
#ifdef USE_IMBALANCE
		  grid[Bi][Bj].imbalance = -1;
#endif
#else
		  grid[Bi][Bj].distance = INFINITE;
#endif
	      
	      Node(Bj,Bi);
	      
	      if( Bj > 0 )
	      {
		for( int u = 0; u < N; ++u )
		{
		  const int i = Bi*N + u;
		  const int j = Bj*N;
		  if( i < image.n && j < image.m ) {
		    grid[Bi][Bj-1].right.difference[u] = image.r(i,j) - image.r(i,j-1);
		  }
		}
		std::sort( grid[Bi][Bj-1].right.difference, grid[Bi][Bj-1].right.difference+N );
	      } else {
		for( int u = 0; u < N; ++u )
		{
		  grid[Bi][width-1].right.difference[u] = 0;
		}
	      }
	      if( Bi > 0 )
	      {
		for( int v = 0; v < N; ++v )
		{
		  const int i = Bi*N;
		  const int j = Bj*N + v;
		  if( i < image.n && j < image.m ) {
		    grid[Bi-1][Bj].down.difference[v] = image.r(i,j) - image.r(i-1,j);
		  }
		}
		std::sort( grid[Bi-1][Bj].down.difference, grid[Bi-1][Bj].down.difference+N );
	      }else {
		for( int u = 0; u < N; ++u )
		{
		  grid[height-1][Bj].down.difference[u] = 0;
		}
	      }
	    }
	  }
// 	  toImage();
// 	  stringstream filename;
// 	  filename << "frame" << framenumber << ".pgm";
// 	  ofstream file( filename.str().c_str() );
// 	  file << image;
// 	  file.close();
// 	  ++framenumber;
// 	  toImage(-1);
	}
	
	void initialize() {
	  framenumber = 100000;
	  for( int Bi = 0; Bi < height; ++Bi )
	  {
	    for( int Bj = 0; Bj < width; ++Bj )
	    {
	      grid[Bi][Bj].deficiency = 0;
	      grid[Bi][Bj].potential = 0;
#ifdef SSP_VARIANT
		  grid[Bi][Bj].visited = false;
		  grid[Bi][Bj].depth = -1;
#ifdef USE_IMBALANCE
		  grid[Bi][Bj].imbalance = -1;
#endif
#else
		  grid[Bi][Bj].distance = INFINITE;
#endif
	    }
	  }
	  apex.deficiency = 0;
	  apex.apexcost = 0;
	  apex.potential = 0;
#ifdef SSP_VARIANT
	  apex.visited = false;
	  apex.depth = -1;
#ifdef USE_IMBALANCE
	  apex.imbalance = -1;
#endif
#else
	  apex.distance = INFINITE;
#endif
	  Apex.x = -1;
	  Apex.y = -1;
	}
	
	~ApexGrid() {
	  if( owner ) {
	  for( int Bi = 0; Bi < height; ++Bi ) {
	    if( blockmin[Bi] ) delete[] blockmin[Bi];
	    if( blockmax[Bi] ) delete[] blockmax[Bi];
	  }
	  if( blockmin )
	    delete[] blockmin;
	  if( blockmax )
	    delete[] blockmax;
	  for( int i = 0; i < height; ++i ) {
	    if( grid[i] ) delete[] grid[i];
	  }
	    if( grid )
			delete[] grid;
	  }
	}

	/**
	 * \brief Reconstructs the image based on the min cost flow solution
	 * \Return true or false, for success, failure respectively
	 */
	bool toImage( int sign = 1 ) {
	  for( int Bi = 0; Bi < height; ++Bi )
	  {
	    for( int Bj = 0; Bj < width; ++Bj )
	    {
	      const int offset = apex.potential - grid[Bi][Bj].potential;
	      for( int u = 0; u < N; ++u )
	      {
		for( int v = 0; v < N; ++v )
		{
		  const int i = Bi*N + u;
		  const int j = Bj*N + v;
		  if( i < image.n && j < image.m ) {
		    image.r(i,j) += sign*offset;
		    if( image.r(i,j) < 0  )
		    {
		      cout<< "Error: image.r(i,j) < 0"<<endl;
		      return false;
		    }
		    if(image.r(i,j) > image.maxval)
		    {
		      cout<< "Error: image.r(i,j) > image.maxval"<<endl;
		      return false;
		    }
		  }
		}
	      }
	    }
	  }
	  return true;
	}
	
	
	/**
	 * \brief Return position with the smallest non-negative arc cost otherwise the number of arcs
	 * \param differece : table with arc costs
	 * \return position
	 */
	int MinAvailableDifference(int difference[]  ) const
	{
		for(int k=0; k<N; k++)
		{
			if(difference[k]>=0)				
			{
				return k;
			}
		}
		return N;
	}
	
	/**
	 * \brief Return position in difference table according to the current flow, in order to be used for forward_length calculation
	 * \param flow : current flow
	 * \return Position
	 */
	int ReturnR(int flow) const
  	{
  		return (int)((int)(N+flow)/(int)2);
  	}
  	
  	/**
	 * \brief Return position in difference table according to the current flow, in order to be used for backward_length calculation
	 * \param flow : current flow
	 * \return Position
	 */
  	int ReturnL(int flow) const
  	{
		
  		return (int)((int)(N-flow)/(int)2);
  	}
	
	/**
	 * \brief Return forward length from a left node to rightwards or from an up node to downwards
	 * \param a : position in difference table
	 * \param u : Visited node
	 * \param adjacent: adjacent node of u
	 * \param difference: table with arc costs
	 * \return Forward length
	 */
	int ForwardLength(int a,Node u, Node adjacent, const int difference[]) const
	{
		assert(((u.x + 1 == adjacent.x) != (u.y + 1 == adjacent.y))
			&& "the adjacent node must be either to the right or below v");
		return difference[a]+grid[u.y][u.x].potential - grid[adjacent.y][adjacent.x].potential;
	}
	
	
	/**
	 * \brief Return backward length from a right node to leftwards or from a down node to upwards
	 * \param b : position in difference table
	 * \param u : Visited node
	 * \param adjacent: adjacent node of u
	 * \param difference: table with arc costs
	 * \return Forward length
	 */
	int BackwardLength(int b,Node u, Node adjacent, const int difference[]) const
	{
		assert(((u.x - 1 == adjacent.x) != (u.y - 1 == adjacent.y))
			&& "the adjacent node must be either to the left or above v");
		return -difference[N - 1 - b] + grid[u.y][u.x].potential - grid[adjacent.y][adjacent.x].potential;
	}
	
	
	/**
	 * \brief Calculates the shortest path from the apex to the rest grid nodes
	 * \param s : Source node
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the shortest path calculation takes place
	 */
	Node Source_ApexNode(const Node &s,const int top, const int bottom, const int left_bound, const int right_bound )
	{
		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
			  int cost;
			  if(grid[i][j].flowtoapex>0)
			    cost = - grid[i][j].apexcost + apex.potential - grid[i][j].potential;
			  else
			    cost= apex.potential - grid[i][j].potential;
			  if( cost == 0 && grid[i][j].deficiency < 0 ) {
			    assert( apex.distance <= image.maxval );
			    Q.push( Node( i, j, apex.distance ) );
			    grid[i][j].distance =  apex.distance;
			    grid[i][j].parent.x = -1;
			    grid[i][j].parent.y = -1;
			    return Node(i,j);
			  }
			}
		}
		for (int i=top; i<bottom; i++)
		{
		  for (int j=left_bound; j<right_bound; j++)
		  {
		    
				int cost;
				if(grid[i][j].flowtoapex>0)
					cost = - grid[i][j].apexcost + apex.potential - grid[i][j].potential;
				else
					cost= apex.potential - grid[i][j].potential;
				
				assert(cost>=0);
				const int newdistance = apex.distance + cost;
				assert( newdistance <= 2*image.maxval );
				if ( newdistance < (grid[i][j].distance) )
				{
					Q.push( Node( i, j, newdistance ) );
					grid[i][j].distance =  newdistance;
					grid[i][j].parent.x = -1;
					grid[i][j].parent.y = -1;
				}
			}
		}
		return Node();
	}
	
	/**
	 * \brief Calculates the shortest path from a node to the rest available adjacent nodes
	 * \param s: Source node
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the shortest path calculation takes place
	 * \return destination, a node with negative deficiency
	 */
	Node shortest_path( const Node &s, const int top, const int bottom, const int left_bound, const int right_bound ) 
	{
		
		int EdgeCost_u_adjacentnode=0;
		visited.clear();
		
		while (!Q.empty()) 
		{
			const Node u = Q.top();
			Q.pop();
			if(u.x==-1)
			{
				if(u.key!=apex.distance)
				{
					//Q.pop();
					continue;
				}
				if( apex.deficiency < 0 ) return Apex;
				const Node t = Source_ApexNode(s,top,bottom,left_bound,right_bound);
				if( t.x != -1 ) return t;
				continue;
			}
			
			if (grid[u.y][u.x].distance != u.key)
			{
				continue;
			}
			
			if( grid[u.y][u.x].deficiency < 0 ) return u;
			
			visited.push_back( u );
			
			if (u.y-1>=top) //Upward adjacent node
			{
				const int l=ReturnL(grid[u.y-1][u.x].down.flow);
				if(l<N )
				{
					EdgeCost_u_adjacentnode = BackwardLength(l,u, Node( u.y-1, u.x), grid[u.y-1][u.x].down.difference);
					assert(EdgeCost_u_adjacentnode>=0);
					grid[u.y-1][u.x].down.backward_length = EdgeCost_u_adjacentnode;
					UpdateDistance(Node(u.y-1,u.x),EdgeCost_u_adjacentnode,u);
				}
			}
			if(u.y+1<bottom) //Downward adjacent node
			{
				const int r=ReturnR(grid[u.y][u.x].down.flow);
				if(r<N )
				{	
					EdgeCost_u_adjacentnode = ForwardLength(r,u, Node(u.y+1,u.x), grid[u.y][u.x].down.difference);
					assert(EdgeCost_u_adjacentnode>=0);
					grid[u.y][u.x].down.forward_length = EdgeCost_u_adjacentnode;
					UpdateDistance(Node(u.y+1,u.x),EdgeCost_u_adjacentnode,u);
				}
			}
			if(u.x-1>=left_bound) //Leftward adjacent node
			{
				const int l=ReturnL(grid[u.y][u.x-1].right.flow);
				if(l<N)
				{
					EdgeCost_u_adjacentnode = BackwardLength(l,u, Node(u.y,u.x-1), grid[u.y][u.x-1].right.difference);
					assert(EdgeCost_u_adjacentnode>=0);
					grid[u.y][u.x-1].right.backward_length = EdgeCost_u_adjacentnode;
					UpdateDistance(Node(u.y,u.x-1),EdgeCost_u_adjacentnode,u);
					
				}
			}
			if(u.x+1<right_bound) //Rightward adjacent node 
			{
				const int r=ReturnR(grid[u.y][u.x].right.flow);
				if(r<N)
				{
					EdgeCost_u_adjacentnode = ForwardLength(r,u, Node(u.y,u.x+1), grid[u.y][u.x].right.difference);
					assert(EdgeCost_u_adjacentnode>=0);
					grid[u.y][u.x].right.forward_length = EdgeCost_u_adjacentnode;
					UpdateDistance(Node(u.y,u.x+1),EdgeCost_u_adjacentnode,u);
					
				}
			}
			
			//Apex adjacent node
			if(grid[u.y][u.x].flowtoapex<0)
				UpdateDistanceApex( grid[u.y][u.x].potential - apex.potential,u);
			else
				UpdateDistanceApex( grid[u.y][u.x].apexcost + grid[u.y][u.x].potential - apex.potential,u);
			
		}
		return s;
	}
	
	/**
	 * \brief Updates the distance of adjacent node from the source
	 * \param EdgeCost_u_adjacentnode: Cost between current visited node and his adjacent node
	 * \param u: Current visited node
	 */
	inline void UpdateDistanceApex( const int& EdgeCost_u_adjacentnode, const Node& u)
	{
		
		const int distance = grid[u.y][u.x].distance + EdgeCost_u_adjacentnode ;
		
		assert(distance>=0);
		
		if ( distance < apex.distance)
			{				
				Apex.key=distance;
				assert( distance <= 2*image.maxval );
				Q.push(Apex);	
				apex.distance =  distance;
				apex.parent.x = u.x;
				apex.parent.y = u.y;
			}
		return;
	}

	/**
	 * \brief Updates the distance of adjacent node from the source
	 * \param adjacent: Adjacent node
	 * \param EdgeCost_u_adjacentnode: Cost between current visited node and his adjacent node
	 * \param u: Current visited node
	 */
	inline void UpdateDistance( const Node& adjacent, const int& EdgeCost_u_adjacentnode, const Node& u)
	{
	  
	  const int distance = grid[u.y][u.x].distance + EdgeCost_u_adjacentnode ;
	  
	  assert(distance>=0);
	  
	  if( distance < grid[adjacent.y][adjacent.x].distance)
	    {
	      assert( distance <= 2*image.maxval );
	      Q.push(Node(adjacent.y,adjacent.x,distance));
	      grid[adjacent.y][adjacent.x].distance =  distance;
	      grid[adjacent.y][adjacent.x].parent.x = u.x;
	      grid[adjacent.y][adjacent.x].parent.y = u.y;
	    }

	    return;
	}
	
	/**
	 * \brief Updates the flow of nodes that consist the shortest path
	 * \param source: source node
	 * \param destination: destination node
	 */
	void UpdateFlowPath(Node source,Node destination)
	{
		Node w = destination;
		
		//Traverse the path in reverse order
		while((w.x!=source.x) || (w.y!=source.y))
		{
			if(w.x==-1)
			{
				const Node v = apex.parent;
				grid[v.y][v.x].flowtoapex++;
				w = v;
			}
			else
			{
				const Node v = grid[w.y][w.x].parent;
				if( v.x == -1 ) {
				  grid[w.y][w.x].flowtoapex--;
				}
				//Rightward
				else if( v.x<w.x && v.y==w.y )
				{
				  grid[v.y][v.x].right.flow++;
				}
				//Downward
				else if(v.x==w.x && v.y<w.y)
				{
				  grid[v.y][v.x].down.flow++;
				}
				//Leftward
				else if(v.x>w.x && v.y==w.y)
				{
				  grid[w.y][w.x].right.flow--;
				}
				//Upward
				else if(v.x==w.x && v.y>w.y)
				{
				  grid[w.y][w.x].down.flow--;
				}
				else
				{
				  cout << "--> Error: wrong coordinates during path traversal <--"<<endl;
				}
				w = v;
			}
		}
		
		return;
	}
	
	
	/**
	 * \brief Initializes the environment for a shortest path calculation, sends a ource to the shortest path calculation, updates information for each grid node concering potentials and deficiencies
	 * \param sources: Vector with available sources
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the sources and destination belong to
	 * \return the maximum apex potential
	 */
	int successive_shortest_path(vector< Node >& sources, const int top, const int bottom, const int left_bound, const int right_bound ) 
	{	
		while(!sources.empty())
		{
			
			Node source = sources.back();
			
			if(source.y>=top && source.y<bottom && source.x>=left_bound && source.x<right_bound)
			{
				if(grid[source.y][source.x].deficiency>0)
				{
					
					grid[source.y][source.x].deficiency-=1;
					grid[source.y][source.x].distance = 0;
					Q.push(Node(source.y,source.x,0));
					
					Node dest = shortest_path(source,top,bottom,left_bound,right_bound);
					
					//Apex node as destination
					if(dest.x==-1)
					{
						assert(apex.deficiency<0);
						apex.deficiency+=1;
						UpdateFlowPath(source,Apex);
						const int distt = apex.distance;
						for( vector< Node >::iterator iter = visited.begin(), end = visited.end(); iter != end; ++iter ) {
							const Node v = *iter;
							if( (grid[v.y][v.x].distance) < INFINITE ) {
								grid[v.y][v.x].potential += grid[v.y][v.x].distance - distt;
							}
							grid[v.y][v.x].distance = INFINITE;
							grid[v.y][v.x].parent.x=v.x;
							grid[v.y][v.x].parent.y=v.y;
							grid[v.y][v.x].parent.key=INFINITE;
							Node(v.y,v.x,INFINITE);
							
						}
						while( !Q.empty() ) {
							const Node v = Q.top();
							grid[v.y][v.x].distance = INFINITE;
							grid[v.y][v.x].parent.x=v.x;
							grid[v.y][v.x].parent.y=v.y;
							grid[v.y][v.x].parent.key=INFINITE;
							Node(v.y,v.x,INFINITE);
							Q.pop();
						}
						
						
					}
					else
					{
						assert(grid[dest.y][dest.x].deficiency<0);
						grid[dest.y][dest.x].deficiency+=1;
						UpdateFlowPath(source,dest);
						const int distt = grid[dest.y][dest.x].distance;
						for( vector< Node >::iterator iter = visited.begin(), end = visited.end(); iter != end; ++iter ) {
							const Node v = *iter;
							if( (grid[v.y][v.x].distance) < INFINITE ) {
								grid[v.y][v.x].potential += grid[v.y][v.x].distance - distt;
							}
							grid[v.y][v.x].distance = INFINITE;
							grid[v.y][v.x].parent.x=v.x;
							grid[v.y][v.x].parent.y=v.y;
							grid[v.y][v.x].parent.key=INFINITE;
							Node(v.y,v.x,INFINITE);
						}
						while( !Q.empty() ) {
							const Node v = Q.top();
							if( v.x != -1 ) {
								grid[v.y][v.x].distance = INFINITE;
								grid[v.y][v.x].parent.x=v.x;
								grid[v.y][v.x].parent.y=v.y;
								grid[v.y][v.x].parent.key=INFINITE;
								Node(v.y,v.x,INFINITE);
							}
							Q.pop();
						}
						if( (apex.distance) < distt ) {
							apex.potential += apex.distance - distt;
						}
						grid[dest.y][dest.x].distance = INFINITE;
						grid[dest.y][dest.x].parent.x=dest.x;
						grid[dest.y][dest.x].parent.y=dest.y;
						grid[dest.y][dest.x].parent.key=INFINITE;
						Node(dest.y,dest.x,INFINITE);
						
					}
					
					apex.distance = INFINITE;
					apex.parent.x = -1;
					apex.parent.y = -1;
					Apex.key = INFINITE;

// 					if( framenumber%framerate == 0 ) {
// 					  toImage();
// 					  stringstream filename;
// 					  filename << "frame" << framenumber << ".pgm";
// 					  ofstream file( filename.str().c_str() );
// 					  file << image;
// 					  file.close();
// 					  toImage(-1);
// 					}
// 					++framenumber;
				}
				else
				{
					sources.pop_back();		
				}
			}
			else
			{
				cout<<"--> Error: Out of bounds source <--" <<endl;
				return apex.potential;
			}
		}
		
		return apex.potential;
	}


	auto is_apex(const Node &v) const -> bool {
		return v.x == -1;
	}

	auto get_nodedata(const Node &v) -> decltype(apex) & {
		return is_apex(v) ? apex : grid[v.y][v.x];
	}

	// v and w must be neighbors and neither v or w may be the apex
	// returns a reference to the arcdata between the two nodes and a bool indicating the direction
	// in which the edge is stored: true if the edge is stored from v to w, false otherwise
	auto get_arcdata(const Node &v, const Node &w) const -> std::pair<decltype(grid[v.y][v.x].right)&, bool> {
		assert(!is_apex(v) && !is_apex(w));
		if (v.x < w.x)
			return {grid[v.y][v.x].right, true};
		if (v.y < w.y)
			return {grid[v.y][v.x].down, true};
		if (v.x > w.x)
			return {grid[w.y][w.x].right, false};
		if (v.y > w.y)
			return {grid[w.y][w.x].down, false};
		assert(false && "unreachable unless the assumptions of this functions are violated in the call");
		exit(1);
	}

	template<typename UnaryFunction>
	void for_each_nodedata(UnaryFunction f, int top, int bottom, int left_bound, int right_bound) {
		for (auto i = top; i < bottom; ++i)
			for (auto j = left_bound; j < right_bound; ++j)
				f(grid[i][j]);
		f(apex);
	}

	template<typename UnaryFunction>
	void for_each_node(UnaryFunction f, int top, int bottom, int left_bound, int right_bound) {
		for (auto i = top; i < bottom; ++i)
			for (auto j = left_bound; j < right_bound; ++j)
				f(Node(i, j));
		f(Apex);
	}

	void update_x_by_tree(const std::vector<std::vector<std::pair<Node, Node>>> &tree, int top, int bottom, int left_bound, int right_bound) {
#ifdef USE_IMBALANCE
		// initially set imbalance to the deficiency of each node
		for_each_nodedata([](decltype(apex) &nodedata) {nodedata.imbalance = nodedata.deficiency; },
			top, bottom, left_bound, right_bound);
#endif

		for (const auto &depth : boost::adaptors::reverse(tree)) {
			for (const auto &edge : depth) {
				Node v, w;
				std::tie(v, w) = edge;
				auto forward = false;

				if (get_nodedata(v).depth < get_nodedata(w).depth) {
					std::swap(v, w);
					forward = true;
				}

				assert(get_nodedata(v).depth > get_nodedata(w).depth && "in the tree w is v's parent");

				auto &data_v = get_nodedata(v);
				auto &data_w = get_nodedata(w);

#ifdef USE_IMBALANCE
				// update imbalance
				data_w.imbalance += data_v.imbalance;
#endif

				if (is_apex(v)) {
#ifdef USE_IMBALANCE
					if (data_v.imbalance < 0 && (forward == (data_w.flowtoapex >= 0))) {
						auto flow_incr = forward ? data_v.imbalance : std::max(data_w.flowtoapex, data_v.imbalance);
#else
					if (data_v.deficiency < 0 && (forward == (data_w.flowtoapex >= 0))) {
						auto flow_incr = forward ? data_v.deficiency : std::max(data_w.flowtoapex, data_v.deficiency);
#endif
						// send flow to apex
						assert(data_w.flowtoapex < 0 ? data_v.potential - data_w.potential == 0 : data_w.apexcost + data_w.potential - data_v.potential == 0);
						assert(flow_incr < 0);
						data_w.flowtoapex -= flow_incr;
						data_v.deficiency -= flow_incr;
						data_w.deficiency += flow_incr;
					}
#ifdef USE_IMBALANCE
					if (data_v.imbalance > 0 && (forward != (data_w.flowtoapex <= 0))) {
						auto flow_incr = !forward ? data_v.imbalance : std::min(data_w.flowtoapex, data_v.imbalance);
#else
					if (data_v.deficiency > 0 && (forward != (data_w.flowtoapex <= 0))) {
						auto flow_incr = !forward ? data_v.deficiency : std::min(data_w.flowtoapex, data_v.deficiency);
#endif
						// send flow away from apex
						assert(data_w.flowtoapex > 0 ? data_w.apexcost + data_w.potential - data_v.potential == 0 : data_v.potential - data_w.potential == 0);
						assert(flow_incr > 0);
						data_w.flowtoapex -= flow_incr;
						data_v.deficiency -= flow_incr;
						data_w.deficiency += flow_incr;
					}
					continue;
				}

				if (is_apex(w)) {
#ifdef USE_IMBALANCE
					if (data_v.imbalance < 0 && (forward == (data_v.flowtoapex <= 0))) {
						auto flow_incr = forward ? data_v.imbalance : std::max(-data_v.flowtoapex, data_v.imbalance);
#else
					if (data_v.deficiency < 0 && (forward == (data_v.flowtoapex <= 0))) {
						auto flow_incr = forward ? data_v.deficiency : std::max(-data_v.flowtoapex, data_v.deficiency);
#endif
						// send flow away from apex
						assert(data_v.flowtoapex > 0 ? data_v.apexcost + data_v.potential - data_w.potential == 0 : data_w.potential - data_v.potential == 0);
						assert(flow_incr < 0);
						data_v.flowtoapex += flow_incr;
						data_v.deficiency -= flow_incr;
						data_w.deficiency += flow_incr;
					}
#ifdef USE_IMBALANCE
					if (data_v.imbalance > 0 && (forward != (data_v.flowtoapex >= 0))) {
						auto flow_incr = !forward ? data_v.imbalance : std::min(-data_v.flowtoapex, data_v.imbalance);
#else
					if (data_v.deficiency > 0 && (forward != (data_v.flowtoapex >= 0))) {
						auto flow_incr = !forward ? data_v.deficiency : std::min(-data_v.flowtoapex, data_v.deficiency);
#endif
						// send flow to apex
						assert(data_v.flowtoapex < 0 ? data_w.potential - data_v.potential == 0 : data_v.apexcost + data_v.potential - data_w.potential == 0);
						assert(flow_incr > 0);
						data_v.flowtoapex += flow_incr;
						data_v.deficiency -= flow_incr;
						data_w.deficiency += flow_incr;
					}
					continue;
				}

				// neither v nor w are the apex
				auto arcdata = get_arcdata(v, w);
				auto x_old = arcdata.first.flow;

				// update flow
				if (arcdata.second) {
					// edge from v to w
#ifdef USE_IMBALANCE
					if (data_v.imbalance < 0 && arcdata.first.flow > -N && !forward) {
#else
					if (data_v.deficiency < 0 && arcdata.first.flow > -N && !forward) {
#endif
						assert(BackwardLength(ReturnL(arcdata.first.flow), w, v, arcdata.first.difference) == 0);
						arcdata.first.flow -= 1;
					}
#ifdef USE_IMBALANCE
					if (data_v.imbalance > 0 && arcdata.first.flow < N && forward) {
#else
					if (data_v.deficiency > 0 && arcdata.first.flow < N && forward) {
#endif
						assert(ForwardLength(ReturnR(arcdata.first.flow), v, w, arcdata.first.difference) == 0);
						arcdata.first.flow += 1;
					}

					// update deficiency
					auto flow_incr = arcdata.first.flow - x_old;
					data_v.deficiency -= flow_incr;
					data_w.deficiency += flow_incr;
				} else {
					// edge from w to v
#ifdef USE_IMBALANCE
					if (data_v.imbalance < 0 && arcdata.first.flow < N && forward) {
#else
					if (data_v.deficiency < 0 && arcdata.first.flow < N && forward) {
#endif
						assert(ForwardLength(ReturnR(arcdata.first.flow), w, v, arcdata.first.difference) == 0);
						arcdata.first.flow += 1;
					}
#ifdef USE_IMBALANCE
					if (data_v.imbalance > 0 && arcdata.first.flow > -N && !forward) {
#else
					if (data_v.deficiency > 0 && arcdata.first.flow > -N && !forward) {
#endif
						assert(BackwardLength(ReturnL(arcdata.first.flow), v, w, arcdata.first.difference) == 0);
						arcdata.first.flow -= 1;
					}

					// update deficiency
					auto flow_incr = arcdata.first.flow - x_old;
					data_v.deficiency += flow_incr;
					data_w.deficiency -= flow_incr;
				}
			}
		}
	}

	template<typename QueueType>
	void push(QueueType &queue, int key, const std::pair<Node, Node> &edge) {
		queue.push({key, edge});
	}

	template<typename ElementType>
	void push(BucketQueue<ElementType> &queue, int key, const std::pair<Node, Node> &edge) {
		assert(std::abs(key) <= image.maxval);
		queue.push({key + (image.maxval), edge});
	}

	template<typename QueueType>
	void update_s_in_s_out(Node v, QueueType &s_in, QueueType &s_out, int top, int bottom, int left_bound, int right_bound) {
		if (is_apex(v)) {
			for (auto i = top; i < bottom; ++i) {
				for (auto j = left_bound; j < right_bound; ++j) {
					const auto &nodedata = grid[i][j];

					if (nodedata.visited)
						continue;

					if (nodedata.flowtoapex >= 0) {
						push(s_in, nodedata.apexcost + nodedata.potential - apex.potential, {Node(i, j), v});
					} else {
						push(s_in, nodedata.potential - apex.potential, {v, Node(i, j)});
					}

					if (nodedata.flowtoapex <= 0) {
						push(s_out, apex.potential - nodedata.potential, {v, Node(i, j)});
					} else {
						push(s_out, -nodedata.apexcost - nodedata.potential + apex.potential, {Node(i, j), v});
					}
				}
			}
			return;
		}

		// up
		if (v.y - 1 >= top && !grid[v.y - 1][v.x].visited) {
			Node w(v.y - 1, v.x); // the node above v
			const auto &arcdata = get_nodedata(w).down;
			if (arcdata.flow > -N) {
				// can send flow upward
				auto length = BackwardLength(ReturnL(arcdata.flow), v, w, arcdata.difference);
				push(s_out, length, {w, v});
			}
			if (arcdata.flow < N) {
				// can send flow downward
				auto length = ForwardLength(ReturnR(arcdata.flow), w, v, arcdata.difference);
				push(s_in, length, {v, w});
			}
		}

		// down
		if (v.y + 1 < bottom && !grid[v.y + 1][v.x].visited) {
			Node w(v.y + 1, v.x); // the node below w
			const auto &arcdata = get_nodedata(v).down;
			if (arcdata.flow > -N) {
				// can send flow upward
				auto length = BackwardLength(ReturnL(arcdata.flow), w, v, arcdata.difference);
				push(s_in, length, {w, v});
			}
			if (arcdata.flow < N) {
				// can send flow downward
				auto length = ForwardLength(ReturnR(arcdata.flow), v, w, arcdata.difference);
				push(s_out, length, {v, w});
			}
		}

		// left
		if (v.x - 1 >= left_bound && !grid[v.y][v.x - 1].visited) {
			Node w(v.y, v.x - 1); // the node to the left of v
			const auto &arcdata = get_nodedata(w).right;
			if (arcdata.flow > -N) {
				// can send flow to the left
				auto length = BackwardLength(ReturnL(arcdata.flow), v, w, arcdata.difference);
				push(s_out, length, {w, v});
			}
			if (arcdata.flow < N) {
				// can send flow to the right
				auto length = ForwardLength(ReturnR(arcdata.flow), w, v, arcdata.difference);
				push(s_in, length, {v, w});
			}
		}

		// right
		if (v.x + 1 < right_bound && !grid[v.y][v.x + 1].visited) {
			Node w(v.y, v.x + 1); // the node to the right of v
			const auto &arcdata = get_nodedata(v).right;
			if (arcdata.flow > -N) {
				// can send flow to the left
				auto length = BackwardLength(ReturnL(arcdata.flow), w, v, arcdata.difference);
				push(s_in, length, {w, v});
			}
			if (arcdata.flow < N) {
				// can send flow to the right
				auto length = ForwardLength(ReturnR(arcdata.flow), v, w, arcdata.difference);
				push(s_out, length, {v, w});
			}
		}

		if (!apex.visited) {
			const auto &nodedata = grid[v.y][v.x];

			if (nodedata.flowtoapex <= 0) {
				push(s_in, apex.potential - nodedata.potential, {Apex, v});
			} else {
				push(s_in, -nodedata.apexcost - nodedata.potential + apex.potential, {v, Apex});
			}

			if (nodedata.flowtoapex >= 0) {
				push(s_out, nodedata.apexcost + nodedata.potential - apex.potential, {v, Apex});
			} else {
				push(s_out, nodedata.potential - apex.potential, {Apex, v});
			}
		}
	}

	/**
	* \brief SSP variant, updates potentials and deficiencies
	* \param sources: Vector with available sources
	* \param bottom,top,left_bound,right_bound : Coordinates of grid in which the sources and destination belong to
	* \return the maximum apex potential
	*/
	int successive_shortest_path_variant(vector<Node> &sources, const int top, const int bottom, const int left_bound, const int right_bound) {
		const auto num_nodes = (bottom - top) * (right_bound - left_bound) + 1;
		
		int step_count = 0;
		while (!sources.empty()) {
			int deficit_1 = 0;
			for_each_nodedata([&deficit_1](const decltype(apex) &nodedata) { deficit_1 += std::abs(nodedata.deficiency); },
				top, bottom, left_bound, right_bound);
			std::cout << "step " << ++step_count << ", deficit_1 = " << deficit_1 << std::endl;

			//auto v_start = *std::max_element(std::begin(sources), std::end(sources) [   this](const Node &lhs, const Node &rhs) { return get_nodedata(lhs).deficiency < get_nodedata(rhs).deficiency || lhs < rhs; });
			auto v_start = Apex;
			auto &v_start_data = get_nodedata(v_start);

			//assert(v_start_data.deficiency > 0);

			auto deficit_s = v_start_data.deficiency;

#ifdef USE_BUCKET_QUEUE
			struct PriorityQueueElementType {
				int key;
				std::pair<Node, Node> edge;
			};

			auto s_in = BucketQueue<PriorityQueueElementType>(2 * image.maxval + 1);
			auto s_out = BucketQueue<PriorityQueueElementType>(2 * image.maxval + 1);
#else
			// (cost, (source, target))
			using PriorityQueueElementType = std::pair<decltype(v_start_data.potential), std::pair<Node, Node>>;
			//auto compare_first = [](const PriorityQueueElementType &lhs, const PriorityQueueElementType &rhs) {
			//	return lhs.first > rhs.first;
			//};
			auto compare_first = [](const PriorityQueueElementType &lhs, const PriorityQueueElementType &rhs) {
				if (lhs.first != rhs.first)
					return lhs.first > rhs.first;
				if (lhs.second.first.x != rhs.second.first.x)
					return lhs.second.first.x > rhs.second.first.x;
				if (lhs.second.first.y != rhs.second.first.y)
					return lhs.second.first.y > rhs.second.first.y;
				if (lhs.second.second.x != rhs.second.second.x)
					return lhs.second.second.x > rhs.second.second.x;
				return lhs.second.second.y > rhs.second.second.y;
			};
			auto s_in = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);
			auto s_out = std::priority_queue<PriorityQueueElementType, std::vector<PriorityQueueElementType>, decltype(compare_first)>(compare_first);
#endif

			auto num_visited = 1;

			v_start_data.depth = 0;
			v_start_data.visited = true;

			auto tree = std::vector<std::vector<std::pair<Node, Node>>>();

			update_s_in_s_out(v_start, s_in, s_out, top, bottom, left_bound, right_bound);

			while (num_visited < num_nodes && !(s_in.empty() && s_out.empty())) {
				decltype(v_start_data.potential) s_a_hat;
				Node v, w;
				int delta;
				bool forward = true;

				if (deficit_s > 0 || s_in.empty()) {
					// flow is supposed to go out of S
					assert(!s_out.empty());
#ifdef USE_BUCKET_QUEUE
					auto queue_element = s_out.top();
					s_a_hat = queue_element.key - (image.maxval);
					std::tie(v, w) = queue_element.edge;
#else
					std::forward_as_tuple(s_a_hat, std::tie(v, w)) = s_out.top();
#endif
					//std::cout << "from s_out: (" << v.x << ", " << v.y << "), (" << w.x << ", " << w.y << ")" << std::endl;
					s_out.pop();
					delta = s_a_hat;
				} else {
					// flow is supposed to go into S
					assert(!s_in.empty());
#ifdef USE_BUCKET_QUEUE
					auto queue_element = s_in.top();
					s_a_hat = queue_element.key - (image.maxval);
					std::tie(v, w) = queue_element.edge;
#else
					std::forward_as_tuple(s_a_hat, std::tie(v, w)) = s_in.top();
#endif
					//std::cout << "from s_in: (" << v.x << ", " << v.y << "), (" << w.x << ", " << w.y << ")" << std::endl;
					s_in.pop();
					delta = -s_a_hat;
				}

				if (!get_nodedata(v).visited) {
					std::swap(v, w);
					forward = false;
				}

				auto &data_v = get_nodedata(v);
				auto &data_w = get_nodedata(w);

				assert(data_v.visited);

				if (data_w.visited)
					// the edge is not a cut edge (anymore)
					continue;

				auto depth = data_v.depth + 1;
				data_w.depth = depth;
				if (tree.size() < depth)
					tree.push_back(std::vector<std::pair<Node, Node>>());
				tree[depth - 1].push_back(forward ? make_pair(v, w) : make_pair(w, v));

				data_w.potential += delta;

#ifdef UPDATE_S_CUTS
				// update x[a_hat]
				auto x_old = N.arcdata[std::abs(a_hat)].xlower;

				update_x_by_s_cuts(N, a_hat, deficit_s, forward);

				auto flow_incr = N.arcdata[std::abs(a_hat)].xlower - x_old;
#endif

				update_s_in_s_out(w, s_in, s_out, top, bottom, left_bound, right_bound);

				// put w in S
				deficit_s += data_w.deficiency;
				data_w.visited = true;
				++num_visited;

#ifdef UPDATE_S_CUTS
				// update deficit vector
				if (out_it_goes == forward) {
					deficit[v] += flow_incr;
					deficit[w] -= flow_incr;
				} else {
					deficit[v] -= flow_incr;
					deficit[w] += flow_incr;
				}
#endif
			}

#ifndef UPDATE_S_CUTS
			update_x_by_tree(tree, top, bottom, left_bound, right_bound);
#endif

			// update A_n
			sources.clear();
			for_each_node([this, &sources](const Node &v) {if (get_nodedata(v).deficiency > 0) sources.push_back(v); },
				top, bottom, left_bound, right_bound);

			// reset visited flag
			for_each_nodedata([](decltype(apex) &nodedata) {nodedata.visited = false; },
				top, bottom, left_bound, right_bound);
		}


		// check flow conservation
		auto check_flow_conservation = [this, top, bottom, left_bound, right_bound](const Node &v) {
			auto flow = 0;
			if (is_apex(v)) {
				for (auto i = top; i < bottom; ++i)
					for (auto j = left_bound; j < right_bound; ++j)
						flow -= grid[i][j].flowtoapex;
			} else {
				if (v.y - 1 >= top)
					flow -= get_nodedata(Node(v.y - 1, v.x)).down.flow;
				if (v.y + 1 < bottom)
					flow += get_nodedata(v).down.flow;
				if (v.x - 1 >= left_bound)
					flow -= get_nodedata(Node(v.y, v.x - 1)).right.flow;
				if (v.x + 1 < right_bound)
					flow += get_nodedata(v).right.flow;
				flow += get_nodedata(v).flowtoapex;
			}
			assert(flow == 0 && "flow conservation check failed");
		};

		for_each_node(check_flow_conservation, top, bottom, left_bound, right_bound);


		// check capacity constraints and complementary slackness
		auto check_capacity_constraints_and_complementary_slackness = [](const ArcData<N> &arcdata, int d_a) {
			assert(-N <= arcdata.flow && arcdata.flow <= N && "capacity constraints");
			for (auto a = 0; a < N; ++a) {
				auto s_a = arcdata.difference[a] + d_a;
				auto flow = 0;
				if (arcdata.flow > 2 * a - N + 1)
					flow = 1;
				if (arcdata.flow < 2 * a - N + 1)
					flow = -1;
				assert((s_a <= 0 || flow == -1) && "complementary slackness");
				assert((s_a >= 0 || flow == 1) && "complementary slackness");
			}
		};

		for (auto i = top; i < bottom; ++i) {
			for (auto j = left_bound; j < right_bound; ++j) {
				const auto &nodedata = grid[i][j];

				// down
				if (i < bottom - 1)
					check_capacity_constraints_and_complementary_slackness(nodedata.down, nodedata.potential - grid[i + 1][j].potential);

				// right
				if (j < right_bound - 1)
					check_capacity_constraints_and_complementary_slackness(nodedata.right, nodedata.potential - grid[i][j + 1].potential);

				// apex
				// no capacity constraints
				auto s_a_1 = nodedata.apexcost + nodedata.potential - apex.potential;
				auto s_a_2 = apex.potential - nodedata.potential;
				assert((s_a_1 <= 0 || nodedata.flowtoapex <= 0) && "complementary slackness");
				assert((s_a_2 <= 0 || nodedata.flowtoapex >= 0) && "complementary slackness");
				assert((nodedata.flowtoapex <= 0 || s_a_1 == 0) && "complementary slackness");
				assert((nodedata.flowtoapex >= 0 || s_a_2 == 0) && "complementary slackness");
			}
		}

		return apex.potential;
	}

	/**
	 * \brief Initializes the flow for each grid node
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the flow initialization takes place
	 */
	void FlowInitialization(const int top, const int bottom, const int left_bound, const int right_bound)
	{
		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
				grid[i][j].flowtoapex = 0;
				
				int r = MinAvailableDifference(grid[i][j].right.difference);
				grid[i][j].right.flow=(2*r)-N;
				
				r =MinAvailableDifference(grid[i][j].down.difference);
				grid[i][j].down.flow=(2*r)-N;
			}
		}
		return;
	}
	
	/**
	 * \brief Initializes the deficiency for each grid node, choosing the available sources
	 * \param Sources: Nodes with positive deficiency
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the deficiency initialization takes place
	 */
	bool DeficiencyInitialization(vector <Node> &Sources,const int top, const int bottom, const int left_bound, const int right_bound)
	{
		int excess = 0;
		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
				int total_deficiency=0;
				//Upward
				if(i-1>=top)
				{
					total_deficiency+= grid[i-1][j].down.flow;
				}
				//Downward
				if(i+1<bottom)
				{
					total_deficiency-= grid[i][j].down.flow;
				}
				//Leftward
				if(j-1>=left_bound)
				{
					total_deficiency+= grid[i][j-1].right.flow;
				}
				//Rightwarad
				if(j+1<right_bound)
				{
					total_deficiency-= grid[i][j].right.flow;
				}
				
				total_deficiency -= grid[i][j].flowtoapex;
				
				grid[i][j].deficiency = total_deficiency;
				if(total_deficiency>0)
				{
					Sources.push_back(Node(i,j));
					excess += total_deficiency;
				}
			}
		}
		//cout << excess << endl;
		return true;
	}
	
	/**
	 * \brief Starts the algorithm to solve the min cost flow problem
	 * \return true for successfull solution of min cost flow problem
	 */
	bool solve_min_cost_flow( const int& levels ) 
	{
		min_cost_circulation(levels,0,height,0,width);

		assert(check(0,height,0,width));
		assert(CheckDeficiencies(0,height,0,width));
		return true;
	}
	
	
	
	
	int success(int top,int bottom, int left_bound, int right_bound)
	{
		cout << "-----"<<endl;
		
		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
				cout<<"("<<i<<","<<j<<") ";
			}
			cout << endl;
		}
		cout << "-----"<<endl;
		return right_bound;
	}
	
	/**
	 * \brief Solves normally or recursively the min cost flow problem
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the min cost flow problem is solved
	 * \param temp_levels : determines the level of recursion
	 * \return the maximum apex potential
	 */
	int min_cost_circulation(int temp_levels,int top, int bottom,  int left_bound, int right_bound)
	{
		if( ( (right_bound-left_bound)>2 || (bottom-top)>2) && (temp_levels >1))
		{
			
			//Width>height, grid division by width dimension
			if((right_bound-left_bound)>(bottom-top))
			{
				int half_width = (left_bound + right_bound)/2;
				
				apex.potential=0;
				int potential_1 = min_cost_circulation(temp_levels-1,top,bottom,left_bound,half_width);
				
				
				apex.potential=0;
				int potential_2 = min_cost_circulation(temp_levels-1,top,bottom,half_width,right_bound);
				
				
				apex.potential=0;	
				UpdatePotentialsWithMax(potential_1,top,bottom,left_bound,half_width);
				UpdatePotentialsWithMax(potential_2,top,bottom,half_width,right_bound);
				
				
				
				vector <Node> Sources;
				assert(!Sources.size());
				
				VerticalLengthBalance(half_width-1,top,bottom,left_bound,right_bound);
				DeficiencyInitialization(Sources,top,bottom,left_bound,right_bound);
				
#ifdef SSP_VARIANT
				auto pot = successive_shortest_path_variant(Sources, top, bottom, left_bound, right_bound);
#else
				int pot = successive_shortest_path(Sources,top,bottom,left_bound,right_bound);
#endif
				assert(!apex.deficiency);
				return pot;
			}
			
			else 
			{
				int half_height = (top + bottom)/2;
				
				apex.potential=0;
				int potential_1 = min_cost_circulation(temp_levels-1,top,half_height,left_bound,right_bound);
				
				
				apex.potential=0;
				int potential_2 = min_cost_circulation(temp_levels-1,half_height,bottom,left_bound,right_bound);
				
				
				apex.potential=0;
				UpdatePotentialsWithMax(potential_1,top,half_height,left_bound,right_bound);
				UpdatePotentialsWithMax(potential_2,half_height,bottom,left_bound,right_bound);
				
				
				
				vector <Node> Sources;
				assert(!Sources.size());
				
				HorizontalLengthBalance(half_height-1,top,bottom,left_bound,right_bound);
				DeficiencyInitialization(Sources,top,bottom,left_bound,right_bound);
				
#ifdef SSP_VARIANT
				auto pot = successive_shortest_path_variant(Sources, top, bottom, left_bound, right_bound);
#else
				int pot = successive_shortest_path(Sources, top, bottom, left_bound, right_bound);
#endif
				assert(!apex.deficiency);
				return pot;
			}
		}
		else
		{
			vector <Node> Sources;
			assert(!Sources.size());
			FlowInitialization(top,bottom,left_bound,right_bound);
			DeficiencyInitialization(Sources,top,bottom,left_bound,right_bound);
			apex.deficiency = -TotalDeficiency(top,bottom,left_bound,right_bound);

#ifdef ONLY_CONVERT_TO_DIMACS
			ofstream dimacs_out("image.dimacs");
			graph_as_dimacs(dimacs_out);
			dimacs_out.close();
			exit(0);
#endif

#ifdef SSP_VARIANT
			return successive_shortest_path_variant(Sources, top, bottom, left_bound, right_bound);
#else
			return successive_shortest_path(Sources,top,bottom,left_bound,right_bound);
#endif
		}
	}


	void graph_as_dimacs(std::ostream& os) {
		auto num_nodes = width * height + 1;
		auto num_edges = 0;
		num_edges += (height - 1) * width * N; // vertical edges
		num_edges += (width - 1) * height * N; // horizontal edges
		num_edges += width * height * 2; // apex edges

		// problem line
		os << "p min " << num_nodes << " " << num_edges << '\n';

		auto node_id = [this](int i, int j) { return 1 + i * width + j; };
		auto apex_id = width * height + 1;

		// node descriptors
		for (auto i = 0; i < height; ++i) {
			for (auto j = 0; j < width; ++j) {
				const auto &nodedata = grid[i][j];
				if (nodedata.deficiency != 0)
					os << "n " << node_id(i, j) << " " << nodedata.deficiency << '\n';
			}
		}

		if (apex.deficiency != 0)
			os << "n " << apex_id << " " << apex.deficiency << '\n';

		// arc descriptors
		for (auto i = 0; i < height; ++i) {
			for (auto j = 0; j < width; ++j) {
				const auto &nodedata = grid[i][j];

				// down
				if (i < height - 1) {
					const auto &arcdata = nodedata.down;
					assert(arcdata.flow % 2 == 0 && "flow must be an even number after the initial circulation");
					const auto r = ReturnR(arcdata.flow);
					for (auto a = 0; a < r; ++a)
						os << "a " << node_id(i + 1, j) << " " << node_id(i, j) << " 0 2 " << -arcdata.difference[a] << '\n';
					for (auto a = r; a < N; ++a)
						os << "a " << node_id(i, j) << " " << node_id(i + 1, j) << " 0 2 " << arcdata.difference[a] << '\n';
				}

				// right
				if (j < width - 1) {
					const auto &arcdata = nodedata.right;
					assert(arcdata.flow % 2 == 0 && "flow must be an even number after the initial circulation");
					const auto r = ReturnR(arcdata.flow);
					for (auto a = 0; a < r; ++a)
						os << "a " << node_id(i, j + 1) << " " << node_id(i, j) << " 0 2 " << -arcdata.difference[a] << '\n';
					for (auto a = r; a < N; ++a)
						os << "a " << node_id(i, j) << " " << node_id(i, j + 1) << " 0 2 " << arcdata.difference[a] << '\n';
				}

				// apex
				os << "a " << node_id(i, j) << " " << apex_id << " 0 " << INFINITE << " " << nodedata.apexcost << '\n';
				os << "a " << apex_id << " " << node_id(i, j) << " 0 " << INFINITE << " 0\n";
			}
		}
		os.flush();
	}

	
	/**
	 * \brief Iterates over all edges that cross the boundary between the 2 parts and modifies the flow such that forward_length and backward_lengthe become non-negative. 
	 * \param half_width: the central point of the grid based on the width
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the boundary balance takes place
	 */
  	void VerticalLengthBalance(int half_width, int top, int bottom,  int left_bound, int right_bound)
  	{
  		const int temp_height= bottom - top;
  		for(int i=0; i<temp_height; i++ )
  		{
  			
  			int r = ReturnR(grid[top + i][half_width].right.flow);
  			grid[top + i][half_width].right.forward_length = grid[top + i][half_width].right.difference[r] +  grid[top + i][half_width].potential - grid[top + i][half_width+1].potential;
  			while( (grid[top + i][half_width].right.forward_length<0 && r<N) )
  			{
  				if(grid[top + i][half_width].right.forward_length<0 && r<N)
  				{
  					grid[top + i][half_width].right.flow++;
					
  				}
				
  				r = ReturnR(grid[top + i][half_width].right.flow);
  				grid[top + i][half_width].right.forward_length = grid[top + i][half_width].right.difference[r] +  grid[top + i][half_width].potential - grid[top + i][half_width+1].potential;
  				if(r==N)
  				{
  					break;
  				}
  				
  				
  			}
  			
  			int l = ReturnL(grid[top + i][half_width].right.flow);
  			grid[top + i][half_width].right.backward_length = -grid[top + i][half_width].right.difference[N-1-l] -  grid[top + i][half_width].potential + grid[top + i][half_width+1].potential;
  			while((grid[top + i][half_width].right.backward_length<0 && l<N))
  			{
  				if (grid[top + i][half_width].right.backward_length<0 && l<N)
  				{
  					grid[top + i][half_width].right.flow--;
  				}
  				
  				l = ReturnL(grid[top + i][half_width].right.flow);
  				grid[top + i][half_width].right.backward_length = -grid[top + i][half_width].right.difference[N-1-l] -  grid[top + i][half_width].potential + grid[top + i][half_width+1].potential;
  				if(l==N)
  				{
  					break;
  				}
  				
  			}
  			
  			assert(grid[top + i][half_width].right.flow>=-N);
  			assert(grid[top + i][half_width].right.flow<=N);
  		}
  		return;
  	}
  	
	
	/**
	 * \brief Iterates over all edges that cross the boundary between the 2 parts and modifies the flow such that forward_length and backward_lengthe become non-negative. 
	 * \param half_height: the central point of the grid based on the height
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the boundary balance takes place
	 */
  	void HorizontalLengthBalance(int half_height, int top, int bottom,  int left_bound, int right_bound)
  	{
  		int temp_width= right_bound - left_bound;
  		for(int i=0; i<temp_width; i++ )
  		{
  			int r = ReturnR(grid[half_height][left_bound + i].down.flow);
  			grid[half_height][left_bound + i].down.forward_length = grid[half_height][left_bound + i].down.difference[r] +  grid[half_height][left_bound + i].potential - grid[half_height+1][left_bound + i].potential;
  			while((grid[half_height][left_bound + i].down.forward_length<0 && r<N) )
  			{
  				if(grid[half_height][left_bound + i].down.forward_length<0 && r<N)
  				{
  					grid[half_height][left_bound + i].down.flow++;
  				}
				
  				r = ReturnR(grid[half_height][left_bound + i].down.flow);
  				grid[half_height][left_bound + i].down.forward_length = grid[half_height][left_bound + i].down.difference[r] +  grid[half_height][left_bound + i].potential - grid[half_height+1][left_bound + i].potential;
  				if(r==N)
  				{
  					break;
  				}
  				
  			}
  			
  			int l = ReturnL(grid[half_height][left_bound + i].down.flow);
  			grid[half_height][left_bound + i].down.backward_length = -grid[half_height][left_bound + i].down.difference[N-1-l] - grid[half_height][left_bound + i].potential + grid[half_height+1][left_bound + i].potential;
  			while((grid[half_height][left_bound + i].down.backward_length<0 && l<N))
  			{
  				if(grid[half_height][left_bound + i].down.backward_length<0 && l<N)
  				{
  					grid[half_height][left_bound + i].down.flow--;
  				}
  				
  				l = ReturnL(grid[half_height][left_bound + i].down.flow);
  				grid[half_height][left_bound + i].down.backward_length = -grid[half_height][left_bound + i].down.difference[N-1-l] - grid[half_height][left_bound + i].potential + grid[half_height+1][left_bound + i].potential;
  				if(l==N)
  				{
  					break;
  				}
  				
				
  			}
  			
  			assert(grid[half_height][left_bound + i].down.flow>=-N);
  			assert(grid[half_height][left_bound + i].down.flow<=N);
  		}
  		return;
  	}
	
	/**
	 * \brief Subtract from each grid node the apex potential
	 * \param apex_potential: current apex potential of the relative grid
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the update of potentials takes place
	 */
  	void UpdatePotentialsWithMax(int apex_potential, int top, int bottom,  int left_bound, int right_bound)
  	{
  		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
				grid[i][j].potential = grid[i][j].potential - apex_potential;
			}
		}
		
  		return;
  	}
  	
  	/**
	 * \brief Calculates the total deficiency of the subgrid
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the calculation of total deficiency takes place
	 * \Return total deficiency
	 */
  	int TotalDeficiency(const int top, const int bottom, const int left_bound, const int right_bound)
  	{
  		int total_deficiency = 0;
  		
  		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
				total_deficiency += grid[i][j].deficiency;
			}
		}
		
		return total_deficiency;
	}	
	
	/**
	 * \brief Checks the forward and bacward lengths for non-negativity
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the non-negativity check takes place
	 * \Return true or false for non-negativity, negatitity respectively.
	 */
	bool check(const int top, const int bottom, const int left_bound, const int right_bound)
	{
		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
				
				int r = ReturnR(grid[i][j].right.flow);
				int l = ReturnL(grid[i][j].right.flow);
				if(j+1<right_bound )
				{	
					if(r<N)
					{
						grid[i][j].right.forward_length		= grid[i][j].right.difference[r] +  grid[i][j].potential - grid[i][j+1].potential;
					}
					else
					{
						grid[i][j].right.forward_length = 0;
					}
					
					if(l<N)
					{
						grid[i][j].right.backward_length 	= -grid[i][j].right.difference[N-1-l] -  grid[i][j].potential + grid[i][j+1].potential;
					}
					else
					{
						grid[i][j].right.backward_length=0;
					}
				}
				else
				{
					grid[i][j].right.forward_length=0;
					grid[i][j].right.backward_length=0;
				}
				
				
				r = ReturnR(grid[i][j].down.flow);
				l = ReturnL(grid[i][j].down.flow);
				if(i+1<bottom)
				{
					if(r<N)
					{
						grid[i][j].down.forward_length		= grid[i][j].down.difference[r] +  grid[i][j].potential - grid[i+1][j].potential;
					}
					else
					{
						grid[i][j].down.forward_length = 0;
					}
					
					if(l<N)
					{
						grid[i][j].down.backward_length 	= -grid[i][j].down.difference[N-1-l] -  grid[i][j].potential + grid[i+1][j].potential;
					}
					else
					{
						grid[i][j].down.backward_length=0;
					}
				}
				else
				{
					grid[i][j].down.forward_length =0;
					grid[i][j].down.backward_length=0;
				}
				
				
				
				
				
				int enable = 0;
				if(	grid[i][j].right.forward_length<0)
				{
					cout << "--> Negative Length RF<--"<<" ("<<i<<","<<j<<")"<<endl;
					enable =1;
				}
				
				if(grid[i][j].right.backward_length<0) {
					cout << "--> Negative Length RB<--"<<" ("<<i<<","<<j<<")"<<endl;
					enable =1;
				}
				
				if(grid[i][j].down.forward_length<0)
				{
					cout << "--> Negative Length DF<--"<<" ("<<i<<","<<j<<")"<<endl;
					enable =1;
				}	
				
				if(grid[i][j].down.backward_length<0)
				{
					cout << "--> Negative Length DB<--"<<" ("<<i<<","<<j<<")"<<endl;
					enable =1;
				}
				const int h = apex.potential - grid[i][j].potential;
				if( h < 0 || h > grid[i][j].apexcost || (grid[i][j].flowtoapex < 0 && h != 0 ) ) {
					cout << "apex potential error @" << i <<"," << j << ": "<< h << " flow: " << grid[i][j].flowtoapex << endl;
				}
				if(enable)
				{	
					return false;
				}
			}
		}
		return true;
	}
	
	
	/**
	 * \brief Checks the deficiencies to be zero after the min cost flow solution
	 * \param bottom,top,left_bound,right_bound : Coordinates of grid in which the deficiency check takes place
	 * \Return true or false for zero, non-zero respectively.
	 */
	bool CheckDeficiencies(const int top, const int bottom, const int left_bound, const int right_bound)
	{
		int total_deficiency=0;
		for (int i=top; i<bottom; i++) 
		{
			for (int j=left_bound; j<right_bound; j++) 
			{
				total_deficiency=0;
				//Upward
				if(i-1>-1)
				{
					total_deficiency+= grid[i-1][j].down.flow;
				}
				//Downward
				if(i+1<height)
				{
					total_deficiency-= grid[i][j].down.flow;
				}
				//Leftward
				if(j-1>-1)
				{
					total_deficiency+= grid[i][j-1].right.flow;
				}
				//Rightwarad
				if(j+1<width)
				{
					total_deficiency-= grid[i][j].right.flow;		
				}
				
				total_deficiency-=grid[i][j].flowtoapex;
				
				if(total_deficiency>0  || total_deficiency<0)
				{
					return false;
				}
			}
		}
		
		assert(!apex.deficiency);
		
		return true;
	}
	
	
	void PrintFinalGrid(void)
	{
		int i,j;
		for ( i=0; i<height; i++)
		{
			for ( j=0; j<width; j++)
			{
				cout<<"("<<i<<","<<j<<")"<< "[Flow:"<< grid[i][j].right.flow << "," <<grid[i][j].down.flow <<"][Pot:"<<  grid[i][j].potential<< "][Dis:" << grid[i][j].distance << "][Par:"<<grid[i][j].parent.y<<","<<grid[i][j].parent.x<<"]\t";
			}
			cout << endl;
		}
		
		cout<< "Apex Pot:"<< apex.potential<<endl;
		
		return;
	}
};




/** \struct PPMImage
 * \brief Initializes the grid
 */
template< int N >
struct InstanceConversion {
	ApexGrid<N> apexgrid;
	
	InstanceConversion( const PPMImage& img ) : apexgrid( img ) {}

	bool toImage() {
	  return apexgrid.toImage();
  }
	    
};


#endif