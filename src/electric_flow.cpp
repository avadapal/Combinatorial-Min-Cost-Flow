#include<iostream>
#include<queue>
#include<vector>
#include<algorithm>

using namespace std;

class Graph

{

	public:
	
	bool visited


}



int main()

{

/* Taking the matrix B as input */
cout<<"Please enter the number of Edges:";
int no_of_edges;
cin>>no_of_edges;

cout<<"Please enter the number of verticies:";
int no_of_verticies;
cin>>no_of_verticies;

int B[no_of_edges][no_of_verticies];

cout<<"Please enter the values in the matrix B:";

int i,j; /* For loop counters */

for(i=0; i<no_of_edges; i++)
	{
	for(j=0; j<no_of_verticies; j++)
		{
		cin>>B[i][j];
		}
	}
/* The matrix B is taken as the input */



/* Taking the vector R (Resistances) as input */
int R[no_of_edges];
cout<<"\n\n";
cout<<"Please enter the Resistance Vector";
for(int i=0; i<no_of_edges; i++)
	{
	cin>>R[i];
	} 
/* The vector R (Resistances) is taken as input */


/* Taking the vector V (Volatages) as input */
int V[no_of_edges];
cout<<"\n\n";
cout<<"Please enter the Voltage Vector";
for(i=0; i<no_of_edges; i++)
	{
	cin>>V[i];
	}
/* Vector V (Voltages) is taken as input */






/* Getting the spanning tree T */

int edges_spanning_tree;
edges_spanning_tree = no_of_verticies -1;

int T[edges_spanning_tree][no_of_verticies];
//int k;
int visited[no_of_verticies];
int visited_counter;

/* Initializing the array - Visited */ 

	for(visited_counter = 0; visited_counter < no_of_verticies; visited_counter++)
			{
				visited[visited_counter]=0;
			}
 /*Initialization of the array - Visited Done! */ 


visited_counter = 0;
int k,l;
int flag;
int t = 0;
int visit;


/**************Sanity Check Starts *****************/
cout<<"\n\n";

for(i=0; i<no_of_edges; i++)
	{
		for(j=0; j<no_of_verticies; j++)
			{
				cout<<B[i][j];
				
			}
		cout<<"\n\n";
	}

cout<<"\n\n\n";
/****************Sanity Check Ends *******************/

int number_non_tree_edges;
number_non_tree_edges = no_of_edges - edges_spanning_tree;

int edges_in_spanning_tree[edges_spanning_tree];
int edges_spanning_tree_count = 0;
int edges_not_in_spanning_tree[number_non_tree_edges];
int number_non_tree_edges_count = 0;

for(j=0; j<no_of_verticies; j++)
	{
		
 		for(i=0; i<no_of_edges; i++)
			{
				if(B[i][j] == 1)
	   
					{
						for(k=0; k<no_of_verticies; k++)
							{
								if(B[i][k]== -1)
								 	{
										visit = k+1;
										
										flag = 0;

										for(l=0; l<no_of_verticies; l++)
											{
												if(visited[l]== visit)
													{
														flag = 1;
													}
											 }
									}
							}

					
	
							if(flag==0)
								{

									for(l=0; l< no_of_verticies; l++)
										{
									T[t][l]=B[i][l];
										}
									
									edges_in_spanning_tree[edges_spanning_tree_count] = i;
									edges_spanning_tree_count++;

									visited[visited_counter] = visit;
									visited_counter++;
									t++;

				      				}

					}

			}
	}	
	
									
																									
 /* The Spanning Tree is obtained */						
																						


/* Finding all the non-tree edges*/

for(i=0; i<no_of_edges; i++)
	{

		flag = 0;

		for(j=0; j<edges_spanning_tree; j++)
			{			
			if(i==edges_in_spanning_tree[j])
				{	
				flag=1;
				}
		
			}
		
		if(flag==0)
			{
			
			edges_not_in_spanning_tree[number_non_tree_edges_count]=i;
			number_non_tree_edges_count++;
			}
	}


/* All the non-tree edges found */




/* Setting the initial feasible flow to an all-zero vector */
int initial_feasible_flow[no_of_edges];

for(i=0; i<no_of_edges; i++)
	{
	initial_feasible_flow[i]=0;
	}
/* We have set the initial feasible flow to be the all-zero vector */ 



/* Initializing the current flow to the initial feasible flow */

int current_flow[no_of_edges];
for(i=0; i<no_of_edges; i++)
	{
	current_flow[i] = initial_feasible_flow[i];
	}
/* Initializing the current flow to the initial feasible flwo done */




	 


/* Find the path between the Start-Vertex and End-Vertex in the Tree */
int paths[no_of_verticies][no_of_edges];
int current_vertex;
//current_vertex = start_vertex;



//std::queue<int> vertex[no_of_verticies];

vector<int> vertex[no_of_verticies];  // An array of vectors which contains all the incident edges of the graph

/* Initializing a two dimensional queue which contains all the edges incident on each vertex on the Graph */

for(j=0; j<no_of_verticies; j++)
	{
		
		for(i=0; i<no_of_edges; i++)
			{
			if(B[i][j]==1||B[i][j]==-1)
				{
			vertex[j].push_back(i);
				}
			}
		cout<<endl;	
	}


vector<int> vertex_tree[no_of_verticies]; //An array which contains all the incident edges of the spanning tree



for(i=0; i<no_of_verticies; i++)
{
	for(j=0; j<vertex[i].size(); j++)
		{
	

			flag = 0;
			for(k=0; k<number_non_tree_edges_count; k++)
				{	
					
					if(vertex[i][j]==edges_not_in_spanning_tree[k])
						{
							flag=1;
						}
					
				}
			if(flag==0)
						{
							vertex_tree[i].push_back(vertex[i][j]);
						}
		}
}
						
				
				




/* Sanity Check */

for(i=0; i<no_of_verticies; i++)
	{
		cout<<"Vertex "; cout<<i; cout<<":";	
		for(j=0; j<vertex[i].size(); j++)
			{
			cout<<vertex[i][j];
			cout<<" ;";
			}
		cout<<endl;
	}


cout<<endl;
cout<<"Another Sanity Check";
cout<<endl; cout<<endl;
for(i=0; i<no_of_verticies; i++)
	{
		cout<<"Vertex "; cout<<i; cout<<":";	
		for(j=0; j<vertex_tree[i].size(); j++)
			{
			cout<<vertex_tree[i][j];
			cout<<" ;";
			}
		cout<<endl;
	}



/* Initializing a two dimensional array which contains the end-points of each edge */

int edges[no_of_edges][2];
int count;



for(i=0; i<no_of_edges; i++)
	{
		count = 0;
		for(j=0; j<no_of_verticies; j++)
			{	
			if(B[i][j]==1||B[i][j]==-1)
				{
					edges[i][count]=j;
					count++;
				}
			}
	}


cout<<"\n\n";
cout<<"The Edge Matrix:";
cout<<"\n\n";
for(i=0; i<no_of_edges; i++)
	{
		for(j=0;j<2;j++)
			{
			cout<<edges[i][j];
			}
	}

 std::queue<int>path;
 while(!path.empty())
	{
	path.pop();
	}

path.push(current_vertex);

int current_edge;
int next_vertex;


/* Spanning Tree Work */
// We randomly pick a root, say the 0th vertex


int parent[no_of_verticies];
parent[0] = -1;  
std::queue<int> bfs;
int incident_edge;
int child_vertex;
int root;
bfs.push(0);
while(!bfs.empty())
{
root = bfs.front();
bfs.pop();
for(i=0; i< vertex_tree[root].size(); i++)
	{

	incident_edge = vertex_tree[root][i];
	cout<<"Incident Edge: "; cout<<incident_edge; cout<<endl;
	cout<<"root:"; cout<<root;
	if(edges[incident_edge][0]==root&&edges[incident_edge][1]!=parent[root])
		{
		child_vertex = edges[incident_edge][1];
		bfs.push(child_vertex);
		cout<<"root: "; cout<<root;
		parent[child_vertex] = root;		
		cout<<"Child Vertex: "; cout<<child_vertex;
				
		}
	
	if(edges[incident_edge][1]==root&&edges[incident_edge][0]!=parent[root])
		{
		
		child_vertex = edges[incident_edge][0];
		bfs.push(child_vertex);
		cout<<"root: "; cout<<root;
		parent[child_vertex]= root;
		cout<<"Child Vertex: "; cout<<child_vertex;
		
		}	

	}

}
cout<<endl; cout<<endl;
cout<<"The Big moment in my code!";
cout<<endl; cout<<endl;
for(i=0; i<no_of_verticies; i++)
	{
	cout<<"Parent of "; cout<<i; cout<<" :"; cout<<parent[i];
	cout<<endl;	
	}




/* Finding the cycle which is formed with the i-th non-tree edge and the spanning tree */


	int edge;
	int start_vertex;
	int end_vertex;
	int count_spanning_tree = 0;
	vector<int> edges_in_cycle[number_no_tree_edges_count];

	while(count_spanning_tree<number_non_tree_edges_count)
	{
	edge = edges_not_in_spanning_tree[count_spanning_tree]; //This is the edge, which along with the spanning tree forms the cycle
	for(i=0; i<no_of_verticies; i++)
		{
		if(B[edge][i]== -1)
			{
			start_vertex = i;
			}
		
		if(B[edge][i] == 1)
			{
			end_vertex = i;
			cout<<"End Vertex:"; cout<<end_vertex;
			}
		
		}


	cout<<"\n\n";
	cout<<"Edge:"; cout<<edge;
	cout<<"Start Vertex:"; cout<<start_vertex;
	cout<<"End Vertex:";   cout<<end_vertex;

	vector<int> verticies_in_path_1;
	vector<int> verticies_in_path_2;
	vector<int> verticies_in_path;
	int v,u;
	v=start_vertex;
	u=end_vertex;
	verticies_in_path_1.push_back(v);
	verticies_in_path_2.push_back(u);
	for(;;)
	{

	if(parent[v]!=-1&&parent[v]!=u)
		{
		v=parent[v];
		verticies_in_path_1.push_back(v);
		}

	
	if(parent[v]==-1)
		{
			
			break;
		}
	if(parent[v]==u)
		{
			break;
		}
	
	}

	
	for(;;)
	{

	if(parent[u]!=-1&&parent[u]!=v)
		{
		u=parent[u];
		cout<<"U =";
		cout<<u;
		verticies_in_path_2.push_back(u);
		}
	if(parent[u]== -1)
		{
			cout<<"Enters here with U = :"; cout<<u;
			break;
		}
	if(parent[u]==v)
		{
		cout<<"Enters here with U = :"; cout<<u;
			
			break;
		}
	u=parent[u];
	verticies_in_path_2.push_back(u);
	}
	
	verticies_in_path.reserve(verticies_in_path_1.size() + verticies_in_path_2.size());
	reverse(verticies_in_path_2.begin(),verticies_in_path_2.end());
	verticies_in_path.insert(verticies_in_path.end(), verticies_in_path_1.begin(), verticies_in_path_1.end());
	verticies_in_path.insert(verticies_in_path.end(), verticies_in_path_2.begin(), verticies_in_path_2.end());
	cout<<endl; cout<<endl;
	cout<<"The Paths:";
	
	for(i=0; i<verticies_in_path.size(); i++)
		{
			
			cout<<verticies_in_path[i];
		}
	cout<<endl<<endl;	
	
/* Getting all the edges in the cycle */

int current_edges;
for(i=0; i<verticies_in_path.size() ; i++)
	{
	
		for(j=0; j< vertex_tree[i].size() ; j++)
			{
			
			current_edges= vertex_tree[i][j];
			if(edges[current_edges][0]==verticies_in_path[i]&&edges[current_edges][1]==verticies_in_path[i+1])
				{
				edges_in_cycle[l].push_back(current_edges);
				}
			if(edges[current_edges][1]==verticies_in_path[i]&&edges[current_edges][0]==verticies_in_path[i+1])
				{
				edges_in_cycle[l].push_back(current_edges);
				}		
	
			}

	}
	
edges_in_cycle.push_back(edge);

cout<<"Edges in the cycle";
cout<<endl<<endl;



count_spanning_tree++;

}



for(i=0; i< edges_in_cycle.size() ; i++)
	{

	cout<<edges_in_cycle[i];	

	}
cout<<endl<<endl;



/* The main iterations of the Simple_Solver */






/********************************************Sanity Check Starts*******************************/
cout<<"Matrix B";
cout<<"\n";
cout<<"====================================================";
cout<<"\n\n";

for(i=0; i<no_of_edges; i++)
	{
		for(j=0; j<no_of_verticies; j++)
			{
				cout<<B[i][j];
				
			}
		cout<<"\n\n";
	}

cout<<"\n\n\n"; 


cout<<"Matrix T";
cout<<"\n";
cout<<"====================================================";
cout<<"\n\n";

for(i=0; i<edges_spanning_tree; i++)
	{
		for(j=0; j<no_of_verticies; j++)
			{
				cout<<T[i][j];
				
			}
		cout<<"\n\n";
	}

cout<<"\n\n\n";

cout<<"The Edges which form the spanning tree";
cout<<"\n\n";
cout<<"==============================================================";
cout<<"\n\n";
for(i=0; i<edges_spanning_tree; i++)
	{
		cout<<edges_in_spanning_tree[i];
		cout<<";";		
	}


cout<<"\n\n\n";

cout<<"The Edges which form the  non-tree edges";
cout<<"\n\n";
cout<<"==============================================================";
cout<<"\n\n";
for(i=0; i<number_non_tree_edges_count; i++)
	{
		cout<<edges_not_in_spanning_tree[i];
		cout<<";";		
	}







/********************************************Sanity Check Ends************************************/



return 0;


}

