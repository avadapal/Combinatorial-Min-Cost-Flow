#include<iostream>
#include<queue>
#include<vector>
#include<algorithm>
#include "graph.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include "simple_solver.h"
#include "random.h"
#include <unistd.h>


using namespace std;

int main()
{
  ofstream file_obj,file_obj1;
  file_obj.open("stretch_iterations.txt");
  struct stat st = {0};
  //cout<<"Enter the number of verticies in the graph:";
  unsigned int n; //n represents the number of verticies and m represents the number of edges
  cin>>n;
  cout<<std::endl<<"Enter the value of beta: "<<std::endl;
  float beta;
  cin>>beta;
  if(n==20)
   {
   if (stat("/home/avadapal/electrical_flow/data/graph_20", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_20", 0700);
    }
  
   std::string path = "/home/avadapal/electrical_flow/data/experiment_20.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl;

   }
  if(n==40)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_40", 0700);
    }
   std::string path = "/home/avadapal/electrical_flow/data/experiment_40.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl;	
   }
  if(n==160&&beta==120)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_160_b14", 0700);
    }
  
   std::string path = "/home/avadapal/electrical_flow/data/experiment_160_b14.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl;
   }

 if(n==160&&beta==240)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_160_b15", 0700);
    }
  
   std::string path = "/home/avadapal/electrical_flow/data/experiment_160_b15.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl;
   }

 if(n==160&&beta==480)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_160_b16", 0700);
    }
  
   std::string path = "/home/avadapal/electrical_flow/data/experiment_160_b16.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl;
   }

 if(n==160&&beta==960)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_160_b17", 0700);
    }
  
   std::string path = "/home/avadapal/electrical_flow/data/experiment_160_b17.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl;
   }

  if(n==80&&beta==0.5)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b1", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b1.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }


   if(n==80&&beta==0.75)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b2", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b2.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

  if(n==80&&beta==1.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b3", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b3.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

   if(n==80&&beta==1.5)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b4", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b4.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

   if(n==80&&beta==2.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b5", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b5.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

   if(n==80&&beta==3.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b6", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b6.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

   if(n==80&&beta==4.5)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b7", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b7.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }


   if(n==80&&beta==7.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b8", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b8.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

   if(n==80&&beta==10.5)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b9", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b9.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

  if(n==80&&beta==15.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b10", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b10.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }
  

   if(n==80&&beta==30.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b11", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b11.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }
  
  if(n==80&&beta==45.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b12", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b12.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

  if(n==80&&beta==60.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b13", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b13.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

 if(n==80&&beta==120.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b14", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b14.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

 if(n==80&&beta==240.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b15", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b15.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

 if(n==80&&beta==480.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b16", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b16.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

 if(n==80&&beta==960.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b17", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b17.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

 if(n==80&&beta==1920.0)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_80_b18", 0700);
    }

   std::string path = "/home/avadapal/electrical_flow/data/experiment_80_b18.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path); 
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
  }

  if(n==320)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_320", 0700);
    }
   std::string path = "/home/avadapal/electrical_flow/data/experiment_320.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path);
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
   }
 
 if(n==640)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_640", 0700);
    }
   std::string path = "/home/avadapal/electrical_flow/data/experiment_640.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path);
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
   }

if(n==1280)
   {
   if (stat("/home/avadapal/electrical_flow/data/folder", &st) == -1) {
    mkdir("/home/avadapal/electrical_flow/data/graph_1280", 0700);
    }
   std::string path = "/home/avadapal/electrical_flow/data/experiment_1280.txt";
   std::ofstream outfile (path);
   outfile.close();
   file_obj1.open(path);
   file_obj1<<"Stretch"<<"	"<<"TCN"<<"	"<<"Number of Iterations"<<"	"<<"Error"<<std::endl; 
   }


  Graph<double> G(n); //,G1(n);
  //G=G1;
  Random rg;
  //Graph<double> G1(n);
  G.create_graph();
  /*
  ofstream fo;
  fo.open("graph.txt");
  for(int i=1; i<=G.count; i++)
   {
    int t = G.tails[i];
    int h = G.heads[i];
    fo<<t<<"	"<<h<<"	"<<i<<std::endl;
    fo<<h<<"	"<<t<<"	"<<-i<<std::endl;
   } */


  G.generate_cycles();
  G.generate_edge_cycle_correspondence();
  int parents_temp[n+1];
  vector<int> tree_edges_temp, non_tree_edges_temp; 
   vector<double> resistances_temp, voltages_temp;
 
  //int str = G.calculate_stretch(); 
  //simple_solver(G, file_obj, str);
  //G.print();
  //G.the_sanity_check();	  
  //cout<<"Stretch: "<<x<<std::endl;

    int count = 0;
    int acceptance  = 0;
    int rejection = 0;
    while(count <1000)
   {  
    float tree_condition_number;    
  
    float stretch = G.calculate_stretch();
    for(unsigned int i=1;i<=n;i++)
     {
      parents_temp[i]=G.parents[i];
     } 

    // parents_temp = G.parents;
     resistances_temp = G.resistances;
     voltages_temp = G.voltages;
     tree_edges_temp = G.tree_edges;
     non_tree_edges_temp = G.non_tree_edges;
    
    tree_condition_number = stretch+G.count -2*G.no_of_verticies + 2;
 
  //  double number = run_experiments(file_obj, stretch, file_obj1, G, beta);
    count++;

    
   
    double e1, e2;

    double k;
    k = stretch*log(stretch*tree_condition_number/0.1);

    e1 = k/number; 
  
    G.change_spanning_tree();

    stretch = G.calculate_stretch();    

    tree_condition_number = stretch+G.count -2*G.no_of_verticies + 2;

    number = run_experiments(file_obj, stretch, file_obj1, G, beta);
    count++;

    k = stretch*log(stretch*tree_condition_number/0.1);
    e2 = k/number;
  
    double delta_e = e2-e1;
    
    if((delta_e) <= 0){
    //Do nothing
    }
    else{
    //float beta = 3.0;
    float prob = exp(-1*beta*delta_e);

    double r = rg.rng();
    if(r>prob){
 
     rejection++;

    for(unsigned int i=1;i<=n;i++)
     {
      G.parents[i] = parents_temp[i];
     }
 
     G.voltages = voltages_temp;
     G.resistances =resistances_temp; 
     G.tree_edges = tree_edges_temp;
     G.non_tree_edges = non_tree_edges_temp;
      }
    else
     {
      acceptance++;    
     }
    }     

  }   

int total = acceptance + rejection;
float rate = acceptance/total;
cout<<std::endl<<std::endl<<"Acceptance Rate: "<<acceptance<<"	"<<rejection<<"	"<<rate<<std::endl;
file_obj1.close();
      

  return 0;
  
}		





