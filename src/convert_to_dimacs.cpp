#include<iostream>
#include<fstream>
#include<string>
#include<queue>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>

using namespace std;

int main(){


  std::fstream file;
  std::fstream file_dimacs;
  
  
  
  
  file.open ("Nodes.txt", std::fstream::in | std::fstream::out);
  file_dimacs.open("dimacs.txt", std::fstream::in | std::fstream::out);
  string id, price, rank, excess, line;  
  
  file_dimacs << "p" << " " << "min" << " " << 10001 << " " << 40400 << endl;
  int j =0;
  while(!file.eof())
{
    if(j > 0){
        file_dimacs << "n" << " ";
    }
     getline(file, id, ',');
     cout << id << " " ; 
     if(j > 0){
     file_dimacs << id;
     }
     getline(file, price, ',') ;
     cout << price;
     //file_dimacs << price << " ";
    
     getline(file, rank, ',') ;
     cout <<  rank << " "  ; 
     //file_dimacs << rank << " "; 
     
     getline(file, excess, '\n') ; 
     if(j > 0){
     file_dimacs << excess  << endl;
     }
     j++;
}
  std::fstream file1;
   file1.open ("Arcs.txt", std::fstream::in | std::fstream::out);
  
  string  head, tail, capacity, cost;  
  int i = 0;
  while(!file1.eof())
{
    if(i > 0){
    file_dimacs << "a ";
    }
     getline(file1, head, ',');
    if(i > 0){
     file_dimacs << head ;
    }
     getline(file1, tail , ',') ;
 
     if(i > 0){
     file_dimacs << tail << " ";
     }
     getline(file1, cost, ',') ;
     if(i > 0)
     {
     file_dimacs << 0 ;
     }
     if(i > 0){
     file_dimacs << cost;
     }
     getline(file1, capacity, '\n') ; 
     
     if(i > 0){
     file_dimacs << capacity  << endl;
     }
  i++;
    
}
file.close();
file1.close();
return 0;
}