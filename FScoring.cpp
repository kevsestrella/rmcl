
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <boost/regex.hpp>
#include <stdio.h>
using namespace std;

#define BUFFSIZE 4096 //ala doc mana! haha

typedef map < string, set<string> >::iterator itermap;
typedef map<int, string>::iterator itertemp;
typedef set<string>::iterator iterset;
typedef vector<string>::iterator itervec;

int main(int argc, char **argv){

  if(argc != 5) {
    printf("Usage: <input CAw output file> <input id-protein mapping file> <input cluster-protein mapping file> <output file>\n");
    return 0;
  }

  FILE *inputclusters = fopen(argv[1], "r");
  FILE *proteinmap = fopen(argv[2], "r");
  FILE *clustermap = fopen(argv[3], "r");
  FILE *outputfile = fopen(argv[4], "w+");

  if(inputclusters == NULL || proteinmap == NULL || clustermap == NULL || outputfile == NULL) {
      perror("Unable to open file!");
      return 0;
  }

	map<string, set<string> > complex_map; //map of complexes from literature. key is the name of the complex
	map<int, pair<string, double> > input_clusters; //map of output. key: cluster number, value: complex name and F-score (pair)
  map<int, string> protein_map; //map of the protein ids to their names
  set<string> protein_members; //set of proteins currently being considered from file

  string temp_name;
  char line[BUFFSIZE];
  int pos, postart;
  int protein_id = 1;
  double cyc_size, clust_size, intersect_size;
  double fscore = 0.0;
  double tempscore = 0.0;
  size_t len;

//READ CLUSTER MAPPING FILE AND STORE IN MAP
  while(fgets(line, sizeof line, clustermap) != NULL){
    
    len = strlen(line);
    if(line[len-1] == '\n'){
      line[len-1] = '\0'; 
    }

   string line_c (line);

   pos = line_c.find("   ");
   string complex_name = line_c.substr(0,pos);
   string protein_list = line_c.substr(pos+3);
   pos = 0;
   postart = 0;
   
   while(pos != (protein_list.size() -1)){
      pos = protein_list.find(" ", pos);
      if(pos != string::npos){
        string proteinim = protein_list.substr(postart, (pos-postart));
        if(proteinim.at(0) == ' '){
          proteinim = proteinim.substr(1);
        }
        complex_map[complex_name].insert(proteinim);
        if(pos != protein_list.size()-1){
          postart = pos;
          pos++;
        } 
      }
    }
  }
//END OF READING CLUSTER MAPPING FILE

//READ PROTEIN MAPPING AND STORE IN MAP
  while(fgets(line, sizeof line, proteinmap) != NULL){
    len = strlen(line);
    if(line[len-1] == '\n'){
      line[len-1] = '\0';
    }
    protein_map[protein_id] = string(line);
    protein_id++;
  }
// END OF READING PROTEIN MAPPING

//READ INPUT FILE AND CALCULATE FSCORE
  while(fgets(line, sizeof line, inputclusters) != NULL){
    char *tok = strtok(line, " ");
    while(tok != NULL){
      protein_members.insert(protein_map[atoi(tok)]);
      tok = strtok(NULL, " ");
    }

    //COMPUTE FSCORE
    int temp = 0;

    for (itermap it = complex_map.begin(); it != complex_map.end(); it++){
      
      vector<string> intersect;
      set_intersection(protein_members.begin(), protein_members.end(),
                        it->second.begin(), it->second.end(), back_inserter(intersect));
          
      cyc_size = it->second.size();
      clust_size = protein_members.size();
      intersect_size = intersect.size();

      if(intersect_size != 0){
        tempscore = (2 * (intersect_size/cyc_size) * (intersect_size/clust_size) ) / 
        ( (intersect_size/cyc_size) + (intersect_size/clust_size) );
      }
      else{
        tempscore = 0.0;
      }  
      if(tempscore > fscore){
        fscore = tempscore;
        temp_name = it->first;
        temp = 1;
      }
    }
    
    if(temp == 0){//(fscore == 0.0){
      temp_name = "NAN!";
    }
  //END COMPUTE FSCORE, OUTPUT TO FILE

  fprintf(outputfile, "%s   %f\n", temp_name.c_str() ,fscore);
  
  //reset to defaults
  fscore = 0.0;
  temp_name = "NAN";
  temp = 0;
  protein_members.clear();
  }

  fclose(inputclusters);
  fclose(proteinmap);
  fclose(clustermap);
  fclose(outputfile);
}
   