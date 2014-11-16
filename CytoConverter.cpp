#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <cstdio>
#include <cstring>
using namespace std;

#define BUFFSIZE 4096 //ala doc mana! haha

typedef map < int, pair<string, int> >::iterator itermap;

int main(int argc, char **argv){


  char line[BUFFSIZE];
  int protein_id = 1;
  int len = 0;
  int index;
  int linenum = 0; 
  map<int, pair<string, int> > protein_map; //map of the protein ids to their names
  string source, dest;

  if(argc != 5) {
    printf("Usage: <input file to convert> <protein mapping file> <graph file> <output file>\n");
    return 0;
  }

  FILE *inputclusters = fopen(argv[1], "r");
  FILE *proteinmap = fopen(argv[2], "r");
  FILE *graphfile = fopen(argv[3], "r");
  FILE *outputfile = fopen(argv[4], "w+");

  if(inputclusters == NULL || proteinmap == NULL ||outputfile == NULL) {
      perror("Unable to open file!");
      return 0;
  }

  printf("%s\n", "Begin reading protein mapping file...");
//END OF READING CLUSTER MAPPING FILE

//READ PROTEIN MAPPING AND STORE IN MAP
  while(fgets(line, sizeof line, proteinmap) != NULL){
    len = strlen(line);
    if(line[len-1] == '\n'){
      line[len-1] = '\0';
    }
    protein_map[protein_id] = make_pair(string(line), 0);
    protein_id++;
  }

  printf("Finished reading mapping. Now reading clustering...\n");
  while(fgets(line, sizeof line, inputclusters) != NULL){
    len = strlen(line);
    if(line[len-1] == '\n'){
      line[len-1] = '\0';
    }

    char *splitter = strtok(line, " ");
    while(splitter != NULL){
      index = atoi(splitter);
      itermap itf = protein_map.find(index);
      if( itf != protein_map.end()){//it's in the map!
        itf->second.second = 1; 
      }
      splitter = strtok(NULL, " ");
    }
  }

  itermap it = protein_map.begin();
  while(it != protein_map.end()){
      if(it->second.second == 0){
        protein_map.erase(it++);
      }
      else{
        it++;
      }
  } 
  
  //for(itermap it2 = protein_map.begin(); it2 != protein_map.end(); it2++){
    //cout << "we have2: " << it2->first << " " << it2->second.first << " " << it2->second.second << "\n";
  //}

  itermap sourcenodes;

  while(fgets(line, sizeof line, graphfile) != NULL){
    //cout << "AT LINENUM: " << linenum << "\n";

    sourcenodes = protein_map.find(linenum);
    if(sourcenodes != protein_map.end()){
      //cout << " YEP!! " << "\n";
      len = strlen(line);
      if(line[len-1] == '\n'){
        line[len-1] = '\0';
      }

      char *splitter2 = strtok(line, " ");

      while(splitter2 != NULL){
        //cout << "split2: " << splitter2 << "\n";
        index = atoi(splitter2);
        itermap itf = protein_map.find(index);
        if( itf != protein_map.end()){
          //cout << "BOOM AT " << itf->first << "\n";
          source = sourcenodes->second.first;
          dest = itf->second.first;
          //cout << linenum << ": " << source << " " << index << ": " << dest << "\n";
          fprintf(outputfile, "%s%s%s\n", source.c_str(), ";" , dest.c_str());    
        }

        splitter2 = strtok(NULL, " ");
        //splitter2 = strtok(NULL, " ");
      }
    }


    linenum++;
    //cout << "\n\n\n";
  }

  fclose(inputclusters);
  fclose(proteinmap);
  fclose(graphfile);
  fclose(outputfile);  

}