#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <cstring>
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
  int sumclusters = 0;
  int flag = 0;
  int protein_id = 1;
  double cyc_size, clust_size, intersect_size;
  double precision = 0.0;
  double recall = 0.0;
  double fscore = 0.0;
  double tempscore = 0.0;
  double avgPrecision = 0.0;
  double avgRecall = 0.0;
  double avgFscore = 0.0;
  size_t len;


  int intersect2 = 0;
//READ CLUSTER MAPPING FILE AND STORE IN MAP
  printf("%s\n", "Begin reading complex (literature) mapping...");
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
   
   int temp = 0;
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
      temp++;
      }
    }
  }
  printf("%s\n","Finished reading complex mapping file.");
  printf("%s\n", "Begin reading protein mapping file...");
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
printf("%s\n","Finished reading protein mapping file.");
//READ INPUT FILE AND CALCULATE FSCORE
printf("%s\n", "Now calculating F-scores of output clusters...");
  while(fgets(line, sizeof line, inputclusters) != NULL){
    char *tok = strtok(line, " ");
    while(tok != NULL){
      protein_members.insert(protein_map[atoi(tok)]);
      tok = strtok(NULL, " ");
    }

    //COMPUTE FSCORE
    clust_size = protein_members.size();
    
    for (itermap it = complex_map.begin(); it != complex_map.end(); it++){
      vector<string> intersect;
      set_intersection(protein_members.begin(), protein_members.end(),
                        it->second.begin(), it->second.end(), back_inserter(intersect));
      intersect_size = intersect.size();
      cyc_size = it->second.size();

      //bdhedmgfkflkf  -- WEALTHY WUZ HERE OHH YEAAAAHHH

      if(intersect_size != 0){
        tempscore = (2 * (intersect_size/cyc_size) * (intersect_size/clust_size) ) / 
        ( (intersect_size/cyc_size) + (intersect_size/clust_size) );
      }
      else{
        tempscore = 0.0;
      }  
      if(tempscore > fscore){
        precision = intersect_size/clust_size;
        recall = intersect_size/cyc_size;
        fscore = tempscore;
        temp_name = it->first;
        flag = 1;
        
      }
      intersect2 = 0;
    }
    
    if(flag == 0){//(fscore == 0.0){
      temp_name = "NAN!";
    }
    else{ //comment out the "else" if you want to include NAN clusters
      sumclusters = sumclusters + clust_size;
    }
  //END COMPUTE FSCORE, OUTPUT TO FILE

  fprintf(outputfile, "%s precision: %f recall: %f fscore: %f\n", temp_name.c_str() ,precision,recall, fscore);
  avgFscore = avgFscore + (fscore * clust_size);
  avgPrecision =  avgPrecision + (precision * clust_size);
  avgRecall = avgRecall + (recall * clust_size);
  //reset to defaults
  fscore = 0.0;
  temp_name = "NAN";
  flag = 0;
  protein_members.clear();
  }

   double finalAvgPrecision = avgPrecision/sumclusters;
   double finalAvgRecall = avgRecall/sumclusters;
  double finalAvgFScore = (avgFscore/sumclusters);
  fprintf(outputfile, "Clusters: %d, Average Precision: %f, Average Recall: %f,Average Fscore: %f\n", sumclusters, finalAvgPrecision, finalAvgRecall, finalAvgFScore);
  printf("%s\n", "Finished computing F-scores.");

  fclose(inputclusters);
  fclose(proteinmap);
  fclose(clustermap);
  fclose(outputfile);
}
   
