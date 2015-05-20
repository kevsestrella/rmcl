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
#include <math.h>
using namespace std;

#define BUFFSIZE 4096 //ala doc mana! haha
typedef map < string, set<string> >::iterator itermap;

int main(int argc, char **argv){

  if(argc != 8) {
    printf("Usage: <input CAw output file> <input id-protein mapping file> <input cluster-protein mapping file> <output file> <threshold> <alpha> <gamma>\n");
    return 0;
  }

  FILE *inputclusters = fopen(argv[1], "r");
  FILE *proteinmap = fopen(argv[2], "r");
  FILE *clustermap = fopen(argv[3], "r");
  FILE *outputfile = fopen(argv[4], "a");
  float threshold = atof(argv[5]);
  float alpha = atof(argv[6]);
  float gamma = atof(argv[7]);
  if(inputclusters == NULL || proteinmap == NULL || clustermap == NULL || outputfile == NULL) {
      perror("Unable to open file!");
      return 0;
  }

	map<string, set<string> > complex_map; //map of complexes from literature. key is the name of the complex
  map<int, string> protein_map; //map of the protein ids to their names
  set<string> protein_members; //set of proteins currently being considered from file

  string temp_name;
  size_t len;
  char line[BUFFSIZE];
  int pos, postart;
  int sumclusters = 0;
  int flag = 0;
  int protein_id = 1;
  int clustcount = 1;
  double cyc_size, clust_size, intersect_size;
  
  double fscore = 0.0;
  double tempscore = 0.0;
  double avgFscore = 0.0;
  
  double precision = 0.0;
  double tempprecis = 0.0;
  double avgPrecis = 0.0;

  double recall = 0.0;
  double temprec = 0.0;
  double avgRec = 0.0;


  int intersect2 = 0;
//READ CLUSTER MAPPING FILE AND STORE IN MAP
  //fprintf(outputfile, "%s %s %s %s %s %s %s %s\n", "threshold", "alpha", "gamma", "avgfscore", "avgprecision", "avgrecall", "geommean", "numclusts");
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

      if(intersect_size != 0){
        tempscore = (2 * (intersect_size/cyc_size) * (intersect_size/clust_size) ) / 
        ( (intersect_size/cyc_size) + (intersect_size/clust_size) );

        tempprecis = (intersect_size/clust_size);
        temprec = (intersect_size/cyc_size);
      }
      else{
        tempscore = 0.0;
        tempprecis = 0.0;
        temprec = 0.0;
      }  

      if(tempscore > fscore){
        fscore = tempscore;
        precision = tempprecis;
        recall = temprec;

        temp_name = it->first;
        flag = 1;
        
      }
    
    }
    
    if(flag == 0){//(fscore == 0.0){
      temp_name = "NAN!";
    }
    else{ //comment out the "else" if you want to include NAN clusters
      sumclusters = sumclusters + clust_size;
    }
  //END COMPUTE FSCORE, OUTPUT TO FILE

  //fprintf(outputfile, "CLUSTER %d: %f %f %f\n", clustcount, precision, recall, fscore);
  //fprintf(outputfile, "Precision: %f\n", precision);
  //fprintf(outputfile, "Recall: %f\n\n", recall);

  avgFscore = avgFscore + (fscore * clust_size);
  avgPrecis = avgPrecis + (precision * clust_size);
  avgRec = avgRec + (recall * clust_size);
  //reset to defaults
  fscore = 0.0;
  precision = 0.0;
  recall = 0.0;
  temp_name = "NAN";
  flag = 0;
  protein_members.clear();
  clustcount++;
  }

  double finalAvgFScore = (avgFscore/sumclusters);
  double finalAvgPrecision = (avgPrecis/sumclusters);
  double finalAvgRecall = (avgRec/sumclusters); 
  double geometricmean = pow((finalAvgFScore * finalAvgRecall * finalAvgPrecision), 1.0/3.0);

  //(finalAvgFScore * finalAvgRecall * finalAvgPrecision)^(1/3); 


  fprintf(outputfile, "%f %f %f %f %f %f %f %d\n", threshold, alpha, gamma, finalAvgFScore, finalAvgPrecision, finalAvgRecall, geometricmean, (clustcount-1));
  
  //fprintf(outputfile, "Average Fscore: %f\n", finalAvgFScore);
  //fprintf(outputfile, "Average Precision: %f\n", finalAvgPrecision);
  //fprintf(outputfile, "Average Recall: %f\n", finalAvgRecall);

  printf("%s\n", "Finished computing F-scores.");

  fclose(inputclusters);
  fclose(proteinmap);
  fclose(clustermap);
  fclose(outputfile);
}
   
