
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
using namespace std;
typedef map < string, set<string> >::iterator itermap;
typedef map < int, pair<string, double> >::iterator iterout;
typedef set<string>::iterator iterset;
typedef vector<string>::iterator itervec;

int main(){
	ifstream p_complexfile("testfile.txt");  //suppl2.xls
	ifstream clusteroutput("testclusteroutput");
	ofstream p_complexlist("outlist.txt");

	map<string, set<string> > complex_map;
	map<int, pair<string, double> > output_clusters;
	set<string> protein_members; 

	string line, protein, complex_name, buff, temp_name;
	size_t n, m, p, end; 
	int clust_num;
	double  cyc_size, clust_size, intersect_size, fscore, tempscore; 
	
	if (p_complexfile.is_open()){
		while(getline(p_complexfile,line)){	 	
		 	n = line.find("  ");
		 	
		 	if(n != string::npos){
		 		protein = line.substr(0, n);
		 		p = protein.find(" ");
		 		protein = protein.substr(p+1);
		 	}	

		 	m = line.find_first_not_of(" ", n+1);
		 	complex_name = line.substr(m);
		 	
		 	end = line.find("  ", m+1);
		 	if(end != string::npos){
		 		complex_name = line.substr(m, (end-m));
		 		
		 	}
		 	
		 	//clean complex_name for possible GO ID's inside it!
		 	
		 	

		 	complex_map[complex_name].insert(protein);
		
    	}
    p_complexfile.close();
    }

    else
    	cout << "Unable to open complex library file";
    
    if(clusteroutput.is_open()){
    	stringstream ss;
    	clust_num = 1;
    	while(getline(clusteroutput,line)){
    		//tokenize line and add to set.
    		//cout << "THE OUTPUT LINE IS: " << line << "\n";
    		ss.str(line);
    		while(ss >> buff){
    			protein_members.insert(buff);
    		}

    		//get max Fscore
    		fscore = 0.0;
    		for (itermap it = complex_map.begin(); it != complex_map.end(); it++) {
  				vector<string> intersect;
  				set_intersection(protein_members.begin(), protein_members.end(),
  											it->second.begin(), it->second.end(), back_inserter(intersect) );
  				cyc_size = it->second.size();
  				clust_size = protein_members.size();
  				intersect_size = intersect.size();

  				tempscore = (2 * (intersect_size/cyc_size) * (intersect_size/clust_size) ) / 
  						( (intersect_size/cyc_size) + (intersect_size/clust_size) );

  				if(tempscore > fscore){
  					fscore = tempscore;
  					temp_name = it->first;
  				}
  			}

  			if(fscore == 0){
  				temp_name = "NAN!";
  			}	
   			
  			pair<string, double> fscore_entry = make_pair (temp_name,  fscore);
  			output_clusters[clust_num] = fscore_entry;

    		//clear string stream and set
    		protein_members.clear();
    		ss.str( std::string() );
			ss.clear();
    		clust_num++;
    	}
    
    	/*print the output clusters and their fscores */

    	for (iterout i  = output_clusters.begin(); i != output_clusters.end(); i++) {
  			cout << "CLUSTER " << i->first << ": " << i->second.first << ", " << i->second.second << "\n";
  		}


  	p_complexlist.close();

    }

    else
    	cout << "Unable to open with output file";
	
}