#include<cstdio>
#include<cstring>
#include<cmath>
#include<iostream>
#include<vector>
#include<map>
#include<algorithm>
using namespace std;
#define edge pair<int, double>
#define EPS 1e-8
#define LEN 4096

void print_help(char *program_name) {
	printf("Usage: %s [options] <FS graph> <MLR-MCL clustering> <output file>\n", program_name);
	printf("Options:\n");
	printf("	-a <float>	alpha parameter (default 1.0)\n");
	printf("	-g <float>	gamma parameter (default 0.75)\n");
}

int main(int argc, char **argv) {
	double alpha = 1.0, gamma = 0.75;
	
	if(argc < 4) {
		print_help(argv[0]);
		return 0;
	}
	int offset = 1;
	for(; argv[offset][0] == '-';) {
		if(strlen(argv[offset]) != 2) {
			print_help(argv[0]);
			return 0;
		}
		switch(argv[offset][1]) {
			case 'a':
				alpha = atof(argv[offset+1]);
				offset+=2;
				if(alpha == 0.0) {
					printf("Invalid value.\n");
					print_help(argv[0]);
					return 0;
				}
				break;
			case 'g':
				gamma = atof(argv[offset+1]);
				offset+=2;
				if(gamma == 0.0) {
					printf("Invalid value.\n");
					print_help(argv[0]);
					return 0;
				}
				break;
			default:
				print_help(argv[0]);
				return 0;
		}
	}
	if(offset+2 != argc-1) {
		print_help(argv[0]);
		return 0;
	}
	
	FILE* wgraph = fopen(argv[offset], "r");
	FILE* clus = fopen(argv[offset+1], "r");
	FILE* complexes = fopen(argv[offset+2], "w+");
	
	if(wgraph == NULL || clus == NULL) {
		perror("");
    	return 0;
	}
	
	int V, Vnew, E;
	int numclusters = 0;
	char str[LEN];
	
	fscanf(wgraph, "%d%d%d", &V,&Vnew,&E);
	fgets(str, LEN, wgraph);
	
	int cluster[Vnew];
	int mapping[Vnew];
	bool isCore[Vnew];
	map<int,int> backmap;
	vector<edge> G[Vnew];
	
	fill(isCore, isCore+Vnew, false);
	
	for(int c = 0; c < Vnew; c++) {
		fscanf(clus, "%d", &cluster[c]);
		numclusters = max(numclusters, cluster[c]);
		
		fgets(str, LEN, wgraph);
		char *splitter = strtok(str, " ");
		mapping[c] = atoi(splitter)-1;
		backmap[mapping[c]] = c;
		splitter = strtok(NULL, " ");
		while(splitter != NULL) {
			int node = atoi(splitter);
			splitter = strtok(NULL, " ");
			double weight = atof(splitter);
			G[c].push_back(edge(node-1,weight));
			splitter = strtok(NULL, " ");
		}
	}
	numclusters++;
	
	double p_winconn[Vnew];                //weighted in-connectivity of protein p with proteins in cluster C
	double p_woutconn[Vnew];               //weighted out-connectivity of protein p wrt cluster C
	double avg_winconn[numclusters];       //average weighted in-connectivity of all proteins in cluster C
	double cluster_size[numclusters];      //cluster size
	double p_winconnC[Vnew];               //weighted in-connectivity of protein p with prelim core proteins of cluster C
	double p_woutconnC[Vnew];              //weighted out-connectivity of protein p wrt to prelim core of cluster C
	double avg_winconnC[numclusters];      //average weighted in-connectivity of all proteins in cluster C wrt prelim core of C
	double cluster_inter[numclusters];     //total interaction score between core proteins of cluster C
	vector<int> members[numclusters];      //cores of complexes
	vector<int> mematts[numclusters]; 	   //attachments of complexes
	
	fill(p_winconn, p_winconn+Vnew, 0.0);
	fill(p_woutconn, p_woutconn+Vnew, 0.0);
	fill(avg_winconn, avg_winconn+numclusters, 0.0);
	fill(cluster_size, cluster_size+numclusters, 0.0);
	fill(p_winconnC, p_winconn+Vnew, 0.0);
	fill(p_woutconnC, p_woutconn+Vnew, 0.0);
	fill(avg_winconnC, avg_winconn+numclusters, 0.0);
	fill(cluster_inter, cluster_inter+numclusters, 0.0);
	
	//[BEGIN] Preliminary Core
	
	for(int c = 0; c < Vnew; c++) {
		cluster_size[cluster[c]] += 1.0;
		for(int d = 0; d < G[c].size(); d++) {
			int j = backmap[G[c][d].first];
			if(cluster[c] == cluster[j])
				p_winconn[c] += G[c][d].second;
			else
				p_woutconn[c] += G[c][d].second;
		}
		avg_winconn[cluster[c]] += p_winconn[c];
	}
	for(int c = 0; c < numclusters; c++)
		avg_winconn[c] /= (2.0*cluster_size[c]);
	
	for(int c = 0; c < Vnew; c++) { /*p_winconn[c] - avg_winconn[cluster[c]] >= EPS && p_winconn[c] - p_woutconn[c] > EPS*/
		if(p_winconn[c] >= avg_winconn[cluster[c]] && p_winconn[c] > p_woutconn[c]) {
			members[cluster[c]].push_back(mapping[c]+1);
			isCore[c] = true;	
		}
	}
	
	//[END] Preliminary Core
	
	//[BEGIN] Extended Core
	
	for(int c = 0; c < Vnew; c++) {
		if(isCore[c]) continue;
		for(int d = 0; d < G[c].size(); d++) {
			int j = backmap[G[c][d].first];
			if(cluster[c] == cluster[j]) {
				if(isCore[j])
					p_winconnC[c] += G[c][d].second;
				else
					p_woutconnC[c] += G[c][d].second;
			}
		}
		avg_winconnC[cluster[c]] += p_winconnC[c];
	}
	
	for(int c = 0; c < numclusters; c++)
		avg_winconnC[c] /= (2.0*cluster_size[c]);
		
	for(int c = 0; c < Vnew; c++) {
		if(!isCore[c]) { /*p_winconnC[c] - avg_winconnC[cluster[c]] >= EPS && p_winconnC[c] - p_woutconnC[c] > EPS*/
			if(p_winconnC[c] >= avg_winconnC[cluster[c]] && p_winconnC[c] > p_woutconnC[c]) {
				members[cluster[c]].push_back(mapping[c]+1);
				isCore[c] = true;
			}
		}
	}
	//[END] Extended Core

	//[BEGIN] Attachments
	int continuedcount = 0;
	int attachcount = 0;
	for(int c = 0; c < Vnew; c++) {
		if(isCore[c]) {
			for(int d = 0; d < G[c].size(); d++) {
				int j = backmap[G[c][d].first];
				if(isCore[j] && cluster[c] == cluster[j])
					cluster_inter[cluster[c]] += G[c][d].second;
			}
		}
	}
	for(int c = 0; c < numclusters; c++)
		cluster_inter[c] /= 2.0;
	
	for(int c = 0; c < Vnew; c++) {
		if(isCore[c]){
			continuedcount++;
		 	continue;
		}
		double p_interC[numclusters];
		fill(p_interC, p_interC+numclusters, 0.0);
		for(int d = 0; d < G[c].size(); d++) {
			int j = backmap[G[c][d].first];
			if(isCore[j])
				p_interC[cluster[j]] += G[c][d].second;
		}
		for(int d = 0; d < numclusters; d++) {
			if(p_interC[d] > 0.0 && cluster_inter[d] > 0.0) {
				if(p_interC[d] >= alpha * cluster_inter[d] * pow(2.0/cluster_size[d], gamma)) {
					//members[d].push_back(mapping[c]+1);
					mematts[d].push_back(mapping[c]+1);
					attachcount++;
				}
			}
			
		}
	}
	//[END] Attachments

	int clustcount = 1;
	for(int c = 0; c < numclusters; c++) {
		if( (members[c].size() + mematts[c].size()) >=4  ){
			fprintf(complexes, "CLUSTER: %d\nCORE: ", clustcount);
		
			sort(members[c].begin(), members[c].end());
			sort(mematts[c].begin(), mematts[c].end());
			for(int d = 0; d < members[c].size(); d++) {
					if(d > 0) {fprintf(complexes, " ");}
					fprintf(complexes, "%d", members[c][d]);
			}
			fprintf(complexes, "\nATTACHMENTS:");
			for(int d = 0; d < mematts[c].size(); d++) {
				
					if(d > 0) {fprintf(complexes, " ");}
					fprintf(complexes, "%d", mematts[c][d]);
			}
			fprintf(complexes, "\n\n");
			clustcount++;
		}
	}
	
	fclose(complexes);
	fclose(wgraph);
	fclose(clus);
}
