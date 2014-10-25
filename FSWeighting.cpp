#include<cstdio>
#include<iostream>
#include<vector>
#include<set>
#include<algorithm>
#include<map>
using namespace std;
#define edge pair<int, double>

/*
[1] Hon Nian Chua, Kang Ning, Wing-Kin Sung, Hon Wai Leong, and Lim-soon Wong. 
    Using indirect protein–protein interactions for protein complex prediction.
    Journal of bioinformatics and computational biology, 6(03):435–466, 2008.
*/

int main(int argc, char **argv) {
	if(argc != 5) {
		printf("Usage: <input filename> <output filename> <mlr filename> <threshold>\n");
		return 0;
	}
	double SFSThresh = atof(argv[4]);
	if(SFSThresh == 0.0) {
		printf("Usage: <input filename> <output filename> <mlr filename> <threshold>\n");
		return 0;
	}
	
	FILE* input = fopen(argv[1], "r");
    
    if(input == NULL) {
    	perror("");
    	return 0;
    }
    
    FILE* output = fopen(argv[2], "w+");
    
    if(output == NULL) {
    	perror("");
    	return 0;
    }
    
    FILE* mlrformat = fopen(argv[3], "w+");
    
    if(mlrformat == NULL) {
    	perror("");
    	return 0;
    }
	
    int V, i;
    int Vnew = 0, Enew = 0;
    double navg = 0.0;
    fscanf(input, "%d", &V);
    set<int> G[V];
    vector<edge > Gw[V];
    vector<int> active;
    map<int,int> mapping;
    for(int c = 0; c < V; c++) {
        while(true) {
            fscanf(input, "%d", &i);
            if(i == -1) break;
            G[c].insert(i-1);
        }
        navg += (double)(G[c].size());
        G[c].insert(c);
    }
    navg /= (double)(V);
    for(int u = 0; u < V; u++) {
		set<int> neighbors;
		for(set<int>:: iterator j = G[u].begin(); j != G[u].end(); j++) {
			int k = *j;
			if(k != u) {
				if(k > u)
				neighbors.insert(k);
				for(set<int>:: iterator l = G[k].begin(); l != G[k].end(); l++) {
					if(*l > u)
					neighbors.insert(*l);
				}
			}
		}
        for(set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
			int v = *it;
            double interuv = 0.0;
            double diffuv = (double)(G[u].size());
            double diffvu = (double)(G[v].size());
            double lambdauv = 0.0;
			double lambdavu = 0.0;
            
           	for(set<int>::iterator it = G[u].begin(); it != G[u].end(); it++) {
                if(G[v].find(*it) != G[v].end()) interuv += 1.0;
            }
            if(interuv == 0.0) continue;
			diffuv -= interuv;
			diffvu -= interuv;
            lambdauv = max(lambdauv, navg - diffuv - interuv);
			lambdavu = max(lambdavu, navg - diffvu - interuv);
            double FSuv = (4.0*interuv*interuv) / ((diffuv+(2.0*interuv)+lambdauv)*(diffvu+(2.0*interuv)+lambdavu));
            if(FSuv >= SFSThresh) {
            	Enew++;
            	if(Gw[u].empty()) { Vnew++; active.push_back(u+1); }
            	if(Gw[v].empty()) { Vnew++; active.push_back(v+1); }
            	Gw[u].push_back(edge(v, FSuv));
            	Gw[v].push_back(edge(u, FSuv));
            }
        }
		if((u+1)%500 == 0)
			printf("Done with node %d.\n", u+1);
    }
    sort(active.begin(), active.end());
    for(int c = 0; c < active.size(); c++)
    	mapping[active[c]] = c+1;
    
    fprintf(output, "%d %d %d\n", V, Vnew, Enew);
    fprintf(mlrformat, "%d %d 1\n", Vnew, Enew);
    
    for(int c = 0; c < V; c++) {
    	if(!Gw[c].empty()) {
    		fprintf(output, "%d", c+1);
    		for(int d = 0; d < Gw[c].size(); d++) {
    			fprintf(output, " %d %lf", Gw[c][d].first+1, Gw[c][d].second);
    			if(d > 0) fprintf(mlrformat, " ");
    			fprintf(mlrformat, "%d 1", mapping[Gw[c][d].first+1]);
    		}
    		fprintf(output, "\n");
    		fprintf(mlrformat, "\n");
    	}
    }
    fclose(input);
    fclose(output);
    fclose(mlrformat);
}

