#include<cstdio>
#include<iostream>
#include<vector>
#include<set>
#include<algorithm>
using namespace std;
#define edge pair<int, double>

/*
[1] Hon Nian Chua, Kang Ning, Wing-Kin Sung, Hon Wai Leong, and Lim-soon Wong. 
    Using indirect protein–protein interactions for protein complex prediction.
    Journal of bioinformatics and computational biology, 6(03):435–466, 2008.
*/

int main(int argc, char **argv) {
	if(argc != 4) {
		printf("Usage: <input filename> <output filename> <threshold>\n");
		return 0;
	}
	double SFSThresh = atof(argv[3]);
	if(SFSThresh == 0.0) {
		printf("Usage: <input filename> <output filename> <threshold>\n");
		return 0;
	}
	
    int V, E, i;
    double navg = 0.0;
	FILE* input = fopen(argv[1], "r");
    FILE* output = fopen(argv[2], "w+");
    fscanf(input, "%d", &V);
    set<int> G[V];
    vector<edge > Gw[V];
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
    	bool found = false;
		set<int> neighbors;
		for(set<int>:: iterator j = G[u].begin(); j != G[u].end(); j++) {
			int k = *j;
			if(k != u) {
				//if(k > u)
				neighbors.insert(k);
				for(set<int>:: iterator l = G[k].begin(); l != G[k].end(); l++) {
					if(*l != u)
					neighbors.insert(*l);
				}
			}
		}
        for(set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
			int v = *it;
            double interuv = 0.0;
            double unionuv = (double)(G[v].size());
            double diffuv = (double)(G[u].size());
            double diffvu = (double)(G[v].size());
            double lambdauv = 0.0;
			double lambdavu = 0.0;
            
           	for(set<int>::iterator it = G[u].begin(); it != G[u].end(); it++) {
                if(G[v].find(*it) != G[v].end()) interuv += 1.0;
            }
            if(interuv == 0.0) continue;
            for(set<int>::iterator it = G[u].begin(); it != G[u].end(); it++) {
                if(G[v].find(*it) == G[v].end()) unionuv += 1.0;
            }
			diffuv -= interuv;
			diffvu -= interuv;
            lambdauv = max(lambdauv, navg - diffuv - interuv);
			lambdavu = max(lambdavu, navg - diffvu - interuv);
            double FSuv = (4.0*interuv*interuv) / ((diffuv+(2.0*unionuv)+lambdauv)*(diffvu+(2.0*unionuv)+lambdavu));
            if(FSuv >= SFSThresh) {
            	if(!found) {
            		fprintf(output, "%d ", u+1);
            		found = true;
            	}
            	fprintf(output, " %d %lf", v+1, FSuv);
            }
            /*if(FSuv >= SFSThresh) { //[1]
            	Gw[u].push_back(edge(v, FSuv));
            	Gw[v].push_back(edge(u, FSuv));
            }*/
        }
        if(found) fprintf(output, "\n");
		if((u+1)%500 == 0)
			printf("Done with node %d.\n", u+1);
    }
    /*for(int c = 0; c < V; c++) {
		if(Gw[c].size() > 0) {
			printf("%d ", c+1);
			for(int d = 0; d < Gw[c].size(); d++) {
				printf(" (%d,%lf)", Gw[c][d].first+1, Gw[c][d].second);
			}
			printf("\n");
		}
    }*/
}

