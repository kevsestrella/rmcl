#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<algorithm>
#include<map>
#include<set>
#include<vector>
using namespace std;
#define BUFFSZ 1024

int main(int argc, char **argv) {
	if(argc != 4) {
		printf("Usage: %s <mapping file> <clusters file> <complexes file>\n", argv[0]);
		return 0;
	}

	FILE* mapping = fopen(argv[1], "r");
    
    if(mapping == NULL) {
    	perror("");
    	return 0;
    }
    
    FILE* clusters = fopen(argv[2], "r");
    
    if(clusters == NULL) {
    	perror("");
    	return 0;
    }
    
    FILE* complexes = fopen(argv[3], "r");
    
    if(complexes == NULL) {
    	perror("");
    	return 0;
    }
    
	char str[BUFFSZ];
	map<string,int> M;
	int assign = 1;
	while(fscanf(mapping, "%s", str)!=-1) {
		string S(str);
		M[str] = assign++;
	}
	vector<set<int> > cyc;
	while(fgets(str, BUFFSZ, complexes) != NULL) {
		cyc.push_back(set<int>());
		char *spl = strtok(str, " ");
		while(spl != NULL) {
			string S(spl);
			if(M.find(S) != M.end())
				cyc[cyc.size()-1].insert(M[S]);
			spl = strtok(NULL, " ");
		}
	}
	vector<double> ppv;
	vector<double> tj;
	while(fgets(str, BUFFSZ, clusters) != NULL) {
		ppv.push_back(0.0);
		tj.push_back(0.0);
		vector<int> cl;
		char *spl = strtok(str, " ");
		while(spl != NULL) {
			cl.push_back(atoi(spl));
			spl = strtok(NULL, " ");
		}
		vector<double> found;
		for(int c = 0; c < cyc.size(); c++) {
			found.push_back(0.0);
			for(int d = 0; d < cl.size(); d++) {
				if(cyc[c].find(cl[d]) != cyc[c].end()) {
					found[found.size()-1] += 1.0;
					tj[tj.size()-1] += 1.0;
				}
			}
		}
		for(int c = 0; c < found.size(); c++)
			ppv[ppv.size()-1] = max(ppv[ppv.size()-1], found[c]/tj[tj.size()-1]);
	}
	double num = 0.0;
	double denum = 0.0;
	for(int c = 0; c < ppv.size(); c++) {
		num += ppv[c]*tj[c];
		denum += tj[c];
	}
	printf("%lf\n", num/denum);
}
