FS = FSWeighting
FSE = FSWeight
CA = CAw
CAE = CAw
FSC = FScoring
FSCE = FScoring
IN = biogrid.in
MAPPING = BioGRID_AGaser_graph_mapping
CLUSTER_MAPPING = clusters_with_proteins_list

all:	compile fscore
	
compile:
	g++ $(FS).cpp -o $(FSE)
	g++ $(CA).cpp -o $(CAE)

fscore:
	(cd ./predicted ; g++ $(FSC).cpp -o $(FSCE) ; )
