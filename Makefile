FS = src/FSWeighting
FSE = FSWeight
CA = src/CAw
CAE = CAw
FSC = FScoring
FSCE = FScoring

all:	compile fscore
	
compile:
	g++ $(FS).cpp -o $(FSE)
	g++ $(CA).cpp -o $(CAE)

fscore:
	(cd scoring; make )

