WORKDIR=$(PWD)
BIN=$(WORKDIR)/bin

UNAME := $(shell uname -s)

ifeq ($(UNAME), Linux)
# do something Linux-y
machine=Linux
endif
ifeq ($(UNAME), Darwin)
# do something Solaris-y
machine=OSX
endif

FS = src/FSWeighting
FSE = FSWeight
CA = src/CAw
CAE = CAw
FSC = FScoring
FSCE = FScoring

all:	compile fscore
	
compile:
	g++ $(FS).cpp -o $(BIN)/$(machine)/$(FSE)
	g++ $(CA).cpp -o $(BIN)/$(machine)/$(CAE)

fscore:
	(cd scoring; make )

