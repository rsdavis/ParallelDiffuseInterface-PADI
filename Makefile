
debug = 2
std = c++11

mpi = /Volumes/DATA/Software/openmpi/1.8.4
hdf5 = /usr/local/Cellar/hdf5/1.8.14

all: bin/spf

bin/spf: bin/main.o
	mpic++ -Wall bin/model.o bin/main.o -o bin/spf -DFILELOG_MAX_LEVEL=$(debug) \
		-L $(mpi)/lib -lmpi \
		-L $(hdf5)/lib -lhdf5 \

bin/main.o: include/* bin/model.o src/main.cpp
	mpic++ -Wall -std=$(std) -c src/main.cpp -o bin/main.o -DFILELOG_MAX_LEVEL=$(debug) \
		-I $(mpi)/include \
		-I $(hdf5)/include \
		-I ./include

bin/model.o: include/* src/model.cpp
	mpic++ -Wall -std=$(std) -c src/model.cpp -o bin/model.o -DFILELOG_MAX_LEVEL=$(debug) \
		-I ./include


clean:
	rm bin/*
