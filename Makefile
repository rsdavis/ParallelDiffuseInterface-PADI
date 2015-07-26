
mpi = /Volumes/DATA/Software/openmpi/1.8.4
hdf5 = /usr/local/Cellar/hdf5/1.8.14

default:
	mpic++ -Wall model.cpp main.cpp -DFILELOG_MAX_LEVEL=2 \
	-I$(mpi)/include -L$(mpi)/lib -lmpi  \
	-I$(hdf5)/include -L$(hdf5)/lib -lhdf5
