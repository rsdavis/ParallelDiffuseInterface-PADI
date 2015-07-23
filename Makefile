
mpi = /Volumes/DATA/Software/openmpi/1.8.4
hdf5 = /usr/local/Cellar/hdf5/1.8.14

default:
	mpic++ -Wall main.cpp \
	-I$(mpi)/include -L$(mpi)/lib -lmpi  \
	-I$(hdf5)/include -L$(hdf5)/lib -lhdf5
