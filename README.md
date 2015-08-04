# SimplePhaseField
A simple and user friendly simulation code designed for use in phase field methods.

## Introduction
SimplePhaseField enables a user to implement their specific phase field model into a larger code-base.
It is generic enough to accomodate the majority of phase field models, yet simple enough for scientific researchers to easily jump in and begin simulating.

The features include:
  - Any number of named order parameters
  - Any number of user-specific model parameters
  - Support for 2d and 3d systems
  - Multidimensional MPI domain decomposition
  - HDF5 data storage
  - Parameter file input
  - Event logging
  - Preprogrammed stencil operations

## How it works
The user modifies the compute kernel to implement the calculation required for their specific model.
This kernel plugs into the larger code base which handles the parallelization, input/output, logging, error handling, etc.

## Install
SimplePhaseField requires that both the MPI and HDF5 libraries be installed on the system. 
Once these are installed, modify the Makefile so that the mpi and hdf5 variables contain the appropriate paths to these libraries. 
Run make.

The install can be tested by cd'ing to the "test" directory and running "mpirun -np 4 ./spf". An sample input and configuration file are included.

## Usage
Three files are required to run a simulation: 
- the input file
- the initial configuration hdf5 file
- the executable

The input file must be named "input.txt". It contains parameter values that will be used in the model implementation as well as a few parameters that are used for the simulation. All input parameters have a key and value separated by an equal sign (e.g. "dx = 1.0"). The user can add as many parameters to the input file as they want. The parameter names must not included spaces or commas. Comment lines can be added by prepending a "#" to the front of the line.

The initial configuration file essentially contains the initial conditions and is the starting point of the simulation. It must be in hdf5 format. Any number of order parameters can be included by adding additional datasets to the root directory of the file. The name of the order parameter, which is chosen by the user, is the same name that is used in the model file to refer to the data.

The executable is called "spf" and is located in the bin directory. In most environments it is executed by running "mpirun -np 4 /path/to/spf" where 4 is the number of processors.

A sucessful simulation output 2 files:
- the log file
- the strand file

The log file includes information relevant to that particular simulation. It contains parameter values, profiling information, and error events.

The strand file is an hdf5 file where all of the order parameter data is stored. Each order parameter has its own directory that takes the name of the order parameter. Within each directory is a list of dataset that were outputted periodically during the simulation.

## Programming
There is only one file that the user is required to edit in order to implement thier model. This file is called model.cpp and is located in the src directory. 
This file already includes code to implement a simple cahn-hilliard type model used to simulate spinodal decomposition.
This can be used as a template and expanded to accomodate more complex models.
Additional details are included as comments in the model.cpp file.

An additional file called "preprocessor.hpp" is in the include directory. It contains several system-wide parameter, such as dimensionality and number of ghost rows, that the user may want to edit.
