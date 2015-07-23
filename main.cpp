
#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <ctime>

#include "log.h"
#include "mpigrid.hpp"
#include "h5grid.hpp"
#include "model.hpp"
#include "preprocessor.hpp"


void readParameters(std::map<std::string, std::string> &params)
{
    std::ifstream input("input.txt");
    std::string line;

    while (std::getline(input, line))
    {
        size_t pos;
        std::string key, value;

        // remove whitespace before parsing
        line.erase(std::remove(line.begin(),line.end(),' '), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\t'), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\n'), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\r'), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\v'), line.end());

        // skip empty lines
        if (line.empty()) continue;

        // skip comments
        if (line[0] == '#') continue;

        // parse line using the equal sign
        pos = line.find("=");
        key = line.substr(0, pos);
        value = line.substr(pos+1, line.size());

        params[key] = value;
    }

    input.close();
}

void logStart()
{
    time_t start_time = time(0);
    FILE * logFile = fopen("log.txt", "w");
    Output2FILE::Stream() = logFile;
    FILE_LOG(logINFO) << "Beginning simulation";
    FILE_LOG(logINFO) << std::ctime(&start_time);
}

void logParameters(std::map<std::string, std::string> params)
{
    std::map<std::string, std::string> :: iterator it;
    for (it = params.begin(); it!=params.end(); it++)
        FILE_LOG(logINFO) << "Parameter:"<< it->first << ":" << it->second << ":";
}

void readGlobalData(std::string filename, 
                    double *& global_data, 
                    int * global_dims, 
                    int &ndims, 
                    std::vector<std::string> &name_list)
{
    int err;
    int nx, ny, nz, vol;

    /* Open HDF5 Initial Configuration File */

    H5Grid h5;
    err = h5.open(filename, "r", nx, ny, nz);

    if (err > 0) {
        FILE_LOG(logERROR) << "Error opening file " << filename;
        FILE_LOG(logERROR) << "H5GRID OPEN CODE: " << err;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    if (ny > 0) ndims = 2;
    if (nz > 0) ndims = 3;

    if (ndims != SPF_NDIMS) {
        FILE_LOG(logERROR) << "Data has dimensionality = " << ndims;
        FILE_LOG(logERROR) << "Simulation has dimensionality = " << SPF_NDIMS;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    /* Get List Of Names Of The Order Parameters From HDF5 File */

    err = h5.list("/", name_list);
    if (err > 0)  {
        FILE_LOG(logERROR) << "Error reading order parameters from file " << filename;
        FILE_LOG(logERROR) << "H5GRID LIST CODE: " << err;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    /* Allocate Global Data */

    global_dims[0] = nx;
    global_dims[1] = ny;
    global_dims[2] = nz;

    vol = 1;
    for (int i=0; i<ndims; i++) vol *= global_dims[i];

    try {
        global_data = new double [vol*name_list.size()];
    } catch (std::bad_alloc &ba) {
        FILE_LOG(logERROR) << "Bad Alloc Global Data: " << ba.what();
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    /* Read Data */

    for (int i=0; i<name_list.size(); i++) 
    {
        int offset = vol*i;
        err = h5.read_dataset(name_list[i], global_data + offset);

        if (err > 0) {
            FILE_LOG(logERROR) << "Error reading datasets from file " << filename;
            FILE_LOG(logERROR) << "H5GRID READ_DATASET CODE: " << err;
            MPI_Abort(MPI_COMM_WORLD, 0);
        };
    }

    err = h5.close();
}

int main(int argc, char ** argv)
{

    MPI_Init(&argc, &argv);

    int rank;
    int np;

    int err;
    int global_dims[3];
    int local_dims[3] = {1,1,1};
    int np_dims[3];

    double * global_data;
    double * local_data;
    double * local_chem_pot;
    double * local_mobility;

    int ndims;
    char * phase_names;
    int nphases;
    std::map<std::string, int> name_index;


    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    logStart();
    std::map<std::string, std::string> params;

    MPI_Status status;

    // the processes take turns reading input file

    if (rank != 0) 
        MPI_Recv(NULL, 0, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);

    readParameters(params);

    if (rank != np-1) 
        MPI_Send(NULL, 0, MPI_INT, rank+1, 1, MPI_COMM_WORLD);

    if (rank == 0) logParameters(params);

    // master reads initial configuration

    if (rank == 0) {
        std::vector<std::string> name_list;
        readGlobalData(params["init_file"], global_data, global_dims, ndims, name_list);

        nphases = name_list.size();
        phase_names = new char [100*nphases];
        for (int i=0; i<nphases; i++)
            strncpy(phase_names+i*100, name_list[i].c_str(), name_list[i].length());
    }


    // broadcast names of order parameters

    MPI_Bcast(&nphases, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank!=0) phase_names = new char [100*nphases];
    MPI_Bcast(phase_names, nphases*100, MPI_CHAR, 0, MPI_COMM_WORLD);

    // setup grid

    MPI_Bcast(&ndims, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_dims, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&np_dims, 3, MPI_INT, 0, MPI_COMM_WORLD);

    MPIGrid grid;
    MPI_Dims_create(np, 2, np_dims);
    err = grid.setup(MPI_COMM_WORLD, global_dims, np_dims, ndims, SPF_NROWS, local_dims);
    if (err > 0) {
        FILE_LOG(logERROR) << "MPIGRID SETUP CODE: " << err;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    // allocate local data

    int local_volume = 1;
    int global_volume = 1;
    for (int i=0; i<ndims; i++)
    {
        local_volume *= local_dims[i];
        global_volume *= global_dims[i];
    }

    local_data = new double [local_volume*nphases];
    local_chem_pot = new double [local_volume*nphases];
    local_mobility = new double [local_volume*nphases];

    // create map of names and index

    for (int i=0; i<nphases; i++)
    {
        std::string name(phase_names+i*100);
        name_index[name] = i;
    }

    // scatter and share data

    for (int i=0; i<nphases; i++)
        grid.scatter(global_data+global_volume*i, local_data+local_volume*i);

    for (int i=0; i<nphases; i++)
        grid.share(local_data + local_volume*i);

    // create user friendly alias for data access

    double ** data_alias = new double * [nphases]; 
    double ** chem_pot_alias = new double * [nphases];
    double ** mobility_alias = new double * [nphases];

    for (int i=0; i<nphases; i++)
    {
        data_alias[i] = local_data + local_volume*i;
        chem_pot_alias[i] = local_chem_pot + local_volume*i;
        mobility_alias[i] = local_mobility + local_volume*i;
    }


    preprocess(local_data, local_dims, params, name_index);

    // begin stepping in time

    for (int istep=1; istep<=std::stoi(params["nsteps"]); istep++)
    {
        integrate(data_alias, chem_pot_alias, mobility_alias, local_dims);

        for (int i=0; i<nphases; i++)
            grid.share(local_data + local_volume*i);
    }

    for (int i=0; i<nphases; i++)
        grid.gather(global_data+global_volume*i, local_data+local_volume*i);

    if (rank == 0) {
        H5Grid h5;
        int nx = global_dims[0];
        int ny = global_dims[1];
        int nz = 0;
        h5.open("out.h5", "w", nx, ny, nz);
        std::string name(phase_names);
        h5.write_dataset(name, global_data);
        h5.close();
    }

    free(data_alias);
    free(chem_pot_alias);
    free(mobility_alias);
    free(local_data);
    free(local_chem_pot);
    free(local_mobility);

    MPI_Finalize();

    return 0;

}

