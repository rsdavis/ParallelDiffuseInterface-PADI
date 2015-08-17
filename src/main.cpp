
#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm> // std::remove

#include "log.h"
#include "MPIGrid.hpp"
#include "H5Grid.hpp"
#include "model.hpp"
#include "preprocessor.hpp"
#include "mpitimer.h"

std::vector<std::string> parse_by_comma(std::string str)
{
    std::stringstream ss(str);
    std::vector<std::string> result;
    while (ss.good())
    {
        std::string substr;
        std::getline( ss, substr, ',' );
        result.push_back( substr );
    }
    return result;
}

void readParameters(std::map<std::string, std::string> &params)
{
    /** 
    This function parses the input file.
    It reads a key-value pair from any line involving an "=" sign.
    All key-value are stored as string in params for use later.
    Comment lines begin with "#"
    */

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
                    double *& global_phase, 
                    int * global_dims, 
                    std::vector<std::string> &name_list)
{
    int err;
    int vol;
    int ndims=1;

    /* Open HDF5 Initial Configuration File */

    FILE_LOG(logDEBUG) << "Reading Initial Configuration file = " << filename;

    H5Grid h5;
    err = h5.open(filename, "r", global_dims, ndims);

    if (err > 0) {
        FILE_LOG(logERROR) << "Error opening file " << filename;
        FILE_LOG(logERROR) << "H5GRID OPEN CODE: " << err;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    FILE_LOG(logDEBUG) << "global_dims[0]=" << global_dims[0];
    FILE_LOG(logDEBUG) << "global_dims[1]=" << global_dims[1];

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

    vol = 1;
    for (int i=0; i<SPF_NDIMS; i++) vol *= global_dims[i];

    FILE_LOG(logDEBUG) << "Allocating global data: " << vol;

    try {
        global_phase = new double [vol*name_list.size()];
    } catch (std::bad_alloc &ba) {
        FILE_LOG(logERROR) << "Bad Alloc Global Data: " << ba.what();
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    /* Read Data */

    for (int i=0; i<name_list.size(); i++) 
    {
        int offset = vol*i;
        err = h5.read_dataset(name_list[i], global_phase + offset);

        if (err > 0) {
            FILE_LOG(logERROR) << "Error reading datasets from file " << filename;
            FILE_LOG(logERROR) << "H5GRID READ_DATASET CODE: " << err;
            MPI_Abort(MPI_COMM_WORLD, 0);
        };
    }

    err = h5.close();
}

std::string output_path(std::string phase_name, int num)
{
    std::ostringstream ss;
    ss << "/" << phase_name << "/" << std::setw(6) << std::setfill('0') << num;
    return ss.str();
}

void updateLog(int istep, int nsteps, double elapsed_time)
{
    char buffer[80];
    time_t current_time, eta;
    struct tm * timeinfo;

    double fraction_completed = ((double)istep) / nsteps;
    double remaining_time = (1.0/fraction_completed -1.0)*elapsed_time;

    time(&current_time);
    eta = current_time + (int) remaining_time;

    timeinfo = localtime(&eta);
    strftime(buffer, 80, "%c", timeinfo);
    FILE_LOG(logINFO) << (int) (fraction_completed*100) << "% Complete, ETA: " << buffer;
    //FILE_LOG(logINFO) << (int) (fraction_completed*100) << "\% Complete, ETA: " << asctime(timeinfo);
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

    double * global_phase;
    double * local_phase;
    double * local_chem_pot;
    double * local_mobility;

    char * phase_names;
    int nphases;

    std::map<std::string, int> name_index;
    std::vector<int> output_phase;
    std::vector<int> output_mobility;
    std::vector<int> output_chem_pot;

    mpitimer * comm_time = mpitimer_new();
    mpitimer * comp_time = mpitimer_new();
    mpitimer * io_time = mpitimer_new();
    mpitimer * total_time = mpitimer_new();

    mpitimer_start(total_time);


    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    logStart();
    //Log::ReportingLevel() = logINFO;
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
        readGlobalData(params["init_file"], global_phase, global_dims, name_list);

        nphases = name_list.size();
        phase_names = new char [100*nphases];
        for (int i=0; i<nphases; i++)
            strncpy(phase_names+i*100, name_list[i].c_str(), 100);
    }

    // broadcast names of order parameters

    MPI_Bcast(&nphases, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank!=0) phase_names = new char [100*nphases];
    MPI_Bcast(phase_names, nphases*100, MPI_CHAR, 0, MPI_COMM_WORLD);

    // setup grid

    MPI_Bcast(&global_dims, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&np_dims, 3, MPI_INT, 0, MPI_COMM_WORLD);

    MPIGrid grid;
    for (int i=0; i<SPF_NDIMS; i++) np_dims[i] = 0;
    MPI_Dims_create(np, SPF_NDIMS, np_dims);
    err = grid.setup(MPI_COMM_WORLD, global_dims, np_dims, SPF_NDIMS, SPF_NROWS, local_dims);
    if (err > 0) {
        FILE_LOG(logERROR) << "MPIGRID SETUP CODE: " << err;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    if (rank == 0) {
        FILE_LOG(logINFO);
        FILE_LOG(logINFO) << "Number of Processors: " << np;
        FILE_LOG(logINFO) << "Number of Dimensions: " << SPF_NDIMS;
        FILE_LOG(logINFO) << "Global Grid Dimensions: " << global_dims[0] <<","<< global_dims[1] <<","<< global_dims[2];
        FILE_LOG(logINFO) << "Processor Dimensions: " << np_dims[0] <<","<< np_dims[1] <<","<< np_dims[2];
        FILE_LOG(logINFO) << "Local Grid Dimensions: " << local_dims[0] <<","<< local_dims[1] <<","<< local_dims[2];
        FILE_LOG(logINFO);
    }

    // allocate local data

    int local_volume = 1;
    int global_volume = 1;
    for (int i=0; i<SPF_NDIMS; i++)
    {
        local_volume *= local_dims[i];
        global_volume *= global_dims[i];
    }

    local_phase = new double [local_volume*nphases];
    local_chem_pot = new double [local_volume*nphases];
    local_mobility = new double [local_volume*nphases];

    // create map of names and index
    std::vector<std::string> out_p = parse_by_comma(params["output_phase"]);
    std::vector<std::string> out_m = parse_by_comma(params["output_mobility"]);
    std::vector<std::string> out_c = parse_by_comma(params["output_chem_pot"]);

    output_phase.resize(nphases);
    output_mobility.resize(nphases);
    output_chem_pot.resize(nphases);


    for (int i=0; i<nphases; i++)
    {
        std::string name(phase_names+i*100);
        name_index[name] = i;

        output_phase[i] = 0;
        output_mobility[i] = 0;
        output_chem_pot[i] = 0;

        if (std::find(out_p.begin(), out_p.end(), name) != out_p.end()) 
            output_phase[i] = 1;
        if (std::find(out_m.begin(), out_m.end(), name) != out_m.end()) 
            output_mobility[i] = 1;
        if (std::find(out_c.begin(), out_c.end(), name) != out_c.end()) 
            output_chem_pot[i] = 1;
    }

    // scatter and share data

    for (int i=0; i<nphases; i++)
        grid.scatter(global_phase+global_volume*i, local_phase+local_volume*i);

    for (int i=0; i<nphases; i++)
        grid.share(local_phase + local_volume*i);

    // create alias for more user-friendly data access

    double ** data_alias = new double * [nphases]; 
    double ** chem_pot_alias = new double * [nphases];
    double ** mobility_alias = new double * [nphases];

    for (int i=0; i<nphases; i++)
    {
        data_alias[i] = local_phase + local_volume*i;
        chem_pot_alias[i] = local_chem_pot + local_volume*i;
        mobility_alias[i] = local_mobility + local_volume*i;
    }

    if (rank==0) {
        H5Grid h5;
        int ndims = SPF_NDIMS;
        h5.open("strand.h5", "w", global_dims, ndims);
        for (int i=0; i<nphases; i++)
        {
            std::string path = output_path(phase_names+i*100, 0);
            h5.write_dataset(path, global_phase + global_volume*i);
        }
        h5.close();
    }

    preprocess(data_alias, local_dims, params, name_index);

    // begin stepping in time

    for (int istep=1; istep<=std::stoi(params["nsteps"]); istep++)
    {
        FILE_LOG(logDEBUG) << "Step: " << istep;

        mpitimer_start(comp_time);
        kernel(data_alias, chem_pot_alias, mobility_alias, local_dims);
        mpitimer_stop(comp_time);

        FILE_LOG(logDEBUG) << "Share data";
        mpitimer_start(comm_time);
        for (int i=0; i<nphases; i++)
            grid.share(local_phase + local_volume*i);
        mpitimer_stop(comm_time);

        // output
        if (istep % std::stoi(params["output_frequency"]) == 0) {

            FILE_LOG(logDEBUG) << "Output";
            
            H5Grid h5;
            int stat, ndims = SPF_NDIMS;
            double * buffer; 
            int frame;
            mpitimer_start(io_time);

            if (rank == 0) h5.open("strand.h5", "a", global_dims, ndims);

            buffer = new double [global_volume];
            frame = istep / std::stoi(params["output_frequency"]);

            // cycle through order parameters
            for (int i=0; i<nphases; i++) {

                std::string name = phase_names+i*100;

                // output phases
                grid.gather(global_phase+global_volume*i, local_phase+local_volume*i);
                if (output_phase[i] && rank == 0) {
                        std::string path = output_path(name, frame);
                        stat = h5.write_dataset(path, global_phase+global_volume*i);
                }

                // output mobility
                if (output_mobility[i]) {
                    grid.gather(buffer, local_mobility+local_volume*i);
                    if (rank==0) {
                        std::string mod = "_mobility";
                        std::string path = output_path(name+mod, frame);
                        stat = h5.write_dataset(path, buffer);
                    }
                }

                // output_chem_pot
                if (output_chem_pot[i]) {
                    grid.gather(buffer, local_chem_pot+local_volume*i);
                    if (rank==0) {
                        std::string mod = "_chem_pot";
                        std::string path = output_path(name+mod, frame);
                        stat = h5.write_dataset(path, buffer);
                    }
                }
            }

            if (rank == 0) h5.close();

            /*
            if (rank == 0) {
                H5Grid chckpt;
                stat = chckpt.open("checkpoint.h5", "w", global_dims, ndims);
                if (stat != 0) {FILE_LOG(logERROR) << "H5Grid open checkpoint: " << stat;}
                for (int i=0; i<nphases; i++)
                {
                    std::string name = phase_names+i*100;
                    std::cout << name << std::endl;
                    stat = chckpt.write_dataset(name, global_phase+global_volume*i);
                    if (stat != 0) {FILE_LOG(logERROR) << "H5Grid write checkpoint: " << stat;}
                }
                stat = chckpt.close();
                    if (stat != 0) {FILE_LOG(logERROR) << "H5Grid close checkpoint: " << stat;}
            }
            */

            mpitimer_stop(io_time);
        } // end output

    } // end timesteping

    if (rank == 0) {
        H5Grid chckpt;
        int ndims = SPF_NDIMS;
        int stat = chckpt.open("checkpoint.h5", "w", global_dims, ndims);
        if (stat != 0) {FILE_LOG(logERROR) << "H5Grid open checkpoint: " << stat;}
        for (int i=0; i<nphases; i++)
        {
            std::string name = phase_names+i*100;
            std::cout << name << std::endl;
            stat = chckpt.write_dataset(name, global_phase+global_volume*i);
            if (stat != 0) {FILE_LOG(logERROR) << "H5Grid write checkpoint: " << stat;}
        }
        stat = chckpt.close();
        if (stat != 0) {FILE_LOG(logERROR) << "H5Grid close checkpoint: " << stat;}
    }

    FILE_LOG(logDEBUG) << "Postprocess";
    postprocess(data_alias, chem_pot_alias, mobility_alias, local_dims);

    if (rank == 0) {
        mpitimer_stop(total_time);

        int total_seconds = (int) mpitimer_get_time(total_time);
        int days = total_seconds / (60*60*24);
        int hours = (total_seconds % (60*60*24)) / 3600;
        int mins = ((total_seconds % (60*60*24)) % 3600) / 60;
        int secs = ((total_seconds % (60*60*24)) % 3600) % 60;

        FILE_LOG(logINFO);
        FILE_LOG(logINFO) << "Computation Time (%):   " << mpitimer_get_time(comp_time) / mpitimer_get_time(total_time) * 100;
        FILE_LOG(logINFO) << "Communication Time (%): " << mpitimer_get_time(comm_time) / mpitimer_get_time(total_time) * 100;
        FILE_LOG(logINFO) << "Input/Output Time (%):  " << mpitimer_get_time(io_time) / mpitimer_get_time(total_time) * 100;
        FILE_LOG(logINFO) << "Total Time (s):     " << mpitimer_get_time(total_time);
        FILE_LOG(logINFO) << "Elapsed Time: " << days << " days, " << hours << " hours, " << mins << " minutes, " << secs << " seconds";
    }


    delete [] data_alias;
    delete [] chem_pot_alias;
    delete [] mobility_alias;

    delete [] local_phase;
    delete [] local_chem_pot;
    delete [] local_mobility;

    if (rank == 0) {
        delete [] global_phase;
    }

    delete [] phase_names;


    MPI_Finalize();

    return 0;

}

