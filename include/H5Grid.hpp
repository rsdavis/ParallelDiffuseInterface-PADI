
#include <iostream>
#include <string>
#include <vector>
#include "hdf5.h"

//! A simple interface for storing/accessing grid data in hdf5 files

class H5Grid {

    private:
        hid_t file_id;
        hsize_t dims[3];

        void create_group(std::string path);
        void parse (std::string attr_name, std::string &path, std::string &name);

    public:
        H5Grid();
        int open(std::string filename, std::string mode, int &nx, int &ny, int &nz);
        int read_dataset(std::string dataset_name, double * dataset);
        int write_dataset(std::string dataset_name, double * dataset);
        int set_attribute(std::string attr_name, int attr_value);
        int get_attribute(std::string attr_name, int &attr_value);
        int list(std::string path, std::vector<std::string> &list);
        int close();
};

H5Grid :: H5Grid () {
    file_id = 0;
}

void H5Grid :: create_group(std::string path)
{
    int start = 1;
    int pos = 0;
    if (path[0]!='/') path = "/" + path;

    while ((pos = path.find("/", start)) != std::string::npos)
    {
        std::string substr = path.substr(0, pos);
        htri_t group_exists = H5Lexists(file_id, substr.c_str(), H5P_DEFAULT);
        if (!group_exists) {
            hid_t grp_id = H5Gcreate(file_id, substr.c_str() , H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Gclose(grp_id);
        }
        start = pos+1;
    }
}


int H5Grid :: open(std::string filename, std::string mode, int &nx, int &ny, int &nz)
{
    /**
    @param filename Name of the hdf5 file
    @param mode of file access "r", "w", or "a"
    @param nx number of grid points in x dimension
    @param ny number of grid points in y dimension
    @param nz number of grid points in z dimension
    **/

    /**
    nx, ny, nz must be variables, and cannot be constant integers (e.g. "nx" instead of "100")
    For 2D data, nz should be set to 0.
    **/

    /** \error ERROR 1: mode must be either "r", "w", or "a" */
    /** \error ERROR 2: a file is already open */
    /** \error ERROR 3: mode is "r" or "a" but the file does not exist */
    /** \error ERROR 4: file is not in hdf5 format */

    unsigned read, write, append;

    read = (mode == "r");
    write = (mode == "w");
    append = (mode == "a");

    // check that one is true
    if ( (read || write || append) != 1) return 1;

    // check that a file is not already open
    if (file_id != 0) return 2;


    if ( write ) {
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        set_attribute("/Nx", nx);
        set_attribute("/Ny", ny);
        set_attribute("/Nz", nz);

    } else if (read || append) {
        FILE * file_exists = fopen(filename.c_str(), "r");
        if (file_exists == NULL) return 3;
        else fclose(file_exists);

        //this doesn't work for some reason
        //if (H5Fis_hdf5(filename.c_str()) < 1) return 4;

        if (read)
            file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (append)
            file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

        get_attribute("/Nx", nx);
        get_attribute("/Ny", ny);
        get_attribute("/Nz", nz);

    }

    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;

    return 0;
}

int H5Grid :: read_dataset(std::string dataset_name, double * dataset)
{
    /**
    * @param dataset_name name of the dataset
    * @param dataset pointer to the first element of data
    **/

    hid_t data_id; 

    //check that dataset exists
    /* requires that you check every group in the path
    htri_t dataset_exists = H5Lexists(file_id, dataset_name.c_str(), H5P_DEFAULT);
    if (dataset_exists <= 0) return 1;
    */

    data_id = H5Dopen(file_id, dataset_name.c_str(), H5P_DEFAULT);
    H5Dread(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset);

    H5Dclose(data_id);

    return 0;
}

int H5Grid :: write_dataset(std::string dataset_name, double * dataset)
{
    /**
    * @param dataset_name name of the dataset
    * @param dataset pointer to the first element of data
    **/

    /** \error ERROR 1: dataset already existed and was overwritten */
    /** \error ERROR 2: internal error when writing the dataset */

    hid_t dcpl_id, space_id, data_id;
    herr_t error;
    int rank;

    if (dims[2] == 0) rank = 2;
    else rank = 3;

    std::string path;
    std::string name;

    // create groups in the path
    create_group(dataset_name);

    // check if dataset already exisits
    htri_t dataset_exists = H5Lexists(file_id, dataset_name.c_str(), H5P_DEFAULT);

    if (dataset_exists==0) { // create dataset

        // chunk for data compression
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl_id, rank, dims);
        H5Pset_shuffle(dcpl_id);
        H5Pset_deflate(dcpl_id, 1);

        space_id = H5Screate_simple(rank, dims, NULL);
        data_id = H5Dcreate(file_id, 
                            dataset_name.c_str(), 
                            H5T_NATIVE_DOUBLE, 
                            space_id, 
                            H5P_DEFAULT, 
                            dcpl_id, 
                            H5P_DEFAULT);

    } else { // dataset_exists, just open and overwrite data
        data_id = H5Dopen(file_id, dataset_name.c_str(), H5P_DEFAULT);
    }

    error = H5Dwrite(data_id,
                     H5T_NATIVE_DOUBLE,
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     dataset);

    if (dataset_exists) return 1;
    if (error < 0) return 2;

    return 0;
}

void H5Grid :: parse (std::string attr_name, std::string &path, std::string &name)
{
    size_t start = 0;
    size_t end = 0;

    while (end != std::string::npos)
    {
        end = attr_name.find("/", start);
        path = attr_name.substr(0,start);
        name = attr_name.substr(start,end-start);
        start = end + 1;
    }
}

int H5Grid :: set_attribute(std::string attr_name, int  attr_value)
{
    /**
    * @param attr_name Name/address of the attribute
    * @param attr_value Value of the attribute
    **/

    /**
    The attr_name must have atleast one backslash.
    For example, setting an attribute named my_attr in the root group
    should be specified as "/my_attr"
    */

    /** \error ERROR 1: attr_name not sufficient, ensure that it begins with a backslash */
    /** \error ERROR 2: attr_name not sufficient, ensure that it does not end in a backslash */

    std::string path;
    std::string name;
    parse(attr_name, path, name);

    if (path.size() == 0) return 1;
    if (name.size() == 0) return 2;

    create_group(attr_name);
    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate_by_name(file_id, 
                              path.c_str(),
                              name.c_str(), 
                              H5T_NATIVE_INT, 
                              space_id, 
                              H5P_DEFAULT, 
                              H5P_DEFAULT,
                              H5P_DEFAULT);

    H5Awrite(attr_id, H5T_NATIVE_INT, &attr_value);
    H5Aclose(attr_id);
    H5Sclose(space_id);

    return 0;
}

int H5Grid :: get_attribute(std::string attr_name, int &attr_value)
{
    /**
    * @param attr_name Name/address of the attribute
    * @param attr_value Value of the attribute
    **/

    /** \error ERROR 1: attribute does not exist */

    hid_t attr_id;
    std::string path;
    std::string name;

    parse(attr_name, path, name);

    if (H5Aexists_by_name(file_id, path.c_str(), name.c_str(), H5P_DEFAULT) <= 0) return 1;

    attr_id = H5Aopen_by_name(file_id, path.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    H5Aread(attr_id, H5T_NATIVE_INT, &attr_value);

    H5Aclose(attr_id);

    return 0;
}

int H5Grid :: list(std::string path, std::vector<std::string> &list)
{

    /**
    @param path list all the dataset in this path
    @param list modified to contain the names of the datasets
    **/

    hid_t group_id;
    H5G_info_t group_info;

    // need to implement group existence check

    group_id = H5Gopen(file_id, path.c_str(), H5P_DEFAULT);
    H5Gget_info(group_id, &group_info);

    list.resize(group_info.nlinks);

    for (int i=0; i<group_info.nlinks; i++)
    {
        char name[256];
        H5Gget_objname_by_idx(group_id, i, name, 256);
        list[i] = name;
    }


    H5Gclose(group_id);

    return 0;
}

int H5Grid :: close()
{
    /** \error ERROR 1: there is no file to close **/

    if (file_id == 0) {
        // no file opened
        return 1;
    } else {
        H5Fclose(file_id);
        file_id = 0;
    }

    return 0;
}
