
#include <string>

#include "mpi.h"

/// MPIGrid is a class that enables the user to easily perform N-dimensional domain decomposition

/**
MPIGrid simplifies distributed memory parallelization for applications
that run on rectilinear grids.
There are four main functions (setup, scatter, share, gather) that take global system data,
distribute it to all of the processes, share the ghost rows with neighboring processes, 
and gather the data for output.
The functions are written to allow decomposition along each dimensions, for systems of arbitrary dimension.
The basic datatypes (int, float, double) are supported.
*/

class MPIGrid {

    private:
                                /// m prefix indicates membor variable
        int m_np;               /// number of processes
        int m_rank;             /// rank of process
        int m_ndims;            /// number of dimensions
        int m_nrows;            /// number of ghost rows
        int * m_local_dims;     /// dimensions of local share, with ghost rows
        int * m_global_dims;    /// dimensions of global system
        int * m_np_dims;        /// number of processes in each dimension

        MPI_Comm topology;      /// MPI Virtual Topology

        template<typename T>
        void pack(T const * const __restrict__ data, T * const __restrict__ packed_data, int count, int block_length, int stride);
        template<typename T>
        void unpack(T * const __restrict__ data, T const * const __restrict__ packed_data, int count, int block_length, int stride);
        template<typename T>
        MPI_Datatype getMPI_Datatype();

    public:

        MPIGrid();
        ~MPIGrid();

        int setup(MPI_Comm comm, int const * const global_dims, int const * const np_dims, int ndims, int nrows, int * const local_dims);

        template <typename T>
        int scatter(T const * const __restrict__ global_data, T * const __restrict__ local_data);

        template <typename T>
        int gather(T * const __restrict__ global_data, T const * const __restrict__ local_data);

        template <typename T>
        int share(T * const local_data);

};

MPIGrid :: MPIGrid()
{
}

MPIGrid :: ~MPIGrid()
{
    delete [] m_local_dims;
    delete [] m_global_dims;;
    delete [] m_np_dims;
}

template <typename T> 
MPI_Datatype MPIGrid :: getMPI_Datatype() {}
template <>
MPI_Datatype MPIGrid :: getMPI_Datatype<double>() { return MPI_DOUBLE; }
template <>
MPI_Datatype MPIGrid :: getMPI_Datatype<int>() { return MPI_INT; }
template <>
MPI_Datatype MPIGrid :: getMPI_Datatype<float>() { return MPI_FLOAT; }


template<typename T>
void MPIGrid :: pack(T const * const __restrict__ data, T * const __restrict__ packed_data, int count, int block_length, int stride)
{
    size_t num = block_length * sizeof(T);

    for (int i=0; i<count; i++)
    {
        void * source = (void *) (data + i*stride);
        void * destination = (void *) (packed_data + i*block_length);
        memcpy(destination, source, num);
    }
}

template<typename T>
void MPIGrid :: unpack(T * const __restrict__ data, T const * const __restrict__ packed_data, int count, int block_length, int stride)
{
    size_t num = block_length * sizeof(T);

    for (int i=0; i<count; i++)
    {
        void * source = (void *) (packed_data + i*block_length);
        void * destination = (void *) (data + i*stride);
        memcpy(destination, source, num);
    }

}

int MPIGrid :: setup(MPI_Comm comm_old, int const * const global_dims, int const * const np_dims, int ndims, int nrows, int * const local_dims)
{

    /// Setup a grid and store all the information needed for communication

    /**
    Setup must be called once (and only once) before any of the other class functions are called.
    Most of the error checking happens here, so the return error value should be checked for 0 (no errors);
    global_dims and local_dims must already be allocated with ndims elements
    */

    /**
    @param [in] comm_old the MPI communicator (usually MPI_COMM_WORLD)
    @param [in] global_dims extents of the global system
    @param [in] np_dims the number of processors in each dimensions (the decomposition)
    @param [in] ndims the number of dimensions
    @param [in] nrows the number of ghost rows needed (1 per laplacian)
    @param [out] local_dims dimensions of local data 
    */

    int periodic[ndims];
    int np_product;
    MPI_Comm_size(comm_old, &m_np);

    /** \error ERROR 1 the number of dimensions must be greater than 0 */
    if (ndims < 1) return 1;

    /** \error ERROR 2 the number of processors in each dimensions must be greater than 0 */
    for (int i=0; i<ndims; i++) 
        if (np_dims[i] < 1) return 2;

    /** \error ERROR 3 the global dimensions must be >= m_nrows */
    for (int i=0; i<ndims; i++) 
        if (global_dims[i] < nrows) return 3;

    /** \error ERROR 4 the number of processors must divide evenly in global dims in each dimension */
    for (int i=0; i<ndims; i++) 
        if (global_dims[i] % np_dims[i] != 0) return 4;

    /** \error ERROR 5 the number of processors in each dimension must equal the total number of processors */
    np_product = 1;
    for (int i=0; i<ndims; i++) 
        np_product *= np_dims[i];

    if (np_product != m_np) return 5;

    // create a cartesian topology
    for (int i=0; i<ndims; i++) periodic[i] = 1;

    MPI_Cart_create(comm_old, ndims, np_dims, periodic, 0, &topology);
    MPI_Comm_rank(topology, &m_rank);

    m_nrows = nrows;
    m_ndims = ndims;
    m_global_dims = new int [m_ndims];
    m_local_dims = new int [m_ndims];
    m_np_dims = new int [m_ndims];

    for (int i=0; i<ndims; i++)
    {
        m_global_dims[i] = global_dims[i];
        m_np_dims[i] = np_dims[i];
        m_local_dims[i] = m_global_dims[i] / m_np_dims[i] + m_nrows*2;
        local_dims[i] = m_local_dims[i];
    }

    return 0;
}

template <typename T>
int MPIGrid :: scatter(T const * const __restrict__ global_data, T * const __restrict__ local_data)
{

    /// Send the global data on the master process to the local data for each process

    /**
    All processes should call this function. 
    Only the master process needs to have allocated resources for global_data
    The contents of global_data will be scattered to local_data, so local_data does not need to be initialized.
    Supported data types are int, float, and double.
    */

    /**
    @param [in] global_data pointer to the first element of the global data on the master process
    @param [out] local_data pointer to the first element of the local data on each process
    */

    MPI_Request request;
    MPI_Status status;

    int source = 0;
    int tag = 1;
    int block_length;
    int stride;
    int count;
    int offset;
    
    int subdomain[m_ndims];
    int subdomain_volume;

    // calculate extents of subdomains
    subdomain_volume = 1;
    for (int i=0; i<m_ndims; i++)
    {
        subdomain[i] = m_global_dims[i] / m_np_dims[i];
        subdomain_volume *= subdomain[i];
    }

    // calculate number of contiguous chunks
    count = 1;
    for (int i=0; i<m_ndims-1; i++) count *= subdomain[i];

    T * packed_data = new T [subdomain_volume];

    /* ============== master sends data ============= */
    if (m_rank == 0) {

        int coord_stride[m_ndims];
        for (int i=0; i<m_ndims; i++)
        {
            coord_stride[i] = subdomain[i];
            for (int j=i+1; j<m_ndims; j++)
                coord_stride[i] *= m_np_dims[j]*subdomain[j];
        }

        block_length = subdomain[m_ndims-1];
        stride = m_global_dims[m_ndims-1];

        for (int id=0; id<m_np; id++) {

            int coords[m_ndims];
            MPI_Cart_coords(topology, id, m_ndims, coords);

            // calculate subdomain offset
            offset = 0;
            for (int i=0; i<m_ndims; i++) offset += coords[i] * coord_stride[i];

            pack(global_data+offset, packed_data, count, block_length, stride);
            MPI_Isend(packed_data, subdomain_volume, getMPI_Datatype<T>(), id, tag, topology, &request);
        }
    }

    /* ============== everyone receives data ============= */

    offset = 0;
    for (int i=0; i<m_ndims; i++)
    {
        int local_stride_i = 1;
        for (int j=i+1; j<m_ndims; j++)
            local_stride_i *= m_local_dims[j];
        offset += local_stride_i*m_nrows;
    }

    block_length = subdomain[m_ndims-1];
    stride = m_local_dims[m_ndims-1];

    MPI_Recv(packed_data, subdomain_volume, getMPI_Datatype<T>(), source, tag, topology, &status);
    unpack(local_data+offset, packed_data, count, block_length, stride);

    delete [] packed_data;

    return 0;
}

template<typename T>
int MPIGrid :: gather(T * const __restrict__ global_data, T const * const __restrict__ local_data)
{

    /// Collect the local data from each local process onto the master process

    /**
    All processes should call this function. 
    Only the master process needs to have allocated resources for global_data
    Supported data types are int, float, and double.
    */

    /**
    @param [out] global_data pointer to first element of global data on master process
    @param [in] local_data pointer to the first element of local data on each process.
    */

    MPI_Request request;
    MPI_Status status;

    int tag = 1;
    int destination = 0;
    int count;
    int block_length;
    int stride;
    int offset;

    int subdomain[m_ndims];
    int subdomain_volume;

    subdomain_volume = 1;
    for (int i=0; i<m_ndims; i++)
    {
        subdomain[i] = m_global_dims[i] / m_np_dims[i];
        subdomain_volume *= subdomain[i];
    }

    count = 1;
    for (int i=0; i<m_ndims-1; i++) count *= subdomain[i];

    block_length = subdomain[m_ndims-1];
    stride = m_local_dims[m_ndims-1];

    offset = 0;
    for (int i=0; i<m_ndims; i++)
    {
        int local_stride_i = 1;
        for (int j=i+1; j<m_ndims; j++)
            local_stride_i *= m_local_dims[j];
        offset += local_stride_i*m_nrows;
    }

    T * packed_send = new T [subdomain_volume];
    T * packed_recv = new T [subdomain_volume];

    /* ====================== EVERY PROCESS SENDS DATA ========================== */
    pack(local_data + offset, packed_send, count, block_length, stride);
    MPI_Isend(packed_send, subdomain_volume, getMPI_Datatype<T>(), destination, tag, topology, &request);

    /* ====================== MASTER RECEIVES DATA ========================== */
    if (m_rank == 0) {

        int coord_stride[m_ndims];

        for (int i=0; i<m_ndims; i++)
        {
            coord_stride[i] = subdomain[i];
            for (int j=i+1; j<m_ndims; j++)
                coord_stride[i] *= m_np_dims[j]*subdomain[j];
        }

        for (int id=0; id<m_np; id++) {
            int coords[m_ndims];
            MPI_Cart_coords(topology, id, m_ndims, coords);

            offset = 0;
            for (int i=0; i<m_ndims; i++) offset += coords[i] * coord_stride[i];

            stride = m_global_dims[m_ndims-1];

            MPI_Recv(packed_recv, subdomain_volume, getMPI_Datatype<T>(), id, tag, topology, &status);
            unpack(global_data+offset, packed_recv, count, block_length, stride);
        }

    }

    // must wait for master to finish recieving before free'ing packaged data
    MPI_Barrier(topology);

    delete [] packed_send;
    delete [] packed_recv;

    return 0;
}

template <typename T>
int MPIGrid :: share(T * const local_data)
{

    /// This function communicates ghost rows to neighboring processes
    
    /**
    Each process calls this function.
    All of the ghost row information is shared, the number of rows was indicated at setup.
    */

    /**
    @param [in,out] local_data pointer to the first element of the local data
    */

    for (int i=0; i<m_ndims; i++)
    {
        int tag = 1;
        int source, destination;
        MPI_Status status;

        int count = 1;
        int block_length = m_nrows;
        int stride = 1;
        int step = 1;

        int send_offset;
        int recv_offset;

        for (int j=0; j<i; j++)
            count *= m_local_dims[j];

        for (int j=i+1; j<m_ndims; j++)
            block_length *= m_local_dims[j];

        for (int j=i; j<m_ndims; j++)
            stride *= m_local_dims[j];

        for (int j=i+1; j<m_ndims; j++)
            step *= m_local_dims[j];

        T * packed_send = new T [count*block_length];
        T * packed_recv = new T [count*block_length];

        /* =========== SENDRECV POSITIVE DIRECTION =================== */

        send_offset = (m_local_dims[i] - 2*m_nrows) * step; 
        recv_offset = 0;

        pack(local_data+send_offset, packed_send, count, block_length, stride);

        MPI_Cart_shift(topology, i, 1, &source, &destination);

        MPI_Sendrecv(packed_send, count*block_length, getMPI_Datatype<T>(), destination, tag,
                     packed_recv, count*block_length, getMPI_Datatype<T>(), source, tag,
                     topology, &status);

        unpack(local_data+recv_offset, packed_recv, count, block_length, stride);

        /* =========== SENDRECV NEGATIVE DIRECTION =================== */

        send_offset = m_nrows * step; 
        recv_offset = (m_local_dims[i] - m_nrows) * step;

        pack(local_data+send_offset, packed_send, count, block_length, stride);

        MPI_Cart_shift(topology, i, -1, &source, &destination);

        MPI_Sendrecv(packed_send, count*block_length, getMPI_Datatype<T>(), destination, tag,
                     packed_recv, count*block_length, getMPI_Datatype<T>(), source, tag,
                     topology, &status);

        unpack(local_data+recv_offset, packed_recv, count, block_length, stride);

        delete [] packed_send;
        delete [] packed_recv;
    }

    return 0;

}
