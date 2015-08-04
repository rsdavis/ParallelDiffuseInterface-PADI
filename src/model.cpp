
#include "model.hpp"

namespace model
{

    /** 
    The model namespace includes variable definitions and data indices.
    It enables the parameters to be used 'globally' within this model file.
    These variables are referenced using model::variable.
    */

    int phi;
    double dx;
    double dt;
    double a2;
    double a4;
    double w;
}


/** model-specific functions can be included here */


void preprocess(double ** phase,  // order parameter data
                int * dims,       // system dimensions
                std::map<std::string, std::string> params, // input file parameters
                std::map<std::string, int> name_index)     // phase indices
{

    /**
    Preprocess is only called once, before the time-stepping begins.
    It is used primarily for unpacking parameters from the input file and
    getting the index to access each order parameter.
    It can also be used to allocate additional data storage and 
    calculate model-specific data that will be constant throughout the simulation.
    */

    // unpack phase index

    model::phi = name_index["phi"];

    // unpack model parameters

    unpack_parameter(model::dx, "dx", params);
    unpack_parameter(model::dt, "dt", params);

    unpack_parameter(model::a2, "a2", params);
    unpack_parameter(model::a4, "a4", params);
    unpack_parameter(model::w,   "w", params);

}

void kernel(double ** phase, double ** chem_pot, double ** mobility, int * dims)
{
    /**

    The kernel is run at every timestep and should include all calculations 
    required by the model in order to integrate the order parameters.

    All of the order parameter data is stored in the phase array.
    It is accessed using the phase index, and the ijk index ("ndx").

    Space is already allocated for storing chemical potential (chem_pot) and
    mobilities (mobility) values for each order parameter. 
    Their use is optional. They are only included for convenience.

    Looping is handled by using the for_loop_ijk(x) macro 
    where "x" is the number of ghost rows to be included in the loop.
    The calc_ijk_index macro calculates an "ndx" which indicates the position
    within the system at every iteration.
    These two macros are required in order to make the implementation independent of 
    system dimensionality.

    The stencil object includes common finite differencing operations. 
    A setup routine must be called in order to determine indexing values.

    Don't modify i,j,k variables, they are used for looping.

    */

    Stencil stencil;
    stencil.setup(dims, model::dx);

    for_loop_ijk(1)
    {
        int ndx = calc_ijk_index();

        double p = phase[model::phi][ndx];

        double laplace = stencil.laplacian_h2(phase[model::phi], ndx);

        chem_pot[model::phi][ndx] = model::a2 * p 
                                  + model::a4 * p*p*p 
                                  - model::w * laplace;

        mobility[model::phi][ndx] = 1.0;
    }

    for_loop_ijk(0)
    {
        int ndx = calc_ijk_index();
        phase[model::phi][ndx] += model::dt * stencil.laplacian_h2(chem_pot[model::phi], ndx);
    }
}


void postprocess(double ** phase, double ** chem_pot, double ** mobility, int * dims)
{

    /**
    This routine is run only once after the simulation has completed.
    It is primarily used for free'ing resources allocated in the preprocess routine.
    */
}

