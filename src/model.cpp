
#include "model.hpp"

namespace model
{
    /** 
    The model namespace includes variable definitions and data indices.
    It enables the parameters to be used 'globally' within this model file.
    These variables are referenced using model::variable.
    */

    int a_ndx, b_ndx, c_ndx;

    double dx;
    double dt;

    double A;
    double B;
    double C;

    double Wa;
    double Wb;
    double Wc;

    double Ma;
    double Mb;
    double Mc;

}


/** model-specific functions can be included here */


void preprocess(double ** phase,  // order parameter data
                int * dims,       // system dimensions
                std::map<std::string, std::string> params, // input file parameters
                std::map<std::string, int> phase_index)     // phase indices
{

    /**
    Preprocess is only called once, before the time-stepping begins.
    It is used primarily for unpacking parameters from the input file and
    getting the index to access each order parameter.
    Preprocess can also be used to allocate additional data storage and 
    calculate model-specific data that will be constant throughout the simulation.
    */

    // unpack phase index

    unpack(phase_index, "a", model::a_ndx);
    unpack(phase_index, "b", model::b_ndx);
    unpack(phase_index, "c", model::c_ndx);

    // unpack model parameters

    unpack(params, "dx", model::dx);
    unpack(params, "dt", model::dt);

    unpack(params, "A", model::A);
    unpack(params, "B", model::B);
    unpack(params, "C", model::C);

    unpack(params, "Wa", model::Wa);
    unpack(params, "Wb", model::Wb);
    unpack(params, "Wc", model::Wc);

    unpack(params, "Ma", model::Ma);
    unpack(params, "Mb", model::Mb);
    unpack(params, "Mc", model::Mc);

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
    Any operations involving neighbor grid points should be implemented in Stencil.

    Don't modify i,j,k variables, they are used for looping.
    */

    Stencil stencil;
    stencil.setup(dims, model::dx);

    for_loop_ijk(1)
    {
        int ndx = calc_ijk_index();

        double a = phase[model::a_ndx][ndx];
        double b = phase[model::b_ndx][ndx];
        double c = phase[model::c_ndx][ndx];

        double asq = a*a;
        double bsq = b*b;
        double csq = c*c;

        double laplace_a = stencil.laplacian_h2(phase[model::a_ndx], ndx);
        double laplace_b = stencil.laplacian_h2(phase[model::b_ndx], ndx);
        double laplace_c = stencil.laplacian_h2(phase[model::c_ndx], ndx);

        chem_pot[model::a_ndx][ndx] = model::A * (2*a-6*asq+4*asq*a) 
                                    + 2*a*bsq + 2*a*csq
                                    - 2*(1-a)*(1-b)*(1-b)*(1-c)*(1-c)
                                    - model::Wa * laplace_a;

        chem_pot[model::b_ndx][ndx] = model::B * (2*b-6*bsq+4*bsq*b)
                                    + 2*asq*b + 2*b*csq
                                    - 2*(1-a)*(1-a)*(1-b)*(1-c)*(1-c)
                                    - model::Wb * laplace_b;


        chem_pot[model::c_ndx][ndx] = model::C * (2*c-6*csq+4*csq*c)
                                    + 2*asq*c + 2*bsq*c
                                    - 2*(1-a)*(1-a)*(1-b)*(1-b)*(1-c)
                                    - model::Wc * laplace_c;
    }

    for_loop_ijk(0)
    {
        int ndx = calc_ijk_index();
        phase[model::a_ndx][ndx] += model::dt*model::Ma*stencil.laplacian_h2(chem_pot[model::a_ndx],ndx);
        phase[model::b_ndx][ndx] += model::dt*model::Mb*stencil.laplacian_h2(chem_pot[model::b_ndx],ndx);
        phase[model::c_ndx][ndx] += model::dt*model::Mc*stencil.laplacian_h2(chem_pot[model::c_ndx],ndx);
    }
}


void postprocess(double ** phase, double ** chem_pot, double ** mobility, int * dims)
{

    /**
    This routine is run only once after the simulation has completed.
    It is primarily used for free'ing resources allocated in the preprocess routine.
    */
}

