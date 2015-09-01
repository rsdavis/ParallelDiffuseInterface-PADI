
#include "model.hpp"

namespace model
{
    /** 
    The model namespace includes variable definitions and data indices.
    It enables the parameters to be used 'globally' within this model file.
    These variables are referenced using model::variable.
    */

    int phi;
    int eta0;
    int eta1;

    double dx;
    double dt;

    double m;
    double Wp;
    double Wn;

    double A;
    double B;
    double C;


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

    unpack(phase_index, "phi", model::phi);
    unpack(phase_index, "eta0", model::eta0);
    unpack(phase_index, "eta1", model::eta1);

    // unpack model parameters

    unpack(params, "dx", model::dx);
    unpack(params, "dt", model::dt);

    unpack(params, "Wp", model::Wp);
    unpack(params, "Wn", model::Wn);

    unpack(params, "A", model::A);
    unpack(params, "B", model::B);
    unpack(params, "C", model::C);

}

double P(double x)
{
    double xsq = x*x;
    return (2*x-6*xsq+4*xsq*x);
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

        double p = phase[model::phi][ndx];
        double n0 = phase[model::eta0][ndx];
        double n1 = phase[model::eta1][ndx];

        double psq = p*p;
        double n0sq = n0*n0;
        double n1sq = n1*n1;

        double laplace_phi = stencil.laplacian_h2(phase[model::phi], ndx);
        double laplace_eta0 = stencil.laplacian_h2(phase[model::eta0], ndx);
        double laplace_eta1 = stencil.laplacian_h2(phase[model::eta1], ndx);

        chem_pot[model::phi][ndx] = model::A*P(p) + model::C*(
                                  - 2*n0sq*(1-p)
                                  - 2*n1sq*(1-p)
                                  + 2*p*(1-n0)*(1-n0)*(1-n1)*(1-n1)
                                  + 2*p*n0sq*n1sq)
                                  - model::Wp * laplace_phi;

        chem_pot[model::eta0][ndx] = model::B*P(n0) + model::C*(
                                   + 2*n0*(1-p)*(1-p)
                                   - 2*psq*(1-n0)*(1-n1)*(1-n1)
                                   + 2*psq*n0*n1sq)
                                   - model::Wn * laplace_eta0;

        chem_pot[model::eta1][ndx] = model::B*P(n1) + model::C*(
                                   + 2*n1*(1-p)*(1-p)
                                   - 2*psq*(1-n0)*(1-n0)*(1-n1)
                                   + 2*psq*n0sq*n1)
                                   - model::Wn * laplace_eta1;
        
    }

    for_loop_ijk(0)
    {
        int ndx = calc_ijk_index();

        phase[model::phi][ndx] += model::dt * stencil.laplacian_h2(chem_pot[model::phi], ndx);
        phase[model::eta0][ndx] += -model::dt * chem_pot[model::eta0][ndx];
        phase[model::eta1][ndx] += -model::dt * chem_pot[model::eta1][ndx];

    }
}


void postprocess(double ** phase, double ** chem_pot, double ** mobility, int * dims)
{

    /**
    This routine is run only once after the simulation has completed.
    It is primarily used for free'ing resources allocated in the preprocess routine.
    */
}

