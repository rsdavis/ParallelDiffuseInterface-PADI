
#include "model.hpp"

namespace model
{
    int phi;
    double dx;
    double dt;
    double a2;
    double a4;
    double w;
    double inv_dx_sq;
}

void preprocess(double ** phase, int * dims,
                std::map<std::string, std::string> params,
                std::map<std::string, int> name_index)
{

    // unpack phase index

    model::phi = name_index["phi"];

    // unpack model parameters

    model::dx = std::stod(params["dx"]);
    model::dt = std::stod(params["dt"]);

    model::a2 = std::stod(params["a2"]);
    model::a4 = std::stod(params["a4"]);
    model::w = std::stod(params["w"]);

    model::inv_dx_sq = 1.0/(model::dx*model::dx);

}

void integrate(double ** phase, double ** chem_pot, double ** mobility, int * dims)
{
    // int i, j, k; // reserved for looping, don't modify

    Derivatives derivatives;
    derivatives.setup(dims);

    for_loop_ijk(1)
    {
        int ndx = calc_ijk_index();

        double p = phase[model::phi][ndx];

        double laplace = derivatives.laplacian_h2(phase[model::phi], ndx, model::inv_dx_sq);

        chem_pot[model::phi][ndx] = model::a2 * p 
                                  + model::a4 * p*p*p 
                                  - model::w * laplace;

        mobility[model::phi][ndx] = 1.0;
    }

    for_loop_ijk(0)
    {
        int ndx = calc_ijk_index();
        phase[model::phi][ndx] += model::dt * derivatives.laplacian_h2(chem_pot[model::phi], ndx, model::inv_dx_sq);
    }
}


void postprocess(double ** phase, double ** chem_pot, double ** mobility, int * dims)
{
}

