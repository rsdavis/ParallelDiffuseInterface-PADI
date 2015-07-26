
#include <iostream>
#include <map>
#include "preprocessor.hpp"
#include "derivatives.h"

#if (SPF_NDIMS==1)
#define for_loop_ijk(x) for (int i=SPF_NROWS-x;i<dims[0]-SPF_NROWS+x;i++)

#elif (SPF_NDIMS==2)
#define for_loop_ijk(x) for (int i=SPF_NROWS-x;i<dims[0]-SPF_NROWS+x;i++) \
                        for (int j=SPF_NROWS-x;j<dims[1]-SPF_NROWS+x;j++)

#elif (SPF_NDIMS==3)
#define for_loop_ijk(x) for (int i=SPF_NROWS-x;i<dims[0]-SPF_NROWS+x;i++) \
                        for (int j=SPF_NROWS-x;j<dims[1]-SPF_NROWS+x;j++) \
                        for (int k=SPF_NROWS-x;k<dims[2]-SPF_NROWS+x;k++)
#endif

#if (SPF_NDIMS==1)
#define calc_ijk_index() i

#elif (SPF_NDIMS==2)
#define calc_ijk_index() i*dims[1] + j

#elif (SPF_NDIMS==3)
#define calc_ijk_index() i*dims[1]*dims[2] + j*dims[2] + k
#endif

void preprocess(double ** phase, int * dims,
                std::map<std::string, std::string> params,
                std::map<std::string, int> name_index);

void integrate(double ** phase, double ** chem_pot, double ** mobility, int * dims);

void postprocess(double ** phase, double ** chem_pot, double ** mobility, int * dims);
