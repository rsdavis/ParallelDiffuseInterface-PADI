
#ifndef SPF_DERIVATIVES_H
#define SPF_DERIVATIVES_H


class Derivatives {

    private:
    
        int X0Y0Z0;
        int XPY0Z0;
        int XMY0Z0;
        int X0YPZ0;
        int X0YMZ0;
        int XPYPZ0;
        int XPYMZ0;
        int XMYPZ0;
        int XMYMZ0;

        int X0Y0ZP;
        int XPY0ZP;
        int XMY0ZP;
        int X0YPZP;
        int X0YMZP;
        int XPYPZP;
        int XPYMZP;
        int XMYPZP;
        int XMYMZP;

        int X0Y0ZM;
        int XPY0ZM;
        int XMY0ZM;
        int X0YPZM;
        int X0YMZM;
        int XPYPZM;
        int XPYMZM;
        int XMYPZM;
        int XMYMZM;

    public:

        inline void setup(int * dims);
        inline double laplacian_h2(double * data, int ndx, double inv_dx_sq);

};

inline void Derivatives :: setup(int * dims)
{
    X0Y0Z0 = ( 0)*dims[1]*dims[2] + ( 0)*dims[2] + ( 0);
    XPY0Z0 = ( 1)*dims[1]*dims[2] + ( 0)*dims[2] + ( 0);
    XMY0Z0 = (-1)*dims[1]*dims[2] + ( 0)*dims[2] + ( 0);
    X0YPZ0 = ( 0)*dims[1]*dims[2] + ( 1)*dims[2] + ( 0);
    X0YMZ0 = ( 0)*dims[1]*dims[2] + (-1)*dims[2] + ( 0);
    XPYPZ0 = ( 1)*dims[1]*dims[2] + ( 1)*dims[2] + ( 0);
    XPYMZ0 = ( 1)*dims[1]*dims[2] + (-1)*dims[2] + ( 0);
    XMYPZ0 = (-1)*dims[1]*dims[2] + ( 1)*dims[2] + ( 0);
    XMYMZ0 = (-1)*dims[1]*dims[2] + (-1)*dims[2] + ( 0);

    X0Y0ZP = ( 0)*dims[1]*dims[2] + ( 0)*dims[2] + ( 1);
    XPY0ZP = ( 1)*dims[1]*dims[2] + ( 0)*dims[2] + ( 1);
    XMY0ZP = (-1)*dims[1]*dims[2] + ( 0)*dims[2] + ( 1);
    X0YPZP = ( 0)*dims[1]*dims[2] + ( 1)*dims[2] + ( 1);
    X0YMZP = ( 0)*dims[1]*dims[2] + (-1)*dims[2] + ( 1);
    XPYPZP = ( 1)*dims[1]*dims[2] + ( 1)*dims[2] + ( 1);
    XPYMZP = ( 1)*dims[1]*dims[2] + (-1)*dims[2] + ( 1);
    XMYPZP = (-1)*dims[1]*dims[2] + ( 1)*dims[2] + ( 1);
    XMYMZP = (-1)*dims[1]*dims[2] + (-1)*dims[2] + ( 1);

    X0Y0ZM = ( 0)*dims[1]*dims[2] + ( 0)*dims[2] + (-1);
    XPY0ZM = ( 1)*dims[1]*dims[2] + ( 0)*dims[2] + (-1);
    XMY0ZM = (-1)*dims[1]*dims[2] + ( 0)*dims[2] + (-1);
    X0YPZM = ( 0)*dims[1]*dims[2] + ( 1)*dims[2] + (-1);
    X0YMZM = ( 0)*dims[1]*dims[2] + (-1)*dims[2] + (-1);
    XPYPZM = ( 1)*dims[1]*dims[2] + ( 1)*dims[2] + (-1);
    XPYMZM = ( 1)*dims[1]*dims[2] + (-1)*dims[2] + (-1);
    XMYPZM = (-1)*dims[1]*dims[2] + ( 1)*dims[2] + (-1);
    XMYMZM = (-1)*dims[1]*dims[2] + (-1)*dims[2] + (-1);
}

#if SPF_NDIMS == 3
inline double Derivatives :: laplacian_h2(double * data, int ndx, double inv_dx_sq)
{
    double * ptr = data + ndx;
    double lap = *(ptr + XPY0Z0)
               + *(ptr + XMY0Z0)
               + *(ptr + X0YPZ0)
               + *(ptr + X0YMZ0)
               + *(ptr + X0Y0ZP)
               + *(ptr + X0Y0ZM)
               - *(ptr + X0Y0Z0)*6;
    return lap*inv_dx_sq;
}

#elif SPF_NDIMS == 2
inline double Derivatives :: laplacian_h2(double * data, int ndx, double inv_dx_sq)
{
    double * ptr = data + ndx;
    double lap = *(ptr + XPY0Z0)
               + *(ptr + XMY0Z0)
               + *(ptr + X0YPZ0)
               + *(ptr + X0YMZ0)
               - *(ptr + X0Y0Z0)*4;
    return lap*inv_dx_sq;
}

#endif

/*
inline double laplacian_h4(double * data, int ndx, double inv_dx_sq)
{
    double * ptr = data + ndx;
    double lap = 
               - *(ptr + X0Y0Z0)*24
               + *(ptr + XPY0Z0)*2
               + *(ptr + XMY0Z0)*2
               + *(ptr + X0YPZ0)*2
               + *(ptr + X0YMZ0)*2
               + *(ptr + XPYPZ0)
               + *(ptr + XPYMZ0)
               + *(ptr + XMYPZ0)
               + *(ptr + XMYMZ0)

               + *(ptr + X0Y0ZP)*2
               + *(ptr + XPY0ZP)
               + *(ptr + XMY0ZP)
               + *(ptr + X0YPZP)
               + *(ptr + X0YMZP)

               + *(ptr + X0Y0ZM)*2
               + *(ptr + XPY0ZM)
               + *(ptr + XMY0ZM)
               + *(ptr + X0YPZM)
               + *(ptr + X0YMZM);

    return lap*inv_dx_sq;
}


#if (SOFC_DIMENSIONALITY==3)
inline double div_M_grad_Mu(double **** mobility,
                            double **** chem_pot,
                            int ndx, int i, int j, int k,
                            double inv_dx_sq)
{

    double mx1 = mobility[ndx][i+1][j][k] + mobility[ndx][i][j][k];
    double mx2 = mobility[ndx][i][j][k] + mobility[ndx][i-1][j][k];

    double my1 = mobility[ndx][i][j+1][k] + mobility[ndx][i][j][k];
    double my2 = mobility[ndx][i][j][k] + mobility[ndx][i][j-1][k];

    double mz1 = mobility[ndx][i][j][k+1] + mobility[ndx][i][j][k];
    double mz2 = mobility[ndx][i][j][k] + mobility[ndx][i][j][k-1];

    double divx = mx1 * (chem_pot[ndx][i+1][j][k] - chem_pot[ndx][i][j][k]) -
                  mx2 * (chem_pot[ndx][i][j][k] - chem_pot[ndx][i-1][j][k]);

    double divy = my1 * (chem_pot[ndx][i][j+1][k] - chem_pot[ndx][i][j][k]) -
                  my2 * (chem_pot[ndx][i][j][k] - chem_pot[ndx][i][j-1][k]);

    double divz = mz1 * (chem_pot[ndx][i][j][k+1] - chem_pot[ndx][i][j][k]) -
                  mz2 * (chem_pot[ndx][i][j][k] - chem_pot[ndx][i][j][k-1]);

    double ddt = inv_dx_sq * 0.5 * (divx + divy + divz);

    return ddt;

}    

#elif (SOFC_DIMENSIONALITY==2)
inline double div_M_grad_Mu(double **** mobility,
                            double **** chem_pot,
                            int ndx, int i, int j, int k,
                            double inv_dx_sq)
{
    double mx1 = mobility[ndx][i+1][j][k] + mobility[ndx][i][j][k];
    double mx2 = mobility[ndx][i][j][k] + mobility[ndx][i-1][j][k];

    double my1 = mobility[ndx][i][j+1][k] + mobility[ndx][i][j][k];
    double my2 = mobility[ndx][i][j][k] + mobility[ndx][i][j-1][k];

    double divx = mx1 * (chem_pot[ndx][i+1][j][k] - chem_pot[ndx][i][j][k]) -
                  mx2 * (chem_pot[ndx][i][j][k] - chem_pot[ndx][i-1][j][k]);

    double divy = my1 * (chem_pot[ndx][i][j+1][k] - chem_pot[ndx][i][j][k]) -
                  my2 * (chem_pot[ndx][i][j][k] - chem_pot[ndx][i][j-1][k]);

    double ddt = inv_dx_sq * 0.5 * (divx + divy);

    return ddt;
}    

#endif
*/

#endif
