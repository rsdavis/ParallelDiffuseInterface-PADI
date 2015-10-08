
#ifndef SPF_DERIVATIVES_H
#define SPF_DERIVATIVES_H

#include <math.h>

class Stencil {

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

        double dx;
        double inv_dx;
        double inv_dx_sq;

    public:

        inline void setup(int * dims, double dx);
        inline SPF_DATATYPE laplacian_h2 (SPF_DATATYPE * data, int ndx);
        inline SPF_DATATYPE laplacian_h4 (SPF_DATATYPE * data, int ndx);
        inline SPF_DATATYPE grad_norm    (SPF_DATATYPE * data, int ndx);
        inline SPF_DATATYPE grad_sq      (SPF_DATATYPE * data, int ndx);
        inline SPF_DATATYPE div_A_grad_B (SPF_DATATYPE * A, SPF_DATATYPE * B, int ndx);

};

inline void Stencil :: setup(int * dims, double dx)
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

    Stencil::dx = dx;
    Stencil::inv_dx = 1.0/dx;
    Stencil::inv_dx_sq = 1.0/(dx*dx);
}

#if SPF_NDIMS == 3
inline SPF_DATATYPE Stencil :: laplacian_h2(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE lap = *(ptr + XPY0Z0)
               + *(ptr + XMY0Z0)
               + *(ptr + X0YPZ0)
               + *(ptr + X0YMZ0)
               + *(ptr + X0Y0ZP)
               + *(ptr + X0Y0ZM)
               - *(ptr + X0Y0Z0)*6;
    return lap*inv_dx_sq;
}

#elif SPF_NDIMS == 2
inline SPF_DATATYPE Stencil :: laplacian_h2(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE lap = *(ptr + XPY0Z0)
               + *(ptr + XMY0Z0)
               + *(ptr + X0YPZ0)
               + *(ptr + X0YMZ0)
               - *(ptr + X0Y0Z0)*4;
    return lap*inv_dx_sq;
}
#endif

#if SPF_NDIMS == 2
inline SPF_DATATYPE Stencil :: laplacian_h4(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE lap = *(ptr + XPY0Z0)*4
               + *(ptr + XMY0Z0)*4
               + *(ptr + X0YPZ0)*4
               + *(ptr + X0YMZ0)*4
               + *(ptr + XPYPZ0)
               + *(ptr + XMYPZ0)
               + *(ptr + XPYMZ0)
               + *(ptr + XMYMZ0)
               - *(ptr + X0Y0Z0)*20;
    return lap*inv_dx_sq/6.0;
}
#elif SPF_NDIMS == 3
inline SPF_DATATYPE Stencil :: laplacian_h4(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE lap = 
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
#endif

#if SPF_NDIMS == 2
inline SPF_DATATYPE Stencil :: grad_norm(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE grad_x = *(ptr + XPY0Z0) - *(ptr + XMY0Z0);
    SPF_DATATYPE grad_y = *(ptr + X0YPZ0) - *(ptr + X0YMZ0);

    return 0.5*inv_dx*sqrt(grad_x*grad_x + grad_y*grad_y);
}
#elif SPF_NDIMS == 3
inline SPF_DATATYPE Stencil :: grad_norm(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE grad_x = *(ptr + XPY0Z0) - *(ptr + XMY0Z0);
    SPF_DATATYPE grad_y = *(ptr + X0YPZ0) - *(ptr + X0YMZ0);
    SPF_DATATYPE grad_z = *(ptr + X0Y0ZP) - *(ptr + X0Y0ZM);

    return 0.5*inv_dx*sqrt(grad_x*grad_x + grad_y*grad_y + grad_z*grad_z);
}
#endif

#if SPF_NDIMS == 2
inline SPF_DATATYPE Stencil :: grad_sq(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE grad_x = *(ptr + XPY0Z0) - *(ptr + XMY0Z0);
    SPF_DATATYPE grad_y = *(ptr + X0YPZ0) - *(ptr + X0YMZ0);

    return 0.25*inv_dx_sq*(grad_x*grad_x + grad_y*grad_y);
}
#elif SPF_NDIMS == 3
inline SPF_DATATYPE Stencil :: grad_sq(SPF_DATATYPE * data, int ndx)
{
    SPF_DATATYPE * ptr = data + ndx;
    SPF_DATATYPE grad_x = *(ptr + XPY0Z0) - *(ptr + XMY0Z0);
    SPF_DATATYPE grad_y = *(ptr + X0YPZ0) - *(ptr + X0YMZ0);
    SPF_DATATYPE grad_z = *(ptr + X0Y0ZP) - *(ptr + X0Y0ZM);

    return 0.25*inv_dx_sq*(grad_x*grad_x + grad_y*grad_y + grad_z*grad_z);
}
#endif

#if SPF_NDIMS==3
inline SPF_DATATYPE Stencil :: div_A_grad_B(SPF_DATATYPE * A, SPF_DATATYPE * B, int ndx)
{

    SPF_DATATYPE * a = A + ndx;
    SPF_DATATYPE * b = B + ndx;

    SPF_DATATYPE axp = *(a+XPY0Z0) + *(a+X0Y0Z0);
    SPF_DATATYPE axm = *(a+X0Y0Z0) + *(a+XMY0Z0);

    SPF_DATATYPE ayp = *(a+X0YPZ0) + *(a+X0Y0Z0);
    SPF_DATATYPE aym = *(a+X0Y0Z0) + *(a+X0YMZ0);

    SPF_DATATYPE azp = *(a+X0Y0ZP) + *(a+X0Y0Z0);
    SPF_DATATYPE azm = *(a+X0Y0Z0) + *(a+X0Y0ZM);

    SPF_DATATYPE divx = axp * ( *(b+XPY0Z0) - *(b + X0Y0Z0) ) -
                  axm * ( *(b+X0Y0Z0) - *(b + XMY0Z0) );

    SPF_DATATYPE divy = ayp * ( *(b+X0YPZ0) - *(b + X0Y0Z0) ) -
                  aym * ( *(b+X0Y0Z0) - *(b + X0YMZ0) );

    SPF_DATATYPE divz = azp * ( *(b+X0Y0ZP) - *(b + X0Y0Z0) ) -
                  azm * ( *(b+X0Y0Z0) - *(b + X0Y0ZM) );

    SPF_DATATYPE ddt = inv_dx_sq * 0.5 * (divx + divy + divz);

    return ddt;

}    

#elif SPF_NDIMS==2
inline SPF_DATATYPE Stencil :: div_A_grad_B(SPF_DATATYPE * A, SPF_DATATYPE * B, int ndx)
{

    SPF_DATATYPE * a = A + ndx;
    SPF_DATATYPE * b = B + ndx;

    SPF_DATATYPE axp = *(a+XPY0Z0) + *(a+X0Y0Z0);
    SPF_DATATYPE axm = *(a+X0Y0Z0) + *(a+XMY0Z0);

    SPF_DATATYPE ayp = *(a+X0YPZ0) + *(a+X0Y0Z0);
    SPF_DATATYPE aym = *(a+X0Y0Z0) + *(a+X0YMZ0);

    SPF_DATATYPE divx = axp * ( *(b+XPY0Z0) - *(b + X0Y0Z0) ) -
                  axm * ( *(b+X0Y0Z0) - *(b + XMY0Z0) );

    SPF_DATATYPE divy = ayp * ( *(b+X0YPZ0) - *(b + X0Y0Z0) ) -
                  aym * ( *(b+X0Y0Z0) - *(b + X0YMZ0) );

    SPF_DATATYPE ddt = inv_dx_sq * 0.5 * (divx + divy);

    return ddt;

}    
#endif

#endif
