#ifndef _Transport_F_H_
#define _Transport_F_H_
#include <AMReX_BLFort.H>

#include <AMReX.H>

extern "C"
{
  void get_transport_coeffs
    (const int &nc, const int *lo, const int *hi,
     const amrex::Real *q, const int *qlo, const int *qhi, 
     const amrex::Real *c, const int *clo, const int *chi);
}

#endif // _Transport_F_H_
