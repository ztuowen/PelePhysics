#include "mechanism.h"
#include "trans_params.H"
#include <AMReX_Arena.H>
#include <cstdio>

using namespace amrex;

namespace trans_params {

AMREX_GPU_DEVICE_MANAGED amrex::Real* wt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* iwt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eps;
AMREX_GPU_DEVICE_MANAGED amrex::Real* sig;
AMREX_GPU_DEVICE_MANAGED amrex::Real* dip;
AMREX_GPU_DEVICE_MANAGED amrex::Real* pol;
AMREX_GPU_DEVICE_MANAGED amrex::Real* zrot;
AMREX_GPU_DEVICE_MANAGED amrex::Real* fitmu;
AMREX_GPU_DEVICE_MANAGED amrex::Real* fitlam;
AMREX_GPU_DEVICE_MANAGED amrex::Real* fitdbin;
AMREX_GPU_DEVICE_MANAGED amrex::intg* nlin;
AMREX_GPU_DEVICE_MANAGED int array_size;
AMREX_GPU_DEVICE_MANAGED int fit_length;

void init ()
{
    array_size = NUM_SPECIES;
    wt   = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));
    iwt  = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));
    eps  = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));
    sig  = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));
    dip  = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));
    pol  = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));
    zrot = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));

    fitmu   = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size*fit_length));
    fitlam  = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size*fit_length));
    fitdbin = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size*array_size*fit_length));

    nlin = static_cast<int*>(The_Managed_Arena()->alloc(sizeof(int)*array_size));

    egtransetWT(wt);
    egtransetEPS(eps);
    egtransetSIG(sig);
    egtransetDIP(dip);
    egtransetPOL(pol);
    egtransetZROT(zrot);
    egtransetNLIN(nlin);
    egtransetCOFETA(fitmu);
    egtransetCOFLAM(fitlam);
    egtransetCOFD(fitdbin);

    for (int i=0; i < array_size; ++i){
        iwt[i] = 1. / wt[i];
    }
}

void finalize ()
{
    The_Managed_Arena()->free(wt);
    The_Managed_Arena()->free(iwt);
    The_Managed_Arena()->free(eps);
    The_Managed_Arena()->free(sig);
    The_Managed_Arena()->free(dip);
    The_Managed_Arena()->free(pol);
    The_Managed_Arena()->free(zrot);
    The_Managed_Arena()->free(fitmu);
    The_Managed_Arena()->free(fitlam);
    The_Managed_Arena()->free(fitdbin);
    The_Managed_Arena()->free(nlin);
}


}
