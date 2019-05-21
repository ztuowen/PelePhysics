#ifndef _actual_sdc_F_H_
#define _actual_sdc_F_H_ 

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

extern "C" {

  int reactor_init(const int* iE_in, const int* nLobato_in, const int* nsdcite_in, const int* iverbose_in);

  int react(amrex::Real *rY_in, amrex::Real *rY_src_in, 
		  amrex::Real *rX_in, amrex::Real *rX_src_in,
		  amrex::Real *P_in,  amrex::Real *dt_react, amrex::Real *time, int *Init);

}

#endif

#endif /*_PeleLM_F_H_*/
