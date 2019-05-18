#ifndef transport_H_
#define transport_H_

void actual_transport(bool wtr_get_xi, bool wtr_get_mu, bool wtr_get_lam, bool wtr_get_Ddiag, amrex::Real& Tloc,
              amrex::Real& rholoc, mrex::GpuArray<amrex::Real,NUM_SPECIES> const& Yloc,
              amrex::Real& mu, amrex::Real& xi,  amrex::Real& lam,
              amrex::Real& rho, amrex::GpuArray<amrex::Real,NUM_SPECIES> const& Ddiag);

 void comp_pure_bulk(amrex::Real Tloc,  amrex::GpuArray<amrex::Real,NUM_SPECIES> const& muloc,
          amrex::GpuArray<amrex::Real,NUM_SPECIES> const& xiloc );


#endif
