module eos_bind_module
  use iso_c_binding
  use amrex_fort_module, only: amrex_real
  implicit none 

#ifdef AMREX_USE_CUDA
  !device interface
  interface 
    AMREX_DEVICE subroutine ckpy_d(rho, T, massfrac, iwrk, rwrk, p) &
      bind(C,name="CKPY")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: rho, T, p
      real(c_double), intent(inout), dimension(9) :: massfrac
      real(c_double), intent(in) :: rwrk
      integer(c_int), intent(in) :: iwrk
    end subroutine ckpy_d

    AMREX_DEVICE subroutine ckums_d(T, iwrk, rwk, ei) &
      bind(C,name="CKUMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: ei(*)
      real(c_double), intent(in) :: T, rwk
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckums_d

    AMREX_DEVICE subroutine ckcpms_d(T, iwrk, rwk, cpi) &
      bind(C,name="CKCPMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(in) :: T, rwk
      real(c_double), intent(inout) :: cpi(*) 
      integer(c_int), intent(in) :: iwrk
    end subroutine ckcpms_d

    AMREX_DEVICE subroutine ckhms_d(T, iwrk, rwk, hi) & 
      bind(C,name="CKHMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(in) :: T, rwk
      real(c_double), intent(inout) :: hi(*)
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckhms_d

    AMREX_DEVICE subroutine get_T_given_ey_d(e, massfrac, iwrk, rwrk, T , lierr) &
      bind(C,name="GET_T_GIVEN_EY")
      use iso_c_binding, only: c_double, c_int 
      real(c_double), intent(in) :: e, rwrk
      real(c_double),  intent(inout) :: T,  massfrac(*)
      integer(c_int), intent(in) :: iwrk, lierr
    end subroutine get_T_given_ey_d
   
    AMREX_DEVICE subroutine ckcvms_d(T, i, r, cvms) &
      bind(C, name="CKCVMS_D")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(in) :: T, r
      real(c_double), intent(inout) :: cvms(*)
      integer(c_int), intent(in) :: i 
    end subroutine ckcvms_d
  end interface  
#endif
  !host interface
  interface 

    AMREX_CUDA_FORT_HOST subroutine ckpy(rho, T, massfrac, iwrk, rwrk, p) &
      bind(C,name="CKPY")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: rho, T, massfrac(*), p
      real(c_double), intent(in) :: rwrk
      integer(c_int), intent(in) :: iwrk
    end subroutine ckpy

    AMREX_CUDA_FORT_HOST subroutine ckums(T, iwrk, rwk, ei) &
      bind(C,name="CKUMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: ei(*)
      real(c_double), intent(in) :: T, rwk
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckums

    AMREX_CUDA_FORT_HOST subroutine ckcpms(T, iwrk, rwk, cpi) &
      bind(C,name="CKCPMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(in) :: T, rwk
      real(c_double), intent(inout) :: cpi(*) 
      integer(c_int), intent(in) :: iwrk
    end subroutine ckcpms

    AMREX_CUDA_FORT_HOST subroutine ckhms(T, iwrk, rwk, hi) & 
      bind(C,name="CKHMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(in) :: T, rwk
      real(c_double), intent(inout) :: hi(*)
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckhms

    AMREX_CUDA_FORT_HOST subroutine get_T_given_ey(e, massfrac, iwrk, rwrk, T , lierr) &
      bind(C,name="GET_T_GIVEN_EY")
      use iso_c_binding, only: c_double, c_int 
      real(c_double), intent(in)  :: e, rwrk
      real(c_double),  intent(inout) :: T, massfrac(*)
      integer(c_int), intent(in) :: iwrk
      integer(c_int), intent(out) :: lierr
    end subroutine get_T_given_ey
 
    AMREX_CUDA_FORT_HOST subroutine ckcvms(T, i, r, cvms) &
      bind(C, name="CKCVMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(in) :: T, r
      real(c_double), intent(inout) :: cvms(*)
      integer(c_int), intent(in) :: i 
    end subroutine ckcvms
  end interface 
 
end module eos_bind_module
