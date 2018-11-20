module eos_bind_module
  use iso_c_binding
  use amrex_fort_module, only: amrex_real
  implicit none 


  !device interface
  interface 
    attributes(device) subroutine ckpy_d(rho, T, massfrac, iwrk, rwrk, p) &
      bind(C,name="CKPY")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: rho, T, massfrac(*), rwrk, p
      integer(c_int), intent(in) :: iwrk
    end subroutine ckpy_d

    attributes(device) subroutine ckums_d(T, iwrk, rwk, ei) &
      bind(C,name="CKUMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: ei(*)
      real(c_double), value :: T, rwk
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckums_d

    attributes(device) subroutine ckcpms_d(T, iwrk, rwk, cpi) &
      bind(C,name="CKCPMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, rwk
      real(c_double), intent(inout) :: cpi(*) 
      integer(c_int), intent(in) :: iwrk
    end subroutine ckcpms_d

    attributes(device) subroutine ckhms_d(T, iwrk, rwk, hi) & 
      bind(C,name="CKHMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, rwk
      real(c_double), intent(inout) :: hi(*)
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckhms_d

    attributes(device) subroutine get_T_given_ey_d(e, massfrac, iwrk, rwrk, T , lierr) &
      bind(C,name="GET_T_GIVEN_EY")
      use iso_c_binding, only: c_double, c_int 
      real(c_double), value :: e, rwrk, T
      real(c_double),  intent(inout) :: massfrac(*)
      integer(c_int), intent(in) :: iwrk, lierr
    end subroutine get_T_given_ey_d
   
    attributes(device) subroutine ckcvms_d(T, i, r, cvms) &
      bind(C, name="CKCVMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, r
      real(c_double), intent(inout) :: cvms(*)
      integer(c_int), intent(in) :: i 
    end subroutine ckcvms_d
  end interface  

  !host interface
  interface 

    attributes(host) subroutine ckpy(rho, T, massfrac, iwrk, rwrk, p) &
      bind(C,name="CKPY")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: rho, T, massfrac(*), rwrk, p
      integer(c_int), intent(in) :: iwrk
    end subroutine ckpy

    attributes(host) subroutine ckums(T, iwrk, rwk, ei) &
      bind(C,name="CKUMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: ei(*)
      real(c_double), value :: T, rwk
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckums

    attributes(host) subroutine ckcpms(T, iwrk, rwk, cpi) &
      bind(C,name="CKCPMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, rwk
      real(c_double), intent(inout) :: cpi(*) 
      integer(c_int), intent(in) :: iwrk
    end subroutine ckcpms

    attributes(host) subroutine ckhms(T, iwrk, rwk, hi) & 
      bind(C,name="CKHMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, rwk
      real(c_double), intent(inout) :: hi(*)
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckhms

    attributes(host) subroutine get_T_given_ey(e, massfrac, iwrk, rwrk, T , lierr) &
      bind(C,name="GET_T_GIVEN_EY")
      use iso_c_binding, only: c_double, c_int 
      real(c_double), value :: e, rwrk, T
      real(c_double),  intent(inout) :: massfrac(*)
      integer(c_int), intent(in) :: iwrk, lierr
    end subroutine get_T_given_ey
 
    attributes(host) subroutine ckcvms(T, i, r, cvms) &
      bind(C, name="CKCVMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, r
      real(c_double), intent(inout) :: cvms(*)
      integer(c_int), intent(in) :: i 
    end subroutine ckcvms
  end interface 
 
end module eos_bind_module
