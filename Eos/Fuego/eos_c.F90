module eos_bind_module
  use iso_c_binding
  use amrex_fort_module, only: amrex_real

  interface 

    attributes(device,host) subroutine ckpy_d(rho, T, massfrac, iwrk, rwrk, p) &
      bind(C,name="CKPY")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: rho, T, massfrac(:), rwrk, p
      integer(c_int), intent(in) :: iwrk
    end subroutine ckpy

    attributes(device,host) subroutine ckums_d(T, iwrk, rwk, ei) &
      bind(C,name="CKUMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), intent(inout) :: ei(:)
      real(c_double), value :: T, rwk
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckums

    attributes(device,host) subroutine ckcpms_d(T, iwrk, rwk, cpi) &
      bind(C,name="CKCPMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, rwk
      real(c_double) :: cpi(:) 
      integer(c_int), value :: iwrk
    end subroutine ckcpms

    attributes(device,host) subroutine ckhms_d(T, iwrk, rwk, hi) & 
      bind(C,name="CKHMS")
      use iso_c_binding, only: c_double, c_int
      real(c_double), value :: T, rwk
      real(c_double), intent(inout) :: hi(:)
      integer(c_int), intent(in) :: iwrk 
    end subroutine ckhms

    attributes(device,host) subroutine get_T_given_ey_d(e, massfrac, iwrk, rwrk, T , lierr) &
      bind(C,name="GET_T_GIVEN_EY")
      use iso_c_binding, only: c_double, c_int 
      real(c_double), value :: e, rwrk, T
      real(c_double),  intent(inout) :: massfrac(:)
      integer(c_int), intent(in) :: iwrk, lierr
    end subroutine get_T_given_ey
    
  end interface  

end module eos_bind_module
