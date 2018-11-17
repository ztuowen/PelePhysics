module eos_bind_module
  use iso_c_binding
  use amrex_fort_module, only: amrex_real
  implicit none 



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
      real(c_double) :: cpi(*) 
      integer(c_int), value :: iwrk
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
    
  end interface  
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
      real(c_double) :: cpi(*) 
      integer(c_int), value :: iwrk
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
    
  end interface 
 
!  contains 
  
!   attributes(device) subroutine ckpy_d(rho, T, y, iwrk, rwrk, p) 
!      implicit none     
!      real(amrex_real), intent(out) :: p
!      real(amrex_real), intent(in) :: rwrk, rho, T, y(*) 
!      integer, intent(in) :: iwrk
!      real(amrex_real) :: YOW, imw(9)
!      imw = (/ 1.0 / 2.015940,1.0 / 31.998800, 1.0 / 18.015340, 1.0 / 1.007970, &
!               1.0 / 15.999400, 1.0 / 17.007370, 1.0 / 33.006770, 1.0 / 34.014740, &
!              1.0 / 28.013400 /)
!      YOW = 0 !/* for computing mean MW */
!      YOW = YOW + y(1)*imw(1)  ! /*H2 */
!      YOW = YOW + y(2)*imw(2) !/*O2 */
!      YOW = YOW + y(3)*imw(3) !/*H2O */
!      YOW = YOW + y(4)*imw(4) !/*H */
!      YOW = YOW + y(5)*imw(5) !/*O */
!      YOW = YOW + y(6)*imw(6) !/*OH */
!      YOW = YOW + y(7)*imw(7) !/*HO2 */
!      YOW = YOW + y(8)*imw(8) !/*H2O2 */
!      YOW = YOW + y(9)*imw(9) !/*N2 */
!      p = rho * 8.31451d+07 * T * YOW
!    end subroutine ckpy_d

end module eos_bind_module
