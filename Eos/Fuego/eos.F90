! This is a constant gamma equation of state
!
! This a simplified version of the more general eos_gamma_general.
!

module eos_module

  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use eos_type_module
  use eos_bind_module
  use chemistry_module, only : nspecies, Ru, inv_mwt, chemistry_init, chemistry_initialized, spec_names, elem_names

  implicit none
  character (len=64) :: eos_name = "fuego"
  logical, save, private :: initialized = .false.

  real(amrex_real), save, managed, public :: smallT = 1.d-50
  integer, managed :: iwrk
  real(amrex_real), managed :: rwrk

  public :: eos_init, eos_xty, eos_ytx, eos_ytx2, eos_ytx_vec, eos_cpi, eos_hi, eos_hi_vec, eos_cv, eos_cp, eos_p_wb, eos_wb, eos_get_activity, eos_rt, eos_tp, eos_rp, eos_re, eos_ps, eos_ph, eos_th, eos_rh, eos_get_transport, eos_h, eos_deriv, eos_mui
  private :: nspecies, Ru, inv_mwt

  interface
     subroutine amrex_array_init_snan (p, nelem) bind(C,name="amrex_array_init_snan")
       use iso_c_binding, only : c_double, c_size_t
       real(c_double),intent(inout) :: p
       integer (kind=c_size_t),intent(in),value :: nelem
     end subroutine amrex_array_init_snan
  end interface

contains

  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module
    use parallel
    use iso_c_binding, only : c_double, c_size_t

    implicit none

    real(amrex_real), optional :: small_temp
    real(amrex_real), optional :: small_dens

    integer (kind=c_size_t) :: nelem

    nelem = 1
    call amrex_array_init_snan(mintemp,nelem)
    call amrex_array_init_snan(maxtemp,nelem)
    call amrex_array_init_snan(mindens,nelem)
    call amrex_array_init_snan(maxdens,nelem)
    call amrex_array_init_snan(minmassfrac,nelem)
    call amrex_array_init_snan(maxmassfrac,nelem)
    call amrex_array_init_snan(mine,nelem)
    call amrex_array_init_snan(maxe,nelem)
    call amrex_array_init_snan(minp,nelem)
    call amrex_array_init_snan(maxp,nelem)
    call amrex_array_init_snan(mins,nelem)
    call amrex_array_init_snan(maxs,nelem)
    call amrex_array_init_snan(minh,nelem)
    call amrex_array_init_snan(maxh,nelem)

    mintemp     = 1.d-200
    maxtemp     = 1.d200
    mindens     = 1.d-200
    maxdens     = 1.d200
    minmassfrac = 1.d-200
    maxmassfrac = 1.d0
    mine        = -1.d200
    maxe        = +1.d200
    minp        = 1.d-200
    maxp        = +1.d200
    mins        = -1.d200
    maxs        = +1.d200
    minh        = -1.d200
    maxh        = +1.d200

    if (.not. chemistry_initialized)  call chemistry_init()

    if (present(small_temp)) then
       if (small_temp < mintemp) then
          small_temp = mintemp
       else
          mintemp = small_temp
       endif
    endif

    if (present(small_dens)) then
       if (small_dens < mindens) then
          small_dens = mindens
       else
          mindens = small_dens
       endif
    endif

    initialized = .true.

  end subroutine eos_init

AMREX_CUDA_FORT_DEVICE  subroutine eos_bottom(state)

    use amrex_constants_module
    use amrex_error_module
    implicit none

    type (eos_t), intent(inout) :: state
    real(amrex_real) :: Cvx

    call ckcvms(state % T, iwrk, rwrk, state % cvi)  ! erg/gi.K
    call ckcpms(state % T, iwrk, rwrk, state % cpi)  ! erg/gi.K
    call ckhms (state % T, iwrk, rwrk, state % hi)    ! erg/gi

    state % cv = sum(state % massfrac(:) * state % cvi(:)) ! erg/g.K
    state % cp = sum(state % massfrac(:) * state % cpi(:)) ! erg/g.K
    state % h  = sum(state % massfrac(:) * state %  hi(:)) ! erg/g

    Cvx = state % wbar  *  state % cv ! erg/mole.K

    state % gam1 = (Cvx + Ru) / Cvx ! -
    state % cs = sqrt(state % gam1 * state % p / state % rho) ! cm/s

    state % dpdr_e = state % p / state % rho
    state % dpde = (state % gam1 - ONE) * state % rho

    ! Try to avoid the expensive log function.  Since we don't need entropy
    ! in hydro solver, set it to an invalid but "nice" value for the plotfile.
    state % s = ONE

    ! Actually not sure what this is...used in composition derivatives
    ! in general system (not yet supported)
    state % dPdr = ZERO

  end subroutine eos_bottom

  attributes(device,host)  subroutine ckcvms(T, i, r, cvms) 
    real(amrex_real), intent(inout) :: T, cvms(9)
    real(amrex_real), intent(in) :: r
    integer, intent(in) :: i 
    real(amrex_real) ::  tc(5)
    
    tc = (/ 0, T, T*T, T*T*T, T*T*T*T /) !temperature cache

    call cv_R(cvms, tc)
    !multiply by R/molecularweight 
    cvms(1) = cvms(1) * 4.124383662212169e+07 ! /*H2 */
    cvms(2) = cvms(2) * 2.598381814318037e+06 ! /*O2 */
    cvms(3) = cvms(3) * 4.615239012974499e+06 ! /*H2O */
    cvms(4) = cvms(4) * 8.248767324424338e+07 ! /*H */
    cvms(5) = cvms(5) * 5.196763628636074e+06 ! /*O */
    cvms(6) = cvms(6) * 4.888768810227566e+06 ! /*OH */
    cvms(7) = cvms(7) * 2.519031701678171e+06 ! /*HO2 */
    cvms(8) = cvms(8) * 2.444384405113783e+06 ! /*H2O2 */
    cvms(9) = cvms(9) * 2.968047434442088e+06 ! /*N2 */
  end subroutine ckcvms




! compute Cv/R at the given temperature */
! tc contains precomputed powers of T, tc[0] = log(T) */
attributes(device,host) subroutine cv_R(species,tc)
    implicit none 
    !)*temperature */
    real(amrex_real), intent(inout) :: species(9), tc(5)
    real(amrex_real) ::  T 

    T = tc(2)

    !)*species with midpoint at T=1000 kelvin */
    if (T < 1000) then 
        !)*species 0: H2 */
        species(1) =            +2.29812431e+00             +8.24944174e-04 * tc(2)            -8.14301529e-07 * tc(3)            -9.47543433e-11 * tc(4)            +4.13487224e-13 * tc(5)
        !)*species 1: O2 */
        species(2) =            +2.21293640e+00             +1.12748635e-03 * tc(2)            -5.75615047e-07 * tc(3)            +1.31387723e-09 * tc(4)            -8.76855392e-13 * tc(5)
        !)*species 2: H2O */
        species(3) =            +2.38684249e+00             +3.47498246e-03 * tc(2)            -6.35469633e-06 * tc(3)            +6.96858127e-09 * tc(4)            -2.50658847e-12 * tc(5)
        !)*species 3: H */
        species(4) =            +1.50000000e+00             +0.00000000e+00 * tc(2)            +0.00000000e+00 * tc(3)            +0.00000000e+00 * tc(4)            +0.00000000e+00 * tc(5)
        !)*species 4: O */
        species(5) =            +1.94642878e+00             -1.63816649e-03 * tc(2)            +2.42103170e-06 * tc(3)            -1.60284319e-09 * tc(4)            +3.89069636e-13 * tc(5)
        !)*species 5: OH */
        species(6) =            +3.12530561e+00             -3.22544939e-03 * tc(2)            +6.52764691e-06 * tc(3)            -5.79853643e-09 * tc(4)            +2.06237379e-12 * tc(5)
        !)*species 6: HO2 */
        species(7) =            +3.30179801e+00             -4.74912051e-03 * tc(2)            +2.11582891e-05 * tc(3)            -2.42763894e-08 * tc(4)            +9.29225124e-12 * tc(5)
        !)*species 7: H2O2 */
        species(8) =            +2.38875365e+00             +6.56922581e-03 * tc(2)            -1.48501258e-07 * tc(3)            -4.62580552e-09 * tc(4)            +2.47151475e-12 * tc(5)
        !)*species 8: N2 */
        species(9) =            +2.29867700e+00             +1.40824000e-03 * tc(2)            -3.96322200e-06 * tc(3)            +5.64151500e-09 * tc(4)            -2.44485500e-12 * tc(5)
     else
        !)*species 0: H2 */
        species(1) =            +1.99142337e+00             +7.00064411e-04 * tc(2)            -5.63382869e-08 * tc(3)            -9.23157818e-12 * tc(4)            +1.58275179e-15 * tc(5)
        !)*species 1: O2 */
        species(2) =            +2.69757819e+00             +6.13519689e-04 * tc(2)            -1.25884199e-07 * tc(3)            +1.77528148e-11 * tc(4)            -1.13643531e-15 * tc(5)
        !)*species 2: H2O */
        species(3) =            +1.67214561e+00             +3.05629289e-03 * tc(2)            -8.73026011e-07 * tc(3)            +1.20099639e-10 * tc(4)            -6.39161787e-15 * tc(5)
        !)*species 3: H */
        species(4) =            +1.50000000e+00             +0.00000000e+00 * tc(2)            +0.00000000e+00 * tc(3)            +0.00000000e+00 * tc(4)            +0.00000000e+00 * tc(5)
        !)*species 4: O */
        species(5) =            +1.54205966e+00             -2.75506191e-05 * tc(2)            -3.10280335e-09 * tc(3)            +4.55106742e-12 * tc(4)            -4.36805150e-16 * tc(5)
        !)*species 5: OH */
        species(6) =            +1.86472886e+00             +1.05650448e-03 * tc(2)            -2.59082758e-07 * tc(3)            +3.05218674e-11 * tc(4)            -1.33195876e-15 * tc(5)
        !)*species 6: HO2 */
        species(7) =            +3.01721090e+00             +2.23982013e-03 * tc(2)            -6.33658150e-07 * tc(3)            +1.14246370e-10 * tc(4)            -1.07908535e-14 * tc(5)
        !)*species 7: H2O2 */
        species(8) = +3.57316685e+00 +4.33613639e-03 * tc(2) -1.47468882e-06 * tc(3) +2.34890357e-10 * tc(4) -1.43165356e-14 * tc(5)
        !)*species 8: N2 */
        species(9) = +1.92664000e+00 + 1.48797700e-03 * tc(2) - 5.68476100e-07 * tc(3) +1.00970400e-10 * tc(4)-6.75335100e-15 * tc(5)
    endif 
  end subroutine cv_R


attributes(device,host)  subroutine cp_R(species,tc)
  implicit none 
    real(amrex_real), intent(inout) :: species(9), tc(5)
    !temperature 
    real(amrex_real) :: T 

    T = tc(2)

    ! species with midpoint at T=1000 kelvin */
    if (T < 1000) then  
        !species 1: H2 */
        species(1) =                &
            +3.29812431e+00         & 
            +8.24944174e-04 * tc(2) &
            -8.14301529e-07 * tc(3) &
            -9.47543433e-11 * tc(4) &
            +4.13487224e-13 * tc(5)
        ! species 2: O2 
        species(2) =                &
            +3.21293640e+00         &
            +1.12748635e-03 * tc(2) &
            -5.75615047e-07 * tc(3) &
            +1.31387723e-09 * tc(4) &
            -8.76855392e-13 * tc(5)
        ! species 3: H2O 
        species(3) =                &
            +3.38684249e+00         &
            +3.47498246e-03 * tc(2) &
            -6.35469633e-06 * tc(3) &
            +6.96858127e-09 * tc(4) &
            -2.50658847e-12 * tc(5) 
        ! species 4: H 
        species(4) =                &
            +2.50000000e+00         &
            +0.00000000e+00 * tc(2) &
            +0.00000000e+00 * tc(3) &
            +0.00000000e+00 * tc(4) &
            +0.00000000e+00 * tc(5)
        ! species 5: O 
        species(5) =                &
            +2.94642878e+00         &
            -1.63816649e-03 * tc(2) &
            +2.42103170e-06 * tc(3) &
            -1.60284319e-09 * tc(4) &
            +3.89069636e-13 * tc(5)
        ! species 6: OH 
        species(6) =                &
            +4.12530561e+00         &
            -3.22544939e-03 * tc(2) &
            +6.52764691e-06 * tc(3) &
            -5.79853643e-09 * tc(4) &
            +2.06237379e-12 * tc(5)
        ! species 7: HO2 
        species(7) =                &
            +4.30179801e+00         &
            -4.74912051e-03 * tc(2) &
            +2.11582891e-05 * tc(3) &
            -2.42763894e-08 * tc(4) &
            +9.29225124e-12 * tc(5)
        ! species 8: H2O2 
        species(8) =                &
            +3.38875365e+00         &
            +6.56922581e-03 * tc(2) &
            -1.48501258e-07 * tc(3) &
            -4.62580552e-09 * tc(4) &
            +2.47151475e-12 * tc(5)
        ! species 9: N2 
        species(9) =                &
            +3.29867700e+00         &
            +1.40824000e-03 * tc(2) &
            -3.96322200e-06 * tc(3) &
            +5.64151500e-09 * tc(4) &
            -2.44485500e-12 * tc(5)
    else 
        ! species 0: H2 
        species(1) =            +2.99142337e+00            +7.00064411e-04 * tc(2)            -5.63382869e-08 * tc(3)            -9.23157818e-12 * tc(4)            +1.58275179e-15 * tc(5)
        !)*species 1: O2 */
        species(2) =            +3.69757819e+00            +6.13519689e-04 * tc(2)            -1.25884199e-07 * tc(3)            +1.77528148e-11 * tc(4)            -1.13643531e-15 * tc(5)
        !)*species 2: H2O */
        species(3) =            +2.67214561e+00            +3.05629289e-03 * tc(2)            -8.73026011e-07 * tc(3)            +1.20099639e-10 * tc(4)            -6.39161787e-15 * tc(5)
        !)*species 3: H */
        species(4) =            +2.50000000e+00            +0.00000000e+00 * tc(2)            +0.00000000e+00 * tc(3)            +0.00000000e+00 * tc(4)            +0.00000000e+00 * tc(5)
        !)*species 4: O */
        species(5) =            +2.54205966e+00            -2.75506191e-05 * tc(2)            -3.10280335e-09 * tc(3)            +4.55106742e-12 * tc(4)            -4.36805150e-16 * tc(5)
        !)*species 5: OH */
        species(6) =            +2.86472886e+00            +1.05650448e-03 * tc(2)            -2.59082758e-07 * tc(3)            +3.05218674e-11 * tc(4)            -1.33195876e-15 * tc(5)
        !)*species 6: HO2 */
        species(7) =            +4.01721090e+00            +2.23982013e-03 * tc(2)            -6.33658150e-07 * tc(3)            +1.14246370e-10 * tc(4)            -1.07908535e-14 * tc(5)
        !)*species 7: H2O2 */
        species(8) =            +4.57316685e+00            +4.33613639e-03 * tc(2)            -1.47468882e-06 * tc(3)            +2.34890357e-10 * tc(4)            -1.43165356e-14 * tc(5)
        !)*species 8: N2 */
        species(9) =            +2.92664000e+00            +1.48797700e-03 * tc(2)            -5.68476100e-07 * tc(3)            +1.00970400e-10 * tc(4)            -6.75335100e-15 * tc(5)
    endif  
 end subroutine cp_R


  subroutine eos_xty(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckxty(state % molefrac,iwrk,rwrk,state % massfrac)

  end subroutine eos_xty

  subroutine eos_ytx(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckytx (state % massfrac,iwrk,rwrk,state % molefrac)

  end subroutine eos_ytx

  subroutine eos_ytx2(Y, X, Nsp)

    implicit none

    double precision, intent(in), dimension(1:Nsp) :: Y
    double precision, intent(out), dimension(1:Nsp) :: X
    integer, intent(in) :: Nsp

    call ckytx(Y(:),iwrk,rwrk,X(:))

  end subroutine eos_ytx2

  subroutine eos_ytx_vec(Y, ylo, yhi, X, xlo, xhi, lo, hi, Nsp)

    implicit none

    integer, intent(in) :: ylo(3), yhi(3)
    integer, intent(in) :: xlo(3), xhi(3)    
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: Nsp
    double precision, intent(in), dimension(ylo(1)-1:yhi(1)+1, ylo(2)-1:yhi(2)+1, ylo(3)-1:yhi(3)+1, 1:Nsp) :: Y
    double precision, intent(out), dimension(xlo(1)-1:xhi(1)+1, xlo(2)-1:xhi(2)+1, xlo(3)-1:xhi(3)+1, 1:Nsp) :: X

    integer :: j, k
    integer :: npts

    npts = (hi(1)+1)-(lo(1)-1)+1
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
         call VCKYTX( npts, Y(lo(1)-1:hi(1)+1, j, k, :), iwrk, rwrk, X( lo(1)-1:hi(1)+1, j, k, :) )
       enddo
    enddo

  end subroutine eos_ytx_vec

  attributes(device,host)  subroutine eos_cpi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckcpms(state % T, iwrk, rwrk, state % cpi)

  end subroutine eos_cpi

  attributes(device,host)  subroutine eos_hi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckhms(state % T, iwrk, rwrk, state % hi)

  end subroutine eos_hi

  attributes(device,host)  subroutine eos_hi2(T, hi, Nsp)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in), dimension(1:Nsp) :: hi
    integer, intent(in) :: Nsp

    call ckhms(T,iwrk,rwrk,hi(:))

  end subroutine eos_hi2

  subroutine eos_hi_vec(mass, masslo, masshi, T, Tlo, Thi, hi, hilo, hihi, low, high, Nsp)

    implicit none

    integer, intent(in) :: masslo(3), masshi(3)
    integer, intent(in) :: Tlo(3), Thi(3)
    integer, intent(in) :: hilo(3), hihi(3)    
    integer, intent(in) :: low(3), high(3)    
    integer, intent(in) :: Nsp

    double precision, intent(in), dimension(masslo(1)-1:masshi(1)+1, masslo(2)-1:masshi(2)+1, masslo(3)-1:masshi(3)+1, 1:Nsp) :: mass
    double precision, intent(in), dimension(Tlo(1)-1:Thi(1)+1, Tlo(2)-1:Thi(2)+1, Tlo(3)-1:Thi(3)+1 ) :: T
    double precision, intent(out), dimension(hilo(1)-1:hihi(1)+1, hilo(2)-1:hihi(2)+1, hilo(3)-1:hihi(3)+1, 1:Nsp) :: hi
    
    integer :: j, k
    integer :: npts

    npts = (high(1)+1)-(low(1)-1)+1
    do k = low(3)-1, high(3)+1
       do j = low(2)-1, high(2)+1
          call VCKHMS( npts, T(low(1)-1:high(1)+1, j, k), iwrk, rwrk, hi( low(1)-1:high(1)+1, j, k, :) )
       enddo
    enddo

  end subroutine eos_hi_vec

  subroutine eos_cv(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckcvbs(state % T, state % massfrac, iwrk, rwrk, state % cv)

  end subroutine eos_cv

  subroutine eos_cp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckcpbs(state % T, state % massfrac, iwrk, rwrk, state % cp)

  end subroutine eos_cp

AMREX_CUDA_FORT_DEVICE  subroutine eos_p_wb(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)
    call ckpy(state % rho, state % T, state % massfrac, iwrk, rwrk, state % p)

  end subroutine eos_p_wb

AMREX_CUDA_FORT_DEVICE  subroutine eos_wb(state)

    implicit none

    type (eos_t), intent(inout) :: state

            state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  end subroutine eos_wb

  subroutine eos_get_activity(state)

    implicit none

    type (eos_t), intent(inout) :: state

    double precision :: Cvx

    call ckytcr(state%rho, state % T, state % massfrac, iwrk, rwrk, state % Acti)
          
    call eos_wb(state)

    call eos_bottom(state)

  end subroutine eos_get_activity

 attributes(device,host)  subroutine eos_rt(state)

    implicit none
    type (eos_t), intent(inout) :: state

    call eos_wb(state)

    call ckpy(state % rho, state % T, state % massfrac, iwrk, rwrk, state % p)
    call ckums(state % T, iwrk, rwrk, state % ei)
    state % e = sum(state % massfrac(:) * state % ei(:))

    call eos_bottom(state)

  end subroutine eos_rt

  subroutine eos_tp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)

    call ckrhoy(state % p,state % T,state % massfrac,iwrk,rwrk,state % rho)
    call ckums(state % T, iwrk, rwrk, state % ei)
    state % e = sum(state % massfrac(:) * state % ei(:))

    call eos_bottom(state)

  end subroutine eos_tp

 attributes(device,host)  subroutine eos_rp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)

    state % T = state % p * state % wbar / (state % rho * Ru)
    call ckums(state % T, iwrk, rwrk, state % ei)
    state % e = sum(state % massfrac(:) * state % ei(:))

    call eos_bottom(state)

  end subroutine eos_rp

AMREX_CUDA_FORT_DEVICE subroutine eos_re(state)

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: lierr

    call eos_wb(state)

    call get_T_given_eY_d(state % e, state % massfrac, iwrk, rwrk, state % T, lierr)
    if (lierr .ne. 0) then
       print *, 'EOS: get_T_given_eY failed, T, e, Y = ', &
            state % T, state % e, state % massfrac
    end if
    state % T = max(state % T, smallT)
    call ckums_d(state % T, iwrk, rwrk, state % ei)
    call ckpy_d(state % rho, state % T, state % massfrac, iwrk, rwrk, state % p)

    call eos_bottom(state)

  end subroutine eos_re

  subroutine eos_ps(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_ps is not supported in this EOS.')

  end subroutine eos_ps

  subroutine eos_ph(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_ph is not supported in this EOS.')

  end subroutine eos_ph

  subroutine eos_th(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_th is not supported in this EOS.')

  end subroutine eos_th

  subroutine eos_rh(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_rh is not supported in this EOS.')

  end subroutine eos_rh

  subroutine eos_get_transport(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_get_transport is not supported in this EOS.')

  end subroutine eos_get_transport

  subroutine eos_h(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_h is not supported in this EOS.')

  end subroutine eos_h

  subroutine eos_deriv(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_deriv is not supported in this EOS.')

  end subroutine eos_deriv

  subroutine eos_mui(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_mui is not supported in this EOS.')

  end subroutine eos_mui

end module eos_module

