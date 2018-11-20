! Ths is a constant gamma equation of state
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
  attributes(device) subroutine eos_re_d(state)

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: lierr

    call eos_wb_d(state)

    call get_T_given_eY_d(state % e, state % massfrac, iwrk, rwrk, state % T, lierr)
!    if (lierr .ne. 0) then
!       print *, 'EOS: get_T_given_eY failed, T, e, Y = ', &
!            state % T, state % e, state % massfrac
!    end if
    state % T = max(state % T, smallT)
    call ckums_d(state % T, iwrk, rwrk, state % ei)
    call ckpy_d(state % rho, state % T, state % massfrac,iwrk, rwrk,  state % p)

    call eos_bottom_d(state)

  end subroutine eos_re_d

 attributes(device) subroutine eos_wb_d(state)

    implicit none

    type (eos_t), intent(inout) :: state

            state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  end subroutine eos_wb_d

  attributes(device)  subroutine eos_bottom_d(state)

    use amrex_constants_module
    use amrex_error_module
    implicit none

    type (eos_t), intent(inout) :: state
    real(amrex_real) :: Cvx

    call ckcvms_d(state % T, iwrk, rwrk, state % cvi)  ! erg/gi.K
    call ckcpms_d(state % T, iwrk, rwrk, state % cpi)  ! erg/gi.K
    call ckhms_d(state % T, iwrk, rwrk, state % hi)    ! erg/gi

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

  end subroutine eos_bottom_d

end module eos_module

