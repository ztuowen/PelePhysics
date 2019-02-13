module eos_module

  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use eos_type_module
  use eos_bind_module
  use chemistry_module, only : nspecies, inv_mwt, Ru, chemistry_init, chemistry_initialized, spec_names

  implicit none
  character (len=64) :: eos_name = "fuego"
  logical, save, private :: initialized = .false.

  real(amrex_real), public, parameter :: smallT = 1.d-50
  integer, parameter :: iwrk = 0
  real(amrex_real), parameter :: rwrk = 1.0d0

  public :: eos_init, eos_xty, eos_ytx, eos_ytx2, eos_ytx_vec, eos_cpi, eos_hi, eos_hi_vec, eos_cv, eos_cp, eos_p_wb, eos_wb, eos_get_activity, eos_rt, eos_tp, eos_rp, eos_re, eos_ps, eos_ph, eos_th, eos_rh, eos_get_transport, eos_h, eos_deriv, eos_mui
#ifdef AMREX_USE_CUDA
  public  :: eos_re_d
#endif

  interface
     subroutine amrex_array_init_snan (p, nelem) bind(C,name="amrex_array_init_snan")
       use iso_c_binding, only : c_double, c_size_t
       real(c_double),intent(inout) :: p
       integer (kind=c_size_t),intent(in),value :: nelem
     end subroutine amrex_array_init_snan
  end interface

contains

  subroutine eos_init()

    implicit none

    if (.not. chemistry_initialized)  call chemistry_init()

    initialized = .true.

  end subroutine eos_init

#ifdef AMREX_USE_CUDA
  AMREX_DEVICE  subroutine eos_bottom_d(state)

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
    state % dpde = (state % gam1 - 1.d0) * state % rho

    ! Try to avoid the expensive log function.  Since we don't need entropy
    ! in hydro solver, set it to an invalid but "nice" value for the plotfile.
    state % s = 1.d0

    ! Actually not sure what this is...used in composition derivatives
    ! in general system (not yet supported)
    state % dPdr = 0.d0

  end subroutine eos_bottom_d
#endif

  AMREX_CUDA_FORT_HOST  subroutine eos_bottom(state)

    use amrex_constants_module
    use amrex_error_module
    implicit none

    type (eos_t), intent(inout) :: state
    real(amrex_real) :: Cvx
    integer :: iwk
    real(amrex_real) :: rwk
    
    call ckcvms(state % T, iwk, rwk, state % cvi)  ! erg/gi.K
    call ckcpms(state % T, iwk, rwk, state % cpi)  ! erg/gi.K
    call ckhms(state % T, iwk, rwk, state % hi)    ! erg/gi
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

  subroutine eos_cpi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckcpms(state % T, iwrk, rwrk, state % cpi)

  end subroutine eos_cpi

  subroutine eos_hi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckhms(state % T, iwrk, rwrk, state % hi)

  end subroutine eos_hi

  subroutine eos_hi2(T, hi, Nsp)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout), dimension(1:Nsp) :: hi
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

  subroutine eos_p_wb(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)
    call ckpy(state % rho, state % T, state % massfrac, iwrk, rwrk, state % p)

  end subroutine eos_p_wb

#ifdef AMREX_USE_CUDA
 AMREX_DEVICE subroutine eos_wb_d(state)

    implicit none

    type (eos_t), intent(inout) :: state

    state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  end subroutine eos_wb_d
#endif

  AMREX_CUDA_FORT_HOST subroutine eos_wb(state)

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

  subroutine eos_rt(state)

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

  subroutine eos_rp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)

    state % T = state % p * state % wbar / (state % rho * Ru)
    call ckums(state % T, iwrk, rwrk, state % ei)
    state % e = sum(state % massfrac(:) * state % ei(:))

    call eos_bottom(state)

  end subroutine eos_rp


#ifdef AMREX_USE_CUDA
  AMREX_DEVICE subroutine eos_re_d(state)

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: lierr
    call eos_wb_d(state)

    call get_T_given_eY_d(state % e, state % massfrac, iwrk, rwrk, state % T, lierr)
    state % T = max(state % T, smallT)
    call ckums_d(state % T, iwrk, rwrk, state % ei)
    call ckpy_d(state % rho, state % T, state % massfrac, iwrk, rwrk, state % p)

    call eos_bottom_d(state)

  end subroutine eos_re_d
#endif

   subroutine eos_re(state)

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: lierr

    call eos_wb(state)

    call get_T_given_eY(state % e, state % massfrac, iwrk, rwrk, state % T, lierr)
    if (lierr .ne. 0) then
       print *, 'EOS: get_T_given_eY failed, T, e, Y = ', &
            state % T, state % e, state % massfrac
    end if
    state % T = max(state % T, smallT)
    call ckums(state % T, iwrk, rwrk, state % ei)
    call ckpy(state % rho, state % T, state % massfrac, iwrk, rwrk, state % p)

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

