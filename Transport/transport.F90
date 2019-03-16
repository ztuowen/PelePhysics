module transport_module

  use egz_module
  use eos_type_module
  use transport_type_module
  use actual_transport_module

  implicit none

contains

  ! Transport initialization routine: read in general transport parameters, then 
  ! call any specific initialization used by the transport

  ! This subroutine should be called outside OMP PARALLEL
  subroutine transport_init() bind(C, name="transport_init")

    use extern_probin_module

    implicit none

    ! Set up any specific parameters or initialization steps required by the transport we are using.
    call actual_transport_init()

  end subroutine transport_init

  subroutine transport_close()

    use extern_probin_module

    implicit none

    ! Clean up any specific parameters or initialization steps required by the transport we are using.
    call actual_transport_close

  end subroutine transport_close


  subroutine transport(which, coeff)

    implicit none

    ! Input arguments

    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff

    ! TODO: May want to check/process input state
    
    ! Call the transport evaluator
    call actual_transport(which, coeff)

  end subroutine transport

  subroutine get_transport_coeffs( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       density,      r_lo,  r_hi, &
       D,            D_lo,  D_hi, &
       mu,          mu_lo, mu_hi, &
       xi,          xi_lo, xi_hi, &
       lam,        lam_lo,lam_hi) &
       bind(C, name="get_transport_coeffs")

    !use network, only: nspec
    !use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   r_lo(3),  r_hi(3)
    integer         , intent(in   ) ::   D_lo(3),  D_hi(3)
    integer         , intent(in   ) ::  mu_lo(3), mu_hi(3)
    integer         , intent(in   ) ::  xi_lo(3), xi_hi(3)
    integer         , intent(in   ) :: lam_lo(3),lam_hi(3)
    real (amrex_real), intent(in   ) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),nspec)
    real (amrex_real), intent(in   ) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real (amrex_real), intent(in   ) :: density(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real (amrex_real), intent(inout) :: D(D_lo(1):D_hi(1),D_lo(2):D_hi(2),D_lo(3):D_hi(3),nspec)
    real (amrex_real), intent(inout) :: mu(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))
    real (amrex_real), intent(inout) :: xi(xi_lo(1):xi_hi(1),xi_lo(2):xi_hi(2),xi_lo(3):xi_hi(3))
    real (amrex_real), intent(inout) :: lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3))

    double precision :: trv_eos_state_massfrac(nspec)
    double precision :: trv_eos_state_molefrac(nspec)
    double precision :: trv_eos_state_cpi(nspec)
    double precision :: trv_ddiag(nspec)
    double precision :: trv_eos_state_t
    double precision :: trv_eos_state_rho
    double precision :: trv_mu
    double precision :: trv_xi
    double precision :: trv_lam
    integer :: trv_npts
    double precision :: wt(nspec)
    double precision :: eps(nspec)
    double precision :: dip(nspec)
    double precision :: pol(nspec)
    double precision :: sig(nspec)
    double precision :: zrot(nspec)
    integer :: nlin(nspec)
    integer, parameter :: no=4
    double precision :: cfe(no,nspec)
    double precision :: cfl(no,nspec)
    double precision :: cfd(no,nspec,nspec)
    integer, parameter :: nfit=7
    double precision :: fita(nfit,nspec,nspec)
    double precision :: fita0(nfit) = (/ &
     .1106910525D+01, -.7065517161D-02, -.1671975393D-01, .1188708609D-01, &
     .7569367323D-03, -.1313998345D-02,  .1720853282D-03 /)
    double precision :: eps2(nspec,nspec)
    double precision :: xtr(nspec)
    double precision :: ytr(nspec)
    double precision :: aux(nspec)
    double precision :: cxi(nspec)
    double precision :: cint(nspec)
    double precision :: dlt(6)
    double precision :: beta(nspec)
    double precision :: eta(nspec)
    double precision :: etalg(nspec)
    double precision :: rn(nspec)
    double precision :: an(nspec)
    double precision :: zn(nspec)
    double precision :: dmi(nspec)
    double precision :: g(nspec,nspec)
    double precision :: bin(nspec,nspec)
    double precision :: a(nspec,nspec)

    !! local variables
    integer      :: i, j, k, n, np
    !type (wtr_t) :: which_trans
    !type (trv_t) :: coeff

    !np = hi(1)-lo(1)+1
    !call build(coeff,np)

    !which_trans % wtr_get_xi    = .true.
    !which_trans % wtr_get_mu    = .true.
    !which_trans % wtr_get_lam   = .true.
    !which_trans % wtr_get_Ddiag = .true.

    !do k = lo(3),hi(3)
    !   do j = lo(2),hi(2)

    !      do i=1,np
    !         coeff%eos_state(i)%massfrac(1:nspec) = massfrac(lo(1)+i-1,j,k,1:nspec)
    !      end do

    !      do i=lo(1),hi(1)
    !         coeff%eos_state(i-lo(1)+1)%T = temperature(i,j,k)
    !         coeff%eos_state(i-lo(1)+1)%rho = density(i,j,k)
    !      end do

    !      call transport(which_trans, coeff)

    !      do i=lo(1),hi(1)
    !         mu(i,j,k)  = coeff %  mu(i-lo(1)+1)
    !         xi(i,j,k)  = coeff %  xi(i-lo(1)+1)
    !         lam(i,j,k) = coeff % lam(i-lo(1)+1)
    !      end do
    !      do n=1,nspec
    !         do i=lo(1),hi(1)
    !            D(i,j,k,n) = coeff % Ddiag(i-lo(1)+1,n)
    !         end do
    !      end do

    !   end do
    !end do

    !call destroy(coeff)

    call egz_init_gpu(wt,eps,sig,dip,pol,zrot,nlin,cfe,cfl,cfd,fita,fita0,eps2)

    !$acc update device(nspec)
    !$acc enter data copyin(hi,lo,massfrac,temperature,density,mu,xi,lam,D) copyin(wt,eps,sig,dip,pol,zrot,nlin,cfe,cfl,cfd,fita,fita0,eps2)
    !$acc parallel loop gang vector collapse(3) private(k,j,i,n,trv_eos_state_t,trv_eos_state_rho,trv_eos_state_massfrac,trv_eos_state_molefrac,trv_eos_state_cpi,trv_ddiag,trv_mu,trv_xi,trv_lam,xtr,ytr,aux,cxi,cint,dlt,beta,eta,etalg,rn,an,zn,dmi,g,bin,a) present(hi,lo,massfrac,temperature,density,mu,xi,lam,D,wt,eps,sig,dip,pol,zrot,nlin,cfe,cfl,cfd,fita,fita0,eps2)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             do n=1,nspec
                trv_eos_state_massfrac(n) = massfrac(i,j,k,n)
             end do
             trv_eos_state_T = temperature(i,j,k)
             trv_eos_state_rho = density(i,j,k)

             call actual_transport_gpu(trv_eos_state_massfrac, trv_eos_state_molefrac, trv_eos_state_T, trv_eos_state_rho, trv_eos_state_cpi, trv_mu, trv_xi, trv_lam, trv_ddiag, no, nspec, nfit, wt, eps, zrot, nlin, cfe, cfl, cfd, fita, fita0, xtr, ytr, aux, cxi, cint, dlt, beta, eta, etalg, rn, an, zn, dmi, g, bin, a)

             mu(i,j,k)  = trv_mu
             xi(i,j,k)  = trv_xi
             lam(i,j,k) = trv_lam

             do n=1,nspec
                D(i,j,k,n) = trv_Ddiag(n)
             end do
          end do
       end do
    end do
    !$acc end parallel loop
    !$acc exit data copyout(mu,xi,lam,d) delete(hi,lo,massfrac,temperature,density,wt,eps,sig,dip,pol,zrot,nlin,cfe,cfl,cfd,fita,fita0,eps2)

  end subroutine get_transport_coeffs

  subroutine get_visco_coeffs( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       mu,          mu_lo, mu_hi) &
       bind(C, name="get_visco_coeffs")

    use network, only: nspec
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::  mu_lo(3), mu_hi(3)
    real (amrex_real), intent(in   ) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),nspec)
    real (amrex_real), intent(in   ) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real (amrex_real), intent(inout) :: mu(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))

    ! local variables
    integer      :: i, j, k, n, np
    type (wtr_t) :: which_trans
    type (trv_t) :: coeff

    np = hi(1)-lo(1)+1
    call build(coeff,np)

    which_trans % wtr_get_xi    = .false.
    which_trans % wtr_get_mu    = .true.
    which_trans % wtr_get_lam   = .false.
    which_trans % wtr_get_Ddiag = .false.

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)

          do i=1,np
             coeff%eos_state(i)%massfrac(1:nspec) = massfrac(lo(1)+i-1,j,k,1:nspec)
          end do

          do i=lo(1),hi(1)
             coeff%eos_state(i-lo(1)+1)%T = temperature(i,j,k)
          end do

          call transport(which_trans, coeff)

          do i=lo(1),hi(1)
             mu(i,j,k)  = coeff %  mu(i-lo(1)+1)
          end do

       end do
    end do

    call destroy(coeff)

  end subroutine get_visco_coeffs

end module transport_module
