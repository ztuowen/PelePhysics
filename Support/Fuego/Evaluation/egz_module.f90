module egz_module
  
  use fuego_chemistry

  implicit none

  private

  double precision, parameter :: Ru = 8.314d7
  double precision, parameter :: Patmos = 1.01325d6
  integer, parameter :: no_gpu = 4
  integer, parameter :: nfit = 7
  integer, parameter :: nspec_gpu = 9
  integer, parameter :: iflag_gpu = 5
  !$acc declare create(ru,patmos,iflag_gpu,no_gpu,nfit,nspec_gpu)

  logical, save :: use_bulk_visc = .true.
  integer, save :: iflag = -1
  integer, save :: np = -1
  integer, save :: ns, no
  !need to be declared, but are unused
  !$acc declare create(use_bulk_visc,iflag,np,ns,no)

  double precision, allocatable, save :: wt(:), iwt(:), eps(:), sig(:), dip(:), pol(:), zrot(:)
  integer, allocatable, save :: nlin(:) 
  double precision, allocatable, save :: cfe(:,:), cfl(:,:), cfd(:,:,:)
  double precision, allocatable, save :: eps2(:,:)
  double precision, allocatable, save :: fita(:,:,:)
  double precision, parameter :: fita0(nfit) = (/ &
       .1106910525D+01, -.7065517161D-02, -.1671975393D-01, .1188708609D-01, &
       .7569367323D-03, -.1313998345D-02,  .1720853282D-03 /)

  !
  ! dimension(np,ns)
  !
  double precision, allocatable, save :: xtr(:,:), ytr(:,:), aux(:,:)
  double precision, allocatable, save :: cxi(:,:), cint(:,:)
  !
  ! dimension(np)
  !
  double precision, allocatable, save :: sumtr(:), wwtr(:)
  !
  double precision, allocatable, save :: dlt(:,:)

  !
  ! dimension(np,ns)
  !
  double precision, allocatable, save :: beta(:,:), eta(:,:), etalg(:,:), &
       rn(:,:), an(:,:), zn(:,:), dmi(:,:)
  ! 
  ! dimension(np,ns,ns)
  !
  double precision, allocatable, save :: G(:,:,:), bin(:,:,:), A(:,:,:)

  !$omp threadprivate(xtr,ytr,aux,cxi,cint,sumtr,wwtr,dlt,beta,eta,etalg)
  !$omp threadprivate(rn,an,zn,dmi,G,bin,A,np)

  public :: iflag
  public :: egz_init, egz_close, EGZINI, EGZPAR, EGZE1, EGZE3, EGZK1, EGZK3, EGZL1, EGZVR1
  public :: egz_init_gpu, EGZPAR_gpu, EGZE3_gpu, EGZK3_gpu, EGZL1_gpu, EGZVR1_gpu
  ! egz_init and egz_close should be called outside OMP PARALLEL,
  ! whereas others are inside

contains

  ! This subroutine should be called outside OMP PARALLEL
  subroutine egz_init(use_bulk_visc_in)
    use fuego_module, only : egtransetKK, egtransetNO, egtransetWT, egtransetEPS, egtransetSIG, egtransetDIP, egtransetPOL, egtransetZROT, egtransetNLIN, egtransetCOFETA, egtransetCOFLAM, egtransetCOFD
    implicit none

    logical, intent(in) :: use_bulk_visc_in

    use_bulk_visc = use_bulk_visc_in

    if (use_bulk_visc) then
       iflag = 5
    else
       iflag = 3
    end if
    
    call egtransetKK(ns)
    call egtransetNO(no)
    
    allocate(wt(ns))
    allocate(iwt(ns))
    allocate(eps(ns))
    allocate(sig(ns))
    allocate(dip(ns))
    allocate(pol(ns))
    allocate(zrot(ns))
    allocate(nlin(ns))
       
    allocate(cfe(no,ns))
    allocate(cfl(no,ns))
    allocate(cfd(no,ns,ns))
    
    call egtransetWT(wt)
    iwt = 1.d0/wt
    call egtransetEPS(eps)
    call egtransetSIG(sig)
    call egtransetDIP(dip)
    call egtransetPOL(pol)
    call egtransetZROT(zrot)
    call egtransetNLIN(nlin)
    call egtransetCOFETA(cfe)
    call egtransetCOFLAM(cfl)
    call egtransetCOFD(cfd)
    
    allocate(eps2(ns,ns))
    allocate(fita(nfit,ns,ns))
    
    call LEVEPS()

    call EGZABC(fita,fita0)

  end subroutine egz_init

  ! This subroutine should be called outside OMP PARALLEL
  subroutine egz_close()
    if (allocated(wt)) deallocate(wt)
    if (allocated(iwt)) deallocate(iwt)
    if (allocated(eps)) deallocate(eps)
    if (allocated(sig)) deallocate(sig)
    if (allocated(dip)) deallocate(dip)
    if (allocated(pol)) deallocate(pol)
    if (allocated(zrot)) deallocate(zrot)
    if (allocated(nlin)) deallocate(nlin)
    if (allocated(cfe)) deallocate(cfe)
    if (allocated(cfl)) deallocate(cfl)
    if (allocated(cfd)) deallocate(cfd)
    if (allocated(eps2)) deallocate(eps2)    
    if (allocated(fita)) deallocate(fita)
    !$omp parallel
    call egz_close_np()
    !$omp end parallel
  end subroutine egz_close


  ! This subroutine can be called inside OMP PARALLEL
  subroutine EGZINI(np_in)
    implicit none
    integer, intent(in) :: np_in
    logical, save :: first_call = .true.
    !$omp threadprivate(first_call)

    if (first_call) then
       np = np_in
    end if

    if (first_call .or. np.ne.np_in) then

       call egz_close_np()

       np = np_in

       allocate(xtr(np,ns))
       allocate(ytr(np,ns))
       allocate(aux(np,ns))
       allocate(sumtr(np))
       allocate(wwtr(np))
       allocate(dlt(np,6))

       allocate(beta(np,ns))
       allocate(eta(np,ns))
       allocate(etalg(np,ns))
       allocate(rn(np,ns))
       allocate(an(np,ns))
       allocate(zn(np,ns))
       allocate(dmi(np,ns))

       if (iflag .gt. 1) then
          allocate(bin(np,ns,ns))
       end if
       if (iflag.eq.3 .or. iflag.eq.5) then
          allocate(G(np,ns,ns))
          allocate(A(np,ns,ns))
       end if

       if (iflag > 3) then
          allocate(cxi(np,ns))
          allocate(cint(np,ns))
       end if
       
    end if

    first_call = .false.

  end subroutine EGZINI


  ! This subroutine can be called inside OMP PARALLEL
  subroutine egz_close_np()
    if (allocated(xtr)) deallocate(xtr)
    if (allocated(ytr)) deallocate(ytr)
    if (allocated(aux)) deallocate(aux)
    if (allocated(cxi)) deallocate(cxi)
    if (allocated(cint)) deallocate(cint)
    if (allocated(sumtr)) deallocate(sumtr)
    if (allocated(wwtr)) deallocate(wwtr)
    if (allocated(dlt)) deallocate(dlt)
    if (allocated(beta)) deallocate(beta)
    if (allocated(eta)) deallocate(eta)
    if (allocated(etalg)) deallocate(etalg)
    if (allocated(rn)) deallocate(rn)
    if (allocated(an)) deallocate(an)
    if (allocated(zn)) deallocate(zn)
    if (allocated(dmi)) deallocate(dmi)
    if (allocated(G)) deallocate(G)
    if (allocated(bin)) deallocate(bin)
    if (allocated(A)) deallocate(A)
  end subroutine egz_close_np


  subroutine LEVEPS()
    implicit none
    double precision, parameter :: pi = 3.1415926535D0, &
         fac = 1.0D-12, dipmin = 1.0D-20, boltz = 1.38056D-16
    integer :: j, k
    double precision :: rooteps(ns)
    do j=1,ns
       rooteps(j) = sqrt(EPS(j))
    end do
    do j=1,ns
       !DEC$ IVDEP
       do k=1,j
          IF((DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN)) THEN
!-----------------------------------------------------------------------
!                K IS POLAR, J IS NONPOLAR
!-----------------------------------------------------------------------
             eps2(K,J) = 1.0D0 + 0.25D0*(POL(J)/SIG(J)**3) * &
                  (FAC/BOLTZ) * &
                  (DIP(K)**2/(EPS(K)*SIG(K)**3)) * &
                  rooteps(k)/rooteps(j)
          ELSE IF((DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN)) THEN
!-----------------------------------------------------------------------
!             J IS POLAR, K IS NONPOLAR
!-----------------------------------------------------------------------
             eps2(K,J) = 1.0D0 + 0.25D0*(POL(K)/SIG(K)**3) * &
                  (FAC/BOLTZ) * &
                  (DIP(J)**2/(EPS(J)*SIG(J)**3)) * &
                  rooteps(j)/rooteps(k)
          ELSE
!-----------------------------------------------------------------------
!              NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
!-----------------------------------------------------------------------
             eps2(K,J) = 1.d0
          ENDIF
          eps2(K,J) = log(rooteps(j)*rooteps(k)* eps2(K,J)*eps2(K,J))
       end do
    end do
    do j=1,ns
       do k=j+1,ns
          eps2(k,j) = eps2(j,k)
       end do
    end do
  end subroutine LEVEPS


  subroutine egzabc(FA, FA0)
    implicit none
    double precision, intent(in) :: FA0(nfit)
    double precision, intent(out) :: FA(nfit,ns,ns)
    integer i,j,k,l,m,mm
    double precision :: SUMA, prod
    DO J = 1, NS
       DO I = J, NS
          do m = 1, nfit
             SUMA = 0.0D0
             mm   = m - 1
             do k = mm, nfit-1
                prod = 1.0d0
                do l = 1, k-mm
                   prod = prod * (-eps2(i,j)) * dble(mm+l) / dble(l)
                enddo
                SUMA = SUMA + FA0(k+1) * PROD
             enddo
             FA(m,I,J) = SUMA
          enddo
       ENDDO
    ENDDO
    DO J = 1, NS
       DO I = 1, J-1
          do m = 1, nfit
             FA(m,I,J) = FA(m,J,I)
          enddo
       ENDDO
    ENDDO
  end subroutine egzabc


  subroutine EGZPAR_gpu(T, X, cpms, wt_gpu, eps_gpu, zrot_gpu, nlin_gpu, cfe_gpu, cfd_gpu, fita_gpu, xtr_gpu, ytr_gpu, aux_gpu, cxi_gpu, cint_gpu, dlt_gpu, eta_gpu, etalg_gpu, bin_gpu, A_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: T, X(nspec_gpu)
    double precision, intent(in) :: cpms(nspec_gpu)
    double precision, intent(in) :: wt_gpu(nspec_gpu)
    double precision, intent(in) :: eps_gpu(nspec_gpu)
    double precision, intent(in) :: zrot_gpu(nspec_gpu)
    integer, intent(in) :: nlin_gpu(nspec_gpu)
    double precision, intent(in) :: cfe_gpu(no_gpu,nspec_gpu)
    double precision, intent(in) :: cfd_gpu(no_gpu,nspec_gpu,nspec_gpu)
    double precision, intent(in) :: fita_gpu(nfit,nspec_gpu,nspec_gpu)
    double precision, intent(out) :: xtr_gpu(nspec_gpu)
    double precision, intent(out) :: ytr_gpu(nspec_gpu)
    double precision, intent(out) :: aux_gpu(nspec_gpu)
    double precision, intent(out) :: cxi_gpu(nspec_gpu)
    double precision, intent(out) :: cint_gpu(nspec_gpu)
    double precision, intent(out) :: dlt_gpu(6)
    double precision, intent(out) :: eta_gpu(nspec_gpu)
    double precision, intent(out) :: etalg_gpu(nspec_gpu)
    double precision, intent(out) :: bin_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(out) :: A_gpu(nspec_gpu,nspec_gpu)

    integer :: n
    double precision :: aaa, sumtr, wwtr
    double precision, parameter :: sss = 1.d-16

    call LZPAR_gpu(T, cpms, wt_gpu, eps_gpu, zrot_gpu, nlin_gpu, cfe_gpu, cfd_gpu, fita_gpu, cxi_gpu, cint_gpu, dlt_gpu, eta_gpu, etalg_gpu, bin_gpu, A_gpu)
!-----------------------------------------------------------------------
!     Add a small constant to the mole and mass fractions
!-----------------------------------------------------------------------

      sumtr = 0.d0
      do n=1,nspec_gpu
         sumtr = sumtr + X(n)
      enddo

      aaa = sumtr / dble(nspec_gpu)

      wwtr = 0.0d0
      do n=1,nspec_gpu
         xtr_gpu(n) = X(n) + sss*(aaa - X(n))
         wwtr = wwtr + xtr_gpu(n) * wt_gpu(n)
      end do
!-----------------------------------------------------------------------
!     AUX(i) = \sum_{j .ne. i} YTR(j)
!-----------------------------------------------------------------------

      aaa = 0.d0
      do n=1,nspec_gpu
         ytr_gpu(n) = xtr_gpu(n) * wt_gpu(n) / wwtr
         aaa = aaa + ytr_gpu(n)
      end do
      do n=1,nspec_gpu
         aux_gpu(n) = aaa - ytr_gpu(n)
      end do

  end subroutine EGZPAR_gpu

  subroutine LZPAR_gpu(T, cpms, wt_gpu, eps_gpu, zrot_gpu, nlin_gpu, cfe_gpu, cfd_gpu, fita_gpu, cxi_gpu, cint_gpu, dlt_gpu, eta_gpu, etalg_gpu, bin_gpu, A_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: T
    double precision, intent(in) :: cpms(nspec_gpu)
    double precision, intent(in) :: wt_gpu(nspec_gpu)
    double precision, intent(in) :: eps_gpu(nspec_gpu)
    double precision, intent(in) :: zrot_gpu(nspec_gpu)
    integer, intent(in) :: nlin_gpu(nspec_gpu)
    double precision, intent(in) :: cfe_gpu(no_gpu,nspec_gpu)
    double precision, intent(in) :: cfd_gpu(no_gpu,nspec_gpu,nspec_gpu)
    double precision, intent(in) :: fita_gpu(nfit,nspec_gpu,nspec_gpu)
    double precision, intent(out) :: cxi_gpu(nspec_gpu)
    double precision, intent(out) :: cint_gpu(nspec_gpu)
    double precision, intent(out) :: dlt_gpu(6)
    double precision, intent(out) :: eta_gpu(nspec_gpu)
    double precision, intent(out) :: etalg_gpu(nspec_gpu)
    double precision, intent(out) :: bin_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(out) :: A_gpu(nspec_gpu,nspec_gpu)

    integer :: m, n
    double precision :: tmp, crot
    double precision :: wru, dr, sqdr, dr32, aaaa1, dd, sqdd, dd32, bbbb
    double precision, parameter :: PI1=1.d0/3.1415926535D0, PI32O2=2.7842D+00, &
         P2O4P2=4.4674D+00, PI32=5.5683D+00

    dlt_gpu(1) = log(T)
    dlt_gpu(2) = dlt_gpu(1) * dlt_gpu(1)
    dlt_gpu(3) = dlt_gpu(2) * dlt_gpu(1)
    dlt_gpu(4) = dlt_gpu(3) * dlt_gpu(1)
    dlt_gpu(5) = dlt_gpu(4) * dlt_gpu(1)
    dlt_gpu(6) = dlt_gpu(5) * dlt_gpu(1)

    do n=1,nspec_gpu
       etalg_gpu(n) = cfe_gpu(1,n) + cfe_gpu(2,n)*dlt_gpu(1) + cfe_gpu(3,n)*dlt_gpu(2) + cfe_gpu(4,n)*dlt_gpu(3)
       eta_gpu(n) = exp(etalg_gpu(n))
    end do

    !if (iflag_gpu .le. 1) return

    do n=1,nspec_gpu
       do m=1,n-1
          tmp = -(cfd_gpu(1,m,n)+cfd_gpu(2,m,n)*dlt_gpu(1)+cfd_gpu(3,m,n)*dlt_gpu(2) &
               + cfd_gpu(4,m,n)*dlt_gpu(3))
          bin_gpu(m,n) = exp(tmp)
          bin_gpu(n,m) = bin_gpu(m,n)
       end do
       bin_gpu(n,n) = 0.d0
    end do

    !if (iflag_gpu .le. 2) return

    !if (iflag_gpu.eq.3 .or. iflag_gpu.eq.5) then
       do n=1,nspec_gpu
          do m=1,n-1
             A_gpu(m,n) = fita_gpu(1,m,n) + fita_gpu(2,m,n)*dlt_gpu(1) + fita_gpu(3,m,n)*dlt_gpu(2) &
                  + fita_gpu(4,m,n)*dlt_gpu(3) + fita_gpu(5,m,n)*dlt_gpu(4) &
                  + fita_gpu(6,m,n)*dlt_gpu(5) + fita_gpu(7,m,n)*dlt_gpu(6)
             A_gpu(n,m) = A_gpu(m,n)
          end do
          A_gpu(n,n) = fita_gpu(1,n,n) + fita_gpu(2,n,n)*dlt_gpu(1) + fita_gpu(3,n,n)*dlt_gpu(2) &
               + fita_gpu(4,n,n)*dlt_gpu(3) + fita_gpu(5,n,n)*dlt_gpu(4) &
               + fita_gpu(6,n,n)*dlt_gpu(5) + fita_gpu(7,n,n)*dlt_gpu(6)
       end do
    !end if

    !if (iflag_gpu .eq. 3) return
!-----------------------------------------------------------------------
!         COMPUTE PARKER CORRECTION FOR ZROT
!         AND ALSO THE ROTATIONAL AND INTERNAL PARTS OF SPECIFIC HEAT
!-----------------------------------------------------------------------

    do n=1,nspec_gpu
       select case(nlin_gpu(n))
       case (0)
          crot = 0.d0
          cint_gpu(n) = 0.d0
       case (1)
          wru = wt_gpu(n) / Ru
          crot = 1.d0
          cint_gpu(n) = cpms(n) * wru - 2.50d0
       case (2)
          wru = wt_gpu(n) / Ru
          crot = 1.5d0
          cint_gpu(n) = cpms(n) * wru - 2.50d0
       !case default
       !   print *, "EFZ: wrong value in nlin"
       !   stop
       end select
       
       dr = eps_gpu(n) / 298.d0
       sqdr = sqrt(dr)
       dr32 = sqdr*dr
       aaaa1 = 1.d0/((1.0d0 + PI32O2*sqdr + P2O4P2*dr + PI32*dr32) * max(1.0d0, zrot_gpu(n)))

       dd = eps_gpu(n) / T
       sqdd = sqrt(dd)
       dd32 = sqdd*dd
       bbbb = (1.0d0 + PI32O2*sqdd + P2O4P2*dd + PI32*dd32) 
       cxi_gpu(n) = crot * PI1 * bbbb * aaaa1
    end do
  end subroutine LZPAR_gpu

  subroutine EGZE3_gpu(T, mu, wt_gpu, xtr_gpu, beta_gpu, eta_gpu, rn_gpu, an_gpu, zn_gpu, dmi_gpu, G_gpu, bin_gpu, A_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: T
    double precision, intent(out) :: mu
    double precision, intent(in) :: wt_gpu(nspec_gpu)
    double precision, intent(in) :: xtr_gpu(nspec_gpu)
    double precision, intent(inout) :: beta_gpu(nspec_gpu)
    double precision, intent(in) :: eta_gpu(nspec_gpu)
    double precision, intent(inout) :: rn_gpu(nspec_gpu)
    double precision, intent(inout) :: an_gpu(nspec_gpu)
    double precision, intent(inout) :: zn_gpu(nspec_gpu)
    double precision, intent(inout) :: dmi_gpu(nspec_gpu)
    double precision, intent(inout) :: G_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(inout) :: bin_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(inout) :: A_gpu(nspec_gpu,nspec_gpu)

    integer :: n

    call EGZEMH_gpu(T, wt_gpu, xtr_gpu, beta_gpu, eta_gpu, G_gpu, bin_gpu, A_gpu)
    
    rn_gpu = beta_gpu

    call EGZCG1_gpu(1, rn_gpu, an_gpu, zn_gpu, dmi_gpu, G_gpu)

    mu = 0.d0
    do n=1,nspec_gpu
       mu = mu + an_gpu(n) * beta_gpu(n)
    end do
  end subroutine EGZE3_gpu

  subroutine EGZEMH_gpu(T, wt_gpu, xtr_gpu, beta_gpu, eta_gpu, G_gpu, bin_gpu, A_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: T
    double precision, intent(in) :: wt_gpu(nspec_gpu)
    double precision, intent(in) :: xtr_gpu(nspec_gpu)
    double precision, intent(inout) :: beta_gpu(nspec_gpu)
    double precision, intent(in) :: eta_gpu(nspec_gpu)
    double precision, intent(inout) :: G_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(in) :: bin_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(in) :: A_gpu(nspec_gpu,nspec_gpu)

    integer ::  m, n
    double precision :: FAC, CCC, wtfac, wtmn, wtnm, aaa

    ! EVALUATE THE MATRIX H
    ! Note:    FAC * BIN = 2 W_{ij} / \eta_{ij} / A_{ij}
    FAC = (6.0D0 * RU / ( 5.0D0 * PATMOS )) * T
    CCC = 5.0D0 / 3.0D0

    do n=1,nspec_gpu
       G_gpu(n,n) = xtr_gpu(n) * xtr_gpu(n) / eta_gpu(n)
       ! EVALUATE THE RHS BETA
       beta_gpu(n) = xtr_gpu(n)
    end do

    do n=1,nspec_gpu
       do m=1,n-1
          wtfac = 1.d0/(wt_gpu(m) + wt_gpu(n))
          wtmn = wt_gpu(m)*(1.d0/wt_gpu(n))
          wtnm = wt_gpu(n)*(1.d0/wt_gpu(m))
          aaa = bin_gpu(m,n) * xtr_gpu(n) * xtr_gpu(m) * FAC * wtfac
          G_gpu(m,m) = G_gpu(m,m) + aaa*(A_gpu(m,n)*wtnm + CCC)
          G_gpu(m,n) = aaa * (A_gpu(m,n) - CCC)
          G_gpu(n,m) = G_gpu(m,n)
          G_gpu(n,n) = G_gpu(n,n) + aaa*(A_gpu(m,n)*wtmn + CCC)
       end do
    end do

  end subroutine EGZEMH_gpu


  subroutine EGZCG1_gpu(itmax, rn_gpu, an_gpu, zn_gpu, dmi_gpu, G_gpu)
    !$acc routine seq
    implicit none
    integer, intent(in) :: itmax
    double precision, intent(inout) :: rn_gpu(nspec_gpu)
    double precision, intent(inout) :: an_gpu(nspec_gpu)
    double precision, intent(inout) :: zn_gpu(nspec_gpu)
    double precision, intent(inout) :: dmi_gpu(nspec_gpu)
    double precision, intent(in) :: G_gpu(nspec_gpu,nspec_gpu)

    integer :: niter,  n
    double precision :: betan, aaa, bbb, ccc, temp(nspec_gpu)

    aaa = 0.d0
    betan = 0.d0

    do n = 1, nspec_gpu
       an_gpu(n) = 0.0d0
       zn_gpu(n) = 0.0d0
       dmi_gpu(n) = 1.0D0 / G_gpu(n,n)
       aaa = aaa + dmi_gpu(n) * rn_gpu(n)*rn_gpu(n)
    enddo

    do niter=1, itmax
       do n=1, nspec_gpu
          zn_gpu(n) = dmi_gpu(n)*rn_gpu(n) + betan*zn_gpu(n)
       end do

       CALL EGZAXS_gpu(G_gpu, zn_gpu, temp)

       bbb = 0.d0
       do n=1,nspec_gpu
          bbb = bbb + zn_gpu(n) * temp(n)
       end do

       do n=1,nspec_gpu
          an_gpu(n) = an_gpu(n) + aaa/bbb*zn_gpu(n)
          rn_gpu(n) = rn_gpu(n) - aaa/bbb*temp(n)
       end do

       if (niter .eq. itmax) exit

       ccc = 0.d0
       do n=1,nspec_gpu
          ccc = ccc + dmi_gpu(n) * rn_gpu(n)*rn_gpu(n)
       end do

       betan = ccc / aaa
       aaa = ccc

    end do

  end subroutine EGZCG1_gpu

  subroutine EGZAXS_gpu(AA, X, B) ! B = AA.X, AA is symmetric.
    !$acc routine seq
    implicit none
    double precision, intent(in) :: AA(nspec_gpu,nspec_gpu)
    double precision, intent(in) :: X(nspec_gpu)
    double precision, intent(out) :: B(nspec_gpu)
    integer ::  m, n
    B = 0.d0
    do n=1,nspec_gpu
       B(n) = 0.d0
       do m=1,nspec_gpu
          B(n) = B(n) + AA(m,n) * X(m)
       end do
    end do
  end subroutine EGZAXS_gpu


  subroutine EGZK3_gpu(T, VV, wt_gpu, xtr_gpu, cxi_gpu, cint_gpu, beta_gpu, eta_gpu, G_gpu, bin_gpu, A_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: T
    double precision, intent(out) :: VV
    double precision, intent(in) :: wt_gpu(nspec_gpu)
    double precision, intent(in) :: xtr_gpu(nspec_gpu)
    double precision, intent(in) :: cxi_gpu(nspec_gpu)
    double precision, intent(in) :: cint_gpu(nspec_gpu)
    double precision, intent(inout) :: beta_gpu(nspec_gpu)
    double precision, intent(in) :: eta_gpu(nspec_gpu)
    double precision, intent(inout) :: G_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(in) :: bin_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(in) :: A_gpu(nspec_gpu,nspec_gpu)
    
    integer ::  m, n
    double precision :: ccc, wtfac, bb
    double precision, parameter :: denfac = 2.4d0 * Ru / Patmos

    !if (.not. use_bulk_visc) then
    !   VV = 0.d0
    !   return
    !end if

    ccc = 0.0d0
    do n=1,nspec_gpu
       ccc = ccc + xtr_gpu(n)*cint_gpu(n)
    end do

    do n=1,nspec_gpu
       G_gpu(n,n) = 4.d0*cxi_gpu(n)/eta_gpu(n)*xtr_gpu(n)*xtr_gpu(n)
       beta_gpu(n) = -xtr_gpu(n) * cint_gpu(n) / (ccc+1.5d0)          
    end do

    do n=1,nspec_gpu
       do m=1,n-1
          wtfac = (1.d0/wt_gpu(m)) + (1.d0/wt_gpu(n))
          bb = xtr_gpu(m)*xtr_gpu(n)*bin_gpu(m,n)*denfac*T*A_gpu(m,n)*wtfac
          G_gpu(m,m) = G_gpu(m,m) + bb*cxi_gpu(m)
!          G_gpu(n,m) = 0.d0
!          G_gpu(m,n) = 0.d0
          G_gpu(n,n) = G_gpu(n,n) + bb*cxi_gpu(n)
       end do
    end do

    VV = 0.d0
    do n=1,nspec_gpu
       if (cxi_gpu(n) .eq. 0.d0) then
          VV = VV + beta_gpu(n) * beta_gpu(n) 
       else
          VV = VV + beta_gpu(n) * beta_gpu(n) / G_gpu(n,n)
       end if
    end do
  end subroutine EGZK3_gpu


  subroutine EGZL1_gpu(alpha, X, con, cfl_gpu, dlt_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: alpha
    double precision, intent(in) :: X(nspec_gpu)
    double precision, intent(out) :: con
    double precision, intent(in) :: cfl_gpu(no_gpu,nspec_gpu)
    double precision, intent(in) :: dlt_gpu(6)

    integer ::  n
    double precision :: asum, alpha1

    asum = 0.d0
    if (alpha .eq. 0.d0) then
       do n=1,nspec_gpu
          asum = asum + X(n)*(cfl_gpu(1,n) + cfl_gpu(2,n)*dlt_gpu(1) &
               + cfl_gpu(3,n)*dlt_gpu(2) + cfl_gpu(4,n)*dlt_gpu(3))
       end do
       con = exp(asum) 
    else
       alpha1 = 1.d0 / alpha
       do n=1,nspec_gpu
          asum = asum + X(n)*exp(alpha*(cfl_gpu(1,n) + cfl_gpu(2,n)*dlt_gpu(1) &
               + cfl_gpu(3,n)*dlt_gpu(2) + cfl_gpu(4,n)*dlt_gpu(3)))
       end do
       if (alpha .eq. 1.d0) then
          con = asum
       else if (alpha .eq. -1.d0) then
          con = 1.d0/asum
       else
          con = asum**alpha1
       end if
    end if
    
  end subroutine EGZL1_gpu
 

  subroutine EGZVR1_gpu(T, D, wt_gpu, xtr_gpu, aux_gpu, bin_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: T
    double precision, intent(out) :: D(nspec_gpu)
    double precision, intent(in) :: wt_gpu(nspec_gpu)
    double precision, intent(in) :: xtr_gpu(nspec_gpu)
    double precision, intent(in) :: aux_gpu(nspec_gpu)
    double precision, intent(in) :: bin_gpu(nspec_gpu,nspec_gpu)

    integer ::  m, n
    double precision :: fac

    D = 0.d0
    do n=1,nspec_gpu
       do m=1,nspec_gpu
          D(m) = D(m) + xtr_gpu(n)*bin_gpu(m,n)
       end do
    end do

    fac = (Patmos/Ru) / T

    do n=1,nspec_gpu
       D(n) = wt_gpu(n) * fac * aux_gpu(n) / D(n)
    end do

  end subroutine EGZVR1_gpu

  subroutine egzabc_gpu(FA, FA0, eps2_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: FA0(nfit)
    double precision, intent(out) :: FA(nfit,nspec_gpu,nspec_gpu)
    double precision, intent(in) :: eps2_gpu(nspec_gpu,nspec_gpu)
    integer i,j,k,l,m,mm
    double precision :: SUMA, prod
    DO J = 1, NS
       DO I = J, NS
          do m = 1, nfit
             SUMA = 0.0D0
             mm   = m - 1
             do k = mm, nfit-1
                prod = 1.0d0
                do l = 1, k-mm
                   prod = prod * (-eps2_gpu(i,j)) * dble(mm+l) / dble(l)
                enddo
                SUMA = SUMA + FA0(k+1) * PROD
             enddo
             FA(m,I,J) = SUMA
          enddo
       ENDDO
    ENDDO
    DO J = 1, NS
       DO I = 1, J-1
          do m = 1, nfit
             FA(m,I,J) = FA(m,J,I)
          enddo
       ENDDO
    ENDDO
  end subroutine egzabc_gpu


  subroutine LEVEPS_gpu(eps_gpu, eps2_gpu, dip_gpu, pol_gpu, sig_gpu)
    !$acc routine seq
    implicit none
    double precision, intent(in) :: eps_gpu(nspec_gpu)
    double precision, intent(inout) :: eps2_gpu(nspec_gpu,nspec_gpu)
    double precision, intent(in) :: dip_gpu(nspec_gpu)
    double precision, intent(in) :: pol_gpu(nspec_gpu)
    double precision, intent(in) :: sig_gpu(nspec_gpu)

    double precision, parameter :: pi = 3.1415926535D0, &
         fac = 1.0D-12, dipmin = 1.0D-20, boltz = 1.38056D-16
    integer :: j, k
    double precision :: rooteps(nspec_gpu)

    do j=1,nspec_gpu
       rooteps(j) = sqrt(EPS_gpu(j))
    end do
    do j=1,nspec_gpu
       do k=1,j
          IF((DIP_gpu(J).LT.DIPMIN .AND. DIP_gpu(K).GT.DIPMIN)) THEN
!-----------------------------------------------------------------------
!                K IS POLAR, J IS NONPOLAR
!-----------------------------------------------------------------------
             eps2_gpu(K,J) = 1.0D0 + 0.25D0*(POL_gpu(J)/SIG_gpu(J)**3) * &
                  (FAC/BOLTZ) * &
                  (DIP_gpu(K)**2/(EPS_gpu(K)*SIG_gpu(K)**3)) * &
                  rooteps(k)/rooteps(j)
          ELSE IF((DIP_gpu(J).GT.DIPMIN .AND. DIP_gpu(K).LT.DIPMIN)) THEN
!-----------------------------------------------------------------------
!             J IS POLAR, K IS NONPOLAR
!-----------------------------------------------------------------------
             eps2_gpu(K,J) = 1.0D0 + 0.25D0*(POL_gpu(K)/SIG_gpu(K)**3) * &
                  (FAC/BOLTZ) * &
                  (DIP_gpu(J)**2/(EPS_gpu(J)*SIG_gpu(J)**3)) * &
                  rooteps(j)/rooteps(k)
          ELSE
!-----------------------------------------------------------------------
!              NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
!-----------------------------------------------------------------------
             eps2_gpu(K,J) = 1.d0
          ENDIF
          eps2_gpu(K,J) = log(rooteps(j)*rooteps(k)* eps2_gpu(K,J)*eps2_gpu(K,J))
       end do
    end do
    do j=1,nspec_gpu
       do k=j+1,nspec_gpu
          eps2_gpu(k,j) = eps2_gpu(j,k)
       end do
    end do
  end subroutine LEVEPS_gpu


  subroutine egz_init_gpu(wt_gpu,eps_gpu,sig_gpu,dip_gpu,pol_gpu,zrot_gpu,nlin_gpu,cfe_gpu,cfl_gpu,cfd_gpu,fita_gpu,fita0_gpu,eps2_gpu)
    !$acc routine(egtransetWT,egtransetEPS, egtransetZROT, egtransetNLIN, egtransetCOFETA, egtransetCOFLAM, egtransetCOFD, egtransetDIP, egtransetSIG, egtransetPOL) seq
    USE fuego_module, ONLY: egtransetWT, egtransetEPS, egtransetZROT, egtransetNLIN, egtransetCOFETA, egtransetCOFLAM, egtransetCOFD, egtransetDIP, egtransetSIG, egtransetPOL
    implicit none
    double precision, intent(out) :: wt_gpu(nspec_gpu)
    double precision, intent(out) :: eps_gpu(nspec_gpu)
    double precision, intent(out) :: sig_gpu(nspec_gpu)
    double precision, intent(out) :: dip_gpu(nspec_gpu)
    double precision, intent(out) :: pol_gpu(nspec_gpu)
    double precision, intent(out) :: zrot_gpu(nspec_gpu)
    integer, intent(out) :: nlin_gpu(nspec_gpu)
    double precision, intent(out) :: cfe_gpu(no_gpu,nspec_gpu)
    double precision, intent(out) :: cfl_gpu(no_gpu,nspec_gpu)
    double precision, intent(out) :: cfd_gpu(no_gpu,nspec_gpu,nspec_gpu)
    double precision, intent(out) :: fita_gpu(nfit,nspec_gpu,nspec_gpu)
    double precision, intent(in) :: fita0_gpu(nfit)
    double precision, intent(out) :: eps2_gpu(nspec_gpu,nspec_gpu)
    !$acc enter data create(wt_gpu,eps_gpu,sig_gpu,dip_gpu,pol_gpu,zrot_gpu,nlin_gpu,cfe_gpu,cfl_gpu,cfd_gpu,eps2_gpu)
    !$acc parallel
    call egtransetWT(wt_gpu)
    call egtransetEPS(eps_gpu)
    call egtransetSIG(sig_gpu)
    call egtransetDIP(dip_gpu)
    call egtransetPOL(pol_gpu)
    call egtransetZROT(zrot_gpu)
    call egtransetNLIN(nlin_gpu)
    call egtransetCOFETA(cfe_gpu)
    call egtransetCOFLAM(cfl_gpu)
    call egtransetCOFD(cfd_gpu)
    call levEPS_gpu(eps_gpu,eps2_gpu,dip_gpu,pol_gpu,sig_gpu)
    !$acc end parallel
    !$acc update host(eps2_gpu)
    call egzABC_gpu(fita_gpu,fita0_gpu,eps2_gpu)
  end subroutine egz_init_gpu

  ! This subroutine can be called inside OMP PARALLEL
  subroutine EGZPAR(T, X, cpms)
    implicit none
    double precision, intent(in) :: T(np), X(np,ns)
    double precision, intent(in), optional :: cpms(np,ns)

    integer :: i, n
    double precision :: aaa(np)
    double precision, parameter :: sss = 1.d-16

    call LZPAR(T, cpms)

!-----------------------------------------------------------------------
!     Add a small constant to the mole and mass fractions
!-----------------------------------------------------------------------
      sumtr = 0.d0
      do n=1,ns
         do i=1,np
            sumtr(i) = sumtr(i) + X(i,n)
         end do
      enddo

      do i=1,np
         aaa(i)  = sumtr(i) / dble(ns)
      end do

      wwtr = 0.0d0
      do n=1,ns
         do i=1,np
            xtr(i,n) = X(i,n) + sss*(aaa(i) - X(i,n))
            wwtr(i) = wwtr(i) + xtr(i,n) * wt(n)
         end do
      end do

!-----------------------------------------------------------------------
!     AUX(i) = \sum_{j .ne. i} YTR(j)
!-----------------------------------------------------------------------
      aaa = 0.d0
      do n=1,ns
         do i=1,np
            ytr(i,n) = xtr(i,n) * wt(n) / wwtr(i)
            aaa(i) = aaa(i) + ytr(i,n)
         end do
      end do
      do n=1,ns
         do i=1,np
            aux(i,n) = aaa(i) - ytr(i,n)
         end do
      end do

  end subroutine EGZPAR


  subroutine LZPAR(T, cpms)
    implicit none
    double precision, intent(in) :: T(np)
    double precision, intent(in), optional :: cpms(np,ns)
    integer :: i, m, n
    double precision :: tmp(np), crot(np) 
    double precision :: wru, dr, sqdr, dr32, aaaa1, dd, sqdd, dd32, bbbb
    double precision, parameter :: PI1=1.d0/3.1415926535D0, PI32O2=2.7842D+00, &
         P2O4P2=4.4674D+00, PI32=5.5683D+00

    !DEC$ SIMD
    do i=1,np
       dlt(i,1) = log(T(i))
       dlt(i,2) = dlt(i,1) * dlt(i,1)
       dlt(i,3) = dlt(i,2) * dlt(i,1)
       dlt(i,4) = dlt(i,3) * dlt(i,1)
       dlt(i,5) = dlt(i,4) * dlt(i,1)
       dlt(i,6) = dlt(i,5) * dlt(i,1)
    end do

    do n=1,ns
       !DEC$ SIMD
       do i=1,np
          etalg(i,n) = cfe(1,n) + cfe(2,n)*dlt(i,1) + cfe(3,n)*dlt(i,2) + cfe(4,n)*dlt(i,3)
          eta(i,n) = exp(etalg(i,n))
       end do
    end do

    if (iflag .le. 1) return

    do n=1,ns
       do m=1,n-1
          !DEC$ SIMD
          do i=1,np
             tmp(i) = -(cfd(1,m,n)+cfd(2,m,n)*dlt(i,1)+cfd(3,m,n)*dlt(i,2) &
                  + cfd(4,m,n)*dlt(i,3))
          end do
          !DEC$ SIMD
          do i=1,np
             bin(i,m,n) = exp(tmp(i))
             bin(i,n,m) = bin(i,m,n)
          end do
       end do
       do i=1,np
          bin(i,n,n) = 0.d0
       end do
    end do

    if (iflag .le. 2) return

    if (iflag.eq.3 .or. iflag.eq.5) then
       do n=1,ns
          do m=1,n-1
             do i=1,np
                A(i,m,n) = fita(1,m,n) + fita(2,m,n)*dlt(i,1) + fita(3,m,n)*dlt(i,2) &
                     + fita(4,m,n)*dlt(i,3) + fita(5,m,n)*dlt(i,4) &
                     + fita(6,m,n)*dlt(i,5) + fita(7,m,n)*dlt(i,6)
                A(i,n,m) = A(i,m,n)
             end do
          end do
          do i=1,np
             A(i,n,n) = fita(1,n,n) + fita(2,n,n)*dlt(i,1) + fita(3,n,n)*dlt(i,2) &
                  + fita(4,n,n)*dlt(i,3) + fita(5,n,n)*dlt(i,4) &
                  + fita(6,n,n)*dlt(i,5) + fita(7,n,n)*dlt(i,6)
          end do
       end do
    end if

    if (iflag .eq. 3) return

!-----------------------------------------------------------------------
!         COMPUTE PARKER CORRECTION FOR ZROT
!         AND ALSO THE ROTATIONAL AND INTERNAL PARTS OF SPECIFIC HEAT
!-----------------------------------------------------------------------
    do n=1,ns
       select case(nlin(n))
       case (0)
          do i=1,np
             crot(i) = 0.d0
             cint(i,n) = 0.d0
          end do
       case (1)
          wru = wt(n) / Ru
          do i=1,np
             crot(i) = 1.d0
             cint(i,n) = cpms(i,n) * wru - 2.50d0
          end do
       case (2)
          wru = wt(n) / Ru
          do i=1,np
             crot(i) = 1.5d0
             cint(i,n) = cpms(i,n) * wru - 2.50d0
          end do
       case default
          print *, "EFZ: wrong value in nlin"
          stop
       end select
       
       dr = eps(n) / 298.d0
       sqdr = sqrt(dr)
       dr32 = sqdr*dr
       aaaa1 = 1.d0/((1.0d0 + PI32O2*sqdr + P2O4P2*dr + PI32*dr32) * max(1.0d0, zrot(n)))

       do i=1,np
          dd = eps(n) / T(i)
          sqdd = sqrt(dd)
          dd32 = sqdd*dd
          bbbb = (1.0d0 + PI32O2*sqdd + P2O4P2*dd + PI32*dd32) 
          cxi(i,n) = crot(i) * PI1 * bbbb * aaaa1
       end do
    end do

    return
  end subroutine LZPAR


  ! shear viscosity
  subroutine EGZE1(alpha, X, mu)
    implicit none
    double precision, intent(in) :: alpha
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: mu(np)

    integer :: i, n
    double precision :: alpha1

    mu = 0.d0
    if (alpha .eq. 0.d0) then
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)*etalg(i,n)
          end do
       end do
       do i=1,np
          mu(i) = exp(mu(i))
       end do
    else if (alpha .eq. 1.d0) then
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)*eta(i,n)
          end do
       end do
    else if (alpha .eq. -1.d0) then
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)/eta(i,n)
          end do
       end do
       do i=1,np
          mu(i) = 1.d0/mu(i)
       end do
    else
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)*exp(alpha*etalg(i,n))
          end do
       end do
       alpha1 = 1.d0/alpha
       do i=1,np
          mu(i) = mu(i)**alpha1
       end do       
    end if

  end subroutine EGZE1


  ! shear viscosity
  subroutine EGZE3(T, mu)
    implicit none
    double precision, intent(in) :: T (np)
    double precision, intent(out) :: mu(np)

    integer :: i, n

    call EGZEMH(T)
    
    rn = beta

    call EGZCG1(1)

    mu = 0.d0
    do n=1,ns
       do i=1,np
          mu(i) = mu(i) + an(i,n) * beta(i,n)
       end do
    end do
  end subroutine EGZE3

  subroutine EGZEMH(T)
    implicit none
    double precision, intent(in) :: T (np)
    integer :: i, m, n
    double precision :: FAC(np), CCC, wtfac, wtmn, wtnm, aaa

    do i = 1, np
       ! EVALUATE THE MATRIX H
       ! Note:    FAC * BIN = 2 W_{ij} / \eta_{ij} / A_{ij}
       FAC(i) = (6.0D0 * RU / ( 5.0D0 * PATMOS )) * T(i)
    end do

    CCC = 5.0D0 / 3.0D0

    do n=1,ns
       do i=1,np
          G(i,n,n) = xtr(i,n) * xtr(i,n) / eta(i,n)
          ! EVALUATE THE RHS BETA
          beta(i,n) = xtr(i,n)
       end do
    end do

    do n=1,ns
       do m=1,n-1
          wtfac = 1.d0/(wt(m) + wt(n))
          wtmn = wt(m)*iwt(n)
          wtnm = wt(n)*iwt(m)
          !DEC$ SIMD PRIVATE(aaa)
          do i=1,np
             aaa = bin(i,m,n) * xtr(i,n) * xtr(i,m) * FAC(i) * wtfac
             G(i,m,m) = G(i,m,m) + aaa*(A(i,m,n)*wtnm + CCC)
             G(i,m,n) = aaa * (A(i,m,n) - CCC)
             G(i,n,m) = G(i,m,n)
             G(i,n,n) = G(i,n,n) + aaa*(A(i,m,n)*wtmn + CCC)
          end do
       end do
    end do

  end subroutine EGZEMH


  subroutine EGZCG1(itmax)
    implicit none
    integer, intent(in) :: itmax

    integer :: niter, i, n
    double precision :: betan(np), aaa(np), bbb(np), ccc(np), temp(np,ns)

    do i=1,np
       aaa(i) = 0.d0
       betan(i) = 0.d0
    end do

    do n = 1, ns
       do i=1,np
          an(i,n) = 0.0d0
          zn(i,n) = 0.0d0
          dmi(i,n) = 1.0D0 / G(i,n,n)
          aaa(i) = aaa(i) + dmi(i,n) * rn(i,n)*rn(i,n)
       end do
    enddo

    do niter=1, itmax
       do n=1, ns
          do i = 1, np
             zn(i,n) = dmi(i,n)*rn(i,n) + betan(i)*zn(i,n)
          enddo
       end do

       CALL EGZAXS(G, zn, temp)

       bbb = 0.d0
       do n=1,ns
          do i=1,np
             bbb(i) = bbb(i) + zn(i,n) * temp(i,n)
          end do
       end do

       do n=1,ns
          do i=1,np
             an(i,n) = an(i,n) + aaa(i)/bbb(i)*zn(i,n)
             rn(i,n) = rn(i,n) - aaa(i)/bbb(i)*temp(i,n)
          end do
       end do

       if (niter .eq. itmax) exit

       ccc = 0.d0
       do n=1,ns
          do i=1,np
             ccc(i) = ccc(i) + dmi(i,n) * rn(i,n)*rn(i,n)
          end do
       end do

       do i=1, np
          betan(i) = ccc(i) / aaa(i)
          aaa(i) = ccc(i)
       end do

    end do

  end subroutine EGZCG1

  subroutine EGZAXS(AA, X, B) ! B = AA.X, AA is symmetric.
    implicit none
    double precision, intent(in) :: AA(np,ns,ns)
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: B(np,ns)
    integer :: i, m, n
    B = 0.d0
    do n=1,ns
       B(:,n) = 0.d0
       do m=1,ns
          do i=1,np
             B(i,n) = B(i,n) + AA(i,m,n) * X(i,m)
          end do
       end do
    end do
  end subroutine EGZAXS


  ! volume viscosity
  subroutine EGZK1(alpha, X, VV)
    implicit none
    double precision, intent(in) :: alpha
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: VV(np)

    integer :: i, n
    double precision :: sxp(np), ccc, vvv, alpha1

    sxp = 0.d0
    do n=1,ns
       do i=1,np
          if (cxi(i,n) .ne. 0.d0) then
             sxp(i) = sxp(i) + X(i,n)
          end if
       end do
    end do
    do i=1,np
       sxp(i) = 1.d0/sxp(i)
       VV(i) = 0.d0
    end do

    if (alpha .eq. 0.d0) then
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = etalg(i,n) + log ( 0.25d0*ccc*ccc/cxi(i,n) )
                VV(i) = VV(i) + sxp(i) * X(i,n) * vvv
             end if
          end do
       end do
       do i=1,np
          VV(i) = exp(VV(i))
       end do
    else if (alpha .eq. 1.d0) then
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = eta(i,n) * 0.25d0*ccc*ccc/cxi(i,n) 
                VV(i) = VV(i) + sxp(i) * X(i,n) * vvv
             end if
          end do
       end do
    else if (alpha .eq. -1.d0) then
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = cxi(i,n) / (eta(i,n) * 0.25d0*ccc*ccc)
                VV(i) = VV(i) + sxp(i) * X(i,n) * vvv
             end if
          end do
       end do
       do i=1,np
          VV(i) = 1.d0/VV(i)
       end do       
    else
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = etalg(i,n) + log ( 0.25d0*ccc*ccc/cxi(i,n) )
                VV(i) = VV(i) + sxp(i) * X(i,n) * exp(alpha*vvv)
             end if
          end do
       end do
       alpha1 = 1.d0/alpha
       do i=1,np
          VV(i) = VV(i)**alpha1
       end do       
    end if

  end subroutine EGZK1

  ! volume viscosity
  subroutine EGZK3(T, VV)
    implicit none
    double precision, intent(in) :: T(np)
    double precision, intent(out) :: VV(np)
    
    integer :: i, m, n
    double precision :: ccc(np), wtfac, bb
    double precision, parameter :: denfac = 2.4d0 * Ru / Patmos

    if (.not. use_bulk_visc) then
       VV = 0.d0
       return
    end if

    ccc = 0.0d0
    do n=1,ns
       do i=1,np
          ccc(i) = ccc(i) + xtr(i,n)*cint(i,n)
       end do
    end do

    do n=1,ns
       do i=1,np
          G(i,n,n) = 4.d0*cxi(i,n)/eta(i,n)*xtr(i,n)*xtr(i,n)
          beta(i,n) = -xtr(i,n) * cint(i,n) / (ccc(i)+1.5d0)          
       end do
    end do

    do n=1,ns
       do m=1,n-1
          wtfac = iwt(m) + iwt(n)
          !DEC$ SIMD PRIVATE(bb)
          do i=1,np
             bb = xtr(i,m)*xtr(i,n)*bin(i,m,n)*denfac*T(i)*A(i,m,n)*wtfac
             G(i,m,m) = G(i,m,m) + bb*cxi(i,m)
!             G(i,n,m) = 0.d0
!             G(i,m,n) = 0.d0
             G(i,n,n) = G(i,n,n) + bb*cxi(i,n)
          end do
       end do
    end do

    VV = 0.d0
    do n=1,ns
       if (cxi(1,n) .eq. 0.d0) then
          do i=1,np
             VV(i) = VV(i) + beta(i,n) * beta(i,n) 
          end do
       else
          do i=1,np
             VV(i) = VV(i) + beta(i,n) * beta(i,n) / G(i,n,n)
          end do
       end if
    end do
  end subroutine EGZK3


  ! therml conductivity
  subroutine EGZL1(alpha, X, con)
    implicit none
    double precision, intent(in) :: alpha
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: con(np)

    integer :: i, n
    double precision :: asum(np), alpha1

    asum = 0.d0
    if (alpha .eq. 0.d0) then
       do n=1,ns
          do i=1,np
             asum(i) = asum(i) + X(i,n)*(cfl(1,n) + cfl(2,n)*dlt(i,1) &
                  + cfl(3,n)*dlt(i,2) + cfl(4,n)*dlt(i,3))
          end do
       end do
       do i=1,np
          con(i) = exp(asum(i)) 
       end do
    else
       alpha1 = 1.d0 / alpha
       do n=1,ns
          !DEC$ SIMD
          do i=1,np
             asum(i) = asum(i) + X(i,n)*exp(alpha*(cfl(1,n) + cfl(2,n)*dlt(i,1) &
                  + cfl(3,n)*dlt(i,2) + cfl(4,n)*dlt(i,3)))
          end do
       end do
       if (alpha .eq. 1.d0) then
          do i=1,np
             con(i) = asum(i)
          end do
       else if (alpha .eq. -1.d0) then
          do i=1,np
             con(i) = 1.d0/asum(i)
          end do
       else
          !DEC$ SIMD
          do i=1,np
             con(i) = asum(i)**alpha1
          end do
       end if
    end if
    
  end subroutine EGZL1


  ! rho * flux diffusion coefficients
  subroutine EGZVR1(T, D)
    implicit none
    double precision, intent(in) :: T(np)
    double precision, intent(out) :: D(np,ns)

    integer :: i, m, n
    double precision :: fac(np)

    D = 0.d0
    do n=1,ns
       do m=1,ns
          !DEC$ SIMD
          do i=1,np
             D(i,m) = D(i,m) + xtr(i,n)*bin(i,m,n)
          end do
       end do
    end do

    do i=1,np
       fac(i) = (Patmos/Ru) / T(i)
    end do

    do n=1,ns
       do i=1,np
          D(i,n) = wt(n) * fac(i) * aux(i,n) / D(i,n)
       end do
    end do

  end subroutine EGZVR1

end module egz_module

