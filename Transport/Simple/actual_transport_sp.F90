#include "mechanism.h"

module actual_transport_module

  use amrex_fort_module, only : amrex_real
  use eos_type_module
  use transport_type_module
  use chemistry_module, only : Ru
  use network, only : nspec

  implicit none

  character (len=64) :: transport_name = "smp"
  logical, save, private :: smp_initialized = .false.
  integer, save, private :: npts_smp = 0
!  integer, parameter :: NUM_SPECIES = 9
!  integer, parameter :: nfit = 4
!  integer, parameter ::  norder = 3
  integer, parameter ::  iter = 1

  real(amrex_real), parameter ::  trace = 1.d-15
  real(amrex_real), parameter ::  Patm = 1.01325d6
  real(amrex_real)                    :: Tloc, rholoc, wbar, pscale
  real(amrex_real)                    :: Yloc(NUM_SPECIES), Xloc(NUM_SPECIES)
  real(amrex_real)                    :: muloc(NUM_SPECIES), lamloc(NUM_SPECIES), dbinloc(NUM_SPECIES,NUM_SPECIES)
  real(amrex_real)                    :: xiloc(NUM_SPECIES)
  real(amrex_real)                    :: logT(NUM_FIT-1)

  real(amrex_real), save :: wt(NUM_SPECIES), iwt(NUM_SPECIES), eps(NUM_SPECIES), sig(NUM_SPECIES)
  real(amrex_real), save :: dip(NUM_SPECIES), pol(NUM_SPECIES), zrot(NUM_SPECIES)
  integer, save          :: nlin(NUM_SPECIES)

  real(amrex_real), save :: fitmu(NUM_FIT,NUM_SPECIES),fitlam(NUM_FIT,NUM_SPECIES),fitdbin(NUM_FIT,NUM_SPECIES,NUM_SPECIES)

  logical, parameter :: use_bulk_viscosity = .true.

contains

  ! This subroutine should be called outside OMP PARALLEL
  subroutine actual_transport_init

    implicit none
 
    call egtransetKK(nspec)
!   call egtransetNO(nfit)

!   allocate(wt(nspec))
!   allocate(iwt(nspec))
!   allocate(eps(nspec))
!   allocate(sig(nspec))
!   allocate(dip(nspec))
!   allocate(pol(nspec))
!   allocate(zrot(nspec))
!   allocate(nlin(nspec))

!   allocate(fitmu(nfit,nspec))
!   allocate(fitlam(nfit,nspec))
!   allocate(fitdbin(nfit,nspec,nspec))

    call egtransetWT(wt)
    call egtransetEPS(eps)
    call egtransetSIG(sig)
    call egtransetDIP(dip)
    call egtransetPOL(pol)
    call egtransetZROT(zrot)
    call egtransetNLIN(nlin)
    call egtransetCOFETA(fitmu)
    call egtransetCOFLAM(fitlam)
    call egtransetCOFD(fitdbin)

    iwt = 1.d0/wt

    smp_initialized = .true.

  end subroutine actual_transport_init


  subroutine actual_transport_close

    implicit none

!   if( allocated(wt)) deallocate(wt)
!   if( allocated(iwt)) deallocate(iwt)
!   if( allocated(eps)) deallocate(eps)
!   if( allocated(sig)) deallocate(sig)
!   if( allocated(dip)) deallocate(dip)
!   if( allocated(pol)) deallocate(pol)
!   if( allocated(zrot)) deallocate(zrot)
!   if( allocated(nlin)) deallocate(nlin)

!   deallocate(fitmu)
!   deallocate(fitlam)
!   deallocate(fitdbin)

    smp_initialized = .false.

  end subroutine actual_transport_close

  subroutine actual_transport(which, coeff)

    use amrex_error_module
    
    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff
    integer :: i, n, nn

    if (.not. smp_initialized) then
       call amrex_error('simple ::actual_transport called before initialized')
    endif

       Yloc(:) = coeff % eos_state(1) % massfrac(:)    
       Tloc    = coeff % eos_state(1) % T
       rholoc  = coeff % eos_state(1) % rho
       logT(1) = log(Tloc)
       logT(2) = logT(1)**2
       logT(3) = logT(1)*logT(2)


! define trace mole and mass fraction to get correct limiting
! behavior of species diffusion as species disappear


!  just to not waste space
     Xloc(1) = sum(Yloc(:))
     wbar = 0.d0


   do n = 1, nspec
       
      Yloc(n) = Yloc(n) + trace*(Xloc(1)/dble(nspec)-Yloc(n)) 

   enddo

   do n = 1, nspec
       
      wbar = wbar + Yloc(n)*iwt(n)

   enddo

   wbar = 1.d0/ wbar
 

   do n = 1, nspec
       
     Xloc(n) = Yloc(n)*wbar*iwt(n)

   enddo

    if (which % wtr_get_mu) then
       do n=1,nspec

            muloc(n) = fitmu(1,n)+fitmu(2,n)*logT(1)+ fitmu(3,n)*logT(2)  &
                        + fitmu(4,n)*logT(3)
            muloc(n) = exp(muloc(n))
           
         enddo

       coeff % mu(1) = 0.d0

       do n=1,nspec

            coeff%mu(1) = coeff%mu(1)+ Xloc(n)*muloc(n)**6.d0

       enddo


       coeff%mu(1) = coeff%mu(1)**(1.d0/6.d0)


!  assumption that we only get bulk viscosity if we are already getting shear viscosity
    if (which % wtr_get_xi) then

       call comp_pure_bulk(coeff)

       coeff % xi(1) = 0.d0

       do n=1,nspec

            coeff%xi(1) = coeff%xi(1)+ Xloc(n)*xiloc(n)**0.75d0

       enddo

       coeff%xi(1) = coeff%xi(1)**(4.d0/3.d0)


    endif

    endif

    if (which % wtr_get_lam) then

       do n=1,nspec

            lamloc(n) = fitlam(1,n)+fitlam(2,n)*logT(1)+ fitlam(3,n)*logT(2)  &
                        + fitlam(4,n)*logT(3)
            lamloc(n) = exp(lamloc(n))
           
       enddo

       coeff % lam(1) = 0.d0

       do n=1,nspec

            coeff%lam(1) = coeff%lam(1)+ Xloc(n)*lamloc(n)**0.25d0

       enddo


       coeff%lam(1) = coeff%lam(1)**4


    endif


!   if (which % wtr_get_Ddiag .or. which % wtr_get_Dmat) then
    if (which % wtr_get_Ddiag ) then

       do n=1,nspec
          do nn=1,n-1

               dbinloc(n,nn) = fitdbin(1,n,nn)+fitdbin(2,n,nn)*logT(1)   &
                   + fitdbin(3,n,nn)*logT(2)+ fitdbin(4,n,nn)*logT(3)
               dbinloc(n,nn) = exp(dbinloc(n,nn))
   
               dbinloc(nn,n) = dbinloc(n,nn)
           
          enddo

          dbinloc(n,n) = 0.d0

       enddo


        pscale = Patm * wbar/ (Ru * coeff % eos_state(1) % T * coeff % eos_state(1) % rho)


       if (which % wtr_get_Ddiag ) then

           call mixture(coeff)

           do n=1,nspec

                 coeff%Ddiag(1,n) = rholoc*pscale*coeff%Ddiag(1,n)

           enddo

!      elseif (.false.) then

!          call matrix(coeff, coeff%npts)


!          do n=1,nspec
!          do nn=1,nspec
!          do i=1,coeff%npts

!                coeff%Dmat(i,nn,n) = rholoc(i)*pscale(i)*coeff%Dmat(i,nn,n)

!          enddo
!          enddo
!          enddo
         
       endif
    endif

  end subroutine actual_transport


  subroutine comp_pure_bulk(coeff)

  implicit none

  type (trv_t), intent(inout) :: coeff

  real(amrex_real) :: cvk(nspec), cvkint(nspec), cvkrot(nspec)
  real(amrex_real) :: FofT(nspec), Fnorm(nspec), epskoverT

  real(amrex_real), parameter :: pi = 3.141592653589793238d0
 

  integer n,i

  call ckcvms( coeff % eos_state(1) % T, coeff % eos_state(1)%cvi )

  do i=1,nspec
 
        if(nlin(i) .eq.0)then

           cvkint(i) = 0.d0
           cvkrot(i) = 0.d0

         elseif(nlin(i).eq.1)then

           cvkint(i) =( coeff% eos_state(1)% cvi(i) ) * wt(i)/Ru - 1.5d0
           cvkrot(i) = 1.d0

         else

           cvkint(i) =( coeff% eos_state(1)% cvi(i) ) * wt(i)/Ru - 1.5d0
           cvkrot(i) = 1.5d0

         endif

   enddo

   do i = 1,nspec

      epskoverT = eps(i)/298.d0

      Fnorm(i) = 1.d0 + 0.5d0*pi**1.5*sqrt(epskoverT) + (2.d0+.5d0*pi**2)*epskoverT &
                  +(pi*epskoverT)**1.5

   enddo

   do i = 1,nspec

         epskoverT = eps(i)/ Tloc

         FofT(i) = 1.d0 + 0.5d0*pi**1.5*sqrt(epskoverT) + (2.d0+.5d0*pi**2)*epskoverT &
                  +(pi*epskoverT)**1.5

   enddo
      
   do i=1,nspec

      if(nlin(i) .ne. 0)then


!   zrot/crot approximately zint / cint by assuming vibrational internal energy is small
!   cvkrot is scaled by wk / Ru = mk / kb relative to standard specific cv

          xiloc(i) = 0.25d0*pi*(cvkint(i)/(cvkint(i)+1.5d0))**2* zrot(i)/cvkrot(i)*  &
                          Fnorm(i)/FofT(i) * muloc(i)


        else

!     no bulk viscosity for monotmic species

          xiloc(i) = 0.d0

      endif

    enddo


  end subroutine comp_pure_bulk

  subroutine mixture(coeff)

  implicit none

  type (trv_t), intent(inout) :: coeff

  real(amrex_real) :: term1,term2
  integer j,k

  do j = 1, nspec
    term1 = 0.d0
    term2 = 0.d0
      do k = 1, nspec
       if(k.ne.j) then

            term1 = term1 + Yloc(k)
            term2 = term2 + Xloc(k)/dbinloc(k,j)

       endif
      enddo

    coeff%Ddiag(1,j) = wt(j)* term1/term2 / wbar

  enddo
   
  end subroutine mixture

! subroutine matrix(coeff,npts)

! implicit none

! type (trv_t), intent(inout) :: coeff

! integer :: npts

! integer ::  i, j, k, jj, n
! real(kind=8) ::  term1(npts), term2(npts)
! real(kind=8) :: Di(1:npts,1:nspec), Diff_ij(1:npts,1:nspec,1:nspec)
! real(kind=8) :: Deltamat(1:npts,1:nspec,1:nspec), Zmat(1:npts,1:nspec,1:nspec)
! real(kind=8), dimension(1:npts,1:nspec,1:nspec) :: Pmat, Jmat
! real(kind=8), dimension(1:npts,1:nspec) :: Minv, Mmat
! real(kind=8), dimension(1:npts,1:nspec,1:nspec) :: PJ, matrix1, matrix2
! real(kind=8) :: scr(npts)


!         ! Find Di matrix 
!         do i = 1, nspec
!          term1 = 0.0d0  
!          term2 = 0.0d0  
!          do j = 1, nspec
!           if(j.ne.i) then
!             do n=1,npts
!               term1(n) = term1(n) + Yloc(n,j)
!               term2(n) = term2(n) + Xloc(n,j)/dbinloc(n,i,j)
!             enddo
!           endif
!          enddo   
!            do n=1,npts
!               Di(n,i) = term1(n)/term2(n) 
!            enddo
!         enddo

!   
!         ! Compute Mmat and Minv
!         do i = 1, nspec
!           do n=1,npts
!          
!            Mmat(n,i) = Xloc(n,i)/Di(n,i)
!            Minv(n,i) = Di(n,i)/Xloc(n,i)

!         enddo
!         enddo
!   

!         ! Compute P matrix
!         Pmat = 0.0d0
!         do i = 1, nspec
!          do j = 1, nspec
!            do n=1,npts
!               Pmat(n,i,j) = - Yloc(n,j)
!            enddo
!               if(i.eq.j) then
!            do n=1,npts
!                Pmat(n,i,j) =  Pmat(n,i,j) + 1.0d0
!            enddo
!               endif
!          enddo
!         enddo


!         ! Compute Deltamat
!         Deltamat = 0.0d0
!         do i = 1, nspec
!          do j = 1, nspec
!            if(i.eq.j) then
!             term1 = 0.0d0
!             do k = 1, nspec
!               if(k.ne.i) then

!                  do n=1,npts

!                  term1(n) = term1(n) + Xloc(n,i)*Xloc(n,k)/dbinloc(n,i,k)

!                  enddo

!               endif
!             enddo

!                do n=1,npts
!                   Deltamat(n,i,i) = term1(n)
!                enddo

!               else
!                do n=1,npts
!                   Deltamat(n,i,j) = -Xloc(n,i)*Xloc(n,j)/dbinloc(n,i,j)
!                enddo
!               endif
!                do n=1,npts
!                   Zmat(n,i,j) = -Deltamat(n,i,j)
!                enddo
!          enddo
!         enddo


!         ! Compute Zmat
!         do i = 1, nspec
!           do n=1,npts
!              Zmat(n,i,i) = Zmat(n,i,i) + Mmat(n,i)
!           enddo
!         enddo

!         ! Compute Jmat
!         Jmat = 0.d0
!         do i = 1, nspec
!          do j = 1, nspec
!            do n=1,npts
!               Jmat(n,i,j) = Minv(n,i)*Zmat(n,i,j)
!            enddo
!           enddo
!          enddo

!         ! Compute PJ
!         PJ = 0.0d0
!         do i = 1, nspec
!          do j = 1, nspec
!           do k = 1, nspec
!             do n=1,npts
!               PJ(n,i,j) = PJ(n,i,j) + Pmat(n,i,k)*Jmat(n,k,j)
!             enddo
!           enddo
!          enddo
!         enddo



          ! Compute P M^-1 Pt; store it in matrix2
!         do i = 1, nspec
!          do j = 1, nspec
!           scr = 0.d0
!           do k = 1, nspec
!               do n=1,npts
!                  scr(n) = scr(n) + Pmat(n,i,k)*Minv(n,k)*Pmat(n,j,k)
!               enddo
!               ! notice the change in indices for Pmat to represent Pmat^t
!           enddo
!               do n=1,npts
!                 matrix2(n,i,j) = scr(n)
!                 Diff_ij(n,i,j) = scr(n)
!               enddo
!          enddo
!         enddo

!         if(iter.gt.0)then

!         do jj = 1,iter


!         matrix1=0
!         do i = 1, nspec
!          do j = 1, nspec
!           scr = 0.d0
!           do k = 1, nspec
!              do n=1,npts
!                 scr(n) = scr(n) + PJ(n,i,k)*Diff_ij(n,k,j)
!              enddo
!           enddo
!              do n=1,npts
!                matrix1(n,i,j) = scr(n)+matrix2(n,i,j)
!              enddo
!          enddo
!         enddo

!         Diff_ij=matrix1

!         enddo

!         endif



!         ! Compute D_tilde
!         do i = 1, nspec
!          do j = 1, nspec
!              do n=1,npts
!                 coeff%Dmat(n,i,j) = Diff_ij(n,i,j)*Yloc(n,i)
!              enddo
!          enddo
!         enddo





! end subroutine  matrix
end module actual_transport_module

