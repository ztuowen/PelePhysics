module fuego_module

  implicit none
  private
  public :: ckcpms
  public :: ckums
  public :: ckrhoy
  public :: ckcvms
  public :: ckxty
  public :: ckytcr
  public :: ckytx
  public :: ckytx_gpu
  public :: ckhms
  public :: ckhms_gpu
  public :: vckytx
  public :: vckhms
  public :: ckcvbs
  public :: ckubms
  public :: ckcpbs
  public :: ckpy
  public :: get_t_given_ey
  public :: cksyme
  public :: cksyms
  public :: ckwt
  public :: ckrp
  public :: ckwc
  public :: ckindx
  public :: ckinit
  public :: ckfinalize
  public :: egtransetCOFTD
  public :: egtransetKTDIF
  public :: egtransetCOFD
  public :: egtransetCOFLAM
  public :: egtransetCOFETA
  public :: egtransetNLIN
  public :: egtransetZROT
  public :: egtransetPOL
  public :: egtransetDIP
  public :: egtransetSIG
  public :: egtransetEPS
  public :: egtransetWT
  public :: egtransetNLITE
  public :: egtransetKK
  public :: egtransetNO
  public :: egtransetLENRMC
  public :: egtransetLENIMC

! Inverse molecular weights
!double precision, parameter :: imw(9) = (/ &
!    1.d0 / 2.015940d0,  & ! H2
!    1.d0 / 31.998800d0,  & ! O2
!    1.d0 / 18.015340d0,  & ! H2O
!    1.d0 / 1.007970d0,  & ! H
!    1.d0 / 15.999400d0,  & ! O
!    1.d0 / 17.007370d0,  & ! OH
!    1.d0 / 33.006770d0,  & ! HO2
!    1.d0 / 34.014740d0,  & ! H2O2
!    1.d0 / 28.013400d0 /) ! N2

type :: nonsquare_matrix_double
   double precision, allocatable :: vector(:)
end type nonsquare_matrix_double

type :: nonsquare_matrix_int
   integer, allocatable :: vector(:)
end type nonsquare_matrix_int

double precision, save :: fwd_A(21), fwd_beta(21), fwd_Ea(21)
double precision, save :: low_A(21), low_beta(21), low_Ea(21)
double precision, save :: rev_A(21), rev_beta(21), rev_Ea(21)
double precision, save :: troe_a(21),troe_Ts(21), troe_Tss(21), troe_Tsss(21)
double precision, save :: sri_a(21), sri_b(21), sri_c(21), sri_d(21), sri_e(21)
double precision, save :: activation_units(21), prefactor_units(21), phase_units(21)
integer, save :: is_PD(21), troe_len(21), sri_len(21), nTB(21)
type(nonsquare_matrix_double) :: TB(21)
type(nonsquare_matrix_int) :: TBid(21)

double precision, save :: fwd_A_DEF(21), fwd_beta_DEF(21), fwd_Ea_DEF(21)
double precision, save :: low_A_DEF(21), low_beta_DEF(21), low_Ea_DEF(21)
double precision, save :: rev_A_DEF(21), rev_beta_DEF(21), rev_Ea_DEF(21)
double precision, save :: troe_a_DEF(21),troe_Ts_DEF(21), troe_Tss_DEF(21), troe_Tsss_DEF(21)
double precision, save :: sri_a_DEF(21), sri_b_DEF(21), sri_c_DEF(21), sri_d_DEF(21), sri_e_DEF(21)
double precision, save :: activation_units_DEF(21), prefactor_units_DEF(21), phase_units_DEF(21)
integer, save :: is_PD_DEF(21), troe_len_DEF(21), sri_len_DEF(21), nTB_DEF(21)
type(nonsquare_matrix_double) :: TB_DEF(21)
type(nonsquare_matrix_int) :: TBid_DEF(21)

! productionRate() static variables
double precision, save :: T_save = -1
double precision, save :: k_f_save(21)
double precision, save :: Kc_save(21)

contains

subroutine SetAllDefaults()

    implicit none

    integer :: i, j

    do i=1, 21
        if (nTB_DEF(i) /= 0) then
            nTB_DEF(i) = 0
            if (allocated(TB_DEF(i) % vector)) deallocate(TB_DEF(i) % vector)
            if (allocated(TBid_DEF(i) % vector)) deallocate(TBid_DEF(i) % vector)
        end if

        fwd_A_DEF(i)    = fwd_A(i)
        fwd_beta_DEF(i) = fwd_beta(i)
        fwd_Ea_DEF(i)   = fwd_Ea(i)

        low_A_DEF(i)    = low_A(i)
        low_beta_DEF(i) = low_beta(i)
        low_Ea_DEF(i)   = low_Ea(i)

        rev_A_DEF(i)    = rev_A(i)
        rev_beta_DEF(i) = rev_beta(i)
        rev_Ea_DEF(i)   = rev_Ea(i)

        troe_a_DEF(i)    = troe_a(i)
        troe_Ts_DEF(i)   = troe_Ts(i)
        troe_Tss_DEF(i)  = troe_Tss(i)
        troe_Tsss_DEF(i) = troe_Tsss(i)

        sri_a_DEF(i) = sri_a(i)
        sri_b_DEF(i) = sri_b(i)
        sri_c_DEF(i) = sri_c(i)
        sri_d_DEF(i) = sri_d(i)
        sri_e_DEF(i) = sri_e(i)

        is_PD_DEF(i)    = is_PD(i)
        troe_len_DEF(i) = troe_len(i)
        sri_len_DEF(i)  = sri_len(i)

        activation_units_DEF(i) = activation_units(i)
        prefactor_units_DEF(i)  = prefactor_units(i)
        phase_units_DEF(i)      = phase_units(i)

        nTB_DEF(i)  = nTB(i)
        if (nTB_DEF(i) /= 0) then
           if (.not. allocated(TB_DEF(i) % vector)) allocate(TB_DEF(i) % vector(nTB_DEF(i)))
           if (.not. allocated(TBid_DEF(i) % vector)) allocate(TBid_DEF(i) % vector(nTB_DEF(i)))
           do j=1, nTB_DEF(i)
             TB_DEF(i) % vector(j) = TB(i) % vector(j)
             TBid_DEF(i) % vector(j) = TBid(i) % vector(j)
           end do
        end if
    end do

end subroutine


! Finalizes parameter database
subroutine ckfinalize()

  implicit none

  integer :: i

  do i=1, 21
    if (allocated(TB(i) % vector)) deallocate(TB(i) % vector)
    !TB(i) = 0
    if (allocated(TBid(i) % vector)) deallocate(TBid(i) % vector)
    !TBid(i) = 0
    nTB(i) = 0

    if (allocated(TB_DEF(i) % vector)) deallocate(TB_DEF(i) % vector)
    !TB_DEF(i) = 0
    if (allocated(TBid_DEF(i) % vector)) deallocate(TBid_DEF(i) % vector)
    !TBid_DEF(i) = 0
    nTB_DEF(i) = 0
  end do

end subroutine

! Initializes parameter database
subroutine ckinit()

    implicit none

    ! (0):  H + O2 <=> O + OH
    fwd_A(7)     = 3547000000000000d0
    fwd_beta(7)  = -0.40600000000000003d0
    fwd_Ea(7)    = 16599d0
    prefactor_units(7)  = 1.0000000000000002d-06
    activation_units(7) = 0.50321666580471969d0
    phase_units(7)      = 1d-12
    is_PD(7) = 0
    nTB(7) = 0

    ! (1):  O + H2 <=> H + OH
    fwd_A(8)     = 50800d0
    fwd_beta(8)  = 2.6699999999999999d0
    fwd_Ea(8)    = 6290d0
    prefactor_units(8)  = 1.0000000000000002d-06
    activation_units(8) = 0.50321666580471969d0
    phase_units(8)      = 1d-12
    is_PD(8) = 0
    nTB(8) = 0

    ! (2):  H2 + OH <=> H2O + H
    fwd_A(9)     = 216000000d0
    fwd_beta(9)  = 1.51d0
    fwd_Ea(9)    = 3430d0
    prefactor_units(9)  = 1.0000000000000002d-06
    activation_units(9) = 0.50321666580471969d0
    phase_units(9)      = 1d-12
    is_PD(9) = 0
    nTB(9) = 0

    ! (3):  O + H2O <=> OH + OH
    fwd_A(10)     = 2970000d0
    fwd_beta(10)  = 2.02d0
    fwd_Ea(10)    = 13400d0
    prefactor_units(10)  = 1.0000000000000002d-06
    activation_units(10) = 0.50321666580471969d0
    phase_units(10)      = 1d-12
    is_PD(10) = 0
    nTB(10) = 0

    ! (4):  H2 + M <=> H + H + M
    fwd_A(3)     = 4.577d+19
    fwd_beta(3)  = -1.3999999999999999d0
    fwd_Ea(3)    = 104380d0
    prefactor_units(3)  = 1.0000000000000002d-06
    activation_units(3) = 0.50321666580471969d0
    phase_units(3)      = 1d-6
    is_PD(3) = 0
    nTB(3) = 2
    if (.not. allocated(TB(3) % vector)) allocate(TB(3) % vector(2))
    if (.not. allocated(TBid(3) % vector)) allocate(TBid(3) % vector(2))
    TBid(3) % vector(1) = 0d0
    TB(3) % vector(1) = 2.5d0 ! H2
    TBid(3) % vector(2) = 2d0
    TB(3) % vector(2) = 12d0 ! H2O

    ! (5):  O + O + M <=> O2 + M
    fwd_A(4)     = 6165000000000000d0
    fwd_beta(4)  = -0.5d0
    fwd_Ea(4)    = 0d0
    prefactor_units(4)  = 1.0000000000000002d-12
    activation_units(4) = 0.50321666580471969d0
    phase_units(4)      = 1d-12
    is_PD(4) = 0
    nTB(4) = 2
    if (.not. allocated(TB(4) % vector)) allocate(TB(4) % vector(2))
    if (.not. allocated(TBid(4) % vector)) allocate(TBid(4) % vector(2))
    TBid(4) % vector(1) = 0d0
    TB(4) % vector(1) = 2.5d0 ! H2
    TBid(4) % vector(2) = 2d0
    TB(4) % vector(2) = 12d0 ! H2O

    ! (6):  O + H + M <=> OH + M
    fwd_A(5)     = 4.714d+18
    fwd_beta(5)  = -1d0
    fwd_Ea(5)    = 0d0
    prefactor_units(5)  = 1.0000000000000002d-12
    activation_units(5) = 0.50321666580471969d0
    phase_units(5)      = 1d-12
    is_PD(5) = 0
    nTB(5) = 2
    if (.not. allocated(TB(5) % vector)) allocate(TB(5) % vector(2))
    if (.not. allocated(TBid(5) % vector)) allocate(TBid(5) % vector(2))
    TBid(5) % vector(1) = 0d0
    TB(5) % vector(1) = 2.5d0 ! H2
    TBid(5) % vector(2) = 2d0
    TB(5) % vector(2) = 12d0 ! H2O

    ! (7):  H + OH + M <=> H2O + M
    fwd_A(6)     = 3.8000000000000004d+22
    fwd_beta(6)  = -2d0
    fwd_Ea(6)    = 0d0
    prefactor_units(6)  = 1.0000000000000002d-12
    activation_units(6) = 0.50321666580471969d0
    phase_units(6)      = 1d-12
    is_PD(6) = 0
    nTB(6) = 2
    if (.not. allocated(TB(6) % vector)) allocate(TB(6) % vector(2))
    if (.not. allocated(TBid(6) % vector)) allocate(TBid(6) % vector(2))
    TBid(6) % vector(1) = 0d0
    TB(6) % vector(1) = 2.5d0 ! H2
    TBid(6) % vector(2) = 2d0
    TB(6) % vector(2) = 12d0 ! H2O

    ! (8):  H + O2 (+M) <=> HO2 (+M)
    fwd_A(1)     = 1475000000000d0
    fwd_beta(1)  = 0.59999999999999998d0
    fwd_Ea(1)    = 0d0
    low_A(1)     = 6.366d+20
    low_beta(1)  = -1.72d0
    low_Ea(1)    = 524.79999999999995d0
    troe_a(1)    = 0.80000000000000004d0
    troe_Tsss(1) = 1.0000000000000001d-30
    troe_Ts(1)   = 1d+30
    troe_len(1)  = 3
    prefactor_units(1)  = 1.0000000000000002d-06
    activation_units(1) = 0.50321666580471969d0
    phase_units(1)      = 1d-12
    is_PD(1) = 1
    nTB(1) = 3
    if (.not. allocated(TB(1) % vector)) allocate(TB(1) % vector(3))
    if (.not. allocated(TBid(1) % vector)) allocate(TBid(1) % vector(3))
    TBid(1) % vector(1) = 0d0
    TB(1) % vector(1) = 2d0 ! H2
    TBid(1) % vector(2) = 2d0
    TB(1) % vector(2) = 11d0 ! H2O
    TBid(1) % vector(3) = 1d0
    TB(1) % vector(3) = 0.78000000000000003d0 ! O2

    ! (9):  HO2 + H <=> H2 + O2
    fwd_A(11)     = 16600000000000d0
    fwd_beta(11)  = 0d0
    fwd_Ea(11)    = 823d0
    prefactor_units(11)  = 1.0000000000000002d-06
    activation_units(11) = 0.50321666580471969d0
    phase_units(11)      = 1d-12
    is_PD(11) = 0
    nTB(11) = 0

    ! (10):  HO2 + H <=> OH + OH
    fwd_A(12)     = 70790000000000d0
    fwd_beta(12)  = 0d0
    fwd_Ea(12)    = 295d0
    prefactor_units(12)  = 1.0000000000000002d-06
    activation_units(12) = 0.50321666580471969d0
    phase_units(12)      = 1d-12
    is_PD(12) = 0
    nTB(12) = 0

    ! (11):  HO2 + O <=> O2 + OH
    fwd_A(13)     = 32500000000000d0
    fwd_beta(13)  = 0d0
    fwd_Ea(13)    = 0d0
    prefactor_units(13)  = 1.0000000000000002d-06
    activation_units(13) = 0.50321666580471969d0
    phase_units(13)      = 1d-12
    is_PD(13) = 0
    nTB(13) = 0

    ! (12):  HO2 + OH <=> H2O + O2
    fwd_A(14)     = 28900000000000d0
    fwd_beta(14)  = 0d0
    fwd_Ea(14)    = -497d0
    prefactor_units(14)  = 1.0000000000000002d-06
    activation_units(14) = 0.50321666580471969d0
    phase_units(14)      = 1d-12
    is_PD(14) = 0
    nTB(14) = 0

    ! (13):  HO2 + HO2 <=> H2O2 + O2
    fwd_A(15)     = 420000000000000d0
    fwd_beta(15)  = 0d0
    fwd_Ea(15)    = 11982d0
    prefactor_units(15)  = 1.0000000000000002d-06
    activation_units(15) = 0.50321666580471969d0
    phase_units(15)      = 1d-12
    is_PD(15) = 0
    nTB(15) = 0

    ! (14):  HO2 + HO2 <=> H2O2 + O2
    fwd_A(16)     = 130000000000d0
    fwd_beta(16)  = 0d0
    fwd_Ea(16)    = -1629.3d0
    prefactor_units(16)  = 1.0000000000000002d-06
    activation_units(16) = 0.50321666580471969d0
    phase_units(16)      = 1d-12
    is_PD(16) = 0
    nTB(16) = 0

    ! (15):  H2O2 (+M) <=> OH + OH (+M)
    fwd_A(2)     = 295100000000000d0
    fwd_beta(2)  = 0d0
    fwd_Ea(2)    = 48430d0
    low_A(2)     = 1.202d+17
    low_beta(2)  = 0d0
    low_Ea(2)    = 45500d0
    troe_a(2)    = 0.5d0
    troe_Tsss(2) = 1.0000000000000001d-30
    troe_Ts(2)   = 1d+30
    troe_len(2)  = 3
    prefactor_units(2)  = 1d0
    activation_units(2) = 0.50321666580471969d0
    phase_units(2)      = 1d-6
    is_PD(2) = 1
    nTB(2) = 2
    if (.not. allocated(TB(2) % vector)) allocate(TB(2) % vector(2))
    if (.not. allocated(TBid(2) % vector)) allocate(TBid(2) % vector(2))
    TBid(2) % vector(1) = 0d0
    TB(2) % vector(1) = 2.5d0 ! H2
    TBid(2) % vector(2) = 2d0
    TB(2) % vector(2) = 12d0 ! H2O

    ! (16):  H2O2 + H <=> H2O + OH
    fwd_A(17)     = 24100000000000d0
    fwd_beta(17)  = 0d0
    fwd_Ea(17)    = 3970d0
    prefactor_units(17)  = 1.0000000000000002d-06
    activation_units(17) = 0.50321666580471969d0
    phase_units(17)      = 1d-12
    is_PD(17) = 0
    nTB(17) = 0

    ! (17):  H2O2 + H <=> HO2 + H2
    fwd_A(18)     = 48200000000000d0
    fwd_beta(18)  = 0d0
    fwd_Ea(18)    = 7950d0
    prefactor_units(18)  = 1.0000000000000002d-06
    activation_units(18) = 0.50321666580471969d0
    phase_units(18)      = 1d-12
    is_PD(18) = 0
    nTB(18) = 0

    ! (18):  H2O2 + O <=> OH + HO2
    fwd_A(19)     = 9550000d0
    fwd_beta(19)  = 2d0
    fwd_Ea(19)    = 3970d0
    prefactor_units(19)  = 1.0000000000000002d-06
    activation_units(19) = 0.50321666580471969d0
    phase_units(19)      = 1d-12
    is_PD(19) = 0
    nTB(19) = 0

    ! (19):  H2O2 + OH <=> HO2 + H2O
    fwd_A(20)     = 1000000000000d0
    fwd_beta(20)  = 0d0
    fwd_Ea(20)    = 0d0
    prefactor_units(20)  = 1.0000000000000002d-06
    activation_units(20) = 0.50321666580471969d0
    phase_units(20)      = 1d-12
    is_PD(20) = 0
    nTB(20) = 0

    ! (20):  H2O2 + OH <=> HO2 + H2O
    fwd_A(21)     = 580000000000000d0
    fwd_beta(21)  = 0d0
    fwd_Ea(21)    = 9557d0
    prefactor_units(21)  = 1.0000000000000002d-06
    activation_units(21) = 0.50321666580471969d0
    phase_units(21)      = 1d-12
    is_PD(21) = 0
    nTB(21) = 0

    call SetAllDefaults()

end subroutine


! A few mechanism parameters
subroutine ckindx(mm, kk, ii, nfit)

    implicit none

    integer, intent(out) :: mm
    integer, intent(out) :: kk
    integer, intent(out) :: ii
    integer, intent(out) :: nfit

    mm = 3
    kk = 9
    ii = 21
    nfit = -1 ! Why do you need this anyway?

end subroutine

! Returns the char strings of element names
subroutine cksyme(kname, plenkname)

    implicit none

    integer, intent(in) :: plenkname
    integer, intent(out) :: kname(plenkname*3)

    integer :: i
    integer :: lenkname

    lenkname = plenkname

    !clear kname
    do i=1, lenkname*3
        kname(i) = ichar(' ')
    end do

    ! H 
    kname(0*lenkname+1) = ichar('H')
    kname(0*lenkname+2) = ichar(' ')
    ! O 
    kname(1*lenkname+1) = ichar('O')
    kname(1*lenkname+2) = ichar(' ')
    ! N 
    kname(2*lenkname+1) = ichar('N')
    kname(2*lenkname+2) = ichar(' ')

end subroutine

! Returns the char strings of species names
subroutine cksyms(kname, plenkname)

    implicit none

    integer, intent(in) :: plenkname
    integer, intent(out) :: kname(plenkname*9)

    integer :: i
    integer :: lenkname

    lenkname = plenkname

    !clear kname
    do i=1, lenkname*9
        kname(i) = ichar(' ')
    end do

    ! H2 
    kname(0*lenkname+1) = ichar('H')
    kname(0*lenkname+2) = ichar('2')
    kname(0*lenkname+3) = ichar(' ')
    ! O2 
    kname(1*lenkname+1) = ichar('O')
    kname(1*lenkname+2) = ichar('2')
    kname(1*lenkname+3) = ichar(' ')
    ! H2O 
    kname(2*lenkname+1) = ichar('H')
    kname(2*lenkname+2) = ichar('2')
    kname(2*lenkname+3) = ichar('O')
    kname(2*lenkname+4) = ichar(' ')
    ! H 
    kname(3*lenkname+1) = ichar('H')
    kname(3*lenkname+2) = ichar(' ')
    ! O 
    kname(4*lenkname+1) = ichar('O')
    kname(4*lenkname+2) = ichar(' ')
    ! OH 
    kname(5*lenkname+1) = ichar('O')
    kname(5*lenkname+2) = ichar('H')
    kname(5*lenkname+3) = ichar(' ')
    ! HO2 
    kname(6*lenkname+1) = ichar('H')
    kname(6*lenkname+2) = ichar('O')
    kname(6*lenkname+3) = ichar('2')
    kname(6*lenkname+4) = ichar(' ')
    ! H2O2 
    kname(7*lenkname+1) = ichar('H')
    kname(7*lenkname+2) = ichar('2')
    kname(7*lenkname+3) = ichar('O')
    kname(7*lenkname+4) = ichar('2')
    kname(7*lenkname+5) = ichar(' ')
    ! N2 
    kname(8*lenkname+1) = ichar('N')
    kname(8*lenkname+2) = ichar('2')
    kname(8*lenkname+3) = ichar(' ')

end subroutine

! Returns R, Rc, Patm
subroutine ckrp(ru, ruc, pa)

    implicit none

    double precision, intent(out) :: ru
    double precision, intent(out) :: ruc
    double precision, intent(out) :: pa

    ru  = 83145100.000000d0 
    ruc = 1.98721558317399615845d0 
    pa  = 1013250.000000d0 

end subroutine

! Compute P = rhoRT/W(y)
subroutine ckpy(rho, T, y, P)

    !$acc routine(ckpy) seq

    implicit none

    double precision, intent(in) :: rho
    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    double precision, intent(out) :: P

    double precision :: YOW ! for computing mean MW
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    YOW = 0.d0

    YOW = YOW + (y(1) * imw(1)) ! H2
    YOW = YOW + (y(2) * imw(2)) ! O2
    YOW = YOW + (y(3) * imw(3)) ! H2O
    YOW = YOW + (y(4) * imw(4)) ! H
    YOW = YOW + (y(5) * imw(5)) ! O
    YOW = YOW + (y(6) * imw(6)) ! OH
    YOW = YOW + (y(7) * imw(7)) ! HO2
    YOW = YOW + (y(8) * imw(8)) ! H2O2
    YOW = YOW + (y(9) * imw(9)) ! N2

    ! YOW holds the reciprocal of the mean molecular wt
    P = rho * 8.31451000d+07 * T * YOW ! P = rho*R*T/W

end subroutine

! Compute rho = P*W(y)/RT
subroutine ckrhoy(P, T, y, rho)

    implicit none

    double precision, intent(in) :: P
    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    double precision, intent(out) :: rho

    double precision :: YOW, tmp(9)
    integer :: i
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    YOW = 0.d0

    do i=1, 9
        tmp(i) = y(i) * imw(i)
    end do
    do i=1, 9
        YOW = YOW + tmp(i)
    end do

    rho = P / ( 8.31451000d+07 * T * YOW) ! rho = P*W/(R*T)

end subroutine

! get molecular weight for all species
subroutine ckwt(wt)

    implicit none

    double precision, intent(out) :: wt(9)

    call molecularWeight(wt)

end subroutine

! convert y[species] (mass fracs) to x[species] (mole fracs)
subroutine ckytx(y, x)

    !$acc routine(ckytx) seq

    implicit none

    double precision, intent(in) :: y(9)
    double precision, intent(out) :: x(9)

    double precision :: YOW, YOWINV
    double precision :: tmp(9)
    integer :: i
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    do i=1, 9
        tmp(i) = y(i) * imw(i)
    end do
    do i=1, 9
        YOW = YOW + tmp(i)
    end do

    YOWINV = 1.d0 / YOW

    do i=1, 9
        x(i) = y(i) * imw(i) * YOWINV
    end do

end subroutine

! convert y[species] (mass fracs) to x[species] (mole fracs)
subroutine ckytx_gpu(q, x, lo, hi, i, j, k, qfs, qvar)

    !$acc routine(ckytx_gpu) seq

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: i, j, k, qfs, qvar
    double precision, intent(in), dimension(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:qvar) :: q
    double precision, intent(out), dimension(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:9) :: x

    double precision :: YOW
    integer :: n
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    YOW = 0.d0

    do n=1, 9
        YOW = YOW + q(i,j,k,qfs+n-1) * imw(n)
    end do

    do n=1, 9
        x(i,j,k,n) = q(i,j,k,qfs+n-1) * imw(n) * 1.d0 / YOW
    end do

end subroutine

! convert y(npoints,species) (mass fracs) to x(npoints,species) (mole fracs)
subroutine vckytx(np, y, x)

    implicit none

    integer, intent(in) :: np
    double precision, intent(in) :: y(np,9)
    double precision, intent(inout) :: x(np,9)

    double precision :: YOW(np)
    integer :: i, n
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    do i=1, np
        YOW(i) = 0.d0
    end do

    do n=1, 9
        do i=1, np
            x(i,n) = y(i,n) * imw(n)
            YOW(i) = YOW(i) + x(i,n)
        end do
    end do

    do i=1, np
        YOW(i) = 1.d0 / YOW(i)
    end do

    do n=1, 9
        do i=1, np
            x(i,n) = x(i,n) * YOW(i)
        end do
    end do

end subroutine

! convert y[species] (mass fracs) to c[species] (molar conc)
subroutine ckytcr(rho, T, y, c)

    implicit none

    double precision, intent(in) :: rho
    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    double precision, intent(out) :: c(9)

    integer :: i
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    do i=1, 9
        c(i) = rho * y(i) * imw(i)
    end do

end subroutine

! convert x[species] (mole fracs) to y[species] (mass fracs)
subroutine ckxty(x, y)

    implicit none

    double precision, intent(in) :: x(9)
    double precision, intent(out) :: y(9)

    double precision :: XW, XWinv

    XW = 0.d0 ! See Eq 4, 9 in CK Manual

    ! Compute mean molecular wt first
    XW = XW + (x(1) * 2.01594000d+00) ! H2
    XW = XW + (x(2) * 3.19988000d+01) ! O2
    XW = XW + (x(3) * 1.80153400d+01) ! H2O
    XW = XW + (x(4) * 1.00797000d+00) ! H
    XW = XW + (x(5) * 1.59994000d+01) ! O
    XW = XW + (x(6) * 1.70073700d+01) ! OH
    XW = XW + (x(7) * 3.30067700d+01) ! HO2
    XW = XW + (x(8) * 3.40147400d+01) ! H2O2
    XW = XW + (x(9) * 2.80134000d+01) ! N2

    ! Now compute conversion
    XWinv = 1.d0 / XW
    y(1) = x(1) * 2.01594000d+00 * XWinv 
    y(2) = x(2) * 3.19988000d+01 * XWinv 
    y(3) = x(3) * 1.80153400d+01 * XWinv 
    y(4) = x(4) * 1.00797000d+00 * XWinv 
    y(5) = x(5) * 1.59994000d+01 * XWinv 
    y(6) = x(6) * 1.70073700d+01 * XWinv 
    y(7) = x(7) * 3.30067700d+01 * XWinv 
    y(8) = x(8) * 3.40147400d+01 * XWinv 
    y(9) = x(9) * 2.80134000d+01 * XWinv 

end subroutine

! Returns the specific heats at constant volume
! in mass units (Eq. 29)
subroutine ckcvms(T, cvms)

    !$acc routine(ckcvms) seq
    !$acc routine(cv_r) seq

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout) :: cvms(9)

    double precision :: tT, tc(5)

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cv_R(cvms, tc)

    ! multiply by R/molecularweight
    cvms(1) = cvms(1) * 4.124383662212169d+07 !H2
    cvms(2) = cvms(2) * 2.598381814318037d+06 !O2
    cvms(3) = cvms(3) * 4.615239012974499d+06 !H2O
    cvms(4) = cvms(4) * 8.248767324424338d+07 !H
    cvms(5) = cvms(5) * 5.196763628636074d+06 !O
    cvms(6) = cvms(6) * 4.888768810227566d+06 !OH
    cvms(7) = cvms(7) * 2.519031701678171d+06 !HO2
    cvms(8) = cvms(8) * 2.444384405113783d+06 !H2O2
    cvms(9) = cvms(9) * 2.968047434442088d+06 !N2

end subroutine

! Returns the specific heats at constant pressure
! in mass units (Eq. 26)
subroutine ckcpms(T, cpms)

    !$acc routine(ckcpms) seq

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout) :: cpms(9)

    double precision :: tT, tc(5)

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cp_R(cpms, tc)

    ! multiply by R/molecularweight
    cpms(1) = cpms(1) * 4.124383662212169d+07 ! H2
    cpms(2) = cpms(2) * 2.598381814318037d+06 ! O2
    cpms(3) = cpms(3) * 4.615239012974499d+06 ! H2O
    cpms(4) = cpms(4) * 8.248767324424338d+07 ! H
    cpms(5) = cpms(5) * 5.196763628636074d+06 ! O
    cpms(6) = cpms(6) * 4.888768810227566d+06 ! OH
    cpms(7) = cpms(7) * 2.519031701678171d+06 ! HO2
    cpms(8) = cpms(8) * 2.444384405113783d+06 ! H2O2
    cpms(9) = cpms(9) * 2.968047434442088d+06 ! N2

end subroutine

! Returns internal energy in mass units (Eq 30.)
subroutine ckums(T, ums)

    !$acc routine(ckums) seq
    !$acc routine(speciesinternalenergy) seq

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout) :: ums(9)

    double precision :: tT, tc(5)
    double precision :: RT
    integer :: i
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesInternalEnergy(ums, tc)

    do i=1, 9
        ums(i) = ums(i) * (RT * imw(i))
    end do

end subroutine ckums

subroutine ckhms(T, hms)

    !$acc routine(ckhms) seq
    !$acc routine(speciesEnthalpy) seq

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout) :: hms(9)

    double precision :: tT, RT
    double precision :: tc(5), h(9)
    integer :: i
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesEnthalpy(hms, tc)

    do i=1, 9
        hms(i) = hms(i) * (RT * imw(i))
    end do

end subroutine

! Returns enthalpy in mass units (Eq 27.)
subroutine ckhms_gpu(q, hii, lo, hi, i, j, k, qvar, qtemp, qfs, nspec)

    !$acc routine(ckhms_gpu) seq
    !$acc routine(speciesEnthalpy_gpu) seq

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: i, j, k, qvar, qtemp, qfs, nspec
    double precision, intent(in), dimension(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:qvar) :: q
    double precision, intent(out), dimension(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nspec) :: hii

    double precision :: tT, RT
    double precision :: tc(5), h(9)
    integer :: n
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    tT = q(i,j,k,qtemp) ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesEnthalpy_gpu(hii, lo, hi, i, j, k, tc, nspec)

    do n=1, 9
        hii(i,j,k,n) = hii(i,j,k,n) * (RT * imw(n))
    end do

end subroutine

! Returns enthalpy in mass units (Eq 27.)
subroutine vckhms(np, T, hms)

    implicit none

    integer, intent(in) :: np
    double precision, intent(in) :: T(np)
    double precision, intent(inout) :: hms(np,9)

    double precision :: tc(5), h(9)
    integer :: i, n
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    do i=1, np
        tc(1) = 0.d0
        tc(2) = T(i)
        tc(3) = T(i)*T(i)
        tc(4) = T(i)*T(i)*T(i)
        tc(5) = T(i)*T(i)*T(i)*T(i)

        call speciesEnthalpy(h, tc)

        hms(i, 1) = h(1)
        hms(i, 2) = h(2)
        hms(i, 3) = h(3)
        hms(i, 4) = h(4)
        hms(i, 5) = h(5)
        hms(i, 6) = h(6)
        hms(i, 7) = h(7)
        hms(i, 8) = h(8)
        hms(i, 9) = h(9)
    end do

    do n=1, 9
        do i=1, np
            hms(i,n) = hms(i,n) * ( 8.31451000d+07 * T(i) * imw(n))
        end do
    end do

end subroutine

! Returns the mean specific heat at CP (Eq. 34)
subroutine ckcpbs(T, y, cpbs)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    double precision, intent(out) :: cpbs

    double precision :: cpor(9)
    double precision :: tresult(9)
    double precision :: tT, tc(5)
    double precision :: res
    integer :: i
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    res = 0.d0

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cp_R(cpor, tc)

    do i=1, 9
        tresult(i) = cpor(i) * y(i) * imw(i)
    end do
    do i=1, 9
        res = res + tresult(i)
    end do

    cpbs = res * 8.31451000d+07

end subroutine

! Returns the mean specific heat at CV (Eq. 36)
subroutine ckcvbs(T, y, cvbs)

    !$acc routine(ckcvbs) seq

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    double precision, intent(out) :: cvbs

    double precision :: cvor(9)
    double precision :: tT, tc(5)
    double precision :: res
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    res = 0.d0
    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cv_R(cvor, tc)

    ! multiply by y/molecularweight
    res = res + (cvor(1) * y(1) * imw(1)) ! H2
    res = res + (cvor(2) * y(2) * imw(2)) ! O2
    res = res + (cvor(3) * y(3) * imw(3)) ! H2O
    res = res + (cvor(4) * y(4) * imw(4)) ! H
    res = res + (cvor(5) * y(5) * imw(5)) ! O
    res = res + (cvor(6) * y(6) * imw(6)) ! OH
    res = res + (cvor(7) * y(7) * imw(7)) ! HO2
    res = res + (cvor(8) * y(8) * imw(8)) ! H2O2
    res = res + (cvor(9) * y(9) * imw(9)) ! N2

    cvbs = res * 8.31451000d+07

end subroutine

! get mean internal energy in mass units
subroutine ckubms(T, y, ubms)

    !$acc routine(ckubms) seq
    !$acc routine(speciesinternalenergy) seq

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    double precision, intent(out) :: ubms

    double precision :: ums(9) ! temporary energy array
    double precision :: res
    double precision :: RT, tT, tc(5)
    double precision, parameter :: imw(9) = (/ &
        1.d0 / 2.015940d0,  & ! H2
        1.d0 / 31.998800d0,  & ! O2
        1.d0 / 18.015340d0,  & ! H2O
        1.d0 / 1.007970d0,  & ! H
        1.d0 / 15.999400d0,  & ! O
        1.d0 / 17.007370d0,  & ! OH
        1.d0 / 33.006770d0,  & ! HO2
        1.d0 / 34.014740d0,  & ! H2O2
        1.d0 / 28.013400d0 /) ! N2

    res = 0.d0

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesInternalEnergy(ums, tc)

    ! perform dot product + scaling by wt
    res = res + (y(1) * ums(1) * imw(1)) ! H2
    res = res + (y(2) * ums(2) * imw(2)) ! O2
    res = res + (y(3) * ums(3) * imw(3)) ! H2O
    res = res + (y(4) * ums(4) * imw(4)) ! H
    res = res + (y(5) * ums(5) * imw(5)) ! O
    res = res + (y(6) * ums(6) * imw(6)) ! OH
    res = res + (y(7) * ums(7) * imw(7)) ! HO2
    res = res + (y(8) * ums(8) * imw(8)) ! H2O2
    res = res + (y(9) * ums(9) * imw(9)) ! N2

    ubms = res * RT

end subroutine

! compute the production rate for each species
subroutine ckwc(T, C, wdot)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout) :: C(9)
    double precision, intent(inout) :: wdot(9)

    integer :: id ! loop counter

    ! convert to SI
    do id = 1, 9
        C(id) = C(id) * 1.0d6
    end do

    ! convert to chemkin units
    call productionRate(wdot, C, T)

    ! convert to chemkin units
    do id=1, 9
        C(id) = C(id) * 1.0d-6
        wdot(id) = wdot(id) * 1.0d-6
    end do

end subroutine

! compute the production rate for each species
subroutine productionRate(wdot, sc, T)

    implicit none

    double precision, intent(inout) :: wdot(9)
    double precision, intent(in) :: sc(9)
    double precision, intent(in) :: T

    double precision :: tc(5)
    double precision :: invT
    double precision :: qdot, q_f(21), q_r(21)
    integer :: i

    tc = (/ log(T), T, T*T, T*T*T, T*T*T*T /)
    invT = 1.d0 / tc(2)

    if (T /= T_save) then
        T_save = T
        call comp_k_f(tc,invT,k_f_save)
        call comp_Kc(tc,invT,Kc_save)
    end if

    call comp_qfqr(q_f, q_r, sc, tc, invT)

    do i=1, 9
        wdot(i) = 0.d0
    end do

    qdot = q_f(1)-q_r(1)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot

    qdot = q_f(2)-q_r(2)
    wdot(6) = wdot(6) + qdot
    wdot(6) = wdot(6) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(3)-q_r(3)
    wdot(1) = wdot(1) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(4) = wdot(4) + qdot

    qdot = q_f(4)-q_r(4)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(5) = wdot(5) - qdot

    qdot = q_f(5)-q_r(5)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(6)-q_r(6)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(6) = wdot(6) - qdot

    qdot = q_f(7)-q_r(7)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(8)-q_r(8)
    wdot(1) = wdot(1) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(9)-q_r(9)
    wdot(1) = wdot(1) - qdot
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) + qdot
    wdot(6) = wdot(6) - qdot

    qdot = q_f(10)-q_r(10)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(11)-q_r(11)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(12)-q_r(12)
    wdot(4) = wdot(4) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(13)-q_r(13)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(14)-q_r(14)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(15)-q_r(15)
    wdot(2) = wdot(2) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(7) = wdot(7) - qdot
    wdot(8) = wdot(8) + qdot

    qdot = q_f(16)-q_r(16)
    wdot(2) = wdot(2) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(7) = wdot(7) - qdot
    wdot(8) = wdot(8) + qdot

    qdot = q_f(17)-q_r(17)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(18)-q_r(18)
    wdot(1) = wdot(1) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(19)-q_r(19)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(20)-q_r(20)
    wdot(3) = wdot(3) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(21)-q_r(21)
    wdot(3) = wdot(3) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

end subroutine

subroutine comp_k_f(tc, invT, k_f)

    implicit none

    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT
    double precision, intent(out) :: k_f(21)

    integer :: i

    do i=1, 21
        k_f(i) = prefactor_units(i)*fwd_A(i)*exp(fwd_beta(i)*tc(1)-activation_units(i)*fwd_Ea(i)*invT)
    end do

end subroutine

subroutine comp_Kc(tc, invT, Kc)

    implicit none

    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT
    double precision, intent(inout) :: Kc(21)

    double precision :: g_RT(9)
    double precision :: refC, refCinv
    integer :: i

    ! compute the Gibbs free energy
    call gibbs(g_RT, tc)

    Kc(1) = g_RT(2) + g_RT(4) - g_RT(7)
    Kc(2) = -g_RT(6) - g_RT(6) + g_RT(8)
    Kc(3) = g_RT(1) - g_RT(4) - g_RT(4)
    Kc(4) = -g_RT(2) + g_RT(5) + g_RT(5)
    Kc(5) = g_RT(4) + g_RT(5) - g_RT(6)
    Kc(6) = -g_RT(3) + g_RT(4) + g_RT(6)
    Kc(7) = g_RT(2) + g_RT(4) - g_RT(5) - g_RT(6)
    Kc(8) = g_RT(1) - g_RT(4) + g_RT(5) - g_RT(6)
    Kc(9) = g_RT(1) - g_RT(3) - g_RT(4) + g_RT(6)
    Kc(10) = g_RT(3) + g_RT(5) - g_RT(6) - g_RT(6)
    Kc(11) = -g_RT(1) - g_RT(2) + g_RT(4) + g_RT(7)
    Kc(12) = g_RT(4) - g_RT(6) - g_RT(6) + g_RT(7)
    Kc(13) = -g_RT(2) + g_RT(5) - g_RT(6) + g_RT(7)
    Kc(14) = -g_RT(2) - g_RT(3) + g_RT(6) + g_RT(7)
    Kc(15) = -g_RT(2) + g_RT(7) + g_RT(7) - g_RT(8)
    Kc(16) = -g_RT(2) + g_RT(7) + g_RT(7) - g_RT(8)
    Kc(17) = -g_RT(3) + g_RT(4) - g_RT(6) + g_RT(8)
    Kc(18) = -g_RT(1) + g_RT(4) - g_RT(7) + g_RT(8)
    Kc(19) = g_RT(5) - g_RT(6) - g_RT(7) + g_RT(8)
    Kc(20) = -g_RT(3) + g_RT(6) - g_RT(7) + g_RT(8)
    Kc(21) = -g_RT(3) + g_RT(6) - g_RT(7) + g_RT(8)

    do i=1, 21
        Kc(i) = exp(Kc(i))
    end do

    ! reference concentration: P_atm / (RT) in inverse mol/m^3
    refC = 101325d0 / 8.31451d0 * invT
    refCinv = 1.d0 / refC

    Kc(1) = Kc(1) * (refCinv)
    Kc(2) = Kc(2) * (refC)
    Kc(3) = Kc(3) * (refC)
    Kc(4) = Kc(4) * (refCinv)
    Kc(5) = Kc(5) * (refCinv)
    Kc(6) = Kc(6) * (refCinv)

end subroutine

subroutine comp_qfqr(qf, qr, sc, tc, invT)

    implicit none

    double precision, intent(out) :: qf(21)
    double precision, intent(out) :: qr(21)
    double precision, intent(in) :: sc(9)
    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT

    double precision :: T
    double precision :: mixture
    double precision :: Corr(21)
    double precision :: alpha_troe(2)
    double precision :: redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe
    double precision :: tmp1, tmp2, tmp3
    integer :: i

    ! reaction 1: H + O2 (+M) <=> HO2 (+M)
    qf(1) = sc(2)*sc(4)
    qr(1) = sc(7)
    ! reaction 2: H2O2 (+M) <=> OH + OH (+M)
    qf(2) = sc(8)
    qr(2) = sc(6)*sc(6)
    ! reaction 3: H2 + M <=> H + H + M
    qf(3) = sc(1)
    qr(3) = sc(4)*sc(4)
    ! reaction 4: O + O + M <=> O2 + M
    qf(4) = sc(5)*sc(5)
    qr(4) = sc(2)
    ! reaction 5: O + H + M <=> OH + M
    qf(5) = sc(4)*sc(5)
    qr(5) = sc(6)
    ! reaction 6: H + OH + M <=> H2O + M
    qf(6) = sc(4)*sc(6)
    qr(6) = sc(3)
    ! reaction 7: H + O2 <=> O + OH
    qf(7) = sc(2)*sc(4)
    qr(7) = sc(5)*sc(6)
    ! reaction 8: O + H2 <=> H + OH
    qf(8) = sc(1)*sc(5)
    qr(8) = sc(4)*sc(6)
    ! reaction 9: H2 + OH <=> H2O + H
    qf(9) = sc(1)*sc(6)
    qr(9) = sc(3)*sc(4)
    ! reaction 10: O + H2O <=> OH + OH
    qf(10) = sc(3)*sc(5)
    qr(10) = sc(6)*sc(6)
    ! reaction 11: HO2 + H <=> H2 + O2
    qf(11) = sc(4)*sc(7)
    qr(11) = sc(1)*sc(2)
    ! reaction 12: HO2 + H <=> OH + OH
    qf(12) = sc(4)*sc(7)
    qr(12) = sc(6)*sc(6)
    ! reaction 13: HO2 + O <=> O2 + OH
    qf(13) = sc(5)*sc(7)
    qr(13) = sc(2)*sc(6)
    ! reaction 14: HO2 + OH <=> H2O + O2
    qf(14) = sc(6)*sc(7)
    qr(14) = sc(2)*sc(3)
    ! reaction 15: HO2 + HO2 <=> H2O2 + O2
    qf(15) = sc(7)*sc(7)
    qr(15) = sc(2)*sc(8)
    ! reaction 16: HO2 + HO2 <=> H2O2 + O2
    qf(16) = sc(7)*sc(7)
    qr(16) = sc(2)*sc(8)
    ! reaction 17: H2O2 + H <=> H2O + OH
    qf(17) = sc(4)*sc(8)
    qr(17) = sc(3)*sc(6)
    ! reaction 18: H2O2 + H <=> HO2 + H2
    qf(18) = sc(4)*sc(8)
    qr(18) = sc(1)*sc(7)
    ! reaction 19: H2O2 + O <=> OH + HO2
    qf(19) = sc(5)*sc(8)
    qr(19) = sc(6)*sc(7)
    ! reaction 20: H2O2 + OH <=> HO2 + H2O
    qf(20) = sc(6)*sc(8)
    qr(20) = sc(3)*sc(7)
    ! reaction 21: H2O2 + OH <=> HO2 + H2O
    qf(21) = sc(6)*sc(8)
    qr(21) = sc(3)*sc(7)

    T = tc(2)

    ! compute the mixture concentration
    mixture = 0.d0
    do i=1, 9
        mixture = mixture + sc(i)
    end do

    do i=1, 21
        Corr(i) = 1.d0
    end do

    ! troe
    alpha_troe(1) = mixture + (TB(1) % vector(1) - 1)*sc(1) + (TB(1) % vector(2) - 1)*sc(3) + (TB(1) % vector(3) - 1)*sc(2)
    alpha_troe(2) = mixture + (TB(2) % vector(1) - 1)*sc(1) + (TB(2) % vector(2) - 1)*sc(3)

    do i=1, 2
        redP = alpha_troe(i-0) / k_f_save(i) * phase_units(i) * low_A(i) * exp(low_beta(i) * tc(1) - activation_units(i) * low_Ea(i) *invT)
        F = redP / (1.d0 + redP)
        logPred = log10(redP)

        if (abs(troe_Tsss(i)) > 1.d-100) then
           tmp1 = (1.d0-troe_a(i))*exp(-T/troe_Tsss(i))
        else
           tmp1 = 0.d0
        end if

        if (abs(troe_Ts(i)) > 1.d-100) then
           tmp2 = troe_a(i) * exp(-T/troe_Ts(i))
        else
           tmp2 = 0.d0
        end if

        if (troe_len(i) == 4) then
           tmp3 = exp(-troe_Tss(i) * invT)
        else
           tmp3 = 0.d0
        end if

        logFcent = log10(tmp1+tmp2+tmp3)
        troe_c = -0.4d0 - 0.67d0 * logFcent
        troe_n = 0.75d0 - 1.27d0 * logFcent
        troe = (troe_c + logPred) / (troe_n - 0.14d0*(troe_c + logPred))
        F_troe = 10.d0 ** (logFcent / (1.d0 + troe*troe))
        Corr(i) = F * F_troe
    end do

    ! simple three-body correction
    Corr(3) = mixture + (TB(3) % vector(1) - 1)*sc(1) + (TB(3) % vector(2) - 1)*sc(3)
    Corr(4) = mixture + (TB(4) % vector(1) - 1)*sc(1) + (TB(4) % vector(2) - 1)*sc(3)
    Corr(5) = mixture + (TB(5) % vector(1) - 1)*sc(1) + (TB(5) % vector(2) - 1)*sc(3)
    Corr(6) = mixture + (TB(6) % vector(1) - 1)*sc(1) + (TB(6) % vector(2) - 1)*sc(3)

    do i=1, 21
        qf(i) = qf(i) * (Corr(i) * k_f_save(i))
        qr(i) = qr(i) * (Corr(i) * k_f_save(i) / Kc_save(i))
    end do

end subroutine

! compute the g/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine gibbs(species, tc)

    implicit none

    double precision, intent(out) :: species(9)
    double precision, intent(in) :: tc(5)

    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            -1.012520870000000d+03 * invT &
            +6.592218400000000d+00 &
            -3.298124310000000d+00 * tc(1) &
            -4.124720870000000d-04 * tc(2) &
            +1.357169215000000d-07 * tc(3) &
            +7.896195275000000d-12 * tc(4) &
            -2.067436120000000d-14 * tc(5)
        ! species 2: O2
        species(2) = &
            -1.005249020000000d+03 * invT &
            -2.821801190000000d+00 &
            -3.212936400000000d+00 * tc(1) &
            -5.637431750000000d-04 * tc(2) &
            +9.593584116666666d-08 * tc(3) &
            -1.094897691666667d-10 * tc(4) &
            +4.384276960000000d-14 * tc(5)
        ! species 3: H2O
        species(3) = &
            -3.020811330000000d+04 * invT &
            +7.966096399999998d-01 &
            -3.386842490000000d+00 * tc(1) &
            -1.737491230000000d-03 * tc(2) &
            +1.059116055000000d-06 * tc(3) &
            -5.807151058333333d-10 * tc(4) &
            +1.253294235000000d-13 * tc(5)
        ! species 4: H
        species(4) = &
            +2.547162700000000d+04 * invT &
            +2.960117608000000d+00 &
            -2.500000000000000d+00 * tc(1) &
            -0.000000000000000d+00 * tc(2) &
            -0.000000000000000d+00 * tc(3) &
            -0.000000000000000d+00 * tc(4) &
            -0.000000000000000d+00 * tc(5)
        ! species 5: O
        species(5) = &
            +2.914764450000000d+04 * invT &
            -1.756619999999964d-02 &
            -2.946428780000000d+00 * tc(1) &
            +8.190832450000000d-04 * tc(2) &
            -4.035052833333333d-07 * tc(3) &
            +1.335702658333333d-10 * tc(4) &
            -1.945348180000000d-14 * tc(5)
        ! species 6: OH
        species(6) = &
            +3.346309130000000d+03 * invT &
            +4.815738570000000d+00 &
            -4.125305610000000d+00 * tc(1) &
            +1.612724695000000d-03 * tc(2) &
            -1.087941151666667d-06 * tc(3) &
            +4.832113691666666d-10 * tc(4) &
            -1.031186895000000d-13 * tc(5)
        ! species 7: HO2
        species(7) = &
            +2.948080400000000d+02 * invT &
            +5.851355599999999d-01 &
            -4.301798010000000d+00 * tc(1) &
            +2.374560255000000d-03 * tc(2) &
            -3.526381516666666d-06 * tc(3) &
            +2.023032450000000d-09 * tc(4) &
            -4.646125620000001d-13 * tc(5)
        ! species 8: H2O2
        species(8) = &
            -1.766314650000000d+04 * invT &
            -3.396609550000000d+00 &
            -3.388753650000000d+00 * tc(1) &
            -3.284612905000000d-03 * tc(2) &
            +2.475020966666667d-08 * tc(3) &
            +3.854837933333333d-10 * tc(4) &
            -1.235757375000000d-13 * tc(5)
        ! species 9: N2
        species(9) = &
            -1.020900000000000d+03 * invT &
            -6.516950000000001d-01 &
            -3.298677000000000d+00 * tc(1) &
            -7.041200000000000d-04 * tc(2) &
            +6.605369999999999d-07 * tc(3) &
            -4.701262500000001d-10 * tc(4) &
            +1.222427500000000d-13 * tc(5)
    else
        !species 1: H2
        species(1) = &
            -8.350339970000000d+02 * invT &
            +4.346533540000000d+00 &
            -2.991423370000000d+00 * tc(1) &
            -3.500322055000000d-04 * tc(2) &
            +9.389714483333333d-09 * tc(3) &
            +7.692981816666667d-13 * tc(4) &
            -7.913758950000000d-17 * tc(5)
        !species 2: O2
        species(2) = &
            -1.233930180000000d+03 * invT &
            +5.084126000000002d-01 &
            -3.697578190000000d+00 * tc(1) &
            -3.067598445000000d-04 * tc(2) &
            +2.098069983333333d-08 * tc(3) &
            -1.479401233333333d-12 * tc(4) &
            +5.682176550000000d-17 * tc(5)
        !species 3: H2O
        species(3) = &
            -2.989920900000000d+04 * invT &
            -4.190671200000001d+00 &
            -2.672145610000000d+00 * tc(1) &
            -1.528146445000000d-03 * tc(2) &
            +1.455043351666667d-07 * tc(3) &
            -1.000830325000000d-11 * tc(4) &
            +3.195808935000000d-16 * tc(5)
        !species 4: H
        species(4) = &
            +2.547162700000000d+04 * invT &
            +2.960117638000000d+00 &
            -2.500000000000000d+00 * tc(1) &
            -0.000000000000000d+00 * tc(2) &
            -0.000000000000000d+00 * tc(3) &
            -0.000000000000000d+00 * tc(4) &
            -0.000000000000000d+00 * tc(5)
        !species 5: O
        species(5) = &
            +2.923080270000000d+04 * invT &
            -2.378248450000000d+00 &
            -2.542059660000000d+00 * tc(1) &
            +1.377530955000000d-05 * tc(2) &
            +5.171338916666667d-10 * tc(3) &
            -3.792556183333333d-13 * tc(4) &
            +2.184025750000000d-17 * tc(5)
        !species 6: OH
        species(6) = &
            +3.683628750000000d+03 * invT &
            -2.836911870000000d+00 &
            -2.864728860000000d+00 * tc(1) &
            -5.282522400000000d-04 * tc(2) &
            +4.318045966666667d-08 * tc(3) &
            -2.543488950000000d-12 * tc(4) &
            +6.659793800000000d-17 * tc(5)
        !species 7: HO2
        species(7) = &
            +1.118567130000000d+02 * invT &
            +2.321087500000001d-01 &
            -4.017210900000000d+00 * tc(1) &
            -1.119910065000000d-03 * tc(2) &
            +1.056096916666667d-07 * tc(3) &
            -9.520530833333334d-12 * tc(4) &
            +5.395426750000000d-16 * tc(5)
        !species 8: H2O2
        species(8) = &
            -1.800696090000000d+04 * invT &
            +4.072029891000000d+00 &
            -4.573166850000000d+00 * tc(1) &
            -2.168068195000000d-03 * tc(2) &
            +2.457814700000000d-07 * tc(3) &
            -1.957419641666667d-11 * tc(4) &
            +7.158267800000000d-16 * tc(5)
        !species 9: N2
        species(9) = &
            -9.227977000000000d+02 * invT &
            -3.053888000000000d+00 &
            -2.926640000000000d+00 * tc(1) &
            -7.439885000000000d-04 * tc(2) &
            +9.474601666666666d-08 * tc(3) &
            -8.414199999999999d-12 * tc(4) &
            +3.376675500000000d-16 * tc(5)
    end if

end subroutine


! compute Cv/R at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine cv_R(species, tc)

    !$acc routine(cv_r) seq

    implicit none

    double precision, intent(out) :: species(9)
    double precision, intent(in) :: tc(5)

    double precision :: T

    T = tc(2)

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +2.29812431d+00 &
            +8.24944174d-04 * tc(2) &
            -8.14301529d-07 * tc(3) &
            -9.47543433d-11 * tc(4) &
            +4.13487224d-13 * tc(5)
        ! species 2: O2
        species(2) = &
            +2.21293640d+00 &
            +1.12748635d-03 * tc(2) &
            -5.75615047d-07 * tc(3) &
            +1.31387723d-09 * tc(4) &
            -8.76855392d-13 * tc(5)
        ! species 3: H2O
        species(3) = &
            +2.38684249d+00 &
            +3.47498246d-03 * tc(2) &
            -6.35469633d-06 * tc(3) &
            +6.96858127d-09 * tc(4) &
            -2.50658847d-12 * tc(5)
        ! species 4: H
        species(4) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
        ! species 5: O
        species(5) = &
            +1.94642878d+00 &
            -1.63816649d-03 * tc(2) &
            +2.42103170d-06 * tc(3) &
            -1.60284319d-09 * tc(4) &
            +3.89069636d-13 * tc(5)
        ! species 6: OH
        species(6) = &
            +3.12530561d+00 &
            -3.22544939d-03 * tc(2) &
            +6.52764691d-06 * tc(3) &
            -5.79853643d-09 * tc(4) &
            +2.06237379d-12 * tc(5)
        ! species 7: HO2
        species(7) = &
            +3.30179801d+00 &
            -4.74912051d-03 * tc(2) &
            +2.11582891d-05 * tc(3) &
            -2.42763894d-08 * tc(4) &
            +9.29225124d-12 * tc(5)
        ! species 8: H2O2
        species(8) = &
            +2.38875365d+00 &
            +6.56922581d-03 * tc(2) &
            -1.48501258d-07 * tc(3) &
            -4.62580552d-09 * tc(4) &
            +2.47151475d-12 * tc(5)
        ! species 9: N2
        species(9) = &
            +2.29867700d+00 &
            +1.40824000d-03 * tc(2) &
            -3.96322200d-06 * tc(3) &
            +5.64151500d-09 * tc(4) &
            -2.44485500d-12 * tc(5)
    else
        !species 1: H2
        species(1) = &
            +1.99142337d+00 &
            +7.00064411d-04 * tc(2) &
            -5.63382869d-08 * tc(3) &
            -9.23157818d-12 * tc(4) &
            +1.58275179d-15 * tc(5)
        !species 2: O2
        species(2) = &
            +2.69757819d+00 &
            +6.13519689d-04 * tc(2) &
            -1.25884199d-07 * tc(3) &
            +1.77528148d-11 * tc(4) &
            -1.13643531d-15 * tc(5)
        !species 3: H2O
        species(3) = &
            +1.67214561d+00 &
            +3.05629289d-03 * tc(2) &
            -8.73026011d-07 * tc(3) &
            +1.20099639d-10 * tc(4) &
            -6.39161787d-15 * tc(5)
        !species 4: H
        species(4) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
        !species 5: O
        species(5) = &
            +1.54205966d+00 &
            -2.75506191d-05 * tc(2) &
            -3.10280335d-09 * tc(3) &
            +4.55106742d-12 * tc(4) &
            -4.36805150d-16 * tc(5)
        !species 6: OH
        species(6) = &
            +1.86472886d+00 &
            +1.05650448d-03 * tc(2) &
            -2.59082758d-07 * tc(3) &
            +3.05218674d-11 * tc(4) &
            -1.33195876d-15 * tc(5)
        !species 7: HO2
        species(7) = &
            +3.01721090d+00 &
            +2.23982013d-03 * tc(2) &
            -6.33658150d-07 * tc(3) &
            +1.14246370d-10 * tc(4) &
            -1.07908535d-14 * tc(5)
        !species 8: H2O2
        species(8) = &
            +3.57316685d+00 &
            +4.33613639d-03 * tc(2) &
            -1.47468882d-06 * tc(3) &
            +2.34890357d-10 * tc(4) &
            -1.43165356d-14 * tc(5)
        !species 9: N2
        species(9) = &
            +1.92664000d+00 &
            +1.48797700d-03 * tc(2) &
            -5.68476100d-07 * tc(3) &
            +1.00970400d-10 * tc(4) &
            -6.75335100d-15 * tc(5)
    end if

end subroutine


! compute Cp/R at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine cp_R(species, tc)

    !$acc routine(cp_r) seq

    implicit none

    double precision, intent(out) :: species(9)
    double precision, intent(in) :: tc(5)

    double precision :: T

    T = tc(2)

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +3.29812431d+00 &
            +8.24944174d-04 * tc(2) &
            -8.14301529d-07 * tc(3) &
            -9.47543433d-11 * tc(4) &
            +4.13487224d-13 * tc(5)
        ! species 2: O2
        species(2) = &
            +3.21293640d+00 &
            +1.12748635d-03 * tc(2) &
            -5.75615047d-07 * tc(3) &
            +1.31387723d-09 * tc(4) &
            -8.76855392d-13 * tc(5)
        ! species 3: H2O
        species(3) = &
            +3.38684249d+00 &
            +3.47498246d-03 * tc(2) &
            -6.35469633d-06 * tc(3) &
            +6.96858127d-09 * tc(4) &
            -2.50658847d-12 * tc(5)
        ! species 4: H
        species(4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
        ! species 5: O
        species(5) = &
            +2.94642878d+00 &
            -1.63816649d-03 * tc(2) &
            +2.42103170d-06 * tc(3) &
            -1.60284319d-09 * tc(4) &
            +3.89069636d-13 * tc(5)
        ! species 6: OH
        species(6) = &
            +4.12530561d+00 &
            -3.22544939d-03 * tc(2) &
            +6.52764691d-06 * tc(3) &
            -5.79853643d-09 * tc(4) &
            +2.06237379d-12 * tc(5)
        ! species 7: HO2
        species(7) = &
            +4.30179801d+00 &
            -4.74912051d-03 * tc(2) &
            +2.11582891d-05 * tc(3) &
            -2.42763894d-08 * tc(4) &
            +9.29225124d-12 * tc(5)
        ! species 8: H2O2
        species(8) = &
            +3.38875365d+00 &
            +6.56922581d-03 * tc(2) &
            -1.48501258d-07 * tc(3) &
            -4.62580552d-09 * tc(4) &
            +2.47151475d-12 * tc(5)
        ! species 9: N2
        species(9) = &
            +3.29867700d+00 &
            +1.40824000d-03 * tc(2) &
            -3.96322200d-06 * tc(3) &
            +5.64151500d-09 * tc(4) &
            -2.44485500d-12 * tc(5)
    else
        !species 1: H2
        species(1) = &
            +2.99142337d+00 &
            +7.00064411d-04 * tc(2) &
            -5.63382869d-08 * tc(3) &
            -9.23157818d-12 * tc(4) &
            +1.58275179d-15 * tc(5)
        !species 2: O2
        species(2) = &
            +3.69757819d+00 &
            +6.13519689d-04 * tc(2) &
            -1.25884199d-07 * tc(3) &
            +1.77528148d-11 * tc(4) &
            -1.13643531d-15 * tc(5)
        !species 3: H2O
        species(3) = &
            +2.67214561d+00 &
            +3.05629289d-03 * tc(2) &
            -8.73026011d-07 * tc(3) &
            +1.20099639d-10 * tc(4) &
            -6.39161787d-15 * tc(5)
        !species 4: H
        species(4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
        !species 5: O
        species(5) = &
            +2.54205966d+00 &
            -2.75506191d-05 * tc(2) &
            -3.10280335d-09 * tc(3) &
            +4.55106742d-12 * tc(4) &
            -4.36805150d-16 * tc(5)
        !species 6: OH
        species(6) = &
            +2.86472886d+00 &
            +1.05650448d-03 * tc(2) &
            -2.59082758d-07 * tc(3) &
            +3.05218674d-11 * tc(4) &
            -1.33195876d-15 * tc(5)
        !species 7: HO2
        species(7) = &
            +4.01721090d+00 &
            +2.23982013d-03 * tc(2) &
            -6.33658150d-07 * tc(3) &
            +1.14246370d-10 * tc(4) &
            -1.07908535d-14 * tc(5)
        !species 8: H2O2
        species(8) = &
            +4.57316685d+00 &
            +4.33613639d-03 * tc(2) &
            -1.47468882d-06 * tc(3) &
            +2.34890357d-10 * tc(4) &
            -1.43165356d-14 * tc(5)
        !species 9: N2
        species(9) = &
            +2.92664000d+00 &
            +1.48797700d-03 * tc(2) &
            -5.68476100d-07 * tc(3) &
            +1.00970400d-10 * tc(4) &
            -6.75335100d-15 * tc(5)
    end if

end subroutine

! compute the e/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine speciesInternalEnergy(species, tc)

    !$acc routine(speciesinternalenergy) seq

    implicit none

    double precision, intent(out) :: species(9)
    double precision, intent(in) :: tc(5)

    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +2.29812431d+00 &
            +4.12472087d-04 * tc(2) &
            -2.71433843d-07 * tc(3) &
            -2.36885858d-11 * tc(4) &
            +8.26974448d-14 * tc(5) &
            -1.01252087d+03 * invT
        ! species 2: O2
        species(2) = &
            +2.21293640d+00 &
            +5.63743175d-04 * tc(2) &
            -1.91871682d-07 * tc(3) &
            +3.28469308d-10 * tc(4) &
            -1.75371078d-13 * tc(5) &
            -1.00524902d+03 * invT
        ! species 3: H2O
        species(3) = &
            +2.38684249d+00 &
            +1.73749123d-03 * tc(2) &
            -2.11823211d-06 * tc(3) &
            +1.74214532d-09 * tc(4) &
            -5.01317694d-13 * tc(5) &
            -3.02081133d+04 * invT
        ! species 4: H
        species(4) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            +2.54716270d+04 * invT
        ! species 5: O
        species(5) = &
            +1.94642878d+00 &
            -8.19083245d-04 * tc(2) &
            +8.07010567d-07 * tc(3) &
            -4.00710797d-10 * tc(4) &
            +7.78139272d-14 * tc(5) &
            +2.91476445d+04 * invT
        ! species 6: OH
        species(6) = &
            +3.12530561d+00 &
            -1.61272470d-03 * tc(2) &
            +2.17588230d-06 * tc(3) &
            -1.44963411d-09 * tc(4) &
            +4.12474758d-13 * tc(5) &
            +3.34630913d+03 * invT
        ! species 7: HO2
        species(7) = &
            +3.30179801d+00 &
            -2.37456025d-03 * tc(2) &
            +7.05276303d-06 * tc(3) &
            -6.06909735d-09 * tc(4) &
            +1.85845025d-12 * tc(5) &
            +2.94808040d+02 * invT
        ! species 8: H2O2
        species(8) = &
            +2.38875365d+00 &
            +3.28461290d-03 * tc(2) &
            -4.95004193d-08 * tc(3) &
            -1.15645138d-09 * tc(4) &
            +4.94302950d-13 * tc(5) &
            -1.76631465d+04 * invT
        ! species 9: N2
        species(9) = &
            +2.29867700d+00 &
            +7.04120000d-04 * tc(2) &
            -1.32107400d-06 * tc(3) &
            +1.41037875d-09 * tc(4) &
            -4.88971000d-13 * tc(5) &
            -1.02090000d+03 * invT
    else
        !species 1: H2
        species(1) = &
            +1.99142337d+00 &
            +3.50032206d-04 * tc(2) &
            -1.87794290d-08 * tc(3) &
            -2.30789455d-12 * tc(4) &
            +3.16550358d-16 * tc(5) &
            -8.35033997d+02 * invT
        !species 2: O2
        species(2) = &
            +2.69757819d+00 &
            +3.06759845d-04 * tc(2) &
            -4.19613997d-08 * tc(3) &
            +4.43820370d-12 * tc(4) &
            -2.27287062d-16 * tc(5) &
            -1.23393018d+03 * invT
        !species 3: H2O
        species(3) = &
            +1.67214561d+00 &
            +1.52814644d-03 * tc(2) &
            -2.91008670d-07 * tc(3) &
            +3.00249098d-11 * tc(4) &
            -1.27832357d-15 * tc(5) &
            -2.98992090d+04 * invT
        !species 4: H
        species(4) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            +2.54716270d+04 * invT
        !species 5: O
        species(5) = &
            +1.54205966d+00 &
            -1.37753096d-05 * tc(2) &
            -1.03426778d-09 * tc(3) &
            +1.13776685d-12 * tc(4) &
            -8.73610300d-17 * tc(5) &
            +2.92308027d+04 * invT
        !species 6: OH
        species(6) = &
            +1.86472886d+00 &
            +5.28252240d-04 * tc(2) &
            -8.63609193d-08 * tc(3) &
            +7.63046685d-12 * tc(4) &
            -2.66391752d-16 * tc(5) &
            +3.68362875d+03 * invT
        !species 7: HO2
        species(7) = &
            +3.01721090d+00 &
            +1.11991006d-03 * tc(2) &
            -2.11219383d-07 * tc(3) &
            +2.85615925d-11 * tc(4) &
            -2.15817070d-15 * tc(5) &
            +1.11856713d+02 * invT
        !species 8: H2O2
        species(8) = &
            +3.57316685d+00 &
            +2.16806820d-03 * tc(2) &
            -4.91562940d-07 * tc(3) &
            +5.87225893d-11 * tc(4) &
            -2.86330712d-15 * tc(5) &
            -1.80069609d+04 * invT
        !species 9: N2
        species(9) = &
            +1.92664000d+00 &
            +7.43988500d-04 * tc(2) &
            -1.89492033d-07 * tc(3) &
            +2.52426000d-11 * tc(4) &
            -1.35067020d-15 * tc(5) &
            -9.22797700d+02 * invT
    end if

end subroutine

! compute the h/(RT) at the given temperature (Eq 20)
! tc contains precomputed powers of T, tc(1) = log(T)
subroutine speciesEnthalpy(species, tc)

    !$acc routine(speciesEnthalpy) seq

    implicit none

    double precision, intent(out) :: species(9)
    double precision, intent(in) :: tc(5)

    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +3.29812431d+00 &
            +4.12472087d-04 * tc(2) &
            -2.71433843d-07 * tc(3) &
            -2.36885858d-11 * tc(4) &
            +8.26974448d-14 * tc(5) &
            -1.01252087d+03 * invT
        ! species 2: O2
        species(2) = &
            +3.21293640d+00 &
            +5.63743175d-04 * tc(2) &
            -1.91871682d-07 * tc(3) &
            +3.28469308d-10 * tc(4) &
            -1.75371078d-13 * tc(5) &
            -1.00524902d+03 * invT
        ! species 3: H2O
        species(3) = &
            +3.38684249d+00 &
            +1.73749123d-03 * tc(2) &
            -2.11823211d-06 * tc(3) &
            +1.74214532d-09 * tc(4) &
            -5.01317694d-13 * tc(5) &
            -3.02081133d+04 * invT
        ! species 4: H
        species(4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            +2.54716270d+04 * invT
        ! species 5: O
        species(5) = &
            +2.94642878d+00 &
            -8.19083245d-04 * tc(2) &
            +8.07010567d-07 * tc(3) &
            -4.00710797d-10 * tc(4) &
            +7.78139272d-14 * tc(5) &
            +2.91476445d+04 * invT
        ! species 6: OH
        species(6) = &
            +4.12530561d+00 &
            -1.61272470d-03 * tc(2) &
            +2.17588230d-06 * tc(3) &
            -1.44963411d-09 * tc(4) &
            +4.12474758d-13 * tc(5) &
            +3.34630913d+03 * invT
        ! species 7: HO2
        species(7) = &
            +4.30179801d+00 &
            -2.37456025d-03 * tc(2) &
            +7.05276303d-06 * tc(3) &
            -6.06909735d-09 * tc(4) &
            +1.85845025d-12 * tc(5) &
            +2.94808040d+02 * invT
        ! species 8: H2O2
        species(8) = &
            +3.38875365d+00 &
            +3.28461290d-03 * tc(2) &
            -4.95004193d-08 * tc(3) &
            -1.15645138d-09 * tc(4) &
            +4.94302950d-13 * tc(5) &
            -1.76631465d+04 * invT
        ! species 9: N2
        species(9) = &
            +3.29867700d+00 &
            +7.04120000d-04 * tc(2) &
            -1.32107400d-06 * tc(3) &
            +1.41037875d-09 * tc(4) &
            -4.88971000d-13 * tc(5) &
            -1.02090000d+03 * invT
    else
        !species 1: H2
        species(1) = &
            +2.99142337d+00 &
            +3.50032206d-04 * tc(2) &
            -1.87794290d-08 * tc(3) &
            -2.30789455d-12 * tc(4) &
            +3.16550358d-16 * tc(5) &
            -8.35033997d+02 * invT
        !species 2: O2
        species(2) = &
            +3.69757819d+00 &
            +3.06759845d-04 * tc(2) &
            -4.19613997d-08 * tc(3) &
            +4.43820370d-12 * tc(4) &
            -2.27287062d-16 * tc(5) &
            -1.23393018d+03 * invT
        !species 3: H2O
        species(3) = &
            +2.67214561d+00 &
            +1.52814644d-03 * tc(2) &
            -2.91008670d-07 * tc(3) &
            +3.00249098d-11 * tc(4) &
            -1.27832357d-15 * tc(5) &
            -2.98992090d+04 * invT
        !species 4: H
        species(4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            +2.54716270d+04 * invT
        !species 5: O
        species(5) = &
            +2.54205966d+00 &
            -1.37753096d-05 * tc(2) &
            -1.03426778d-09 * tc(3) &
            +1.13776685d-12 * tc(4) &
            -8.73610300d-17 * tc(5) &
            +2.92308027d+04 * invT
        !species 6: OH
        species(6) = &
            +2.86472886d+00 &
            +5.28252240d-04 * tc(2) &
            -8.63609193d-08 * tc(3) &
            +7.63046685d-12 * tc(4) &
            -2.66391752d-16 * tc(5) &
            +3.68362875d+03 * invT
        !species 7: HO2
        species(7) = &
            +4.01721090d+00 &
            +1.11991006d-03 * tc(2) &
            -2.11219383d-07 * tc(3) &
            +2.85615925d-11 * tc(4) &
            -2.15817070d-15 * tc(5) &
            +1.11856713d+02 * invT
        !species 8: H2O2
        species(8) = &
            +4.57316685d+00 &
            +2.16806820d-03 * tc(2) &
            -4.91562940d-07 * tc(3) &
            +5.87225893d-11 * tc(4) &
            -2.86330712d-15 * tc(5) &
            -1.80069609d+04 * invT
        !species 9: N2
        species(9) = &
            +2.92664000d+00 &
            +7.43988500d-04 * tc(2) &
            -1.89492033d-07 * tc(3) &
            +2.52426000d-11 * tc(4) &
            -1.35067020d-15 * tc(5) &
            -9.22797700d+02 * invT
    end if

end subroutine

! compute the h/(RT) at the given temperature (Eq 20)
! tc contains precomputed powers of T, tc(1) = log(T)
subroutine speciesEnthalpy_gpu(hii, lo, hi, i, j, k, tc, nspec)

    !$acc routine(speciesEnthalpy_gpu) seq

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    double precision, intent(in) :: tc(5)
    integer, intent(in) :: i, j, k, nspec
    double precision, intent(out), dimension(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nspec) :: hii

    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        hii(i, j, k, 1) = &
            +3.29812431d+00 &
            +4.12472087d-04 * tc(2) &
            -2.71433843d-07 * tc(3) &
            -2.36885858d-11 * tc(4) &
            +8.26974448d-14 * tc(5) &
            -1.01252087d+03 * invT
        ! species 2: O2
        hii(i, j, k, 2) = &
            +3.21293640d+00 &
            +5.63743175d-04 * tc(2) &
            -1.91871682d-07 * tc(3) &
            +3.28469308d-10 * tc(4) &
            -1.75371078d-13 * tc(5) &
            -1.00524902d+03 * invT
        ! species 3: H2O
        hii(i, j, k, 3) = &
            +3.38684249d+00 &
            +1.73749123d-03 * tc(2) &
            -2.11823211d-06 * tc(3) &
            +1.74214532d-09 * tc(4) &
            -5.01317694d-13 * tc(5) &
            -3.02081133d+04 * invT
        ! species 4: H
        hii(i, j, k, 4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            +2.54716270d+04 * invT
        ! species 5: O
        hii(i, j, k, 5) = &
            +2.94642878d+00 &
            -8.19083245d-04 * tc(2) &
            +8.07010567d-07 * tc(3) &
            -4.00710797d-10 * tc(4) &
            +7.78139272d-14 * tc(5) &
            +2.91476445d+04 * invT
        ! species 6: OH
        hii(i, j, k, 6) = &
            +4.12530561d+00 &
            -1.61272470d-03 * tc(2) &
            +2.17588230d-06 * tc(3) &
            -1.44963411d-09 * tc(4) &
            +4.12474758d-13 * tc(5) &
            +3.34630913d+03 * invT
        ! species 7: HO2
        hii(i, j, k, 7) = &
            +4.30179801d+00 &
            -2.37456025d-03 * tc(2) &
            +7.05276303d-06 * tc(3) &
            -6.06909735d-09 * tc(4) &
            +1.85845025d-12 * tc(5) &
            +2.94808040d+02 * invT
        ! species 8: H2O2
        hii(i, j, k, 8) = &
            +3.38875365d+00 &
            +3.28461290d-03 * tc(2) &
            -4.95004193d-08 * tc(3) &
            -1.15645138d-09 * tc(4) &
            +4.94302950d-13 * tc(5) &
            -1.76631465d+04 * invT
        ! species 9: N2
        hii(i, j, k, 9) = &
            +3.29867700d+00 &
            +7.04120000d-04 * tc(2) &
            -1.32107400d-06 * tc(3) &
            +1.41037875d-09 * tc(4) &
            -4.88971000d-13 * tc(5) &
            -1.02090000d+03 * invT
    else
        !species 1: H2
        hii(i, j, k, 1) = &
            +2.99142337d+00 &
            +3.50032206d-04 * tc(2) &
            -1.87794290d-08 * tc(3) &
            -2.30789455d-12 * tc(4) &
            +3.16550358d-16 * tc(5) &
            -8.35033997d+02 * invT
        !species 2: O2
        hii(i, j, k, 2) = &
            +3.69757819d+00 &
            +3.06759845d-04 * tc(2) &
            -4.19613997d-08 * tc(3) &
            +4.43820370d-12 * tc(4) &
            -2.27287062d-16 * tc(5) &
            -1.23393018d+03 * invT
        !species 3: H2O
        hii(i, j, k, 3) = &
            +2.67214561d+00 &
            +1.52814644d-03 * tc(2) &
            -2.91008670d-07 * tc(3) &
            +3.00249098d-11 * tc(4) &
            -1.27832357d-15 * tc(5) &
            -2.98992090d+04 * invT
        !species 4: H
        hii(i, j, k, 4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            +2.54716270d+04 * invT
        !species 5: O
        hii(i, j, k, 5) = &
            +2.54205966d+00 &
            -1.37753096d-05 * tc(2) &
            -1.03426778d-09 * tc(3) &
            +1.13776685d-12 * tc(4) &
            -8.73610300d-17 * tc(5) &
            +2.92308027d+04 * invT
        !species 6: OH
        hii(i, j, k, 6) = &
            +2.86472886d+00 &
            +5.28252240d-04 * tc(2) &
            -8.63609193d-08 * tc(3) &
            +7.63046685d-12 * tc(4) &
            -2.66391752d-16 * tc(5) &
            +3.68362875d+03 * invT
        !species 7: HO2
        hii(i, j, k, 7) = &
            +4.01721090d+00 &
            +1.11991006d-03 * tc(2) &
            -2.11219383d-07 * tc(3) &
            +2.85615925d-11 * tc(4) &
            -2.15817070d-15 * tc(5) &
            +1.11856713d+02 * invT
        !species 8: H2O2
        hii(i, j, k, 8) = &
            +4.57316685d+00 &
            +2.16806820d-03 * tc(2) &
            -4.91562940d-07 * tc(3) &
            +5.87225893d-11 * tc(4) &
            -2.86330712d-15 * tc(5) &
            -1.80069609d+04 * invT
        !species 9: N2
        hii(i, j, k, 9) = &
            +2.92664000d+00 &
            +7.43988500d-04 * tc(2) &
            -1.89492033d-07 * tc(3) &
            +2.52426000d-11 * tc(4) &
            -1.35067020d-15 * tc(5) &
            -9.22797700d+02 * invT
    end if

end subroutine

! save molecular weights into array
subroutine molecularWeight(wt)

    implicit none

    double precision, intent(out) :: wt(9)

    wt(1) = 2.015940d0 ! H2
    wt(2) = 31.998800d0 ! O2
    wt(3) = 18.015340d0 ! H2O
    wt(4) = 1.007970d0 ! H
    wt(5) = 15.999400d0 ! O
    wt(6) = 17.007370d0 ! OH
    wt(7) = 33.006770d0 ! HO2
    wt(8) = 34.014740d0 ! H2O2
    wt(9) = 28.013400d0 ! N2

end subroutine

! get temperature given internal energy in mass units and mass fracs
subroutine get_t_given_ey(e, y, t, ierr)

    !$acc routine(get_t_given_ey) seq
    !$acc routine(ckubms) seq
    !$acc routine(ckcvbs) seq

    implicit none

    double precision, intent(in) :: e
    double precision, intent(in) :: y(9)
    double precision, intent(out) :: t
    integer, intent(out) :: ierr

#ifdef CONVERGENCE
    integer, parameter :: maxiter = 5000
    double precision, parameter tol = 1.d-12
#else
    integer, parameter :: maxiter = 200
    double precision, parameter :: tol = 1.d-6
#endif

    double precision :: ein
    double precision, parameter :: tmin = 90.d0 ! max lower bound for thermo def
    double precision, parameter :: tmax = 4000.d0 ! min upper bound for thermo def
    double precision :: e1,emin,emax,cv,t1,dt
    integer :: i ! loop counter

    ein = e

    call ckubms(tmin, y, emin)
    call ckubms(tmax, y, emax)

    if (ein < emin) then
        ! Linear Extrapolation below tmin
        call ckcvbs(tmin, y, cv)
        t = tmin - (emin-ein) / cv
        ierr = 1
        return
    end if
    if (ein > emax) then
        ! Linear Extrapolation above tmax
        call ckcvbs(tmax, y, cv)
        t = tmax - (emax - ein) / cv
        ierr = 1
        return
    end if
    t1 = t
    if (t1 < tmin .or. t1 > tmax) then
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin)
    end if
    do i=1, maxiter
        call ckubms(t1,y,e1)
        call ckcvbs(t1,y,cv)
        dt = (ein - e1) / cv
        if (dt > 100.d0) then
            dt = 100.d0
        else if (dt < -100.d0) then
            dt = -100.d0
        else if (abs(dt) < tol) then
            exit
        else if (t1+dt == t1) then
            exit
        end if
        t1 = t1 + dt
    end do

    t = t1
    ierr = 0

end subroutine

subroutine egtransetLENIMC(LENIMC)

    implicit none

    integer, intent(out) :: LENIMC

    LENIMC = 38

end subroutine

subroutine egtransetLENRMC(LENRMC)

    implicit none

    integer, intent(out) :: LENRMC

    LENRMC = 1854

end subroutine

subroutine egtransetNO(NO)

    implicit none

    integer, intent(out) :: NO

    NO = 4

end subroutine

subroutine egtransetKK(KK)

    implicit none

    integer, intent(out) :: KK

    KK = 9

end subroutine

subroutine egtransetNLITE(NLITE)

    implicit none

    integer, intent(out) :: NLITE

    NLITE = 2

end subroutine

subroutine egtransetPATM(PATM)

    implicit none

    double precision, intent(out) :: PATM

    PATM = 0.1013250000000000d+07

end subroutine

! the molecular weights in g/mol
subroutine egtransetWT(WT)

    implicit none

    double precision, intent(out) :: WT(9)

    WT(1) = 2.01594000d+00
    WT(2) = 3.19988000d+01
    WT(3) = 1.80153400d+01
    WT(4) = 1.00797000d+00
    WT(5) = 1.59994000d+01
    WT(6) = 1.70073700d+01
    WT(7) = 3.30067700d+01
    WT(8) = 3.40147400d+01
    WT(9) = 2.80134000d+01

end subroutine

! the lennard-jones potential well depth eps/kb in K
subroutine egtransetEPS(EPS)

    implicit none

    double precision, intent(out) :: EPS(9)

    EPS(4) = 1.45000000d+02
    EPS(5) = 8.00000000d+01
    EPS(6) = 8.00000000d+01
    EPS(7) = 1.07400000d+02
    EPS(8) = 1.07400000d+02
    EPS(1) = 3.80000000d+01
    EPS(2) = 1.07400000d+02
    EPS(3) = 5.72400000d+02
    EPS(9) = 9.75300000d+01

end subroutine

! the lennard-jones collision diameter in Angstroms
subroutine egtransetSIG(SIG)

    implicit none

    double precision, intent(out) :: SIG(9)

    SIG(4) = 2.05000000d+00
    SIG(5) = 2.75000000d+00
    SIG(6) = 2.75000000d+00
    SIG(7) = 3.45800000d+00
    SIG(8) = 3.45800000d+00
    SIG(1) = 2.92000000d+00
    SIG(2) = 3.45800000d+00
    SIG(3) = 2.60500000d+00
    SIG(9) = 3.62100000d+00

end subroutine

! the dipole moment in Debye
subroutine egtransetDIP(DIP)

    implicit none

    double precision, intent(out) :: DIP(9)

    DIP(4) = 0.00000000d+00
    DIP(5) = 0.00000000d+00
    DIP(6) = 0.00000000d+00
    DIP(7) = 0.00000000d+00
    DIP(8) = 0.00000000d+00
    DIP(1) = 0.00000000d+00
    DIP(2) = 0.00000000d+00
    DIP(3) = 1.84400000d+00
    DIP(9) = 0.00000000d+00

end subroutine

! the polarizability in cubic Angstroms
subroutine egtransetPOL(POL)

    implicit none

    double precision, intent(out) :: POL(9)

    POL(4) = 0.00000000d+00
    POL(5) = 0.00000000d+00
    POL(6) = 0.00000000d+00
    POL(7) = 0.00000000d+00
    POL(8) = 0.00000000d+00
    POL(1) = 7.90000000d-01
    POL(2) = 1.60000000d+00
    POL(3) = 0.00000000d+00
    POL(9) = 1.76000000d+00

end subroutine

! the rotational relaxation collision number at 298 K
subroutine egtransetZROT(ZROT)

    implicit none

    double precision, intent(out) :: ZROT(9)

    ZROT(4) = 0.00000000d+00
    ZROT(5) = 0.00000000d+00
    ZROT(6) = 0.00000000d+00
    ZROT(7) = 1.00000000d+00
    ZROT(8) = 3.80000000d+00
    ZROT(1) = 2.80000000d+02
    ZROT(2) = 3.80000000d+00
    ZROT(3) = 4.00000000d+00
    ZROT(9) = 4.00000000d+00

end subroutine

! 0: monoatomic, 1: linear, 2: nonlinear
subroutine egtransetNLIN(NLIN)

    implicit none

    integer, intent(out) :: NLIN(9)

    NLIN(4) = 0
    NLIN(5) = 0
    NLIN(6) = 1
    NLIN(7) = 2
    NLIN(8) = 2
    NLIN(1) = 1
    NLIN(2) = 1
    NLIN(3) = 2
    NLIN(9) = 1

end subroutine


! Poly fits for the viscosities, dim NO*KK
subroutine egtransetCOFETA(COFETA)

    implicit none

    double precision, intent(out) :: COFETA(36)

    COFETA(1) = -1.37549435d+01
    COFETA(2) = 9.65530587d-01
    COFETA(3) = -4.45720114d-02
    COFETA(4) = 2.05871810d-03
    COFETA(5) = -1.68118868d+01
    COFETA(6) = 2.52362554d+00
    COFETA(7) = -2.49309128d-01
    COFETA(8) = 1.10211025d-02
    COFETA(9) = -1.17770937d+01
    COFETA(10) = -8.26742721d-01
    COFETA(11) = 3.39009079d-01
    COFETA(12) = -2.00674327d-02
    COFETA(13) = -1.98744496d+01
    COFETA(14) = 3.41660514d+00
    COFETA(15) = -3.63206306d-01
    COFETA(16) = 1.58671021d-02
    COFETA(17) = -1.48001581d+01
    COFETA(18) = 1.79491990d+00
    COFETA(19) = -1.54008440d-01
    COFETA(20) = 6.86719439d-03
    COFETA(21) = -1.47696103d+01
    COFETA(22) = 1.79491990d+00
    COFETA(23) = -1.54008440d-01
    COFETA(24) = 6.86719439d-03
    COFETA(25) = -1.67963797d+01
    COFETA(26) = 2.52362554d+00
    COFETA(27) = -2.49309128d-01
    COFETA(28) = 1.10211025d-02
    COFETA(29) = -1.67813391d+01
    COFETA(30) = 2.52362554d+00
    COFETA(31) = -2.49309128d-01
    COFETA(32) = 1.10211025d-02
    COFETA(33) = -1.62526779d+01
    COFETA(34) = 2.24839597d+00
    COFETA(35) = -2.13428438d-01
    COFETA(36) = 9.46192413d-03

end subroutine


! Poly fits for the conductivities, dim NO*KK
subroutine egtransetCOFLAM(COFLAM)

    implicit none

    double precision, intent(out) :: COFLAM(36)

    COFLAM(1) = 1.11035551d+01
    COFLAM(2) = -1.31884089d+00
    COFLAM(3) = 2.44042735d-01
    COFLAM(4) = -8.99837633d-03
    COFLAM(5) = -2.51299170d+00
    COFLAM(6) = 3.15166792d+00
    COFLAM(7) = -3.10009273d-01
    COFLAM(8) = 1.34523096d-02
    COFLAM(9) = 2.21730076d+01
    COFLAM(10) = -8.46933644d+00
    COFLAM(11) = 1.46153515d+00
    COFLAM(12) = -7.29500936d-02
    COFLAM(13) = -3.24539191d-01
    COFLAM(14) = 3.41660514d+00
    COFLAM(15) = -3.63206306d-01
    COFLAM(16) = 1.58671021d-02
    COFLAM(17) = 1.98513952d+00
    COFLAM(18) = 1.79491990d+00
    COFLAM(19) = -1.54008440d-01
    COFLAM(20) = 6.86719439d-03
    COFLAM(21) = 1.60618734d+01
    COFLAM(22) = -4.10626869d+00
    COFLAM(23) = 6.63571339d-01
    COFLAM(24) = -2.97906324d-02
    COFLAM(25) = 5.56023763d-01
    COFLAM(26) = 1.59073590d+00
    COFLAM(27) = -5.28053839d-02
    COFLAM(28) = 4.07601571d-04
    COFLAM(29) = 1.48801857d+00
    COFLAM(30) = 1.06175905d+00
    COFLAM(31) = 5.72200266d-02
    COFLAM(32) = -6.38393531d-03
    COFLAM(33) = 1.15507419d+01
    COFLAM(34) = -2.91453917d+00
    COFLAM(35) = 5.55045765d-01
    COFLAM(36) = -2.75173485d-02

end subroutine

! Poly fits for the diffusion coefficients, dim NO*KK*KK
subroutine egtransetCOFD(COFD)

    implicit none

    double precision, intent(out) :: COFD(324)

    COFD(1) = -1.02395222d+01
    COFD(2) = 2.15403244d+00
    COFD(3) = -6.97480266d-02
    COFD(4) = 3.23666871d-03
    COFD(5) = -1.15797750d+01
    COFD(6) = 2.43235504d+00
    COFD(7) = -1.02890179d-01
    COFD(8) = 4.52903603d-03
    COFD(9) = -1.68758926d+01
    COFD(10) = 4.49460303d+00
    COFD(11) = -3.64766132d-01
    COFD(12) = 1.56457153d-02
    COFD(13) = -1.11808682d+01
    COFD(14) = 2.66936727d+00
    COFD(15) = -1.34411514d-01
    COFD(16) = 5.92957488d-03
    COFD(17) = -1.06250182d+01
    COFD(18) = 2.15849701d+00
    COFD(19) = -6.53886401d-02
    COFD(20) = 2.81453370d-03
    COFD(21) = -1.06283453d+01
    COFD(22) = 2.15849701d+00
    COFD(23) = -6.53886401d-02
    COFD(24) = 2.81453370d-03
    COFD(25) = -1.15806808d+01
    COFD(26) = 2.43235504d+00
    COFD(27) = -1.02890179d-01
    COFD(28) = 4.52903603d-03
    COFD(29) = -1.15815344d+01
    COFD(30) = 2.43235504d+00
    COFD(31) = -1.02890179d-01
    COFD(32) = 4.52903603d-03
    COFD(33) = -1.13253458d+01
    COFD(34) = 2.31195095d+00
    COFD(35) = -8.63988037d-02
    COFD(36) = 3.77573452d-03
    COFD(37) = -1.15797750d+01
    COFD(38) = 2.43235504d+00
    COFD(39) = -1.02890179d-01
    COFD(40) = 4.52903603d-03
    COFD(41) = -1.53110708d+01
    COFD(42) = 3.37317428d+00
    COFD(43) = -2.24900439d-01
    COFD(44) = 9.81228151d-03
    COFD(45) = -2.10640014d+01
    COFD(46) = 5.50980695d+00
    COFD(47) = -4.78335488d-01
    COFD(48) = 1.98515434d-02
    COFD(49) = -1.43712864d+01
    COFD(50) = 3.70920439d+00
    COFD(51) = -2.67274113d-01
    COFD(52) = 1.15967481d-02
    COFD(53) = -1.40864894d+01
    COFD(54) = 3.07458927d+00
    COFD(55) = -1.86899591d-01
    COFD(56) = 8.19829781d-03
    COFD(57) = -1.41066459d+01
    COFD(58) = 3.07458927d+00
    COFD(59) = -1.86899591d-01
    COFD(60) = 8.19829781d-03
    COFD(61) = -1.53187643d+01
    COFD(62) = 3.37317428d+00
    COFD(63) = -2.24900439d-01
    COFD(64) = 9.81228151d-03
    COFD(65) = -1.53261114d+01
    COFD(66) = 3.37317428d+00
    COFD(67) = -2.24900439d-01
    COFD(68) = 9.81228151d-03
    COFD(69) = -1.50096240d+01
    COFD(70) = 3.25515933d+00
    COFD(71) = -2.09710110d-01
    COFD(72) = 9.15941830d-03
    COFD(73) = -1.68758926d+01
    COFD(74) = 4.49460303d+00
    COFD(75) = -3.64766132d-01
    COFD(76) = 1.56457153d-02
    COFD(77) = -2.10640014d+01
    COFD(78) = 5.50980695d+00
    COFD(79) = -4.78335488d-01
    COFD(80) = 1.98515434d-02
    COFD(81) = -1.31492641d+01
    COFD(82) = 1.48004311d+00
    COFD(83) = 1.60499553d-01
    COFD(84) = -1.19765679d-02
    COFD(85) = -1.93611051d+01
    COFD(86) = 5.51579726d+00
    COFD(87) = -4.76061961d-01
    COFD(88) = 1.96329391d-02
    COFD(89) = -1.91096797d+01
    COFD(90) = 5.02608697d+00
    COFD(91) = -4.26959993d-01
    COFD(92) = 1.80709910d-02
    COFD(93) = -1.91256261d+01
    COFD(94) = 5.02608697d+00
    COFD(95) = -4.26959993d-01
    COFD(96) = 1.80709910d-02
    COFD(97) = -2.04177482d+01
    COFD(98) = 5.31457079d+00
    COFD(99) = -4.58216496d-01
    COFD(100) = 1.91825910d-02
    COFD(101) = -2.04230073d+01
    COFD(102) = 5.31457079d+00
    COFD(103) = -4.58216496d-01
    COFD(104) = 1.91825910d-02
    COFD(105) = -2.08123325d+01
    COFD(106) = 5.42470154d+00
    COFD(107) = -4.69700416d-01
    COFD(108) = 1.95706904d-02
    COFD(109) = -1.11808682d+01
    COFD(110) = 2.66936727d+00
    COFD(111) = -1.34411514d-01
    COFD(112) = 5.92957488d-03
    COFD(113) = -1.43712864d+01
    COFD(114) = 3.70920439d+00
    COFD(115) = -2.67274113d-01
    COFD(116) = 1.15967481d-02
    COFD(117) = -1.93611051d+01
    COFD(118) = 5.51579726d+00
    COFD(119) = -4.76061961d-01
    COFD(120) = 1.96329391d-02
    COFD(121) = -1.43693056d+01
    COFD(122) = 4.03992999d+00
    COFD(123) = -3.08044800d-01
    COFD(124) = 1.32757775d-02
    COFD(125) = -1.31860117d+01
    COFD(126) = 3.38003453d+00
    COFD(127) = -2.25783856d-01
    COFD(128) = 9.85028660d-03
    COFD(129) = -1.31877711d+01
    COFD(130) = 3.38003453d+00
    COFD(131) = -2.25783856d-01
    COFD(132) = 9.85028660d-03
    COFD(133) = -1.43717529d+01
    COFD(134) = 3.70920439d+00
    COFD(135) = -2.67274113d-01
    COFD(136) = 1.15967481d-02
    COFD(137) = -1.43721922d+01
    COFD(138) = 3.70920439d+00
    COFD(139) = -2.67274113d-01
    COFD(140) = 1.15967481d-02
    COFD(141) = -1.40298830d+01
    COFD(142) = 3.55837688d+00
    COFD(143) = -2.47785790d-01
    COFD(144) = 1.07555332d-02
    COFD(145) = -1.06250182d+01
    COFD(146) = 2.15849701d+00
    COFD(147) = -6.53886401d-02
    COFD(148) = 2.81453370d-03
    COFD(149) = -1.40864894d+01
    COFD(150) = 3.07458927d+00
    COFD(151) = -1.86899591d-01
    COFD(152) = 8.19829781d-03
    COFD(153) = -1.91096797d+01
    COFD(154) = 5.02608697d+00
    COFD(155) = -4.26959993d-01
    COFD(156) = 1.80709910d-02
    COFD(157) = -1.31860117d+01
    COFD(158) = 3.38003453d+00
    COFD(159) = -2.25783856d-01
    COFD(160) = 9.85028660d-03
    COFD(161) = -1.29877365d+01
    COFD(162) = 2.80841511d+00
    COFD(163) = -1.52629888d-01
    COFD(164) = 6.72604927d-03
    COFD(165) = -1.30027772d+01
    COFD(166) = 2.80841511d+00
    COFD(167) = -1.52629888d-01
    COFD(168) = 6.72604927d-03
    COFD(169) = -1.40916052d+01
    COFD(170) = 3.07458927d+00
    COFD(171) = -1.86899591d-01
    COFD(172) = 8.19829781d-03
    COFD(173) = -1.40964661d+01
    COFD(174) = 3.07458927d+00
    COFD(175) = -1.86899591d-01
    COFD(176) = 8.19829781d-03
    COFD(177) = -1.38756407d+01
    COFD(178) = 2.98558426d+00
    COFD(179) = -1.75507216d-01
    COFD(180) = 7.71173691d-03
    COFD(181) = -1.06283453d+01
    COFD(182) = 2.15849701d+00
    COFD(183) = -6.53886401d-02
    COFD(184) = 2.81453370d-03
    COFD(185) = -1.41066459d+01
    COFD(186) = 3.07458927d+00
    COFD(187) = -1.86899591d-01
    COFD(188) = 8.19829781d-03
    COFD(189) = -1.91256261d+01
    COFD(190) = 5.02608697d+00
    COFD(191) = -4.26959993d-01
    COFD(192) = 1.80709910d-02
    COFD(193) = -1.31877711d+01
    COFD(194) = 3.38003453d+00
    COFD(195) = -2.25783856d-01
    COFD(196) = 9.85028660d-03
    COFD(197) = -1.30027772d+01
    COFD(198) = 2.80841511d+00
    COFD(199) = -1.52629888d-01
    COFD(200) = 6.72604927d-03
    COFD(201) = -1.30182843d+01
    COFD(202) = 2.80841511d+00
    COFD(203) = -1.52629888d-01
    COFD(204) = 6.72604927d-03
    COFD(205) = -1.41119732d+01
    COFD(206) = 3.07458927d+00
    COFD(207) = -1.86899591d-01
    COFD(208) = 8.19829781d-03
    COFD(209) = -1.41170372d+01
    COFD(210) = 3.07458927d+00
    COFD(211) = -1.86899591d-01
    COFD(212) = 8.19829781d-03
    COFD(213) = -1.38948667d+01
    COFD(214) = 2.98558426d+00
    COFD(215) = -1.75507216d-01
    COFD(216) = 7.71173691d-03
    COFD(217) = -1.15806808d+01
    COFD(218) = 2.43235504d+00
    COFD(219) = -1.02890179d-01
    COFD(220) = 4.52903603d-03
    COFD(221) = -1.53187643d+01
    COFD(222) = 3.37317428d+00
    COFD(223) = -2.24900439d-01
    COFD(224) = 9.81228151d-03
    COFD(225) = -2.04177482d+01
    COFD(226) = 5.31457079d+00
    COFD(227) = -4.58216496d-01
    COFD(228) = 1.91825910d-02
    COFD(229) = -1.43717529d+01
    COFD(230) = 3.70920439d+00
    COFD(231) = -2.67274113d-01
    COFD(232) = 1.15967481d-02
    COFD(233) = -1.40916052d+01
    COFD(234) = 3.07458927d+00
    COFD(235) = -1.86899591d-01
    COFD(236) = 8.19829781d-03
    COFD(237) = -1.41119732d+01
    COFD(238) = 3.07458927d+00
    COFD(239) = -1.86899591d-01
    COFD(240) = 8.19829781d-03
    COFD(241) = -1.53265780d+01
    COFD(242) = 3.37317428d+00
    COFD(243) = -2.24900439d-01
    COFD(244) = 9.81228151d-03
    COFD(245) = -1.53340417d+01
    COFD(246) = 3.37317428d+00
    COFD(247) = -2.24900439d-01
    COFD(248) = 9.81228151d-03
    COFD(249) = -1.50168028d+01
    COFD(250) = 3.25515933d+00
    COFD(251) = -2.09710110d-01
    COFD(252) = 9.15941830d-03
    COFD(253) = -1.15815344d+01
    COFD(254) = 2.43235504d+00
    COFD(255) = -1.02890179d-01
    COFD(256) = 4.52903603d-03
    COFD(257) = -1.53261114d+01
    COFD(258) = 3.37317428d+00
    COFD(259) = -2.24900439d-01
    COFD(260) = 9.81228151d-03
    COFD(261) = -2.04230073d+01
    COFD(262) = 5.31457079d+00
    COFD(263) = -4.58216496d-01
    COFD(264) = 1.91825910d-02
    COFD(265) = -1.43721922d+01
    COFD(266) = 3.70920439d+00
    COFD(267) = -2.67274113d-01
    COFD(268) = 1.15967481d-02
    COFD(269) = -1.40964661d+01
    COFD(270) = 3.07458927d+00
    COFD(271) = -1.86899591d-01
    COFD(272) = 8.19829781d-03
    COFD(273) = -1.41170372d+01
    COFD(274) = 3.07458927d+00
    COFD(275) = -1.86899591d-01
    COFD(276) = 8.19829781d-03
    COFD(277) = -1.53340417d+01
    COFD(278) = 3.37317428d+00
    COFD(279) = -2.24900439d-01
    COFD(280) = 9.81228151d-03
    COFD(281) = -1.53416186d+01
    COFD(282) = 3.37317428d+00
    COFD(283) = -2.24900439d-01
    COFD(284) = 9.81228151d-03
    COFD(285) = -1.50236516d+01
    COFD(286) = 3.25515933d+00
    COFD(287) = -2.09710110d-01
    COFD(288) = 9.15941830d-03
    COFD(289) = -1.13253458d+01
    COFD(290) = 2.31195095d+00
    COFD(291) = -8.63988037d-02
    COFD(292) = 3.77573452d-03
    COFD(293) = -1.50096240d+01
    COFD(294) = 3.25515933d+00
    COFD(295) = -2.09710110d-01
    COFD(296) = 9.15941830d-03
    COFD(297) = -2.08123325d+01
    COFD(298) = 5.42470154d+00
    COFD(299) = -4.69700416d-01
    COFD(300) = 1.95706904d-02
    COFD(301) = -1.40298830d+01
    COFD(302) = 3.55837688d+00
    COFD(303) = -2.47785790d-01
    COFD(304) = 1.07555332d-02
    COFD(305) = -1.38756407d+01
    COFD(306) = 2.98558426d+00
    COFD(307) = -1.75507216d-01
    COFD(308) = 7.71173691d-03
    COFD(309) = -1.38948667d+01
    COFD(310) = 2.98558426d+00
    COFD(311) = -1.75507216d-01
    COFD(312) = 7.71173691d-03
    COFD(313) = -1.50168028d+01
    COFD(314) = 3.25515933d+00
    COFD(315) = -2.09710110d-01
    COFD(316) = 9.15941830d-03
    COFD(317) = -1.50236516d+01
    COFD(318) = 3.25515933d+00
    COFD(319) = -2.09710110d-01
    COFD(320) = 9.15941830d-03
    COFD(321) = -1.47639290d+01
    COFD(322) = 3.15955654d+00
    COFD(323) = -1.97590757d-01
    COFD(324) = 8.64692156d-03

end subroutine

! List of specs with small weight, dim NLITE
subroutine egtransetKTDIF(KTDIF)

    implicit none

    integer, intent(out) :: KTDIF(2)

    KTDIF(1) = 1
    KTDIF(2) = 4

end subroutine


! Poly fits for thermal diff ratios, dim NO*NLITE*KK
subroutine egtransetCOFTD(COFTD)

    implicit none

    double precision, intent(out) :: COFTD(72)

    COFTD(1) = 0.00000000d+00
    COFTD(2) = 0.00000000d+00
    COFTD(3) = 0.00000000d+00
    COFTD(4) = 0.00000000d+00
    COFTD(5) = 4.42739084d-01
    COFTD(6) = 7.11770818d-05
    COFTD(7) = -3.84768062d-08
    COFTD(8) = 6.86323437d-12
    COFTD(9) = 6.02028221d-02
    COFTD(10) = 5.61561867d-04
    COFTD(11) = -2.55372862d-07
    COFTD(12) = 3.63389913d-11
    COFTD(13) = -1.52534742d-01
    COFTD(14) = -5.46404022d-05
    COFTD(15) = 2.93412470d-08
    COFTD(16) = -4.87091914d-12
    COFTD(17) = 4.15583337d-01
    COFTD(18) = 1.09738399d-05
    COFTD(19) = -3.96021963d-09
    COFTD(20) = 1.14414443d-12
    COFTD(21) = 4.21932443d-01
    COFTD(22) = 1.11414935d-05
    COFTD(23) = -4.02072219d-09
    COFTD(24) = 1.16162418d-12
    COFTD(25) = 4.44452569d-01
    COFTD(26) = 7.14525507d-05
    COFTD(27) = -3.86257187d-08
    COFTD(28) = 6.88979640d-12
    COFTD(29) = 4.46070183d-01
    COFTD(30) = 7.17126069d-05
    COFTD(31) = -3.87662996d-08
    COFTD(32) = 6.91487226d-12
    COFTD(33) = 4.45261966d-01
    COFTD(34) = 4.94697174d-05
    COFTD(35) = -2.63023442d-08
    COFTD(36) = 4.90306217d-12
    COFTD(37) = 1.52534742d-01
    COFTD(38) = 5.46404022d-05
    COFTD(39) = -2.93412470d-08
    COFTD(40) = 4.87091914d-12
    COFTD(41) = 2.20482843d-01
    COFTD(42) = 4.80164288d-04
    COFTD(43) = -2.32927944d-07
    COFTD(44) = 3.46470436d-11
    COFTD(45) = -1.41883744d-01
    COFTD(46) = 7.66558810d-04
    COFTD(47) = -3.06550003d-07
    COFTD(48) = 4.02959502d-11
    COFTD(49) = 0.00000000d+00
    COFTD(50) = 0.00000000d+00
    COFTD(51) = 0.00000000d+00
    COFTD(52) = 0.00000000d+00
    COFTD(53) = 2.70010150d-01
    COFTD(54) = 3.61555093d-04
    COFTD(55) = -1.80744752d-07
    COFTD(56) = 2.75321248d-11
    COFTD(57) = 2.72041664d-01
    COFTD(58) = 3.64275376d-04
    COFTD(59) = -1.82104647d-07
    COFTD(60) = 2.77392722d-11
    COFTD(61) = 2.20907853d-01
    COFTD(62) = 4.81089870d-04
    COFTD(63) = -2.33376944d-07
    COFTD(64) = 3.47138305d-11
    COFTD(65) = 2.21308399d-01
    COFTD(66) = 4.81962174d-04
    COFTD(67) = -2.33800100d-07
    COFTD(68) = 3.47767730d-11
    COFTD(69) = 2.40744421d-01
    COFTD(70) = 4.45343451d-04
    COFTD(71) = -2.18173874d-07
    COFTD(72) = 3.26958506d-11

end subroutine

end module fuego_module



