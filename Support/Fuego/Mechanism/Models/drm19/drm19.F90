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
  public :: ckhms
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
double precision, parameter :: imw(21) = (/ &
    1.d0 / 2.015940d0,  & ! H2
    1.d0 / 1.007970d0,  & ! H
    1.d0 / 15.999400d0,  & ! O
    1.d0 / 31.998800d0,  & ! O2
    1.d0 / 17.007370d0,  & ! OH
    1.d0 / 18.015340d0,  & ! H2O
    1.d0 / 33.006770d0,  & ! HO2
    1.d0 / 14.027090d0,  & ! CH2
    1.d0 / 14.027090d0,  & ! CH2(S)
    1.d0 / 15.035060d0,  & ! CH3
    1.d0 / 16.043030d0,  & ! CH4
    1.d0 / 28.010550d0,  & ! CO
    1.d0 / 44.009950d0,  & ! CO2
    1.d0 / 29.018520d0,  & ! HCO
    1.d0 / 30.026490d0,  & ! CH2O
    1.d0 / 31.034460d0,  & ! CH3O
    1.d0 / 28.054180d0,  & ! C2H4
    1.d0 / 29.062150d0,  & ! C2H5
    1.d0 / 30.070120d0,  & ! C2H6
    1.d0 / 28.013400d0,  & ! N2
    1.d0 / 39.948000d0 /) ! AR

type :: nonsquare_matrix_double
   double precision, allocatable :: vector(:)
end type nonsquare_matrix_double

type :: nonsquare_matrix_int
   integer, allocatable :: vector(:)
end type nonsquare_matrix_int

double precision, save :: fwd_A(84), fwd_beta(84), fwd_Ea(84)
double precision, save :: low_A(84), low_beta(84), low_Ea(84)
double precision, save :: rev_A(84), rev_beta(84), rev_Ea(84)
double precision, save :: troe_a(84),troe_Ts(84), troe_Tss(84), troe_Tsss(84)
double precision, save :: sri_a(84), sri_b(84), sri_c(84), sri_d(84), sri_e(84)
double precision, save :: activation_units(84), prefactor_units(84), phase_units(84)
integer, save :: is_PD(84), troe_len(84), sri_len(84), nTB(84)
type(nonsquare_matrix_double) :: TB(84)
type(nonsquare_matrix_int) :: TBid(84)

double precision, save :: fwd_A_DEF(84), fwd_beta_DEF(84), fwd_Ea_DEF(84)
double precision, save :: low_A_DEF(84), low_beta_DEF(84), low_Ea_DEF(84)
double precision, save :: rev_A_DEF(84), rev_beta_DEF(84), rev_Ea_DEF(84)
double precision, save :: troe_a_DEF(84),troe_Ts_DEF(84), troe_Tss_DEF(84), troe_Tsss_DEF(84)
double precision, save :: sri_a_DEF(84), sri_b_DEF(84), sri_c_DEF(84), sri_d_DEF(84), sri_e_DEF(84)
double precision, save :: activation_units_DEF(84), prefactor_units_DEF(84), phase_units_DEF(84)
integer, save :: is_PD_DEF(84), troe_len_DEF(84), sri_len_DEF(84), nTB_DEF(84)
type(nonsquare_matrix_double) :: TB_DEF(84)
type(nonsquare_matrix_int) :: TBid_DEF(84)

! productionRate() static variables
double precision, save :: T_save = -1
double precision, save :: k_f_save(84)
double precision, save :: Kc_save(84)

contains

subroutine SetAllDefaults()

    implicit none

    integer :: i, j

    do i=1, 84
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

  do i=1, 84
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

    ! (0):  O + H + M <=> OH + M
    fwd_A(9)     = 5d+17
    fwd_beta(9)  = -1d0
    fwd_Ea(9)    = 0d0
    prefactor_units(9)  = 1.0000000000000002d-12
    activation_units(9) = 0.50321666580471969d0
    phase_units(9)      = 1d-12
    is_PD(9) = 0
    nTB(9) = 7
    if (.not. allocated(TB(9) % vector)) allocate(TB(9) % vector(7))
    if (.not. allocated(TBid(9) % vector)) allocate(TBid(9) % vector(7))
    TBid(9) % vector(1) = 0d0
    TB(9) % vector(1) = 2d0 ! H2
    TBid(9) % vector(2) = 5d0
    TB(9) % vector(2) = 6d0 ! H2O
    TBid(9) % vector(3) = 10d0
    TB(9) % vector(3) = 2d0 ! CH4
    TBid(9) % vector(4) = 11d0
    TB(9) % vector(4) = 1.5d0 ! CO
    TBid(9) % vector(5) = 12d0
    TB(9) % vector(5) = 2d0 ! CO2
    TBid(9) % vector(6) = 18d0
    TB(9) % vector(6) = 3d0 ! C2H6
    TBid(9) % vector(7) = 20d0
    TB(9) % vector(7) = 0.69999999999999996d0 ! AR

    ! (1):  O + H2 <=> H + OH
    fwd_A(15)     = 50000d0
    fwd_beta(15)  = 2.6699999999999999d0
    fwd_Ea(15)    = 6290d0
    prefactor_units(15)  = 1.0000000000000002d-06
    activation_units(15) = 0.50321666580471969d0
    phase_units(15)      = 1d-12
    is_PD(15) = 0
    nTB(15) = 0

    ! (2):  O + HO2 <=> OH + O2
    fwd_A(16)     = 20000000000000d0
    fwd_beta(16)  = 0d0
    fwd_Ea(16)    = 0d0
    prefactor_units(16)  = 1.0000000000000002d-06
    activation_units(16) = 0.50321666580471969d0
    phase_units(16)      = 1d-12
    is_PD(16) = 0
    nTB(16) = 0

    ! (3):  O + CH2 <=> H + HCO
    fwd_A(17)     = 80000000000000d0
    fwd_beta(17)  = 0d0
    fwd_Ea(17)    = 0d0
    prefactor_units(17)  = 1.0000000000000002d-06
    activation_units(17) = 0.50321666580471969d0
    phase_units(17)      = 1d-12
    is_PD(17) = 0
    nTB(17) = 0

    ! (4):  O + CH2(S) <=> H + HCO
    fwd_A(18)     = 15000000000000d0
    fwd_beta(18)  = 0d0
    fwd_Ea(18)    = 0d0
    prefactor_units(18)  = 1.0000000000000002d-06
    activation_units(18) = 0.50321666580471969d0
    phase_units(18)      = 1d-12
    is_PD(18) = 0
    nTB(18) = 0

    ! (5):  O + CH3 <=> H + CH2O
    fwd_A(19)     = 84300000000000d0
    fwd_beta(19)  = 0d0
    fwd_Ea(19)    = 0d0
    prefactor_units(19)  = 1.0000000000000002d-06
    activation_units(19) = 0.50321666580471969d0
    phase_units(19)      = 1d-12
    is_PD(19) = 0
    nTB(19) = 0

    ! (6):  O + CH4 <=> OH + CH3
    fwd_A(20)     = 1020000000d0
    fwd_beta(20)  = 1.5d0
    fwd_Ea(20)    = 8600d0
    prefactor_units(20)  = 1.0000000000000002d-06
    activation_units(20) = 0.50321666580471969d0
    phase_units(20)      = 1d-12
    is_PD(20) = 0
    nTB(20) = 0

    ! (7):  O + CO + M <=> CO2 + M
    fwd_A(10)     = 602000000000000d0
    fwd_beta(10)  = 0d0
    fwd_Ea(10)    = 3000d0
    prefactor_units(10)  = 1.0000000000000002d-12
    activation_units(10) = 0.50321666580471969d0
    phase_units(10)      = 1d-12
    is_PD(10) = 0
    nTB(10) = 8
    if (.not. allocated(TB(10) % vector)) allocate(TB(10) % vector(8))
    if (.not. allocated(TBid(10) % vector)) allocate(TBid(10) % vector(8))
    TBid(10) % vector(1) = 0d0
    TB(10) % vector(1) = 2d0 ! H2
    TBid(10) % vector(2) = 3d0
    TB(10) % vector(2) = 6d0 ! O2
    TBid(10) % vector(3) = 5d0
    TB(10) % vector(3) = 6d0 ! H2O
    TBid(10) % vector(4) = 10d0
    TB(10) % vector(4) = 2d0 ! CH4
    TBid(10) % vector(5) = 11d0
    TB(10) % vector(5) = 1.5d0 ! CO
    TBid(10) % vector(6) = 12d0
    TB(10) % vector(6) = 3.5d0 ! CO2
    TBid(10) % vector(7) = 18d0
    TB(10) % vector(7) = 3d0 ! C2H6
    TBid(10) % vector(8) = 20d0
    TB(10) % vector(8) = 0.5d0 ! AR

    ! (8):  O + HCO <=> OH + CO
    fwd_A(21)     = 30000000000000d0
    fwd_beta(21)  = 0d0
    fwd_Ea(21)    = 0d0
    prefactor_units(21)  = 1.0000000000000002d-06
    activation_units(21) = 0.50321666580471969d0
    phase_units(21)      = 1d-12
    is_PD(21) = 0
    nTB(21) = 0

    ! (9):  O + HCO <=> H + CO2
    fwd_A(22)     = 30000000000000d0
    fwd_beta(22)  = 0d0
    fwd_Ea(22)    = 0d0
    prefactor_units(22)  = 1.0000000000000002d-06
    activation_units(22) = 0.50321666580471969d0
    phase_units(22)      = 1d-12
    is_PD(22) = 0
    nTB(22) = 0

    ! (10):  O + CH2O <=> OH + HCO
    fwd_A(23)     = 39000000000000d0
    fwd_beta(23)  = 0d0
    fwd_Ea(23)    = 3540d0
    prefactor_units(23)  = 1.0000000000000002d-06
    activation_units(23) = 0.50321666580471969d0
    phase_units(23)      = 1d-12
    is_PD(23) = 0
    nTB(23) = 0

    ! (11):  O + C2H4 <=> CH3 + HCO
    fwd_A(24)     = 19200000d0
    fwd_beta(24)  = 1.8300000000000001d0
    fwd_Ea(24)    = 220d0
    prefactor_units(24)  = 1.0000000000000002d-06
    activation_units(24) = 0.50321666580471969d0
    phase_units(24)      = 1d-12
    is_PD(24) = 0
    nTB(24) = 0

    ! (12):  O + C2H5 <=> CH3 + CH2O
    fwd_A(25)     = 132000000000000d0
    fwd_beta(25)  = 0d0
    fwd_Ea(25)    = 0d0
    prefactor_units(25)  = 1.0000000000000002d-06
    activation_units(25) = 0.50321666580471969d0
    phase_units(25)      = 1d-12
    is_PD(25) = 0
    nTB(25) = 0

    ! (13):  O + C2H6 <=> OH + C2H5
    fwd_A(26)     = 89800000d0
    fwd_beta(26)  = 1.9199999999999999d0
    fwd_Ea(26)    = 5690d0
    prefactor_units(26)  = 1.0000000000000002d-06
    activation_units(26) = 0.50321666580471969d0
    phase_units(26)      = 1d-12
    is_PD(26) = 0
    nTB(26) = 0

    ! (14):  O2 + CO <=> O + CO2
    fwd_A(27)     = 2500000000000d0
    fwd_beta(27)  = 0d0
    fwd_Ea(27)    = 47800d0
    prefactor_units(27)  = 1.0000000000000002d-06
    activation_units(27) = 0.50321666580471969d0
    phase_units(27)      = 1d-12
    is_PD(27) = 0
    nTB(27) = 0

    ! (15):  O2 + CH2O <=> HO2 + HCO
    fwd_A(28)     = 100000000000000d0
    fwd_beta(28)  = 0d0
    fwd_Ea(28)    = 40000d0
    prefactor_units(28)  = 1.0000000000000002d-06
    activation_units(28) = 0.50321666580471969d0
    phase_units(28)      = 1d-12
    is_PD(28) = 0
    nTB(28) = 0

    ! (16):  H + O2 + M <=> HO2 + M
    fwd_A(11)     = 2.8d+18
    fwd_beta(11)  = -0.85999999999999999d0
    fwd_Ea(11)    = 0d0
    prefactor_units(11)  = 1.0000000000000002d-12
    activation_units(11) = 0.50321666580471969d0
    phase_units(11)      = 1d-12
    is_PD(11) = 0
    nTB(11) = 7
    if (.not. allocated(TB(11) % vector)) allocate(TB(11) % vector(7))
    if (.not. allocated(TBid(11) % vector)) allocate(TBid(11) % vector(7))
    TBid(11) % vector(1) = 3d0
    TB(11) % vector(1) = 0d0 ! O2
    TBid(11) % vector(2) = 5d0
    TB(11) % vector(2) = 0d0 ! H2O
    TBid(11) % vector(3) = 11d0
    TB(11) % vector(3) = 0.75d0 ! CO
    TBid(11) % vector(4) = 12d0
    TB(11) % vector(4) = 1.5d0 ! CO2
    TBid(11) % vector(5) = 18d0
    TB(11) % vector(5) = 1.5d0 ! C2H6
    TBid(11) % vector(6) = 19d0
    TB(11) % vector(6) = 0d0 ! N2
    TBid(11) % vector(7) = 20d0
    TB(11) % vector(7) = 0d0 ! AR

    ! (17):  H + 2 O2 <=> HO2 + O2
    fwd_A(29)     = 3d+20
    fwd_beta(29)  = -1.72d0
    fwd_Ea(29)    = 0d0
    prefactor_units(29)  = 1.0000000000000002d-12
    activation_units(29) = 0.50321666580471969d0
    phase_units(29)      = 1d-18
    is_PD(29) = 0
    nTB(29) = 0

    ! (18):  H + O2 + H2O <=> HO2 + H2O
    fwd_A(30)     = 9.38d+18
    fwd_beta(30)  = -0.76000000000000001d0
    fwd_Ea(30)    = 0d0
    prefactor_units(30)  = 1.0000000000000002d-12
    activation_units(30) = 0.50321666580471969d0
    phase_units(30)      = 1d-18
    is_PD(30) = 0
    nTB(30) = 0

    ! (19):  H + O2 + N2 <=> HO2 + N2
    fwd_A(31)     = 3.75d+20
    fwd_beta(31)  = -1.72d0
    fwd_Ea(31)    = 0d0
    prefactor_units(31)  = 1.0000000000000002d-12
    activation_units(31) = 0.50321666580471969d0
    phase_units(31)      = 1d-18
    is_PD(31) = 0
    nTB(31) = 0

    ! (20):  H + O2 + AR <=> HO2 + AR
    fwd_A(32)     = 7d+17
    fwd_beta(32)  = -0.80000000000000004d0
    fwd_Ea(32)    = 0d0
    prefactor_units(32)  = 1.0000000000000002d-12
    activation_units(32) = 0.50321666580471969d0
    phase_units(32)      = 1d-18
    is_PD(32) = 0
    nTB(32) = 0

    ! (21):  H + O2 <=> O + OH
    fwd_A(33)     = 83000000000000d0
    fwd_beta(33)  = 0d0
    fwd_Ea(33)    = 14413d0
    prefactor_units(33)  = 1.0000000000000002d-06
    activation_units(33) = 0.50321666580471969d0
    phase_units(33)      = 1d-12
    is_PD(33) = 0
    nTB(33) = 0

    ! (22):  2 H + M <=> H2 + M
    fwd_A(12)     = 1d+18
    fwd_beta(12)  = -1d0
    fwd_Ea(12)    = 0d0
    prefactor_units(12)  = 1.0000000000000002d-12
    activation_units(12) = 0.50321666580471969d0
    phase_units(12)      = 1d-12
    is_PD(12) = 0
    nTB(12) = 6
    if (.not. allocated(TB(12) % vector)) allocate(TB(12) % vector(6))
    if (.not. allocated(TBid(12) % vector)) allocate(TBid(12) % vector(6))
    TBid(12) % vector(1) = 0d0
    TB(12) % vector(1) = 0d0 ! H2
    TBid(12) % vector(2) = 5d0
    TB(12) % vector(2) = 0d0 ! H2O
    TBid(12) % vector(3) = 10d0
    TB(12) % vector(3) = 2d0 ! CH4
    TBid(12) % vector(4) = 12d0
    TB(12) % vector(4) = 0d0 ! CO2
    TBid(12) % vector(5) = 18d0
    TB(12) % vector(5) = 3d0 ! C2H6
    TBid(12) % vector(6) = 20d0
    TB(12) % vector(6) = 0.63d0 ! AR

    ! (23):  2 H + H2 <=> 2 H2
    fwd_A(34)     = 90000000000000000d0
    fwd_beta(34)  = -0.59999999999999998d0
    fwd_Ea(34)    = 0d0
    prefactor_units(34)  = 1.0000000000000002d-12
    activation_units(34) = 0.50321666580471969d0
    phase_units(34)      = 1d-18
    is_PD(34) = 0
    nTB(34) = 0

    ! (24):  2 H + H2O <=> H2 + H2O
    fwd_A(35)     = 6d+19
    fwd_beta(35)  = -1.25d0
    fwd_Ea(35)    = 0d0
    prefactor_units(35)  = 1.0000000000000002d-12
    activation_units(35) = 0.50321666580471969d0
    phase_units(35)      = 1d-18
    is_PD(35) = 0
    nTB(35) = 0

    ! (25):  2 H + CO2 <=> H2 + CO2
    fwd_A(36)     = 5.5d+20
    fwd_beta(36)  = -2d0
    fwd_Ea(36)    = 0d0
    prefactor_units(36)  = 1.0000000000000002d-12
    activation_units(36) = 0.50321666580471969d0
    phase_units(36)      = 1d-18
    is_PD(36) = 0
    nTB(36) = 0

    ! (26):  H + OH + M <=> H2O + M
    fwd_A(13)     = 2.2d+22
    fwd_beta(13)  = -2d0
    fwd_Ea(13)    = 0d0
    prefactor_units(13)  = 1.0000000000000002d-12
    activation_units(13) = 0.50321666580471969d0
    phase_units(13)      = 1d-12
    is_PD(13) = 0
    nTB(13) = 5
    if (.not. allocated(TB(13) % vector)) allocate(TB(13) % vector(5))
    if (.not. allocated(TBid(13) % vector)) allocate(TBid(13) % vector(5))
    TBid(13) % vector(1) = 0d0
    TB(13) % vector(1) = 0.72999999999999998d0 ! H2
    TBid(13) % vector(2) = 5d0
    TB(13) % vector(2) = 3.6499999999999999d0 ! H2O
    TBid(13) % vector(3) = 10d0
    TB(13) % vector(3) = 2d0 ! CH4
    TBid(13) % vector(4) = 18d0
    TB(13) % vector(4) = 3d0 ! C2H6
    TBid(13) % vector(5) = 20d0
    TB(13) % vector(5) = 0.38d0 ! AR

    ! (27):  H + HO2 <=> O2 + H2
    fwd_A(37)     = 28000000000000d0
    fwd_beta(37)  = 0d0
    fwd_Ea(37)    = 1068d0
    prefactor_units(37)  = 1.0000000000000002d-06
    activation_units(37) = 0.50321666580471969d0
    phase_units(37)      = 1d-12
    is_PD(37) = 0
    nTB(37) = 0

    ! (28):  H + HO2 <=> 2 OH
    fwd_A(38)     = 134000000000000d0
    fwd_beta(38)  = 0d0
    fwd_Ea(38)    = 635d0
    prefactor_units(38)  = 1.0000000000000002d-06
    activation_units(38) = 0.50321666580471969d0
    phase_units(38)      = 1d-12
    is_PD(38) = 0
    nTB(38) = 0

    ! (29):  H + CH2 (+M) <=> CH3 (+M)
    fwd_A(1)     = 25000000000000000d0
    fwd_beta(1)  = -0.80000000000000004d0
    fwd_Ea(1)    = 0d0
    low_A(1)     = 3.2000000000000002d+27
    low_beta(1)  = -3.1400000000000001d0
    low_Ea(1)    = 1230d0
    troe_a(1)    = 0.68000000000000005d0
    troe_Tsss(1) = 78d0
    troe_Ts(1)   = 1995d0
    troe_Tss(1)  = 5590d0
    troe_len(1)  = 4
    prefactor_units(1)  = 1.0000000000000002d-06
    activation_units(1) = 0.50321666580471969d0
    phase_units(1)      = 1d-12
    is_PD(1) = 1
    nTB(1) = 7
    if (.not. allocated(TB(1) % vector)) allocate(TB(1) % vector(7))
    if (.not. allocated(TBid(1) % vector)) allocate(TBid(1) % vector(7))
    TBid(1) % vector(1) = 0d0
    TB(1) % vector(1) = 2d0 ! H2
    TBid(1) % vector(2) = 5d0
    TB(1) % vector(2) = 6d0 ! H2O
    TBid(1) % vector(3) = 10d0
    TB(1) % vector(3) = 2d0 ! CH4
    TBid(1) % vector(4) = 11d0
    TB(1) % vector(4) = 1.5d0 ! CO
    TBid(1) % vector(5) = 12d0
    TB(1) % vector(5) = 2d0 ! CO2
    TBid(1) % vector(6) = 18d0
    TB(1) % vector(6) = 3d0 ! C2H6
    TBid(1) % vector(7) = 20d0
    TB(1) % vector(7) = 0.69999999999999996d0 ! AR

    ! (30):  H + CH3 (+M) <=> CH4 (+M)
    fwd_A(2)     = 12700000000000000d0
    fwd_beta(2)  = -0.63d0
    fwd_Ea(2)    = 383d0
    low_A(2)     = 2.4769999999999999d+33
    low_beta(2)  = -4.7599999999999998d0
    low_Ea(2)    = 2440d0
    troe_a(2)    = 0.78300000000000003d0
    troe_Tsss(2) = 74d0
    troe_Ts(2)   = 2941d0
    troe_Tss(2)  = 6964d0
    troe_len(2)  = 4
    prefactor_units(2)  = 1.0000000000000002d-06
    activation_units(2) = 0.50321666580471969d0
    phase_units(2)      = 1d-12
    is_PD(2) = 1
    nTB(2) = 7
    if (.not. allocated(TB(2) % vector)) allocate(TB(2) % vector(7))
    if (.not. allocated(TBid(2) % vector)) allocate(TBid(2) % vector(7))
    TBid(2) % vector(1) = 0d0
    TB(2) % vector(1) = 2d0 ! H2
    TBid(2) % vector(2) = 5d0
    TB(2) % vector(2) = 6d0 ! H2O
    TBid(2) % vector(3) = 10d0
    TB(2) % vector(3) = 2d0 ! CH4
    TBid(2) % vector(4) = 11d0
    TB(2) % vector(4) = 1.5d0 ! CO
    TBid(2) % vector(5) = 12d0
    TB(2) % vector(5) = 2d0 ! CO2
    TBid(2) % vector(6) = 18d0
    TB(2) % vector(6) = 3d0 ! C2H6
    TBid(2) % vector(7) = 20d0
    TB(2) % vector(7) = 0.69999999999999996d0 ! AR

    ! (31):  H + CH4 <=> CH3 + H2
    fwd_A(39)     = 660000000d0
    fwd_beta(39)  = 1.6200000000000001d0
    fwd_Ea(39)    = 10840d0
    prefactor_units(39)  = 1.0000000000000002d-06
    activation_units(39) = 0.50321666580471969d0
    phase_units(39)      = 1d-12
    is_PD(39) = 0
    nTB(39) = 0

    ! (32):  H + HCO (+M) <=> CH2O (+M)
    fwd_A(3)     = 1090000000000d0
    fwd_beta(3)  = 0.47999999999999998d0
    fwd_Ea(3)    = -260d0
    low_A(3)     = 1.35d+24
    low_beta(3)  = -2.5699999999999998d0
    low_Ea(3)    = 1425d0
    troe_a(3)    = 0.78239999999999998d0
    troe_Tsss(3) = 271d0
    troe_Ts(3)   = 2755d0
    troe_Tss(3)  = 6570d0
    troe_len(3)  = 4
    prefactor_units(3)  = 1.0000000000000002d-06
    activation_units(3) = 0.50321666580471969d0
    phase_units(3)      = 1d-12
    is_PD(3) = 1
    nTB(3) = 7
    if (.not. allocated(TB(3) % vector)) allocate(TB(3) % vector(7))
    if (.not. allocated(TBid(3) % vector)) allocate(TBid(3) % vector(7))
    TBid(3) % vector(1) = 0d0
    TB(3) % vector(1) = 2d0 ! H2
    TBid(3) % vector(2) = 5d0
    TB(3) % vector(2) = 6d0 ! H2O
    TBid(3) % vector(3) = 10d0
    TB(3) % vector(3) = 2d0 ! CH4
    TBid(3) % vector(4) = 11d0
    TB(3) % vector(4) = 1.5d0 ! CO
    TBid(3) % vector(5) = 12d0
    TB(3) % vector(5) = 2d0 ! CO2
    TBid(3) % vector(6) = 18d0
    TB(3) % vector(6) = 3d0 ! C2H6
    TBid(3) % vector(7) = 20d0
    TB(3) % vector(7) = 0.69999999999999996d0 ! AR

    ! (33):  H + HCO <=> H2 + CO
    fwd_A(40)     = 73400000000000d0
    fwd_beta(40)  = 0d0
    fwd_Ea(40)    = 0d0
    prefactor_units(40)  = 1.0000000000000002d-06
    activation_units(40) = 0.50321666580471969d0
    phase_units(40)      = 1d-12
    is_PD(40) = 0
    nTB(40) = 0

    ! (34):  H + CH2O (+M) <=> CH3O (+M)
    fwd_A(4)     = 540000000000d0
    fwd_beta(4)  = 0.45400000000000001d0
    fwd_Ea(4)    = 2600d0
    low_A(4)     = 2.2d+30
    low_beta(4)  = -4.7999999999999998d0
    low_Ea(4)    = 5560d0
    troe_a(4)    = 0.75800000000000001d0
    troe_Tsss(4) = 94d0
    troe_Ts(4)   = 1555d0
    troe_Tss(4)  = 4200d0
    troe_len(4)  = 4
    prefactor_units(4)  = 1.0000000000000002d-06
    activation_units(4) = 0.50321666580471969d0
    phase_units(4)      = 1d-12
    is_PD(4) = 1
    nTB(4) = 6
    if (.not. allocated(TB(4) % vector)) allocate(TB(4) % vector(6))
    if (.not. allocated(TBid(4) % vector)) allocate(TBid(4) % vector(6))
    TBid(4) % vector(1) = 0d0
    TB(4) % vector(1) = 2d0 ! H2
    TBid(4) % vector(2) = 5d0
    TB(4) % vector(2) = 6d0 ! H2O
    TBid(4) % vector(3) = 10d0
    TB(4) % vector(3) = 2d0 ! CH4
    TBid(4) % vector(4) = 11d0
    TB(4) % vector(4) = 1.5d0 ! CO
    TBid(4) % vector(5) = 12d0
    TB(4) % vector(5) = 2d0 ! CO2
    TBid(4) % vector(6) = 18d0
    TB(4) % vector(6) = 3d0 ! C2H6

    ! (35):  H + CH2O <=> HCO + H2
    fwd_A(41)     = 23000000000d0
    fwd_beta(41)  = 1.05d0
    fwd_Ea(41)    = 3275d0
    prefactor_units(41)  = 1.0000000000000002d-06
    activation_units(41) = 0.50321666580471969d0
    phase_units(41)      = 1d-12
    is_PD(41) = 0
    nTB(41) = 0

    ! (36):  H + CH3O <=> OH + CH3
    fwd_A(42)     = 32000000000000d0
    fwd_beta(42)  = 0d0
    fwd_Ea(42)    = 0d0
    prefactor_units(42)  = 1.0000000000000002d-06
    activation_units(42) = 0.50321666580471969d0
    phase_units(42)      = 1d-12
    is_PD(42) = 0
    nTB(42) = 0

    ! (37):  H + C2H4 (+M) <=> C2H5 (+M)
    fwd_A(5)     = 1080000000000d0
    fwd_beta(5)  = 0.45400000000000001d0
    fwd_Ea(5)    = 1820d0
    low_A(5)     = 1.1999999999999999d+42
    low_beta(5)  = -7.6200000000000001d0
    low_Ea(5)    = 6970d0
    troe_a(5)    = 0.97529999999999994d0
    troe_Tsss(5) = 210d0
    troe_Ts(5)   = 984d0
    troe_Tss(5)  = 4374d0
    troe_len(5)  = 4
    prefactor_units(5)  = 1.0000000000000002d-06
    activation_units(5) = 0.50321666580471969d0
    phase_units(5)      = 1d-12
    is_PD(5) = 1
    nTB(5) = 7
    if (.not. allocated(TB(5) % vector)) allocate(TB(5) % vector(7))
    if (.not. allocated(TBid(5) % vector)) allocate(TBid(5) % vector(7))
    TBid(5) % vector(1) = 0d0
    TB(5) % vector(1) = 2d0 ! H2
    TBid(5) % vector(2) = 5d0
    TB(5) % vector(2) = 6d0 ! H2O
    TBid(5) % vector(3) = 10d0
    TB(5) % vector(3) = 2d0 ! CH4
    TBid(5) % vector(4) = 11d0
    TB(5) % vector(4) = 1.5d0 ! CO
    TBid(5) % vector(5) = 12d0
    TB(5) % vector(5) = 2d0 ! CO2
    TBid(5) % vector(6) = 18d0
    TB(5) % vector(6) = 3d0 ! C2H6
    TBid(5) % vector(7) = 20d0
    TB(5) % vector(7) = 0.69999999999999996d0 ! AR

    ! (38):  H + C2H5 (+M) <=> C2H6 (+M)
    fwd_A(6)     = 5.21d+17
    fwd_beta(6)  = -0.98999999999999999d0
    fwd_Ea(6)    = 1580d0
    low_A(6)     = 1.9900000000000001d+41
    low_beta(6)  = -7.0800000000000001d0
    low_Ea(6)    = 6685d0
    troe_a(6)    = 0.84219999999999995d0
    troe_Tsss(6) = 125d0
    troe_Ts(6)   = 2219d0
    troe_Tss(6)  = 6882d0
    troe_len(6)  = 4
    prefactor_units(6)  = 1.0000000000000002d-06
    activation_units(6) = 0.50321666580471969d0
    phase_units(6)      = 1d-12
    is_PD(6) = 1
    nTB(6) = 7
    if (.not. allocated(TB(6) % vector)) allocate(TB(6) % vector(7))
    if (.not. allocated(TBid(6) % vector)) allocate(TBid(6) % vector(7))
    TBid(6) % vector(1) = 0d0
    TB(6) % vector(1) = 2d0 ! H2
    TBid(6) % vector(2) = 5d0
    TB(6) % vector(2) = 6d0 ! H2O
    TBid(6) % vector(3) = 10d0
    TB(6) % vector(3) = 2d0 ! CH4
    TBid(6) % vector(4) = 11d0
    TB(6) % vector(4) = 1.5d0 ! CO
    TBid(6) % vector(5) = 12d0
    TB(6) % vector(5) = 2d0 ! CO2
    TBid(6) % vector(6) = 18d0
    TB(6) % vector(6) = 3d0 ! C2H6
    TBid(6) % vector(7) = 20d0
    TB(6) % vector(7) = 0.69999999999999996d0 ! AR

    ! (39):  H + C2H6 <=> C2H5 + H2
    fwd_A(43)     = 115000000d0
    fwd_beta(43)  = 1.8999999999999999d0
    fwd_Ea(43)    = 7530d0
    prefactor_units(43)  = 1.0000000000000002d-06
    activation_units(43) = 0.50321666580471969d0
    phase_units(43)      = 1d-12
    is_PD(43) = 0
    nTB(43) = 0

    ! (40):  H2 + CO (+M) <=> CH2O (+M)
    fwd_A(7)     = 43000000d0
    fwd_beta(7)  = 1.5d0
    fwd_Ea(7)    = 79600d0
    low_A(7)     = 5.0699999999999998d+27
    low_beta(7)  = -3.4199999999999999d0
    low_Ea(7)    = 84350d0
    troe_a(7)    = 0.93200000000000005d0
    troe_Tsss(7) = 197d0
    troe_Ts(7)   = 1540d0
    troe_Tss(7)  = 10300d0
    troe_len(7)  = 4
    prefactor_units(7)  = 1.0000000000000002d-06
    activation_units(7) = 0.50321666580471969d0
    phase_units(7)      = 1d-12
    is_PD(7) = 1
    nTB(7) = 7
    if (.not. allocated(TB(7) % vector)) allocate(TB(7) % vector(7))
    if (.not. allocated(TBid(7) % vector)) allocate(TBid(7) % vector(7))
    TBid(7) % vector(1) = 0d0
    TB(7) % vector(1) = 2d0 ! H2
    TBid(7) % vector(2) = 5d0
    TB(7) % vector(2) = 6d0 ! H2O
    TBid(7) % vector(3) = 10d0
    TB(7) % vector(3) = 2d0 ! CH4
    TBid(7) % vector(4) = 11d0
    TB(7) % vector(4) = 1.5d0 ! CO
    TBid(7) % vector(5) = 12d0
    TB(7) % vector(5) = 2d0 ! CO2
    TBid(7) % vector(6) = 18d0
    TB(7) % vector(6) = 3d0 ! C2H6
    TBid(7) % vector(7) = 20d0
    TB(7) % vector(7) = 0.69999999999999996d0 ! AR

    ! (41):  OH + H2 <=> H + H2O
    fwd_A(44)     = 216000000d0
    fwd_beta(44)  = 1.51d0
    fwd_Ea(44)    = 3430d0
    prefactor_units(44)  = 1.0000000000000002d-06
    activation_units(44) = 0.50321666580471969d0
    phase_units(44)      = 1d-12
    is_PD(44) = 0
    nTB(44) = 0

    ! (42):  2 OH <=> O + H2O
    fwd_A(45)     = 35700d0
    fwd_beta(45)  = 2.3999999999999999d0
    fwd_Ea(45)    = -2110d0
    prefactor_units(45)  = 1.0000000000000002d-06
    activation_units(45) = 0.50321666580471969d0
    phase_units(45)      = 1d-12
    is_PD(45) = 0
    nTB(45) = 0

    ! (43):  OH + HO2 <=> O2 + H2O
    fwd_A(46)     = 29000000000000d0
    fwd_beta(46)  = 0d0
    fwd_Ea(46)    = -500d0
    prefactor_units(46)  = 1.0000000000000002d-06
    activation_units(46) = 0.50321666580471969d0
    phase_units(46)      = 1d-12
    is_PD(46) = 0
    nTB(46) = 0

    ! (44):  OH + CH2 <=> H + CH2O
    fwd_A(47)     = 20000000000000d0
    fwd_beta(47)  = 0d0
    fwd_Ea(47)    = 0d0
    prefactor_units(47)  = 1.0000000000000002d-06
    activation_units(47) = 0.50321666580471969d0
    phase_units(47)      = 1d-12
    is_PD(47) = 0
    nTB(47) = 0

    ! (45):  OH + CH2(S) <=> H + CH2O
    fwd_A(48)     = 30000000000000d0
    fwd_beta(48)  = 0d0
    fwd_Ea(48)    = 0d0
    prefactor_units(48)  = 1.0000000000000002d-06
    activation_units(48) = 0.50321666580471969d0
    phase_units(48)      = 1d-12
    is_PD(48) = 0
    nTB(48) = 0

    ! (46):  OH + CH3 <=> CH2 + H2O
    fwd_A(49)     = 56000000d0
    fwd_beta(49)  = 1.6000000000000001d0
    fwd_Ea(49)    = 5420d0
    prefactor_units(49)  = 1.0000000000000002d-06
    activation_units(49) = 0.50321666580471969d0
    phase_units(49)      = 1d-12
    is_PD(49) = 0
    nTB(49) = 0

    ! (47):  OH + CH3 <=> CH2(S) + H2O
    fwd_A(50)     = 25010000000000d0
    fwd_beta(50)  = 0d0
    fwd_Ea(50)    = 0d0
    prefactor_units(50)  = 1.0000000000000002d-06
    activation_units(50) = 0.50321666580471969d0
    phase_units(50)      = 1d-12
    is_PD(50) = 0
    nTB(50) = 0

    ! (48):  OH + CH4 <=> CH3 + H2O
    fwd_A(51)     = 100000000d0
    fwd_beta(51)  = 1.6000000000000001d0
    fwd_Ea(51)    = 3120d0
    prefactor_units(51)  = 1.0000000000000002d-06
    activation_units(51) = 0.50321666580471969d0
    phase_units(51)      = 1d-12
    is_PD(51) = 0
    nTB(51) = 0

    ! (49):  OH + CO <=> H + CO2
    fwd_A(52)     = 47600000d0
    fwd_beta(52)  = 1.228d0
    fwd_Ea(52)    = 70d0
    prefactor_units(52)  = 1.0000000000000002d-06
    activation_units(52) = 0.50321666580471969d0
    phase_units(52)      = 1d-12
    is_PD(52) = 0
    nTB(52) = 0

    ! (50):  OH + HCO <=> H2O + CO
    fwd_A(53)     = 50000000000000d0
    fwd_beta(53)  = 0d0
    fwd_Ea(53)    = 0d0
    prefactor_units(53)  = 1.0000000000000002d-06
    activation_units(53) = 0.50321666580471969d0
    phase_units(53)      = 1d-12
    is_PD(53) = 0
    nTB(53) = 0

    ! (51):  OH + CH2O <=> HCO + H2O
    fwd_A(54)     = 3430000000d0
    fwd_beta(54)  = 1.1799999999999999d0
    fwd_Ea(54)    = -447d0
    prefactor_units(54)  = 1.0000000000000002d-06
    activation_units(54) = 0.50321666580471969d0
    phase_units(54)      = 1d-12
    is_PD(54) = 0
    nTB(54) = 0

    ! (52):  OH + C2H6 <=> C2H5 + H2O
    fwd_A(55)     = 3540000d0
    fwd_beta(55)  = 2.1200000000000001d0
    fwd_Ea(55)    = 870d0
    prefactor_units(55)  = 1.0000000000000002d-06
    activation_units(55) = 0.50321666580471969d0
    phase_units(55)      = 1d-12
    is_PD(55) = 0
    nTB(55) = 0

    ! (53):  HO2 + CH2 <=> OH + CH2O
    fwd_A(56)     = 20000000000000d0
    fwd_beta(56)  = 0d0
    fwd_Ea(56)    = 0d0
    prefactor_units(56)  = 1.0000000000000002d-06
    activation_units(56) = 0.50321666580471969d0
    phase_units(56)      = 1d-12
    is_PD(56) = 0
    nTB(56) = 0

    ! (54):  HO2 + CH3 <=> O2 + CH4
    fwd_A(57)     = 1000000000000d0
    fwd_beta(57)  = 0d0
    fwd_Ea(57)    = 0d0
    prefactor_units(57)  = 1.0000000000000002d-06
    activation_units(57) = 0.50321666580471969d0
    phase_units(57)      = 1d-12
    is_PD(57) = 0
    nTB(57) = 0

    ! (55):  HO2 + CH3 <=> OH + CH3O
    fwd_A(58)     = 20000000000000d0
    fwd_beta(58)  = 0d0
    fwd_Ea(58)    = 0d0
    prefactor_units(58)  = 1.0000000000000002d-06
    activation_units(58) = 0.50321666580471969d0
    phase_units(58)      = 1d-12
    is_PD(58) = 0
    nTB(58) = 0

    ! (56):  HO2 + CO <=> OH + CO2
    fwd_A(59)     = 150000000000000d0
    fwd_beta(59)  = 0d0
    fwd_Ea(59)    = 23600d0
    prefactor_units(59)  = 1.0000000000000002d-06
    activation_units(59) = 0.50321666580471969d0
    phase_units(59)      = 1d-12
    is_PD(59) = 0
    nTB(59) = 0

    ! (57):  CH2 + O2 <=> OH + HCO
    fwd_A(60)     = 13200000000000d0
    fwd_beta(60)  = 0d0
    fwd_Ea(60)    = 1500d0
    prefactor_units(60)  = 1.0000000000000002d-06
    activation_units(60) = 0.50321666580471969d0
    phase_units(60)      = 1d-12
    is_PD(60) = 0
    nTB(60) = 0

    ! (58):  CH2 + H2 <=> H + CH3
    fwd_A(61)     = 500000d0
    fwd_beta(61)  = 2d0
    fwd_Ea(61)    = 7230d0
    prefactor_units(61)  = 1.0000000000000002d-06
    activation_units(61) = 0.50321666580471969d0
    phase_units(61)      = 1d-12
    is_PD(61) = 0
    nTB(61) = 0

    ! (59):  CH2 + CH3 <=> H + C2H4
    fwd_A(62)     = 40000000000000d0
    fwd_beta(62)  = 0d0
    fwd_Ea(62)    = 0d0
    prefactor_units(62)  = 1.0000000000000002d-06
    activation_units(62) = 0.50321666580471969d0
    phase_units(62)      = 1d-12
    is_PD(62) = 0
    nTB(62) = 0

    ! (60):  CH2 + CH4 <=> 2 CH3
    fwd_A(63)     = 2460000d0
    fwd_beta(63)  = 2d0
    fwd_Ea(63)    = 8270d0
    prefactor_units(63)  = 1.0000000000000002d-06
    activation_units(63) = 0.50321666580471969d0
    phase_units(63)      = 1d-12
    is_PD(63) = 0
    nTB(63) = 0

    ! (61):  CH2(S) + N2 <=> CH2 + N2
    fwd_A(64)     = 15000000000000d0
    fwd_beta(64)  = 0d0
    fwd_Ea(64)    = 600d0
    prefactor_units(64)  = 1.0000000000000002d-06
    activation_units(64) = 0.50321666580471969d0
    phase_units(64)      = 1d-12
    is_PD(64) = 0
    nTB(64) = 0

    ! (62):  CH2(S) + AR <=> CH2 + AR
    fwd_A(65)     = 9000000000000d0
    fwd_beta(65)  = 0d0
    fwd_Ea(65)    = 600d0
    prefactor_units(65)  = 1.0000000000000002d-06
    activation_units(65) = 0.50321666580471969d0
    phase_units(65)      = 1d-12
    is_PD(65) = 0
    nTB(65) = 0

    ! (63):  CH2(S) + O2 <=> H + OH + CO
    fwd_A(66)     = 28000000000000d0
    fwd_beta(66)  = 0d0
    fwd_Ea(66)    = 0d0
    prefactor_units(66)  = 1.0000000000000002d-06
    activation_units(66) = 0.50321666580471969d0
    phase_units(66)      = 1d-12
    is_PD(66) = 0
    nTB(66) = 0

    ! (64):  CH2(S) + O2 <=> CO + H2O
    fwd_A(67)     = 12000000000000d0
    fwd_beta(67)  = 0d0
    fwd_Ea(67)    = 0d0
    prefactor_units(67)  = 1.0000000000000002d-06
    activation_units(67) = 0.50321666580471969d0
    phase_units(67)      = 1d-12
    is_PD(67) = 0
    nTB(67) = 0

    ! (65):  CH2(S) + H2 <=> CH3 + H
    fwd_A(68)     = 70000000000000d0
    fwd_beta(68)  = 0d0
    fwd_Ea(68)    = 0d0
    prefactor_units(68)  = 1.0000000000000002d-06
    activation_units(68) = 0.50321666580471969d0
    phase_units(68)      = 1d-12
    is_PD(68) = 0
    nTB(68) = 0

    ! (66):  CH2(S) + H2O <=> CH2 + H2O
    fwd_A(69)     = 30000000000000d0
    fwd_beta(69)  = 0d0
    fwd_Ea(69)    = 0d0
    prefactor_units(69)  = 1.0000000000000002d-06
    activation_units(69) = 0.50321666580471969d0
    phase_units(69)      = 1d-12
    is_PD(69) = 0
    nTB(69) = 0

    ! (67):  CH2(S) + CH3 <=> H + C2H4
    fwd_A(70)     = 12000000000000d0
    fwd_beta(70)  = 0d0
    fwd_Ea(70)    = -570d0
    prefactor_units(70)  = 1.0000000000000002d-06
    activation_units(70) = 0.50321666580471969d0
    phase_units(70)      = 1d-12
    is_PD(70) = 0
    nTB(70) = 0

    ! (68):  CH2(S) + CH4 <=> 2 CH3
    fwd_A(71)     = 16000000000000d0
    fwd_beta(71)  = 0d0
    fwd_Ea(71)    = -570d0
    prefactor_units(71)  = 1.0000000000000002d-06
    activation_units(71) = 0.50321666580471969d0
    phase_units(71)      = 1d-12
    is_PD(71) = 0
    nTB(71) = 0

    ! (69):  CH2(S) + CO <=> CH2 + CO
    fwd_A(72)     = 9000000000000d0
    fwd_beta(72)  = 0d0
    fwd_Ea(72)    = 0d0
    prefactor_units(72)  = 1.0000000000000002d-06
    activation_units(72) = 0.50321666580471969d0
    phase_units(72)      = 1d-12
    is_PD(72) = 0
    nTB(72) = 0

    ! (70):  CH2(S) + CO2 <=> CH2 + CO2
    fwd_A(73)     = 7000000000000d0
    fwd_beta(73)  = 0d0
    fwd_Ea(73)    = 0d0
    prefactor_units(73)  = 1.0000000000000002d-06
    activation_units(73) = 0.50321666580471969d0
    phase_units(73)      = 1d-12
    is_PD(73) = 0
    nTB(73) = 0

    ! (71):  CH2(S) + CO2 <=> CO + CH2O
    fwd_A(74)     = 14000000000000d0
    fwd_beta(74)  = 0d0
    fwd_Ea(74)    = 0d0
    prefactor_units(74)  = 1.0000000000000002d-06
    activation_units(74) = 0.50321666580471969d0
    phase_units(74)      = 1d-12
    is_PD(74) = 0
    nTB(74) = 0

    ! (72):  CH3 + O2 <=> O + CH3O
    fwd_A(75)     = 26750000000000d0
    fwd_beta(75)  = 0d0
    fwd_Ea(75)    = 28800d0
    prefactor_units(75)  = 1.0000000000000002d-06
    activation_units(75) = 0.50321666580471969d0
    phase_units(75)      = 1d-12
    is_PD(75) = 0
    nTB(75) = 0

    ! (73):  CH3 + O2 <=> OH + CH2O
    fwd_A(76)     = 36000000000d0
    fwd_beta(76)  = 0d0
    fwd_Ea(76)    = 8940d0
    prefactor_units(76)  = 1.0000000000000002d-06
    activation_units(76) = 0.50321666580471969d0
    phase_units(76)      = 1d-12
    is_PD(76) = 0
    nTB(76) = 0

    ! (74):  2 CH3 (+M) <=> C2H6 (+M)
    fwd_A(8)     = 21200000000000000d0
    fwd_beta(8)  = -0.96999999999999997d0
    fwd_Ea(8)    = 620d0
    low_A(8)     = 1.7700000000000001d+50
    low_beta(8)  = -9.6699999999999999d0
    low_Ea(8)    = 6220d0
    troe_a(8)    = 0.53249999999999997d0
    troe_Tsss(8) = 151d0
    troe_Ts(8)   = 1038d0
    troe_Tss(8)  = 4970d0
    troe_len(8)  = 4
    prefactor_units(8)  = 1.0000000000000002d-06
    activation_units(8) = 0.50321666580471969d0
    phase_units(8)      = 1d-12
    is_PD(8) = 1
    nTB(8) = 7
    if (.not. allocated(TB(8) % vector)) allocate(TB(8) % vector(7))
    if (.not. allocated(TBid(8) % vector)) allocate(TBid(8) % vector(7))
    TBid(8) % vector(1) = 0d0
    TB(8) % vector(1) = 2d0 ! H2
    TBid(8) % vector(2) = 5d0
    TB(8) % vector(2) = 6d0 ! H2O
    TBid(8) % vector(3) = 10d0
    TB(8) % vector(3) = 2d0 ! CH4
    TBid(8) % vector(4) = 11d0
    TB(8) % vector(4) = 1.5d0 ! CO
    TBid(8) % vector(5) = 12d0
    TB(8) % vector(5) = 2d0 ! CO2
    TBid(8) % vector(6) = 18d0
    TB(8) % vector(6) = 3d0 ! C2H6
    TBid(8) % vector(7) = 20d0
    TB(8) % vector(7) = 0.69999999999999996d0 ! AR

    ! (75):  2 CH3 <=> H + C2H5
    fwd_A(77)     = 4990000000000d0
    fwd_beta(77)  = 0.10000000000000001d0
    fwd_Ea(77)    = 10600d0
    prefactor_units(77)  = 1.0000000000000002d-06
    activation_units(77) = 0.50321666580471969d0
    phase_units(77)      = 1d-12
    is_PD(77) = 0
    nTB(77) = 0

    ! (76):  CH3 + HCO <=> CH4 + CO
    fwd_A(78)     = 26480000000000d0
    fwd_beta(78)  = 0d0
    fwd_Ea(78)    = 0d0
    prefactor_units(78)  = 1.0000000000000002d-06
    activation_units(78) = 0.50321666580471969d0
    phase_units(78)      = 1d-12
    is_PD(78) = 0
    nTB(78) = 0

    ! (77):  CH3 + CH2O <=> HCO + CH4
    fwd_A(79)     = 3320d0
    fwd_beta(79)  = 2.8100000000000001d0
    fwd_Ea(79)    = 5860d0
    prefactor_units(79)  = 1.0000000000000002d-06
    activation_units(79) = 0.50321666580471969d0
    phase_units(79)      = 1d-12
    is_PD(79) = 0
    nTB(79) = 0

    ! (78):  CH3 + C2H6 <=> C2H5 + CH4
    fwd_A(80)     = 6140000d0
    fwd_beta(80)  = 1.74d0
    fwd_Ea(80)    = 10450d0
    prefactor_units(80)  = 1.0000000000000002d-06
    activation_units(80) = 0.50321666580471969d0
    phase_units(80)      = 1d-12
    is_PD(80) = 0
    nTB(80) = 0

    ! (79):  HCO + H2O <=> H + CO + H2O
    fwd_A(81)     = 2.244d+18
    fwd_beta(81)  = -1d0
    fwd_Ea(81)    = 17000d0
    prefactor_units(81)  = 1.0000000000000002d-06
    activation_units(81) = 0.50321666580471969d0
    phase_units(81)      = 1d-12
    is_PD(81) = 0
    nTB(81) = 0

    ! (80):  HCO + M <=> H + CO + M
    fwd_A(14)     = 1.87d+17
    fwd_beta(14)  = -1d0
    fwd_Ea(14)    = 17000d0
    prefactor_units(14)  = 1.0000000000000002d-06
    activation_units(14) = 0.50321666580471969d0
    phase_units(14)      = 1d-6
    is_PD(14) = 0
    nTB(14) = 6
    if (.not. allocated(TB(14) % vector)) allocate(TB(14) % vector(6))
    if (.not. allocated(TBid(14) % vector)) allocate(TBid(14) % vector(6))
    TBid(14) % vector(1) = 0d0
    TB(14) % vector(1) = 2d0 ! H2
    TBid(14) % vector(2) = 5d0
    TB(14) % vector(2) = 0d0 ! H2O
    TBid(14) % vector(3) = 10d0
    TB(14) % vector(3) = 2d0 ! CH4
    TBid(14) % vector(4) = 11d0
    TB(14) % vector(4) = 1.5d0 ! CO
    TBid(14) % vector(5) = 12d0
    TB(14) % vector(5) = 2d0 ! CO2
    TBid(14) % vector(6) = 18d0
    TB(14) % vector(6) = 3d0 ! C2H6

    ! (81):  HCO + O2 <=> HO2 + CO
    fwd_A(82)     = 7600000000000d0
    fwd_beta(82)  = 0d0
    fwd_Ea(82)    = 400d0
    prefactor_units(82)  = 1.0000000000000002d-06
    activation_units(82) = 0.50321666580471969d0
    phase_units(82)      = 1d-12
    is_PD(82) = 0
    nTB(82) = 0

    ! (82):  CH3O + O2 <=> HO2 + CH2O
    fwd_A(83)     = 4.2799999999999999d-13
    fwd_beta(83)  = 7.5999999999999996d0
    fwd_Ea(83)    = -3530d0
    prefactor_units(83)  = 1.0000000000000002d-06
    activation_units(83) = 0.50321666580471969d0
    phase_units(83)      = 1d-12
    is_PD(83) = 0
    nTB(83) = 0

    ! (83):  C2H5 + O2 <=> HO2 + C2H4
    fwd_A(84)     = 840000000000d0
    fwd_beta(84)  = 0d0
    fwd_Ea(84)    = 3875d0
    prefactor_units(84)  = 1.0000000000000002d-06
    activation_units(84) = 0.50321666580471969d0
    phase_units(84)      = 1d-12
    is_PD(84) = 0
    nTB(84) = 0

    call SetAllDefaults()

end subroutine


! A few mechanism parameters
subroutine ckindx(iwrk, rwrk, mm, kk, ii, nfit)

    implicit none

    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    integer, intent(out) :: mm
    integer, intent(out) :: kk
    integer, intent(out) :: ii
    integer, intent(out) :: nfit

    mm = 5
    kk = 21
    ii = 84
    nfit = -1 ! Why do you need this anyway?

end subroutine

! Returns the char strings of element names
subroutine cksyme(kname, plenkname)

    implicit none

    integer, intent(out) :: kname(plenkname*5)
    integer, intent(in) :: plenkname

    integer :: i
    integer :: lenkname

    lenkname = plenkname

    !clear kname
    do i=1, lenkname*5
        kname(i) = ichar(' ')
    end do

    ! O 
    kname(0*lenkname+1) = ichar('O')
    kname(0*lenkname+2) = ichar(' ')
    ! H 
    kname(1*lenkname+1) = ichar('H')
    kname(1*lenkname+2) = ichar(' ')
    ! C 
    kname(2*lenkname+1) = ichar('C')
    kname(2*lenkname+2) = ichar(' ')
    ! N 
    kname(3*lenkname+1) = ichar('N')
    kname(3*lenkname+2) = ichar(' ')
    ! AR 
    kname(4*lenkname+1) = ichar('A')
    kname(4*lenkname+2) = ichar('R')
    kname(4*lenkname+3) = ichar(' ')

end subroutine

! Returns the char strings of species names
subroutine cksyms(kname, plenkname)

    implicit none

    integer, intent(out) :: kname(plenkname*21)
    integer, intent(in) :: plenkname

    integer :: i
    integer :: lenkname

    lenkname = plenkname

    !clear kname
    do i=1, lenkname*21
        kname(i) = ichar(' ')
    end do

    ! H2 
    kname(0*lenkname+1) = ichar('H')
    kname(0*lenkname+2) = ichar('2')
    kname(0*lenkname+3) = ichar(' ')
    ! H 
    kname(1*lenkname+1) = ichar('H')
    kname(1*lenkname+2) = ichar(' ')
    ! O 
    kname(2*lenkname+1) = ichar('O')
    kname(2*lenkname+2) = ichar(' ')
    ! O2 
    kname(3*lenkname+1) = ichar('O')
    kname(3*lenkname+2) = ichar('2')
    kname(3*lenkname+3) = ichar(' ')
    ! OH 
    kname(4*lenkname+1) = ichar('O')
    kname(4*lenkname+2) = ichar('H')
    kname(4*lenkname+3) = ichar(' ')
    ! H2O 
    kname(5*lenkname+1) = ichar('H')
    kname(5*lenkname+2) = ichar('2')
    kname(5*lenkname+3) = ichar('O')
    kname(5*lenkname+4) = ichar(' ')
    ! HO2 
    kname(6*lenkname+1) = ichar('H')
    kname(6*lenkname+2) = ichar('O')
    kname(6*lenkname+3) = ichar('2')
    kname(6*lenkname+4) = ichar(' ')
    ! CH2 
    kname(7*lenkname+1) = ichar('C')
    kname(7*lenkname+2) = ichar('H')
    kname(7*lenkname+3) = ichar('2')
    kname(7*lenkname+4) = ichar(' ')
    ! CH2(S) 
    kname(8*lenkname+1) = ichar('C')
    kname(8*lenkname+2) = ichar('H')
    kname(8*lenkname+3) = ichar('2')
    kname(8*lenkname+4) = ichar('(')
    kname(8*lenkname+5) = ichar('S')
    kname(8*lenkname+6) = ichar(')')
    kname(8*lenkname+7) = ichar(' ')
    ! CH3 
    kname(9*lenkname+1) = ichar('C')
    kname(9*lenkname+2) = ichar('H')
    kname(9*lenkname+3) = ichar('3')
    kname(9*lenkname+4) = ichar(' ')
    ! CH4 
    kname(10*lenkname+1) = ichar('C')
    kname(10*lenkname+2) = ichar('H')
    kname(10*lenkname+3) = ichar('4')
    kname(10*lenkname+4) = ichar(' ')
    ! CO 
    kname(11*lenkname+1) = ichar('C')
    kname(11*lenkname+2) = ichar('O')
    kname(11*lenkname+3) = ichar(' ')
    ! CO2 
    kname(12*lenkname+1) = ichar('C')
    kname(12*lenkname+2) = ichar('O')
    kname(12*lenkname+3) = ichar('2')
    kname(12*lenkname+4) = ichar(' ')
    ! HCO 
    kname(13*lenkname+1) = ichar('H')
    kname(13*lenkname+2) = ichar('C')
    kname(13*lenkname+3) = ichar('O')
    kname(13*lenkname+4) = ichar(' ')
    ! CH2O 
    kname(14*lenkname+1) = ichar('C')
    kname(14*lenkname+2) = ichar('H')
    kname(14*lenkname+3) = ichar('2')
    kname(14*lenkname+4) = ichar('O')
    kname(14*lenkname+5) = ichar(' ')
    ! CH3O 
    kname(15*lenkname+1) = ichar('C')
    kname(15*lenkname+2) = ichar('H')
    kname(15*lenkname+3) = ichar('3')
    kname(15*lenkname+4) = ichar('O')
    kname(15*lenkname+5) = ichar(' ')
    ! C2H4 
    kname(16*lenkname+1) = ichar('C')
    kname(16*lenkname+2) = ichar('2')
    kname(16*lenkname+3) = ichar('H')
    kname(16*lenkname+4) = ichar('4')
    kname(16*lenkname+5) = ichar(' ')
    ! C2H5 
    kname(17*lenkname+1) = ichar('C')
    kname(17*lenkname+2) = ichar('2')
    kname(17*lenkname+3) = ichar('H')
    kname(17*lenkname+4) = ichar('5')
    kname(17*lenkname+5) = ichar(' ')
    ! C2H6 
    kname(18*lenkname+1) = ichar('C')
    kname(18*lenkname+2) = ichar('2')
    kname(18*lenkname+3) = ichar('H')
    kname(18*lenkname+4) = ichar('6')
    kname(18*lenkname+5) = ichar(' ')
    ! N2 
    kname(19*lenkname+1) = ichar('N')
    kname(19*lenkname+2) = ichar('2')
    kname(19*lenkname+3) = ichar(' ')
    ! AR 
    kname(20*lenkname+1) = ichar('A')
    kname(20*lenkname+2) = ichar('R')
    kname(20*lenkname+3) = ichar(' ')

end subroutine

! Returns R, Rc, Patm
subroutine ckrp(ickwrk, rckwrk, ru, ruc, pa)

    implicit none

    integer, intent(in) :: ickwrk
    double precision, intent(in) :: rckwrk
    double precision, intent(out) :: ru
    double precision, intent(out) :: ruc
    double precision, intent(out) :: pa

    ru  = 83145100.000000d0 
    ruc = 1.98721558317399615845d0 
    pa  = 1013250.000000d0 

end subroutine

! Compute P = rhoRT/W(y)
subroutine ckpy(rho, T, y, iwrk, rwrk, P)

    implicit none

    double precision, intent(in) :: rho
    double precision, intent(in) :: T
    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: P

    double precision :: YOW ! for computing mean MW

    YOW = 0.d0

    YOW = YOW + (y(1) * imw(1)) ! H2
    YOW = YOW + (y(2) * imw(2)) ! H
    YOW = YOW + (y(3) * imw(3)) ! O
    YOW = YOW + (y(4) * imw(4)) ! O2
    YOW = YOW + (y(5) * imw(5)) ! OH
    YOW = YOW + (y(6) * imw(6)) ! H2O
    YOW = YOW + (y(7) * imw(7)) ! HO2
    YOW = YOW + (y(8) * imw(8)) ! CH2
    YOW = YOW + (y(9) * imw(9)) ! CH2(S)
    YOW = YOW + (y(10) * imw(10)) ! CH3
    YOW = YOW + (y(11) * imw(11)) ! CH4
    YOW = YOW + (y(12) * imw(12)) ! CO
    YOW = YOW + (y(13) * imw(13)) ! CO2
    YOW = YOW + (y(14) * imw(14)) ! HCO
    YOW = YOW + (y(15) * imw(15)) ! CH2O
    YOW = YOW + (y(16) * imw(16)) ! CH3O
    YOW = YOW + (y(17) * imw(17)) ! C2H4
    YOW = YOW + (y(18) * imw(18)) ! C2H5
    YOW = YOW + (y(19) * imw(19)) ! C2H6
    YOW = YOW + (y(20) * imw(20)) ! N2
    YOW = YOW + (y(21) * imw(21)) ! AR

    ! YOW holds the reciprocal of the mean molecular wt
    P = rho * 8.31451000d+07 * T * YOW ! P = rho*R*T/W

end subroutine

! Compute rho = P*W(y)/RT
subroutine ckrhoy(P, T, y, iwrk, rwrk, rho)

    implicit none

    double precision, intent(in) :: P
    double precision, intent(in) :: T
    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: rho

    double precision :: YOW, tmp(21)
    integer :: i

    YOW = 0.d0

    do i=1, 21
        tmp(i) = y(i) * imw(i)
    end do
    do i=1, 21
        YOW = YOW + tmp(i)
    end do

    rho = P / ( 8.31451000d+07 * T * YOW) ! rho = P*W/(R*T)

end subroutine

! get molecular weight for all species
subroutine ckwt(iwrk, rwrk, wt)

    implicit none

    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: wt(21)

    call molecularWeight(wt)

end subroutine

! convert y[species] (mass fracs) to x[species] (mole fracs)
subroutine ckytx(y, iwrk, rwrk, x)

    implicit none

    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: x(21)

    double precision :: YOW, YOWINV
    integer :: i

    YOW = 0.d0

    do i=1, 21
        YOW = YOW + y(i) * imw(i)
    end do

    YOWINV = 1.d0 / YOW

    do i=1, 21
        x(i) = y(i) * imw(i) * YOWINV
    end do

end subroutine

! convert y(npoints,species) (mass fracs) to x(npoints,species) (mole fracs)
subroutine vckytx(np, y, iwrk, rwrk, x)

    implicit none

    integer, intent(in) :: np
    double precision, intent(in) :: y(np,21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: x(np,21)

    double precision :: YOW(np)
    integer :: i, n

    do i=1, np
        YOW(i) = 0.d0
    end do

    do n=1, 21
        do i=1, np
            x(i,n) = y(i,n) * imw(n)
            YOW(i) = YOW(i) + x(i,n)
        end do
    end do

    do i=1, np
        YOW(i) = 1.d0 / YOW(i)
    end do

    do n=1, 21
        do i=1, np
            x(i,n) = x(i,n) * YOW(i)
        end do
    end do

end subroutine

! convert y[species] (mass fracs) to c[species] (molar conc)
subroutine ckytcr(rho, T, y, iwrk, rwrk, c)

    implicit none

    double precision, intent(in) :: rho
    double precision, intent(in) :: T
    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: c(21)

    integer :: i

    do i=1, 21
        c(i) = rho * y(i) * imw(i)
    end do

end subroutine

! convert x[species] (mole fracs) to y[species] (mass fracs)
subroutine ckxty(x, iwrk, rwrk, y)

    implicit none

    double precision, intent(in) :: x(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: y(21)

    double precision :: XW, XWinv

    XW = 0.d0 ! See Eq 4, 9 in CK Manual

    ! Compute mean molecular wt first
    XW = XW + (x(1) * 2.01594000d+00) ! H2
    XW = XW + (x(2) * 1.00797000d+00) ! H
    XW = XW + (x(3) * 1.59994000d+01) ! O
    XW = XW + (x(4) * 3.19988000d+01) ! O2
    XW = XW + (x(5) * 1.70073700d+01) ! OH
    XW = XW + (x(6) * 1.80153400d+01) ! H2O
    XW = XW + (x(7) * 3.30067700d+01) ! HO2
    XW = XW + (x(8) * 1.40270900d+01) ! CH2
    XW = XW + (x(9) * 1.40270900d+01) ! CH2(S)
    XW = XW + (x(10) * 1.50350600d+01) ! CH3
    XW = XW + (x(11) * 1.60430300d+01) ! CH4
    XW = XW + (x(12) * 2.80105500d+01) ! CO
    XW = XW + (x(13) * 4.40099500d+01) ! CO2
    XW = XW + (x(14) * 2.90185200d+01) ! HCO
    XW = XW + (x(15) * 3.00264900d+01) ! CH2O
    XW = XW + (x(16) * 3.10344600d+01) ! CH3O
    XW = XW + (x(17) * 2.80541800d+01) ! C2H4
    XW = XW + (x(18) * 2.90621500d+01) ! C2H5
    XW = XW + (x(19) * 3.00701200d+01) ! C2H6
    XW = XW + (x(20) * 2.80134000d+01) ! N2
    XW = XW + (x(21) * 3.99480000d+01) ! AR

    ! Now compute conversion
    XWinv = 1.d0 / XW
    y(1) = x(1) * 2.01594000d+00 * XWinv 
    y(2) = x(2) * 1.00797000d+00 * XWinv 
    y(3) = x(3) * 1.59994000d+01 * XWinv 
    y(4) = x(4) * 3.19988000d+01 * XWinv 
    y(5) = x(5) * 1.70073700d+01 * XWinv 
    y(6) = x(6) * 1.80153400d+01 * XWinv 
    y(7) = x(7) * 3.30067700d+01 * XWinv 
    y(8) = x(8) * 1.40270900d+01 * XWinv 
    y(9) = x(9) * 1.40270900d+01 * XWinv 
    y(10) = x(10) * 1.50350600d+01 * XWinv 
    y(11) = x(11) * 1.60430300d+01 * XWinv 
    y(12) = x(12) * 2.80105500d+01 * XWinv 
    y(13) = x(13) * 4.40099500d+01 * XWinv 
    y(14) = x(14) * 2.90185200d+01 * XWinv 
    y(15) = x(15) * 3.00264900d+01 * XWinv 
    y(16) = x(16) * 3.10344600d+01 * XWinv 
    y(17) = x(17) * 2.80541800d+01 * XWinv 
    y(18) = x(18) * 2.90621500d+01 * XWinv 
    y(19) = x(19) * 3.00701200d+01 * XWinv 
    y(20) = x(20) * 2.80134000d+01 * XWinv 
    y(21) = x(21) * 3.99480000d+01 * XWinv 

end subroutine

! Returns the specific heats at constant volume
! in mass units (Eq. 29)
subroutine ckcvms(T, iwrk, rwrk, cvms)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: cvms(21)

    double precision :: tT, tc(5)

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cv_R(cvms, tc)

    ! multiply by R/molecularweight
    cvms(1) = cvms(1) * 4.124383662212169d+07 !H2
    cvms(2) = cvms(2) * 8.248767324424338d+07 !H
    cvms(3) = cvms(3) * 5.196763628636074d+06 !O
    cvms(4) = cvms(4) * 2.598381814318037d+06 !O2
    cvms(5) = cvms(5) * 4.888768810227566d+06 !OH
    cvms(6) = cvms(6) * 4.615239012974499d+06 !H2O
    cvms(7) = cvms(7) * 2.519031701678171d+06 !HO2
    cvms(8) = cvms(8) * 5.927466067445207d+06 !CH2
    cvms(9) = cvms(9) * 5.927466067445207d+06 !CH2(S)
    cvms(10) = cvms(10) * 5.530081023953346d+06 !CH3
    cvms(11) = cvms(11) * 5.182630712527496d+06 !CH4
    cvms(12) = cvms(12) * 2.968349425484326d+06 !CO
    cvms(13) = cvms(13) * 1.889234139098090d+06 !CO2
    cvms(14) = cvms(14) * 2.865242610581105d+06 !HCO
    cvms(15) = cvms(15) * 2.769058254894261d+06 !CH2O
    cvms(16) = cvms(16) * 2.679121853578248d+06 !CH3O
    cvms(17) = cvms(17) * 2.963733033722604d+06 !C2H4
    cvms(18) = cvms(18) * 2.860941121011349d+06 !C2H5
    cvms(19) = cvms(19) * 2.765040511976673d+06 !C2H6
    cvms(20) = cvms(20) * 2.968047434442088d+06 !N2
    cvms(21) = cvms(21) * 2.081333233203164d+06 !AR

end subroutine

! Returns the specific heats at constant pressure
! in mass units (Eq. 26)
subroutine ckcpms(T, iwrk, rwrk, cpms)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: cpms(21)

    double precision :: tT, tc(5)

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cp_R(cpms, tc)

    ! multiply by R/molecularweight
    cpms(1) = cpms(1) * 4.124383662212169d+07 ! H2
    cpms(2) = cpms(2) * 8.248767324424338d+07 ! H
    cpms(3) = cpms(3) * 5.196763628636074d+06 ! O
    cpms(4) = cpms(4) * 2.598381814318037d+06 ! O2
    cpms(5) = cpms(5) * 4.888768810227566d+06 ! OH
    cpms(6) = cpms(6) * 4.615239012974499d+06 ! H2O
    cpms(7) = cpms(7) * 2.519031701678171d+06 ! HO2
    cpms(8) = cpms(8) * 5.927466067445207d+06 ! CH2
    cpms(9) = cpms(9) * 5.927466067445207d+06 ! CH2(S)
    cpms(10) = cpms(10) * 5.530081023953346d+06 ! CH3
    cpms(11) = cpms(11) * 5.182630712527496d+06 ! CH4
    cpms(12) = cpms(12) * 2.968349425484326d+06 ! CO
    cpms(13) = cpms(13) * 1.889234139098090d+06 ! CO2
    cpms(14) = cpms(14) * 2.865242610581105d+06 ! HCO
    cpms(15) = cpms(15) * 2.769058254894261d+06 ! CH2O
    cpms(16) = cpms(16) * 2.679121853578248d+06 ! CH3O
    cpms(17) = cpms(17) * 2.963733033722604d+06 ! C2H4
    cpms(18) = cpms(18) * 2.860941121011349d+06 ! C2H5
    cpms(19) = cpms(19) * 2.765040511976673d+06 ! C2H6
    cpms(20) = cpms(20) * 2.968047434442088d+06 ! N2
    cpms(21) = cpms(21) * 2.081333233203164d+06 ! AR

end subroutine

! Returns internal energy in mass units (Eq 30.)
subroutine ckums(T, iwrk, rwrk, ums)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: ums(21)

    double precision :: tT, tc(5)
    double precision :: RT
    integer :: i

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesInternalEnergy(ums, tc)

    do i=1, 21
        ums(i) = ums(i) * (RT * imw(i))
    end do

end subroutine

! Returns enthalpy in mass units (Eq 27.)
subroutine ckhms(T, iwrk, rwrk, hms)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: hms(21)

    double precision :: tT, RT
    double precision :: tc(5), h(21)
    integer :: i

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesEnthalpy(hms, tc)

    do i=1, 21
        hms(i) = hms(i) * (RT * imw(i))
    end do

end subroutine

! Returns enthalpy in mass units (Eq 27.)
subroutine vckhms(np, T, iwrk, rwrk, hms)

    implicit none

    integer, intent(in) :: np
    double precision, intent(in) :: T(np)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: hms(np,21)

    double precision :: tc(5), h(21)
    integer :: i, n

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
        hms(i, 10) = h(10)
        hms(i, 11) = h(11)
        hms(i, 12) = h(12)
        hms(i, 13) = h(13)
        hms(i, 14) = h(14)
        hms(i, 15) = h(15)
        hms(i, 16) = h(16)
        hms(i, 17) = h(17)
        hms(i, 18) = h(18)
        hms(i, 19) = h(19)
        hms(i, 20) = h(20)
        hms(i, 21) = h(21)
    end do

    do n=1, 21
        do i=1, np
            hms(i,n) = hms(i,n) * ( 8.31451000d+07 * T(i) * imw(n))
        end do
    end do

end subroutine

! Returns the mean specific heat at CP (Eq. 34)
subroutine ckcpbs(T, y, iwrk, rwrk, cpbs)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: cpbs

    double precision :: cpor(21)
    double precision :: tresult(21)
    double precision :: tT, tc(5)
    double precision :: res
    integer :: i

    res = 0.d0

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cp_R(cpor, tc)

    do i=1, 21
        tresult(i) = cpor(i) * y(i) * imw(i)
    end do
    do i=1, 21
        res = res + tresult(i)
    end do

    cpbs = res * 8.31451000d+07

end subroutine

! Returns the mean specific heat at CV (Eq. 36)
subroutine ckcvbs(T, y, iwrk, rwrk, cvbs)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: cvbs

    double precision :: cvor(21)
    double precision :: tT, tc(5)
    double precision :: res

    res = 0.d0
    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cv_R(cvor, tc)

    ! multiply by y/molecularweight
    res = res + (cvor(1) * y(1) * imw(1)) ! H2
    res = res + (cvor(2) * y(2) * imw(2)) ! H
    res = res + (cvor(3) * y(3) * imw(3)) ! O
    res = res + (cvor(4) * y(4) * imw(4)) ! O2
    res = res + (cvor(5) * y(5) * imw(5)) ! OH
    res = res + (cvor(6) * y(6) * imw(6)) ! H2O
    res = res + (cvor(7) * y(7) * imw(7)) ! HO2
    res = res + (cvor(8) * y(8) * imw(8)) ! CH2
    res = res + (cvor(9) * y(9) * imw(9)) ! CH2(S)
    res = res + (cvor(10) * y(10) * imw(10)) ! CH3
    res = res + (cvor(11) * y(11) * imw(11)) ! CH4
    res = res + (cvor(12) * y(12) * imw(12)) ! CO
    res = res + (cvor(13) * y(13) * imw(13)) ! CO2
    res = res + (cvor(14) * y(14) * imw(14)) ! HCO
    res = res + (cvor(15) * y(15) * imw(15)) ! CH2O
    res = res + (cvor(16) * y(16) * imw(16)) ! CH3O
    res = res + (cvor(17) * y(17) * imw(17)) ! C2H4
    res = res + (cvor(18) * y(18) * imw(18)) ! C2H5
    res = res + (cvor(19) * y(19) * imw(19)) ! C2H6
    res = res + (cvor(20) * y(20) * imw(20)) ! N2
    res = res + (cvor(21) * y(21) * imw(21)) ! AR

    cvbs = res * 8.31451000d+07

end subroutine

! get mean internal energy in mass units
subroutine ckubms(T, y, iwrk, rwrk, ubms)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: ubms

    double precision :: ums(21) ! temporary energy array
    double precision :: res
    double precision :: RT, tT, tc(5)

    res = 0.d0

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesInternalEnergy(ums, tc)

    ! perform dot product + scaling by wt
    res = res + (y(1) * ums(1) * imw(1)) ! H2
    res = res + (y(2) * ums(2) * imw(2)) ! H
    res = res + (y(3) * ums(3) * imw(3)) ! O
    res = res + (y(4) * ums(4) * imw(4)) ! O2
    res = res + (y(5) * ums(5) * imw(5)) ! OH
    res = res + (y(6) * ums(6) * imw(6)) ! H2O
    res = res + (y(7) * ums(7) * imw(7)) ! HO2
    res = res + (y(8) * ums(8) * imw(8)) ! CH2
    res = res + (y(9) * ums(9) * imw(9)) ! CH2(S)
    res = res + (y(10) * ums(10) * imw(10)) ! CH3
    res = res + (y(11) * ums(11) * imw(11)) ! CH4
    res = res + (y(12) * ums(12) * imw(12)) ! CO
    res = res + (y(13) * ums(13) * imw(13)) ! CO2
    res = res + (y(14) * ums(14) * imw(14)) ! HCO
    res = res + (y(15) * ums(15) * imw(15)) ! CH2O
    res = res + (y(16) * ums(16) * imw(16)) ! CH3O
    res = res + (y(17) * ums(17) * imw(17)) ! C2H4
    res = res + (y(18) * ums(18) * imw(18)) ! C2H5
    res = res + (y(19) * ums(19) * imw(19)) ! C2H6
    res = res + (y(20) * ums(20) * imw(20)) ! N2
    res = res + (y(21) * ums(21) * imw(21)) ! AR

    ubms = res * RT

end subroutine

! compute the production rate for each species
subroutine ckwc(T, C, iwrk, rwrk, wdot)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout) :: C(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: wdot(21)

    integer :: id ! loop counter

    ! convert to SI
    do id = 1, 21
        C(id) = C(id) * 1.0d6
    end do

    ! convert to chemkin units
    call productionRate(wdot, C, T)

    ! convert to chemkin units
    do id=1, 21
        C(id) = C(id) * 1.0d-6
        wdot(id) = wdot(id) * 1.0d-6
    end do

end subroutine

! compute the production rate for each species
subroutine productionRate(wdot, sc, T)

    implicit none

    double precision, intent(inout) :: wdot(21)
    double precision, intent(in) :: sc(21)
    double precision, intent(in) :: T

    double precision :: tc(5)
    double precision :: invT
    double precision :: qdot, q_f(84), q_r(84)
    integer :: i

    tc = (/ log(T), T, T*T, T*T*T, T*T*T*T /)
    invT = 1.d0 / tc(2)

    if (T /= T_save) then
        T_save = T
        call comp_k_f(tc,invT,k_f_save)
        call comp_Kc(tc,invT,Kc_save)
    end if

    call comp_qfqr(q_f, q_r, sc, tc, invT)

    do i=1, 21
        wdot(i) = 0.d0
    end do

    qdot = q_f(1)-q_r(1)
    wdot(2) = wdot(2) - qdot
    wdot(8) = wdot(8) - qdot
    wdot(10) = wdot(10) + qdot

    qdot = q_f(2)-q_r(2)
    wdot(2) = wdot(2) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(11) = wdot(11) + qdot

    qdot = q_f(3)-q_r(3)
    wdot(2) = wdot(2) - qdot
    wdot(14) = wdot(14) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(4)-q_r(4)
    wdot(2) = wdot(2) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(5)-q_r(5)
    wdot(2) = wdot(2) - qdot
    wdot(17) = wdot(17) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(6)-q_r(6)
    wdot(2) = wdot(2) - qdot
    wdot(18) = wdot(18) - qdot
    wdot(19) = wdot(19) + qdot

    qdot = q_f(7)-q_r(7)
    wdot(1) = wdot(1) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(8)-q_r(8)
    wdot(10) = wdot(10) - (2 * qdot)
    wdot(19) = wdot(19) + qdot

    qdot = q_f(9)-q_r(9)
    wdot(2) = wdot(2) - qdot
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot

    qdot = q_f(10)-q_r(10)
    wdot(3) = wdot(3) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(11)-q_r(11)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot

    qdot = q_f(12)-q_r(12)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - (2 * qdot)

    qdot = q_f(13)-q_r(13)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(14)-q_r(14)
    wdot(2) = wdot(2) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(15)-q_r(15)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot

    qdot = q_f(16)-q_r(16)
    wdot(3) = wdot(3) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(17)-q_r(17)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(8) = wdot(8) - qdot
    wdot(14) = wdot(14) + qdot

    qdot = q_f(18)-q_r(18)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(9) = wdot(9) - qdot
    wdot(14) = wdot(14) + qdot

    qdot = q_f(19)-q_r(19)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(20)-q_r(20)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(10) = wdot(10) + qdot
    wdot(11) = wdot(11) - qdot

    qdot = q_f(21)-q_r(21)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(22)-q_r(22)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(23)-q_r(23)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(14) = wdot(14) + qdot
    wdot(15) = wdot(15) - qdot

    qdot = q_f(24)-q_r(24)
    wdot(3) = wdot(3) - qdot
    wdot(10) = wdot(10) + qdot
    wdot(14) = wdot(14) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(25)-q_r(25)
    wdot(3) = wdot(3) - qdot
    wdot(10) = wdot(10) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(18) = wdot(18) - qdot

    qdot = q_f(26)-q_r(26)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(27)-q_r(27)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(28)-q_r(28)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(14) = wdot(14) + qdot
    wdot(15) = wdot(15) - qdot

    qdot = q_f(29)-q_r(29)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - (2 * qdot)
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) + qdot

    qdot = q_f(30)-q_r(30)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) + qdot

    qdot = q_f(31)-q_r(31)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(20) = wdot(20) - qdot
    wdot(20) = wdot(20) + qdot

    qdot = q_f(32)-q_r(32)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(21) = wdot(21) - qdot
    wdot(21) = wdot(21) + qdot

    qdot = q_f(33)-q_r(33)
    wdot(2) = wdot(2) - qdot
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot

    qdot = q_f(34)-q_r(34)
    wdot(1) = wdot(1) - qdot
    wdot(1) = wdot(1) + (2 * qdot)
    wdot(2) = wdot(2) - (2 * qdot)

    qdot = q_f(35)-q_r(35)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - (2 * qdot)
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(36)-q_r(36)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - (2 * qdot)
    wdot(13) = wdot(13) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(37)-q_r(37)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(38)-q_r(38)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + (2 * qdot)
    wdot(7) = wdot(7) - qdot

    qdot = q_f(39)-q_r(39)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(10) = wdot(10) + qdot
    wdot(11) = wdot(11) - qdot

    qdot = q_f(40)-q_r(40)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(12) = wdot(12) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(41)-q_r(41)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(15) = wdot(15) - qdot

    qdot = q_f(42)-q_r(42)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(10) = wdot(10) + qdot
    wdot(16) = wdot(16) - qdot

    qdot = q_f(43)-q_r(43)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(44)-q_r(44)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(45)-q_r(45)
    wdot(3) = wdot(3) + qdot
    wdot(5) = wdot(5) - (2 * qdot)
    wdot(6) = wdot(6) + qdot

    qdot = q_f(46)-q_r(46)
    wdot(4) = wdot(4) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(47)-q_r(47)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(8) = wdot(8) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(48)-q_r(48)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(9) = wdot(9) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(49)-q_r(49)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(8) = wdot(8) + qdot
    wdot(10) = wdot(10) - qdot

    qdot = q_f(50)-q_r(50)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(9) = wdot(9) + qdot
    wdot(10) = wdot(10) - qdot

    qdot = q_f(51)-q_r(51)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(10) = wdot(10) + qdot
    wdot(11) = wdot(11) - qdot

    qdot = q_f(52)-q_r(52)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(53)-q_r(53)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(54)-q_r(54)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(14) = wdot(14) + qdot
    wdot(15) = wdot(15) - qdot

    qdot = q_f(55)-q_r(55)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(56)-q_r(56)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(8) = wdot(8) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(57)-q_r(57)
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(11) = wdot(11) + qdot

    qdot = q_f(58)-q_r(58)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(59)-q_r(59)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(60)-q_r(60)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(8) = wdot(8) - qdot
    wdot(14) = wdot(14) + qdot

    qdot = q_f(61)-q_r(61)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(8) = wdot(8) - qdot
    wdot(10) = wdot(10) + qdot

    qdot = q_f(62)-q_r(62)
    wdot(2) = wdot(2) + qdot
    wdot(8) = wdot(8) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(17) = wdot(17) + qdot

    qdot = q_f(63)-q_r(63)
    wdot(8) = wdot(8) - qdot
    wdot(10) = wdot(10) + (2 * qdot)
    wdot(11) = wdot(11) - qdot

    qdot = q_f(64)-q_r(64)
    wdot(8) = wdot(8) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(20) = wdot(20) - qdot
    wdot(20) = wdot(20) + qdot

    qdot = q_f(65)-q_r(65)
    wdot(8) = wdot(8) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(21) = wdot(21) - qdot
    wdot(21) = wdot(21) + qdot

    qdot = q_f(66)-q_r(66)
    wdot(2) = wdot(2) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(12) = wdot(12) + qdot

    qdot = q_f(67)-q_r(67)
    wdot(4) = wdot(4) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(12) = wdot(12) + qdot

    qdot = q_f(68)-q_r(68)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(10) = wdot(10) + qdot

    qdot = q_f(69)-q_r(69)
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(8) = wdot(8) + qdot
    wdot(9) = wdot(9) - qdot

    qdot = q_f(70)-q_r(70)
    wdot(2) = wdot(2) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(17) = wdot(17) + qdot

    qdot = q_f(71)-q_r(71)
    wdot(9) = wdot(9) - qdot
    wdot(10) = wdot(10) + (2 * qdot)
    wdot(11) = wdot(11) - qdot

    qdot = q_f(72)-q_r(72)
    wdot(8) = wdot(8) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(12) = wdot(12) + qdot

    qdot = q_f(73)-q_r(73)
    wdot(8) = wdot(8) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(74)-q_r(74)
    wdot(9) = wdot(9) - qdot
    wdot(12) = wdot(12) + qdot
    wdot(13) = wdot(13) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(75)-q_r(75)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(76)-q_r(76)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(77)-q_r(77)
    wdot(2) = wdot(2) + qdot
    wdot(10) = wdot(10) - (2 * qdot)
    wdot(18) = wdot(18) + qdot

    qdot = q_f(78)-q_r(78)
    wdot(10) = wdot(10) - qdot
    wdot(11) = wdot(11) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(79)-q_r(79)
    wdot(10) = wdot(10) - qdot
    wdot(11) = wdot(11) + qdot
    wdot(14) = wdot(14) + qdot
    wdot(15) = wdot(15) - qdot

    qdot = q_f(80)-q_r(80)
    wdot(10) = wdot(10) - qdot
    wdot(11) = wdot(11) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(81)-q_r(81)
    wdot(2) = wdot(2) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(82)-q_r(82)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(83)-q_r(83)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(16) = wdot(16) - qdot

    qdot = q_f(84)-q_r(84)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) - qdot

end subroutine

subroutine comp_k_f(tc, invT, k_f)

    implicit none

    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT
    double precision, intent(out) :: k_f(84)

    integer :: i

    do i=1, 84
        k_f(i) = prefactor_units(i)*fwd_A(i)*exp(fwd_beta(i)*tc(1)-activation_units(i)*fwd_Ea(i)*invT)
    end do

end subroutine

subroutine comp_Kc(tc, invT, Kc)

    implicit none

    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT
    double precision, intent(inout) :: Kc(84)

    double precision :: g_RT(21)
    double precision :: refC, refCinv
    integer :: i

    ! compute the Gibbs free energy
    call gibbs(g_RT, tc)

    Kc(1) = g_RT(2) + g_RT(8) - g_RT(10)
    Kc(2) = g_RT(2) + g_RT(10) - g_RT(11)
    Kc(3) = g_RT(2) + g_RT(14) - g_RT(15)
    Kc(4) = g_RT(2) + g_RT(15) - g_RT(16)
    Kc(5) = g_RT(2) + g_RT(17) - g_RT(18)
    Kc(6) = g_RT(2) + g_RT(18) - g_RT(19)
    Kc(7) = g_RT(1) + g_RT(12) - g_RT(15)
    Kc(8) = 2*g_RT(10) - g_RT(19)
    Kc(9) = g_RT(2) + g_RT(3) - g_RT(5)
    Kc(10) = g_RT(3) + g_RT(12) - g_RT(13)
    Kc(11) = g_RT(2) + g_RT(4) - g_RT(7)
    Kc(12) = -g_RT(1) + 2*g_RT(2)
    Kc(13) = g_RT(2) + g_RT(5) - g_RT(6)
    Kc(14) = -g_RT(2) - g_RT(12) + g_RT(14)
    Kc(15) = g_RT(1) - g_RT(2) + g_RT(3) - g_RT(5)
    Kc(16) = g_RT(3) - g_RT(4) - g_RT(5) + g_RT(7)
    Kc(17) = -g_RT(2) + g_RT(3) + g_RT(8) - g_RT(14)
    Kc(18) = -g_RT(2) + g_RT(3) + g_RT(9) - g_RT(14)
    Kc(19) = -g_RT(2) + g_RT(3) + g_RT(10) - g_RT(15)
    Kc(20) = g_RT(3) - g_RT(5) - g_RT(10) + g_RT(11)
    Kc(21) = g_RT(3) - g_RT(5) - g_RT(12) + g_RT(14)
    Kc(22) = -g_RT(2) + g_RT(3) - g_RT(13) + g_RT(14)
    Kc(23) = g_RT(3) - g_RT(5) - g_RT(14) + g_RT(15)
    Kc(24) = g_RT(3) - g_RT(10) - g_RT(14) + g_RT(17)
    Kc(25) = g_RT(3) - g_RT(10) - g_RT(15) + g_RT(18)
    Kc(26) = g_RT(3) - g_RT(5) - g_RT(18) + g_RT(19)
    Kc(27) = -g_RT(3) + g_RT(4) + g_RT(12) - g_RT(13)
    Kc(28) = g_RT(4) - g_RT(7) - g_RT(14) + g_RT(15)
    Kc(29) = g_RT(2) + 2*g_RT(4) - g_RT(4) - g_RT(7)
    Kc(30) = g_RT(2) + g_RT(4) + g_RT(6) - g_RT(6) - g_RT(7)
    Kc(31) = g_RT(2) + g_RT(4) - g_RT(7) + g_RT(20) - g_RT(20)
    Kc(32) = g_RT(2) + g_RT(4) - g_RT(7) + g_RT(21) - g_RT(21)
    Kc(33) = g_RT(2) - g_RT(3) + g_RT(4) - g_RT(5)
    Kc(34) = g_RT(1) - 2*g_RT(1) + 2*g_RT(2)
    Kc(35) = -g_RT(1) + 2*g_RT(2) + g_RT(6) - g_RT(6)
    Kc(36) = -g_RT(1) + 2*g_RT(2) + g_RT(13) - g_RT(13)
    Kc(37) = -g_RT(1) + g_RT(2) - g_RT(4) + g_RT(7)
    Kc(38) = g_RT(2) - 2*g_RT(5) + g_RT(7)
    Kc(39) = -g_RT(1) + g_RT(2) - g_RT(10) + g_RT(11)
    Kc(40) = -g_RT(1) + g_RT(2) - g_RT(12) + g_RT(14)
    Kc(41) = -g_RT(1) + g_RT(2) - g_RT(14) + g_RT(15)
    Kc(42) = g_RT(2) - g_RT(5) - g_RT(10) + g_RT(16)
    Kc(43) = -g_RT(1) + g_RT(2) - g_RT(18) + g_RT(19)
    Kc(44) = g_RT(1) - g_RT(2) + g_RT(5) - g_RT(6)
    Kc(45) = -g_RT(3) + 2*g_RT(5) - g_RT(6)
    Kc(46) = -g_RT(4) + g_RT(5) - g_RT(6) + g_RT(7)
    Kc(47) = -g_RT(2) + g_RT(5) + g_RT(8) - g_RT(15)
    Kc(48) = -g_RT(2) + g_RT(5) + g_RT(9) - g_RT(15)
    Kc(49) = g_RT(5) - g_RT(6) - g_RT(8) + g_RT(10)
    Kc(50) = g_RT(5) - g_RT(6) - g_RT(9) + g_RT(10)
    Kc(51) = g_RT(5) - g_RT(6) - g_RT(10) + g_RT(11)
    Kc(52) = -g_RT(2) + g_RT(5) + g_RT(12) - g_RT(13)
    Kc(53) = g_RT(5) - g_RT(6) - g_RT(12) + g_RT(14)
    Kc(54) = g_RT(5) - g_RT(6) - g_RT(14) + g_RT(15)
    Kc(55) = g_RT(5) - g_RT(6) - g_RT(18) + g_RT(19)
    Kc(56) = -g_RT(5) + g_RT(7) + g_RT(8) - g_RT(15)
    Kc(57) = -g_RT(4) + g_RT(7) + g_RT(10) - g_RT(11)
    Kc(58) = -g_RT(5) + g_RT(7) + g_RT(10) - g_RT(16)
    Kc(59) = -g_RT(5) + g_RT(7) + g_RT(12) - g_RT(13)
    Kc(60) = g_RT(4) - g_RT(5) + g_RT(8) - g_RT(14)
    Kc(61) = g_RT(1) - g_RT(2) + g_RT(8) - g_RT(10)
    Kc(62) = -g_RT(2) + g_RT(8) + g_RT(10) - g_RT(17)
    Kc(63) = g_RT(8) - 2*g_RT(10) + g_RT(11)
    Kc(64) = -g_RT(8) + g_RT(9) + g_RT(20) - g_RT(20)
    Kc(65) = -g_RT(8) + g_RT(9) + g_RT(21) - g_RT(21)
    Kc(66) = -g_RT(2) + g_RT(4) - g_RT(5) + g_RT(9) - g_RT(12)
    Kc(67) = g_RT(4) - g_RT(6) + g_RT(9) - g_RT(12)
    Kc(68) = g_RT(1) - g_RT(2) + g_RT(9) - g_RT(10)
    Kc(69) = g_RT(6) - g_RT(6) - g_RT(8) + g_RT(9)
    Kc(70) = -g_RT(2) + g_RT(9) + g_RT(10) - g_RT(17)
    Kc(71) = g_RT(9) - 2*g_RT(10) + g_RT(11)
    Kc(72) = -g_RT(8) + g_RT(9) + g_RT(12) - g_RT(12)
    Kc(73) = -g_RT(8) + g_RT(9) + g_RT(13) - g_RT(13)
    Kc(74) = g_RT(9) - g_RT(12) + g_RT(13) - g_RT(15)
    Kc(75) = -g_RT(3) + g_RT(4) + g_RT(10) - g_RT(16)
    Kc(76) = g_RT(4) - g_RT(5) + g_RT(10) - g_RT(15)
    Kc(77) = -g_RT(2) + 2*g_RT(10) - g_RT(18)
    Kc(78) = g_RT(10) - g_RT(11) - g_RT(12) + g_RT(14)
    Kc(79) = g_RT(10) - g_RT(11) - g_RT(14) + g_RT(15)
    Kc(80) = g_RT(10) - g_RT(11) - g_RT(18) + g_RT(19)
    Kc(81) = -g_RT(2) + g_RT(6) - g_RT(6) - g_RT(12) + g_RT(14)
    Kc(82) = g_RT(4) - g_RT(7) - g_RT(12) + g_RT(14)
    Kc(83) = g_RT(4) - g_RT(7) - g_RT(15) + g_RT(16)
    Kc(84) = g_RT(4) - g_RT(7) - g_RT(17) + g_RT(18)

    do i=1, 84
        Kc(i) = exp(Kc(i))
    end do

    ! reference concentration: P_atm / (RT) in inverse mol/m^3
    refC = 101325d0 / 8.31451d0 * invT
    refCinv = 1.d0 / refC

    Kc(1) = Kc(1) * (refCinv)
    Kc(2) = Kc(2) * (refCinv)
    Kc(3) = Kc(3) * (refCinv)
    Kc(4) = Kc(4) * (refCinv)
    Kc(5) = Kc(5) * (refCinv)
    Kc(6) = Kc(6) * (refCinv)
    Kc(7) = Kc(7) * (refCinv)
    Kc(8) = Kc(8) * (refCinv)
    Kc(9) = Kc(9) * (refCinv)
    Kc(10) = Kc(10) * (refCinv)
    Kc(11) = Kc(11) * (refCinv)
    Kc(12) = Kc(12) * (refCinv)
    Kc(13) = Kc(13) * (refCinv)
    Kc(14) = Kc(14) * (refC)
    Kc(29) = Kc(29) * (refCinv)
    Kc(30) = Kc(30) * (refCinv)
    Kc(31) = Kc(31) * (refCinv)
    Kc(32) = Kc(32) * (refCinv)
    Kc(34) = Kc(34) * (refCinv)
    Kc(35) = Kc(35) * (refCinv)
    Kc(36) = Kc(36) * (refCinv)
    Kc(66) = Kc(66) * (refC)
    Kc(81) = Kc(81) * (refC)

end subroutine

subroutine comp_qfqr(qf, qr, sc, tc, invT)

    implicit none

    double precision, intent(out) :: qf(84)
    double precision, intent(out) :: qr(84)
    double precision, intent(in) :: sc(21)
    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT

    double precision :: T
    double precision :: mixture
    double precision :: Corr(84)
    double precision :: alpha_troe(8)
    double precision :: redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe
    double precision :: tmp1, tmp2, tmp3
    integer :: i

    ! reaction 1: H + CH2 (+M) <=> CH3 (+M)
    qf(1) = sc(2)*sc(8)
    qr(1) = sc(10)
    ! reaction 2: H + CH3 (+M) <=> CH4 (+M)
    qf(2) = sc(2)*sc(10)
    qr(2) = sc(11)
    ! reaction 3: H + HCO (+M) <=> CH2O (+M)
    qf(3) = sc(2)*sc(14)
    qr(3) = sc(15)
    ! reaction 4: H + CH2O (+M) <=> CH3O (+M)
    qf(4) = sc(2)*sc(15)
    qr(4) = sc(16)
    ! reaction 5: H + C2H4 (+M) <=> C2H5 (+M)
    qf(5) = sc(2)*sc(17)
    qr(5) = sc(18)
    ! reaction 6: H + C2H5 (+M) <=> C2H6 (+M)
    qf(6) = sc(2)*sc(18)
    qr(6) = sc(19)
    ! reaction 7: H2 + CO (+M) <=> CH2O (+M)
    qf(7) = sc(1)*sc(12)
    qr(7) = sc(15)
    ! reaction 8: 2 CH3 (+M) <=> C2H6 (+M)
    qf(8) = sc(10)*sc(10)
    qr(8) = sc(19)
    ! reaction 9: O + H + M <=> OH + M
    qf(9) = sc(2)*sc(3)
    qr(9) = sc(5)
    ! reaction 10: O + CO + M <=> CO2 + M
    qf(10) = sc(3)*sc(12)
    qr(10) = sc(13)
    ! reaction 11: H + O2 + M <=> HO2 + M
    qf(11) = sc(2)*sc(4)
    qr(11) = sc(7)
    ! reaction 12: 2 H + M <=> H2 + M
    qf(12) = sc(2)*sc(2)
    qr(12) = sc(1)
    ! reaction 13: H + OH + M <=> H2O + M
    qf(13) = sc(2)*sc(5)
    qr(13) = sc(6)
    ! reaction 14: HCO + M <=> H + CO + M
    qf(14) = sc(14)
    qr(14) = sc(2)*sc(12)
    ! reaction 15: O + H2 <=> H + OH
    qf(15) = sc(1)*sc(3)
    qr(15) = sc(2)*sc(5)
    ! reaction 16: O + HO2 <=> OH + O2
    qf(16) = sc(3)*sc(7)
    qr(16) = sc(4)*sc(5)
    ! reaction 17: O + CH2 <=> H + HCO
    qf(17) = sc(3)*sc(8)
    qr(17) = sc(2)*sc(14)
    ! reaction 18: O + CH2(S) <=> H + HCO
    qf(18) = sc(3)*sc(9)
    qr(18) = sc(2)*sc(14)
    ! reaction 19: O + CH3 <=> H + CH2O
    qf(19) = sc(3)*sc(10)
    qr(19) = sc(2)*sc(15)
    ! reaction 20: O + CH4 <=> OH + CH3
    qf(20) = sc(3)*sc(11)
    qr(20) = sc(5)*sc(10)
    ! reaction 21: O + HCO <=> OH + CO
    qf(21) = sc(3)*sc(14)
    qr(21) = sc(5)*sc(12)
    ! reaction 22: O + HCO <=> H + CO2
    qf(22) = sc(3)*sc(14)
    qr(22) = sc(2)*sc(13)
    ! reaction 23: O + CH2O <=> OH + HCO
    qf(23) = sc(3)*sc(15)
    qr(23) = sc(5)*sc(14)
    ! reaction 24: O + C2H4 <=> CH3 + HCO
    qf(24) = sc(3)*sc(17)
    qr(24) = sc(10)*sc(14)
    ! reaction 25: O + C2H5 <=> CH3 + CH2O
    qf(25) = sc(3)*sc(18)
    qr(25) = sc(10)*sc(15)
    ! reaction 26: O + C2H6 <=> OH + C2H5
    qf(26) = sc(3)*sc(19)
    qr(26) = sc(5)*sc(18)
    ! reaction 27: O2 + CO <=> O + CO2
    qf(27) = sc(4)*sc(12)
    qr(27) = sc(3)*sc(13)
    ! reaction 28: O2 + CH2O <=> HO2 + HCO
    qf(28) = sc(4)*sc(15)
    qr(28) = sc(7)*sc(14)
    ! reaction 29: H + 2 O2 <=> HO2 + O2
    qf(29) = sc(2)*sc(4)*sc(4)
    qr(29) = sc(4)*sc(7)
    ! reaction 30: H + O2 + H2O <=> HO2 + H2O
    qf(30) = sc(2)*sc(4)*sc(6)
    qr(30) = sc(6)*sc(7)
    ! reaction 31: H + O2 + N2 <=> HO2 + N2
    qf(31) = sc(2)*sc(4)*sc(20)
    qr(31) = sc(7)*sc(20)
    ! reaction 32: H + O2 + AR <=> HO2 + AR
    qf(32) = sc(2)*sc(4)*sc(21)
    qr(32) = sc(7)*sc(21)
    ! reaction 33: H + O2 <=> O + OH
    qf(33) = sc(2)*sc(4)
    qr(33) = sc(3)*sc(5)
    ! reaction 34: 2 H + H2 <=> 2 H2
    qf(34) = sc(1)*sc(2)*sc(2)
    qr(34) = sc(1)*sc(1)
    ! reaction 35: 2 H + H2O <=> H2 + H2O
    qf(35) = sc(2)*sc(2)*sc(6)
    qr(35) = sc(1)*sc(6)
    ! reaction 36: 2 H + CO2 <=> H2 + CO2
    qf(36) = sc(2)*sc(2)*sc(13)
    qr(36) = sc(1)*sc(13)
    ! reaction 37: H + HO2 <=> O2 + H2
    qf(37) = sc(2)*sc(7)
    qr(37) = sc(1)*sc(4)
    ! reaction 38: H + HO2 <=> 2 OH
    qf(38) = sc(2)*sc(7)
    qr(38) = sc(5)*sc(5)
    ! reaction 39: H + CH4 <=> CH3 + H2
    qf(39) = sc(2)*sc(11)
    qr(39) = sc(1)*sc(10)
    ! reaction 40: H + HCO <=> H2 + CO
    qf(40) = sc(2)*sc(14)
    qr(40) = sc(1)*sc(12)
    ! reaction 41: H + CH2O <=> HCO + H2
    qf(41) = sc(2)*sc(15)
    qr(41) = sc(1)*sc(14)
    ! reaction 42: H + CH3O <=> OH + CH3
    qf(42) = sc(2)*sc(16)
    qr(42) = sc(5)*sc(10)
    ! reaction 43: H + C2H6 <=> C2H5 + H2
    qf(43) = sc(2)*sc(19)
    qr(43) = sc(1)*sc(18)
    ! reaction 44: OH + H2 <=> H + H2O
    qf(44) = sc(1)*sc(5)
    qr(44) = sc(2)*sc(6)
    ! reaction 45: 2 OH <=> O + H2O
    qf(45) = sc(5)*sc(5)
    qr(45) = sc(3)*sc(6)
    ! reaction 46: OH + HO2 <=> O2 + H2O
    qf(46) = sc(5)*sc(7)
    qr(46) = sc(4)*sc(6)
    ! reaction 47: OH + CH2 <=> H + CH2O
    qf(47) = sc(5)*sc(8)
    qr(47) = sc(2)*sc(15)
    ! reaction 48: OH + CH2(S) <=> H + CH2O
    qf(48) = sc(5)*sc(9)
    qr(48) = sc(2)*sc(15)
    ! reaction 49: OH + CH3 <=> CH2 + H2O
    qf(49) = sc(5)*sc(10)
    qr(49) = sc(6)*sc(8)
    ! reaction 50: OH + CH3 <=> CH2(S) + H2O
    qf(50) = sc(5)*sc(10)
    qr(50) = sc(6)*sc(9)
    ! reaction 51: OH + CH4 <=> CH3 + H2O
    qf(51) = sc(5)*sc(11)
    qr(51) = sc(6)*sc(10)
    ! reaction 52: OH + CO <=> H + CO2
    qf(52) = sc(5)*sc(12)
    qr(52) = sc(2)*sc(13)
    ! reaction 53: OH + HCO <=> H2O + CO
    qf(53) = sc(5)*sc(14)
    qr(53) = sc(6)*sc(12)
    ! reaction 54: OH + CH2O <=> HCO + H2O
    qf(54) = sc(5)*sc(15)
    qr(54) = sc(6)*sc(14)
    ! reaction 55: OH + C2H6 <=> C2H5 + H2O
    qf(55) = sc(5)*sc(19)
    qr(55) = sc(6)*sc(18)
    ! reaction 56: HO2 + CH2 <=> OH + CH2O
    qf(56) = sc(7)*sc(8)
    qr(56) = sc(5)*sc(15)
    ! reaction 57: HO2 + CH3 <=> O2 + CH4
    qf(57) = sc(7)*sc(10)
    qr(57) = sc(4)*sc(11)
    ! reaction 58: HO2 + CH3 <=> OH + CH3O
    qf(58) = sc(7)*sc(10)
    qr(58) = sc(5)*sc(16)
    ! reaction 59: HO2 + CO <=> OH + CO2
    qf(59) = sc(7)*sc(12)
    qr(59) = sc(5)*sc(13)
    ! reaction 60: CH2 + O2 <=> OH + HCO
    qf(60) = sc(4)*sc(8)
    qr(60) = sc(5)*sc(14)
    ! reaction 61: CH2 + H2 <=> H + CH3
    qf(61) = sc(1)*sc(8)
    qr(61) = sc(2)*sc(10)
    ! reaction 62: CH2 + CH3 <=> H + C2H4
    qf(62) = sc(8)*sc(10)
    qr(62) = sc(2)*sc(17)
    ! reaction 63: CH2 + CH4 <=> 2 CH3
    qf(63) = sc(8)*sc(11)
    qr(63) = sc(10)*sc(10)
    ! reaction 64: CH2(S) + N2 <=> CH2 + N2
    qf(64) = sc(9)*sc(20)
    qr(64) = sc(8)*sc(20)
    ! reaction 65: CH2(S) + AR <=> CH2 + AR
    qf(65) = sc(9)*sc(21)
    qr(65) = sc(8)*sc(21)
    ! reaction 66: CH2(S) + O2 <=> H + OH + CO
    qf(66) = sc(4)*sc(9)
    qr(66) = sc(2)*sc(5)*sc(12)
    ! reaction 67: CH2(S) + O2 <=> CO + H2O
    qf(67) = sc(4)*sc(9)
    qr(67) = sc(6)*sc(12)
    ! reaction 68: CH2(S) + H2 <=> CH3 + H
    qf(68) = sc(1)*sc(9)
    qr(68) = sc(2)*sc(10)
    ! reaction 69: CH2(S) + H2O <=> CH2 + H2O
    qf(69) = sc(6)*sc(9)
    qr(69) = sc(6)*sc(8)
    ! reaction 70: CH2(S) + CH3 <=> H + C2H4
    qf(70) = sc(9)*sc(10)
    qr(70) = sc(2)*sc(17)
    ! reaction 71: CH2(S) + CH4 <=> 2 CH3
    qf(71) = sc(9)*sc(11)
    qr(71) = sc(10)*sc(10)
    ! reaction 72: CH2(S) + CO <=> CH2 + CO
    qf(72) = sc(9)*sc(12)
    qr(72) = sc(8)*sc(12)
    ! reaction 73: CH2(S) + CO2 <=> CH2 + CO2
    qf(73) = sc(9)*sc(13)
    qr(73) = sc(8)*sc(13)
    ! reaction 74: CH2(S) + CO2 <=> CO + CH2O
    qf(74) = sc(9)*sc(13)
    qr(74) = sc(12)*sc(15)
    ! reaction 75: CH3 + O2 <=> O + CH3O
    qf(75) = sc(4)*sc(10)
    qr(75) = sc(3)*sc(16)
    ! reaction 76: CH3 + O2 <=> OH + CH2O
    qf(76) = sc(4)*sc(10)
    qr(76) = sc(5)*sc(15)
    ! reaction 77: 2 CH3 <=> H + C2H5
    qf(77) = sc(10)*sc(10)
    qr(77) = sc(2)*sc(18)
    ! reaction 78: CH3 + HCO <=> CH4 + CO
    qf(78) = sc(10)*sc(14)
    qr(78) = sc(11)*sc(12)
    ! reaction 79: CH3 + CH2O <=> HCO + CH4
    qf(79) = sc(10)*sc(15)
    qr(79) = sc(11)*sc(14)
    ! reaction 80: CH3 + C2H6 <=> C2H5 + CH4
    qf(80) = sc(10)*sc(19)
    qr(80) = sc(11)*sc(18)
    ! reaction 81: HCO + H2O <=> H + CO + H2O
    qf(81) = sc(6)*sc(14)
    qr(81) = sc(2)*sc(6)*sc(12)
    ! reaction 82: HCO + O2 <=> HO2 + CO
    qf(82) = sc(4)*sc(14)
    qr(82) = sc(7)*sc(12)
    ! reaction 83: CH3O + O2 <=> HO2 + CH2O
    qf(83) = sc(4)*sc(16)
    qr(83) = sc(7)*sc(15)
    ! reaction 84: C2H5 + O2 <=> HO2 + C2H4
    qf(84) = sc(4)*sc(18)
    qr(84) = sc(7)*sc(17)

    T = tc(2)

    ! compute the mixture concentration
    mixture = 0.d0
    do i=1, 21
        mixture = mixture + sc(i)
    end do

    do i=1, 84
        Corr(i) = 1.d0
    end do

    ! troe
    alpha_troe(1) = mixture + (TB(1) % vector(1) - 1)*sc(1) + (TB(1) % vector(2) - 1)*sc(6) + (TB(1) % vector(3) - 1)*sc(11) + (TB(1) % vector(4) - 1)*sc(12) + (TB(1) % vector(5) - 1)*sc(13) + (TB(1) % vector(6) - 1)*sc(19) + (TB(1) % vector(7) - 1)*sc(21)
    alpha_troe(2) = mixture + (TB(2) % vector(1) - 1)*sc(1) + (TB(2) % vector(2) - 1)*sc(6) + (TB(2) % vector(3) - 1)*sc(11) + (TB(2) % vector(4) - 1)*sc(12) + (TB(2) % vector(5) - 1)*sc(13) + (TB(2) % vector(6) - 1)*sc(19) + (TB(2) % vector(7) - 1)*sc(21)
    alpha_troe(3) = mixture + (TB(3) % vector(1) - 1)*sc(1) + (TB(3) % vector(2) - 1)*sc(6) + (TB(3) % vector(3) - 1)*sc(11) + (TB(3) % vector(4) - 1)*sc(12) + (TB(3) % vector(5) - 1)*sc(13) + (TB(3) % vector(6) - 1)*sc(19) + (TB(3) % vector(7) - 1)*sc(21)
    alpha_troe(4) = mixture + (TB(4) % vector(1) - 1)*sc(1) + (TB(4) % vector(2) - 1)*sc(6) + (TB(4) % vector(3) - 1)*sc(11) + (TB(4) % vector(4) - 1)*sc(12) + (TB(4) % vector(5) - 1)*sc(13) + (TB(4) % vector(6) - 1)*sc(19)
    alpha_troe(5) = mixture + (TB(5) % vector(1) - 1)*sc(1) + (TB(5) % vector(2) - 1)*sc(6) + (TB(5) % vector(3) - 1)*sc(11) + (TB(5) % vector(4) - 1)*sc(12) + (TB(5) % vector(5) - 1)*sc(13) + (TB(5) % vector(6) - 1)*sc(19) + (TB(5) % vector(7) - 1)*sc(21)
    alpha_troe(6) = mixture + (TB(6) % vector(1) - 1)*sc(1) + (TB(6) % vector(2) - 1)*sc(6) + (TB(6) % vector(3) - 1)*sc(11) + (TB(6) % vector(4) - 1)*sc(12) + (TB(6) % vector(5) - 1)*sc(13) + (TB(6) % vector(6) - 1)*sc(19) + (TB(6) % vector(7) - 1)*sc(21)
    alpha_troe(7) = mixture + (TB(7) % vector(1) - 1)*sc(1) + (TB(7) % vector(2) - 1)*sc(6) + (TB(7) % vector(3) - 1)*sc(11) + (TB(7) % vector(4) - 1)*sc(12) + (TB(7) % vector(5) - 1)*sc(13) + (TB(7) % vector(6) - 1)*sc(19) + (TB(7) % vector(7) - 1)*sc(21)
    alpha_troe(8) = mixture + (TB(8) % vector(1) - 1)*sc(1) + (TB(8) % vector(2) - 1)*sc(6) + (TB(8) % vector(3) - 1)*sc(11) + (TB(8) % vector(4) - 1)*sc(12) + (TB(8) % vector(5) - 1)*sc(13) + (TB(8) % vector(6) - 1)*sc(19) + (TB(8) % vector(7) - 1)*sc(21)

    do i=1, 8
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
    Corr(9) = mixture + (TB(9) % vector(1) - 1)*sc(1) + (TB(9) % vector(2) - 1)*sc(6) + (TB(9) % vector(3) - 1)*sc(11) + (TB(9) % vector(4) - 1)*sc(12) + (TB(9) % vector(5) - 1)*sc(13) + (TB(9) % vector(6) - 1)*sc(19) + (TB(9) % vector(7) - 1)*sc(21)
    Corr(10) = mixture + (TB(10) % vector(1) - 1)*sc(1) + (TB(10) % vector(2) - 1)*sc(4) + (TB(10) % vector(3) - 1)*sc(6) + (TB(10) % vector(4) - 1)*sc(11) + (TB(10) % vector(5) - 1)*sc(12) + (TB(10) % vector(6) - 1)*sc(13) + (TB(10) % vector(7) - 1)*sc(19) + (TB(10) % vector(8) - 1)*sc(21)
    Corr(11) = mixture + (TB(11) % vector(1) - 1)*sc(4) + (TB(11) % vector(2) - 1)*sc(6) + (TB(11) % vector(3) - 1)*sc(12) + (TB(11) % vector(4) - 1)*sc(13) + (TB(11) % vector(5) - 1)*sc(19) + (TB(11) % vector(6) - 1)*sc(20) + (TB(11) % vector(7) - 1)*sc(21)
    Corr(12) = mixture + (TB(12) % vector(1) - 1)*sc(1) + (TB(12) % vector(2) - 1)*sc(6) + (TB(12) % vector(3) - 1)*sc(11) + (TB(12) % vector(4) - 1)*sc(13) + (TB(12) % vector(5) - 1)*sc(19) + (TB(12) % vector(6) - 1)*sc(21)
    Corr(13) = mixture + (TB(13) % vector(1) - 1)*sc(1) + (TB(13) % vector(2) - 1)*sc(6) + (TB(13) % vector(3) - 1)*sc(11) + (TB(13) % vector(4) - 1)*sc(19) + (TB(13) % vector(5) - 1)*sc(21)
    Corr(14) = mixture + (TB(14) % vector(1) - 1)*sc(1) + (TB(14) % vector(2) - 1)*sc(6) + (TB(14) % vector(3) - 1)*sc(11) + (TB(14) % vector(4) - 1)*sc(12) + (TB(14) % vector(5) - 1)*sc(13) + (TB(14) % vector(6) - 1)*sc(19)

    do i=1, 84
        qf(i) = qf(i) * (Corr(i) * k_f_save(i))
        qr(i) = qr(i) * (Corr(i) * k_f_save(i) / Kc_save(i))
    end do

end subroutine

! compute the g/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine gibbs(species, tc)

    implicit none

    double precision, intent(out) :: species(21)
    double precision, intent(in) :: tc(5)

    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            -9.179351730000000d+02 * invT &
            +1.661320882000000d+00 &
            -2.344331120000000d+00 * tc(1) &
            -3.990260375000000d-03 * tc(2) &
            +3.246358500000000d-06 * tc(3) &
            -1.679767450000000d-09 * tc(4) &
            +3.688058805000000d-13 * tc(5)
        ! species 2: H
        species(2) = &
            +2.547365990000000d+04 * invT &
            +2.946682853000000d+00 &
            -2.500000000000000d+00 * tc(1) &
            -3.526664095000000d-13 * tc(2) &
            +3.326532733333333d-16 * tc(3) &
            -1.917346933333333d-19 * tc(4) &
            +4.638661660000000d-23 * tc(5)
        ! species 3: O
        species(3) = &
            +2.912225920000000d+04 * invT &
            +1.116333640000000d+00 &
            -3.168267100000000d+00 * tc(1) &
            +1.639659420000000d-03 * tc(2) &
            -1.107177326666667d-06 * tc(3) &
            +5.106721866666666d-10 * tc(4) &
            -1.056329855000000d-13 * tc(5)
        ! species 4: O2
        species(4) = &
            -1.063943560000000d+03 * invT &
            +1.247806300000001d-01 &
            -3.782456360000000d+00 * tc(1) &
            +1.498367080000000d-03 * tc(2) &
            -1.641217001666667d-06 * tc(3) &
            +8.067745908333334d-10 * tc(4) &
            -1.621864185000000d-13 * tc(5)
        ! species 5: OH
        species(5) = &
            +3.615080560000000d+03 * invT &
            +4.095940888000000d+00 &
            -3.992015430000000d+00 * tc(1) &
            +1.200658760000000d-03 * tc(2) &
            -7.696564016666666d-07 * tc(3) &
            +3.234277775000000d-10 * tc(4) &
            -6.820573500000000d-14 * tc(5)
        ! species 6: H2O
        species(6) = &
            -3.029372670000000d+04 * invT &
            +5.047672768000000d+00 &
            -4.198640560000000d+00 * tc(1) &
            +1.018217050000000d-03 * tc(2) &
            -1.086733685000000d-06 * tc(3) &
            +4.573308850000000d-10 * tc(4) &
            -8.859890850000000d-14 * tc(5)
        ! species 7: HO2
        species(7) = &
            +2.948080400000000d+02 * invT &
            +5.851355599999999d-01 &
            -4.301798010000000d+00 * tc(1) &
            +2.374560255000000d-03 * tc(2) &
            -3.526381516666666d-06 * tc(3) &
            +2.023032450000000d-09 * tc(4) &
            -4.646125620000001d-13 * tc(5)
        ! species 8: CH2
        species(8) = &
            +4.600404010000000d+04 * invT &
            +2.200146820000000d+00 &
            -3.762678670000000d+00 * tc(1) &
            -4.844360715000000d-04 * tc(2) &
            -4.658164016666667d-07 * tc(3) &
            +3.209092941666667d-10 * tc(4) &
            -8.437085950000000d-14 * tc(5)
        ! species 9: CH2(S)
        species(9) = &
            +5.049681630000000d+04 * invT &
            +4.967723077000000d+00 &
            -4.198604110000000d+00 * tc(1) &
            +1.183307095000000d-03 * tc(2) &
            -1.372160366666667d-06 * tc(3) &
            +5.573466508333334d-10 * tc(4) &
            -9.715736850000000d-14 * tc(5)
        ! species 10: CH3
        species(10) = &
            +1.644499880000000d+04 * invT &
            +2.069026070000000d+00 &
            -3.673590400000000d+00 * tc(1) &
            -1.005475875000000d-03 * tc(2) &
            -9.550364266666668d-07 * tc(3) &
            +5.725978541666666d-10 * tc(4) &
            -1.271928670000000d-13 * tc(5)
        ! species 11: CH4
        species(11) = &
            -1.024664760000000d+04 * invT &
            +9.791179889999999d+00 &
            -5.149876130000000d+00 * tc(1) &
            +6.835489400000000d-03 * tc(2) &
            -8.196676650000000d-06 * tc(3) &
            +4.039525216666667d-09 * tc(4) &
            -8.334697800000000d-13 * tc(5)
        ! species 12: CO
        species(12) = &
            -1.434408600000000d+04 * invT &
            +7.112418999999992d-02 &
            -3.579533470000000d+00 * tc(1) &
            +3.051768400000000d-04 * tc(2) &
            -1.694690550000000d-07 * tc(3) &
            -7.558382366666667d-11 * tc(4) &
            +4.522122495000000d-14 * tc(5)
        ! species 13: CO2
        species(13) = &
            -4.837196970000000d+04 * invT &
            -7.544278700000000d+00 &
            -2.356773520000000d+00 * tc(1) &
            -4.492298385000000d-03 * tc(2) &
            +1.187260448333333d-06 * tc(3) &
            -2.049325183333333d-10 * tc(4) &
            +7.184977399999999d-15 * tc(5)
        ! species 14: HCO
        species(14) = &
            +3.839564960000000d+03 * invT &
            +8.268134100000002d-01 &
            -4.221185840000000d+00 * tc(1) &
            +1.621962660000000d-03 * tc(2) &
            -2.296657433333333d-06 * tc(3) &
            +1.109534108333333d-09 * tc(4) &
            -2.168844325000000d-13 * tc(5)
        ! species 15: CH2O
        species(15) = &
            -1.430895670000000d+04 * invT &
            +4.190910250000000d+00 &
            -4.793723150000000d+00 * tc(1) &
            +4.954166845000000d-03 * tc(2) &
            -6.220333466666666d-06 * tc(3) &
            +3.160710508333333d-09 * tc(4) &
            -6.588632600000000d-13 * tc(5)
        ! species 16: CH3O
        species(16) = &
            +9.786011000000000d+02 * invT &
            -1.104597300000000d+01 &
            -2.106204000000000d+00 * tc(1) &
            -3.608297500000000d-03 * tc(2) &
            -8.897453333333333d-07 * tc(3) &
            +6.148030000000000d-10 * tc(4) &
            -1.037805000000000d-13 * tc(5)
        ! species 17: C2H4
        species(17) = &
            +5.089775930000000d+03 * invT &
            -1.381294799999999d-01 &
            -3.959201480000000d+00 * tc(1) &
            +3.785261235000000d-03 * tc(2) &
            -9.516504866666667d-06 * tc(3) &
            +5.763239608333333d-09 * tc(4) &
            -1.349421865000000d-12 * tc(5)
        ! species 18: C2H5
        species(18) = &
            +1.284162650000000d+04 * invT &
            -4.007435600000004d-01 &
            -4.306465680000000d+00 * tc(1) &
            +2.093294460000000d-03 * tc(2) &
            -8.285713450000000d-06 * tc(3) &
            +4.992721716666666d-09 * tc(4) &
            -1.152545020000000d-12 * tc(5)
        ! species 19: C2H6
        species(19) = &
            -1.152220550000000d+04 * invT &
            +1.624601760000000d+00 &
            -4.291424920000000d+00 * tc(1) &
            +2.750771350000000d-03 * tc(2) &
            -9.990638133333334d-06 * tc(3) &
            +5.903885708333334d-09 * tc(4) &
            -1.343428855000000d-12 * tc(5)
        ! species 20: N2
        species(20) = &
            -1.020899900000000d+03 * invT &
            -6.516950000000001d-01 &
            -3.298677000000000d+00 * tc(1) &
            -7.041202000000000d-04 * tc(2) &
            +6.605369999999999d-07 * tc(3) &
            -4.701262500000001d-10 * tc(4) &
            +1.222427000000000d-13 * tc(5)
        ! species 21: AR
        species(21) = &
            -7.453750000000000d+02 * invT &
            -1.866000000000000d+00 &
            -2.500000000000000d+00 * tc(1) &
            -0.000000000000000d+00 * tc(2) &
            -0.000000000000000d+00 * tc(3) &
            -0.000000000000000d+00 * tc(4) &
            -0.000000000000000d+00 * tc(5)
    else
        !species 1: H2
        species(1) = &
            -9.501589220000000d+02 * invT &
            +6.542302510000000d+00 &
            -3.337279200000000d+00 * tc(1) &
            +2.470123655000000d-05 * tc(2) &
            -8.324279633333333d-08 * tc(3) &
            +1.496386616666667d-11 * tc(4) &
            -1.001276880000000d-15 * tc(5)
        !species 2: H
        species(2) = &
            +2.547365990000000d+04 * invT &
            +2.946682924000000d+00 &
            -2.500000010000000d+00 * tc(1) &
            +1.154214865000000d-11 * tc(2) &
            -2.692699133333334d-15 * tc(3) &
            +3.945960291666667d-19 * tc(4) &
            -2.490986785000000d-23 * tc(5)
        !species 3: O
        species(3) = &
            +2.921757910000000d+04 * invT &
            -2.214917859999999d+00 &
            -2.569420780000000d+00 * tc(1) &
            +4.298705685000000d-05 * tc(2) &
            -6.991409816666667d-09 * tc(3) &
            +8.348149916666666d-13 * tc(4) &
            -6.141684549999999d-17 * tc(5)
        !species 4: O2
        species(4) = &
            -1.088457720000000d+03 * invT &
            -2.170693450000000d+00 &
            -3.282537840000000d+00 * tc(1) &
            -7.415437700000000d-04 * tc(2) &
            +1.263277781666667d-07 * tc(3) &
            -1.745587958333333d-11 * tc(4) &
            +1.083588970000000d-15 * tc(5)
        !species 5: OH
        species(5) = &
            +3.858657000000000d+03 * invT &
            -1.383808430000000d+00 &
            -3.092887670000000d+00 * tc(1) &
            -2.742148580000000d-04 * tc(2) &
            -2.108420466666667d-08 * tc(3) &
            +7.328846300000000d-12 * tc(4) &
            -5.870618800000000d-16 * tc(5)
        !species 6: H2O
        species(6) = &
            -3.000429710000000d+04 * invT &
            -1.932777610000000d+00 &
            -3.033992490000000d+00 * tc(1) &
            -1.088459020000000d-03 * tc(2) &
            +2.734541966666666d-08 * tc(3) &
            +8.086832250000000d-12 * tc(4) &
            -8.410049600000000d-16 * tc(5)
        !species 7: HO2
        species(7) = &
            +1.118567130000000d+02 * invT &
            +2.321087500000001d-01 &
            -4.017210900000000d+00 * tc(1) &
            -1.119910065000000d-03 * tc(2) &
            +1.056096916666667d-07 * tc(3) &
            -9.520530833333334d-12 * tc(4) &
            +5.395426750000000d-16 * tc(5)
        !species 8: CH2
        species(8) = &
            +4.626360400000000d+04 * invT &
            -3.297092110000000d+00 &
            -2.874101130000000d+00 * tc(1) &
            -1.828196460000000d-03 * tc(2) &
            +2.348243283333333d-07 * tc(3) &
            -2.168162908333333d-11 * tc(4) &
            +9.386378350000000d-16 * tc(5)
        !species 9: CH2(S)
        species(9) = &
            +5.092599970000000d+04 * invT &
            -6.334463270000000d+00 &
            -2.292038420000000d+00 * tc(1) &
            -2.327943185000000d-03 * tc(2) &
            +3.353199116666667d-07 * tc(3) &
            -3.482550000000000d-11 * tc(4) &
            +1.698581825000000d-15 * tc(5)
        !species 10: CH3
        species(10) = &
            +1.677558430000000d+04 * invT &
            -6.194354070000000d+00 &
            -2.285717720000000d+00 * tc(1) &
            -3.619950185000000d-03 * tc(2) &
            +4.978572466666667d-07 * tc(3) &
            -4.964038700000000d-11 * tc(4) &
            +2.335771970000000d-15 * tc(5)
        !species 11: CH4
        species(11) = &
            -9.468344590000001d+03 * invT &
            -1.836246650500000d+01 &
            -7.485149500000000d-02 * tc(1) &
            -6.695473350000000d-03 * tc(2) &
            +9.554763483333333d-07 * tc(3) &
            -1.019104458333333d-10 * tc(4) &
            +5.090761500000000d-15 * tc(5)
        !species 12: CO
        species(12) = &
            -1.415187240000000d+04 * invT &
            -5.103502110000000d+00 &
            -2.715185610000000d+00 * tc(1) &
            -1.031263715000000d-03 * tc(2) &
            +1.664709618333334d-07 * tc(3) &
            -1.917108400000000d-11 * tc(4) &
            +1.018238580000000d-15 * tc(5)
        !species 13: CO2
        species(13) = &
            -4.875916600000000d+04 * invT &
            +1.585822230000000d+00 &
            -3.857460290000000d+00 * tc(1) &
            -2.207185130000000d-03 * tc(2) &
            +3.691356733333334d-07 * tc(3) &
            -4.362418233333334d-11 * tc(4) &
            +2.360420820000000d-15 * tc(5)
        !species 14: HCO
        species(14) = &
            +4.011918150000000d+03 * invT &
            -7.026170540000000d+00 &
            -2.772174380000000d+00 * tc(1) &
            -2.478477630000000d-03 * tc(2) &
            +4.140760216666667d-07 * tc(3) &
            -4.909681483333334d-11 * tc(4) &
            +2.667543555000000d-15 * tc(5)
        !species 15: CH2O
        species(15) = &
            -1.399583230000000d+04 * invT &
            -1.189563292000000d+01 &
            -1.760690080000000d+00 * tc(1) &
            -4.600000410000000d-03 * tc(2) &
            +7.370980216666666d-07 * tc(3) &
            -8.386767666666666d-11 * tc(4) &
            +4.419278200000001d-15 * tc(5)
        !species 16: CH3O
        species(16) = &
            +1.278325200000000d+02 * invT &
            +8.412240000000000d-01 &
            -3.770799000000000d+00 * tc(1) &
            -3.935748500000000d-03 * tc(2) &
            +4.427306666666667d-07 * tc(3) &
            -3.287025833333333d-11 * tc(4) &
            +1.056308000000000d-15 * tc(5)
        !species 17: C2H4
        species(17) = &
            +4.939886140000000d+03 * invT &
            -8.269258140000002d+00 &
            -2.036111160000000d+00 * tc(1) &
            -7.322707550000000d-03 * tc(2) &
            +1.118463191666667d-06 * tc(3) &
            -1.226857691666667d-10 * tc(4) &
            +6.285303050000000d-15 * tc(5)
        !species 18: C2H5
        species(18) = &
            +1.285752000000000d+04 * invT &
            -1.150777788000000d+01 &
            -1.954656420000000d+00 * tc(1) &
            -8.698636100000001d-03 * tc(2) &
            +1.330344446666667d-06 * tc(3) &
            -1.460147408333333d-10 * tc(4) &
            +7.482078800000000d-15 * tc(5)
        !species 19: C2H6
        species(19) = &
            -1.142639320000000d+04 * invT &
            -1.404372920000000d+01 &
            -1.071881500000000d+00 * tc(1) &
            -1.084263385000000d-02 * tc(2) &
            +1.670934450000000d-06 * tc(3) &
            -1.845100008333333d-10 * tc(4) &
            +9.500144500000000d-15 * tc(5)
        !species 20: N2
        species(20) = &
            -9.227977000000000d+02 * invT &
            -3.053888000000000d+00 &
            -2.926640000000000d+00 * tc(1) &
            -7.439884000000000d-04 * tc(2) &
            +9.474600000000001d-08 * tc(3) &
            -8.414198333333333d-12 * tc(4) &
            +3.376675500000000d-16 * tc(5)
        !species 21: AR
        species(21) = &
            -7.453750000000000d+02 * invT &
            -1.866000000000000d+00 &
            -2.500000000000000d+00 * tc(1) &
            -0.000000000000000d+00 * tc(2) &
            -0.000000000000000d+00 * tc(3) &
            -0.000000000000000d+00 * tc(4) &
            -0.000000000000000d+00 * tc(5)
    end if

end subroutine


! compute Cv/R at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine cv_R(species, tc)

    implicit none

    double precision, intent(out) :: species(21)
    double precision, intent(in) :: tc(5)

    double precision :: T

    T = tc(2)

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +1.34433112d+00 &
            +7.98052075d-03 * tc(2) &
            -1.94781510d-05 * tc(3) &
            +2.01572094d-08 * tc(4) &
            -7.37611761d-12 * tc(5)
        ! species 2: H
        species(2) = &
            +1.50000000d+00 &
            +7.05332819d-13 * tc(2) &
            -1.99591964d-15 * tc(3) &
            +2.30081632d-18 * tc(4) &
            -9.27732332d-22 * tc(5)
        ! species 3: O
        species(3) = &
            +2.16826710d+00 &
            -3.27931884d-03 * tc(2) &
            +6.64306396d-06 * tc(3) &
            -6.12806624d-09 * tc(4) &
            +2.11265971d-12 * tc(5)
        ! species 4: O2
        species(4) = &
            +2.78245636d+00 &
            -2.99673416d-03 * tc(2) &
            +9.84730201d-06 * tc(3) &
            -9.68129509d-09 * tc(4) &
            +3.24372837d-12 * tc(5)
        ! species 5: OH
        species(5) = &
            +2.99201543d+00 &
            -2.40131752d-03 * tc(2) &
            +4.61793841d-06 * tc(3) &
            -3.88113333d-09 * tc(4) &
            +1.36411470d-12 * tc(5)
        ! species 6: H2O
        species(6) = &
            +3.19864056d+00 &
            -2.03643410d-03 * tc(2) &
            +6.52040211d-06 * tc(3) &
            -5.48797062d-09 * tc(4) &
            +1.77197817d-12 * tc(5)
        ! species 7: HO2
        species(7) = &
            +3.30179801d+00 &
            -4.74912051d-03 * tc(2) &
            +2.11582891d-05 * tc(3) &
            -2.42763894d-08 * tc(4) &
            +9.29225124d-12 * tc(5)
        ! species 8: CH2
        species(8) = &
            +2.76267867d+00 &
            +9.68872143d-04 * tc(2) &
            +2.79489841d-06 * tc(3) &
            -3.85091153d-09 * tc(4) &
            +1.68741719d-12 * tc(5)
        ! species 9: CH2(S)
        species(9) = &
            +3.19860411d+00 &
            -2.36661419d-03 * tc(2) &
            +8.23296220d-06 * tc(3) &
            -6.68815981d-09 * tc(4) &
            +1.94314737d-12 * tc(5)
        ! species 10: CH3
        species(10) = &
            +2.67359040d+00 &
            +2.01095175d-03 * tc(2) &
            +5.73021856d-06 * tc(3) &
            -6.87117425d-09 * tc(4) &
            +2.54385734d-12 * tc(5)
        ! species 11: CH4
        species(11) = &
            +4.14987613d+00 &
            -1.36709788d-02 * tc(2) &
            +4.91800599d-05 * tc(3) &
            -4.84743026d-08 * tc(4) &
            +1.66693956d-11 * tc(5)
        ! species 12: CO
        species(12) = &
            +2.57953347d+00 &
            -6.10353680d-04 * tc(2) &
            +1.01681433d-06 * tc(3) &
            +9.07005884d-10 * tc(4) &
            -9.04424499d-13 * tc(5)
        ! species 13: CO2
        species(13) = &
            +1.35677352d+00 &
            +8.98459677d-03 * tc(2) &
            -7.12356269d-06 * tc(3) &
            +2.45919022d-09 * tc(4) &
            -1.43699548d-13 * tc(5)
        ! species 14: HCO
        species(14) = &
            +3.22118584d+00 &
            -3.24392532d-03 * tc(2) &
            +1.37799446d-05 * tc(3) &
            -1.33144093d-08 * tc(4) &
            +4.33768865d-12 * tc(5)
        ! species 15: CH2O
        species(15) = &
            +3.79372315d+00 &
            -9.90833369d-03 * tc(2) &
            +3.73220008d-05 * tc(3) &
            -3.79285261d-08 * tc(4) &
            +1.31772652d-11 * tc(5)
        ! species 16: CH3O
        species(16) = &
            +1.10620400d+00 &
            +7.21659500d-03 * tc(2) &
            +5.33847200d-06 * tc(3) &
            -7.37763600d-09 * tc(4) &
            +2.07561000d-12 * tc(5)
        ! species 17: C2H4
        species(17) = &
            +2.95920148d+00 &
            -7.57052247d-03 * tc(2) &
            +5.70990292d-05 * tc(3) &
            -6.91588753d-08 * tc(4) &
            +2.69884373d-11 * tc(5)
        ! species 18: C2H5
        species(18) = &
            +3.30646568d+00 &
            -4.18658892d-03 * tc(2) &
            +4.97142807d-05 * tc(3) &
            -5.99126606d-08 * tc(4) &
            +2.30509004d-11 * tc(5)
        ! species 19: C2H6
        species(19) = &
            +3.29142492d+00 &
            -5.50154270d-03 * tc(2) &
            +5.99438288d-05 * tc(3) &
            -7.08466285d-08 * tc(4) &
            +2.68685771d-11 * tc(5)
        ! species 20: N2
        species(20) = &
            +2.29867700d+00 &
            +1.40824040d-03 * tc(2) &
            -3.96322200d-06 * tc(3) &
            +5.64151500d-09 * tc(4) &
            -2.44485400d-12 * tc(5)
        ! species 21: AR
        species(21) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
    else
        !species 1: H2
        species(1) = &
            +2.33727920d+00 &
            -4.94024731d-05 * tc(2) &
            +4.99456778d-07 * tc(3) &
            -1.79566394d-10 * tc(4) &
            +2.00255376d-14 * tc(5)
        !species 2: H
        species(2) = &
            +1.50000001d+00 &
            -2.30842973d-11 * tc(2) &
            +1.61561948d-14 * tc(3) &
            -4.73515235d-18 * tc(4) &
            +4.98197357d-22 * tc(5)
        !species 3: O
        species(3) = &
            +1.56942078d+00 &
            -8.59741137d-05 * tc(2) &
            +4.19484589d-08 * tc(3) &
            -1.00177799d-11 * tc(4) &
            +1.22833691d-15 * tc(5)
        !species 4: O2
        species(4) = &
            +2.28253784d+00 &
            +1.48308754d-03 * tc(2) &
            -7.57966669d-07 * tc(3) &
            +2.09470555d-10 * tc(4) &
            -2.16717794d-14 * tc(5)
        !species 5: OH
        species(5) = &
            +2.09288767d+00 &
            +5.48429716d-04 * tc(2) &
            +1.26505228d-07 * tc(3) &
            -8.79461556d-11 * tc(4) &
            +1.17412376d-14 * tc(5)
        !species 6: H2O
        species(6) = &
            +2.03399249d+00 &
            +2.17691804d-03 * tc(2) &
            -1.64072518d-07 * tc(3) &
            -9.70419870d-11 * tc(4) &
            +1.68200992d-14 * tc(5)
        !species 7: HO2
        species(7) = &
            +3.01721090d+00 &
            +2.23982013d-03 * tc(2) &
            -6.33658150d-07 * tc(3) &
            +1.14246370d-10 * tc(4) &
            -1.07908535d-14 * tc(5)
        !species 8: CH2
        species(8) = &
            +1.87410113d+00 &
            +3.65639292d-03 * tc(2) &
            -1.40894597d-06 * tc(3) &
            +2.60179549d-10 * tc(4) &
            -1.87727567d-14 * tc(5)
        !species 9: CH2(S)
        species(9) = &
            +1.29203842d+00 &
            +4.65588637d-03 * tc(2) &
            -2.01191947d-06 * tc(3) &
            +4.17906000d-10 * tc(4) &
            -3.39716365d-14 * tc(5)
        !species 10: CH3
        species(10) = &
            +1.28571772d+00 &
            +7.23990037d-03 * tc(2) &
            -2.98714348d-06 * tc(3) &
            +5.95684644d-10 * tc(4) &
            -4.67154394d-14 * tc(5)
        !species 11: CH4
        species(11) = &
            -9.25148505d-01 &
            +1.33909467d-02 * tc(2) &
            -5.73285809d-06 * tc(3) &
            +1.22292535d-09 * tc(4) &
            -1.01815230d-13 * tc(5)
        !species 12: CO
        species(12) = &
            +1.71518561d+00 &
            +2.06252743d-03 * tc(2) &
            -9.98825771d-07 * tc(3) &
            +2.30053008d-10 * tc(4) &
            -2.03647716d-14 * tc(5)
        !species 13: CO2
        species(13) = &
            +2.85746029d+00 &
            +4.41437026d-03 * tc(2) &
            -2.21481404d-06 * tc(3) &
            +5.23490188d-10 * tc(4) &
            -4.72084164d-14 * tc(5)
        !species 14: HCO
        species(14) = &
            +1.77217438d+00 &
            +4.95695526d-03 * tc(2) &
            -2.48445613d-06 * tc(3) &
            +5.89161778d-10 * tc(4) &
            -5.33508711d-14 * tc(5)
        !species 15: CH2O
        species(15) = &
            +7.60690080d-01 &
            +9.20000082d-03 * tc(2) &
            -4.42258813d-06 * tc(3) &
            +1.00641212d-09 * tc(4) &
            -8.83855640d-14 * tc(5)
        !species 16: CH3O
        species(16) = &
            +2.77079900d+00 &
            +7.87149700d-03 * tc(2) &
            -2.65638400d-06 * tc(3) &
            +3.94443100d-10 * tc(4) &
            -2.11261600d-14 * tc(5)
        !species 17: C2H4
        species(17) = &
            +1.03611116d+00 &
            +1.46454151d-02 * tc(2) &
            -6.71077915d-06 * tc(3) &
            +1.47222923d-09 * tc(4) &
            -1.25706061d-13 * tc(5)
        !species 18: C2H5
        species(18) = &
            +9.54656420d-01 &
            +1.73972722d-02 * tc(2) &
            -7.98206668d-06 * tc(3) &
            +1.75217689d-09 * tc(4) &
            -1.49641576d-13 * tc(5)
        !species 19: C2H6
        species(19) = &
            +7.18815000d-02 &
            +2.16852677d-02 * tc(2) &
            -1.00256067d-05 * tc(3) &
            +2.21412001d-09 * tc(4) &
            -1.90002890d-13 * tc(5)
        !species 20: N2
        species(20) = &
            +1.92664000d+00 &
            +1.48797680d-03 * tc(2) &
            -5.68476000d-07 * tc(3) &
            +1.00970380d-10 * tc(4) &
            -6.75335100d-15 * tc(5)
        !species 21: AR
        species(21) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
    end if

end subroutine


! compute Cp/R at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine cp_R(species, tc)

    implicit none

    double precision, intent(out) :: species(21)
    double precision, intent(in) :: tc(5)

    double precision :: T

    T = tc(2)

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +2.34433112d+00 &
            +7.98052075d-03 * tc(2) &
            -1.94781510d-05 * tc(3) &
            +2.01572094d-08 * tc(4) &
            -7.37611761d-12 * tc(5)
        ! species 2: H
        species(2) = &
            +2.50000000d+00 &
            +7.05332819d-13 * tc(2) &
            -1.99591964d-15 * tc(3) &
            +2.30081632d-18 * tc(4) &
            -9.27732332d-22 * tc(5)
        ! species 3: O
        species(3) = &
            +3.16826710d+00 &
            -3.27931884d-03 * tc(2) &
            +6.64306396d-06 * tc(3) &
            -6.12806624d-09 * tc(4) &
            +2.11265971d-12 * tc(5)
        ! species 4: O2
        species(4) = &
            +3.78245636d+00 &
            -2.99673416d-03 * tc(2) &
            +9.84730201d-06 * tc(3) &
            -9.68129509d-09 * tc(4) &
            +3.24372837d-12 * tc(5)
        ! species 5: OH
        species(5) = &
            +3.99201543d+00 &
            -2.40131752d-03 * tc(2) &
            +4.61793841d-06 * tc(3) &
            -3.88113333d-09 * tc(4) &
            +1.36411470d-12 * tc(5)
        ! species 6: H2O
        species(6) = &
            +4.19864056d+00 &
            -2.03643410d-03 * tc(2) &
            +6.52040211d-06 * tc(3) &
            -5.48797062d-09 * tc(4) &
            +1.77197817d-12 * tc(5)
        ! species 7: HO2
        species(7) = &
            +4.30179801d+00 &
            -4.74912051d-03 * tc(2) &
            +2.11582891d-05 * tc(3) &
            -2.42763894d-08 * tc(4) &
            +9.29225124d-12 * tc(5)
        ! species 8: CH2
        species(8) = &
            +3.76267867d+00 &
            +9.68872143d-04 * tc(2) &
            +2.79489841d-06 * tc(3) &
            -3.85091153d-09 * tc(4) &
            +1.68741719d-12 * tc(5)
        ! species 9: CH2(S)
        species(9) = &
            +4.19860411d+00 &
            -2.36661419d-03 * tc(2) &
            +8.23296220d-06 * tc(3) &
            -6.68815981d-09 * tc(4) &
            +1.94314737d-12 * tc(5)
        ! species 10: CH3
        species(10) = &
            +3.67359040d+00 &
            +2.01095175d-03 * tc(2) &
            +5.73021856d-06 * tc(3) &
            -6.87117425d-09 * tc(4) &
            +2.54385734d-12 * tc(5)
        ! species 11: CH4
        species(11) = &
            +5.14987613d+00 &
            -1.36709788d-02 * tc(2) &
            +4.91800599d-05 * tc(3) &
            -4.84743026d-08 * tc(4) &
            +1.66693956d-11 * tc(5)
        ! species 12: CO
        species(12) = &
            +3.57953347d+00 &
            -6.10353680d-04 * tc(2) &
            +1.01681433d-06 * tc(3) &
            +9.07005884d-10 * tc(4) &
            -9.04424499d-13 * tc(5)
        ! species 13: CO2
        species(13) = &
            +2.35677352d+00 &
            +8.98459677d-03 * tc(2) &
            -7.12356269d-06 * tc(3) &
            +2.45919022d-09 * tc(4) &
            -1.43699548d-13 * tc(5)
        ! species 14: HCO
        species(14) = &
            +4.22118584d+00 &
            -3.24392532d-03 * tc(2) &
            +1.37799446d-05 * tc(3) &
            -1.33144093d-08 * tc(4) &
            +4.33768865d-12 * tc(5)
        ! species 15: CH2O
        species(15) = &
            +4.79372315d+00 &
            -9.90833369d-03 * tc(2) &
            +3.73220008d-05 * tc(3) &
            -3.79285261d-08 * tc(4) &
            +1.31772652d-11 * tc(5)
        ! species 16: CH3O
        species(16) = &
            +2.10620400d+00 &
            +7.21659500d-03 * tc(2) &
            +5.33847200d-06 * tc(3) &
            -7.37763600d-09 * tc(4) &
            +2.07561000d-12 * tc(5)
        ! species 17: C2H4
        species(17) = &
            +3.95920148d+00 &
            -7.57052247d-03 * tc(2) &
            +5.70990292d-05 * tc(3) &
            -6.91588753d-08 * tc(4) &
            +2.69884373d-11 * tc(5)
        ! species 18: C2H5
        species(18) = &
            +4.30646568d+00 &
            -4.18658892d-03 * tc(2) &
            +4.97142807d-05 * tc(3) &
            -5.99126606d-08 * tc(4) &
            +2.30509004d-11 * tc(5)
        ! species 19: C2H6
        species(19) = &
            +4.29142492d+00 &
            -5.50154270d-03 * tc(2) &
            +5.99438288d-05 * tc(3) &
            -7.08466285d-08 * tc(4) &
            +2.68685771d-11 * tc(5)
        ! species 20: N2
        species(20) = &
            +3.29867700d+00 &
            +1.40824040d-03 * tc(2) &
            -3.96322200d-06 * tc(3) &
            +5.64151500d-09 * tc(4) &
            -2.44485400d-12 * tc(5)
        ! species 21: AR
        species(21) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
    else
        !species 1: H2
        species(1) = &
            +3.33727920d+00 &
            -4.94024731d-05 * tc(2) &
            +4.99456778d-07 * tc(3) &
            -1.79566394d-10 * tc(4) &
            +2.00255376d-14 * tc(5)
        !species 2: H
        species(2) = &
            +2.50000001d+00 &
            -2.30842973d-11 * tc(2) &
            +1.61561948d-14 * tc(3) &
            -4.73515235d-18 * tc(4) &
            +4.98197357d-22 * tc(5)
        !species 3: O
        species(3) = &
            +2.56942078d+00 &
            -8.59741137d-05 * tc(2) &
            +4.19484589d-08 * tc(3) &
            -1.00177799d-11 * tc(4) &
            +1.22833691d-15 * tc(5)
        !species 4: O2
        species(4) = &
            +3.28253784d+00 &
            +1.48308754d-03 * tc(2) &
            -7.57966669d-07 * tc(3) &
            +2.09470555d-10 * tc(4) &
            -2.16717794d-14 * tc(5)
        !species 5: OH
        species(5) = &
            +3.09288767d+00 &
            +5.48429716d-04 * tc(2) &
            +1.26505228d-07 * tc(3) &
            -8.79461556d-11 * tc(4) &
            +1.17412376d-14 * tc(5)
        !species 6: H2O
        species(6) = &
            +3.03399249d+00 &
            +2.17691804d-03 * tc(2) &
            -1.64072518d-07 * tc(3) &
            -9.70419870d-11 * tc(4) &
            +1.68200992d-14 * tc(5)
        !species 7: HO2
        species(7) = &
            +4.01721090d+00 &
            +2.23982013d-03 * tc(2) &
            -6.33658150d-07 * tc(3) &
            +1.14246370d-10 * tc(4) &
            -1.07908535d-14 * tc(5)
        !species 8: CH2
        species(8) = &
            +2.87410113d+00 &
            +3.65639292d-03 * tc(2) &
            -1.40894597d-06 * tc(3) &
            +2.60179549d-10 * tc(4) &
            -1.87727567d-14 * tc(5)
        !species 9: CH2(S)
        species(9) = &
            +2.29203842d+00 &
            +4.65588637d-03 * tc(2) &
            -2.01191947d-06 * tc(3) &
            +4.17906000d-10 * tc(4) &
            -3.39716365d-14 * tc(5)
        !species 10: CH3
        species(10) = &
            +2.28571772d+00 &
            +7.23990037d-03 * tc(2) &
            -2.98714348d-06 * tc(3) &
            +5.95684644d-10 * tc(4) &
            -4.67154394d-14 * tc(5)
        !species 11: CH4
        species(11) = &
            +7.48514950d-02 &
            +1.33909467d-02 * tc(2) &
            -5.73285809d-06 * tc(3) &
            +1.22292535d-09 * tc(4) &
            -1.01815230d-13 * tc(5)
        !species 12: CO
        species(12) = &
            +2.71518561d+00 &
            +2.06252743d-03 * tc(2) &
            -9.98825771d-07 * tc(3) &
            +2.30053008d-10 * tc(4) &
            -2.03647716d-14 * tc(5)
        !species 13: CO2
        species(13) = &
            +3.85746029d+00 &
            +4.41437026d-03 * tc(2) &
            -2.21481404d-06 * tc(3) &
            +5.23490188d-10 * tc(4) &
            -4.72084164d-14 * tc(5)
        !species 14: HCO
        species(14) = &
            +2.77217438d+00 &
            +4.95695526d-03 * tc(2) &
            -2.48445613d-06 * tc(3) &
            +5.89161778d-10 * tc(4) &
            -5.33508711d-14 * tc(5)
        !species 15: CH2O
        species(15) = &
            +1.76069008d+00 &
            +9.20000082d-03 * tc(2) &
            -4.42258813d-06 * tc(3) &
            +1.00641212d-09 * tc(4) &
            -8.83855640d-14 * tc(5)
        !species 16: CH3O
        species(16) = &
            +3.77079900d+00 &
            +7.87149700d-03 * tc(2) &
            -2.65638400d-06 * tc(3) &
            +3.94443100d-10 * tc(4) &
            -2.11261600d-14 * tc(5)
        !species 17: C2H4
        species(17) = &
            +2.03611116d+00 &
            +1.46454151d-02 * tc(2) &
            -6.71077915d-06 * tc(3) &
            +1.47222923d-09 * tc(4) &
            -1.25706061d-13 * tc(5)
        !species 18: C2H5
        species(18) = &
            +1.95465642d+00 &
            +1.73972722d-02 * tc(2) &
            -7.98206668d-06 * tc(3) &
            +1.75217689d-09 * tc(4) &
            -1.49641576d-13 * tc(5)
        !species 19: C2H6
        species(19) = &
            +1.07188150d+00 &
            +2.16852677d-02 * tc(2) &
            -1.00256067d-05 * tc(3) &
            +2.21412001d-09 * tc(4) &
            -1.90002890d-13 * tc(5)
        !species 20: N2
        species(20) = &
            +2.92664000d+00 &
            +1.48797680d-03 * tc(2) &
            -5.68476000d-07 * tc(3) &
            +1.00970380d-10 * tc(4) &
            -6.75335100d-15 * tc(5)
        !species 21: AR
        species(21) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5)
    end if

end subroutine

! compute the e/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine speciesInternalEnergy(species, tc)

    implicit none

    double precision, intent(out) :: species(21)
    double precision, intent(in) :: tc(5)

    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +1.34433112d+00 &
            +3.99026037d-03 * tc(2) &
            -6.49271700d-06 * tc(3) &
            +5.03930235d-09 * tc(4) &
            -1.47522352d-12 * tc(5) &
            -9.17935173d+02 * invT
        ! species 2: H
        species(2) = &
            +1.50000000d+00 &
            +3.52666409d-13 * tc(2) &
            -6.65306547d-16 * tc(3) &
            +5.75204080d-19 * tc(4) &
            -1.85546466d-22 * tc(5) &
            +2.54736599d+04 * invT
        ! species 3: O
        species(3) = &
            +2.16826710d+00 &
            -1.63965942d-03 * tc(2) &
            +2.21435465d-06 * tc(3) &
            -1.53201656d-09 * tc(4) &
            +4.22531942d-13 * tc(5) &
            +2.91222592d+04 * invT
        ! species 4: O2
        species(4) = &
            +2.78245636d+00 &
            -1.49836708d-03 * tc(2) &
            +3.28243400d-06 * tc(3) &
            -2.42032377d-09 * tc(4) &
            +6.48745674d-13 * tc(5) &
            -1.06394356d+03 * invT
        ! species 5: OH
        species(5) = &
            +2.99201543d+00 &
            -1.20065876d-03 * tc(2) &
            +1.53931280d-06 * tc(3) &
            -9.70283332d-10 * tc(4) &
            +2.72822940d-13 * tc(5) &
            +3.61508056d+03 * invT
        ! species 6: H2O
        species(6) = &
            +3.19864056d+00 &
            -1.01821705d-03 * tc(2) &
            +2.17346737d-06 * tc(3) &
            -1.37199266d-09 * tc(4) &
            +3.54395634d-13 * tc(5) &
            -3.02937267d+04 * invT
        ! species 7: HO2
        species(7) = &
            +3.30179801d+00 &
            -2.37456025d-03 * tc(2) &
            +7.05276303d-06 * tc(3) &
            -6.06909735d-09 * tc(4) &
            +1.85845025d-12 * tc(5) &
            +2.94808040d+02 * invT
        ! species 8: CH2
        species(8) = &
            +2.76267867d+00 &
            +4.84436072d-04 * tc(2) &
            +9.31632803d-07 * tc(3) &
            -9.62727883d-10 * tc(4) &
            +3.37483438d-13 * tc(5) &
            +4.60040401d+04 * invT
        ! species 9: CH2(S)
        species(9) = &
            +3.19860411d+00 &
            -1.18330710d-03 * tc(2) &
            +2.74432073d-06 * tc(3) &
            -1.67203995d-09 * tc(4) &
            +3.88629474d-13 * tc(5) &
            +5.04968163d+04 * invT
        ! species 10: CH3
        species(10) = &
            +2.67359040d+00 &
            +1.00547588d-03 * tc(2) &
            +1.91007285d-06 * tc(3) &
            -1.71779356d-09 * tc(4) &
            +5.08771468d-13 * tc(5) &
            +1.64449988d+04 * invT
        ! species 11: CH4
        species(11) = &
            +4.14987613d+00 &
            -6.83548940d-03 * tc(2) &
            +1.63933533d-05 * tc(3) &
            -1.21185757d-08 * tc(4) &
            +3.33387912d-12 * tc(5) &
            -1.02466476d+04 * invT
        ! species 12: CO
        species(12) = &
            +2.57953347d+00 &
            -3.05176840d-04 * tc(2) &
            +3.38938110d-07 * tc(3) &
            +2.26751471d-10 * tc(4) &
            -1.80884900d-13 * tc(5) &
            -1.43440860d+04 * invT
        ! species 13: CO2
        species(13) = &
            +1.35677352d+00 &
            +4.49229839d-03 * tc(2) &
            -2.37452090d-06 * tc(3) &
            +6.14797555d-10 * tc(4) &
            -2.87399096d-14 * tc(5) &
            -4.83719697d+04 * invT
        ! species 14: HCO
        species(14) = &
            +3.22118584d+00 &
            -1.62196266d-03 * tc(2) &
            +4.59331487d-06 * tc(3) &
            -3.32860233d-09 * tc(4) &
            +8.67537730d-13 * tc(5) &
            +3.83956496d+03 * invT
        ! species 15: CH2O
        species(15) = &
            +3.79372315d+00 &
            -4.95416684d-03 * tc(2) &
            +1.24406669d-05 * tc(3) &
            -9.48213152d-09 * tc(4) &
            +2.63545304d-12 * tc(5) &
            -1.43089567d+04 * invT
        ! species 16: CH3O
        species(16) = &
            +1.10620400d+00 &
            +3.60829750d-03 * tc(2) &
            +1.77949067d-06 * tc(3) &
            -1.84440900d-09 * tc(4) &
            +4.15122000d-13 * tc(5) &
            +9.78601100d+02 * invT
        ! species 17: C2H4
        species(17) = &
            +2.95920148d+00 &
            -3.78526124d-03 * tc(2) &
            +1.90330097d-05 * tc(3) &
            -1.72897188d-08 * tc(4) &
            +5.39768746d-12 * tc(5) &
            +5.08977593d+03 * invT
        ! species 18: C2H5
        species(18) = &
            +3.30646568d+00 &
            -2.09329446d-03 * tc(2) &
            +1.65714269d-05 * tc(3) &
            -1.49781651d-08 * tc(4) &
            +4.61018008d-12 * tc(5) &
            +1.28416265d+04 * invT
        ! species 19: C2H6
        species(19) = &
            +3.29142492d+00 &
            -2.75077135d-03 * tc(2) &
            +1.99812763d-05 * tc(3) &
            -1.77116571d-08 * tc(4) &
            +5.37371542d-12 * tc(5) &
            -1.15222055d+04 * invT
        ! species 20: N2
        species(20) = &
            +2.29867700d+00 &
            +7.04120200d-04 * tc(2) &
            -1.32107400d-06 * tc(3) &
            +1.41037875d-09 * tc(4) &
            -4.88970800d-13 * tc(5) &
            -1.02089990d+03 * invT
        ! species 21: AR
        species(21) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            -7.45375000d+02 * invT
    else
        !species 1: H2
        species(1) = &
            +2.33727920d+00 &
            -2.47012365d-05 * tc(2) &
            +1.66485593d-07 * tc(3) &
            -4.48915985d-11 * tc(4) &
            +4.00510752d-15 * tc(5) &
            -9.50158922d+02 * invT
        !species 2: H
        species(2) = &
            +1.50000001d+00 &
            -1.15421486d-11 * tc(2) &
            +5.38539827d-15 * tc(3) &
            -1.18378809d-18 * tc(4) &
            +9.96394714d-23 * tc(5) &
            +2.54736599d+04 * invT
        !species 3: O
        species(3) = &
            +1.56942078d+00 &
            -4.29870569d-05 * tc(2) &
            +1.39828196d-08 * tc(3) &
            -2.50444497d-12 * tc(4) &
            +2.45667382d-16 * tc(5) &
            +2.92175791d+04 * invT
        !species 4: O2
        species(4) = &
            +2.28253784d+00 &
            +7.41543770d-04 * tc(2) &
            -2.52655556d-07 * tc(3) &
            +5.23676387d-11 * tc(4) &
            -4.33435588d-15 * tc(5) &
            -1.08845772d+03 * invT
        !species 5: OH
        species(5) = &
            +2.09288767d+00 &
            +2.74214858d-04 * tc(2) &
            +4.21684093d-08 * tc(3) &
            -2.19865389d-11 * tc(4) &
            +2.34824752d-15 * tc(5) &
            +3.85865700d+03 * invT
        !species 6: H2O
        species(6) = &
            +2.03399249d+00 &
            +1.08845902d-03 * tc(2) &
            -5.46908393d-08 * tc(3) &
            -2.42604967d-11 * tc(4) &
            +3.36401984d-15 * tc(5) &
            -3.00042971d+04 * invT
        !species 7: HO2
        species(7) = &
            +3.01721090d+00 &
            +1.11991006d-03 * tc(2) &
            -2.11219383d-07 * tc(3) &
            +2.85615925d-11 * tc(4) &
            -2.15817070d-15 * tc(5) &
            +1.11856713d+02 * invT
        !species 8: CH2
        species(8) = &
            +1.87410113d+00 &
            +1.82819646d-03 * tc(2) &
            -4.69648657d-07 * tc(3) &
            +6.50448872d-11 * tc(4) &
            -3.75455134d-15 * tc(5) &
            +4.62636040d+04 * invT
        !species 9: CH2(S)
        species(9) = &
            +1.29203842d+00 &
            +2.32794318d-03 * tc(2) &
            -6.70639823d-07 * tc(3) &
            +1.04476500d-10 * tc(4) &
            -6.79432730d-15 * tc(5) &
            +5.09259997d+04 * invT
        !species 10: CH3
        species(10) = &
            +1.28571772d+00 &
            +3.61995018d-03 * tc(2) &
            -9.95714493d-07 * tc(3) &
            +1.48921161d-10 * tc(4) &
            -9.34308788d-15 * tc(5) &
            +1.67755843d+04 * invT
        !species 11: CH4
        species(11) = &
            -9.25148505d-01 &
            +6.69547335d-03 * tc(2) &
            -1.91095270d-06 * tc(3) &
            +3.05731338d-10 * tc(4) &
            -2.03630460d-14 * tc(5) &
            -9.46834459d+03 * invT
        !species 12: CO
        species(12) = &
            +1.71518561d+00 &
            +1.03126372d-03 * tc(2) &
            -3.32941924d-07 * tc(3) &
            +5.75132520d-11 * tc(4) &
            -4.07295432d-15 * tc(5) &
            -1.41518724d+04 * invT
        !species 13: CO2
        species(13) = &
            +2.85746029d+00 &
            +2.20718513d-03 * tc(2) &
            -7.38271347d-07 * tc(3) &
            +1.30872547d-10 * tc(4) &
            -9.44168328d-15 * tc(5) &
            -4.87591660d+04 * invT
        !species 14: HCO
        species(14) = &
            +1.77217438d+00 &
            +2.47847763d-03 * tc(2) &
            -8.28152043d-07 * tc(3) &
            +1.47290445d-10 * tc(4) &
            -1.06701742d-14 * tc(5) &
            +4.01191815d+03 * invT
        !species 15: CH2O
        species(15) = &
            +7.60690080d-01 &
            +4.60000041d-03 * tc(2) &
            -1.47419604d-06 * tc(3) &
            +2.51603030d-10 * tc(4) &
            -1.76771128d-14 * tc(5) &
            -1.39958323d+04 * invT
        !species 16: CH3O
        species(16) = &
            +2.77079900d+00 &
            +3.93574850d-03 * tc(2) &
            -8.85461333d-07 * tc(3) &
            +9.86107750d-11 * tc(4) &
            -4.22523200d-15 * tc(5) &
            +1.27832520d+02 * invT
        !species 17: C2H4
        species(17) = &
            +1.03611116d+00 &
            +7.32270755d-03 * tc(2) &
            -2.23692638d-06 * tc(3) &
            +3.68057308d-10 * tc(4) &
            -2.51412122d-14 * tc(5) &
            +4.93988614d+03 * invT
        !species 18: C2H5
        species(18) = &
            +9.54656420d-01 &
            +8.69863610d-03 * tc(2) &
            -2.66068889d-06 * tc(3) &
            +4.38044223d-10 * tc(4) &
            -2.99283152d-14 * tc(5) &
            +1.28575200d+04 * invT
        !species 19: C2H6
        species(19) = &
            +7.18815000d-02 &
            +1.08426339d-02 * tc(2) &
            -3.34186890d-06 * tc(3) &
            +5.53530003d-10 * tc(4) &
            -3.80005780d-14 * tc(5) &
            -1.14263932d+04 * invT
        !species 20: N2
        species(20) = &
            +1.92664000d+00 &
            +7.43988400d-04 * tc(2) &
            -1.89492000d-07 * tc(3) &
            +2.52425950d-11 * tc(4) &
            -1.35067020d-15 * tc(5) &
            -9.22797700d+02 * invT
        !species 21: AR
        species(21) = &
            +1.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            -7.45375000d+02 * invT
    end if

end subroutine

! compute the h/(RT) at the given temperature (Eq 20)
! tc contains precomputed powers of T, tc(1) = log(T)
subroutine speciesEnthalpy(species, tc)

    implicit none

    double precision, intent(out) :: species(21)
    double precision, intent(in) :: tc(5)

    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
        ! species 1: H2
        species(1) = &
            +2.34433112d+00 &
            +3.99026037d-03 * tc(2) &
            -6.49271700d-06 * tc(3) &
            +5.03930235d-09 * tc(4) &
            -1.47522352d-12 * tc(5) &
            -9.17935173d+02 * invT
        ! species 2: H
        species(2) = &
            +2.50000000d+00 &
            +3.52666409d-13 * tc(2) &
            -6.65306547d-16 * tc(3) &
            +5.75204080d-19 * tc(4) &
            -1.85546466d-22 * tc(5) &
            +2.54736599d+04 * invT
        ! species 3: O
        species(3) = &
            +3.16826710d+00 &
            -1.63965942d-03 * tc(2) &
            +2.21435465d-06 * tc(3) &
            -1.53201656d-09 * tc(4) &
            +4.22531942d-13 * tc(5) &
            +2.91222592d+04 * invT
        ! species 4: O2
        species(4) = &
            +3.78245636d+00 &
            -1.49836708d-03 * tc(2) &
            +3.28243400d-06 * tc(3) &
            -2.42032377d-09 * tc(4) &
            +6.48745674d-13 * tc(5) &
            -1.06394356d+03 * invT
        ! species 5: OH
        species(5) = &
            +3.99201543d+00 &
            -1.20065876d-03 * tc(2) &
            +1.53931280d-06 * tc(3) &
            -9.70283332d-10 * tc(4) &
            +2.72822940d-13 * tc(5) &
            +3.61508056d+03 * invT
        ! species 6: H2O
        species(6) = &
            +4.19864056d+00 &
            -1.01821705d-03 * tc(2) &
            +2.17346737d-06 * tc(3) &
            -1.37199266d-09 * tc(4) &
            +3.54395634d-13 * tc(5) &
            -3.02937267d+04 * invT
        ! species 7: HO2
        species(7) = &
            +4.30179801d+00 &
            -2.37456025d-03 * tc(2) &
            +7.05276303d-06 * tc(3) &
            -6.06909735d-09 * tc(4) &
            +1.85845025d-12 * tc(5) &
            +2.94808040d+02 * invT
        ! species 8: CH2
        species(8) = &
            +3.76267867d+00 &
            +4.84436072d-04 * tc(2) &
            +9.31632803d-07 * tc(3) &
            -9.62727883d-10 * tc(4) &
            +3.37483438d-13 * tc(5) &
            +4.60040401d+04 * invT
        ! species 9: CH2(S)
        species(9) = &
            +4.19860411d+00 &
            -1.18330710d-03 * tc(2) &
            +2.74432073d-06 * tc(3) &
            -1.67203995d-09 * tc(4) &
            +3.88629474d-13 * tc(5) &
            +5.04968163d+04 * invT
        ! species 10: CH3
        species(10) = &
            +3.67359040d+00 &
            +1.00547588d-03 * tc(2) &
            +1.91007285d-06 * tc(3) &
            -1.71779356d-09 * tc(4) &
            +5.08771468d-13 * tc(5) &
            +1.64449988d+04 * invT
        ! species 11: CH4
        species(11) = &
            +5.14987613d+00 &
            -6.83548940d-03 * tc(2) &
            +1.63933533d-05 * tc(3) &
            -1.21185757d-08 * tc(4) &
            +3.33387912d-12 * tc(5) &
            -1.02466476d+04 * invT
        ! species 12: CO
        species(12) = &
            +3.57953347d+00 &
            -3.05176840d-04 * tc(2) &
            +3.38938110d-07 * tc(3) &
            +2.26751471d-10 * tc(4) &
            -1.80884900d-13 * tc(5) &
            -1.43440860d+04 * invT
        ! species 13: CO2
        species(13) = &
            +2.35677352d+00 &
            +4.49229839d-03 * tc(2) &
            -2.37452090d-06 * tc(3) &
            +6.14797555d-10 * tc(4) &
            -2.87399096d-14 * tc(5) &
            -4.83719697d+04 * invT
        ! species 14: HCO
        species(14) = &
            +4.22118584d+00 &
            -1.62196266d-03 * tc(2) &
            +4.59331487d-06 * tc(3) &
            -3.32860233d-09 * tc(4) &
            +8.67537730d-13 * tc(5) &
            +3.83956496d+03 * invT
        ! species 15: CH2O
        species(15) = &
            +4.79372315d+00 &
            -4.95416684d-03 * tc(2) &
            +1.24406669d-05 * tc(3) &
            -9.48213152d-09 * tc(4) &
            +2.63545304d-12 * tc(5) &
            -1.43089567d+04 * invT
        ! species 16: CH3O
        species(16) = &
            +2.10620400d+00 &
            +3.60829750d-03 * tc(2) &
            +1.77949067d-06 * tc(3) &
            -1.84440900d-09 * tc(4) &
            +4.15122000d-13 * tc(5) &
            +9.78601100d+02 * invT
        ! species 17: C2H4
        species(17) = &
            +3.95920148d+00 &
            -3.78526124d-03 * tc(2) &
            +1.90330097d-05 * tc(3) &
            -1.72897188d-08 * tc(4) &
            +5.39768746d-12 * tc(5) &
            +5.08977593d+03 * invT
        ! species 18: C2H5
        species(18) = &
            +4.30646568d+00 &
            -2.09329446d-03 * tc(2) &
            +1.65714269d-05 * tc(3) &
            -1.49781651d-08 * tc(4) &
            +4.61018008d-12 * tc(5) &
            +1.28416265d+04 * invT
        ! species 19: C2H6
        species(19) = &
            +4.29142492d+00 &
            -2.75077135d-03 * tc(2) &
            +1.99812763d-05 * tc(3) &
            -1.77116571d-08 * tc(4) &
            +5.37371542d-12 * tc(5) &
            -1.15222055d+04 * invT
        ! species 20: N2
        species(20) = &
            +3.29867700d+00 &
            +7.04120200d-04 * tc(2) &
            -1.32107400d-06 * tc(3) &
            +1.41037875d-09 * tc(4) &
            -4.88970800d-13 * tc(5) &
            -1.02089990d+03 * invT
        ! species 21: AR
        species(21) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            -7.45375000d+02 * invT
    else
        !species 1: H2
        species(1) = &
            +3.33727920d+00 &
            -2.47012365d-05 * tc(2) &
            +1.66485593d-07 * tc(3) &
            -4.48915985d-11 * tc(4) &
            +4.00510752d-15 * tc(5) &
            -9.50158922d+02 * invT
        !species 2: H
        species(2) = &
            +2.50000001d+00 &
            -1.15421486d-11 * tc(2) &
            +5.38539827d-15 * tc(3) &
            -1.18378809d-18 * tc(4) &
            +9.96394714d-23 * tc(5) &
            +2.54736599d+04 * invT
        !species 3: O
        species(3) = &
            +2.56942078d+00 &
            -4.29870569d-05 * tc(2) &
            +1.39828196d-08 * tc(3) &
            -2.50444497d-12 * tc(4) &
            +2.45667382d-16 * tc(5) &
            +2.92175791d+04 * invT
        !species 4: O2
        species(4) = &
            +3.28253784d+00 &
            +7.41543770d-04 * tc(2) &
            -2.52655556d-07 * tc(3) &
            +5.23676387d-11 * tc(4) &
            -4.33435588d-15 * tc(5) &
            -1.08845772d+03 * invT
        !species 5: OH
        species(5) = &
            +3.09288767d+00 &
            +2.74214858d-04 * tc(2) &
            +4.21684093d-08 * tc(3) &
            -2.19865389d-11 * tc(4) &
            +2.34824752d-15 * tc(5) &
            +3.85865700d+03 * invT
        !species 6: H2O
        species(6) = &
            +3.03399249d+00 &
            +1.08845902d-03 * tc(2) &
            -5.46908393d-08 * tc(3) &
            -2.42604967d-11 * tc(4) &
            +3.36401984d-15 * tc(5) &
            -3.00042971d+04 * invT
        !species 7: HO2
        species(7) = &
            +4.01721090d+00 &
            +1.11991006d-03 * tc(2) &
            -2.11219383d-07 * tc(3) &
            +2.85615925d-11 * tc(4) &
            -2.15817070d-15 * tc(5) &
            +1.11856713d+02 * invT
        !species 8: CH2
        species(8) = &
            +2.87410113d+00 &
            +1.82819646d-03 * tc(2) &
            -4.69648657d-07 * tc(3) &
            +6.50448872d-11 * tc(4) &
            -3.75455134d-15 * tc(5) &
            +4.62636040d+04 * invT
        !species 9: CH2(S)
        species(9) = &
            +2.29203842d+00 &
            +2.32794318d-03 * tc(2) &
            -6.70639823d-07 * tc(3) &
            +1.04476500d-10 * tc(4) &
            -6.79432730d-15 * tc(5) &
            +5.09259997d+04 * invT
        !species 10: CH3
        species(10) = &
            +2.28571772d+00 &
            +3.61995018d-03 * tc(2) &
            -9.95714493d-07 * tc(3) &
            +1.48921161d-10 * tc(4) &
            -9.34308788d-15 * tc(5) &
            +1.67755843d+04 * invT
        !species 11: CH4
        species(11) = &
            +7.48514950d-02 &
            +6.69547335d-03 * tc(2) &
            -1.91095270d-06 * tc(3) &
            +3.05731338d-10 * tc(4) &
            -2.03630460d-14 * tc(5) &
            -9.46834459d+03 * invT
        !species 12: CO
        species(12) = &
            +2.71518561d+00 &
            +1.03126372d-03 * tc(2) &
            -3.32941924d-07 * tc(3) &
            +5.75132520d-11 * tc(4) &
            -4.07295432d-15 * tc(5) &
            -1.41518724d+04 * invT
        !species 13: CO2
        species(13) = &
            +3.85746029d+00 &
            +2.20718513d-03 * tc(2) &
            -7.38271347d-07 * tc(3) &
            +1.30872547d-10 * tc(4) &
            -9.44168328d-15 * tc(5) &
            -4.87591660d+04 * invT
        !species 14: HCO
        species(14) = &
            +2.77217438d+00 &
            +2.47847763d-03 * tc(2) &
            -8.28152043d-07 * tc(3) &
            +1.47290445d-10 * tc(4) &
            -1.06701742d-14 * tc(5) &
            +4.01191815d+03 * invT
        !species 15: CH2O
        species(15) = &
            +1.76069008d+00 &
            +4.60000041d-03 * tc(2) &
            -1.47419604d-06 * tc(3) &
            +2.51603030d-10 * tc(4) &
            -1.76771128d-14 * tc(5) &
            -1.39958323d+04 * invT
        !species 16: CH3O
        species(16) = &
            +3.77079900d+00 &
            +3.93574850d-03 * tc(2) &
            -8.85461333d-07 * tc(3) &
            +9.86107750d-11 * tc(4) &
            -4.22523200d-15 * tc(5) &
            +1.27832520d+02 * invT
        !species 17: C2H4
        species(17) = &
            +2.03611116d+00 &
            +7.32270755d-03 * tc(2) &
            -2.23692638d-06 * tc(3) &
            +3.68057308d-10 * tc(4) &
            -2.51412122d-14 * tc(5) &
            +4.93988614d+03 * invT
        !species 18: C2H5
        species(18) = &
            +1.95465642d+00 &
            +8.69863610d-03 * tc(2) &
            -2.66068889d-06 * tc(3) &
            +4.38044223d-10 * tc(4) &
            -2.99283152d-14 * tc(5) &
            +1.28575200d+04 * invT
        !species 19: C2H6
        species(19) = &
            +1.07188150d+00 &
            +1.08426339d-02 * tc(2) &
            -3.34186890d-06 * tc(3) &
            +5.53530003d-10 * tc(4) &
            -3.80005780d-14 * tc(5) &
            -1.14263932d+04 * invT
        !species 20: N2
        species(20) = &
            +2.92664000d+00 &
            +7.43988400d-04 * tc(2) &
            -1.89492000d-07 * tc(3) &
            +2.52425950d-11 * tc(4) &
            -1.35067020d-15 * tc(5) &
            -9.22797700d+02 * invT
        !species 21: AR
        species(21) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4) &
            +0.00000000d+00 * tc(5) &
            -7.45375000d+02 * invT
    end if

end subroutine

! save molecular weights into array
subroutine molecularWeight(wt)

    implicit none

    double precision, intent(out) :: wt(21)

    wt(1) = 2.015940d0 ! H2
    wt(2) = 1.007970d0 ! H
    wt(3) = 15.999400d0 ! O
    wt(4) = 31.998800d0 ! O2
    wt(5) = 17.007370d0 ! OH
    wt(6) = 18.015340d0 ! H2O
    wt(7) = 33.006770d0 ! HO2
    wt(8) = 14.027090d0 ! CH2
    wt(9) = 14.027090d0 ! CH2(S)
    wt(10) = 15.035060d0 ! CH3
    wt(11) = 16.043030d0 ! CH4
    wt(12) = 28.010550d0 ! CO
    wt(13) = 44.009950d0 ! CO2
    wt(14) = 29.018520d0 ! HCO
    wt(15) = 30.026490d0 ! CH2O
    wt(16) = 31.034460d0 ! CH3O
    wt(17) = 28.054180d0 ! C2H4
    wt(18) = 29.062150d0 ! C2H5
    wt(19) = 30.070120d0 ! C2H6
    wt(20) = 28.013400d0 ! N2
    wt(21) = 39.948000d0 ! AR

end subroutine

! get temperature given internal energy in mass units and mass fracs
subroutine get_t_given_ey(e, y, iwrk, rwrk, t, ierr)

    implicit none

    double precision, intent(in) :: e
    double precision, intent(in) :: y(21)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
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

    call ckubms(tmin, y, iwrk, rwrk, emin)
    call ckubms(tmax, y, iwrk, rwrk, emax)

    if (ein < emin) then
        ! Linear Extrapolation below tmin
        call ckcvbs(tmin, y, iwrk, rwrk, cv)
        t = tmin - (emin-ein) / cv
        ierr = 1
        return
    end if
    if (ein > emax) then
        ! Linear Extrapolation above tmax
        call ckcvbs(tmax, y, iwrk, rwrk, cv)
        t = tmax - (emax - ein) / cv
        ierr = 1
        return
    end if
    t1 = t
    if (t1 < tmin .or. t1 > tmax) then
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin)
    end if
    do i=1, maxiter
        call ckubms(t1,y,iwrk,rwrk,e1)
        call ckcvbs(t1,y,iwrk,rwrk,cv)
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

    LENIMC = 86

end subroutine

subroutine egtransetLENRMC(LENRMC)

    implicit none

    integer, intent(out) :: LENRMC

    LENRMC = 9114

end subroutine

subroutine egtransetNO(NO)

    implicit none

    integer, intent(out) :: NO

    NO = 4

end subroutine

subroutine egtransetKK(KK)

    implicit none

    integer, intent(out) :: KK

    KK = 21

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

    double precision, intent(out) :: WT(21)

    WT(1) = 2.01594000d+00
    WT(2) = 1.00797000d+00
    WT(3) = 1.59994000d+01
    WT(4) = 3.19988000d+01
    WT(5) = 1.70073700d+01
    WT(6) = 1.80153400d+01
    WT(7) = 3.30067700d+01
    WT(8) = 1.40270900d+01
    WT(9) = 1.40270900d+01
    WT(10) = 1.50350600d+01
    WT(11) = 1.60430300d+01
    WT(12) = 2.80105500d+01
    WT(13) = 4.40099500d+01
    WT(14) = 2.90185200d+01
    WT(15) = 3.00264900d+01
    WT(16) = 3.10344600d+01
    WT(17) = 2.80541800d+01
    WT(18) = 2.90621500d+01
    WT(19) = 3.00701200d+01
    WT(20) = 2.80134000d+01
    WT(21) = 3.99480000d+01

end subroutine

! the lennard-jones potential well depth eps/kb in K
subroutine egtransetEPS(EPS)

    implicit none

    double precision, intent(out) :: EPS(21)

    EPS(12) = 9.81000000d+01
    EPS(20) = 9.75300000d+01
    EPS(21) = 1.36500000d+02
    EPS(8) = 1.44000000d+02
    EPS(13) = 2.44000000d+02
    EPS(9) = 1.44000000d+02
    EPS(14) = 4.98000000d+02
    EPS(11) = 1.41400000d+02
    EPS(16) = 4.17000000d+02
    EPS(2) = 1.45000000d+02
    EPS(17) = 2.80800000d+02
    EPS(7) = 1.07400000d+02
    EPS(19) = 2.52300000d+02
    EPS(1) = 3.80000000d+01
    EPS(5) = 8.00000000d+01
    EPS(4) = 1.07400000d+02
    EPS(18) = 2.52300000d+02
    EPS(6) = 5.72400000d+02
    EPS(15) = 4.98000000d+02
    EPS(3) = 8.00000000d+01
    EPS(10) = 1.44000000d+02

end subroutine

! the lennard-jones collision diameter in Angstroms
subroutine egtransetSIG(SIG)

    implicit none

    double precision, intent(out) :: SIG(21)

    SIG(12) = 3.65000000d+00
    SIG(20) = 3.62100000d+00
    SIG(21) = 3.33000000d+00
    SIG(8) = 3.80000000d+00
    SIG(13) = 3.76300000d+00
    SIG(9) = 3.80000000d+00
    SIG(14) = 3.59000000d+00
    SIG(11) = 3.74600000d+00
    SIG(16) = 3.69000000d+00
    SIG(2) = 2.05000000d+00
    SIG(17) = 3.97100000d+00
    SIG(7) = 3.45800000d+00
    SIG(19) = 4.30200000d+00
    SIG(1) = 2.92000000d+00
    SIG(5) = 2.75000000d+00
    SIG(4) = 3.45800000d+00
    SIG(18) = 4.30200000d+00
    SIG(6) = 2.60500000d+00
    SIG(15) = 3.59000000d+00
    SIG(3) = 2.75000000d+00
    SIG(10) = 3.80000000d+00

end subroutine

! the dipole moment in Debye
subroutine egtransetDIP(DIP)

    implicit none

    double precision, intent(out) :: DIP(21)

    DIP(12) = 0.00000000d+00
    DIP(20) = 0.00000000d+00
    DIP(21) = 0.00000000d+00
    DIP(8) = 0.00000000d+00
    DIP(13) = 0.00000000d+00
    DIP(9) = 0.00000000d+00
    DIP(14) = 0.00000000d+00
    DIP(11) = 0.00000000d+00
    DIP(16) = 1.70000000d+00
    DIP(2) = 0.00000000d+00
    DIP(17) = 0.00000000d+00
    DIP(7) = 0.00000000d+00
    DIP(19) = 0.00000000d+00
    DIP(1) = 0.00000000d+00
    DIP(5) = 0.00000000d+00
    DIP(4) = 0.00000000d+00
    DIP(18) = 0.00000000d+00
    DIP(6) = 1.84400000d+00
    DIP(15) = 0.00000000d+00
    DIP(3) = 0.00000000d+00
    DIP(10) = 0.00000000d+00

end subroutine

! the polarizability in cubic Angstroms
subroutine egtransetPOL(POL)

    implicit none

    double precision, intent(out) :: POL(21)

    POL(12) = 1.95000000d+00
    POL(20) = 1.76000000d+00
    POL(21) = 0.00000000d+00
    POL(8) = 0.00000000d+00
    POL(13) = 2.65000000d+00
    POL(9) = 0.00000000d+00
    POL(14) = 0.00000000d+00
    POL(11) = 2.60000000d+00
    POL(16) = 0.00000000d+00
    POL(2) = 0.00000000d+00
    POL(17) = 0.00000000d+00
    POL(7) = 0.00000000d+00
    POL(19) = 0.00000000d+00
    POL(1) = 7.90000000d-01
    POL(5) = 0.00000000d+00
    POL(4) = 1.60000000d+00
    POL(18) = 0.00000000d+00
    POL(6) = 0.00000000d+00
    POL(15) = 0.00000000d+00
    POL(3) = 0.00000000d+00
    POL(10) = 0.00000000d+00

end subroutine

! the rotational relaxation collision number at 298 K
subroutine egtransetZROT(ZROT)

    implicit none

    double precision, intent(out) :: ZROT(21)

    ZROT(12) = 1.80000000d+00
    ZROT(20) = 4.00000000d+00
    ZROT(21) = 0.00000000d+00
    ZROT(8) = 0.00000000d+00
    ZROT(13) = 2.10000000d+00
    ZROT(9) = 0.00000000d+00
    ZROT(14) = 0.00000000d+00
    ZROT(11) = 1.30000000d+01
    ZROT(16) = 2.00000000d+00
    ZROT(2) = 0.00000000d+00
    ZROT(17) = 1.50000000d+00
    ZROT(7) = 1.00000000d+00
    ZROT(19) = 1.50000000d+00
    ZROT(1) = 2.80000000d+02
    ZROT(5) = 0.00000000d+00
    ZROT(4) = 3.80000000d+00
    ZROT(18) = 1.50000000d+00
    ZROT(6) = 4.00000000d+00
    ZROT(15) = 2.00000000d+00
    ZROT(3) = 0.00000000d+00
    ZROT(10) = 0.00000000d+00

end subroutine

! 0: monoatomic, 1: linear, 2: nonlinear
subroutine egtransetNLIN(NLIN)

    implicit none

    integer, intent(out) :: NLIN(21)

    NLIN(12) = 1
    NLIN(20) = 1
    NLIN(21) = 0
    NLIN(8) = 1
    NLIN(13) = 1
    NLIN(9) = 1
    NLIN(14) = 2
    NLIN(11) = 2
    NLIN(16) = 2
    NLIN(2) = 0
    NLIN(17) = 2
    NLIN(7) = 2
    NLIN(19) = 2
    NLIN(1) = 1
    NLIN(5) = 1
    NLIN(4) = 1
    NLIN(18) = 2
    NLIN(6) = 2
    NLIN(15) = 2
    NLIN(3) = 0
    NLIN(10) = 1

end subroutine


! Poly fits for the viscosities, dim NO*KK
subroutine egtransetCOFETA(COFETA)

    implicit none

    double precision, intent(out) :: COFETA(84)

    COFETA(1) = -1.38347699d+01
    COFETA(2) = 1.00106621d+00
    COFETA(3) = -4.98105694d-02
    COFETA(4) = 2.31450475d-03
    COFETA(5) = -2.04078397d+01
    COFETA(6) = 3.65436395d+00
    COFETA(7) = -3.98339635d-01
    COFETA(8) = 1.75883009d-02
    COFETA(9) = -1.50926240d+01
    COFETA(10) = 1.92606504d+00
    COFETA(11) = -1.73487476d-01
    COFETA(12) = 7.82572931d-03
    COFETA(13) = -1.71618309d+01
    COFETA(14) = 2.68036374d+00
    COFETA(15) = -2.72570227d-01
    COFETA(16) = 1.21650964d-02
    COFETA(17) = -1.50620763d+01
    COFETA(18) = 1.92606504d+00
    COFETA(19) = -1.73487476d-01
    COFETA(20) = 7.82572931d-03
    COFETA(21) = -1.05420863d+01
    COFETA(22) = -1.37777096d+00
    COFETA(23) = 4.20502308d-01
    COFETA(24) = -2.40627230d-02
    COFETA(25) = -1.71463238d+01
    COFETA(26) = 2.68036374d+00
    COFETA(27) = -2.72570227d-01
    COFETA(28) = 1.21650964d-02
    COFETA(29) = -2.02663469d+01
    COFETA(30) = 3.63241793d+00
    COFETA(31) = -3.95581049d-01
    COFETA(32) = 1.74725495d-02
    COFETA(33) = -2.02663469d+01
    COFETA(34) = 3.63241793d+00
    COFETA(35) = -3.95581049d-01
    COFETA(36) = 1.74725495d-02
    COFETA(37) = -2.02316497d+01
    COFETA(38) = 3.63241793d+00
    COFETA(39) = -3.95581049d-01
    COFETA(40) = 1.74725495d-02
    COFETA(41) = -2.00094664d+01
    COFETA(42) = 3.57220167d+00
    COFETA(43) = -3.87936446d-01
    COFETA(44) = 1.71483254d-02
    COFETA(45) = -1.66188336d+01
    COFETA(46) = 2.40307799d+00
    COFETA(47) = -2.36167638d-01
    COFETA(48) = 1.05714061d-02
    COFETA(49) = -2.40014975d+01
    COFETA(50) = 5.14359547d+00
    COFETA(51) = -5.74269731d-01
    COFETA(52) = 2.44937679d-02
    COFETA(53) = -1.98501306d+01
    COFETA(54) = 2.69480162d+00
    COFETA(55) = -1.65880845d-01
    COFETA(56) = 3.14504769d-03
    COFETA(57) = -1.98330577d+01
    COFETA(58) = 2.69480162d+00
    COFETA(59) = -1.65880845d-01
    COFETA(60) = 3.14504769d-03
    COFETA(61) = -1.99945919d+01
    COFETA(62) = 2.86923313d+00
    COFETA(63) = -2.03325661d-01
    COFETA(64) = 5.39056989d-03
    COFETA(65) = -2.50655444d+01
    COFETA(66) = 5.33982977d+00
    COFETA(67) = -5.89982992d-01
    COFETA(68) = 2.47780650d-02
    COFETA(69) = -2.46581414d+01
    COFETA(70) = 5.19497183d+00
    COFETA(71) = -5.78827339d-01
    COFETA(72) = 2.46050921d-02
    COFETA(73) = -2.46410937d+01
    COFETA(74) = 5.19497183d+00
    COFETA(75) = -5.78827339d-01
    COFETA(76) = 2.46050921d-02
    COFETA(77) = -1.65695594d+01
    COFETA(78) = 2.39056562d+00
    COFETA(79) = -2.34558144d-01
    COFETA(80) = 1.05024037d-02
    COFETA(81) = -1.90422907d+01
    COFETA(82) = 3.47025711d+00
    COFETA(83) = -3.75102111d-01
    COFETA(84) = 1.66086076d-02

end subroutine


! Poly fits for the conductivities, dim NO*KK
subroutine egtransetCOFLAM(COFLAM)

    implicit none

    double precision, intent(out) :: COFLAM(84)

    COFLAM(1) = 9.24084392d+00
    COFLAM(2) = -4.69568931d-01
    COFLAM(3) = 1.15980279d-01
    COFLAM(4) = -2.61033830d-03
    COFLAM(5) = -8.57929284d-01
    COFLAM(6) = 3.65436395d+00
    COFLAM(7) = -3.98339635d-01
    COFLAM(8) = 1.75883009d-02
    COFLAM(9) = 1.69267361d+00
    COFLAM(10) = 1.92606504d+00
    COFLAM(11) = -1.73487476d-01
    COFLAM(12) = 7.82572931d-03
    COFLAM(13) = -1.93718739d+00
    COFLAM(14) = 2.89110219d+00
    COFLAM(15) = -2.71096923d-01
    COFLAM(16) = 1.15344907d-02
    COFLAM(17) = 1.41666185d+01
    COFLAM(18) = -3.24630496d+00
    COFLAM(19) = 5.33878183d-01
    COFLAM(20) = -2.32905165d-02
    COFLAM(21) = 2.33729817d+01
    COFLAM(22) = -8.96536433d+00
    COFLAM(23) = 1.52828665d+00
    COFLAM(24) = -7.58551778d-02
    COFLAM(25) = -1.12960913d+00
    COFLAM(26) = 2.34014305d+00
    COFLAM(27) = -1.63245030d-01
    COFLAM(28) = 5.80319600d-03
    COFLAM(29) = 1.29177902d+01
    COFLAM(30) = -3.73745535d+00
    COFLAM(31) = 7.15831021d-01
    COFLAM(32) = -3.63846910d-02
    COFLAM(33) = 1.89383266d+01
    COFLAM(34) = -6.51018128d+00
    COFLAM(35) = 1.13292060d+00
    COFLAM(36) = -5.69603226d-02
    COFLAM(37) = 1.39937901d+01
    COFLAM(38) = -4.64256494d+00
    COFLAM(39) = 9.07728674d-01
    COFLAM(40) = -4.77274469d-02
    COFLAM(41) = 1.33091602d+01
    COFLAM(42) = -4.96140261d+00
    COFLAM(43) = 1.03295505d+00
    COFLAM(44) = -5.63420090d-02
    COFLAM(45) = 1.18777264d+01
    COFLAM(46) = -3.15463949d+00
    COFLAM(47) = 6.01973268d-01
    COFLAM(48) = -3.03211261d-02
    COFLAM(49) = -1.13649314d+01
    COFLAM(50) = 5.88177395d+00
    COFLAM(51) = -5.68651819d-01
    COFLAM(52) = 2.03561485d-02
    COFLAM(53) = 6.30243508d+00
    COFLAM(54) = -2.22810801d+00
    COFLAM(55) = 6.37340514d-01
    COFLAM(56) = -3.81056018d-02
    COFLAM(57) = 5.39305086d+00
    COFLAM(58) = -2.39312375d+00
    COFLAM(59) = 7.39585221d-01
    COFLAM(60) = -4.58435589d-02
    COFLAM(61) = -6.14588187d+00
    COFLAM(62) = 2.47428873d+00
    COFLAM(63) = 6.43999571d-02
    COFLAM(64) = -1.45368336d-02
    COFLAM(65) = -1.46152839d+01
    COFLAM(66) = 6.36251406d+00
    COFLAM(67) = -5.03832130d-01
    COFLAM(68) = 1.26121050d-02
    COFLAM(69) = -8.95009705d+00
    COFLAM(70) = 4.02515080d+00
    COFLAM(71) = -1.84063946d-01
    COFLAM(72) = -1.94054752d-03
    COFLAM(73) = -1.09902209d+01
    COFLAM(74) = 4.70647707d+00
    COFLAM(75) = -2.52272495d-01
    COFLAM(76) = 1.75193258d-04
    COFLAM(77) = 1.29306158d+01
    COFLAM(78) = -3.52817362d+00
    COFLAM(79) = 6.45499013d-01
    COFLAM(80) = -3.19375299d-02
    COFLAM(81) = -3.17202048d+00
    COFLAM(82) = 3.47025711d+00
    COFLAM(83) = -3.75102111d-01
    COFLAM(84) = 1.66086076d-02

end subroutine

! Poly fits for the diffusion coefficients, dim NO*KK*KK
subroutine egtransetCOFD(COFD)

    implicit none

    double precision, intent(out) :: COFD(1764)

    COFD(1) = -1.03270606d+01
    COFD(2) = 2.19285409d+00
    COFD(3) = -7.54492786d-02
    COFD(4) = 3.51398213d-03
    COFD(5) = -1.14366381d+01
    COFD(6) = 2.78323501d+00
    COFD(7) = -1.51214064d-01
    COFD(8) = 6.75150012d-03
    COFD(9) = -1.09595712d+01
    COFD(10) = 2.30836460d+00
    COFD(11) = -8.76339315d-02
    COFD(12) = 3.90878445d-03
    COFD(13) = -1.18988955d+01
    COFD(14) = 2.57507000d+00
    COFD(15) = -1.24033737d-01
    COFD(16) = 5.56694959d-03
    COFD(17) = -1.09628982d+01
    COFD(18) = 2.30836460d+00
    COFD(19) = -8.76339315d-02
    COFD(20) = 3.90878445d-03
    COFD(21) = -1.71982995d+01
    COFD(22) = 4.63881404d+00
    COFD(23) = -3.86139633d-01
    COFD(24) = 1.66955081d-02
    COFD(25) = -1.18998012d+01
    COFD(26) = 2.57507000d+00
    COFD(27) = -1.24033737d-01
    COFD(28) = 5.56694959d-03
    COFD(29) = -1.25098960d+01
    COFD(30) = 2.77873601d+00
    COFD(31) = -1.50637360d-01
    COFD(32) = 6.72684281d-03
    COFD(33) = -1.25098960d+01
    COFD(34) = 2.77873601d+00
    COFD(35) = -1.50637360d-01
    COFD(36) = 6.72684281d-03
    COFD(37) = -1.25141260d+01
    COFD(38) = 2.77873601d+00
    COFD(39) = -1.50637360d-01
    COFD(40) = 6.72684281d-03
    COFD(41) = -1.24693568d+01
    COFD(42) = 2.76686648d+00
    COFD(43) = -1.49120141d-01
    COFD(44) = 6.66220432d-03
    COFD(45) = -1.17159737d+01
    COFD(46) = 2.48123210d+00
    COFD(47) = -1.11322604d-01
    COFD(48) = 4.99282389d-03
    COFD(49) = -1.37794315d+01
    COFD(50) = 3.23973858d+00
    COFD(51) = -2.09989036d-01
    COFD(52) = 9.27667906d-03
    COFD(53) = -1.60517370d+01
    COFD(54) = 4.11188603d+00
    COFD(55) = -3.21540884d-01
    COFD(56) = 1.40482564d-02
    COFD(57) = -1.60528285d+01
    COFD(58) = 4.11188603d+00
    COFD(59) = -3.21540884d-01
    COFD(60) = 1.40482564d-02
    COFD(61) = -1.58456300d+01
    COFD(62) = 4.02074783d+00
    COFD(63) = -3.10018522d-01
    COFD(64) = 1.35599552d-02
    COFD(65) = -1.42229194d+01
    COFD(66) = 3.38669384d+00
    COFD(67) = -2.28784122d-01
    COFD(68) = 1.00790953d-02
    COFD(69) = -1.39913897d+01
    COFD(70) = 3.26384506d+00
    COFD(71) = -2.12947087d-01
    COFD(72) = 9.39743888d-03
    COFD(73) = -1.39924781d+01
    COFD(74) = 3.26384506d+00
    COFD(75) = -2.12947087d-01
    COFD(76) = 9.39743888d-03
    COFD(77) = -1.16906297d+01
    COFD(78) = 2.47469981d+00
    COFD(79) = -1.10436257d-01
    COFD(80) = 4.95273813d-03
    COFD(81) = -1.23130152d+01
    COFD(82) = 2.74418790d+00
    COFD(83) = -1.46230156d-01
    COFD(84) = 6.53948886d-03
    COFD(85) = -1.14366381d+01
    COFD(86) = 2.78323501d+00
    COFD(87) = -1.51214064d-01
    COFD(88) = 6.75150012d-03
    COFD(89) = -1.47968712d+01
    COFD(90) = 4.23027636d+00
    COFD(91) = -3.36139991d-01
    COFD(92) = 1.46507621d-02
    COFD(93) = -1.34230272d+01
    COFD(94) = 3.48624238d+00
    COFD(95) = -2.41554467d-01
    COFD(96) = 1.06263545d-02
    COFD(97) = -1.46550083d+01
    COFD(98) = 3.83606243d+00
    COFD(99) = -2.86076532d-01
    COFD(100) = 1.25205829d-02
    COFD(101) = -1.34247866d+01
    COFD(102) = 3.48624238d+00
    COFD(103) = -2.41554467d-01
    COFD(104) = 1.06263545d-02
    COFD(105) = -1.95739570d+01
    COFD(106) = 5.61113230d+00
    COFD(107) = -4.90190187d-01
    COFD(108) = 2.03260675d-02
    COFD(109) = -1.46554748d+01
    COFD(110) = 3.83606243d+00
    COFD(111) = -2.86076532d-01
    COFD(112) = 1.25205829d-02
    COFD(113) = -1.57972369d+01
    COFD(114) = 4.22225052d+00
    COFD(115) = -3.35156428d-01
    COFD(116) = 1.46104855d-02
    COFD(117) = -1.57972369d+01
    COFD(118) = 4.22225052d+00
    COFD(119) = -3.35156428d-01
    COFD(120) = 1.46104855d-02
    COFD(121) = -1.57994893d+01
    COFD(122) = 4.22225052d+00
    COFD(123) = -3.35156428d-01
    COFD(124) = 1.46104855d-02
    COFD(125) = -1.57199037d+01
    COFD(126) = 4.19936335d+00
    COFD(127) = -3.32311009d-01
    COFD(128) = 1.44921003d-02
    COFD(129) = -1.43151174d+01
    COFD(130) = 3.68038508d+00
    COFD(131) = -2.65779346d-01
    COFD(132) = 1.16360771d-02
    COFD(133) = -1.76147026d+01
    COFD(134) = 4.86049500d+00
    COFD(135) = -4.12200578d-01
    COFD(136) = 1.77160971d-02
    COFD(137) = -1.97544450d+01
    COFD(138) = 5.56931926d+00
    COFD(139) = -4.89105511d-01
    COFD(140) = 2.04493129d-02
    COFD(141) = -1.97550088d+01
    COFD(142) = 5.56931926d+00
    COFD(143) = -4.89105511d-01
    COFD(144) = 2.04493129d-02
    COFD(145) = -1.92718582d+01
    COFD(146) = 5.41172124d+00
    COFD(147) = -4.73213887d-01
    COFD(148) = 1.99405473d-02
    COFD(149) = -1.82251914d+01
    COFD(150) = 5.05237312d+00
    COFD(151) = -4.35182396d-01
    COFD(152) = 1.86363074d-02
    COFD(153) = -1.79339327d+01
    COFD(154) = 4.91373893d+00
    COFD(155) = -4.18747629d-01
    COFD(156) = 1.79856610d-02
    COFD(157) = -1.79344949d+01
    COFD(158) = 4.91373893d+00
    COFD(159) = -4.18747629d-01
    COFD(160) = 1.79856610d-02
    COFD(161) = -1.42894441d+01
    COFD(162) = 3.67490723d+00
    COFD(163) = -2.65114792d-01
    COFD(164) = 1.16092671d-02
    COFD(165) = -1.54738604d+01
    COFD(166) = 4.15765300d+00
    COFD(167) = -3.27126237d-01
    COFD(168) = 1.42762611d-02
    COFD(169) = -1.09595712d+01
    COFD(170) = 2.30836460d+00
    COFD(171) = -8.76339315d-02
    COFD(172) = 3.90878445d-03
    COFD(173) = -1.34230272d+01
    COFD(174) = 3.48624238d+00
    COFD(175) = -2.41554467d-01
    COFD(176) = 1.06263545d-02
    COFD(177) = -1.32093628d+01
    COFD(178) = 2.90778936d+00
    COFD(179) = -1.67388544d-01
    COFD(180) = 7.45220609d-03
    COFD(181) = -1.43139231d+01
    COFD(182) = 3.17651319d+00
    COFD(183) = -2.02028974d-01
    COFD(184) = 8.94232502d-03
    COFD(185) = -1.32244035d+01
    COFD(186) = 2.90778936d+00
    COFD(187) = -1.67388544d-01
    COFD(188) = 7.45220609d-03
    COFD(189) = -1.94093572d+01
    COFD(190) = 5.16013126d+00
    COFD(191) = -4.46824543d-01
    COFD(192) = 1.90464887d-02
    COFD(193) = -1.43190389d+01
    COFD(194) = 3.17651319d+00
    COFD(195) = -2.02028974d-01
    COFD(196) = 8.94232502d-03
    COFD(197) = -1.50584249d+01
    COFD(198) = 3.47945612d+00
    COFD(199) = -2.40703722d-01
    COFD(200) = 1.05907441d-02
    COFD(201) = -1.50584249d+01
    COFD(202) = 3.47945612d+00
    COFD(203) = -2.40703722d-01
    COFD(204) = 1.05907441d-02
    COFD(205) = -1.50766130d+01
    COFD(206) = 3.47945612d+00
    COFD(207) = -2.40703722d-01
    COFD(208) = 1.05907441d-02
    COFD(209) = -1.50270339d+01
    COFD(210) = 3.46140064d+00
    COFD(211) = -2.38440092d-01
    COFD(212) = 1.04960087d-02
    COFD(213) = -1.40999008d+01
    COFD(214) = 3.08120012d+00
    COFD(215) = -1.89629903d-01
    COFD(216) = 8.40361952d-03
    COFD(217) = -1.70534856d+01
    COFD(218) = 4.14240922d+00
    COFD(219) = -3.25239774d-01
    COFD(220) = 1.41980687d-02
    COFD(221) = -1.94313116d+01
    COFD(222) = 5.02567894d+00
    COFD(223) = -4.32045169d-01
    COFD(224) = 1.85132214d-02
    COFD(225) = -1.94373127d+01
    COFD(226) = 5.02567894d+00
    COFD(227) = -4.32045169d-01
    COFD(228) = 1.85132214d-02
    COFD(229) = -1.88179418d+01
    COFD(230) = 4.79683898d+00
    COFD(231) = -4.04829719d-01
    COFD(232) = 1.74325475d-02
    COFD(233) = -1.74792112d+01
    COFD(234) = 4.29676909d+00
    COFD(235) = -3.44085306d-01
    COFD(236) = 1.49671135d-02
    COFD(237) = -1.72496634d+01
    COFD(238) = 4.17889917d+00
    COFD(239) = -3.29752510d-01
    COFD(240) = 1.43850275d-02
    COFD(241) = -1.72556499d+01
    COFD(242) = 4.17889917d+00
    COFD(243) = -3.29752510d-01
    COFD(244) = 1.43850275d-02
    COFD(245) = -1.40756935d+01
    COFD(246) = 3.07549274d+00
    COFD(247) = -1.88889344d-01
    COFD(248) = 8.37152866d-03
    COFD(249) = -1.49610527d+01
    COFD(250) = 3.41988961d+00
    COFD(251) = -2.33128386d-01
    COFD(252) = 1.02689994d-02
    COFD(253) = -1.18988955d+01
    COFD(254) = 2.57507000d+00
    COFD(255) = -1.24033737d-01
    COFD(256) = 5.56694959d-03
    COFD(257) = -1.46550083d+01
    COFD(258) = 3.83606243d+00
    COFD(259) = -2.86076532d-01
    COFD(260) = 1.25205829d-02
    COFD(261) = -1.43139231d+01
    COFD(262) = 3.17651319d+00
    COFD(263) = -2.02028974d-01
    COFD(264) = 8.94232502d-03
    COFD(265) = -1.55511344d+01
    COFD(266) = 3.48070094d+00
    COFD(267) = -2.40859499d-01
    COFD(268) = 1.05972514d-02
    COFD(269) = -1.43340796d+01
    COFD(270) = 3.17651319d+00
    COFD(271) = -2.02028974d-01
    COFD(272) = 8.94232502d-03
    COFD(273) = -2.12652533d+01
    COFD(274) = 5.59961818d+00
    COFD(275) = -4.91624858d-01
    COFD(276) = 2.05035550d-02
    COFD(277) = -1.55588279d+01
    COFD(278) = 3.48070094d+00
    COFD(279) = -2.40859499d-01
    COFD(280) = 1.05972514d-02
    COFD(281) = -1.63254691d+01
    COFD(282) = 3.82388595d+00
    COFD(283) = -2.84480724d-01
    COFD(284) = 1.24506311d-02
    COFD(285) = -1.63254691d+01
    COFD(286) = 3.82388595d+00
    COFD(287) = -2.84480724d-01
    COFD(288) = 1.24506311d-02
    COFD(289) = -1.63493345d+01
    COFD(290) = 3.82388595d+00
    COFD(291) = -2.84480724d-01
    COFD(292) = 1.24506311d-02
    COFD(293) = -1.62724462d+01
    COFD(294) = 3.79163564d+00
    COFD(295) = -2.80257365d-01
    COFD(296) = 1.22656902d-02
    COFD(297) = -1.52721107d+01
    COFD(298) = 3.36790500d+00
    COFD(299) = -2.26321740d-01
    COFD(300) = 9.97135055d-03
    COFD(301) = -1.84688406d+01
    COFD(302) = 4.49330851d+00
    COFD(303) = -3.68208715d-01
    COFD(304) = 1.59565402d-02
    COFD(305) = -2.08204449d+01
    COFD(306) = 5.35267674d+00
    COFD(307) = -4.69010505d-01
    COFD(308) = 1.98979152d-02
    COFD(309) = -2.08293255d+01
    COFD(310) = 5.35267674d+00
    COFD(311) = -4.69010505d-01
    COFD(312) = 1.98979152d-02
    COFD(313) = -2.04928958d+01
    COFD(314) = 5.22397933d+00
    COFD(315) = -4.54138171d-01
    COFD(316) = 1.93249285d-02
    COFD(317) = -1.89544778d+01
    COFD(318) = 4.68595732d+00
    COFD(319) = -3.91842840d-01
    COFD(320) = 1.69262542d-02
    COFD(321) = -1.86335932d+01
    COFD(322) = 4.53572533d+00
    COFD(323) = -3.73386925d-01
    COFD(324) = 1.61678881d-02
    COFD(325) = -1.86424545d+01
    COFD(326) = 4.53572533d+00
    COFD(327) = -3.73386925d-01
    COFD(328) = 1.61678881d-02
    COFD(329) = -1.52414485d+01
    COFD(330) = 3.35922578d+00
    COFD(331) = -2.25181399d-01
    COFD(332) = 9.92132878d-03
    COFD(333) = -1.62380075d+01
    COFD(334) = 3.72612300d+00
    COFD(335) = -2.71663673d-01
    COFD(336) = 1.18889643d-02
    COFD(337) = -1.09628982d+01
    COFD(338) = 2.30836460d+00
    COFD(339) = -8.76339315d-02
    COFD(340) = 3.90878445d-03
    COFD(341) = -1.34247866d+01
    COFD(342) = 3.48624238d+00
    COFD(343) = -2.41554467d-01
    COFD(344) = 1.06263545d-02
    COFD(345) = -1.32244035d+01
    COFD(346) = 2.90778936d+00
    COFD(347) = -1.67388544d-01
    COFD(348) = 7.45220609d-03
    COFD(349) = -1.43340796d+01
    COFD(350) = 3.17651319d+00
    COFD(351) = -2.02028974d-01
    COFD(352) = 8.94232502d-03
    COFD(353) = -1.32399106d+01
    COFD(354) = 2.90778936d+00
    COFD(355) = -1.67388544d-01
    COFD(356) = 7.45220609d-03
    COFD(357) = -1.94253036d+01
    COFD(358) = 5.16013126d+00
    COFD(359) = -4.46824543d-01
    COFD(360) = 1.90464887d-02
    COFD(361) = -1.43394069d+01
    COFD(362) = 3.17651319d+00
    COFD(363) = -2.02028974d-01
    COFD(364) = 8.94232502d-03
    COFD(365) = -1.50724636d+01
    COFD(366) = 3.47945612d+00
    COFD(367) = -2.40703722d-01
    COFD(368) = 1.05907441d-02
    COFD(369) = -1.50724636d+01
    COFD(370) = 3.47945612d+00
    COFD(371) = -2.40703722d-01
    COFD(372) = 1.05907441d-02
    COFD(373) = -1.50911794d+01
    COFD(374) = 3.47945612d+00
    COFD(375) = -2.40703722d-01
    COFD(376) = 1.05907441d-02
    COFD(377) = -1.50420953d+01
    COFD(378) = 3.46140064d+00
    COFD(379) = -2.38440092d-01
    COFD(380) = 1.04960087d-02
    COFD(381) = -1.41191261d+01
    COFD(382) = 3.08120012d+00
    COFD(383) = -1.89629903d-01
    COFD(384) = 8.40361952d-03
    COFD(385) = -1.70757047d+01
    COFD(386) = 4.14240922d+00
    COFD(387) = -3.25239774d-01
    COFD(388) = 1.41980687d-02
    COFD(389) = -1.94507876d+01
    COFD(390) = 5.02567894d+00
    COFD(391) = -4.32045169d-01
    COFD(392) = 1.85132214d-02
    COFD(393) = -1.94570287d+01
    COFD(394) = 5.02567894d+00
    COFD(395) = -4.32045169d-01
    COFD(396) = 1.85132214d-02
    COFD(397) = -1.88378874d+01
    COFD(398) = 4.79683898d+00
    COFD(399) = -4.04829719d-01
    COFD(400) = 1.74325475d-02
    COFD(401) = -1.74984476d+01
    COFD(402) = 4.29676909d+00
    COFD(403) = -3.44085306d-01
    COFD(404) = 1.49671135d-02
    COFD(405) = -1.72691500d+01
    COFD(406) = 4.17889917d+00
    COFD(407) = -3.29752510d-01
    COFD(408) = 1.43850275d-02
    COFD(409) = -1.72753760d+01
    COFD(410) = 4.17889917d+00
    COFD(411) = -3.29752510d-01
    COFD(412) = 1.43850275d-02
    COFD(413) = -1.40949196d+01
    COFD(414) = 3.07549274d+00
    COFD(415) = -1.88889344d-01
    COFD(416) = 8.37152866d-03
    COFD(417) = -1.49826725d+01
    COFD(418) = 3.41988961d+00
    COFD(419) = -2.33128386d-01
    COFD(420) = 1.02689994d-02
    COFD(421) = -1.71982995d+01
    COFD(422) = 4.63881404d+00
    COFD(423) = -3.86139633d-01
    COFD(424) = 1.66955081d-02
    COFD(425) = -1.95739570d+01
    COFD(426) = 5.61113230d+00
    COFD(427) = -4.90190187d-01
    COFD(428) = 2.03260675d-02
    COFD(429) = -1.94093572d+01
    COFD(430) = 5.16013126d+00
    COFD(431) = -4.46824543d-01
    COFD(432) = 1.90464887d-02
    COFD(433) = -2.12652533d+01
    COFD(434) = 5.59961818d+00
    COFD(435) = -4.91624858d-01
    COFD(436) = 2.05035550d-02
    COFD(437) = -1.94253036d+01
    COFD(438) = 5.16013126d+00
    COFD(439) = -4.46824543d-01
    COFD(440) = 1.90464887d-02
    COFD(441) = -1.19157919d+01
    COFD(442) = 9.28955130d-01
    COFD(443) = 2.42107090d-01
    COFD(444) = -1.59823963d-02
    COFD(445) = -2.06463744d+01
    COFD(446) = 5.41688482d+00
    COFD(447) = -4.73387188d-01
    COFD(448) = 1.99280175d-02
    COFD(449) = -2.12639214d+01
    COFD(450) = 5.61184117d+00
    COFD(451) = -4.90532156d-01
    COFD(452) = 2.03507922d-02
    COFD(453) = -2.12639214d+01
    COFD(454) = 5.61184117d+00
    COFD(455) = -4.90532156d-01
    COFD(456) = 2.03507922d-02
    COFD(457) = -2.12831323d+01
    COFD(458) = 5.61184117d+00
    COFD(459) = -4.90532156d-01
    COFD(460) = 2.03507922d-02
    COFD(461) = -2.14087397d+01
    COFD(462) = 5.57282008d+00
    COFD(463) = -4.76690890d-01
    COFD(464) = 1.94000719d-02
    COFD(465) = -2.11388331d+01
    COFD(466) = 5.55529675d+00
    COFD(467) = -4.87942518d-01
    COFD(468) = 2.04249054d-02
    COFD(469) = -2.07653719d+01
    COFD(470) = 5.01092022d+00
    COFD(471) = -3.77985635d-01
    COFD(472) = 1.40968645d-02
    COFD(473) = -1.77498543d+01
    COFD(474) = 3.57475686d+00
    COFD(475) = -1.56396297d-01
    COFD(476) = 3.12157721d-03
    COFD(477) = -1.77563250d+01
    COFD(478) = 3.57475686d+00
    COFD(479) = -1.56396297d-01
    COFD(480) = 3.12157721d-03
    COFD(481) = -1.65295288d+01
    COFD(482) = 2.97569206d+00
    COFD(483) = -6.75652842d-02
    COFD(484) = -1.08648422d-03
    COFD(485) = -2.08812333d+01
    COFD(486) = 5.08859217d+00
    COFD(487) = -3.90525428d-01
    COFD(488) = 1.47376395d-02
    COFD(489) = -2.12597312d+01
    COFD(490) = 5.24930667d+00
    COFD(491) = -4.17435088d-01
    COFD(492) = 1.61434424d-02
    COFD(493) = -2.12661865d+01
    COFD(494) = 5.24930667d+00
    COFD(495) = -4.17435088d-01
    COFD(496) = 1.61434424d-02
    COFD(497) = -2.10643259d+01
    COFD(498) = 5.53614847d+00
    COFD(499) = -4.86046736d-01
    COFD(500) = 2.03659188d-02
    COFD(501) = -2.12755888d+01
    COFD(502) = 5.60381989d+00
    COFD(503) = -4.91225459d-01
    COFD(504) = 2.04487844d-02
    COFD(505) = -1.18998012d+01
    COFD(506) = 2.57507000d+00
    COFD(507) = -1.24033737d-01
    COFD(508) = 5.56694959d-03
    COFD(509) = -1.46554748d+01
    COFD(510) = 3.83606243d+00
    COFD(511) = -2.86076532d-01
    COFD(512) = 1.25205829d-02
    COFD(513) = -1.43190389d+01
    COFD(514) = 3.17651319d+00
    COFD(515) = -2.02028974d-01
    COFD(516) = 8.94232502d-03
    COFD(517) = -1.55588279d+01
    COFD(518) = 3.48070094d+00
    COFD(519) = -2.40859499d-01
    COFD(520) = 1.05972514d-02
    COFD(521) = -1.43394069d+01
    COFD(522) = 3.17651319d+00
    COFD(523) = -2.02028974d-01
    COFD(524) = 8.94232502d-03
    COFD(525) = -2.06463744d+01
    COFD(526) = 5.41688482d+00
    COFD(527) = -4.73387188d-01
    COFD(528) = 1.99280175d-02
    COFD(529) = -1.55666415d+01
    COFD(530) = 3.48070094d+00
    COFD(531) = -2.40859499d-01
    COFD(532) = 1.05972514d-02
    COFD(533) = -1.63301444d+01
    COFD(534) = 3.82388595d+00
    COFD(535) = -2.84480724d-01
    COFD(536) = 1.24506311d-02
    COFD(537) = -1.63301444d+01
    COFD(538) = 3.82388595d+00
    COFD(539) = -2.84480724d-01
    COFD(540) = 1.24506311d-02
    COFD(541) = -1.63542394d+01
    COFD(542) = 3.82388595d+00
    COFD(543) = -2.84480724d-01
    COFD(544) = 1.24506311d-02
    COFD(545) = -1.62775714d+01
    COFD(546) = 3.79163564d+00
    COFD(547) = -2.80257365d-01
    COFD(548) = 1.22656902d-02
    COFD(549) = -1.52792891d+01
    COFD(550) = 3.36790500d+00
    COFD(551) = -2.26321740d-01
    COFD(552) = 9.97135055d-03
    COFD(553) = -1.84777607d+01
    COFD(554) = 4.49330851d+00
    COFD(555) = -3.68208715d-01
    COFD(556) = 1.59565402d-02
    COFD(557) = -2.08277598d+01
    COFD(558) = 5.35267674d+00
    COFD(559) = -4.69010505d-01
    COFD(560) = 1.98979152d-02
    COFD(561) = -2.08367725d+01
    COFD(562) = 5.35267674d+00
    COFD(563) = -4.69010505d-01
    COFD(564) = 1.98979152d-02
    COFD(565) = -2.02637994d+01
    COFD(566) = 5.14984081d+00
    COFD(567) = -4.46093018d-01
    COFD(568) = 1.90396647d-02
    COFD(569) = -1.89616623d+01
    COFD(570) = 4.68595732d+00
    COFD(571) = -3.91842840d-01
    COFD(572) = 1.69262542d-02
    COFD(573) = -1.86409139d+01
    COFD(574) = 4.53572533d+00
    COFD(575) = -3.73386925d-01
    COFD(576) = 1.61678881d-02
    COFD(577) = -1.86499071d+01
    COFD(578) = 4.53572533d+00
    COFD(579) = -3.73386925d-01
    COFD(580) = 1.61678881d-02
    COFD(581) = -1.52486273d+01
    COFD(582) = 3.35922578d+00
    COFD(583) = -2.25181399d-01
    COFD(584) = 9.92132878d-03
    COFD(585) = -1.62465583d+01
    COFD(586) = 3.72612300d+00
    COFD(587) = -2.71663673d-01
    COFD(588) = 1.18889643d-02
    COFD(589) = -1.25098960d+01
    COFD(590) = 2.77873601d+00
    COFD(591) = -1.50637360d-01
    COFD(592) = 6.72684281d-03
    COFD(593) = -1.57972369d+01
    COFD(594) = 4.22225052d+00
    COFD(595) = -3.35156428d-01
    COFD(596) = 1.46104855d-02
    COFD(597) = -1.50584249d+01
    COFD(598) = 3.47945612d+00
    COFD(599) = -2.40703722d-01
    COFD(600) = 1.05907441d-02
    COFD(601) = -1.63254691d+01
    COFD(602) = 3.82388595d+00
    COFD(603) = -2.84480724d-01
    COFD(604) = 1.24506311d-02
    COFD(605) = -1.50724636d+01
    COFD(606) = 3.47945612d+00
    COFD(607) = -2.40703722d-01
    COFD(608) = 1.05907441d-02
    COFD(609) = -2.12639214d+01
    COFD(610) = 5.61184117d+00
    COFD(611) = -4.90532156d-01
    COFD(612) = 2.03507922d-02
    COFD(613) = -1.63301444d+01
    COFD(614) = 3.82388595d+00
    COFD(615) = -2.84480724d-01
    COFD(616) = 1.24506311d-02
    COFD(617) = -1.73027557d+01
    COFD(618) = 4.21416723d+00
    COFD(619) = -3.34163932d-01
    COFD(620) = 1.45697432d-02
    COFD(621) = -1.73027557d+01
    COFD(622) = 4.21416723d+00
    COFD(623) = -3.34163932d-01
    COFD(624) = 1.45697432d-02
    COFD(625) = -1.73198034d+01
    COFD(626) = 4.21416723d+00
    COFD(627) = -3.34163932d-01
    COFD(628) = 1.45697432d-02
    COFD(629) = -1.72556729d+01
    COFD(630) = 4.19029808d+00
    COFD(631) = -3.31177076d-01
    COFD(632) = 1.44446234d-02
    COFD(633) = -1.59634533d+01
    COFD(634) = 3.67388294d+00
    COFD(635) = -2.64990709d-01
    COFD(636) = 1.16042706d-02
    COFD(637) = -1.93015555d+01
    COFD(638) = 4.85015581d+00
    COFD(639) = -4.10945109d-01
    COFD(640) = 1.76651398d-02
    COFD(641) = -2.14160703d+01
    COFD(642) = 5.56531152d+00
    COFD(643) = -4.88789821d-01
    COFD(644) = 2.04437116d-02
    COFD(645) = -2.14215700d+01
    COFD(646) = 5.56531152d+00
    COFD(647) = -4.88789821d-01
    COFD(648) = 2.04437116d-02
    COFD(649) = -2.09376196d+01
    COFD(650) = 5.40870099d+00
    COFD(651) = -4.73017610d-01
    COFD(652) = 1.99399066d-02
    COFD(653) = -1.98418115d+01
    COFD(654) = 5.04367502d+00
    COFD(655) = -4.34153325d-01
    COFD(656) = 1.85956055d-02
    COFD(657) = -1.95263312d+01
    COFD(658) = 4.90255048d+00
    COFD(659) = -4.17368501d-01
    COFD(660) = 1.79287358d-02
    COFD(661) = -1.95318173d+01
    COFD(662) = 4.90255048d+00
    COFD(663) = -4.17368501d-01
    COFD(664) = 1.79287358d-02
    COFD(665) = -1.59404882d+01
    COFD(666) = 3.66853818d+00
    COFD(667) = -2.64346221d-01
    COFD(668) = 1.15784613d-02
    COFD(669) = -1.71942502d+01
    COFD(670) = 4.14993355d+00
    COFD(671) = -3.26168062d-01
    COFD(672) = 1.42364115d-02
    COFD(673) = -1.25098960d+01
    COFD(674) = 2.77873601d+00
    COFD(675) = -1.50637360d-01
    COFD(676) = 6.72684281d-03
    COFD(677) = -1.57972369d+01
    COFD(678) = 4.22225052d+00
    COFD(679) = -3.35156428d-01
    COFD(680) = 1.46104855d-02
    COFD(681) = -1.50584249d+01
    COFD(682) = 3.47945612d+00
    COFD(683) = -2.40703722d-01
    COFD(684) = 1.05907441d-02
    COFD(685) = -1.63254691d+01
    COFD(686) = 3.82388595d+00
    COFD(687) = -2.84480724d-01
    COFD(688) = 1.24506311d-02
    COFD(689) = -1.50724636d+01
    COFD(690) = 3.47945612d+00
    COFD(691) = -2.40703722d-01
    COFD(692) = 1.05907441d-02
    COFD(693) = -2.12639214d+01
    COFD(694) = 5.61184117d+00
    COFD(695) = -4.90532156d-01
    COFD(696) = 2.03507922d-02
    COFD(697) = -1.63301444d+01
    COFD(698) = 3.82388595d+00
    COFD(699) = -2.84480724d-01
    COFD(700) = 1.24506311d-02
    COFD(701) = -1.73027557d+01
    COFD(702) = 4.21416723d+00
    COFD(703) = -3.34163932d-01
    COFD(704) = 1.45697432d-02
    COFD(705) = -1.73027557d+01
    COFD(706) = 4.21416723d+00
    COFD(707) = -3.34163932d-01
    COFD(708) = 1.45697432d-02
    COFD(709) = -1.73198034d+01
    COFD(710) = 4.21416723d+00
    COFD(711) = -3.34163932d-01
    COFD(712) = 1.45697432d-02
    COFD(713) = -1.72556729d+01
    COFD(714) = 4.19029808d+00
    COFD(715) = -3.31177076d-01
    COFD(716) = 1.44446234d-02
    COFD(717) = -1.59634533d+01
    COFD(718) = 3.67388294d+00
    COFD(719) = -2.64990709d-01
    COFD(720) = 1.16042706d-02
    COFD(721) = -1.93015555d+01
    COFD(722) = 4.85015581d+00
    COFD(723) = -4.10945109d-01
    COFD(724) = 1.76651398d-02
    COFD(725) = -2.14160703d+01
    COFD(726) = 5.56531152d+00
    COFD(727) = -4.88789821d-01
    COFD(728) = 2.04437116d-02
    COFD(729) = -2.14215700d+01
    COFD(730) = 5.56531152d+00
    COFD(731) = -4.88789821d-01
    COFD(732) = 2.04437116d-02
    COFD(733) = -2.09376196d+01
    COFD(734) = 5.40870099d+00
    COFD(735) = -4.73017610d-01
    COFD(736) = 1.99399066d-02
    COFD(737) = -1.98418115d+01
    COFD(738) = 5.04367502d+00
    COFD(739) = -4.34153325d-01
    COFD(740) = 1.85956055d-02
    COFD(741) = -1.95263312d+01
    COFD(742) = 4.90255048d+00
    COFD(743) = -4.17368501d-01
    COFD(744) = 1.79287358d-02
    COFD(745) = -1.95318173d+01
    COFD(746) = 4.90255048d+00
    COFD(747) = -4.17368501d-01
    COFD(748) = 1.79287358d-02
    COFD(749) = -1.59404882d+01
    COFD(750) = 3.66853818d+00
    COFD(751) = -2.64346221d-01
    COFD(752) = 1.15784613d-02
    COFD(753) = -1.71942502d+01
    COFD(754) = 4.14993355d+00
    COFD(755) = -3.26168062d-01
    COFD(756) = 1.42364115d-02
    COFD(757) = -1.25141260d+01
    COFD(758) = 2.77873601d+00
    COFD(759) = -1.50637360d-01
    COFD(760) = 6.72684281d-03
    COFD(761) = -1.57994893d+01
    COFD(762) = 4.22225052d+00
    COFD(763) = -3.35156428d-01
    COFD(764) = 1.46104855d-02
    COFD(765) = -1.50766130d+01
    COFD(766) = 3.47945612d+00
    COFD(767) = -2.40703722d-01
    COFD(768) = 1.05907441d-02
    COFD(769) = -1.63493345d+01
    COFD(770) = 3.82388595d+00
    COFD(771) = -2.84480724d-01
    COFD(772) = 1.24506311d-02
    COFD(773) = -1.50911794d+01
    COFD(774) = 3.47945612d+00
    COFD(775) = -2.40703722d-01
    COFD(776) = 1.05907441d-02
    COFD(777) = -2.12831323d+01
    COFD(778) = 5.61184117d+00
    COFD(779) = -4.90532156d-01
    COFD(780) = 2.03507922d-02
    COFD(781) = -1.63542394d+01
    COFD(782) = 3.82388595d+00
    COFD(783) = -2.84480724d-01
    COFD(784) = 1.24506311d-02
    COFD(785) = -1.73198034d+01
    COFD(786) = 4.21416723d+00
    COFD(787) = -3.34163932d-01
    COFD(788) = 1.45697432d-02
    COFD(789) = -1.73198034d+01
    COFD(790) = 4.21416723d+00
    COFD(791) = -3.34163932d-01
    COFD(792) = 1.45697432d-02
    COFD(793) = -1.73374529d+01
    COFD(794) = 4.21416723d+00
    COFD(795) = -3.34163932d-01
    COFD(796) = 1.45697432d-02
    COFD(797) = -1.72738845d+01
    COFD(798) = 4.19029808d+00
    COFD(799) = -3.31177076d-01
    COFD(800) = 1.44446234d-02
    COFD(801) = -1.59863030d+01
    COFD(802) = 3.67388294d+00
    COFD(803) = -2.64990709d-01
    COFD(804) = 1.16042706d-02
    COFD(805) = -1.93276434d+01
    COFD(806) = 4.85015581d+00
    COFD(807) = -4.10945109d-01
    COFD(808) = 1.76651398d-02
    COFD(809) = -2.14391943d+01
    COFD(810) = 5.56531152d+00
    COFD(811) = -4.88789821d-01
    COFD(812) = 2.04437116d-02
    COFD(813) = -2.14449559d+01
    COFD(814) = 5.56531152d+00
    COFD(815) = -4.88789821d-01
    COFD(816) = 2.04437116d-02
    COFD(817) = -2.09612557d+01
    COFD(818) = 5.40870099d+00
    COFD(819) = -4.73017610d-01
    COFD(820) = 1.99399066d-02
    COFD(821) = -1.98646734d+01
    COFD(822) = 5.04367502d+00
    COFD(823) = -4.34153325d-01
    COFD(824) = 1.85956055d-02
    COFD(825) = -1.95494668d+01
    COFD(826) = 4.90255048d+00
    COFD(827) = -4.17368501d-01
    COFD(828) = 1.79287358d-02
    COFD(829) = -1.95552142d+01
    COFD(830) = 4.90255048d+00
    COFD(831) = -4.17368501d-01
    COFD(832) = 1.79287358d-02
    COFD(833) = -1.59633387d+01
    COFD(834) = 3.66853818d+00
    COFD(835) = -2.64346221d-01
    COFD(836) = 1.15784613d-02
    COFD(837) = -1.72196961d+01
    COFD(838) = 4.14993355d+00
    COFD(839) = -3.26168062d-01
    COFD(840) = 1.42364115d-02
    COFD(841) = -1.24693568d+01
    COFD(842) = 2.76686648d+00
    COFD(843) = -1.49120141d-01
    COFD(844) = 6.66220432d-03
    COFD(845) = -1.57199037d+01
    COFD(846) = 4.19936335d+00
    COFD(847) = -3.32311009d-01
    COFD(848) = 1.44921003d-02
    COFD(849) = -1.50270339d+01
    COFD(850) = 3.46140064d+00
    COFD(851) = -2.38440092d-01
    COFD(852) = 1.04960087d-02
    COFD(853) = -1.62724462d+01
    COFD(854) = 3.79163564d+00
    COFD(855) = -2.80257365d-01
    COFD(856) = 1.22656902d-02
    COFD(857) = -1.50420953d+01
    COFD(858) = 3.46140064d+00
    COFD(859) = -2.38440092d-01
    COFD(860) = 1.04960087d-02
    COFD(861) = -2.14087397d+01
    COFD(862) = 5.57282008d+00
    COFD(863) = -4.76690890d-01
    COFD(864) = 1.94000719d-02
    COFD(865) = -1.62775714d+01
    COFD(866) = 3.79163564d+00
    COFD(867) = -2.80257365d-01
    COFD(868) = 1.22656902d-02
    COFD(869) = -1.72556729d+01
    COFD(870) = 4.19029808d+00
    COFD(871) = -3.31177076d-01
    COFD(872) = 1.44446234d-02
    COFD(873) = -1.72556729d+01
    COFD(874) = 4.19029808d+00
    COFD(875) = -3.31177076d-01
    COFD(876) = 1.44446234d-02
    COFD(877) = -1.72738845d+01
    COFD(878) = 4.19029808d+00
    COFD(879) = -3.31177076d-01
    COFD(880) = 1.44446234d-02
    COFD(881) = -1.72167708d+01
    COFD(882) = 4.16886779d+00
    COFD(883) = -3.28518156d-01
    COFD(884) = 1.43341626d-02
    COFD(885) = -1.59525102d+01
    COFD(886) = 3.66023858d+00
    COFD(887) = -2.63401043d-01
    COFD(888) = 1.15432000d-02
    COFD(889) = -1.92867554d+01
    COFD(890) = 4.83375900d+00
    COFD(891) = -4.09146560d-01
    COFD(892) = 1.76006599d-02
    COFD(893) = -2.14022336d+01
    COFD(894) = 5.55346617d+00
    COFD(895) = -4.87783156d-01
    COFD(896) = 2.04210886d-02
    COFD(897) = -2.14082453d+01
    COFD(898) = 5.55346617d+00
    COFD(899) = -4.87783156d-01
    COFD(900) = 2.04210886d-02
    COFD(901) = -2.11381508d+01
    COFD(902) = 5.45574440d+00
    COFD(903) = -4.77436155d-01
    COFD(904) = 2.00644596d-02
    COFD(905) = -1.98075055d+01
    COFD(906) = 5.02169524d+00
    COFD(907) = -4.31582804d-01
    COFD(908) = 1.84953568d-02
    COFD(909) = -1.94763688d+01
    COFD(910) = 4.87333294d+00
    COFD(911) = -4.13769241d-01
    COFD(912) = 1.77802244d-02
    COFD(913) = -1.94823660d+01
    COFD(914) = 4.87333294d+00
    COFD(915) = -4.13769241d-01
    COFD(916) = 1.77802244d-02
    COFD(917) = -1.59327297d+01
    COFD(918) = 3.65620899d+00
    COFD(919) = -2.62933804d-01
    COFD(920) = 1.15253223d-02
    COFD(921) = -1.71754154d+01
    COFD(922) = 4.13131681d+00
    COFD(923) = -3.23897559d-01
    COFD(924) = 1.41438222d-02
    COFD(925) = -1.17159737d+01
    COFD(926) = 2.48123210d+00
    COFD(927) = -1.11322604d-01
    COFD(928) = 4.99282389d-03
    COFD(929) = -1.43151174d+01
    COFD(930) = 3.68038508d+00
    COFD(931) = -2.65779346d-01
    COFD(932) = 1.16360771d-02
    COFD(933) = -1.40999008d+01
    COFD(934) = 3.08120012d+00
    COFD(935) = -1.89629903d-01
    COFD(936) = 8.40361952d-03
    COFD(937) = -1.52721107d+01
    COFD(938) = 3.36790500d+00
    COFD(939) = -2.26321740d-01
    COFD(940) = 9.97135055d-03
    COFD(941) = -1.41191261d+01
    COFD(942) = 3.08120012d+00
    COFD(943) = -1.89629903d-01
    COFD(944) = 8.40361952d-03
    COFD(945) = -2.11388331d+01
    COFD(946) = 5.55529675d+00
    COFD(947) = -4.87942518d-01
    COFD(948) = 2.04249054d-02
    COFD(949) = -1.52792891d+01
    COFD(950) = 3.36790500d+00
    COFD(951) = -2.26321740d-01
    COFD(952) = 9.97135055d-03
    COFD(953) = -1.59634533d+01
    COFD(954) = 3.67388294d+00
    COFD(955) = -2.64990709d-01
    COFD(956) = 1.16042706d-02
    COFD(957) = -1.59634533d+01
    COFD(958) = 3.67388294d+00
    COFD(959) = -2.64990709d-01
    COFD(960) = 1.16042706d-02
    COFD(961) = -1.59863030d+01
    COFD(962) = 3.67388294d+00
    COFD(963) = -2.64990709d-01
    COFD(964) = 1.16042706d-02
    COFD(965) = -1.59525102d+01
    COFD(966) = 3.66023858d+00
    COFD(967) = -2.63401043d-01
    COFD(968) = 1.15432000d-02
    COFD(969) = -1.50233475d+01
    COFD(970) = 3.26660767d+00
    COFD(971) = -2.13287177d-01
    COFD(972) = 9.41137857d-03
    COFD(973) = -1.81735763d+01
    COFD(974) = 4.38391495d+00
    COFD(975) = -3.54941287d-01
    COFD(976) = 1.54195107d-02
    COFD(977) = -2.05045578d+01
    COFD(978) = 5.23843909d+00
    COFD(979) = -4.55815614d-01
    COFD(980) = 1.93898040d-02
    COFD(981) = -2.05128705d+01
    COFD(982) = 5.23843909d+00
    COFD(983) = -4.55815614d-01
    COFD(984) = 1.93898040d-02
    COFD(985) = -2.02642227d+01
    COFD(986) = 5.14499740d+00
    COFD(987) = -4.45694430d-01
    COFD(988) = 1.90318646d-02
    COFD(989) = -1.86157761d+01
    COFD(990) = 4.55689508d+00
    COFD(991) = -3.75937921d-01
    COFD(992) = 1.62703488d-02
    COFD(993) = -1.83455435d+01
    COFD(994) = 4.42828044d+00
    COFD(995) = -3.60417833d-01
    COFD(996) = 1.56455103d-02
    COFD(997) = -1.83538377d+01
    COFD(998) = 4.42828044d+00
    COFD(999) = -3.60417833d-01
    COFD(1000) = 1.56455103d-02
    COFD(1001) = -1.50031687d+01
    COFD(1002) = 3.26223357d+00
    COFD(1003) = -2.12746642d-01
    COFD(1004) = 9.38912883d-03
    COFD(1005) = -1.60074211d+01
    COFD(1006) = 3.63723937d+00
    COFD(1007) = -2.60754222d-01
    COFD(1008) = 1.14428814d-02
    COFD(1009) = -1.37794315d+01
    COFD(1010) = 3.23973858d+00
    COFD(1011) = -2.09989036d-01
    COFD(1012) = 9.27667906d-03
    COFD(1013) = -1.76147026d+01
    COFD(1014) = 4.86049500d+00
    COFD(1015) = -4.12200578d-01
    COFD(1016) = 1.77160971d-02
    COFD(1017) = -1.70534856d+01
    COFD(1018) = 4.14240922d+00
    COFD(1019) = -3.25239774d-01
    COFD(1020) = 1.41980687d-02
    COFD(1021) = -1.84688406d+01
    COFD(1022) = 4.49330851d+00
    COFD(1023) = -3.68208715d-01
    COFD(1024) = 1.59565402d-02
    COFD(1025) = -1.70757047d+01
    COFD(1026) = 4.14240922d+00
    COFD(1027) = -3.25239774d-01
    COFD(1028) = 1.41980687d-02
    COFD(1029) = -2.07653719d+01
    COFD(1030) = 5.01092022d+00
    COFD(1031) = -3.77985635d-01
    COFD(1032) = 1.40968645d-02
    COFD(1033) = -1.84777607d+01
    COFD(1034) = 4.49330851d+00
    COFD(1035) = -3.68208715d-01
    COFD(1036) = 1.59565402d-02
    COFD(1037) = -1.93015555d+01
    COFD(1038) = 4.85015581d+00
    COFD(1039) = -4.10945109d-01
    COFD(1040) = 1.76651398d-02
    COFD(1041) = -1.93015555d+01
    COFD(1042) = 4.85015581d+00
    COFD(1043) = -4.10945109d-01
    COFD(1044) = 1.76651398d-02
    COFD(1045) = -1.93276434d+01
    COFD(1046) = 4.85015581d+00
    COFD(1047) = -4.10945109d-01
    COFD(1048) = 1.76651398d-02
    COFD(1049) = -1.92867554d+01
    COFD(1050) = 4.83375900d+00
    COFD(1051) = -4.09146560d-01
    COFD(1052) = 1.76006599d-02
    COFD(1053) = -1.81735763d+01
    COFD(1054) = 4.38391495d+00
    COFD(1055) = -3.54941287d-01
    COFD(1056) = 1.54195107d-02
    COFD(1057) = -2.13425698d+01
    COFD(1058) = 5.40460130d+00
    COFD(1059) = -4.72718910d-01
    COFD(1060) = 1.99362717d-02
    COFD(1061) = -2.19215555d+01
    COFD(1062) = 5.45216133d+00
    COFD(1063) = -4.52916925d-01
    COFD(1064) = 1.80456400d-02
    COFD(1065) = -2.19317743d+01
    COFD(1066) = 5.45216133d+00
    COFD(1067) = -4.52916925d-01
    COFD(1068) = 1.80456400d-02
    COFD(1069) = -2.20421041d+01
    COFD(1070) = 5.52708332d+00
    COFD(1071) = -4.68000808d-01
    COFD(1072) = 1.89131908d-02
    COFD(1073) = -2.16802612d+01
    COFD(1074) = 5.52918296d+00
    COFD(1075) = -4.85360709d-01
    COFD(1076) = 2.03448006d-02
    COFD(1077) = -2.14224484d+01
    COFD(1078) = 5.41729961d+00
    COFD(1079) = -4.73400281d-01
    COFD(1080) = 1.99269567d-02
    COFD(1081) = -2.14326461d+01
    COFD(1082) = 5.41729961d+00
    COFD(1083) = -4.73400281d-01
    COFD(1084) = 1.99269567d-02
    COFD(1085) = -1.81432461d+01
    COFD(1086) = 4.37565431d+00
    COFD(1087) = -3.53906025d-01
    COFD(1088) = 1.53760786d-02
    COFD(1089) = -1.93483692d+01
    COFD(1090) = 4.79506290d+00
    COFD(1091) = -4.04621659d-01
    COFD(1092) = 1.74244230d-02
    COFD(1093) = -1.60517370d+01
    COFD(1094) = 4.11188603d+00
    COFD(1095) = -3.21540884d-01
    COFD(1096) = 1.40482564d-02
    COFD(1097) = -1.97544450d+01
    COFD(1098) = 5.56931926d+00
    COFD(1099) = -4.89105511d-01
    COFD(1100) = 2.04493129d-02
    COFD(1101) = -1.94313116d+01
    COFD(1102) = 5.02567894d+00
    COFD(1103) = -4.32045169d-01
    COFD(1104) = 1.85132214d-02
    COFD(1105) = -2.08204449d+01
    COFD(1106) = 5.35267674d+00
    COFD(1107) = -4.69010505d-01
    COFD(1108) = 1.98979152d-02
    COFD(1109) = -1.94507876d+01
    COFD(1110) = 5.02567894d+00
    COFD(1111) = -4.32045169d-01
    COFD(1112) = 1.85132214d-02
    COFD(1113) = -1.77498543d+01
    COFD(1114) = 3.57475686d+00
    COFD(1115) = -1.56396297d-01
    COFD(1116) = 3.12157721d-03
    COFD(1117) = -2.08277598d+01
    COFD(1118) = 5.35267674d+00
    COFD(1119) = -4.69010505d-01
    COFD(1120) = 1.98979152d-02
    COFD(1121) = -2.14160703d+01
    COFD(1122) = 5.56531152d+00
    COFD(1123) = -4.88789821d-01
    COFD(1124) = 2.04437116d-02
    COFD(1125) = -2.14160703d+01
    COFD(1126) = 5.56531152d+00
    COFD(1127) = -4.88789821d-01
    COFD(1128) = 2.04437116d-02
    COFD(1129) = -2.14391943d+01
    COFD(1130) = 5.56531152d+00
    COFD(1131) = -4.88789821d-01
    COFD(1132) = 2.04437116d-02
    COFD(1133) = -2.14022336d+01
    COFD(1134) = 5.55346617d+00
    COFD(1135) = -4.87783156d-01
    COFD(1136) = 2.04210886d-02
    COFD(1137) = -2.05045578d+01
    COFD(1138) = 5.23843909d+00
    COFD(1139) = -4.55815614d-01
    COFD(1140) = 1.93898040d-02
    COFD(1141) = -2.19215555d+01
    COFD(1142) = 5.45216133d+00
    COFD(1143) = -4.52916925d-01
    COFD(1144) = 1.80456400d-02
    COFD(1145) = -1.90328712d+01
    COFD(1146) = 3.99221757d+00
    COFD(1147) = -2.19854880d-01
    COFD(1148) = 6.22736279d-03
    COFD(1149) = -1.90413348d+01
    COFD(1150) = 3.99221757d+00
    COFD(1151) = -2.19854880d-01
    COFD(1152) = 6.22736279d-03
    COFD(1153) = -2.01801667d+01
    COFD(1154) = 4.53183330d+00
    COFD(1155) = -3.02186760d-01
    COFD(1156) = 1.02756490d-02
    COFD(1157) = -2.16296373d+01
    COFD(1158) = 5.29019717d+00
    COFD(1159) = -4.24502606d-01
    COFD(1160) = 1.65197343d-02
    COFD(1161) = -2.19229190d+01
    COFD(1162) = 5.41841631d+00
    COFD(1163) = -4.46818971d-01
    COFD(1164) = 1.77127652d-02
    COFD(1165) = -2.19313638d+01
    COFD(1166) = 5.41841631d+00
    COFD(1167) = -4.46818971d-01
    COFD(1168) = 1.77127652d-02
    COFD(1169) = -2.04750581d+01
    COFD(1170) = 5.23112374d+00
    COFD(1171) = -4.54967682d-01
    COFD(1172) = 1.93570423d-02
    COFD(1173) = -2.14255087d+01
    COFD(1174) = 5.52240865d+00
    COFD(1175) = -4.84699537d-01
    COFD(1176) = 2.03247833d-02
    COFD(1177) = -1.60528285d+01
    COFD(1178) = 4.11188603d+00
    COFD(1179) = -3.21540884d-01
    COFD(1180) = 1.40482564d-02
    COFD(1181) = -1.97550088d+01
    COFD(1182) = 5.56931926d+00
    COFD(1183) = -4.89105511d-01
    COFD(1184) = 2.04493129d-02
    COFD(1185) = -1.94373127d+01
    COFD(1186) = 5.02567894d+00
    COFD(1187) = -4.32045169d-01
    COFD(1188) = 1.85132214d-02
    COFD(1189) = -2.08293255d+01
    COFD(1190) = 5.35267674d+00
    COFD(1191) = -4.69010505d-01
    COFD(1192) = 1.98979152d-02
    COFD(1193) = -1.94570287d+01
    COFD(1194) = 5.02567894d+00
    COFD(1195) = -4.32045169d-01
    COFD(1196) = 1.85132214d-02
    COFD(1197) = -1.77563250d+01
    COFD(1198) = 3.57475686d+00
    COFD(1199) = -1.56396297d-01
    COFD(1200) = 3.12157721d-03
    COFD(1201) = -2.08367725d+01
    COFD(1202) = 5.35267674d+00
    COFD(1203) = -4.69010505d-01
    COFD(1204) = 1.98979152d-02
    COFD(1205) = -2.14215700d+01
    COFD(1206) = 5.56531152d+00
    COFD(1207) = -4.88789821d-01
    COFD(1208) = 2.04437116d-02
    COFD(1209) = -2.14215700d+01
    COFD(1210) = 5.56531152d+00
    COFD(1211) = -4.88789821d-01
    COFD(1212) = 2.04437116d-02
    COFD(1213) = -2.14449559d+01
    COFD(1214) = 5.56531152d+00
    COFD(1215) = -4.88789821d-01
    COFD(1216) = 2.04437116d-02
    COFD(1217) = -2.14082453d+01
    COFD(1218) = 5.55346617d+00
    COFD(1219) = -4.87783156d-01
    COFD(1220) = 2.04210886d-02
    COFD(1221) = -2.05128705d+01
    COFD(1222) = 5.23843909d+00
    COFD(1223) = -4.55815614d-01
    COFD(1224) = 1.93898040d-02
    COFD(1225) = -2.19317743d+01
    COFD(1226) = 5.45216133d+00
    COFD(1227) = -4.52916925d-01
    COFD(1228) = 1.80456400d-02
    COFD(1229) = -1.90413348d+01
    COFD(1230) = 3.99221757d+00
    COFD(1231) = -2.19854880d-01
    COFD(1232) = 6.22736279d-03
    COFD(1233) = -1.90499441d+01
    COFD(1234) = 3.99221757d+00
    COFD(1235) = -2.19854880d-01
    COFD(1236) = 6.22736279d-03
    COFD(1237) = -2.01889168d+01
    COFD(1238) = 4.53183330d+00
    COFD(1239) = -3.02186760d-01
    COFD(1240) = 1.02756490d-02
    COFD(1241) = -2.16379567d+01
    COFD(1242) = 5.29019717d+00
    COFD(1243) = -4.24502606d-01
    COFD(1244) = 1.65197343d-02
    COFD(1245) = -2.19313890d+01
    COFD(1246) = 5.41841631d+00
    COFD(1247) = -4.46818971d-01
    COFD(1248) = 1.77127652d-02
    COFD(1249) = -2.19399793d+01
    COFD(1250) = 5.41841631d+00
    COFD(1251) = -4.46818971d-01
    COFD(1252) = 1.77127652d-02
    COFD(1253) = -2.04833713d+01
    COFD(1254) = 5.23112374d+00
    COFD(1255) = -4.54967682d-01
    COFD(1256) = 1.93570423d-02
    COFD(1257) = -2.14353267d+01
    COFD(1258) = 5.52240865d+00
    COFD(1259) = -4.84699537d-01
    COFD(1260) = 2.03247833d-02
    COFD(1261) = -1.58456300d+01
    COFD(1262) = 4.02074783d+00
    COFD(1263) = -3.10018522d-01
    COFD(1264) = 1.35599552d-02
    COFD(1265) = -1.92718582d+01
    COFD(1266) = 5.41172124d+00
    COFD(1267) = -4.73213887d-01
    COFD(1268) = 1.99405473d-02
    COFD(1269) = -1.88179418d+01
    COFD(1270) = 4.79683898d+00
    COFD(1271) = -4.04829719d-01
    COFD(1272) = 1.74325475d-02
    COFD(1273) = -2.04928958d+01
    COFD(1274) = 5.22397933d+00
    COFD(1275) = -4.54138171d-01
    COFD(1276) = 1.93249285d-02
    COFD(1277) = -1.88378874d+01
    COFD(1278) = 4.79683898d+00
    COFD(1279) = -4.04829719d-01
    COFD(1280) = 1.74325475d-02
    COFD(1281) = -1.65295288d+01
    COFD(1282) = 2.97569206d+00
    COFD(1283) = -6.75652842d-02
    COFD(1284) = -1.08648422d-03
    COFD(1285) = -2.02637994d+01
    COFD(1286) = 5.14984081d+00
    COFD(1287) = -4.46093018d-01
    COFD(1288) = 1.90396647d-02
    COFD(1289) = -2.09376196d+01
    COFD(1290) = 5.40870099d+00
    COFD(1291) = -4.73017610d-01
    COFD(1292) = 1.99399066d-02
    COFD(1293) = -2.09376196d+01
    COFD(1294) = 5.40870099d+00
    COFD(1295) = -4.73017610d-01
    COFD(1296) = 1.99399066d-02
    COFD(1297) = -2.09612557d+01
    COFD(1298) = 5.40870099d+00
    COFD(1299) = -4.73017610d-01
    COFD(1300) = 1.99399066d-02
    COFD(1301) = -2.11381508d+01
    COFD(1302) = 5.45574440d+00
    COFD(1303) = -4.77436155d-01
    COFD(1304) = 2.00644596d-02
    COFD(1305) = -2.02642227d+01
    COFD(1306) = 5.14499740d+00
    COFD(1307) = -4.45694430d-01
    COFD(1308) = 1.90318646d-02
    COFD(1309) = -2.20421041d+01
    COFD(1310) = 5.52708332d+00
    COFD(1311) = -4.68000808d-01
    COFD(1312) = 1.89131908d-02
    COFD(1313) = -2.01801667d+01
    COFD(1314) = 4.53183330d+00
    COFD(1315) = -3.02186760d-01
    COFD(1316) = 1.02756490d-02
    COFD(1317) = -2.01889168d+01
    COFD(1318) = 4.53183330d+00
    COFD(1319) = -3.02186760d-01
    COFD(1320) = 1.02756490d-02
    COFD(1321) = -1.95877017d+01
    COFD(1322) = 4.27643051d+00
    COFD(1323) = -2.68040901d-01
    COFD(1324) = 8.77650113d-03
    COFD(1325) = -2.19670848d+01
    COFD(1326) = 5.48847873d+00
    COFD(1327) = -4.59558930d-01
    COFD(1328) = 1.84107961d-02
    COFD(1329) = -2.21070030d+01
    COFD(1330) = 5.55072945d+00
    COFD(1331) = -4.72525345d-01
    COFD(1332) = 1.91674202d-02
    COFD(1333) = -2.21157340d+01
    COFD(1334) = 5.55072945d+00
    COFD(1335) = -4.72525345d-01
    COFD(1336) = 1.91674202d-02
    COFD(1337) = -2.02268902d+01
    COFD(1338) = 5.13632093d+00
    COFD(1339) = -4.44839124d-01
    COFD(1340) = 1.90058354d-02
    COFD(1341) = -2.10026861d+01
    COFD(1342) = 5.38326647d+00
    COFD(1343) = -4.71201048d-01
    COFD(1344) = 1.99207516d-02
    COFD(1345) = -1.42229194d+01
    COFD(1346) = 3.38669384d+00
    COFD(1347) = -2.28784122d-01
    COFD(1348) = 1.00790953d-02
    COFD(1349) = -1.82251914d+01
    COFD(1350) = 5.05237312d+00
    COFD(1351) = -4.35182396d-01
    COFD(1352) = 1.86363074d-02
    COFD(1353) = -1.74792112d+01
    COFD(1354) = 4.29676909d+00
    COFD(1355) = -3.44085306d-01
    COFD(1356) = 1.49671135d-02
    COFD(1357) = -1.89544778d+01
    COFD(1358) = 4.68595732d+00
    COFD(1359) = -3.91842840d-01
    COFD(1360) = 1.69262542d-02
    COFD(1361) = -1.74984476d+01
    COFD(1362) = 4.29676909d+00
    COFD(1363) = -3.44085306d-01
    COFD(1364) = 1.49671135d-02
    COFD(1365) = -2.08812333d+01
    COFD(1366) = 5.08859217d+00
    COFD(1367) = -3.90525428d-01
    COFD(1368) = 1.47376395d-02
    COFD(1369) = -1.89616623d+01
    COFD(1370) = 4.68595732d+00
    COFD(1371) = -3.91842840d-01
    COFD(1372) = 1.69262542d-02
    COFD(1373) = -1.98418115d+01
    COFD(1374) = 5.04367502d+00
    COFD(1375) = -4.34153325d-01
    COFD(1376) = 1.85956055d-02
    COFD(1377) = -1.98418115d+01
    COFD(1378) = 5.04367502d+00
    COFD(1379) = -4.34153325d-01
    COFD(1380) = 1.85956055d-02
    COFD(1381) = -1.98646734d+01
    COFD(1382) = 5.04367502d+00
    COFD(1383) = -4.34153325d-01
    COFD(1384) = 1.85956055d-02
    COFD(1385) = -1.98075055d+01
    COFD(1386) = 5.02169524d+00
    COFD(1387) = -4.31582804d-01
    COFD(1388) = 1.84953568d-02
    COFD(1389) = -1.86157761d+01
    COFD(1390) = 4.55689508d+00
    COFD(1391) = -3.75937921d-01
    COFD(1392) = 1.62703488d-02
    COFD(1393) = -2.16802612d+01
    COFD(1394) = 5.52918296d+00
    COFD(1395) = -4.85360709d-01
    COFD(1396) = 2.03448006d-02
    COFD(1397) = -2.16296373d+01
    COFD(1398) = 5.29019717d+00
    COFD(1399) = -4.24502606d-01
    COFD(1400) = 1.65197343d-02
    COFD(1401) = -2.16379567d+01
    COFD(1402) = 5.29019717d+00
    COFD(1403) = -4.24502606d-01
    COFD(1404) = 1.65197343d-02
    COFD(1405) = -2.19670848d+01
    COFD(1406) = 5.48847873d+00
    COFD(1407) = -4.59558930d-01
    COFD(1408) = 1.84107961d-02
    COFD(1409) = -2.19327397d+01
    COFD(1410) = 5.60638188d+00
    COFD(1411) = -4.91272522d-01
    COFD(1412) = 2.04396264d-02
    COFD(1413) = -2.18190539d+01
    COFD(1414) = 5.55753905d+00
    COFD(1415) = -4.88136714d-01
    COFD(1416) = 2.04294957d-02
    COFD(1417) = -2.18273547d+01
    COFD(1418) = 5.55753905d+00
    COFD(1419) = -4.88136714d-01
    COFD(1420) = 2.04294957d-02
    COFD(1421) = -1.85864144d+01
    COFD(1422) = 4.54915847d+00
    COFD(1423) = -3.75000738d-01
    COFD(1424) = 1.62324821d-02
    COFD(1425) = -1.98040322d+01
    COFD(1426) = 4.97569695d+00
    COFD(1427) = -4.26123307d-01
    COFD(1428) = 1.82788664d-02
    COFD(1429) = -1.39913897d+01
    COFD(1430) = 3.26384506d+00
    COFD(1431) = -2.12947087d-01
    COFD(1432) = 9.39743888d-03
    COFD(1433) = -1.79339327d+01
    COFD(1434) = 4.91373893d+00
    COFD(1435) = -4.18747629d-01
    COFD(1436) = 1.79856610d-02
    COFD(1437) = -1.72496634d+01
    COFD(1438) = 4.17889917d+00
    COFD(1439) = -3.29752510d-01
    COFD(1440) = 1.43850275d-02
    COFD(1441) = -1.86335932d+01
    COFD(1442) = 4.53572533d+00
    COFD(1443) = -3.73386925d-01
    COFD(1444) = 1.61678881d-02
    COFD(1445) = -1.72691500d+01
    COFD(1446) = 4.17889917d+00
    COFD(1447) = -3.29752510d-01
    COFD(1448) = 1.43850275d-02
    COFD(1449) = -2.12597312d+01
    COFD(1450) = 5.24930667d+00
    COFD(1451) = -4.17435088d-01
    COFD(1452) = 1.61434424d-02
    COFD(1453) = -1.86409139d+01
    COFD(1454) = 4.53572533d+00
    COFD(1455) = -3.73386925d-01
    COFD(1456) = 1.61678881d-02
    COFD(1457) = -1.95263312d+01
    COFD(1458) = 4.90255048d+00
    COFD(1459) = -4.17368501d-01
    COFD(1460) = 1.79287358d-02
    COFD(1461) = -1.95263312d+01
    COFD(1462) = 4.90255048d+00
    COFD(1463) = -4.17368501d-01
    COFD(1464) = 1.79287358d-02
    COFD(1465) = -1.95494668d+01
    COFD(1466) = 4.90255048d+00
    COFD(1467) = -4.17368501d-01
    COFD(1468) = 1.79287358d-02
    COFD(1469) = -1.94763688d+01
    COFD(1470) = 4.87333294d+00
    COFD(1471) = -4.13769241d-01
    COFD(1472) = 1.77802244d-02
    COFD(1473) = -1.83455435d+01
    COFD(1474) = 4.42828044d+00
    COFD(1475) = -3.60417833d-01
    COFD(1476) = 1.56455103d-02
    COFD(1477) = -2.14224484d+01
    COFD(1478) = 5.41729961d+00
    COFD(1479) = -4.73400281d-01
    COFD(1480) = 1.99269567d-02
    COFD(1481) = -2.19229190d+01
    COFD(1482) = 5.41841631d+00
    COFD(1483) = -4.46818971d-01
    COFD(1484) = 1.77127652d-02
    COFD(1485) = -2.19313890d+01
    COFD(1486) = 5.41841631d+00
    COFD(1487) = -4.46818971d-01
    COFD(1488) = 1.77127652d-02
    COFD(1489) = -2.21070030d+01
    COFD(1490) = 5.55072945d+00
    COFD(1491) = -4.72525345d-01
    COFD(1492) = 1.91674202d-02
    COFD(1493) = -2.18190539d+01
    COFD(1494) = 5.55753905d+00
    COFD(1495) = -4.88136714d-01
    COFD(1496) = 2.04294957d-02
    COFD(1497) = -2.15575659d+01
    COFD(1498) = 5.44803850d+00
    COFD(1499) = -4.76610560d-01
    COFD(1500) = 2.00355294d-02
    COFD(1501) = -2.15660171d+01
    COFD(1502) = 5.44803850d+00
    COFD(1503) = -4.76610560d-01
    COFD(1504) = 2.00355294d-02
    COFD(1505) = -1.83166353d+01
    COFD(1506) = 4.42045763d+00
    COFD(1507) = -3.59451578d-01
    COFD(1508) = 1.56056164d-02
    COFD(1509) = -1.94928815d+01
    COFD(1510) = 4.83189721d+00
    COFD(1511) = -4.08932249d-01
    COFD(1512) = 1.75924650d-02
    COFD(1513) = -1.39924781d+01
    COFD(1514) = 3.26384506d+00
    COFD(1515) = -2.12947087d-01
    COFD(1516) = 9.39743888d-03
    COFD(1517) = -1.79344949d+01
    COFD(1518) = 4.91373893d+00
    COFD(1519) = -4.18747629d-01
    COFD(1520) = 1.79856610d-02
    COFD(1521) = -1.72556499d+01
    COFD(1522) = 4.17889917d+00
    COFD(1523) = -3.29752510d-01
    COFD(1524) = 1.43850275d-02
    COFD(1525) = -1.86424545d+01
    COFD(1526) = 4.53572533d+00
    COFD(1527) = -3.73386925d-01
    COFD(1528) = 1.61678881d-02
    COFD(1529) = -1.72753760d+01
    COFD(1530) = 4.17889917d+00
    COFD(1531) = -3.29752510d-01
    COFD(1532) = 1.43850275d-02
    COFD(1533) = -2.12661865d+01
    COFD(1534) = 5.24930667d+00
    COFD(1535) = -4.17435088d-01
    COFD(1536) = 1.61434424d-02
    COFD(1537) = -1.86499071d+01
    COFD(1538) = 4.53572533d+00
    COFD(1539) = -3.73386925d-01
    COFD(1540) = 1.61678881d-02
    COFD(1541) = -1.95318173d+01
    COFD(1542) = 4.90255048d+00
    COFD(1543) = -4.17368501d-01
    COFD(1544) = 1.79287358d-02
    COFD(1545) = -1.95318173d+01
    COFD(1546) = 4.90255048d+00
    COFD(1547) = -4.17368501d-01
    COFD(1548) = 1.79287358d-02
    COFD(1549) = -1.95552142d+01
    COFD(1550) = 4.90255048d+00
    COFD(1551) = -4.17368501d-01
    COFD(1552) = 1.79287358d-02
    COFD(1553) = -1.94823660d+01
    COFD(1554) = 4.87333294d+00
    COFD(1555) = -4.13769241d-01
    COFD(1556) = 1.77802244d-02
    COFD(1557) = -1.83538377d+01
    COFD(1558) = 4.42828044d+00
    COFD(1559) = -3.60417833d-01
    COFD(1560) = 1.56455103d-02
    COFD(1561) = -2.14326461d+01
    COFD(1562) = 5.41729961d+00
    COFD(1563) = -4.73400281d-01
    COFD(1564) = 1.99269567d-02
    COFD(1565) = -2.19313638d+01
    COFD(1566) = 5.41841631d+00
    COFD(1567) = -4.46818971d-01
    COFD(1568) = 1.77127652d-02
    COFD(1569) = -2.19399793d+01
    COFD(1570) = 5.41841631d+00
    COFD(1571) = -4.46818971d-01
    COFD(1572) = 1.77127652d-02
    COFD(1573) = -2.21157340d+01
    COFD(1574) = 5.55072945d+00
    COFD(1575) = -4.72525345d-01
    COFD(1576) = 1.91674202d-02
    COFD(1577) = -2.18273547d+01
    COFD(1578) = 5.55753905d+00
    COFD(1579) = -4.88136714d-01
    COFD(1580) = 2.04294957d-02
    COFD(1581) = -2.15660171d+01
    COFD(1582) = 5.44803850d+00
    COFD(1583) = -4.76610560d-01
    COFD(1584) = 2.00355294d-02
    COFD(1585) = -2.15746136d+01
    COFD(1586) = 5.44803850d+00
    COFD(1587) = -4.76610560d-01
    COFD(1588) = 2.00355294d-02
    COFD(1589) = -1.83249299d+01
    COFD(1590) = 4.42045763d+00
    COFD(1591) = -3.59451578d-01
    COFD(1592) = 1.56056164d-02
    COFD(1593) = -1.95026789d+01
    COFD(1594) = 4.83189721d+00
    COFD(1595) = -4.08932249d-01
    COFD(1596) = 1.75924650d-02
    COFD(1597) = -1.16906297d+01
    COFD(1598) = 2.47469981d+00
    COFD(1599) = -1.10436257d-01
    COFD(1600) = 4.95273813d-03
    COFD(1601) = -1.42894441d+01
    COFD(1602) = 3.67490723d+00
    COFD(1603) = -2.65114792d-01
    COFD(1604) = 1.16092671d-02
    COFD(1605) = -1.40756935d+01
    COFD(1606) = 3.07549274d+00
    COFD(1607) = -1.88889344d-01
    COFD(1608) = 8.37152866d-03
    COFD(1609) = -1.52414485d+01
    COFD(1610) = 3.35922578d+00
    COFD(1611) = -2.25181399d-01
    COFD(1612) = 9.92132878d-03
    COFD(1613) = -1.40949196d+01
    COFD(1614) = 3.07549274d+00
    COFD(1615) = -1.88889344d-01
    COFD(1616) = 8.37152866d-03
    COFD(1617) = -2.10643259d+01
    COFD(1618) = 5.53614847d+00
    COFD(1619) = -4.86046736d-01
    COFD(1620) = 2.03659188d-02
    COFD(1621) = -1.52486273d+01
    COFD(1622) = 3.35922578d+00
    COFD(1623) = -2.25181399d-01
    COFD(1624) = 9.92132878d-03
    COFD(1625) = -1.59404882d+01
    COFD(1626) = 3.66853818d+00
    COFD(1627) = -2.64346221d-01
    COFD(1628) = 1.15784613d-02
    COFD(1629) = -1.59404882d+01
    COFD(1630) = 3.66853818d+00
    COFD(1631) = -2.64346221d-01
    COFD(1632) = 1.15784613d-02
    COFD(1633) = -1.59633387d+01
    COFD(1634) = 3.66853818d+00
    COFD(1635) = -2.64346221d-01
    COFD(1636) = 1.15784613d-02
    COFD(1637) = -1.59327297d+01
    COFD(1638) = 3.65620899d+00
    COFD(1639) = -2.62933804d-01
    COFD(1640) = 1.15253223d-02
    COFD(1641) = -1.50031687d+01
    COFD(1642) = 3.26223357d+00
    COFD(1643) = -2.12746642d-01
    COFD(1644) = 9.38912883d-03
    COFD(1645) = -1.81432461d+01
    COFD(1646) = 4.37565431d+00
    COFD(1647) = -3.53906025d-01
    COFD(1648) = 1.53760786d-02
    COFD(1649) = -2.04750581d+01
    COFD(1650) = 5.23112374d+00
    COFD(1651) = -4.54967682d-01
    COFD(1652) = 1.93570423d-02
    COFD(1653) = -2.04833713d+01
    COFD(1654) = 5.23112374d+00
    COFD(1655) = -4.54967682d-01
    COFD(1656) = 1.93570423d-02
    COFD(1657) = -2.02268902d+01
    COFD(1658) = 5.13632093d+00
    COFD(1659) = -4.44839124d-01
    COFD(1660) = 1.90058354d-02
    COFD(1661) = -1.85864144d+01
    COFD(1662) = 4.54915847d+00
    COFD(1663) = -3.75000738d-01
    COFD(1664) = 1.62324821d-02
    COFD(1665) = -1.83166353d+01
    COFD(1666) = 4.42045763d+00
    COFD(1667) = -3.59451578d-01
    COFD(1668) = 1.56056164d-02
    COFD(1669) = -1.83249299d+01
    COFD(1670) = 4.42045763d+00
    COFD(1671) = -3.59451578d-01
    COFD(1672) = 1.56056164d-02
    COFD(1673) = -1.49828430d+01
    COFD(1674) = 3.25781069d+00
    COFD(1675) = -2.12199367d-01
    COFD(1676) = 9.36657283d-03
    COFD(1677) = -1.59877782d+01
    COFD(1678) = 3.63340763d+00
    COFD(1679) = -2.60307961d-01
    COFD(1680) = 1.14256954d-02
    COFD(1681) = -1.23130152d+01
    COFD(1682) = 2.74418790d+00
    COFD(1683) = -1.46230156d-01
    COFD(1684) = 6.53948886d-03
    COFD(1685) = -1.54738604d+01
    COFD(1686) = 4.15765300d+00
    COFD(1687) = -3.27126237d-01
    COFD(1688) = 1.42762611d-02
    COFD(1689) = -1.49610527d+01
    COFD(1690) = 3.41988961d+00
    COFD(1691) = -2.33128386d-01
    COFD(1692) = 1.02689994d-02
    COFD(1693) = -1.62380075d+01
    COFD(1694) = 3.72612300d+00
    COFD(1695) = -2.71663673d-01
    COFD(1696) = 1.18889643d-02
    COFD(1697) = -1.49826725d+01
    COFD(1698) = 3.41988961d+00
    COFD(1699) = -2.33128386d-01
    COFD(1700) = 1.02689994d-02
    COFD(1701) = -2.12755888d+01
    COFD(1702) = 5.60381989d+00
    COFD(1703) = -4.91225459d-01
    COFD(1704) = 2.04487844d-02
    COFD(1705) = -1.62465583d+01
    COFD(1706) = 3.72612300d+00
    COFD(1707) = -2.71663673d-01
    COFD(1708) = 1.18889643d-02
    COFD(1709) = -1.71942502d+01
    COFD(1710) = 4.14993355d+00
    COFD(1711) = -3.26168062d-01
    COFD(1712) = 1.42364115d-02
    COFD(1713) = -1.71942502d+01
    COFD(1714) = 4.14993355d+00
    COFD(1715) = -3.26168062d-01
    COFD(1716) = 1.42364115d-02
    COFD(1717) = -1.72196961d+01
    COFD(1718) = 4.14993355d+00
    COFD(1719) = -3.26168062d-01
    COFD(1720) = 1.42364115d-02
    COFD(1721) = -1.71754154d+01
    COFD(1722) = 4.13131681d+00
    COFD(1723) = -3.23897559d-01
    COFD(1724) = 1.41438222d-02
    COFD(1725) = -1.60074211d+01
    COFD(1726) = 3.63723937d+00
    COFD(1727) = -2.60754222d-01
    COFD(1728) = 1.14428814d-02
    COFD(1729) = -1.93483692d+01
    COFD(1730) = 4.79506290d+00
    COFD(1731) = -4.04621659d-01
    COFD(1732) = 1.74244230d-02
    COFD(1733) = -2.14255087d+01
    COFD(1734) = 5.52240865d+00
    COFD(1735) = -4.84699537d-01
    COFD(1736) = 2.03247833d-02
    COFD(1737) = -2.14353267d+01
    COFD(1738) = 5.52240865d+00
    COFD(1739) = -4.84699537d-01
    COFD(1740) = 2.03247833d-02
    COFD(1741) = -2.10026861d+01
    COFD(1742) = 5.38326647d+00
    COFD(1743) = -4.71201048d-01
    COFD(1744) = 1.99207516d-02
    COFD(1745) = -1.98040322d+01
    COFD(1746) = 4.97569695d+00
    COFD(1747) = -4.26123307d-01
    COFD(1748) = 1.82788664d-02
    COFD(1749) = -1.94928815d+01
    COFD(1750) = 4.83189721d+00
    COFD(1751) = -4.08932249d-01
    COFD(1752) = 1.75924650d-02
    COFD(1753) = -1.95026789d+01
    COFD(1754) = 4.83189721d+00
    COFD(1755) = -4.08932249d-01
    COFD(1756) = 1.75924650d-02
    COFD(1757) = -1.59877782d+01
    COFD(1758) = 3.63340763d+00
    COFD(1759) = -2.60307961d-01
    COFD(1760) = 1.14256954d-02
    COFD(1761) = -1.72273911d+01
    COFD(1762) = 4.09361913d+00
    COFD(1763) = -3.19258125d-01
    COFD(1764) = 1.39526981d-02

end subroutine

! List of specs with small weight, dim NLITE
subroutine egtransetKTDIF(KTDIF)

    implicit none

    integer, intent(out) :: KTDIF(2)

    KTDIF(1) = 1
    KTDIF(2) = 2

end subroutine


! Poly fits for thermal diff ratios, dim NO*NLITE*KK
subroutine egtransetCOFTD(COFTD)

    implicit none

    double precision, intent(out) :: COFTD(168)

    COFTD(1) = 0.00000000d+00
    COFTD(2) = 0.00000000d+00
    COFTD(3) = 0.00000000d+00
    COFTD(4) = 0.00000000d+00
    COFTD(5) = -1.44152190d-01
    COFTD(6) = -7.99993584d-05
    COFTD(7) = 4.89707442d-08
    COFTD(8) = -9.14277269d-12
    COFTD(9) = 4.06682492d-01
    COFTD(10) = 3.84705248d-05
    COFTD(11) = -2.54846868d-08
    COFTD(12) = 5.86302354d-12
    COFTD(13) = 4.26579943d-01
    COFTD(14) = 1.20407274d-04
    COFTD(15) = -7.67298757d-08
    COFTD(16) = 1.52090336d-11
    COFTD(17) = 4.12895615d-01
    COFTD(18) = 3.90582612d-05
    COFTD(19) = -2.58740310d-08
    COFTD(20) = 5.95259633d-12
    COFTD(21) = 2.27469146d-02
    COFTD(22) = 6.73078907d-04
    COFTD(23) = -3.40935843d-07
    COFTD(24) = 5.48499211d-11
    COFTD(25) = 4.28230888d-01
    COFTD(26) = 1.20873273d-04
    COFTD(27) = -7.70268349d-08
    COFTD(28) = 1.52678954d-11
    COFTD(29) = 3.24747031d-01
    COFTD(30) = 1.77798548d-04
    COFTD(31) = -1.08934732d-07
    COFTD(32) = 2.03595881d-11
    COFTD(33) = 3.24747031d-01
    COFTD(34) = 1.77798548d-04
    COFTD(35) = -1.08934732d-07
    COFTD(36) = 2.03595881d-11
    COFTD(37) = 3.31191185d-01
    COFTD(38) = 1.81326714d-04
    COFTD(39) = -1.11096391d-07
    COFTD(40) = 2.07635959d-11
    COFTD(41) = 3.39557243d-01
    COFTD(42) = 1.79335036d-04
    COFTD(43) = -1.10135705d-07
    COFTD(44) = 2.06427239d-11
    COFTD(45) = 4.30605547d-01
    COFTD(46) = 9.35961902d-05
    COFTD(47) = -6.03983623d-08
    COFTD(48) = 1.23115170d-11
    COFTD(49) = 2.93191523d-01
    COFTD(50) = 4.01430006d-04
    COFTD(51) = -2.30705763d-07
    COFTD(52) = 4.05176586d-11
    COFTD(53) = 1.22119780d-01
    COFTD(54) = 6.18373616d-04
    COFTD(55) = -3.28422593d-07
    COFTD(56) = 5.44603522d-11
    COFTD(57) = 1.22693382d-01
    COFTD(58) = 6.21278143d-04
    COFTD(59) = -3.29965208d-07
    COFTD(60) = 5.47161548d-11
    COFTD(61) = 1.40314191d-01
    COFTD(62) = 6.01266129d-04
    COFTD(63) = -3.21915137d-07
    COFTD(64) = 5.36679068d-11
    COFTD(65) = 2.49017478d-01
    COFTD(66) = 4.29036573d-04
    COFTD(67) = -2.42668617d-07
    COFTD(68) = 4.20801371d-11
    COFTD(69) = 2.72759599d-01
    COFTD(70) = 3.94402719d-04
    COFTD(71) = -2.25800520d-07
    COFTD(72) = 3.95325634d-11
    COFTD(73) = 2.74036956d-01
    COFTD(74) = 3.96249742d-04
    COFTD(75) = -2.26857964d-07
    COFTD(76) = 3.97176979d-11
    COFTD(77) = 4.31331269d-01
    COFTD(78) = 9.20536800d-05
    COFTD(79) = -5.94509616d-08
    COFTD(80) = 1.21437993d-11
    COFTD(81) = 4.01012808d-01
    COFTD(82) = 1.97252826d-04
    COFTD(83) = -1.21698146d-07
    COFTD(84) = 2.29408847d-11
    COFTD(85) = 1.44152190d-01
    COFTD(86) = 7.99993584d-05
    COFTD(87) = -4.89707442d-08
    COFTD(88) = 9.14277269d-12
    COFTD(89) = 0.00000000d+00
    COFTD(90) = 0.00000000d+00
    COFTD(91) = 0.00000000d+00
    COFTD(92) = 0.00000000d+00
    COFTD(93) = 2.35283119d-01
    COFTD(94) = 4.65670599d-04
    COFTD(95) = -2.60939824d-07
    COFTD(96) = 4.49271822d-11
    COFTD(97) = 1.79840299d-01
    COFTD(98) = 6.01722902d-04
    COFTD(99) = -3.26433894d-07
    COFTD(100) = 5.49112302d-11
    COFTD(101) = 2.37053352d-01
    COFTD(102) = 4.69174231d-04
    COFTD(103) = -2.62903094d-07
    COFTD(104) = 4.52652072d-11
    COFTD(105) = -1.74352698d-01
    COFTD(106) = 8.62246873d-04
    COFTD(107) = -3.79545489d-07
    COFTD(108) = 5.60262093d-11
    COFTD(109) = 1.80186965d-01
    COFTD(110) = 6.02882805d-04
    COFTD(111) = -3.27063140d-07
    COFTD(112) = 5.50170790d-11
    COFTD(113) = 9.90752318d-02
    COFTD(114) = 6.44201384d-04
    COFTD(115) = -3.38485953d-07
    COFTD(116) = 5.57356746d-11
    COFTD(117) = 9.90752318d-02
    COFTD(118) = 6.44201384d-04
    COFTD(119) = -3.38485953d-07
    COFTD(120) = 5.57356746d-11
    COFTD(121) = 1.00039110d-01
    COFTD(122) = 6.50468660d-04
    COFTD(123) = -3.41778999d-07
    COFTD(124) = 5.62779132d-11
    COFTD(125) = 1.05124122d-01
    COFTD(126) = 6.50665957d-04
    COFTD(127) = -3.42564538d-07
    COFTD(128) = 5.64804120d-11
    COFTD(129) = 2.00119897d-01
    COFTD(130) = 5.64793704d-04
    COFTD(131) = -3.09445484d-07
    COFTD(132) = 5.24139335d-11
    COFTD(133) = -2.00309448d-02
    COFTD(134) = 8.50440115d-04
    COFTD(135) = -4.21064468d-07
    COFTD(136) = 6.67959710d-11
    COFTD(137) = -1.60981264d-01
    COFTD(138) = 9.03807572d-04
    COFTD(139) = -4.06927941d-07
    COFTD(140) = 6.09202254d-11
    COFTD(141) = -1.61357564d-01
    COFTD(142) = 9.05920260d-04
    COFTD(143) = -4.07879153d-07
    COFTD(144) = 6.10626290d-11
    COFTD(145) = -1.31244519d-01
    COFTD(146) = 9.03901384d-04
    COFTD(147) = -4.17831507d-07
    COFTD(148) = 6.35725667d-11
    COFTD(149) = -5.08744745d-02
    COFTD(150) = 8.54342586d-04
    COFTD(151) = -4.15926453d-07
    COFTD(152) = 6.53063261d-11
    COFTD(153) = -2.71690558d-02
    COFTD(154) = 8.37233133d-04
    COFTD(155) = -4.12887636d-07
    COFTD(156) = 6.53405197d-11
    COFTD(157) = -2.72323768d-02
    COFTD(158) = 8.39184413d-04
    COFTD(159) = -4.13849924d-07
    COFTD(160) = 6.54928043d-11
    COFTD(161) = 2.01521643d-01
    COFTD(162) = 5.62744089d-04
    COFTD(163) = -3.08519239d-07
    COFTD(164) = 5.22805986d-11
    COFTD(165) = 1.22193921d-01
    COFTD(166) = 6.90321128d-04
    COFTD(167) = -3.64844875d-07
    COFTD(168) = 6.03054876d-11

end subroutine

end module fuego_module



