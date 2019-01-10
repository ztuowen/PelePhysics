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
double precision, parameter :: imw(53) = (/ &
    1.d0 / 2.015940d0,  & ! H2
    1.d0 / 1.007970d0,  & ! H
    1.d0 / 15.999400d0,  & ! O
    1.d0 / 31.998800d0,  & ! O2
    1.d0 / 17.007370d0,  & ! OH
    1.d0 / 18.015340d0,  & ! H2O
    1.d0 / 33.006770d0,  & ! HO2
    1.d0 / 34.014740d0,  & ! H2O2
    1.d0 / 12.011150d0,  & ! C
    1.d0 / 13.019120d0,  & ! CH
    1.d0 / 14.027090d0,  & ! CH2
    1.d0 / 14.027090d0,  & ! CH2(S)
    1.d0 / 15.035060d0,  & ! CH3
    1.d0 / 16.043030d0,  & ! CH4
    1.d0 / 28.010550d0,  & ! CO
    1.d0 / 44.009950d0,  & ! CO2
    1.d0 / 29.018520d0,  & ! HCO
    1.d0 / 30.026490d0,  & ! CH2O
    1.d0 / 31.034460d0,  & ! CH2OH
    1.d0 / 31.034460d0,  & ! CH3O
    1.d0 / 32.042430d0,  & ! CH3OH
    1.d0 / 25.030270d0,  & ! C2H
    1.d0 / 26.038240d0,  & ! C2H2
    1.d0 / 27.046210d0,  & ! C2H3
    1.d0 / 28.054180d0,  & ! C2H4
    1.d0 / 29.062150d0,  & ! C2H5
    1.d0 / 30.070120d0,  & ! C2H6
    1.d0 / 41.029670d0,  & ! HCCO
    1.d0 / 42.037640d0,  & ! CH2CO
    1.d0 / 42.037640d0,  & ! HCCOH
    1.d0 / 14.006700d0,  & ! N
    1.d0 / 15.014670d0,  & ! NH
    1.d0 / 16.022640d0,  & ! NH2
    1.d0 / 17.030610d0,  & ! NH3
    1.d0 / 29.021370d0,  & ! NNH
    1.d0 / 30.006100d0,  & ! NO
    1.d0 / 46.005500d0,  & ! NO2
    1.d0 / 44.012800d0,  & ! N2O
    1.d0 / 31.014070d0,  & ! HNO
    1.d0 / 26.017850d0,  & ! CN
    1.d0 / 27.025820d0,  & ! HCN
    1.d0 / 28.033790d0,  & ! H2CN
    1.d0 / 41.032520d0,  & ! HCNN
    1.d0 / 43.025220d0,  & ! HCNO
    1.d0 / 43.025220d0,  & ! HOCN
    1.d0 / 43.025220d0,  & ! HNCO
    1.d0 / 42.017250d0,  & ! NCO
    1.d0 / 28.013400d0,  & ! N2
    1.d0 / 39.948000d0,  & ! AR
    1.d0 / 43.089240d0,  & ! C3H7
    1.d0 / 44.097210d0,  & ! C3H8
    1.d0 / 43.045610d0,  & ! CH2CHO
    1.d0 / 44.053580d0 /) ! CH3CHO

type :: nonsquare_matrix_double
   double precision, allocatable :: vector(:)
end type nonsquare_matrix_double

type :: nonsquare_matrix_int
   integer, allocatable :: vector(:)
end type nonsquare_matrix_int

double precision, save :: fwd_A(325), fwd_beta(325), fwd_Ea(325)
double precision, save :: low_A(325), low_beta(325), low_Ea(325)
double precision, save :: rev_A(325), rev_beta(325), rev_Ea(325)
double precision, save :: troe_a(325),troe_Ts(325), troe_Tss(325), troe_Tsss(325)
double precision, save :: sri_a(325), sri_b(325), sri_c(325), sri_d(325), sri_e(325)
double precision, save :: activation_units(325), prefactor_units(325), phase_units(325)
integer, save :: is_PD(325), troe_len(325), sri_len(325), nTB(325)
type(nonsquare_matrix_double) :: TB(325)
type(nonsquare_matrix_int) :: TBid(325)

double precision, save :: fwd_A_DEF(325), fwd_beta_DEF(325), fwd_Ea_DEF(325)
double precision, save :: low_A_DEF(325), low_beta_DEF(325), low_Ea_DEF(325)
double precision, save :: rev_A_DEF(325), rev_beta_DEF(325), rev_Ea_DEF(325)
double precision, save :: troe_a_DEF(325),troe_Ts_DEF(325), troe_Tss_DEF(325), troe_Tsss_DEF(325)
double precision, save :: sri_a_DEF(325), sri_b_DEF(325), sri_c_DEF(325), sri_d_DEF(325), sri_e_DEF(325)
double precision, save :: activation_units_DEF(325), prefactor_units_DEF(325), phase_units_DEF(325)
integer, save :: is_PD_DEF(325), troe_len_DEF(325), sri_len_DEF(325), nTB_DEF(325)
type(nonsquare_matrix_double) :: TB_DEF(325)
type(nonsquare_matrix_int) :: TBid_DEF(325)

! productionRate() static variables
double precision, save :: T_save = -1
double precision, save :: k_f_save(325)
double precision, save :: Kc_save(325)

contains

subroutine SetAllDefaults()

    implicit none

    integer :: i, j

    do i=1, 325
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

  do i=1, 325
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

    ! (0):  2 O + M <=> O2 + M
    fwd_A(30)     = 1.2d+17
    fwd_beta(30)  = -1d0
    fwd_Ea(30)    = 0d0
    prefactor_units(30)  = 1.0000000000000002d-12
    activation_units(30) = 0.50321666580471969d0
    phase_units(30)      = 1d-12
    is_PD(30) = 0
    nTB(30) = 7
    if (.not. allocated(TB(30) % vector)) allocate(TB(30) % vector(7))
    if (.not. allocated(TBid(30) % vector)) allocate(TBid(30) % vector(7))
    TBid(30) % vector(1) = 0d0
    TB(30) % vector(1) = 2.3999999999999999d0 ! H2
    TBid(30) % vector(2) = 5d0
    TB(30) % vector(2) = 15.4d0 ! H2O
    TBid(30) % vector(3) = 13d0
    TB(30) % vector(3) = 2d0 ! CH4
    TBid(30) % vector(4) = 14d0
    TB(30) % vector(4) = 1.75d0 ! CO
    TBid(30) % vector(5) = 15d0
    TB(30) % vector(5) = 3.6000000000000001d0 ! CO2
    TBid(30) % vector(6) = 26d0
    TB(30) % vector(6) = 3d0 ! C2H6
    TBid(30) % vector(7) = 48d0
    TB(30) % vector(7) = 0.82999999999999996d0 ! AR

    ! (1):  O + H + M <=> OH + M
    fwd_A(31)     = 5d+17
    fwd_beta(31)  = -1d0
    fwd_Ea(31)    = 0d0
    prefactor_units(31)  = 1.0000000000000002d-12
    activation_units(31) = 0.50321666580471969d0
    phase_units(31)      = 1d-12
    is_PD(31) = 0
    nTB(31) = 7
    if (.not. allocated(TB(31) % vector)) allocate(TB(31) % vector(7))
    if (.not. allocated(TBid(31) % vector)) allocate(TBid(31) % vector(7))
    TBid(31) % vector(1) = 0d0
    TB(31) % vector(1) = 2d0 ! H2
    TBid(31) % vector(2) = 5d0
    TB(31) % vector(2) = 6d0 ! H2O
    TBid(31) % vector(3) = 13d0
    TB(31) % vector(3) = 2d0 ! CH4
    TBid(31) % vector(4) = 14d0
    TB(31) % vector(4) = 1.5d0 ! CO
    TBid(31) % vector(5) = 15d0
    TB(31) % vector(5) = 2d0 ! CO2
    TBid(31) % vector(6) = 26d0
    TB(31) % vector(6) = 3d0 ! C2H6
    TBid(31) % vector(7) = 48d0
    TB(31) % vector(7) = 0.69999999999999996d0 ! AR

    ! (2):  O + H2 <=> H + OH
    fwd_A(42)     = 38700d0
    fwd_beta(42)  = 2.7000000000000002d0
    fwd_Ea(42)    = 6260d0
    prefactor_units(42)  = 1.0000000000000002d-06
    activation_units(42) = 0.50321666580471969d0
    phase_units(42)      = 1d-12
    is_PD(42) = 0
    nTB(42) = 0

    ! (3):  O + HO2 <=> OH + O2
    fwd_A(43)     = 20000000000000d0
    fwd_beta(43)  = 0d0
    fwd_Ea(43)    = 0d0
    prefactor_units(43)  = 1.0000000000000002d-06
    activation_units(43) = 0.50321666580471969d0
    phase_units(43)      = 1d-12
    is_PD(43) = 0
    nTB(43) = 0

    ! (4):  O + H2O2 <=> OH + HO2
    fwd_A(44)     = 9630000d0
    fwd_beta(44)  = 2d0
    fwd_Ea(44)    = 4000d0
    prefactor_units(44)  = 1.0000000000000002d-06
    activation_units(44) = 0.50321666580471969d0
    phase_units(44)      = 1d-12
    is_PD(44) = 0
    nTB(44) = 0

    ! (5):  O + CH <=> H + CO
    fwd_A(45)     = 57000000000000d0
    fwd_beta(45)  = 0d0
    fwd_Ea(45)    = 0d0
    prefactor_units(45)  = 1.0000000000000002d-06
    activation_units(45) = 0.50321666580471969d0
    phase_units(45)      = 1d-12
    is_PD(45) = 0
    nTB(45) = 0

    ! (6):  O + CH2 <=> H + HCO
    fwd_A(46)     = 80000000000000d0
    fwd_beta(46)  = 0d0
    fwd_Ea(46)    = 0d0
    prefactor_units(46)  = 1.0000000000000002d-06
    activation_units(46) = 0.50321666580471969d0
    phase_units(46)      = 1d-12
    is_PD(46) = 0
    nTB(46) = 0

    ! (7):  O + CH2(S) <=> H2 + CO
    fwd_A(47)     = 15000000000000d0
    fwd_beta(47)  = 0d0
    fwd_Ea(47)    = 0d0
    prefactor_units(47)  = 1.0000000000000002d-06
    activation_units(47) = 0.50321666580471969d0
    phase_units(47)      = 1d-12
    is_PD(47) = 0
    nTB(47) = 0

    ! (8):  O + CH2(S) <=> H + HCO
    fwd_A(48)     = 15000000000000d0
    fwd_beta(48)  = 0d0
    fwd_Ea(48)    = 0d0
    prefactor_units(48)  = 1.0000000000000002d-06
    activation_units(48) = 0.50321666580471969d0
    phase_units(48)      = 1d-12
    is_PD(48) = 0
    nTB(48) = 0

    ! (9):  O + CH3 <=> H + CH2O
    fwd_A(49)     = 50600000000000d0
    fwd_beta(49)  = 0d0
    fwd_Ea(49)    = 0d0
    prefactor_units(49)  = 1.0000000000000002d-06
    activation_units(49) = 0.50321666580471969d0
    phase_units(49)      = 1d-12
    is_PD(49) = 0
    nTB(49) = 0

    ! (10):  O + CH4 <=> OH + CH3
    fwd_A(50)     = 1020000000d0
    fwd_beta(50)  = 1.5d0
    fwd_Ea(50)    = 8600d0
    prefactor_units(50)  = 1.0000000000000002d-06
    activation_units(50) = 0.50321666580471969d0
    phase_units(50)      = 1d-12
    is_PD(50) = 0
    nTB(50) = 0

    ! (11):  O + CO (+M) <=> CO2 (+M)
    fwd_A(27)     = 18000000000d0
    fwd_beta(27)  = 0d0
    fwd_Ea(27)    = 2385d0
    low_A(27)     = 602000000000000d0
    low_beta(27)  = 0d0
    low_Ea(27)    = 3000d0
    prefactor_units(27)  = 1.0000000000000002d-06
    activation_units(27) = 0.50321666580471969d0
    phase_units(27)      = 1d-12
    is_PD(27) = 1
    nTB(27) = 8
    if (.not. allocated(TB(27) % vector)) allocate(TB(27) % vector(8))
    if (.not. allocated(TBid(27) % vector)) allocate(TBid(27) % vector(8))
    TBid(27) % vector(1) = 0d0
    TB(27) % vector(1) = 2d0 ! H2
    TBid(27) % vector(2) = 3d0
    TB(27) % vector(2) = 6d0 ! O2
    TBid(27) % vector(3) = 5d0
    TB(27) % vector(3) = 6d0 ! H2O
    TBid(27) % vector(4) = 13d0
    TB(27) % vector(4) = 2d0 ! CH4
    TBid(27) % vector(5) = 14d0
    TB(27) % vector(5) = 1.5d0 ! CO
    TBid(27) % vector(6) = 15d0
    TB(27) % vector(6) = 3.5d0 ! CO2
    TBid(27) % vector(7) = 26d0
    TB(27) % vector(7) = 3d0 ! C2H6
    TBid(27) % vector(8) = 48d0
    TB(27) % vector(8) = 0.5d0 ! AR

    ! (12):  O + HCO <=> OH + CO
    fwd_A(51)     = 30000000000000d0
    fwd_beta(51)  = 0d0
    fwd_Ea(51)    = 0d0
    prefactor_units(51)  = 1.0000000000000002d-06
    activation_units(51) = 0.50321666580471969d0
    phase_units(51)      = 1d-12
    is_PD(51) = 0
    nTB(51) = 0

    ! (13):  O + HCO <=> H + CO2
    fwd_A(52)     = 30000000000000d0
    fwd_beta(52)  = 0d0
    fwd_Ea(52)    = 0d0
    prefactor_units(52)  = 1.0000000000000002d-06
    activation_units(52) = 0.50321666580471969d0
    phase_units(52)      = 1d-12
    is_PD(52) = 0
    nTB(52) = 0

    ! (14):  O + CH2O <=> OH + HCO
    fwd_A(53)     = 39000000000000d0
    fwd_beta(53)  = 0d0
    fwd_Ea(53)    = 3540d0
    prefactor_units(53)  = 1.0000000000000002d-06
    activation_units(53) = 0.50321666580471969d0
    phase_units(53)      = 1d-12
    is_PD(53) = 0
    nTB(53) = 0

    ! (15):  O + CH2OH <=> OH + CH2O
    fwd_A(54)     = 10000000000000d0
    fwd_beta(54)  = 0d0
    fwd_Ea(54)    = 0d0
    prefactor_units(54)  = 1.0000000000000002d-06
    activation_units(54) = 0.50321666580471969d0
    phase_units(54)      = 1d-12
    is_PD(54) = 0
    nTB(54) = 0

    ! (16):  O + CH3O <=> OH + CH2O
    fwd_A(55)     = 10000000000000d0
    fwd_beta(55)  = 0d0
    fwd_Ea(55)    = 0d0
    prefactor_units(55)  = 1.0000000000000002d-06
    activation_units(55) = 0.50321666580471969d0
    phase_units(55)      = 1d-12
    is_PD(55) = 0
    nTB(55) = 0

    ! (17):  O + CH3OH <=> OH + CH2OH
    fwd_A(56)     = 388000d0
    fwd_beta(56)  = 2.5d0
    fwd_Ea(56)    = 3100d0
    prefactor_units(56)  = 1.0000000000000002d-06
    activation_units(56) = 0.50321666580471969d0
    phase_units(56)      = 1d-12
    is_PD(56) = 0
    nTB(56) = 0

    ! (18):  O + CH3OH <=> OH + CH3O
    fwd_A(57)     = 130000d0
    fwd_beta(57)  = 2.5d0
    fwd_Ea(57)    = 5000d0
    prefactor_units(57)  = 1.0000000000000002d-06
    activation_units(57) = 0.50321666580471969d0
    phase_units(57)      = 1d-12
    is_PD(57) = 0
    nTB(57) = 0

    ! (19):  O + C2H <=> CH + CO
    fwd_A(58)     = 50000000000000d0
    fwd_beta(58)  = 0d0
    fwd_Ea(58)    = 0d0
    prefactor_units(58)  = 1.0000000000000002d-06
    activation_units(58) = 0.50321666580471969d0
    phase_units(58)      = 1d-12
    is_PD(58) = 0
    nTB(58) = 0

    ! (20):  O + C2H2 <=> H + HCCO
    fwd_A(59)     = 13500000d0
    fwd_beta(59)  = 2d0
    fwd_Ea(59)    = 1900d0
    prefactor_units(59)  = 1.0000000000000002d-06
    activation_units(59) = 0.50321666580471969d0
    phase_units(59)      = 1d-12
    is_PD(59) = 0
    nTB(59) = 0

    ! (21):  O + C2H2 <=> OH + C2H
    fwd_A(60)     = 4.6d+19
    fwd_beta(60)  = -1.4099999999999999d0
    fwd_Ea(60)    = 28950d0
    prefactor_units(60)  = 1.0000000000000002d-06
    activation_units(60) = 0.50321666580471969d0
    phase_units(60)      = 1d-12
    is_PD(60) = 0
    nTB(60) = 0

    ! (22):  O + C2H2 <=> CO + CH2
    fwd_A(61)     = 6940000d0
    fwd_beta(61)  = 2d0
    fwd_Ea(61)    = 1900d0
    prefactor_units(61)  = 1.0000000000000002d-06
    activation_units(61) = 0.50321666580471969d0
    phase_units(61)      = 1d-12
    is_PD(61) = 0
    nTB(61) = 0

    ! (23):  O + C2H3 <=> H + CH2CO
    fwd_A(62)     = 30000000000000d0
    fwd_beta(62)  = 0d0
    fwd_Ea(62)    = 0d0
    prefactor_units(62)  = 1.0000000000000002d-06
    activation_units(62) = 0.50321666580471969d0
    phase_units(62)      = 1d-12
    is_PD(62) = 0
    nTB(62) = 0

    ! (24):  O + C2H4 <=> CH3 + HCO
    fwd_A(63)     = 12500000d0
    fwd_beta(63)  = 1.8300000000000001d0
    fwd_Ea(63)    = 220d0
    prefactor_units(63)  = 1.0000000000000002d-06
    activation_units(63) = 0.50321666580471969d0
    phase_units(63)      = 1d-12
    is_PD(63) = 0
    nTB(63) = 0

    ! (25):  O + C2H5 <=> CH3 + CH2O
    fwd_A(64)     = 22400000000000d0
    fwd_beta(64)  = 0d0
    fwd_Ea(64)    = 0d0
    prefactor_units(64)  = 1.0000000000000002d-06
    activation_units(64) = 0.50321666580471969d0
    phase_units(64)      = 1d-12
    is_PD(64) = 0
    nTB(64) = 0

    ! (26):  O + C2H6 <=> OH + C2H5
    fwd_A(65)     = 89800000d0
    fwd_beta(65)  = 1.9199999999999999d0
    fwd_Ea(65)    = 5690d0
    prefactor_units(65)  = 1.0000000000000002d-06
    activation_units(65) = 0.50321666580471969d0
    phase_units(65)      = 1d-12
    is_PD(65) = 0
    nTB(65) = 0

    ! (27):  O + HCCO <=> H + 2 CO
    fwd_A(66)     = 100000000000000d0
    fwd_beta(66)  = 0d0
    fwd_Ea(66)    = 0d0
    prefactor_units(66)  = 1.0000000000000002d-06
    activation_units(66) = 0.50321666580471969d0
    phase_units(66)      = 1d-12
    is_PD(66) = 0
    nTB(66) = 0

    ! (28):  O + CH2CO <=> OH + HCCO
    fwd_A(67)     = 10000000000000d0
    fwd_beta(67)  = 0d0
    fwd_Ea(67)    = 8000d0
    prefactor_units(67)  = 1.0000000000000002d-06
    activation_units(67) = 0.50321666580471969d0
    phase_units(67)      = 1d-12
    is_PD(67) = 0
    nTB(67) = 0

    ! (29):  O + CH2CO <=> CH2 + CO2
    fwd_A(68)     = 1750000000000d0
    fwd_beta(68)  = 0d0
    fwd_Ea(68)    = 1350d0
    prefactor_units(68)  = 1.0000000000000002d-06
    activation_units(68) = 0.50321666580471969d0
    phase_units(68)      = 1d-12
    is_PD(68) = 0
    nTB(68) = 0

    ! (30):  O2 + CO <=> O + CO2
    fwd_A(69)     = 2500000000000d0
    fwd_beta(69)  = 0d0
    fwd_Ea(69)    = 47800d0
    prefactor_units(69)  = 1.0000000000000002d-06
    activation_units(69) = 0.50321666580471969d0
    phase_units(69)      = 1d-12
    is_PD(69) = 0
    nTB(69) = 0

    ! (31):  O2 + CH2O <=> HO2 + HCO
    fwd_A(70)     = 100000000000000d0
    fwd_beta(70)  = 0d0
    fwd_Ea(70)    = 40000d0
    prefactor_units(70)  = 1.0000000000000002d-06
    activation_units(70) = 0.50321666580471969d0
    phase_units(70)      = 1d-12
    is_PD(70) = 0
    nTB(70) = 0

    ! (32):  H + O2 + M <=> HO2 + M
    fwd_A(32)     = 2.8d+18
    fwd_beta(32)  = -0.85999999999999999d0
    fwd_Ea(32)    = 0d0
    prefactor_units(32)  = 1.0000000000000002d-12
    activation_units(32) = 0.50321666580471969d0
    phase_units(32)      = 1d-12
    is_PD(32) = 0
    nTB(32) = 7
    if (.not. allocated(TB(32) % vector)) allocate(TB(32) % vector(7))
    if (.not. allocated(TBid(32) % vector)) allocate(TBid(32) % vector(7))
    TBid(32) % vector(1) = 3d0
    TB(32) % vector(1) = 0d0 ! O2
    TBid(32) % vector(2) = 5d0
    TB(32) % vector(2) = 0d0 ! H2O
    TBid(32) % vector(3) = 14d0
    TB(32) % vector(3) = 0.75d0 ! CO
    TBid(32) % vector(4) = 15d0
    TB(32) % vector(4) = 1.5d0 ! CO2
    TBid(32) % vector(5) = 26d0
    TB(32) % vector(5) = 1.5d0 ! C2H6
    TBid(32) % vector(6) = 47d0
    TB(32) % vector(6) = 0d0 ! N2
    TBid(32) % vector(7) = 48d0
    TB(32) % vector(7) = 0d0 ! AR

    ! (33):  H + 2 O2 <=> HO2 + O2
    fwd_A(71)     = 2.08d+19
    fwd_beta(71)  = -1.24d0
    fwd_Ea(71)    = 0d0
    prefactor_units(71)  = 1.0000000000000002d-12
    activation_units(71) = 0.50321666580471969d0
    phase_units(71)      = 1d-18
    is_PD(71) = 0
    nTB(71) = 0

    ! (34):  H + O2 + H2O <=> HO2 + H2O
    fwd_A(72)     = 1.126d+19
    fwd_beta(72)  = -0.76000000000000001d0
    fwd_Ea(72)    = 0d0
    prefactor_units(72)  = 1.0000000000000002d-12
    activation_units(72) = 0.50321666580471969d0
    phase_units(72)      = 1d-18
    is_PD(72) = 0
    nTB(72) = 0

    ! (35):  H + O2 + N2 <=> HO2 + N2
    fwd_A(73)     = 2.6d+19
    fwd_beta(73)  = -1.24d0
    fwd_Ea(73)    = 0d0
    prefactor_units(73)  = 1.0000000000000002d-12
    activation_units(73) = 0.50321666580471969d0
    phase_units(73)      = 1d-18
    is_PD(73) = 0
    nTB(73) = 0

    ! (36):  H + O2 + AR <=> HO2 + AR
    fwd_A(74)     = 7d+17
    fwd_beta(74)  = -0.80000000000000004d0
    fwd_Ea(74)    = 0d0
    prefactor_units(74)  = 1.0000000000000002d-12
    activation_units(74) = 0.50321666580471969d0
    phase_units(74)      = 1d-18
    is_PD(74) = 0
    nTB(74) = 0

    ! (37):  H + O2 <=> O + OH
    fwd_A(75)     = 26500000000000000d0
    fwd_beta(75)  = -0.67069999999999996d0
    fwd_Ea(75)    = 17041d0
    prefactor_units(75)  = 1.0000000000000002d-06
    activation_units(75) = 0.50321666580471969d0
    phase_units(75)      = 1d-12
    is_PD(75) = 0
    nTB(75) = 0

    ! (38):  2 H + M <=> H2 + M
    fwd_A(33)     = 1d+18
    fwd_beta(33)  = -1d0
    fwd_Ea(33)    = 0d0
    prefactor_units(33)  = 1.0000000000000002d-12
    activation_units(33) = 0.50321666580471969d0
    phase_units(33)      = 1d-12
    is_PD(33) = 0
    nTB(33) = 6
    if (.not. allocated(TB(33) % vector)) allocate(TB(33) % vector(6))
    if (.not. allocated(TBid(33) % vector)) allocate(TBid(33) % vector(6))
    TBid(33) % vector(1) = 0d0
    TB(33) % vector(1) = 0d0 ! H2
    TBid(33) % vector(2) = 5d0
    TB(33) % vector(2) = 0d0 ! H2O
    TBid(33) % vector(3) = 13d0
    TB(33) % vector(3) = 2d0 ! CH4
    TBid(33) % vector(4) = 15d0
    TB(33) % vector(4) = 0d0 ! CO2
    TBid(33) % vector(5) = 26d0
    TB(33) % vector(5) = 3d0 ! C2H6
    TBid(33) % vector(6) = 48d0
    TB(33) % vector(6) = 0.63d0 ! AR

    ! (39):  2 H + H2 <=> 2 H2
    fwd_A(76)     = 90000000000000000d0
    fwd_beta(76)  = -0.59999999999999998d0
    fwd_Ea(76)    = 0d0
    prefactor_units(76)  = 1.0000000000000002d-12
    activation_units(76) = 0.50321666580471969d0
    phase_units(76)      = 1d-18
    is_PD(76) = 0
    nTB(76) = 0

    ! (40):  2 H + H2O <=> H2 + H2O
    fwd_A(77)     = 6d+19
    fwd_beta(77)  = -1.25d0
    fwd_Ea(77)    = 0d0
    prefactor_units(77)  = 1.0000000000000002d-12
    activation_units(77) = 0.50321666580471969d0
    phase_units(77)      = 1d-18
    is_PD(77) = 0
    nTB(77) = 0

    ! (41):  2 H + CO2 <=> H2 + CO2
    fwd_A(78)     = 5.5d+20
    fwd_beta(78)  = -2d0
    fwd_Ea(78)    = 0d0
    prefactor_units(78)  = 1.0000000000000002d-12
    activation_units(78) = 0.50321666580471969d0
    phase_units(78)      = 1d-18
    is_PD(78) = 0
    nTB(78) = 0

    ! (42):  H + OH + M <=> H2O + M
    fwd_A(34)     = 2.2d+22
    fwd_beta(34)  = -2d0
    fwd_Ea(34)    = 0d0
    prefactor_units(34)  = 1.0000000000000002d-12
    activation_units(34) = 0.50321666580471969d0
    phase_units(34)      = 1d-12
    is_PD(34) = 0
    nTB(34) = 5
    if (.not. allocated(TB(34) % vector)) allocate(TB(34) % vector(5))
    if (.not. allocated(TBid(34) % vector)) allocate(TBid(34) % vector(5))
    TBid(34) % vector(1) = 0d0
    TB(34) % vector(1) = 0.72999999999999998d0 ! H2
    TBid(34) % vector(2) = 5d0
    TB(34) % vector(2) = 3.6499999999999999d0 ! H2O
    TBid(34) % vector(3) = 13d0
    TB(34) % vector(3) = 2d0 ! CH4
    TBid(34) % vector(4) = 26d0
    TB(34) % vector(4) = 3d0 ! C2H6
    TBid(34) % vector(5) = 48d0
    TB(34) % vector(5) = 0.38d0 ! AR

    ! (43):  H + HO2 <=> O + H2O
    fwd_A(79)     = 3970000000000d0
    fwd_beta(79)  = 0d0
    fwd_Ea(79)    = 671d0
    prefactor_units(79)  = 1.0000000000000002d-06
    activation_units(79) = 0.50321666580471969d0
    phase_units(79)      = 1d-12
    is_PD(79) = 0
    nTB(79) = 0

    ! (44):  H + HO2 <=> O2 + H2
    fwd_A(80)     = 44800000000000d0
    fwd_beta(80)  = 0d0
    fwd_Ea(80)    = 1068d0
    prefactor_units(80)  = 1.0000000000000002d-06
    activation_units(80) = 0.50321666580471969d0
    phase_units(80)      = 1d-12
    is_PD(80) = 0
    nTB(80) = 0

    ! (45):  H + HO2 <=> 2 OH
    fwd_A(81)     = 84000000000000d0
    fwd_beta(81)  = 0d0
    fwd_Ea(81)    = 635d0
    prefactor_units(81)  = 1.0000000000000002d-06
    activation_units(81) = 0.50321666580471969d0
    phase_units(81)      = 1d-12
    is_PD(81) = 0
    nTB(81) = 0

    ! (46):  H + H2O2 <=> HO2 + H2
    fwd_A(82)     = 12100000d0
    fwd_beta(82)  = 2d0
    fwd_Ea(82)    = 5200d0
    prefactor_units(82)  = 1.0000000000000002d-06
    activation_units(82) = 0.50321666580471969d0
    phase_units(82)      = 1d-12
    is_PD(82) = 0
    nTB(82) = 0

    ! (47):  H + H2O2 <=> OH + H2O
    fwd_A(83)     = 10000000000000d0
    fwd_beta(83)  = 0d0
    fwd_Ea(83)    = 3600d0
    prefactor_units(83)  = 1.0000000000000002d-06
    activation_units(83) = 0.50321666580471969d0
    phase_units(83)      = 1d-12
    is_PD(83) = 0
    nTB(83) = 0

    ! (48):  H + CH <=> C + H2
    fwd_A(84)     = 165000000000000d0
    fwd_beta(84)  = 0d0
    fwd_Ea(84)    = 0d0
    prefactor_units(84)  = 1.0000000000000002d-06
    activation_units(84) = 0.50321666580471969d0
    phase_units(84)      = 1d-12
    is_PD(84) = 0
    nTB(84) = 0

    ! (49):  H + CH2 (+M) <=> CH3 (+M)
    fwd_A(1)     = 600000000000000d0
    fwd_beta(1)  = 0d0
    fwd_Ea(1)    = 0d0
    low_A(1)     = 1.0399999999999999d+26
    low_beta(1)  = -2.7599999999999998d0
    low_Ea(1)    = 1600d0
    troe_a(1)    = 0.56200000000000006d0
    troe_Tsss(1) = 91d0
    troe_Ts(1)   = 5836d0
    troe_Tss(1)  = 8552d0
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
    TBid(1) % vector(3) = 13d0
    TB(1) % vector(3) = 2d0 ! CH4
    TBid(1) % vector(4) = 14d0
    TB(1) % vector(4) = 1.5d0 ! CO
    TBid(1) % vector(5) = 15d0
    TB(1) % vector(5) = 2d0 ! CO2
    TBid(1) % vector(6) = 26d0
    TB(1) % vector(6) = 3d0 ! C2H6
    TBid(1) % vector(7) = 48d0
    TB(1) % vector(7) = 0.69999999999999996d0 ! AR

    ! (50):  H + CH2(S) <=> CH + H2
    fwd_A(85)     = 30000000000000d0
    fwd_beta(85)  = 0d0
    fwd_Ea(85)    = 0d0
    prefactor_units(85)  = 1.0000000000000002d-06
    activation_units(85) = 0.50321666580471969d0
    phase_units(85)      = 1d-12
    is_PD(85) = 0
    nTB(85) = 0

    ! (51):  H + CH3 (+M) <=> CH4 (+M)
    fwd_A(2)     = 13900000000000000d0
    fwd_beta(2)  = -0.53400000000000003d0
    fwd_Ea(2)    = 536d0
    low_A(2)     = 2.6200000000000002d+33
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
    TBid(2) % vector(3) = 13d0
    TB(2) % vector(3) = 3d0 ! CH4
    TBid(2) % vector(4) = 14d0
    TB(2) % vector(4) = 1.5d0 ! CO
    TBid(2) % vector(5) = 15d0
    TB(2) % vector(5) = 2d0 ! CO2
    TBid(2) % vector(6) = 26d0
    TB(2) % vector(6) = 3d0 ! C2H6
    TBid(2) % vector(7) = 48d0
    TB(2) % vector(7) = 0.69999999999999996d0 ! AR

    ! (52):  H + CH4 <=> CH3 + H2
    fwd_A(86)     = 660000000d0
    fwd_beta(86)  = 1.6200000000000001d0
    fwd_Ea(86)    = 10840d0
    prefactor_units(86)  = 1.0000000000000002d-06
    activation_units(86) = 0.50321666580471969d0
    phase_units(86)      = 1d-12
    is_PD(86) = 0
    nTB(86) = 0

    ! (53):  H + HCO (+M) <=> CH2O (+M)
    fwd_A(3)     = 1090000000000d0
    fwd_beta(3)  = 0.47999999999999998d0
    fwd_Ea(3)    = -260d0
    low_A(3)     = 2.4700000000000001d+24
    low_beta(3)  = -2.5699999999999998d0
    low_Ea(3)    = 425d0
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
    TBid(3) % vector(3) = 13d0
    TB(3) % vector(3) = 2d0 ! CH4
    TBid(3) % vector(4) = 14d0
    TB(3) % vector(4) = 1.5d0 ! CO
    TBid(3) % vector(5) = 15d0
    TB(3) % vector(5) = 2d0 ! CO2
    TBid(3) % vector(6) = 26d0
    TB(3) % vector(6) = 3d0 ! C2H6
    TBid(3) % vector(7) = 48d0
    TB(3) % vector(7) = 0.69999999999999996d0 ! AR

    ! (54):  H + HCO <=> H2 + CO
    fwd_A(87)     = 73400000000000d0
    fwd_beta(87)  = 0d0
    fwd_Ea(87)    = 0d0
    prefactor_units(87)  = 1.0000000000000002d-06
    activation_units(87) = 0.50321666580471969d0
    phase_units(87)      = 1d-12
    is_PD(87) = 0
    nTB(87) = 0

    ! (55):  H + CH2O (+M) <=> CH2OH (+M)
    fwd_A(4)     = 540000000000d0
    fwd_beta(4)  = 0.45400000000000001d0
    fwd_Ea(4)    = 3600d0
    low_A(4)     = 1.27d+32
    low_beta(4)  = -4.8200000000000003d0
    low_Ea(4)    = 6530d0
    troe_a(4)    = 0.71870000000000001d0
    troe_Tsss(4) = 103d0
    troe_Ts(4)   = 1291d0
    troe_Tss(4)  = 4160d0
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
    TBid(4) % vector(3) = 13d0
    TB(4) % vector(3) = 2d0 ! CH4
    TBid(4) % vector(4) = 14d0
    TB(4) % vector(4) = 1.5d0 ! CO
    TBid(4) % vector(5) = 15d0
    TB(4) % vector(5) = 2d0 ! CO2
    TBid(4) % vector(6) = 26d0
    TB(4) % vector(6) = 3d0 ! C2H6

    ! (56):  H + CH2O (+M) <=> CH3O (+M)
    fwd_A(5)     = 540000000000d0
    fwd_beta(5)  = 0.45400000000000001d0
    fwd_Ea(5)    = 2600d0
    low_A(5)     = 2.2d+30
    low_beta(5)  = -4.7999999999999998d0
    low_Ea(5)    = 5560d0
    troe_a(5)    = 0.75800000000000001d0
    troe_Tsss(5) = 94d0
    troe_Ts(5)   = 1555d0
    troe_Tss(5)  = 4200d0
    troe_len(5)  = 4
    prefactor_units(5)  = 1.0000000000000002d-06
    activation_units(5) = 0.50321666580471969d0
    phase_units(5)      = 1d-12
    is_PD(5) = 1
    nTB(5) = 6
    if (.not. allocated(TB(5) % vector)) allocate(TB(5) % vector(6))
    if (.not. allocated(TBid(5) % vector)) allocate(TBid(5) % vector(6))
    TBid(5) % vector(1) = 0d0
    TB(5) % vector(1) = 2d0 ! H2
    TBid(5) % vector(2) = 5d0
    TB(5) % vector(2) = 6d0 ! H2O
    TBid(5) % vector(3) = 13d0
    TB(5) % vector(3) = 2d0 ! CH4
    TBid(5) % vector(4) = 14d0
    TB(5) % vector(4) = 1.5d0 ! CO
    TBid(5) % vector(5) = 15d0
    TB(5) % vector(5) = 2d0 ! CO2
    TBid(5) % vector(6) = 26d0
    TB(5) % vector(6) = 3d0 ! C2H6

    ! (57):  H + CH2O <=> HCO + H2
    fwd_A(88)     = 57400000d0
    fwd_beta(88)  = 1.8999999999999999d0
    fwd_Ea(88)    = 2742d0
    prefactor_units(88)  = 1.0000000000000002d-06
    activation_units(88) = 0.50321666580471969d0
    phase_units(88)      = 1d-12
    is_PD(88) = 0
    nTB(88) = 0

    ! (58):  H + CH2OH (+M) <=> CH3OH (+M)
    fwd_A(6)     = 1055000000000d0
    fwd_beta(6)  = 0.5d0
    fwd_Ea(6)    = 86d0
    low_A(6)     = 4.3600000000000004d+31
    low_beta(6)  = -4.6500000000000004d0
    low_Ea(6)    = 5080d0
    troe_a(6)    = 0.59999999999999998d0
    troe_Tsss(6) = 100d0
    troe_Ts(6)   = 90000d0
    troe_Tss(6)  = 10000d0
    troe_len(6)  = 4
    prefactor_units(6)  = 1.0000000000000002d-06
    activation_units(6) = 0.50321666580471969d0
    phase_units(6)      = 1d-12
    is_PD(6) = 1
    nTB(6) = 6
    if (.not. allocated(TB(6) % vector)) allocate(TB(6) % vector(6))
    if (.not. allocated(TBid(6) % vector)) allocate(TBid(6) % vector(6))
    TBid(6) % vector(1) = 0d0
    TB(6) % vector(1) = 2d0 ! H2
    TBid(6) % vector(2) = 5d0
    TB(6) % vector(2) = 6d0 ! H2O
    TBid(6) % vector(3) = 13d0
    TB(6) % vector(3) = 2d0 ! CH4
    TBid(6) % vector(4) = 14d0
    TB(6) % vector(4) = 1.5d0 ! CO
    TBid(6) % vector(5) = 15d0
    TB(6) % vector(5) = 2d0 ! CO2
    TBid(6) % vector(6) = 26d0
    TB(6) % vector(6) = 3d0 ! C2H6

    ! (59):  H + CH2OH <=> H2 + CH2O
    fwd_A(89)     = 20000000000000d0
    fwd_beta(89)  = 0d0
    fwd_Ea(89)    = 0d0
    prefactor_units(89)  = 1.0000000000000002d-06
    activation_units(89) = 0.50321666580471969d0
    phase_units(89)      = 1d-12
    is_PD(89) = 0
    nTB(89) = 0

    ! (60):  H + CH2OH <=> OH + CH3
    fwd_A(90)     = 165000000000d0
    fwd_beta(90)  = 0.65000000000000002d0
    fwd_Ea(90)    = -284d0
    prefactor_units(90)  = 1.0000000000000002d-06
    activation_units(90) = 0.50321666580471969d0
    phase_units(90)      = 1d-12
    is_PD(90) = 0
    nTB(90) = 0

    ! (61):  H + CH2OH <=> CH2(S) + H2O
    fwd_A(91)     = 32800000000000d0
    fwd_beta(91)  = -0.089999999999999997d0
    fwd_Ea(91)    = 610d0
    prefactor_units(91)  = 1.0000000000000002d-06
    activation_units(91) = 0.50321666580471969d0
    phase_units(91)      = 1d-12
    is_PD(91) = 0
    nTB(91) = 0

    ! (62):  H + CH3O (+M) <=> CH3OH (+M)
    fwd_A(7)     = 2430000000000d0
    fwd_beta(7)  = 0.51500000000000001d0
    fwd_Ea(7)    = 50d0
    low_A(7)     = 4.6600000000000002d+41
    low_beta(7)  = -7.4400000000000004d0
    low_Ea(7)    = 14080d0
    troe_a(7)    = 0.69999999999999996d0
    troe_Tsss(7) = 100d0
    troe_Ts(7)   = 90000d0
    troe_Tss(7)  = 10000d0
    troe_len(7)  = 4
    prefactor_units(7)  = 1.0000000000000002d-06
    activation_units(7) = 0.50321666580471969d0
    phase_units(7)      = 1d-12
    is_PD(7) = 1
    nTB(7) = 6
    if (.not. allocated(TB(7) % vector)) allocate(TB(7) % vector(6))
    if (.not. allocated(TBid(7) % vector)) allocate(TBid(7) % vector(6))
    TBid(7) % vector(1) = 0d0
    TB(7) % vector(1) = 2d0 ! H2
    TBid(7) % vector(2) = 5d0
    TB(7) % vector(2) = 6d0 ! H2O
    TBid(7) % vector(3) = 13d0
    TB(7) % vector(3) = 2d0 ! CH4
    TBid(7) % vector(4) = 14d0
    TB(7) % vector(4) = 1.5d0 ! CO
    TBid(7) % vector(5) = 15d0
    TB(7) % vector(5) = 2d0 ! CO2
    TBid(7) % vector(6) = 26d0
    TB(7) % vector(6) = 3d0 ! C2H6

    ! (63):  H + CH3O <=> H + CH2OH
    fwd_A(92)     = 41500000d0
    fwd_beta(92)  = 1.6299999999999999d0
    fwd_Ea(92)    = 1924d0
    prefactor_units(92)  = 1.0000000000000002d-06
    activation_units(92) = 0.50321666580471969d0
    phase_units(92)      = 1d-12
    is_PD(92) = 0
    nTB(92) = 0

    ! (64):  H + CH3O <=> H2 + CH2O
    fwd_A(93)     = 20000000000000d0
    fwd_beta(93)  = 0d0
    fwd_Ea(93)    = 0d0
    prefactor_units(93)  = 1.0000000000000002d-06
    activation_units(93) = 0.50321666580471969d0
    phase_units(93)      = 1d-12
    is_PD(93) = 0
    nTB(93) = 0

    ! (65):  H + CH3O <=> OH + CH3
    fwd_A(94)     = 1500000000000d0
    fwd_beta(94)  = 0.5d0
    fwd_Ea(94)    = -110d0
    prefactor_units(94)  = 1.0000000000000002d-06
    activation_units(94) = 0.50321666580471969d0
    phase_units(94)      = 1d-12
    is_PD(94) = 0
    nTB(94) = 0

    ! (66):  H + CH3O <=> CH2(S) + H2O
    fwd_A(95)     = 262000000000000d0
    fwd_beta(95)  = -0.23000000000000001d0
    fwd_Ea(95)    = 1070d0
    prefactor_units(95)  = 1.0000000000000002d-06
    activation_units(95) = 0.50321666580471969d0
    phase_units(95)      = 1d-12
    is_PD(95) = 0
    nTB(95) = 0

    ! (67):  H + CH3OH <=> CH2OH + H2
    fwd_A(96)     = 17000000d0
    fwd_beta(96)  = 2.1000000000000001d0
    fwd_Ea(96)    = 4870d0
    prefactor_units(96)  = 1.0000000000000002d-06
    activation_units(96) = 0.50321666580471969d0
    phase_units(96)      = 1d-12
    is_PD(96) = 0
    nTB(96) = 0

    ! (68):  H + CH3OH <=> CH3O + H2
    fwd_A(97)     = 4200000d0
    fwd_beta(97)  = 2.1000000000000001d0
    fwd_Ea(97)    = 4870d0
    prefactor_units(97)  = 1.0000000000000002d-06
    activation_units(97) = 0.50321666580471969d0
    phase_units(97)      = 1d-12
    is_PD(97) = 0
    nTB(97) = 0

    ! (69):  H + C2H (+M) <=> C2H2 (+M)
    fwd_A(8)     = 1d+17
    fwd_beta(8)  = -1d0
    fwd_Ea(8)    = 0d0
    low_A(8)     = 3.7500000000000002d+33
    low_beta(8)  = -4.7999999999999998d0
    low_Ea(8)    = 1900d0
    troe_a(8)    = 0.64639999999999997d0
    troe_Tsss(8) = 132d0
    troe_Ts(8)   = 1315d0
    troe_Tss(8)  = 5566d0
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
    TBid(8) % vector(3) = 13d0
    TB(8) % vector(3) = 2d0 ! CH4
    TBid(8) % vector(4) = 14d0
    TB(8) % vector(4) = 1.5d0 ! CO
    TBid(8) % vector(5) = 15d0
    TB(8) % vector(5) = 2d0 ! CO2
    TBid(8) % vector(6) = 26d0
    TB(8) % vector(6) = 3d0 ! C2H6
    TBid(8) % vector(7) = 48d0
    TB(8) % vector(7) = 0.69999999999999996d0 ! AR

    ! (70):  H + C2H2 (+M) <=> C2H3 (+M)
    fwd_A(9)     = 5600000000000d0
    fwd_beta(9)  = 0d0
    fwd_Ea(9)    = 2400d0
    low_A(9)     = 3.8d+40
    low_beta(9)  = -7.2699999999999996d0
    low_Ea(9)    = 7220d0
    troe_a(9)    = 0.75070000000000003d0
    troe_Tsss(9) = 98.5d0
    troe_Ts(9)   = 1302d0
    troe_Tss(9)  = 4167d0
    troe_len(9)  = 4
    prefactor_units(9)  = 1.0000000000000002d-06
    activation_units(9) = 0.50321666580471969d0
    phase_units(9)      = 1d-12
    is_PD(9) = 1
    nTB(9) = 7
    if (.not. allocated(TB(9) % vector)) allocate(TB(9) % vector(7))
    if (.not. allocated(TBid(9) % vector)) allocate(TBid(9) % vector(7))
    TBid(9) % vector(1) = 0d0
    TB(9) % vector(1) = 2d0 ! H2
    TBid(9) % vector(2) = 5d0
    TB(9) % vector(2) = 6d0 ! H2O
    TBid(9) % vector(3) = 13d0
    TB(9) % vector(3) = 2d0 ! CH4
    TBid(9) % vector(4) = 14d0
    TB(9) % vector(4) = 1.5d0 ! CO
    TBid(9) % vector(5) = 15d0
    TB(9) % vector(5) = 2d0 ! CO2
    TBid(9) % vector(6) = 26d0
    TB(9) % vector(6) = 3d0 ! C2H6
    TBid(9) % vector(7) = 48d0
    TB(9) % vector(7) = 0.69999999999999996d0 ! AR

    ! (71):  H + C2H3 (+M) <=> C2H4 (+M)
    fwd_A(10)     = 6080000000000d0
    fwd_beta(10)  = 0.27000000000000002d0
    fwd_Ea(10)    = 280d0
    low_A(10)     = 1.3999999999999999d+30
    low_beta(10)  = -3.8599999999999999d0
    low_Ea(10)    = 3320d0
    troe_a(10)    = 0.78200000000000003d0
    troe_Tsss(10) = 207.5d0
    troe_Ts(10)   = 2663d0
    troe_Tss(10)  = 6095d0
    troe_len(10)  = 4
    prefactor_units(10)  = 1.0000000000000002d-06
    activation_units(10) = 0.50321666580471969d0
    phase_units(10)      = 1d-12
    is_PD(10) = 1
    nTB(10) = 7
    if (.not. allocated(TB(10) % vector)) allocate(TB(10) % vector(7))
    if (.not. allocated(TBid(10) % vector)) allocate(TBid(10) % vector(7))
    TBid(10) % vector(1) = 0d0
    TB(10) % vector(1) = 2d0 ! H2
    TBid(10) % vector(2) = 5d0
    TB(10) % vector(2) = 6d0 ! H2O
    TBid(10) % vector(3) = 13d0
    TB(10) % vector(3) = 2d0 ! CH4
    TBid(10) % vector(4) = 14d0
    TB(10) % vector(4) = 1.5d0 ! CO
    TBid(10) % vector(5) = 15d0
    TB(10) % vector(5) = 2d0 ! CO2
    TBid(10) % vector(6) = 26d0
    TB(10) % vector(6) = 3d0 ! C2H6
    TBid(10) % vector(7) = 48d0
    TB(10) % vector(7) = 0.69999999999999996d0 ! AR

    ! (72):  H + C2H3 <=> H2 + C2H2
    fwd_A(98)     = 30000000000000d0
    fwd_beta(98)  = 0d0
    fwd_Ea(98)    = 0d0
    prefactor_units(98)  = 1.0000000000000002d-06
    activation_units(98) = 0.50321666580471969d0
    phase_units(98)      = 1d-12
    is_PD(98) = 0
    nTB(98) = 0

    ! (73):  H + C2H4 (+M) <=> C2H5 (+M)
    fwd_A(11)     = 540000000000d0
    fwd_beta(11)  = 0.45400000000000001d0
    fwd_Ea(11)    = 1820d0
    low_A(11)     = 5.9999999999999997d+41
    low_beta(11)  = -7.6200000000000001d0
    low_Ea(11)    = 6970d0
    troe_a(11)    = 0.97529999999999994d0
    troe_Tsss(11) = 210d0
    troe_Ts(11)   = 984d0
    troe_Tss(11)  = 4374d0
    troe_len(11)  = 4
    prefactor_units(11)  = 1.0000000000000002d-06
    activation_units(11) = 0.50321666580471969d0
    phase_units(11)      = 1d-12
    is_PD(11) = 1
    nTB(11) = 7
    if (.not. allocated(TB(11) % vector)) allocate(TB(11) % vector(7))
    if (.not. allocated(TBid(11) % vector)) allocate(TBid(11) % vector(7))
    TBid(11) % vector(1) = 0d0
    TB(11) % vector(1) = 2d0 ! H2
    TBid(11) % vector(2) = 5d0
    TB(11) % vector(2) = 6d0 ! H2O
    TBid(11) % vector(3) = 13d0
    TB(11) % vector(3) = 2d0 ! CH4
    TBid(11) % vector(4) = 14d0
    TB(11) % vector(4) = 1.5d0 ! CO
    TBid(11) % vector(5) = 15d0
    TB(11) % vector(5) = 2d0 ! CO2
    TBid(11) % vector(6) = 26d0
    TB(11) % vector(6) = 3d0 ! C2H6
    TBid(11) % vector(7) = 48d0
    TB(11) % vector(7) = 0.69999999999999996d0 ! AR

    ! (74):  H + C2H4 <=> C2H3 + H2
    fwd_A(99)     = 1325000d0
    fwd_beta(99)  = 2.5299999999999998d0
    fwd_Ea(99)    = 12240d0
    prefactor_units(99)  = 1.0000000000000002d-06
    activation_units(99) = 0.50321666580471969d0
    phase_units(99)      = 1d-12
    is_PD(99) = 0
    nTB(99) = 0

    ! (75):  H + C2H5 (+M) <=> C2H6 (+M)
    fwd_A(12)     = 5.21d+17
    fwd_beta(12)  = -0.98999999999999999d0
    fwd_Ea(12)    = 1580d0
    low_A(12)     = 1.9900000000000001d+41
    low_beta(12)  = -7.0800000000000001d0
    low_Ea(12)    = 6685d0
    troe_a(12)    = 0.84219999999999995d0
    troe_Tsss(12) = 125d0
    troe_Ts(12)   = 2219d0
    troe_Tss(12)  = 6882d0
    troe_len(12)  = 4
    prefactor_units(12)  = 1.0000000000000002d-06
    activation_units(12) = 0.50321666580471969d0
    phase_units(12)      = 1d-12
    is_PD(12) = 1
    nTB(12) = 7
    if (.not. allocated(TB(12) % vector)) allocate(TB(12) % vector(7))
    if (.not. allocated(TBid(12) % vector)) allocate(TBid(12) % vector(7))
    TBid(12) % vector(1) = 0d0
    TB(12) % vector(1) = 2d0 ! H2
    TBid(12) % vector(2) = 5d0
    TB(12) % vector(2) = 6d0 ! H2O
    TBid(12) % vector(3) = 13d0
    TB(12) % vector(3) = 2d0 ! CH4
    TBid(12) % vector(4) = 14d0
    TB(12) % vector(4) = 1.5d0 ! CO
    TBid(12) % vector(5) = 15d0
    TB(12) % vector(5) = 2d0 ! CO2
    TBid(12) % vector(6) = 26d0
    TB(12) % vector(6) = 3d0 ! C2H6
    TBid(12) % vector(7) = 48d0
    TB(12) % vector(7) = 0.69999999999999996d0 ! AR

    ! (76):  H + C2H5 <=> H2 + C2H4
    fwd_A(100)     = 2000000000000d0
    fwd_beta(100)  = 0d0
    fwd_Ea(100)    = 0d0
    prefactor_units(100)  = 1.0000000000000002d-06
    activation_units(100) = 0.50321666580471969d0
    phase_units(100)      = 1d-12
    is_PD(100) = 0
    nTB(100) = 0

    ! (77):  H + C2H6 <=> C2H5 + H2
    fwd_A(101)     = 115000000d0
    fwd_beta(101)  = 1.8999999999999999d0
    fwd_Ea(101)    = 7530d0
    prefactor_units(101)  = 1.0000000000000002d-06
    activation_units(101) = 0.50321666580471969d0
    phase_units(101)      = 1d-12
    is_PD(101) = 0
    nTB(101) = 0

    ! (78):  H + HCCO <=> CH2(S) + CO
    fwd_A(102)     = 100000000000000d0
    fwd_beta(102)  = 0d0
    fwd_Ea(102)    = 0d0
    prefactor_units(102)  = 1.0000000000000002d-06
    activation_units(102) = 0.50321666580471969d0
    phase_units(102)      = 1d-12
    is_PD(102) = 0
    nTB(102) = 0

    ! (79):  H + CH2CO <=> HCCO + H2
    fwd_A(103)     = 50000000000000d0
    fwd_beta(103)  = 0d0
    fwd_Ea(103)    = 8000d0
    prefactor_units(103)  = 1.0000000000000002d-06
    activation_units(103) = 0.50321666580471969d0
    phase_units(103)      = 1d-12
    is_PD(103) = 0
    nTB(103) = 0

    ! (80):  H + CH2CO <=> CH3 + CO
    fwd_A(104)     = 11300000000000d0
    fwd_beta(104)  = 0d0
    fwd_Ea(104)    = 3428d0
    prefactor_units(104)  = 1.0000000000000002d-06
    activation_units(104) = 0.50321666580471969d0
    phase_units(104)      = 1d-12
    is_PD(104) = 0
    nTB(104) = 0

    ! (81):  H + HCCOH <=> H + CH2CO
    fwd_A(105)     = 10000000000000d0
    fwd_beta(105)  = 0d0
    fwd_Ea(105)    = 0d0
    prefactor_units(105)  = 1.0000000000000002d-06
    activation_units(105) = 0.50321666580471969d0
    phase_units(105)      = 1d-12
    is_PD(105) = 0
    nTB(105) = 0

    ! (82):  H2 + CO (+M) <=> CH2O (+M)
    fwd_A(13)     = 43000000d0
    fwd_beta(13)  = 1.5d0
    fwd_Ea(13)    = 79600d0
    low_A(13)     = 5.0699999999999998d+27
    low_beta(13)  = -3.4199999999999999d0
    low_Ea(13)    = 84350d0
    troe_a(13)    = 0.93200000000000005d0
    troe_Tsss(13) = 197d0
    troe_Ts(13)   = 1540d0
    troe_Tss(13)  = 10300d0
    troe_len(13)  = 4
    prefactor_units(13)  = 1.0000000000000002d-06
    activation_units(13) = 0.50321666580471969d0
    phase_units(13)      = 1d-12
    is_PD(13) = 1
    nTB(13) = 7
    if (.not. allocated(TB(13) % vector)) allocate(TB(13) % vector(7))
    if (.not. allocated(TBid(13) % vector)) allocate(TBid(13) % vector(7))
    TBid(13) % vector(1) = 0d0
    TB(13) % vector(1) = 2d0 ! H2
    TBid(13) % vector(2) = 5d0
    TB(13) % vector(2) = 6d0 ! H2O
    TBid(13) % vector(3) = 13d0
    TB(13) % vector(3) = 2d0 ! CH4
    TBid(13) % vector(4) = 14d0
    TB(13) % vector(4) = 1.5d0 ! CO
    TBid(13) % vector(5) = 15d0
    TB(13) % vector(5) = 2d0 ! CO2
    TBid(13) % vector(6) = 26d0
    TB(13) % vector(6) = 3d0 ! C2H6
    TBid(13) % vector(7) = 48d0
    TB(13) % vector(7) = 0.69999999999999996d0 ! AR

    ! (83):  OH + H2 <=> H + H2O
    fwd_A(106)     = 216000000d0
    fwd_beta(106)  = 1.51d0
    fwd_Ea(106)    = 3430d0
    prefactor_units(106)  = 1.0000000000000002d-06
    activation_units(106) = 0.50321666580471969d0
    phase_units(106)      = 1d-12
    is_PD(106) = 0
    nTB(106) = 0

    ! (84):  2 OH (+M) <=> H2O2 (+M)
    fwd_A(14)     = 74000000000000d0
    fwd_beta(14)  = -0.37d0
    fwd_Ea(14)    = 0d0
    low_A(14)     = 2.3d+18
    low_beta(14)  = -0.90000000000000002d0
    low_Ea(14)    = -1700d0
    troe_a(14)    = 0.73460000000000003d0
    troe_Tsss(14) = 94d0
    troe_Ts(14)   = 1756d0
    troe_Tss(14)  = 5182d0
    troe_len(14)  = 4
    prefactor_units(14)  = 1.0000000000000002d-06
    activation_units(14) = 0.50321666580471969d0
    phase_units(14)      = 1d-12
    is_PD(14) = 1
    nTB(14) = 7
    if (.not. allocated(TB(14) % vector)) allocate(TB(14) % vector(7))
    if (.not. allocated(TBid(14) % vector)) allocate(TBid(14) % vector(7))
    TBid(14) % vector(1) = 0d0
    TB(14) % vector(1) = 2d0 ! H2
    TBid(14) % vector(2) = 5d0
    TB(14) % vector(2) = 6d0 ! H2O
    TBid(14) % vector(3) = 13d0
    TB(14) % vector(3) = 2d0 ! CH4
    TBid(14) % vector(4) = 14d0
    TB(14) % vector(4) = 1.5d0 ! CO
    TBid(14) % vector(5) = 15d0
    TB(14) % vector(5) = 2d0 ! CO2
    TBid(14) % vector(6) = 26d0
    TB(14) % vector(6) = 3d0 ! C2H6
    TBid(14) % vector(7) = 48d0
    TB(14) % vector(7) = 0.69999999999999996d0 ! AR

    ! (85):  2 OH <=> O + H2O
    fwd_A(107)     = 35700d0
    fwd_beta(107)  = 2.3999999999999999d0
    fwd_Ea(107)    = -2110d0
    prefactor_units(107)  = 1.0000000000000002d-06
    activation_units(107) = 0.50321666580471969d0
    phase_units(107)      = 1d-12
    is_PD(107) = 0
    nTB(107) = 0

    ! (86):  OH + HO2 <=> O2 + H2O
    fwd_A(108)     = 14500000000000d0
    fwd_beta(108)  = 0d0
    fwd_Ea(108)    = -500d0
    prefactor_units(108)  = 1.0000000000000002d-06
    activation_units(108) = 0.50321666580471969d0
    phase_units(108)      = 1d-12
    is_PD(108) = 0
    nTB(108) = 0

    ! (87):  OH + H2O2 <=> HO2 + H2O
    fwd_A(109)     = 2000000000000d0
    fwd_beta(109)  = 0d0
    fwd_Ea(109)    = 427d0
    prefactor_units(109)  = 1.0000000000000002d-06
    activation_units(109) = 0.50321666580471969d0
    phase_units(109)      = 1d-12
    is_PD(109) = 0
    nTB(109) = 0

    ! (88):  OH + H2O2 <=> HO2 + H2O
    fwd_A(110)     = 1.7d+18
    fwd_beta(110)  = 0d0
    fwd_Ea(110)    = 29410d0
    prefactor_units(110)  = 1.0000000000000002d-06
    activation_units(110) = 0.50321666580471969d0
    phase_units(110)      = 1d-12
    is_PD(110) = 0
    nTB(110) = 0

    ! (89):  OH + C <=> H + CO
    fwd_A(111)     = 50000000000000d0
    fwd_beta(111)  = 0d0
    fwd_Ea(111)    = 0d0
    prefactor_units(111)  = 1.0000000000000002d-06
    activation_units(111) = 0.50321666580471969d0
    phase_units(111)      = 1d-12
    is_PD(111) = 0
    nTB(111) = 0

    ! (90):  OH + CH <=> H + HCO
    fwd_A(112)     = 30000000000000d0
    fwd_beta(112)  = 0d0
    fwd_Ea(112)    = 0d0
    prefactor_units(112)  = 1.0000000000000002d-06
    activation_units(112) = 0.50321666580471969d0
    phase_units(112)      = 1d-12
    is_PD(112) = 0
    nTB(112) = 0

    ! (91):  OH + CH2 <=> H + CH2O
    fwd_A(113)     = 20000000000000d0
    fwd_beta(113)  = 0d0
    fwd_Ea(113)    = 0d0
    prefactor_units(113)  = 1.0000000000000002d-06
    activation_units(113) = 0.50321666580471969d0
    phase_units(113)      = 1d-12
    is_PD(113) = 0
    nTB(113) = 0

    ! (92):  OH + CH2 <=> CH + H2O
    fwd_A(114)     = 11300000d0
    fwd_beta(114)  = 2d0
    fwd_Ea(114)    = 3000d0
    prefactor_units(114)  = 1.0000000000000002d-06
    activation_units(114) = 0.50321666580471969d0
    phase_units(114)      = 1d-12
    is_PD(114) = 0
    nTB(114) = 0

    ! (93):  OH + CH2(S) <=> H + CH2O
    fwd_A(115)     = 30000000000000d0
    fwd_beta(115)  = 0d0
    fwd_Ea(115)    = 0d0
    prefactor_units(115)  = 1.0000000000000002d-06
    activation_units(115) = 0.50321666580471969d0
    phase_units(115)      = 1d-12
    is_PD(115) = 0
    nTB(115) = 0

    ! (94):  OH + CH3 (+M) <=> CH3OH (+M)
    fwd_A(15)     = 2.79d+18
    fwd_beta(15)  = -1.4299999999999999d0
    fwd_Ea(15)    = 1330d0
    low_A(15)     = 4.0000000000000002d+36
    low_beta(15)  = -5.9199999999999999d0
    low_Ea(15)    = 3140d0
    troe_a(15)    = 0.41199999999999998d0
    troe_Tsss(15) = 195d0
    troe_Ts(15)   = 5900d0
    troe_Tss(15)  = 6394d0
    troe_len(15)  = 4
    prefactor_units(15)  = 1.0000000000000002d-06
    activation_units(15) = 0.50321666580471969d0
    phase_units(15)      = 1d-12
    is_PD(15) = 1
    nTB(15) = 6
    if (.not. allocated(TB(15) % vector)) allocate(TB(15) % vector(6))
    if (.not. allocated(TBid(15) % vector)) allocate(TBid(15) % vector(6))
    TBid(15) % vector(1) = 0d0
    TB(15) % vector(1) = 2d0 ! H2
    TBid(15) % vector(2) = 5d0
    TB(15) % vector(2) = 6d0 ! H2O
    TBid(15) % vector(3) = 13d0
    TB(15) % vector(3) = 2d0 ! CH4
    TBid(15) % vector(4) = 14d0
    TB(15) % vector(4) = 1.5d0 ! CO
    TBid(15) % vector(5) = 15d0
    TB(15) % vector(5) = 2d0 ! CO2
    TBid(15) % vector(6) = 26d0
    TB(15) % vector(6) = 3d0 ! C2H6

    ! (95):  OH + CH3 <=> CH2 + H2O
    fwd_A(116)     = 56000000d0
    fwd_beta(116)  = 1.6000000000000001d0
    fwd_Ea(116)    = 5420d0
    prefactor_units(116)  = 1.0000000000000002d-06
    activation_units(116) = 0.50321666580471969d0
    phase_units(116)      = 1d-12
    is_PD(116) = 0
    nTB(116) = 0

    ! (96):  OH + CH3 <=> CH2(S) + H2O
    fwd_A(117)     = 6.44d+17
    fwd_beta(117)  = -1.3400000000000001d0
    fwd_Ea(117)    = 1417d0
    prefactor_units(117)  = 1.0000000000000002d-06
    activation_units(117) = 0.50321666580471969d0
    phase_units(117)      = 1d-12
    is_PD(117) = 0
    nTB(117) = 0

    ! (97):  OH + CH4 <=> CH3 + H2O
    fwd_A(118)     = 100000000d0
    fwd_beta(118)  = 1.6000000000000001d0
    fwd_Ea(118)    = 3120d0
    prefactor_units(118)  = 1.0000000000000002d-06
    activation_units(118) = 0.50321666580471969d0
    phase_units(118)      = 1d-12
    is_PD(118) = 0
    nTB(118) = 0

    ! (98):  OH + CO <=> H + CO2
    fwd_A(119)     = 47600000d0
    fwd_beta(119)  = 1.228d0
    fwd_Ea(119)    = 70d0
    prefactor_units(119)  = 1.0000000000000002d-06
    activation_units(119) = 0.50321666580471969d0
    phase_units(119)      = 1d-12
    is_PD(119) = 0
    nTB(119) = 0

    ! (99):  OH + HCO <=> H2O + CO
    fwd_A(120)     = 50000000000000d0
    fwd_beta(120)  = 0d0
    fwd_Ea(120)    = 0d0
    prefactor_units(120)  = 1.0000000000000002d-06
    activation_units(120) = 0.50321666580471969d0
    phase_units(120)      = 1d-12
    is_PD(120) = 0
    nTB(120) = 0

    ! (100):  OH + CH2O <=> HCO + H2O
    fwd_A(121)     = 3430000000d0
    fwd_beta(121)  = 1.1799999999999999d0
    fwd_Ea(121)    = -447d0
    prefactor_units(121)  = 1.0000000000000002d-06
    activation_units(121) = 0.50321666580471969d0
    phase_units(121)      = 1d-12
    is_PD(121) = 0
    nTB(121) = 0

    ! (101):  OH + CH2OH <=> H2O + CH2O
    fwd_A(122)     = 5000000000000d0
    fwd_beta(122)  = 0d0
    fwd_Ea(122)    = 0d0
    prefactor_units(122)  = 1.0000000000000002d-06
    activation_units(122) = 0.50321666580471969d0
    phase_units(122)      = 1d-12
    is_PD(122) = 0
    nTB(122) = 0

    ! (102):  OH + CH3O <=> H2O + CH2O
    fwd_A(123)     = 5000000000000d0
    fwd_beta(123)  = 0d0
    fwd_Ea(123)    = 0d0
    prefactor_units(123)  = 1.0000000000000002d-06
    activation_units(123) = 0.50321666580471969d0
    phase_units(123)      = 1d-12
    is_PD(123) = 0
    nTB(123) = 0

    ! (103):  OH + CH3OH <=> CH2OH + H2O
    fwd_A(124)     = 1440000d0
    fwd_beta(124)  = 2d0
    fwd_Ea(124)    = -840d0
    prefactor_units(124)  = 1.0000000000000002d-06
    activation_units(124) = 0.50321666580471969d0
    phase_units(124)      = 1d-12
    is_PD(124) = 0
    nTB(124) = 0

    ! (104):  OH + CH3OH <=> CH3O + H2O
    fwd_A(125)     = 6300000d0
    fwd_beta(125)  = 2d0
    fwd_Ea(125)    = 1500d0
    prefactor_units(125)  = 1.0000000000000002d-06
    activation_units(125) = 0.50321666580471969d0
    phase_units(125)      = 1d-12
    is_PD(125) = 0
    nTB(125) = 0

    ! (105):  OH + C2H <=> H + HCCO
    fwd_A(126)     = 20000000000000d0
    fwd_beta(126)  = 0d0
    fwd_Ea(126)    = 0d0
    prefactor_units(126)  = 1.0000000000000002d-06
    activation_units(126) = 0.50321666580471969d0
    phase_units(126)      = 1d-12
    is_PD(126) = 0
    nTB(126) = 0

    ! (106):  OH + C2H2 <=> H + CH2CO
    fwd_A(127)     = 0.00021800000000000001d0
    fwd_beta(127)  = 4.5d0
    fwd_Ea(127)    = -1000d0
    prefactor_units(127)  = 1.0000000000000002d-06
    activation_units(127) = 0.50321666580471969d0
    phase_units(127)      = 1d-12
    is_PD(127) = 0
    nTB(127) = 0

    ! (107):  OH + C2H2 <=> H + HCCOH
    fwd_A(128)     = 504000d0
    fwd_beta(128)  = 2.2999999999999998d0
    fwd_Ea(128)    = 13500d0
    prefactor_units(128)  = 1.0000000000000002d-06
    activation_units(128) = 0.50321666580471969d0
    phase_units(128)      = 1d-12
    is_PD(128) = 0
    nTB(128) = 0

    ! (108):  OH + C2H2 <=> C2H + H2O
    fwd_A(129)     = 33700000d0
    fwd_beta(129)  = 2d0
    fwd_Ea(129)    = 14000d0
    prefactor_units(129)  = 1.0000000000000002d-06
    activation_units(129) = 0.50321666580471969d0
    phase_units(129)      = 1d-12
    is_PD(129) = 0
    nTB(129) = 0

    ! (109):  OH + C2H2 <=> CH3 + CO
    fwd_A(130)     = 0.00048299999999999998d0
    fwd_beta(130)  = 4d0
    fwd_Ea(130)    = -2000d0
    prefactor_units(130)  = 1.0000000000000002d-06
    activation_units(130) = 0.50321666580471969d0
    phase_units(130)      = 1d-12
    is_PD(130) = 0
    nTB(130) = 0

    ! (110):  OH + C2H3 <=> H2O + C2H2
    fwd_A(131)     = 5000000000000d0
    fwd_beta(131)  = 0d0
    fwd_Ea(131)    = 0d0
    prefactor_units(131)  = 1.0000000000000002d-06
    activation_units(131) = 0.50321666580471969d0
    phase_units(131)      = 1d-12
    is_PD(131) = 0
    nTB(131) = 0

    ! (111):  OH + C2H4 <=> C2H3 + H2O
    fwd_A(132)     = 3600000d0
    fwd_beta(132)  = 2d0
    fwd_Ea(132)    = 2500d0
    prefactor_units(132)  = 1.0000000000000002d-06
    activation_units(132) = 0.50321666580471969d0
    phase_units(132)      = 1d-12
    is_PD(132) = 0
    nTB(132) = 0

    ! (112):  OH + C2H6 <=> C2H5 + H2O
    fwd_A(133)     = 3540000d0
    fwd_beta(133)  = 2.1200000000000001d0
    fwd_Ea(133)    = 870d0
    prefactor_units(133)  = 1.0000000000000002d-06
    activation_units(133) = 0.50321666580471969d0
    phase_units(133)      = 1d-12
    is_PD(133) = 0
    nTB(133) = 0

    ! (113):  OH + CH2CO <=> HCCO + H2O
    fwd_A(134)     = 7500000000000d0
    fwd_beta(134)  = 0d0
    fwd_Ea(134)    = 2000d0
    prefactor_units(134)  = 1.0000000000000002d-06
    activation_units(134) = 0.50321666580471969d0
    phase_units(134)      = 1d-12
    is_PD(134) = 0
    nTB(134) = 0

    ! (114):  2 HO2 <=> O2 + H2O2
    fwd_A(135)     = 130000000000d0
    fwd_beta(135)  = 0d0
    fwd_Ea(135)    = -1630d0
    prefactor_units(135)  = 1.0000000000000002d-06
    activation_units(135) = 0.50321666580471969d0
    phase_units(135)      = 1d-12
    is_PD(135) = 0
    nTB(135) = 0

    ! (115):  2 HO2 <=> O2 + H2O2
    fwd_A(136)     = 420000000000000d0
    fwd_beta(136)  = 0d0
    fwd_Ea(136)    = 12000d0
    prefactor_units(136)  = 1.0000000000000002d-06
    activation_units(136) = 0.50321666580471969d0
    phase_units(136)      = 1d-12
    is_PD(136) = 0
    nTB(136) = 0

    ! (116):  HO2 + CH2 <=> OH + CH2O
    fwd_A(137)     = 20000000000000d0
    fwd_beta(137)  = 0d0
    fwd_Ea(137)    = 0d0
    prefactor_units(137)  = 1.0000000000000002d-06
    activation_units(137) = 0.50321666580471969d0
    phase_units(137)      = 1d-12
    is_PD(137) = 0
    nTB(137) = 0

    ! (117):  HO2 + CH3 <=> O2 + CH4
    fwd_A(138)     = 1000000000000d0
    fwd_beta(138)  = 0d0
    fwd_Ea(138)    = 0d0
    prefactor_units(138)  = 1.0000000000000002d-06
    activation_units(138) = 0.50321666580471969d0
    phase_units(138)      = 1d-12
    is_PD(138) = 0
    nTB(138) = 0

    ! (118):  HO2 + CH3 <=> OH + CH3O
    fwd_A(139)     = 37800000000000d0
    fwd_beta(139)  = 0d0
    fwd_Ea(139)    = 0d0
    prefactor_units(139)  = 1.0000000000000002d-06
    activation_units(139) = 0.50321666580471969d0
    phase_units(139)      = 1d-12
    is_PD(139) = 0
    nTB(139) = 0

    ! (119):  HO2 + CO <=> OH + CO2
    fwd_A(140)     = 150000000000000d0
    fwd_beta(140)  = 0d0
    fwd_Ea(140)    = 23600d0
    prefactor_units(140)  = 1.0000000000000002d-06
    activation_units(140) = 0.50321666580471969d0
    phase_units(140)      = 1d-12
    is_PD(140) = 0
    nTB(140) = 0

    ! (120):  HO2 + CH2O <=> HCO + H2O2
    fwd_A(141)     = 5600000d0
    fwd_beta(141)  = 2d0
    fwd_Ea(141)    = 12000d0
    prefactor_units(141)  = 1.0000000000000002d-06
    activation_units(141) = 0.50321666580471969d0
    phase_units(141)      = 1d-12
    is_PD(141) = 0
    nTB(141) = 0

    ! (121):  C + O2 <=> O + CO
    fwd_A(142)     = 58000000000000d0
    fwd_beta(142)  = 0d0
    fwd_Ea(142)    = 576d0
    prefactor_units(142)  = 1.0000000000000002d-06
    activation_units(142) = 0.50321666580471969d0
    phase_units(142)      = 1d-12
    is_PD(142) = 0
    nTB(142) = 0

    ! (122):  C + CH2 <=> H + C2H
    fwd_A(143)     = 50000000000000d0
    fwd_beta(143)  = 0d0
    fwd_Ea(143)    = 0d0
    prefactor_units(143)  = 1.0000000000000002d-06
    activation_units(143) = 0.50321666580471969d0
    phase_units(143)      = 1d-12
    is_PD(143) = 0
    nTB(143) = 0

    ! (123):  C + CH3 <=> H + C2H2
    fwd_A(144)     = 50000000000000d0
    fwd_beta(144)  = 0d0
    fwd_Ea(144)    = 0d0
    prefactor_units(144)  = 1.0000000000000002d-06
    activation_units(144) = 0.50321666580471969d0
    phase_units(144)      = 1d-12
    is_PD(144) = 0
    nTB(144) = 0

    ! (124):  CH + O2 <=> O + HCO
    fwd_A(145)     = 67100000000000d0
    fwd_beta(145)  = 0d0
    fwd_Ea(145)    = 0d0
    prefactor_units(145)  = 1.0000000000000002d-06
    activation_units(145) = 0.50321666580471969d0
    phase_units(145)      = 1d-12
    is_PD(145) = 0
    nTB(145) = 0

    ! (125):  CH + H2 <=> H + CH2
    fwd_A(146)     = 108000000000000d0
    fwd_beta(146)  = 0d0
    fwd_Ea(146)    = 3110d0
    prefactor_units(146)  = 1.0000000000000002d-06
    activation_units(146) = 0.50321666580471969d0
    phase_units(146)      = 1d-12
    is_PD(146) = 0
    nTB(146) = 0

    ! (126):  CH + H2O <=> H + CH2O
    fwd_A(147)     = 5710000000000d0
    fwd_beta(147)  = 0d0
    fwd_Ea(147)    = -755d0
    prefactor_units(147)  = 1.0000000000000002d-06
    activation_units(147) = 0.50321666580471969d0
    phase_units(147)      = 1d-12
    is_PD(147) = 0
    nTB(147) = 0

    ! (127):  CH + CH2 <=> H + C2H2
    fwd_A(148)     = 40000000000000d0
    fwd_beta(148)  = 0d0
    fwd_Ea(148)    = 0d0
    prefactor_units(148)  = 1.0000000000000002d-06
    activation_units(148) = 0.50321666580471969d0
    phase_units(148)      = 1d-12
    is_PD(148) = 0
    nTB(148) = 0

    ! (128):  CH + CH3 <=> H + C2H3
    fwd_A(149)     = 30000000000000d0
    fwd_beta(149)  = 0d0
    fwd_Ea(149)    = 0d0
    prefactor_units(149)  = 1.0000000000000002d-06
    activation_units(149) = 0.50321666580471969d0
    phase_units(149)      = 1d-12
    is_PD(149) = 0
    nTB(149) = 0

    ! (129):  CH + CH4 <=> H + C2H4
    fwd_A(150)     = 60000000000000d0
    fwd_beta(150)  = 0d0
    fwd_Ea(150)    = 0d0
    prefactor_units(150)  = 1.0000000000000002d-06
    activation_units(150) = 0.50321666580471969d0
    phase_units(150)      = 1d-12
    is_PD(150) = 0
    nTB(150) = 0

    ! (130):  CH + CO (+M) <=> HCCO (+M)
    fwd_A(16)     = 50000000000000d0
    fwd_beta(16)  = 0d0
    fwd_Ea(16)    = 0d0
    low_A(16)     = 2.6899999999999998d+28
    low_beta(16)  = -3.7400000000000002d0
    low_Ea(16)    = 1936d0
    troe_a(16)    = 0.57569999999999999d0
    troe_Tsss(16) = 237d0
    troe_Ts(16)   = 1652d0
    troe_Tss(16)  = 5069d0
    troe_len(16)  = 4
    prefactor_units(16)  = 1.0000000000000002d-06
    activation_units(16) = 0.50321666580471969d0
    phase_units(16)      = 1d-12
    is_PD(16) = 1
    nTB(16) = 7
    if (.not. allocated(TB(16) % vector)) allocate(TB(16) % vector(7))
    if (.not. allocated(TBid(16) % vector)) allocate(TBid(16) % vector(7))
    TBid(16) % vector(1) = 0d0
    TB(16) % vector(1) = 2d0 ! H2
    TBid(16) % vector(2) = 5d0
    TB(16) % vector(2) = 6d0 ! H2O
    TBid(16) % vector(3) = 13d0
    TB(16) % vector(3) = 2d0 ! CH4
    TBid(16) % vector(4) = 14d0
    TB(16) % vector(4) = 1.5d0 ! CO
    TBid(16) % vector(5) = 15d0
    TB(16) % vector(5) = 2d0 ! CO2
    TBid(16) % vector(6) = 26d0
    TB(16) % vector(6) = 3d0 ! C2H6
    TBid(16) % vector(7) = 48d0
    TB(16) % vector(7) = 0.69999999999999996d0 ! AR

    ! (131):  CH + CO2 <=> HCO + CO
    fwd_A(151)     = 190000000000000d0
    fwd_beta(151)  = 0d0
    fwd_Ea(151)    = 15792d0
    prefactor_units(151)  = 1.0000000000000002d-06
    activation_units(151) = 0.50321666580471969d0
    phase_units(151)      = 1d-12
    is_PD(151) = 0
    nTB(151) = 0

    ! (132):  CH + CH2O <=> H + CH2CO
    fwd_A(152)     = 94600000000000d0
    fwd_beta(152)  = 0d0
    fwd_Ea(152)    = -515d0
    prefactor_units(152)  = 1.0000000000000002d-06
    activation_units(152) = 0.50321666580471969d0
    phase_units(152)      = 1d-12
    is_PD(152) = 0
    nTB(152) = 0

    ! (133):  CH + HCCO <=> CO + C2H2
    fwd_A(153)     = 50000000000000d0
    fwd_beta(153)  = 0d0
    fwd_Ea(153)    = 0d0
    prefactor_units(153)  = 1.0000000000000002d-06
    activation_units(153) = 0.50321666580471969d0
    phase_units(153)      = 1d-12
    is_PD(153) = 0
    nTB(153) = 0

    ! (134):  CH2 + O2 => OH + H + CO
    fwd_A(154)     = 5000000000000d0
    fwd_beta(154)  = 0d0
    fwd_Ea(154)    = 1500d0
    prefactor_units(154)  = 1.0000000000000002d-06
    activation_units(154) = 0.50321666580471969d0
    phase_units(154)      = 1d-12
    is_PD(154) = 0
    nTB(154) = 0

    ! (135):  CH2 + H2 <=> H + CH3
    fwd_A(155)     = 500000d0
    fwd_beta(155)  = 2d0
    fwd_Ea(155)    = 7230d0
    prefactor_units(155)  = 1.0000000000000002d-06
    activation_units(155) = 0.50321666580471969d0
    phase_units(155)      = 1d-12
    is_PD(155) = 0
    nTB(155) = 0

    ! (136):  2 CH2 <=> H2 + C2H2
    fwd_A(156)     = 1600000000000000d0
    fwd_beta(156)  = 0d0
    fwd_Ea(156)    = 11944d0
    prefactor_units(156)  = 1.0000000000000002d-06
    activation_units(156) = 0.50321666580471969d0
    phase_units(156)      = 1d-12
    is_PD(156) = 0
    nTB(156) = 0

    ! (137):  CH2 + CH3 <=> H + C2H4
    fwd_A(157)     = 40000000000000d0
    fwd_beta(157)  = 0d0
    fwd_Ea(157)    = 0d0
    prefactor_units(157)  = 1.0000000000000002d-06
    activation_units(157) = 0.50321666580471969d0
    phase_units(157)      = 1d-12
    is_PD(157) = 0
    nTB(157) = 0

    ! (138):  CH2 + CH4 <=> 2 CH3
    fwd_A(158)     = 2460000d0
    fwd_beta(158)  = 2d0
    fwd_Ea(158)    = 8270d0
    prefactor_units(158)  = 1.0000000000000002d-06
    activation_units(158) = 0.50321666580471969d0
    phase_units(158)      = 1d-12
    is_PD(158) = 0
    nTB(158) = 0

    ! (139):  CH2 + CO (+M) <=> CH2CO (+M)
    fwd_A(17)     = 810000000000d0
    fwd_beta(17)  = 0.5d0
    fwd_Ea(17)    = 4510d0
    low_A(17)     = 2.69d+33
    low_beta(17)  = -5.1100000000000003d0
    low_Ea(17)    = 7095d0
    troe_a(17)    = 0.5907d0
    troe_Tsss(17) = 275d0
    troe_Ts(17)   = 1226d0
    troe_Tss(17)  = 5185d0
    troe_len(17)  = 4
    prefactor_units(17)  = 1.0000000000000002d-06
    activation_units(17) = 0.50321666580471969d0
    phase_units(17)      = 1d-12
    is_PD(17) = 1
    nTB(17) = 7
    if (.not. allocated(TB(17) % vector)) allocate(TB(17) % vector(7))
    if (.not. allocated(TBid(17) % vector)) allocate(TBid(17) % vector(7))
    TBid(17) % vector(1) = 0d0
    TB(17) % vector(1) = 2d0 ! H2
    TBid(17) % vector(2) = 5d0
    TB(17) % vector(2) = 6d0 ! H2O
    TBid(17) % vector(3) = 13d0
    TB(17) % vector(3) = 2d0 ! CH4
    TBid(17) % vector(4) = 14d0
    TB(17) % vector(4) = 1.5d0 ! CO
    TBid(17) % vector(5) = 15d0
    TB(17) % vector(5) = 2d0 ! CO2
    TBid(17) % vector(6) = 26d0
    TB(17) % vector(6) = 3d0 ! C2H6
    TBid(17) % vector(7) = 48d0
    TB(17) % vector(7) = 0.69999999999999996d0 ! AR

    ! (140):  CH2 + HCCO <=> C2H3 + CO
    fwd_A(159)     = 30000000000000d0
    fwd_beta(159)  = 0d0
    fwd_Ea(159)    = 0d0
    prefactor_units(159)  = 1.0000000000000002d-06
    activation_units(159) = 0.50321666580471969d0
    phase_units(159)      = 1d-12
    is_PD(159) = 0
    nTB(159) = 0

    ! (141):  CH2(S) + N2 <=> CH2 + N2
    fwd_A(160)     = 15000000000000d0
    fwd_beta(160)  = 0d0
    fwd_Ea(160)    = 600d0
    prefactor_units(160)  = 1.0000000000000002d-06
    activation_units(160) = 0.50321666580471969d0
    phase_units(160)      = 1d-12
    is_PD(160) = 0
    nTB(160) = 0

    ! (142):  CH2(S) + AR <=> CH2 + AR
    fwd_A(161)     = 9000000000000d0
    fwd_beta(161)  = 0d0
    fwd_Ea(161)    = 600d0
    prefactor_units(161)  = 1.0000000000000002d-06
    activation_units(161) = 0.50321666580471969d0
    phase_units(161)      = 1d-12
    is_PD(161) = 0
    nTB(161) = 0

    ! (143):  CH2(S) + O2 <=> H + OH + CO
    fwd_A(162)     = 28000000000000d0
    fwd_beta(162)  = 0d0
    fwd_Ea(162)    = 0d0
    prefactor_units(162)  = 1.0000000000000002d-06
    activation_units(162) = 0.50321666580471969d0
    phase_units(162)      = 1d-12
    is_PD(162) = 0
    nTB(162) = 0

    ! (144):  CH2(S) + O2 <=> CO + H2O
    fwd_A(163)     = 12000000000000d0
    fwd_beta(163)  = 0d0
    fwd_Ea(163)    = 0d0
    prefactor_units(163)  = 1.0000000000000002d-06
    activation_units(163) = 0.50321666580471969d0
    phase_units(163)      = 1d-12
    is_PD(163) = 0
    nTB(163) = 0

    ! (145):  CH2(S) + H2 <=> CH3 + H
    fwd_A(164)     = 70000000000000d0
    fwd_beta(164)  = 0d0
    fwd_Ea(164)    = 0d0
    prefactor_units(164)  = 1.0000000000000002d-06
    activation_units(164) = 0.50321666580471969d0
    phase_units(164)      = 1d-12
    is_PD(164) = 0
    nTB(164) = 0

    ! (146):  CH2(S) + H2O (+M) <=> CH3OH (+M)
    fwd_A(18)     = 4.82d+17
    fwd_beta(18)  = -1.1599999999999999d0
    fwd_Ea(18)    = 1145d0
    low_A(18)     = 1.8799999999999998d+38
    low_beta(18)  = -6.3600000000000003d0
    low_Ea(18)    = 5040d0
    troe_a(18)    = 0.60270000000000001d0
    troe_Tsss(18) = 208d0
    troe_Ts(18)   = 3922d0
    troe_Tss(18)  = 10180d0
    troe_len(18)  = 4
    prefactor_units(18)  = 1.0000000000000002d-06
    activation_units(18) = 0.50321666580471969d0
    phase_units(18)      = 1d-12
    is_PD(18) = 1
    nTB(18) = 6
    if (.not. allocated(TB(18) % vector)) allocate(TB(18) % vector(6))
    if (.not. allocated(TBid(18) % vector)) allocate(TBid(18) % vector(6))
    TBid(18) % vector(1) = 0d0
    TB(18) % vector(1) = 2d0 ! H2
    TBid(18) % vector(2) = 5d0
    TB(18) % vector(2) = 6d0 ! H2O
    TBid(18) % vector(3) = 13d0
    TB(18) % vector(3) = 2d0 ! CH4
    TBid(18) % vector(4) = 14d0
    TB(18) % vector(4) = 1.5d0 ! CO
    TBid(18) % vector(5) = 15d0
    TB(18) % vector(5) = 2d0 ! CO2
    TBid(18) % vector(6) = 26d0
    TB(18) % vector(6) = 3d0 ! C2H6

    ! (147):  CH2(S) + H2O <=> CH2 + H2O
    fwd_A(165)     = 30000000000000d0
    fwd_beta(165)  = 0d0
    fwd_Ea(165)    = 0d0
    prefactor_units(165)  = 1.0000000000000002d-06
    activation_units(165) = 0.50321666580471969d0
    phase_units(165)      = 1d-12
    is_PD(165) = 0
    nTB(165) = 0

    ! (148):  CH2(S) + CH3 <=> H + C2H4
    fwd_A(166)     = 12000000000000d0
    fwd_beta(166)  = 0d0
    fwd_Ea(166)    = -570d0
    prefactor_units(166)  = 1.0000000000000002d-06
    activation_units(166) = 0.50321666580471969d0
    phase_units(166)      = 1d-12
    is_PD(166) = 0
    nTB(166) = 0

    ! (149):  CH2(S) + CH4 <=> 2 CH3
    fwd_A(167)     = 16000000000000d0
    fwd_beta(167)  = 0d0
    fwd_Ea(167)    = -570d0
    prefactor_units(167)  = 1.0000000000000002d-06
    activation_units(167) = 0.50321666580471969d0
    phase_units(167)      = 1d-12
    is_PD(167) = 0
    nTB(167) = 0

    ! (150):  CH2(S) + CO <=> CH2 + CO
    fwd_A(168)     = 9000000000000d0
    fwd_beta(168)  = 0d0
    fwd_Ea(168)    = 0d0
    prefactor_units(168)  = 1.0000000000000002d-06
    activation_units(168) = 0.50321666580471969d0
    phase_units(168)      = 1d-12
    is_PD(168) = 0
    nTB(168) = 0

    ! (151):  CH2(S) + CO2 <=> CH2 + CO2
    fwd_A(169)     = 7000000000000d0
    fwd_beta(169)  = 0d0
    fwd_Ea(169)    = 0d0
    prefactor_units(169)  = 1.0000000000000002d-06
    activation_units(169) = 0.50321666580471969d0
    phase_units(169)      = 1d-12
    is_PD(169) = 0
    nTB(169) = 0

    ! (152):  CH2(S) + CO2 <=> CO + CH2O
    fwd_A(170)     = 14000000000000d0
    fwd_beta(170)  = 0d0
    fwd_Ea(170)    = 0d0
    prefactor_units(170)  = 1.0000000000000002d-06
    activation_units(170) = 0.50321666580471969d0
    phase_units(170)      = 1d-12
    is_PD(170) = 0
    nTB(170) = 0

    ! (153):  CH2(S) + C2H6 <=> CH3 + C2H5
    fwd_A(171)     = 40000000000000d0
    fwd_beta(171)  = 0d0
    fwd_Ea(171)    = -550d0
    prefactor_units(171)  = 1.0000000000000002d-06
    activation_units(171) = 0.50321666580471969d0
    phase_units(171)      = 1d-12
    is_PD(171) = 0
    nTB(171) = 0

    ! (154):  CH3 + O2 <=> O + CH3O
    fwd_A(172)     = 35600000000000d0
    fwd_beta(172)  = 0d0
    fwd_Ea(172)    = 30480d0
    prefactor_units(172)  = 1.0000000000000002d-06
    activation_units(172) = 0.50321666580471969d0
    phase_units(172)      = 1d-12
    is_PD(172) = 0
    nTB(172) = 0

    ! (155):  CH3 + O2 <=> OH + CH2O
    fwd_A(173)     = 2310000000000d0
    fwd_beta(173)  = 0d0
    fwd_Ea(173)    = 20315d0
    prefactor_units(173)  = 1.0000000000000002d-06
    activation_units(173) = 0.50321666580471969d0
    phase_units(173)      = 1d-12
    is_PD(173) = 0
    nTB(173) = 0

    ! (156):  CH3 + H2O2 <=> HO2 + CH4
    fwd_A(174)     = 24500d0
    fwd_beta(174)  = 2.4700000000000002d0
    fwd_Ea(174)    = 5180d0
    prefactor_units(174)  = 1.0000000000000002d-06
    activation_units(174) = 0.50321666580471969d0
    phase_units(174)      = 1d-12
    is_PD(174) = 0
    nTB(174) = 0

    ! (157):  2 CH3 (+M) <=> C2H6 (+M)
    fwd_A(19)     = 67700000000000000d0
    fwd_beta(19)  = -1.1799999999999999d0
    fwd_Ea(19)    = 654d0
    low_A(19)     = 3.4d+41
    low_beta(19)  = -7.0300000000000002d0
    low_Ea(19)    = 2762d0
    troe_a(19)    = 0.61899999999999999d0
    troe_Tsss(19) = 73.200000000000003d0
    troe_Ts(19)   = 1180d0
    troe_Tss(19)  = 9999d0
    troe_len(19)  = 4
    prefactor_units(19)  = 1.0000000000000002d-06
    activation_units(19) = 0.50321666580471969d0
    phase_units(19)      = 1d-12
    is_PD(19) = 1
    nTB(19) = 7
    if (.not. allocated(TB(19) % vector)) allocate(TB(19) % vector(7))
    if (.not. allocated(TBid(19) % vector)) allocate(TBid(19) % vector(7))
    TBid(19) % vector(1) = 0d0
    TB(19) % vector(1) = 2d0 ! H2
    TBid(19) % vector(2) = 5d0
    TB(19) % vector(2) = 6d0 ! H2O
    TBid(19) % vector(3) = 13d0
    TB(19) % vector(3) = 2d0 ! CH4
    TBid(19) % vector(4) = 14d0
    TB(19) % vector(4) = 1.5d0 ! CO
    TBid(19) % vector(5) = 15d0
    TB(19) % vector(5) = 2d0 ! CO2
    TBid(19) % vector(6) = 26d0
    TB(19) % vector(6) = 3d0 ! C2H6
    TBid(19) % vector(7) = 48d0
    TB(19) % vector(7) = 0.69999999999999996d0 ! AR

    ! (158):  2 CH3 <=> H + C2H5
    fwd_A(175)     = 6840000000000d0
    fwd_beta(175)  = 0.10000000000000001d0
    fwd_Ea(175)    = 10600d0
    prefactor_units(175)  = 1.0000000000000002d-06
    activation_units(175) = 0.50321666580471969d0
    phase_units(175)      = 1d-12
    is_PD(175) = 0
    nTB(175) = 0

    ! (159):  CH3 + HCO <=> CH4 + CO
    fwd_A(176)     = 26480000000000d0
    fwd_beta(176)  = 0d0
    fwd_Ea(176)    = 0d0
    prefactor_units(176)  = 1.0000000000000002d-06
    activation_units(176) = 0.50321666580471969d0
    phase_units(176)      = 1d-12
    is_PD(176) = 0
    nTB(176) = 0

    ! (160):  CH3 + CH2O <=> HCO + CH4
    fwd_A(177)     = 3320d0
    fwd_beta(177)  = 2.8100000000000001d0
    fwd_Ea(177)    = 5860d0
    prefactor_units(177)  = 1.0000000000000002d-06
    activation_units(177) = 0.50321666580471969d0
    phase_units(177)      = 1d-12
    is_PD(177) = 0
    nTB(177) = 0

    ! (161):  CH3 + CH3OH <=> CH2OH + CH4
    fwd_A(178)     = 30000000d0
    fwd_beta(178)  = 1.5d0
    fwd_Ea(178)    = 9940d0
    prefactor_units(178)  = 1.0000000000000002d-06
    activation_units(178) = 0.50321666580471969d0
    phase_units(178)      = 1d-12
    is_PD(178) = 0
    nTB(178) = 0

    ! (162):  CH3 + CH3OH <=> CH3O + CH4
    fwd_A(179)     = 10000000d0
    fwd_beta(179)  = 1.5d0
    fwd_Ea(179)    = 9940d0
    prefactor_units(179)  = 1.0000000000000002d-06
    activation_units(179) = 0.50321666580471969d0
    phase_units(179)      = 1d-12
    is_PD(179) = 0
    nTB(179) = 0

    ! (163):  CH3 + C2H4 <=> C2H3 + CH4
    fwd_A(180)     = 227000d0
    fwd_beta(180)  = 2d0
    fwd_Ea(180)    = 9200d0
    prefactor_units(180)  = 1.0000000000000002d-06
    activation_units(180) = 0.50321666580471969d0
    phase_units(180)      = 1d-12
    is_PD(180) = 0
    nTB(180) = 0

    ! (164):  CH3 + C2H6 <=> C2H5 + CH4
    fwd_A(181)     = 6140000d0
    fwd_beta(181)  = 1.74d0
    fwd_Ea(181)    = 10450d0
    prefactor_units(181)  = 1.0000000000000002d-06
    activation_units(181) = 0.50321666580471969d0
    phase_units(181)      = 1d-12
    is_PD(181) = 0
    nTB(181) = 0

    ! (165):  HCO + H2O <=> H + CO + H2O
    fwd_A(182)     = 1.5d+18
    fwd_beta(182)  = -1d0
    fwd_Ea(182)    = 17000d0
    prefactor_units(182)  = 1.0000000000000002d-06
    activation_units(182) = 0.50321666580471969d0
    phase_units(182)      = 1d-12
    is_PD(182) = 0
    nTB(182) = 0

    ! (166):  HCO + M <=> H + CO + M
    fwd_A(35)     = 1.87d+17
    fwd_beta(35)  = -1d0
    fwd_Ea(35)    = 17000d0
    prefactor_units(35)  = 1.0000000000000002d-06
    activation_units(35) = 0.50321666580471969d0
    phase_units(35)      = 1d-6
    is_PD(35) = 0
    nTB(35) = 6
    if (.not. allocated(TB(35) % vector)) allocate(TB(35) % vector(6))
    if (.not. allocated(TBid(35) % vector)) allocate(TBid(35) % vector(6))
    TBid(35) % vector(1) = 0d0
    TB(35) % vector(1) = 2d0 ! H2
    TBid(35) % vector(2) = 5d0
    TB(35) % vector(2) = 0d0 ! H2O
    TBid(35) % vector(3) = 13d0
    TB(35) % vector(3) = 2d0 ! CH4
    TBid(35) % vector(4) = 14d0
    TB(35) % vector(4) = 1.5d0 ! CO
    TBid(35) % vector(5) = 15d0
    TB(35) % vector(5) = 2d0 ! CO2
    TBid(35) % vector(6) = 26d0
    TB(35) % vector(6) = 3d0 ! C2H6

    ! (167):  HCO + O2 <=> HO2 + CO
    fwd_A(183)     = 13450000000000d0
    fwd_beta(183)  = 0d0
    fwd_Ea(183)    = 400d0
    prefactor_units(183)  = 1.0000000000000002d-06
    activation_units(183) = 0.50321666580471969d0
    phase_units(183)      = 1d-12
    is_PD(183) = 0
    nTB(183) = 0

    ! (168):  CH2OH + O2 <=> HO2 + CH2O
    fwd_A(184)     = 18000000000000d0
    fwd_beta(184)  = 0d0
    fwd_Ea(184)    = 900d0
    prefactor_units(184)  = 1.0000000000000002d-06
    activation_units(184) = 0.50321666580471969d0
    phase_units(184)      = 1d-12
    is_PD(184) = 0
    nTB(184) = 0

    ! (169):  CH3O + O2 <=> HO2 + CH2O
    fwd_A(185)     = 4.2799999999999999d-13
    fwd_beta(185)  = 7.5999999999999996d0
    fwd_Ea(185)    = -3530d0
    prefactor_units(185)  = 1.0000000000000002d-06
    activation_units(185) = 0.50321666580471969d0
    phase_units(185)      = 1d-12
    is_PD(185) = 0
    nTB(185) = 0

    ! (170):  C2H + O2 <=> HCO + CO
    fwd_A(186)     = 10000000000000d0
    fwd_beta(186)  = 0d0
    fwd_Ea(186)    = -755d0
    prefactor_units(186)  = 1.0000000000000002d-06
    activation_units(186) = 0.50321666580471969d0
    phase_units(186)      = 1d-12
    is_PD(186) = 0
    nTB(186) = 0

    ! (171):  C2H + H2 <=> H + C2H2
    fwd_A(187)     = 56800000000d0
    fwd_beta(187)  = 0.90000000000000002d0
    fwd_Ea(187)    = 1993d0
    prefactor_units(187)  = 1.0000000000000002d-06
    activation_units(187) = 0.50321666580471969d0
    phase_units(187)      = 1d-12
    is_PD(187) = 0
    nTB(187) = 0

    ! (172):  C2H3 + O2 <=> HCO + CH2O
    fwd_A(188)     = 45800000000000000d0
    fwd_beta(188)  = -1.3899999999999999d0
    fwd_Ea(188)    = 1015d0
    prefactor_units(188)  = 1.0000000000000002d-06
    activation_units(188) = 0.50321666580471969d0
    phase_units(188)      = 1d-12
    is_PD(188) = 0
    nTB(188) = 0

    ! (173):  C2H4 (+M) <=> H2 + C2H2 (+M)
    fwd_A(20)     = 8000000000000d0
    fwd_beta(20)  = 0.44d0
    fwd_Ea(20)    = 86770d0
    low_A(20)     = 1.5800000000000002d+51
    low_beta(20)  = -9.3000000000000007d0
    low_Ea(20)    = 97800d0
    troe_a(20)    = 0.73450000000000004d0
    troe_Tsss(20) = 180d0
    troe_Ts(20)   = 1035d0
    troe_Tss(20)  = 5417d0
    troe_len(20)  = 4
    prefactor_units(20)  = 1d0
    activation_units(20) = 0.50321666580471969d0
    phase_units(20)      = 1d-6
    is_PD(20) = 1
    nTB(20) = 7
    if (.not. allocated(TB(20) % vector)) allocate(TB(20) % vector(7))
    if (.not. allocated(TBid(20) % vector)) allocate(TBid(20) % vector(7))
    TBid(20) % vector(1) = 0d0
    TB(20) % vector(1) = 2d0 ! H2
    TBid(20) % vector(2) = 5d0
    TB(20) % vector(2) = 6d0 ! H2O
    TBid(20) % vector(3) = 13d0
    TB(20) % vector(3) = 2d0 ! CH4
    TBid(20) % vector(4) = 14d0
    TB(20) % vector(4) = 1.5d0 ! CO
    TBid(20) % vector(5) = 15d0
    TB(20) % vector(5) = 2d0 ! CO2
    TBid(20) % vector(6) = 26d0
    TB(20) % vector(6) = 3d0 ! C2H6
    TBid(20) % vector(7) = 48d0
    TB(20) % vector(7) = 0.69999999999999996d0 ! AR

    ! (174):  C2H5 + O2 <=> HO2 + C2H4
    fwd_A(189)     = 840000000000d0
    fwd_beta(189)  = 0d0
    fwd_Ea(189)    = 3875d0
    prefactor_units(189)  = 1.0000000000000002d-06
    activation_units(189) = 0.50321666580471969d0
    phase_units(189)      = 1d-12
    is_PD(189) = 0
    nTB(189) = 0

    ! (175):  HCCO + O2 <=> OH + 2 CO
    fwd_A(190)     = 3200000000000d0
    fwd_beta(190)  = 0d0
    fwd_Ea(190)    = 854d0
    prefactor_units(190)  = 1.0000000000000002d-06
    activation_units(190) = 0.50321666580471969d0
    phase_units(190)      = 1d-12
    is_PD(190) = 0
    nTB(190) = 0

    ! (176):  2 HCCO <=> 2 CO + C2H2
    fwd_A(191)     = 10000000000000d0
    fwd_beta(191)  = 0d0
    fwd_Ea(191)    = 0d0
    prefactor_units(191)  = 1.0000000000000002d-06
    activation_units(191) = 0.50321666580471969d0
    phase_units(191)      = 1d-12
    is_PD(191) = 0
    nTB(191) = 0

    ! (177):  N + NO <=> N2 + O
    fwd_A(192)     = 27000000000000d0
    fwd_beta(192)  = 0d0
    fwd_Ea(192)    = 355d0
    prefactor_units(192)  = 1.0000000000000002d-06
    activation_units(192) = 0.50321666580471969d0
    phase_units(192)      = 1d-12
    is_PD(192) = 0
    nTB(192) = 0

    ! (178):  N + O2 <=> NO + O
    fwd_A(193)     = 9000000000d0
    fwd_beta(193)  = 1d0
    fwd_Ea(193)    = 6500d0
    prefactor_units(193)  = 1.0000000000000002d-06
    activation_units(193) = 0.50321666580471969d0
    phase_units(193)      = 1d-12
    is_PD(193) = 0
    nTB(193) = 0

    ! (179):  N + OH <=> NO + H
    fwd_A(194)     = 33600000000000d0
    fwd_beta(194)  = 0d0
    fwd_Ea(194)    = 385d0
    prefactor_units(194)  = 1.0000000000000002d-06
    activation_units(194) = 0.50321666580471969d0
    phase_units(194)      = 1d-12
    is_PD(194) = 0
    nTB(194) = 0

    ! (180):  N2O + O <=> N2 + O2
    fwd_A(195)     = 1400000000000d0
    fwd_beta(195)  = 0d0
    fwd_Ea(195)    = 10810d0
    prefactor_units(195)  = 1.0000000000000002d-06
    activation_units(195) = 0.50321666580471969d0
    phase_units(195)      = 1d-12
    is_PD(195) = 0
    nTB(195) = 0

    ! (181):  N2O + O <=> 2 NO
    fwd_A(196)     = 29000000000000d0
    fwd_beta(196)  = 0d0
    fwd_Ea(196)    = 23150d0
    prefactor_units(196)  = 1.0000000000000002d-06
    activation_units(196) = 0.50321666580471969d0
    phase_units(196)      = 1d-12
    is_PD(196) = 0
    nTB(196) = 0

    ! (182):  N2O + H <=> N2 + OH
    fwd_A(197)     = 387000000000000d0
    fwd_beta(197)  = 0d0
    fwd_Ea(197)    = 18880d0
    prefactor_units(197)  = 1.0000000000000002d-06
    activation_units(197) = 0.50321666580471969d0
    phase_units(197)      = 1d-12
    is_PD(197) = 0
    nTB(197) = 0

    ! (183):  N2O + OH <=> N2 + HO2
    fwd_A(198)     = 2000000000000d0
    fwd_beta(198)  = 0d0
    fwd_Ea(198)    = 21060d0
    prefactor_units(198)  = 1.0000000000000002d-06
    activation_units(198) = 0.50321666580471969d0
    phase_units(198)      = 1d-12
    is_PD(198) = 0
    nTB(198) = 0

    ! (184):  N2O (+M) <=> N2 + O (+M)
    fwd_A(28)     = 79100000000d0
    fwd_beta(28)  = 0d0
    fwd_Ea(28)    = 56020d0
    low_A(28)     = 637000000000000d0
    low_beta(28)  = 0d0
    low_Ea(28)    = 56640d0
    prefactor_units(28)  = 1d0
    activation_units(28) = 0.50321666580471969d0
    phase_units(28)      = 1d-6
    is_PD(28) = 1
    nTB(28) = 7
    if (.not. allocated(TB(28) % vector)) allocate(TB(28) % vector(7))
    if (.not. allocated(TBid(28) % vector)) allocate(TBid(28) % vector(7))
    TBid(28) % vector(1) = 0d0
    TB(28) % vector(1) = 2d0 ! H2
    TBid(28) % vector(2) = 5d0
    TB(28) % vector(2) = 6d0 ! H2O
    TBid(28) % vector(3) = 13d0
    TB(28) % vector(3) = 2d0 ! CH4
    TBid(28) % vector(4) = 14d0
    TB(28) % vector(4) = 1.5d0 ! CO
    TBid(28) % vector(5) = 15d0
    TB(28) % vector(5) = 2d0 ! CO2
    TBid(28) % vector(6) = 26d0
    TB(28) % vector(6) = 3d0 ! C2H6
    TBid(28) % vector(7) = 48d0
    TB(28) % vector(7) = 0.625d0 ! AR

    ! (185):  HO2 + NO <=> NO2 + OH
    fwd_A(199)     = 2110000000000d0
    fwd_beta(199)  = 0d0
    fwd_Ea(199)    = -480d0
    prefactor_units(199)  = 1.0000000000000002d-06
    activation_units(199) = 0.50321666580471969d0
    phase_units(199)      = 1d-12
    is_PD(199) = 0
    nTB(199) = 0

    ! (186):  NO + O + M <=> NO2 + M
    fwd_A(36)     = 1.06d+20
    fwd_beta(36)  = -1.4099999999999999d0
    fwd_Ea(36)    = 0d0
    prefactor_units(36)  = 1.0000000000000002d-12
    activation_units(36) = 0.50321666580471969d0
    phase_units(36)      = 1d-12
    is_PD(36) = 0
    nTB(36) = 7
    if (.not. allocated(TB(36) % vector)) allocate(TB(36) % vector(7))
    if (.not. allocated(TBid(36) % vector)) allocate(TBid(36) % vector(7))
    TBid(36) % vector(1) = 0d0
    TB(36) % vector(1) = 2d0 ! H2
    TBid(36) % vector(2) = 5d0
    TB(36) % vector(2) = 6d0 ! H2O
    TBid(36) % vector(3) = 13d0
    TB(36) % vector(3) = 2d0 ! CH4
    TBid(36) % vector(4) = 14d0
    TB(36) % vector(4) = 1.5d0 ! CO
    TBid(36) % vector(5) = 15d0
    TB(36) % vector(5) = 2d0 ! CO2
    TBid(36) % vector(6) = 26d0
    TB(36) % vector(6) = 3d0 ! C2H6
    TBid(36) % vector(7) = 48d0
    TB(36) % vector(7) = 0.69999999999999996d0 ! AR

    ! (187):  NO2 + O <=> NO + O2
    fwd_A(200)     = 3900000000000d0
    fwd_beta(200)  = 0d0
    fwd_Ea(200)    = -240d0
    prefactor_units(200)  = 1.0000000000000002d-06
    activation_units(200) = 0.50321666580471969d0
    phase_units(200)      = 1d-12
    is_PD(200) = 0
    nTB(200) = 0

    ! (188):  NO2 + H <=> NO + OH
    fwd_A(201)     = 132000000000000d0
    fwd_beta(201)  = 0d0
    fwd_Ea(201)    = 360d0
    prefactor_units(201)  = 1.0000000000000002d-06
    activation_units(201) = 0.50321666580471969d0
    phase_units(201)      = 1d-12
    is_PD(201) = 0
    nTB(201) = 0

    ! (189):  NH + O <=> NO + H
    fwd_A(202)     = 40000000000000d0
    fwd_beta(202)  = 0d0
    fwd_Ea(202)    = 0d0
    prefactor_units(202)  = 1.0000000000000002d-06
    activation_units(202) = 0.50321666580471969d0
    phase_units(202)      = 1d-12
    is_PD(202) = 0
    nTB(202) = 0

    ! (190):  NH + H <=> N + H2
    fwd_A(203)     = 32000000000000d0
    fwd_beta(203)  = 0d0
    fwd_Ea(203)    = 330d0
    prefactor_units(203)  = 1.0000000000000002d-06
    activation_units(203) = 0.50321666580471969d0
    phase_units(203)      = 1d-12
    is_PD(203) = 0
    nTB(203) = 0

    ! (191):  NH + OH <=> HNO + H
    fwd_A(204)     = 20000000000000d0
    fwd_beta(204)  = 0d0
    fwd_Ea(204)    = 0d0
    prefactor_units(204)  = 1.0000000000000002d-06
    activation_units(204) = 0.50321666580471969d0
    phase_units(204)      = 1d-12
    is_PD(204) = 0
    nTB(204) = 0

    ! (192):  NH + OH <=> N + H2O
    fwd_A(205)     = 2000000000d0
    fwd_beta(205)  = 1.2d0
    fwd_Ea(205)    = 0d0
    prefactor_units(205)  = 1.0000000000000002d-06
    activation_units(205) = 0.50321666580471969d0
    phase_units(205)      = 1d-12
    is_PD(205) = 0
    nTB(205) = 0

    ! (193):  NH + O2 <=> HNO + O
    fwd_A(206)     = 461000d0
    fwd_beta(206)  = 2d0
    fwd_Ea(206)    = 6500d0
    prefactor_units(206)  = 1.0000000000000002d-06
    activation_units(206) = 0.50321666580471969d0
    phase_units(206)      = 1d-12
    is_PD(206) = 0
    nTB(206) = 0

    ! (194):  NH + O2 <=> NO + OH
    fwd_A(207)     = 1280000d0
    fwd_beta(207)  = 1.5d0
    fwd_Ea(207)    = 100d0
    prefactor_units(207)  = 1.0000000000000002d-06
    activation_units(207) = 0.50321666580471969d0
    phase_units(207)      = 1d-12
    is_PD(207) = 0
    nTB(207) = 0

    ! (195):  NH + N <=> N2 + H
    fwd_A(208)     = 15000000000000d0
    fwd_beta(208)  = 0d0
    fwd_Ea(208)    = 0d0
    prefactor_units(208)  = 1.0000000000000002d-06
    activation_units(208) = 0.50321666580471969d0
    phase_units(208)      = 1d-12
    is_PD(208) = 0
    nTB(208) = 0

    ! (196):  NH + H2O <=> HNO + H2
    fwd_A(209)     = 20000000000000d0
    fwd_beta(209)  = 0d0
    fwd_Ea(209)    = 13850d0
    prefactor_units(209)  = 1.0000000000000002d-06
    activation_units(209) = 0.50321666580471969d0
    phase_units(209)      = 1d-12
    is_PD(209) = 0
    nTB(209) = 0

    ! (197):  NH + NO <=> N2 + OH
    fwd_A(210)     = 21600000000000d0
    fwd_beta(210)  = -0.23000000000000001d0
    fwd_Ea(210)    = 0d0
    prefactor_units(210)  = 1.0000000000000002d-06
    activation_units(210) = 0.50321666580471969d0
    phase_units(210)      = 1d-12
    is_PD(210) = 0
    nTB(210) = 0

    ! (198):  NH + NO <=> N2O + H
    fwd_A(211)     = 365000000000000d0
    fwd_beta(211)  = -0.45000000000000001d0
    fwd_Ea(211)    = 0d0
    prefactor_units(211)  = 1.0000000000000002d-06
    activation_units(211) = 0.50321666580471969d0
    phase_units(211)      = 1d-12
    is_PD(211) = 0
    nTB(211) = 0

    ! (199):  NH2 + O <=> OH + NH
    fwd_A(212)     = 3000000000000d0
    fwd_beta(212)  = 0d0
    fwd_Ea(212)    = 0d0
    prefactor_units(212)  = 1.0000000000000002d-06
    activation_units(212) = 0.50321666580471969d0
    phase_units(212)      = 1d-12
    is_PD(212) = 0
    nTB(212) = 0

    ! (200):  NH2 + O <=> H + HNO
    fwd_A(213)     = 39000000000000d0
    fwd_beta(213)  = 0d0
    fwd_Ea(213)    = 0d0
    prefactor_units(213)  = 1.0000000000000002d-06
    activation_units(213) = 0.50321666580471969d0
    phase_units(213)      = 1d-12
    is_PD(213) = 0
    nTB(213) = 0

    ! (201):  NH2 + H <=> NH + H2
    fwd_A(214)     = 40000000000000d0
    fwd_beta(214)  = 0d0
    fwd_Ea(214)    = 3650d0
    prefactor_units(214)  = 1.0000000000000002d-06
    activation_units(214) = 0.50321666580471969d0
    phase_units(214)      = 1d-12
    is_PD(214) = 0
    nTB(214) = 0

    ! (202):  NH2 + OH <=> NH + H2O
    fwd_A(215)     = 90000000d0
    fwd_beta(215)  = 1.5d0
    fwd_Ea(215)    = -460d0
    prefactor_units(215)  = 1.0000000000000002d-06
    activation_units(215) = 0.50321666580471969d0
    phase_units(215)      = 1d-12
    is_PD(215) = 0
    nTB(215) = 0

    ! (203):  NNH <=> N2 + H
    fwd_A(216)     = 330000000d0
    fwd_beta(216)  = 0d0
    fwd_Ea(216)    = 0d0
    prefactor_units(216)  = 1d0
    activation_units(216) = 0.50321666580471969d0
    phase_units(216)      = 1d-6
    is_PD(216) = 0
    nTB(216) = 0

    ! (204):  NNH + M <=> N2 + H + M
    fwd_A(37)     = 130000000000000d0
    fwd_beta(37)  = -0.11d0
    fwd_Ea(37)    = 4980d0
    prefactor_units(37)  = 1.0000000000000002d-06
    activation_units(37) = 0.50321666580471969d0
    phase_units(37)      = 1d-6
    is_PD(37) = 0
    nTB(37) = 7
    if (.not. allocated(TB(37) % vector)) allocate(TB(37) % vector(7))
    if (.not. allocated(TBid(37) % vector)) allocate(TBid(37) % vector(7))
    TBid(37) % vector(1) = 0d0
    TB(37) % vector(1) = 2d0 ! H2
    TBid(37) % vector(2) = 5d0
    TB(37) % vector(2) = 6d0 ! H2O
    TBid(37) % vector(3) = 13d0
    TB(37) % vector(3) = 2d0 ! CH4
    TBid(37) % vector(4) = 14d0
    TB(37) % vector(4) = 1.5d0 ! CO
    TBid(37) % vector(5) = 15d0
    TB(37) % vector(5) = 2d0 ! CO2
    TBid(37) % vector(6) = 26d0
    TB(37) % vector(6) = 3d0 ! C2H6
    TBid(37) % vector(7) = 48d0
    TB(37) % vector(7) = 0.69999999999999996d0 ! AR

    ! (205):  NNH + O2 <=> HO2 + N2
    fwd_A(217)     = 5000000000000d0
    fwd_beta(217)  = 0d0
    fwd_Ea(217)    = 0d0
    prefactor_units(217)  = 1.0000000000000002d-06
    activation_units(217) = 0.50321666580471969d0
    phase_units(217)      = 1d-12
    is_PD(217) = 0
    nTB(217) = 0

    ! (206):  NNH + O <=> OH + N2
    fwd_A(218)     = 25000000000000d0
    fwd_beta(218)  = 0d0
    fwd_Ea(218)    = 0d0
    prefactor_units(218)  = 1.0000000000000002d-06
    activation_units(218) = 0.50321666580471969d0
    phase_units(218)      = 1d-12
    is_PD(218) = 0
    nTB(218) = 0

    ! (207):  NNH + O <=> NH + NO
    fwd_A(219)     = 70000000000000d0
    fwd_beta(219)  = 0d0
    fwd_Ea(219)    = 0d0
    prefactor_units(219)  = 1.0000000000000002d-06
    activation_units(219) = 0.50321666580471969d0
    phase_units(219)      = 1d-12
    is_PD(219) = 0
    nTB(219) = 0

    ! (208):  NNH + H <=> H2 + N2
    fwd_A(220)     = 50000000000000d0
    fwd_beta(220)  = 0d0
    fwd_Ea(220)    = 0d0
    prefactor_units(220)  = 1.0000000000000002d-06
    activation_units(220) = 0.50321666580471969d0
    phase_units(220)      = 1d-12
    is_PD(220) = 0
    nTB(220) = 0

    ! (209):  NNH + OH <=> H2O + N2
    fwd_A(221)     = 20000000000000d0
    fwd_beta(221)  = 0d0
    fwd_Ea(221)    = 0d0
    prefactor_units(221)  = 1.0000000000000002d-06
    activation_units(221) = 0.50321666580471969d0
    phase_units(221)      = 1d-12
    is_PD(221) = 0
    nTB(221) = 0

    ! (210):  NNH + CH3 <=> CH4 + N2
    fwd_A(222)     = 25000000000000d0
    fwd_beta(222)  = 0d0
    fwd_Ea(222)    = 0d0
    prefactor_units(222)  = 1.0000000000000002d-06
    activation_units(222) = 0.50321666580471969d0
    phase_units(222)      = 1d-12
    is_PD(222) = 0
    nTB(222) = 0

    ! (211):  H + NO + M <=> HNO + M
    fwd_A(38)     = 4.48d+19
    fwd_beta(38)  = -1.3200000000000001d0
    fwd_Ea(38)    = 740d0
    prefactor_units(38)  = 1.0000000000000002d-12
    activation_units(38) = 0.50321666580471969d0
    phase_units(38)      = 1d-12
    is_PD(38) = 0
    nTB(38) = 7
    if (.not. allocated(TB(38) % vector)) allocate(TB(38) % vector(7))
    if (.not. allocated(TBid(38) % vector)) allocate(TBid(38) % vector(7))
    TBid(38) % vector(1) = 0d0
    TB(38) % vector(1) = 2d0 ! H2
    TBid(38) % vector(2) = 5d0
    TB(38) % vector(2) = 6d0 ! H2O
    TBid(38) % vector(3) = 13d0
    TB(38) % vector(3) = 2d0 ! CH4
    TBid(38) % vector(4) = 14d0
    TB(38) % vector(4) = 1.5d0 ! CO
    TBid(38) % vector(5) = 15d0
    TB(38) % vector(5) = 2d0 ! CO2
    TBid(38) % vector(6) = 26d0
    TB(38) % vector(6) = 3d0 ! C2H6
    TBid(38) % vector(7) = 48d0
    TB(38) % vector(7) = 0.69999999999999996d0 ! AR

    ! (212):  HNO + O <=> NO + OH
    fwd_A(223)     = 25000000000000d0
    fwd_beta(223)  = 0d0
    fwd_Ea(223)    = 0d0
    prefactor_units(223)  = 1.0000000000000002d-06
    activation_units(223) = 0.50321666580471969d0
    phase_units(223)      = 1d-12
    is_PD(223) = 0
    nTB(223) = 0

    ! (213):  HNO + H <=> H2 + NO
    fwd_A(224)     = 900000000000d0
    fwd_beta(224)  = 0.71999999999999997d0
    fwd_Ea(224)    = 660d0
    prefactor_units(224)  = 1.0000000000000002d-06
    activation_units(224) = 0.50321666580471969d0
    phase_units(224)      = 1d-12
    is_PD(224) = 0
    nTB(224) = 0

    ! (214):  HNO + OH <=> NO + H2O
    fwd_A(225)     = 13000000d0
    fwd_beta(225)  = 1.8999999999999999d0
    fwd_Ea(225)    = -950d0
    prefactor_units(225)  = 1.0000000000000002d-06
    activation_units(225) = 0.50321666580471969d0
    phase_units(225)      = 1d-12
    is_PD(225) = 0
    nTB(225) = 0

    ! (215):  HNO + O2 <=> HO2 + NO
    fwd_A(226)     = 10000000000000d0
    fwd_beta(226)  = 0d0
    fwd_Ea(226)    = 13000d0
    prefactor_units(226)  = 1.0000000000000002d-06
    activation_units(226) = 0.50321666580471969d0
    phase_units(226)      = 1d-12
    is_PD(226) = 0
    nTB(226) = 0

    ! (216):  CN + O <=> CO + N
    fwd_A(227)     = 77000000000000d0
    fwd_beta(227)  = 0d0
    fwd_Ea(227)    = 0d0
    prefactor_units(227)  = 1.0000000000000002d-06
    activation_units(227) = 0.50321666580471969d0
    phase_units(227)      = 1d-12
    is_PD(227) = 0
    nTB(227) = 0

    ! (217):  CN + OH <=> NCO + H
    fwd_A(228)     = 40000000000000d0
    fwd_beta(228)  = 0d0
    fwd_Ea(228)    = 0d0
    prefactor_units(228)  = 1.0000000000000002d-06
    activation_units(228) = 0.50321666580471969d0
    phase_units(228)      = 1d-12
    is_PD(228) = 0
    nTB(228) = 0

    ! (218):  CN + H2O <=> HCN + OH
    fwd_A(229)     = 8000000000000d0
    fwd_beta(229)  = 0d0
    fwd_Ea(229)    = 7460d0
    prefactor_units(229)  = 1.0000000000000002d-06
    activation_units(229) = 0.50321666580471969d0
    phase_units(229)      = 1d-12
    is_PD(229) = 0
    nTB(229) = 0

    ! (219):  CN + O2 <=> NCO + O
    fwd_A(230)     = 6140000000000d0
    fwd_beta(230)  = 0d0
    fwd_Ea(230)    = -440d0
    prefactor_units(230)  = 1.0000000000000002d-06
    activation_units(230) = 0.50321666580471969d0
    phase_units(230)      = 1d-12
    is_PD(230) = 0
    nTB(230) = 0

    ! (220):  CN + H2 <=> HCN + H
    fwd_A(231)     = 295000d0
    fwd_beta(231)  = 2.4500000000000002d0
    fwd_Ea(231)    = 2240d0
    prefactor_units(231)  = 1.0000000000000002d-06
    activation_units(231) = 0.50321666580471969d0
    phase_units(231)      = 1d-12
    is_PD(231) = 0
    nTB(231) = 0

    ! (221):  NCO + O <=> NO + CO
    fwd_A(232)     = 23500000000000d0
    fwd_beta(232)  = 0d0
    fwd_Ea(232)    = 0d0
    prefactor_units(232)  = 1.0000000000000002d-06
    activation_units(232) = 0.50321666580471969d0
    phase_units(232)      = 1d-12
    is_PD(232) = 0
    nTB(232) = 0

    ! (222):  NCO + H <=> NH + CO
    fwd_A(233)     = 54000000000000d0
    fwd_beta(233)  = 0d0
    fwd_Ea(233)    = 0d0
    prefactor_units(233)  = 1.0000000000000002d-06
    activation_units(233) = 0.50321666580471969d0
    phase_units(233)      = 1d-12
    is_PD(233) = 0
    nTB(233) = 0

    ! (223):  NCO + OH <=> NO + H + CO
    fwd_A(234)     = 2500000000000d0
    fwd_beta(234)  = 0d0
    fwd_Ea(234)    = 0d0
    prefactor_units(234)  = 1.0000000000000002d-06
    activation_units(234) = 0.50321666580471969d0
    phase_units(234)      = 1d-12
    is_PD(234) = 0
    nTB(234) = 0

    ! (224):  NCO + N <=> N2 + CO
    fwd_A(235)     = 20000000000000d0
    fwd_beta(235)  = 0d0
    fwd_Ea(235)    = 0d0
    prefactor_units(235)  = 1.0000000000000002d-06
    activation_units(235) = 0.50321666580471969d0
    phase_units(235)      = 1d-12
    is_PD(235) = 0
    nTB(235) = 0

    ! (225):  NCO + O2 <=> NO + CO2
    fwd_A(236)     = 2000000000000d0
    fwd_beta(236)  = 0d0
    fwd_Ea(236)    = 20000d0
    prefactor_units(236)  = 1.0000000000000002d-06
    activation_units(236) = 0.50321666580471969d0
    phase_units(236)      = 1d-12
    is_PD(236) = 0
    nTB(236) = 0

    ! (226):  NCO + M <=> N + CO + M
    fwd_A(39)     = 310000000000000d0
    fwd_beta(39)  = 0d0
    fwd_Ea(39)    = 54050d0
    prefactor_units(39)  = 1.0000000000000002d-06
    activation_units(39) = 0.50321666580471969d0
    phase_units(39)      = 1d-6
    is_PD(39) = 0
    nTB(39) = 7
    if (.not. allocated(TB(39) % vector)) allocate(TB(39) % vector(7))
    if (.not. allocated(TBid(39) % vector)) allocate(TBid(39) % vector(7))
    TBid(39) % vector(1) = 0d0
    TB(39) % vector(1) = 2d0 ! H2
    TBid(39) % vector(2) = 5d0
    TB(39) % vector(2) = 6d0 ! H2O
    TBid(39) % vector(3) = 13d0
    TB(39) % vector(3) = 2d0 ! CH4
    TBid(39) % vector(4) = 14d0
    TB(39) % vector(4) = 1.5d0 ! CO
    TBid(39) % vector(5) = 15d0
    TB(39) % vector(5) = 2d0 ! CO2
    TBid(39) % vector(6) = 26d0
    TB(39) % vector(6) = 3d0 ! C2H6
    TBid(39) % vector(7) = 48d0
    TB(39) % vector(7) = 0.69999999999999996d0 ! AR

    ! (227):  NCO + NO <=> N2O + CO
    fwd_A(237)     = 1.9d+17
    fwd_beta(237)  = -1.52d0
    fwd_Ea(237)    = 740d0
    prefactor_units(237)  = 1.0000000000000002d-06
    activation_units(237) = 0.50321666580471969d0
    phase_units(237)      = 1d-12
    is_PD(237) = 0
    nTB(237) = 0

    ! (228):  NCO + NO <=> N2 + CO2
    fwd_A(238)     = 3.8d+18
    fwd_beta(238)  = -2d0
    fwd_Ea(238)    = 800d0
    prefactor_units(238)  = 1.0000000000000002d-06
    activation_units(238) = 0.50321666580471969d0
    phase_units(238)      = 1d-12
    is_PD(238) = 0
    nTB(238) = 0

    ! (229):  HCN + M <=> H + CN + M
    fwd_A(40)     = 1.0400000000000001d+29
    fwd_beta(40)  = -3.2999999999999998d0
    fwd_Ea(40)    = 126600d0
    prefactor_units(40)  = 1.0000000000000002d-06
    activation_units(40) = 0.50321666580471969d0
    phase_units(40)      = 1d-6
    is_PD(40) = 0
    nTB(40) = 7
    if (.not. allocated(TB(40) % vector)) allocate(TB(40) % vector(7))
    if (.not. allocated(TBid(40) % vector)) allocate(TBid(40) % vector(7))
    TBid(40) % vector(1) = 0d0
    TB(40) % vector(1) = 2d0 ! H2
    TBid(40) % vector(2) = 5d0
    TB(40) % vector(2) = 6d0 ! H2O
    TBid(40) % vector(3) = 13d0
    TB(40) % vector(3) = 2d0 ! CH4
    TBid(40) % vector(4) = 14d0
    TB(40) % vector(4) = 1.5d0 ! CO
    TBid(40) % vector(5) = 15d0
    TB(40) % vector(5) = 2d0 ! CO2
    TBid(40) % vector(6) = 26d0
    TB(40) % vector(6) = 3d0 ! C2H6
    TBid(40) % vector(7) = 48d0
    TB(40) % vector(7) = 0.69999999999999996d0 ! AR

    ! (230):  HCN + O <=> NCO + H
    fwd_A(239)     = 20300d0
    fwd_beta(239)  = 2.6400000000000001d0
    fwd_Ea(239)    = 4980d0
    prefactor_units(239)  = 1.0000000000000002d-06
    activation_units(239) = 0.50321666580471969d0
    phase_units(239)      = 1d-12
    is_PD(239) = 0
    nTB(239) = 0

    ! (231):  HCN + O <=> NH + CO
    fwd_A(240)     = 5070d0
    fwd_beta(240)  = 2.6400000000000001d0
    fwd_Ea(240)    = 4980d0
    prefactor_units(240)  = 1.0000000000000002d-06
    activation_units(240) = 0.50321666580471969d0
    phase_units(240)      = 1d-12
    is_PD(240) = 0
    nTB(240) = 0

    ! (232):  HCN + O <=> CN + OH
    fwd_A(241)     = 3910000000d0
    fwd_beta(241)  = 1.5800000000000001d0
    fwd_Ea(241)    = 26600d0
    prefactor_units(241)  = 1.0000000000000002d-06
    activation_units(241) = 0.50321666580471969d0
    phase_units(241)      = 1d-12
    is_PD(241) = 0
    nTB(241) = 0

    ! (233):  HCN + OH <=> HOCN + H
    fwd_A(242)     = 1100000d0
    fwd_beta(242)  = 2.0299999999999998d0
    fwd_Ea(242)    = 13370d0
    prefactor_units(242)  = 1.0000000000000002d-06
    activation_units(242) = 0.50321666580471969d0
    phase_units(242)      = 1d-12
    is_PD(242) = 0
    nTB(242) = 0

    ! (234):  HCN + OH <=> HNCO + H
    fwd_A(243)     = 4400d0
    fwd_beta(243)  = 2.2599999999999998d0
    fwd_Ea(243)    = 6400d0
    prefactor_units(243)  = 1.0000000000000002d-06
    activation_units(243) = 0.50321666580471969d0
    phase_units(243)      = 1d-12
    is_PD(243) = 0
    nTB(243) = 0

    ! (235):  HCN + OH <=> NH2 + CO
    fwd_A(244)     = 160d0
    fwd_beta(244)  = 2.5600000000000001d0
    fwd_Ea(244)    = 9000d0
    prefactor_units(244)  = 1.0000000000000002d-06
    activation_units(244) = 0.50321666580471969d0
    phase_units(244)      = 1d-12
    is_PD(244) = 0
    nTB(244) = 0

    ! (236):  H + HCN (+M) <=> H2CN (+M)
    fwd_A(29)     = 33000000000000d0
    fwd_beta(29)  = 0d0
    fwd_Ea(29)    = 0d0
    low_A(29)     = 1.4d+26
    low_beta(29)  = -3.3999999999999999d0
    low_Ea(29)    = 1900d0
    prefactor_units(29)  = 1.0000000000000002d-06
    activation_units(29) = 0.50321666580471969d0
    phase_units(29)      = 1d-12
    is_PD(29) = 1
    nTB(29) = 7
    if (.not. allocated(TB(29) % vector)) allocate(TB(29) % vector(7))
    if (.not. allocated(TBid(29) % vector)) allocate(TBid(29) % vector(7))
    TBid(29) % vector(1) = 0d0
    TB(29) % vector(1) = 2d0 ! H2
    TBid(29) % vector(2) = 5d0
    TB(29) % vector(2) = 6d0 ! H2O
    TBid(29) % vector(3) = 13d0
    TB(29) % vector(3) = 2d0 ! CH4
    TBid(29) % vector(4) = 14d0
    TB(29) % vector(4) = 1.5d0 ! CO
    TBid(29) % vector(5) = 15d0
    TB(29) % vector(5) = 2d0 ! CO2
    TBid(29) % vector(6) = 26d0
    TB(29) % vector(6) = 3d0 ! C2H6
    TBid(29) % vector(7) = 48d0
    TB(29) % vector(7) = 0.69999999999999996d0 ! AR

    ! (237):  H2CN + N <=> N2 + CH2
    fwd_A(245)     = 60000000000000d0
    fwd_beta(245)  = 0d0
    fwd_Ea(245)    = 400d0
    prefactor_units(245)  = 1.0000000000000002d-06
    activation_units(245) = 0.50321666580471969d0
    phase_units(245)      = 1d-12
    is_PD(245) = 0
    nTB(245) = 0

    ! (238):  C + N2 <=> CN + N
    fwd_A(246)     = 63000000000000d0
    fwd_beta(246)  = 0d0
    fwd_Ea(246)    = 46020d0
    prefactor_units(246)  = 1.0000000000000002d-06
    activation_units(246) = 0.50321666580471969d0
    phase_units(246)      = 1d-12
    is_PD(246) = 0
    nTB(246) = 0

    ! (239):  CH + N2 <=> HCN + N
    fwd_A(247)     = 3120000000d0
    fwd_beta(247)  = 0.88d0
    fwd_Ea(247)    = 20130d0
    prefactor_units(247)  = 1.0000000000000002d-06
    activation_units(247) = 0.50321666580471969d0
    phase_units(247)      = 1d-12
    is_PD(247) = 0
    nTB(247) = 0

    ! (240):  CH + N2 (+M) <=> HCNN (+M)
    fwd_A(21)     = 3100000000000d0
    fwd_beta(21)  = 0.14999999999999999d0
    fwd_Ea(21)    = 0d0
    low_A(21)     = 1.2999999999999999d+25
    low_beta(21)  = -3.1600000000000001d0
    low_Ea(21)    = 740d0
    troe_a(21)    = 0.66700000000000004d0
    troe_Tsss(21) = 235d0
    troe_Ts(21)   = 2117d0
    troe_Tss(21)  = 4536d0
    troe_len(21)  = 4
    prefactor_units(21)  = 1.0000000000000002d-06
    activation_units(21) = 0.50321666580471969d0
    phase_units(21)      = 1d-12
    is_PD(21) = 1
    nTB(21) = 7
    if (.not. allocated(TB(21) % vector)) allocate(TB(21) % vector(7))
    if (.not. allocated(TBid(21) % vector)) allocate(TBid(21) % vector(7))
    TBid(21) % vector(1) = 0d0
    TB(21) % vector(1) = 2d0 ! H2
    TBid(21) % vector(2) = 5d0
    TB(21) % vector(2) = 6d0 ! H2O
    TBid(21) % vector(3) = 13d0
    TB(21) % vector(3) = 2d0 ! CH4
    TBid(21) % vector(4) = 14d0
    TB(21) % vector(4) = 1.5d0 ! CO
    TBid(21) % vector(5) = 15d0
    TB(21) % vector(5) = 2d0 ! CO2
    TBid(21) % vector(6) = 26d0
    TB(21) % vector(6) = 3d0 ! C2H6
    TBid(21) % vector(7) = 48d0
    TB(21) % vector(7) = 1d0 ! AR

    ! (241):  CH2 + N2 <=> HCN + NH
    fwd_A(248)     = 10000000000000d0
    fwd_beta(248)  = 0d0
    fwd_Ea(248)    = 74000d0
    prefactor_units(248)  = 1.0000000000000002d-06
    activation_units(248) = 0.50321666580471969d0
    phase_units(248)      = 1d-12
    is_PD(248) = 0
    nTB(248) = 0

    ! (242):  CH2(S) + N2 <=> NH + HCN
    fwd_A(249)     = 100000000000d0
    fwd_beta(249)  = 0d0
    fwd_Ea(249)    = 65000d0
    prefactor_units(249)  = 1.0000000000000002d-06
    activation_units(249) = 0.50321666580471969d0
    phase_units(249)      = 1d-12
    is_PD(249) = 0
    nTB(249) = 0

    ! (243):  C + NO <=> CN + O
    fwd_A(250)     = 19000000000000d0
    fwd_beta(250)  = 0d0
    fwd_Ea(250)    = 0d0
    prefactor_units(250)  = 1.0000000000000002d-06
    activation_units(250) = 0.50321666580471969d0
    phase_units(250)      = 1d-12
    is_PD(250) = 0
    nTB(250) = 0

    ! (244):  C + NO <=> CO + N
    fwd_A(251)     = 29000000000000d0
    fwd_beta(251)  = 0d0
    fwd_Ea(251)    = 0d0
    prefactor_units(251)  = 1.0000000000000002d-06
    activation_units(251) = 0.50321666580471969d0
    phase_units(251)      = 1d-12
    is_PD(251) = 0
    nTB(251) = 0

    ! (245):  CH + NO <=> HCN + O
    fwd_A(252)     = 41000000000000d0
    fwd_beta(252)  = 0d0
    fwd_Ea(252)    = 0d0
    prefactor_units(252)  = 1.0000000000000002d-06
    activation_units(252) = 0.50321666580471969d0
    phase_units(252)      = 1d-12
    is_PD(252) = 0
    nTB(252) = 0

    ! (246):  CH + NO <=> H + NCO
    fwd_A(253)     = 16200000000000d0
    fwd_beta(253)  = 0d0
    fwd_Ea(253)    = 0d0
    prefactor_units(253)  = 1.0000000000000002d-06
    activation_units(253) = 0.50321666580471969d0
    phase_units(253)      = 1d-12
    is_PD(253) = 0
    nTB(253) = 0

    ! (247):  CH + NO <=> N + HCO
    fwd_A(254)     = 24600000000000d0
    fwd_beta(254)  = 0d0
    fwd_Ea(254)    = 0d0
    prefactor_units(254)  = 1.0000000000000002d-06
    activation_units(254) = 0.50321666580471969d0
    phase_units(254)      = 1d-12
    is_PD(254) = 0
    nTB(254) = 0

    ! (248):  CH2 + NO <=> H + HNCO
    fwd_A(255)     = 3.1d+17
    fwd_beta(255)  = -1.3799999999999999d0
    fwd_Ea(255)    = 1270d0
    prefactor_units(255)  = 1.0000000000000002d-06
    activation_units(255) = 0.50321666580471969d0
    phase_units(255)      = 1d-12
    is_PD(255) = 0
    nTB(255) = 0

    ! (249):  CH2 + NO <=> OH + HCN
    fwd_A(256)     = 290000000000000d0
    fwd_beta(256)  = -0.68999999999999995d0
    fwd_Ea(256)    = 760d0
    prefactor_units(256)  = 1.0000000000000002d-06
    activation_units(256) = 0.50321666580471969d0
    phase_units(256)      = 1d-12
    is_PD(256) = 0
    nTB(256) = 0

    ! (250):  CH2 + NO <=> H + HCNO
    fwd_A(257)     = 38000000000000d0
    fwd_beta(257)  = -0.35999999999999999d0
    fwd_Ea(257)    = 580d0
    prefactor_units(257)  = 1.0000000000000002d-06
    activation_units(257) = 0.50321666580471969d0
    phase_units(257)      = 1d-12
    is_PD(257) = 0
    nTB(257) = 0

    ! (251):  CH2(S) + NO <=> H + HNCO
    fwd_A(258)     = 3.1d+17
    fwd_beta(258)  = -1.3799999999999999d0
    fwd_Ea(258)    = 1270d0
    prefactor_units(258)  = 1.0000000000000002d-06
    activation_units(258) = 0.50321666580471969d0
    phase_units(258)      = 1d-12
    is_PD(258) = 0
    nTB(258) = 0

    ! (252):  CH2(S) + NO <=> OH + HCN
    fwd_A(259)     = 290000000000000d0
    fwd_beta(259)  = -0.68999999999999995d0
    fwd_Ea(259)    = 760d0
    prefactor_units(259)  = 1.0000000000000002d-06
    activation_units(259) = 0.50321666580471969d0
    phase_units(259)      = 1d-12
    is_PD(259) = 0
    nTB(259) = 0

    ! (253):  CH2(S) + NO <=> H + HCNO
    fwd_A(260)     = 38000000000000d0
    fwd_beta(260)  = -0.35999999999999999d0
    fwd_Ea(260)    = 580d0
    prefactor_units(260)  = 1.0000000000000002d-06
    activation_units(260) = 0.50321666580471969d0
    phase_units(260)      = 1d-12
    is_PD(260) = 0
    nTB(260) = 0

    ! (254):  CH3 + NO <=> HCN + H2O
    fwd_A(261)     = 96000000000000d0
    fwd_beta(261)  = 0d0
    fwd_Ea(261)    = 28800d0
    prefactor_units(261)  = 1.0000000000000002d-06
    activation_units(261) = 0.50321666580471969d0
    phase_units(261)      = 1d-12
    is_PD(261) = 0
    nTB(261) = 0

    ! (255):  CH3 + NO <=> H2CN + OH
    fwd_A(262)     = 1000000000000d0
    fwd_beta(262)  = 0d0
    fwd_Ea(262)    = 21750d0
    prefactor_units(262)  = 1.0000000000000002d-06
    activation_units(262) = 0.50321666580471969d0
    phase_units(262)      = 1d-12
    is_PD(262) = 0
    nTB(262) = 0

    ! (256):  HCNN + O <=> CO + H + N2
    fwd_A(263)     = 22000000000000d0
    fwd_beta(263)  = 0d0
    fwd_Ea(263)    = 0d0
    prefactor_units(263)  = 1.0000000000000002d-06
    activation_units(263) = 0.50321666580471969d0
    phase_units(263)      = 1d-12
    is_PD(263) = 0
    nTB(263) = 0

    ! (257):  HCNN + O <=> HCN + NO
    fwd_A(264)     = 2000000000000d0
    fwd_beta(264)  = 0d0
    fwd_Ea(264)    = 0d0
    prefactor_units(264)  = 1.0000000000000002d-06
    activation_units(264) = 0.50321666580471969d0
    phase_units(264)      = 1d-12
    is_PD(264) = 0
    nTB(264) = 0

    ! (258):  HCNN + O2 <=> O + HCO + N2
    fwd_A(265)     = 12000000000000d0
    fwd_beta(265)  = 0d0
    fwd_Ea(265)    = 0d0
    prefactor_units(265)  = 1.0000000000000002d-06
    activation_units(265) = 0.50321666580471969d0
    phase_units(265)      = 1d-12
    is_PD(265) = 0
    nTB(265) = 0

    ! (259):  HCNN + OH <=> H + HCO + N2
    fwd_A(266)     = 12000000000000d0
    fwd_beta(266)  = 0d0
    fwd_Ea(266)    = 0d0
    prefactor_units(266)  = 1.0000000000000002d-06
    activation_units(266) = 0.50321666580471969d0
    phase_units(266)      = 1d-12
    is_PD(266) = 0
    nTB(266) = 0

    ! (260):  HCNN + H <=> CH2 + N2
    fwd_A(267)     = 100000000000000d0
    fwd_beta(267)  = 0d0
    fwd_Ea(267)    = 0d0
    prefactor_units(267)  = 1.0000000000000002d-06
    activation_units(267) = 0.50321666580471969d0
    phase_units(267)      = 1d-12
    is_PD(267) = 0
    nTB(267) = 0

    ! (261):  HNCO + O <=> NH + CO2
    fwd_A(268)     = 98000000d0
    fwd_beta(268)  = 1.4099999999999999d0
    fwd_Ea(268)    = 8500d0
    prefactor_units(268)  = 1.0000000000000002d-06
    activation_units(268) = 0.50321666580471969d0
    phase_units(268)      = 1d-12
    is_PD(268) = 0
    nTB(268) = 0

    ! (262):  HNCO + O <=> HNO + CO
    fwd_A(269)     = 150000000d0
    fwd_beta(269)  = 1.5700000000000001d0
    fwd_Ea(269)    = 44000d0
    prefactor_units(269)  = 1.0000000000000002d-06
    activation_units(269) = 0.50321666580471969d0
    phase_units(269)      = 1d-12
    is_PD(269) = 0
    nTB(269) = 0

    ! (263):  HNCO + O <=> NCO + OH
    fwd_A(270)     = 2200000d0
    fwd_beta(270)  = 2.1099999999999999d0
    fwd_Ea(270)    = 11400d0
    prefactor_units(270)  = 1.0000000000000002d-06
    activation_units(270) = 0.50321666580471969d0
    phase_units(270)      = 1d-12
    is_PD(270) = 0
    nTB(270) = 0

    ! (264):  HNCO + H <=> NH2 + CO
    fwd_A(271)     = 22500000d0
    fwd_beta(271)  = 1.7d0
    fwd_Ea(271)    = 3800d0
    prefactor_units(271)  = 1.0000000000000002d-06
    activation_units(271) = 0.50321666580471969d0
    phase_units(271)      = 1d-12
    is_PD(271) = 0
    nTB(271) = 0

    ! (265):  HNCO + H <=> H2 + NCO
    fwd_A(272)     = 105000d0
    fwd_beta(272)  = 2.5d0
    fwd_Ea(272)    = 13300d0
    prefactor_units(272)  = 1.0000000000000002d-06
    activation_units(272) = 0.50321666580471969d0
    phase_units(272)      = 1d-12
    is_PD(272) = 0
    nTB(272) = 0

    ! (266):  HNCO + OH <=> NCO + H2O
    fwd_A(273)     = 33000000d0
    fwd_beta(273)  = 1.5d0
    fwd_Ea(273)    = 3600d0
    prefactor_units(273)  = 1.0000000000000002d-06
    activation_units(273) = 0.50321666580471969d0
    phase_units(273)      = 1d-12
    is_PD(273) = 0
    nTB(273) = 0

    ! (267):  HNCO + OH <=> NH2 + CO2
    fwd_A(274)     = 3300000d0
    fwd_beta(274)  = 1.5d0
    fwd_Ea(274)    = 3600d0
    prefactor_units(274)  = 1.0000000000000002d-06
    activation_units(274) = 0.50321666580471969d0
    phase_units(274)      = 1d-12
    is_PD(274) = 0
    nTB(274) = 0

    ! (268):  HNCO + M <=> NH + CO + M
    fwd_A(41)     = 11800000000000000d0
    fwd_beta(41)  = 0d0
    fwd_Ea(41)    = 84720d0
    prefactor_units(41)  = 1.0000000000000002d-06
    activation_units(41) = 0.50321666580471969d0
    phase_units(41)      = 1d-6
    is_PD(41) = 0
    nTB(41) = 7
    if (.not. allocated(TB(41) % vector)) allocate(TB(41) % vector(7))
    if (.not. allocated(TBid(41) % vector)) allocate(TBid(41) % vector(7))
    TBid(41) % vector(1) = 0d0
    TB(41) % vector(1) = 2d0 ! H2
    TBid(41) % vector(2) = 5d0
    TB(41) % vector(2) = 6d0 ! H2O
    TBid(41) % vector(3) = 13d0
    TB(41) % vector(3) = 2d0 ! CH4
    TBid(41) % vector(4) = 14d0
    TB(41) % vector(4) = 1.5d0 ! CO
    TBid(41) % vector(5) = 15d0
    TB(41) % vector(5) = 2d0 ! CO2
    TBid(41) % vector(6) = 26d0
    TB(41) % vector(6) = 3d0 ! C2H6
    TBid(41) % vector(7) = 48d0
    TB(41) % vector(7) = 0.69999999999999996d0 ! AR

    ! (269):  HCNO + H <=> H + HNCO
    fwd_A(275)     = 2100000000000000d0
    fwd_beta(275)  = -0.68999999999999995d0
    fwd_Ea(275)    = 2850d0
    prefactor_units(275)  = 1.0000000000000002d-06
    activation_units(275) = 0.50321666580471969d0
    phase_units(275)      = 1d-12
    is_PD(275) = 0
    nTB(275) = 0

    ! (270):  HCNO + H <=> OH + HCN
    fwd_A(276)     = 270000000000d0
    fwd_beta(276)  = 0.17999999999999999d0
    fwd_Ea(276)    = 2120d0
    prefactor_units(276)  = 1.0000000000000002d-06
    activation_units(276) = 0.50321666580471969d0
    phase_units(276)      = 1d-12
    is_PD(276) = 0
    nTB(276) = 0

    ! (271):  HCNO + H <=> NH2 + CO
    fwd_A(277)     = 170000000000000d0
    fwd_beta(277)  = -0.75d0
    fwd_Ea(277)    = 2890d0
    prefactor_units(277)  = 1.0000000000000002d-06
    activation_units(277) = 0.50321666580471969d0
    phase_units(277)      = 1d-12
    is_PD(277) = 0
    nTB(277) = 0

    ! (272):  HOCN + H <=> H + HNCO
    fwd_A(278)     = 20000000d0
    fwd_beta(278)  = 2d0
    fwd_Ea(278)    = 2000d0
    prefactor_units(278)  = 1.0000000000000002d-06
    activation_units(278) = 0.50321666580471969d0
    phase_units(278)      = 1d-12
    is_PD(278) = 0
    nTB(278) = 0

    ! (273):  HCCO + NO <=> HCNO + CO
    fwd_A(279)     = 9000000000000d0
    fwd_beta(279)  = 0d0
    fwd_Ea(279)    = 0d0
    prefactor_units(279)  = 1.0000000000000002d-06
    activation_units(279) = 0.50321666580471969d0
    phase_units(279)      = 1d-12
    is_PD(279) = 0
    nTB(279) = 0

    ! (274):  CH3 + N <=> H2CN + H
    fwd_A(280)     = 610000000000000d0
    fwd_beta(280)  = -0.31d0
    fwd_Ea(280)    = 290d0
    prefactor_units(280)  = 1.0000000000000002d-06
    activation_units(280) = 0.50321666580471969d0
    phase_units(280)      = 1d-12
    is_PD(280) = 0
    nTB(280) = 0

    ! (275):  CH3 + N <=> HCN + H2
    fwd_A(281)     = 3700000000000d0
    fwd_beta(281)  = 0.14999999999999999d0
    fwd_Ea(281)    = -90d0
    prefactor_units(281)  = 1.0000000000000002d-06
    activation_units(281) = 0.50321666580471969d0
    phase_units(281)      = 1d-12
    is_PD(281) = 0
    nTB(281) = 0

    ! (276):  NH3 + H <=> NH2 + H2
    fwd_A(282)     = 540000d0
    fwd_beta(282)  = 2.3999999999999999d0
    fwd_Ea(282)    = 9915d0
    prefactor_units(282)  = 1.0000000000000002d-06
    activation_units(282) = 0.50321666580471969d0
    phase_units(282)      = 1d-12
    is_PD(282) = 0
    nTB(282) = 0

    ! (277):  NH3 + OH <=> NH2 + H2O
    fwd_A(283)     = 50000000d0
    fwd_beta(283)  = 1.6000000000000001d0
    fwd_Ea(283)    = 955d0
    prefactor_units(283)  = 1.0000000000000002d-06
    activation_units(283) = 0.50321666580471969d0
    phase_units(283)      = 1d-12
    is_PD(283) = 0
    nTB(283) = 0

    ! (278):  NH3 + O <=> NH2 + OH
    fwd_A(284)     = 9400000d0
    fwd_beta(284)  = 1.9399999999999999d0
    fwd_Ea(284)    = 6460d0
    prefactor_units(284)  = 1.0000000000000002d-06
    activation_units(284) = 0.50321666580471969d0
    phase_units(284)      = 1d-12
    is_PD(284) = 0
    nTB(284) = 0

    ! (279):  NH + CO2 <=> HNO + CO
    fwd_A(285)     = 10000000000000d0
    fwd_beta(285)  = 0d0
    fwd_Ea(285)    = 14350d0
    prefactor_units(285)  = 1.0000000000000002d-06
    activation_units(285) = 0.50321666580471969d0
    phase_units(285)      = 1d-12
    is_PD(285) = 0
    nTB(285) = 0

    ! (280):  CN + NO2 <=> NCO + NO
    fwd_A(286)     = 6160000000000000d0
    fwd_beta(286)  = -0.752d0
    fwd_Ea(286)    = 345d0
    prefactor_units(286)  = 1.0000000000000002d-06
    activation_units(286) = 0.50321666580471969d0
    phase_units(286)      = 1d-12
    is_PD(286) = 0
    nTB(286) = 0

    ! (281):  NCO + NO2 <=> N2O + CO2
    fwd_A(287)     = 3250000000000d0
    fwd_beta(287)  = 0d0
    fwd_Ea(287)    = -705d0
    prefactor_units(287)  = 1.0000000000000002d-06
    activation_units(287) = 0.50321666580471969d0
    phase_units(287)      = 1d-12
    is_PD(287) = 0
    nTB(287) = 0

    ! (282):  N + CO2 <=> NO + CO
    fwd_A(288)     = 3000000000000d0
    fwd_beta(288)  = 0d0
    fwd_Ea(288)    = 11300d0
    prefactor_units(288)  = 1.0000000000000002d-06
    activation_units(288) = 0.50321666580471969d0
    phase_units(288)      = 1d-12
    is_PD(288) = 0
    nTB(288) = 0

    ! (283):  O + CH3 => H + H2 + CO
    fwd_A(289)     = 33700000000000d0
    fwd_beta(289)  = 0d0
    fwd_Ea(289)    = 0d0
    prefactor_units(289)  = 1.0000000000000002d-06
    activation_units(289) = 0.50321666580471969d0
    phase_units(289)      = 1d-12
    is_PD(289) = 0
    nTB(289) = 0

    ! (284):  O + C2H4 <=> H + CH2CHO
    fwd_A(290)     = 6700000d0
    fwd_beta(290)  = 1.8300000000000001d0
    fwd_Ea(290)    = 220d0
    prefactor_units(290)  = 1.0000000000000002d-06
    activation_units(290) = 0.50321666580471969d0
    phase_units(290)      = 1d-12
    is_PD(290) = 0
    nTB(290) = 0

    ! (285):  O + C2H5 <=> H + CH3CHO
    fwd_A(291)     = 109600000000000d0
    fwd_beta(291)  = 0d0
    fwd_Ea(291)    = 0d0
    prefactor_units(291)  = 1.0000000000000002d-06
    activation_units(291) = 0.50321666580471969d0
    phase_units(291)      = 1d-12
    is_PD(291) = 0
    nTB(291) = 0

    ! (286):  OH + HO2 <=> O2 + H2O
    fwd_A(292)     = 5000000000000000d0
    fwd_beta(292)  = 0d0
    fwd_Ea(292)    = 17330d0
    prefactor_units(292)  = 1.0000000000000002d-06
    activation_units(292) = 0.50321666580471969d0
    phase_units(292)      = 1d-12
    is_PD(292) = 0
    nTB(292) = 0

    ! (287):  OH + CH3 => H2 + CH2O
    fwd_A(293)     = 8000000000d0
    fwd_beta(293)  = 0.5d0
    fwd_Ea(293)    = -1755d0
    prefactor_units(293)  = 1.0000000000000002d-06
    activation_units(293) = 0.50321666580471969d0
    phase_units(293)      = 1d-12
    is_PD(293) = 0
    nTB(293) = 0

    ! (288):  CH + H2 (+M) <=> CH3 (+M)
    fwd_A(22)     = 1970000000000d0
    fwd_beta(22)  = 0.42999999999999999d0
    fwd_Ea(22)    = -370d0
    low_A(22)     = 4.82d+25
    low_beta(22)  = -2.7999999999999998d0
    low_Ea(22)    = 590d0
    troe_a(22)    = 0.57799999999999996d0
    troe_Tsss(22) = 122d0
    troe_Ts(22)   = 2535d0
    troe_Tss(22)  = 9365d0
    troe_len(22)  = 4
    prefactor_units(22)  = 1.0000000000000002d-06
    activation_units(22) = 0.50321666580471969d0
    phase_units(22)      = 1d-12
    is_PD(22) = 1
    nTB(22) = 7
    if (.not. allocated(TB(22) % vector)) allocate(TB(22) % vector(7))
    if (.not. allocated(TBid(22) % vector)) allocate(TBid(22) % vector(7))
    TBid(22) % vector(1) = 0d0
    TB(22) % vector(1) = 2d0 ! H2
    TBid(22) % vector(2) = 5d0
    TB(22) % vector(2) = 6d0 ! H2O
    TBid(22) % vector(3) = 13d0
    TB(22) % vector(3) = 2d0 ! CH4
    TBid(22) % vector(4) = 14d0
    TB(22) % vector(4) = 1.5d0 ! CO
    TBid(22) % vector(5) = 15d0
    TB(22) % vector(5) = 2d0 ! CO2
    TBid(22) % vector(6) = 26d0
    TB(22) % vector(6) = 3d0 ! C2H6
    TBid(22) % vector(7) = 48d0
    TB(22) % vector(7) = 0.69999999999999996d0 ! AR

    ! (289):  CH2 + O2 => 2 H + CO2
    fwd_A(294)     = 5800000000000d0
    fwd_beta(294)  = 0d0
    fwd_Ea(294)    = 1500d0
    prefactor_units(294)  = 1.0000000000000002d-06
    activation_units(294) = 0.50321666580471969d0
    phase_units(294)      = 1d-12
    is_PD(294) = 0
    nTB(294) = 0

    ! (290):  CH2 + O2 <=> O + CH2O
    fwd_A(295)     = 2400000000000d0
    fwd_beta(295)  = 0d0
    fwd_Ea(295)    = 1500d0
    prefactor_units(295)  = 1.0000000000000002d-06
    activation_units(295) = 0.50321666580471969d0
    phase_units(295)      = 1d-12
    is_PD(295) = 0
    nTB(295) = 0

    ! (291):  CH2 + CH2 => 2 H + C2H2
    fwd_A(296)     = 200000000000000d0
    fwd_beta(296)  = 0d0
    fwd_Ea(296)    = 10989d0
    prefactor_units(296)  = 1.0000000000000002d-06
    activation_units(296) = 0.50321666580471969d0
    phase_units(296)      = 1d-12
    is_PD(296) = 0
    nTB(296) = 0

    ! (292):  CH2(S) + H2O => H2 + CH2O
    fwd_A(297)     = 68200000000d0
    fwd_beta(297)  = 0.25d0
    fwd_Ea(297)    = -935d0
    prefactor_units(297)  = 1.0000000000000002d-06
    activation_units(297) = 0.50321666580471969d0
    phase_units(297)      = 1d-12
    is_PD(297) = 0
    nTB(297) = 0

    ! (293):  C2H3 + O2 <=> O + CH2CHO
    fwd_A(298)     = 303000000000d0
    fwd_beta(298)  = 0.28999999999999998d0
    fwd_Ea(298)    = 11d0
    prefactor_units(298)  = 1.0000000000000002d-06
    activation_units(298) = 0.50321666580471969d0
    phase_units(298)      = 1d-12
    is_PD(298) = 0
    nTB(298) = 0

    ! (294):  C2H3 + O2 <=> HO2 + C2H2
    fwd_A(299)     = 1337000d0
    fwd_beta(299)  = 1.6100000000000001d0
    fwd_Ea(299)    = -384d0
    prefactor_units(299)  = 1.0000000000000002d-06
    activation_units(299) = 0.50321666580471969d0
    phase_units(299)      = 1d-12
    is_PD(299) = 0
    nTB(299) = 0

    ! (295):  O + CH3CHO <=> OH + CH2CHO
    fwd_A(300)     = 2920000000000d0
    fwd_beta(300)  = 0d0
    fwd_Ea(300)    = 1808d0
    prefactor_units(300)  = 1.0000000000000002d-06
    activation_units(300) = 0.50321666580471969d0
    phase_units(300)      = 1d-12
    is_PD(300) = 0
    nTB(300) = 0

    ! (296):  O + CH3CHO => OH + CH3 + CO
    fwd_A(301)     = 2920000000000d0
    fwd_beta(301)  = 0d0
    fwd_Ea(301)    = 1808d0
    prefactor_units(301)  = 1.0000000000000002d-06
    activation_units(301) = 0.50321666580471969d0
    phase_units(301)      = 1d-12
    is_PD(301) = 0
    nTB(301) = 0

    ! (297):  O2 + CH3CHO => HO2 + CH3 + CO
    fwd_A(302)     = 30100000000000d0
    fwd_beta(302)  = 0d0
    fwd_Ea(302)    = 39150d0
    prefactor_units(302)  = 1.0000000000000002d-06
    activation_units(302) = 0.50321666580471969d0
    phase_units(302)      = 1d-12
    is_PD(302) = 0
    nTB(302) = 0

    ! (298):  H + CH3CHO <=> CH2CHO + H2
    fwd_A(303)     = 2050000000d0
    fwd_beta(303)  = 1.1599999999999999d0
    fwd_Ea(303)    = 2405d0
    prefactor_units(303)  = 1.0000000000000002d-06
    activation_units(303) = 0.50321666580471969d0
    phase_units(303)      = 1d-12
    is_PD(303) = 0
    nTB(303) = 0

    ! (299):  H + CH3CHO => CH3 + H2 + CO
    fwd_A(304)     = 2050000000d0
    fwd_beta(304)  = 1.1599999999999999d0
    fwd_Ea(304)    = 2405d0
    prefactor_units(304)  = 1.0000000000000002d-06
    activation_units(304) = 0.50321666580471969d0
    phase_units(304)      = 1d-12
    is_PD(304) = 0
    nTB(304) = 0

    ! (300):  OH + CH3CHO => CH3 + H2O + CO
    fwd_A(305)     = 23430000000d0
    fwd_beta(305)  = 0.72999999999999998d0
    fwd_Ea(305)    = -1113d0
    prefactor_units(305)  = 1.0000000000000002d-06
    activation_units(305) = 0.50321666580471969d0
    phase_units(305)      = 1d-12
    is_PD(305) = 0
    nTB(305) = 0

    ! (301):  HO2 + CH3CHO => CH3 + H2O2 + CO
    fwd_A(306)     = 3010000000000d0
    fwd_beta(306)  = 0d0
    fwd_Ea(306)    = 11923d0
    prefactor_units(306)  = 1.0000000000000002d-06
    activation_units(306) = 0.50321666580471969d0
    phase_units(306)      = 1d-12
    is_PD(306) = 0
    nTB(306) = 0

    ! (302):  CH3 + CH3CHO => CH3 + CH4 + CO
    fwd_A(307)     = 2720000d0
    fwd_beta(307)  = 1.77d0
    fwd_Ea(307)    = 5920d0
    prefactor_units(307)  = 1.0000000000000002d-06
    activation_units(307) = 0.50321666580471969d0
    phase_units(307)      = 1d-12
    is_PD(307) = 0
    nTB(307) = 0

    ! (303):  H + CH2CO (+M) <=> CH2CHO (+M)
    fwd_A(23)     = 486500000000d0
    fwd_beta(23)  = 0.42199999999999999d0
    fwd_Ea(23)    = -1755d0
    low_A(23)     = 1.012d+42
    low_beta(23)  = -7.6299999999999999d0
    low_Ea(23)    = 3854d0
    troe_a(23)    = 0.46500000000000002d0
    troe_Tsss(23) = 201d0
    troe_Ts(23)   = 1773d0
    troe_Tss(23)  = 5333d0
    troe_len(23)  = 4
    prefactor_units(23)  = 1.0000000000000002d-06
    activation_units(23) = 0.50321666580471969d0
    phase_units(23)      = 1d-12
    is_PD(23) = 1
    nTB(23) = 7
    if (.not. allocated(TB(23) % vector)) allocate(TB(23) % vector(7))
    if (.not. allocated(TBid(23) % vector)) allocate(TBid(23) % vector(7))
    TBid(23) % vector(1) = 0d0
    TB(23) % vector(1) = 2d0 ! H2
    TBid(23) % vector(2) = 5d0
    TB(23) % vector(2) = 6d0 ! H2O
    TBid(23) % vector(3) = 13d0
    TB(23) % vector(3) = 2d0 ! CH4
    TBid(23) % vector(4) = 14d0
    TB(23) % vector(4) = 1.5d0 ! CO
    TBid(23) % vector(5) = 15d0
    TB(23) % vector(5) = 2d0 ! CO2
    TBid(23) % vector(6) = 26d0
    TB(23) % vector(6) = 3d0 ! C2H6
    TBid(23) % vector(7) = 48d0
    TB(23) % vector(7) = 0.69999999999999996d0 ! AR

    ! (304):  O + CH2CHO => H + CH2 + CO2
    fwd_A(308)     = 150000000000000d0
    fwd_beta(308)  = 0d0
    fwd_Ea(308)    = 0d0
    prefactor_units(308)  = 1.0000000000000002d-06
    activation_units(308) = 0.50321666580471969d0
    phase_units(308)      = 1d-12
    is_PD(308) = 0
    nTB(308) = 0

    ! (305):  O2 + CH2CHO => OH + CO + CH2O
    fwd_A(309)     = 18100000000d0
    fwd_beta(309)  = 0d0
    fwd_Ea(309)    = 0d0
    prefactor_units(309)  = 1.0000000000000002d-06
    activation_units(309) = 0.50321666580471969d0
    phase_units(309)      = 1d-12
    is_PD(309) = 0
    nTB(309) = 0

    ! (306):  O2 + CH2CHO => OH + 2 HCO
    fwd_A(310)     = 23500000000d0
    fwd_beta(310)  = 0d0
    fwd_Ea(310)    = 0d0
    prefactor_units(310)  = 1.0000000000000002d-06
    activation_units(310) = 0.50321666580471969d0
    phase_units(310)      = 1d-12
    is_PD(310) = 0
    nTB(310) = 0

    ! (307):  H + CH2CHO <=> CH3 + HCO
    fwd_A(311)     = 22000000000000d0
    fwd_beta(311)  = 0d0
    fwd_Ea(311)    = 0d0
    prefactor_units(311)  = 1.0000000000000002d-06
    activation_units(311) = 0.50321666580471969d0
    phase_units(311)      = 1d-12
    is_PD(311) = 0
    nTB(311) = 0

    ! (308):  H + CH2CHO <=> CH2CO + H2
    fwd_A(312)     = 11000000000000d0
    fwd_beta(312)  = 0d0
    fwd_Ea(312)    = 0d0
    prefactor_units(312)  = 1.0000000000000002d-06
    activation_units(312) = 0.50321666580471969d0
    phase_units(312)      = 1d-12
    is_PD(312) = 0
    nTB(312) = 0

    ! (309):  OH + CH2CHO <=> H2O + CH2CO
    fwd_A(313)     = 12000000000000d0
    fwd_beta(313)  = 0d0
    fwd_Ea(313)    = 0d0
    prefactor_units(313)  = 1.0000000000000002d-06
    activation_units(313) = 0.50321666580471969d0
    phase_units(313)      = 1d-12
    is_PD(313) = 0
    nTB(313) = 0

    ! (310):  OH + CH2CHO <=> HCO + CH2OH
    fwd_A(314)     = 30100000000000d0
    fwd_beta(314)  = 0d0
    fwd_Ea(314)    = 0d0
    prefactor_units(314)  = 1.0000000000000002d-06
    activation_units(314) = 0.50321666580471969d0
    phase_units(314)      = 1d-12
    is_PD(314) = 0
    nTB(314) = 0

    ! (311):  CH3 + C2H5 (+M) <=> C3H8 (+M)
    fwd_A(24)     = 9430000000000d0
    fwd_beta(24)  = 0d0
    fwd_Ea(24)    = 0d0
    low_A(24)     = 2.71d+74
    low_beta(24)  = -16.82d0
    low_Ea(24)    = 13065d0
    troe_a(24)    = 0.1527d0
    troe_Tsss(24) = 291d0
    troe_Ts(24)   = 2742d0
    troe_Tss(24)  = 7748d0
    troe_len(24)  = 4
    prefactor_units(24)  = 1.0000000000000002d-06
    activation_units(24) = 0.50321666580471969d0
    phase_units(24)      = 1d-12
    is_PD(24) = 1
    nTB(24) = 7
    if (.not. allocated(TB(24) % vector)) allocate(TB(24) % vector(7))
    if (.not. allocated(TBid(24) % vector)) allocate(TBid(24) % vector(7))
    TBid(24) % vector(1) = 0d0
    TB(24) % vector(1) = 2d0 ! H2
    TBid(24) % vector(2) = 5d0
    TB(24) % vector(2) = 6d0 ! H2O
    TBid(24) % vector(3) = 13d0
    TB(24) % vector(3) = 2d0 ! CH4
    TBid(24) % vector(4) = 14d0
    TB(24) % vector(4) = 1.5d0 ! CO
    TBid(24) % vector(5) = 15d0
    TB(24) % vector(5) = 2d0 ! CO2
    TBid(24) % vector(6) = 26d0
    TB(24) % vector(6) = 3d0 ! C2H6
    TBid(24) % vector(7) = 48d0
    TB(24) % vector(7) = 0.69999999999999996d0 ! AR

    ! (312):  O + C3H8 <=> OH + C3H7
    fwd_A(315)     = 193000d0
    fwd_beta(315)  = 2.6800000000000002d0
    fwd_Ea(315)    = 3716d0
    prefactor_units(315)  = 1.0000000000000002d-06
    activation_units(315) = 0.50321666580471969d0
    phase_units(315)      = 1d-12
    is_PD(315) = 0
    nTB(315) = 0

    ! (313):  H + C3H8 <=> C3H7 + H2
    fwd_A(316)     = 1320000d0
    fwd_beta(316)  = 2.54d0
    fwd_Ea(316)    = 6756d0
    prefactor_units(316)  = 1.0000000000000002d-06
    activation_units(316) = 0.50321666580471969d0
    phase_units(316)      = 1d-12
    is_PD(316) = 0
    nTB(316) = 0

    ! (314):  OH + C3H8 <=> C3H7 + H2O
    fwd_A(317)     = 31600000d0
    fwd_beta(317)  = 1.8d0
    fwd_Ea(317)    = 934d0
    prefactor_units(317)  = 1.0000000000000002d-06
    activation_units(317) = 0.50321666580471969d0
    phase_units(317)      = 1d-12
    is_PD(317) = 0
    nTB(317) = 0

    ! (315):  C3H7 + H2O2 <=> HO2 + C3H8
    fwd_A(318)     = 378d0
    fwd_beta(318)  = 2.7200000000000002d0
    fwd_Ea(318)    = 1500d0
    prefactor_units(318)  = 1.0000000000000002d-06
    activation_units(318) = 0.50321666580471969d0
    phase_units(318)      = 1d-12
    is_PD(318) = 0
    nTB(318) = 0

    ! (316):  CH3 + C3H8 <=> C3H7 + CH4
    fwd_A(319)     = 0.90300000000000002d0
    fwd_beta(319)  = 3.6499999999999999d0
    fwd_Ea(319)    = 7154d0
    prefactor_units(319)  = 1.0000000000000002d-06
    activation_units(319) = 0.50321666580471969d0
    phase_units(319)      = 1d-12
    is_PD(319) = 0
    nTB(319) = 0

    ! (317):  CH3 + C2H4 (+M) <=> C3H7 (+M)
    fwd_A(25)     = 2550000d0
    fwd_beta(25)  = 1.6000000000000001d0
    fwd_Ea(25)    = 5700d0
    low_A(25)     = 3d+63
    low_beta(25)  = -14.6d0
    low_Ea(25)    = 18170d0
    troe_a(25)    = 0.18940000000000001d0
    troe_Tsss(25) = 277d0
    troe_Ts(25)   = 8748d0
    troe_Tss(25)  = 7891d0
    troe_len(25)  = 4
    prefactor_units(25)  = 1.0000000000000002d-06
    activation_units(25) = 0.50321666580471969d0
    phase_units(25)      = 1d-12
    is_PD(25) = 1
    nTB(25) = 7
    if (.not. allocated(TB(25) % vector)) allocate(TB(25) % vector(7))
    if (.not. allocated(TBid(25) % vector)) allocate(TBid(25) % vector(7))
    TBid(25) % vector(1) = 0d0
    TB(25) % vector(1) = 2d0 ! H2
    TBid(25) % vector(2) = 5d0
    TB(25) % vector(2) = 6d0 ! H2O
    TBid(25) % vector(3) = 13d0
    TB(25) % vector(3) = 2d0 ! CH4
    TBid(25) % vector(4) = 14d0
    TB(25) % vector(4) = 1.5d0 ! CO
    TBid(25) % vector(5) = 15d0
    TB(25) % vector(5) = 2d0 ! CO2
    TBid(25) % vector(6) = 26d0
    TB(25) % vector(6) = 3d0 ! C2H6
    TBid(25) % vector(7) = 48d0
    TB(25) % vector(7) = 0.69999999999999996d0 ! AR

    ! (318):  O + C3H7 <=> C2H5 + CH2O
    fwd_A(320)     = 96400000000000d0
    fwd_beta(320)  = 0d0
    fwd_Ea(320)    = 0d0
    prefactor_units(320)  = 1.0000000000000002d-06
    activation_units(320) = 0.50321666580471969d0
    phase_units(320)      = 1d-12
    is_PD(320) = 0
    nTB(320) = 0

    ! (319):  H + C3H7 (+M) <=> C3H8 (+M)
    fwd_A(26)     = 36130000000000d0
    fwd_beta(26)  = 0d0
    fwd_Ea(26)    = 0d0
    low_A(26)     = 4.4199999999999998d+61
    low_beta(26)  = -13.545d0
    low_Ea(26)    = 11357d0
    troe_a(26)    = 0.315d0
    troe_Tsss(26) = 369d0
    troe_Ts(26)   = 3285d0
    troe_Tss(26)  = 6667d0
    troe_len(26)  = 4
    prefactor_units(26)  = 1.0000000000000002d-06
    activation_units(26) = 0.50321666580471969d0
    phase_units(26)      = 1d-12
    is_PD(26) = 1
    nTB(26) = 7
    if (.not. allocated(TB(26) % vector)) allocate(TB(26) % vector(7))
    if (.not. allocated(TBid(26) % vector)) allocate(TBid(26) % vector(7))
    TBid(26) % vector(1) = 0d0
    TB(26) % vector(1) = 2d0 ! H2
    TBid(26) % vector(2) = 5d0
    TB(26) % vector(2) = 6d0 ! H2O
    TBid(26) % vector(3) = 13d0
    TB(26) % vector(3) = 2d0 ! CH4
    TBid(26) % vector(4) = 14d0
    TB(26) % vector(4) = 1.5d0 ! CO
    TBid(26) % vector(5) = 15d0
    TB(26) % vector(5) = 2d0 ! CO2
    TBid(26) % vector(6) = 26d0
    TB(26) % vector(6) = 3d0 ! C2H6
    TBid(26) % vector(7) = 48d0
    TB(26) % vector(7) = 0.69999999999999996d0 ! AR

    ! (320):  H + C3H7 <=> CH3 + C2H5
    fwd_A(321)     = 4060000d0
    fwd_beta(321)  = 2.1899999999999999d0
    fwd_Ea(321)    = 890d0
    prefactor_units(321)  = 1.0000000000000002d-06
    activation_units(321) = 0.50321666580471969d0
    phase_units(321)      = 1d-12
    is_PD(321) = 0
    nTB(321) = 0

    ! (321):  OH + C3H7 <=> C2H5 + CH2OH
    fwd_A(322)     = 24100000000000d0
    fwd_beta(322)  = 0d0
    fwd_Ea(322)    = 0d0
    prefactor_units(322)  = 1.0000000000000002d-06
    activation_units(322) = 0.50321666580471969d0
    phase_units(322)      = 1d-12
    is_PD(322) = 0
    nTB(322) = 0

    ! (322):  HO2 + C3H7 <=> O2 + C3H8
    fwd_A(323)     = 25500000000d0
    fwd_beta(323)  = 0.255d0
    fwd_Ea(323)    = -943d0
    prefactor_units(323)  = 1.0000000000000002d-06
    activation_units(323) = 0.50321666580471969d0
    phase_units(323)      = 1d-12
    is_PD(323) = 0
    nTB(323) = 0

    ! (323):  HO2 + C3H7 => OH + C2H5 + CH2O
    fwd_A(324)     = 24100000000000d0
    fwd_beta(324)  = 0d0
    fwd_Ea(324)    = 0d0
    prefactor_units(324)  = 1.0000000000000002d-06
    activation_units(324) = 0.50321666580471969d0
    phase_units(324)      = 1d-12
    is_PD(324) = 0
    nTB(324) = 0

    ! (324):  CH3 + C3H7 <=> 2 C2H5
    fwd_A(325)     = 19270000000000d0
    fwd_beta(325)  = -0.32000000000000001d0
    fwd_Ea(325)    = 0d0
    prefactor_units(325)  = 1.0000000000000002d-06
    activation_units(325) = 0.50321666580471969d0
    phase_units(325)      = 1d-12
    is_PD(325) = 0
    nTB(325) = 0

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
    kk = 53
    ii = 325
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

    integer, intent(out) :: kname(plenkname*53)
    integer, intent(in) :: plenkname

    integer :: i
    integer :: lenkname

    lenkname = plenkname

    !clear kname
    do i=1, lenkname*53
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
    ! H2O2 
    kname(7*lenkname+1) = ichar('H')
    kname(7*lenkname+2) = ichar('2')
    kname(7*lenkname+3) = ichar('O')
    kname(7*lenkname+4) = ichar('2')
    kname(7*lenkname+5) = ichar(' ')
    ! C 
    kname(8*lenkname+1) = ichar('C')
    kname(8*lenkname+2) = ichar(' ')
    ! CH 
    kname(9*lenkname+1) = ichar('C')
    kname(9*lenkname+2) = ichar('H')
    kname(9*lenkname+3) = ichar(' ')
    ! CH2 
    kname(10*lenkname+1) = ichar('C')
    kname(10*lenkname+2) = ichar('H')
    kname(10*lenkname+3) = ichar('2')
    kname(10*lenkname+4) = ichar(' ')
    ! CH2(S) 
    kname(11*lenkname+1) = ichar('C')
    kname(11*lenkname+2) = ichar('H')
    kname(11*lenkname+3) = ichar('2')
    kname(11*lenkname+4) = ichar('(')
    kname(11*lenkname+5) = ichar('S')
    kname(11*lenkname+6) = ichar(')')
    kname(11*lenkname+7) = ichar(' ')
    ! CH3 
    kname(12*lenkname+1) = ichar('C')
    kname(12*lenkname+2) = ichar('H')
    kname(12*lenkname+3) = ichar('3')
    kname(12*lenkname+4) = ichar(' ')
    ! CH4 
    kname(13*lenkname+1) = ichar('C')
    kname(13*lenkname+2) = ichar('H')
    kname(13*lenkname+3) = ichar('4')
    kname(13*lenkname+4) = ichar(' ')
    ! CO 
    kname(14*lenkname+1) = ichar('C')
    kname(14*lenkname+2) = ichar('O')
    kname(14*lenkname+3) = ichar(' ')
    ! CO2 
    kname(15*lenkname+1) = ichar('C')
    kname(15*lenkname+2) = ichar('O')
    kname(15*lenkname+3) = ichar('2')
    kname(15*lenkname+4) = ichar(' ')
    ! HCO 
    kname(16*lenkname+1) = ichar('H')
    kname(16*lenkname+2) = ichar('C')
    kname(16*lenkname+3) = ichar('O')
    kname(16*lenkname+4) = ichar(' ')
    ! CH2O 
    kname(17*lenkname+1) = ichar('C')
    kname(17*lenkname+2) = ichar('H')
    kname(17*lenkname+3) = ichar('2')
    kname(17*lenkname+4) = ichar('O')
    kname(17*lenkname+5) = ichar(' ')
    ! CH2OH 
    kname(18*lenkname+1) = ichar('C')
    kname(18*lenkname+2) = ichar('H')
    kname(18*lenkname+3) = ichar('2')
    kname(18*lenkname+4) = ichar('O')
    kname(18*lenkname+5) = ichar('H')
    kname(18*lenkname+6) = ichar(' ')
    ! CH3O 
    kname(19*lenkname+1) = ichar('C')
    kname(19*lenkname+2) = ichar('H')
    kname(19*lenkname+3) = ichar('3')
    kname(19*lenkname+4) = ichar('O')
    kname(19*lenkname+5) = ichar(' ')
    ! CH3OH 
    kname(20*lenkname+1) = ichar('C')
    kname(20*lenkname+2) = ichar('H')
    kname(20*lenkname+3) = ichar('3')
    kname(20*lenkname+4) = ichar('O')
    kname(20*lenkname+5) = ichar('H')
    kname(20*lenkname+6) = ichar(' ')
    ! C2H 
    kname(21*lenkname+1) = ichar('C')
    kname(21*lenkname+2) = ichar('2')
    kname(21*lenkname+3) = ichar('H')
    kname(21*lenkname+4) = ichar(' ')
    ! C2H2 
    kname(22*lenkname+1) = ichar('C')
    kname(22*lenkname+2) = ichar('2')
    kname(22*lenkname+3) = ichar('H')
    kname(22*lenkname+4) = ichar('2')
    kname(22*lenkname+5) = ichar(' ')
    ! C2H3 
    kname(23*lenkname+1) = ichar('C')
    kname(23*lenkname+2) = ichar('2')
    kname(23*lenkname+3) = ichar('H')
    kname(23*lenkname+4) = ichar('3')
    kname(23*lenkname+5) = ichar(' ')
    ! C2H4 
    kname(24*lenkname+1) = ichar('C')
    kname(24*lenkname+2) = ichar('2')
    kname(24*lenkname+3) = ichar('H')
    kname(24*lenkname+4) = ichar('4')
    kname(24*lenkname+5) = ichar(' ')
    ! C2H5 
    kname(25*lenkname+1) = ichar('C')
    kname(25*lenkname+2) = ichar('2')
    kname(25*lenkname+3) = ichar('H')
    kname(25*lenkname+4) = ichar('5')
    kname(25*lenkname+5) = ichar(' ')
    ! C2H6 
    kname(26*lenkname+1) = ichar('C')
    kname(26*lenkname+2) = ichar('2')
    kname(26*lenkname+3) = ichar('H')
    kname(26*lenkname+4) = ichar('6')
    kname(26*lenkname+5) = ichar(' ')
    ! HCCO 
    kname(27*lenkname+1) = ichar('H')
    kname(27*lenkname+2) = ichar('C')
    kname(27*lenkname+3) = ichar('C')
    kname(27*lenkname+4) = ichar('O')
    kname(27*lenkname+5) = ichar(' ')
    ! CH2CO 
    kname(28*lenkname+1) = ichar('C')
    kname(28*lenkname+2) = ichar('H')
    kname(28*lenkname+3) = ichar('2')
    kname(28*lenkname+4) = ichar('C')
    kname(28*lenkname+5) = ichar('O')
    kname(28*lenkname+6) = ichar(' ')
    ! HCCOH 
    kname(29*lenkname+1) = ichar('H')
    kname(29*lenkname+2) = ichar('C')
    kname(29*lenkname+3) = ichar('C')
    kname(29*lenkname+4) = ichar('O')
    kname(29*lenkname+5) = ichar('H')
    kname(29*lenkname+6) = ichar(' ')
    ! N 
    kname(30*lenkname+1) = ichar('N')
    kname(30*lenkname+2) = ichar(' ')
    ! NH 
    kname(31*lenkname+1) = ichar('N')
    kname(31*lenkname+2) = ichar('H')
    kname(31*lenkname+3) = ichar(' ')
    ! NH2 
    kname(32*lenkname+1) = ichar('N')
    kname(32*lenkname+2) = ichar('H')
    kname(32*lenkname+3) = ichar('2')
    kname(32*lenkname+4) = ichar(' ')
    ! NH3 
    kname(33*lenkname+1) = ichar('N')
    kname(33*lenkname+2) = ichar('H')
    kname(33*lenkname+3) = ichar('3')
    kname(33*lenkname+4) = ichar(' ')
    ! NNH 
    kname(34*lenkname+1) = ichar('N')
    kname(34*lenkname+2) = ichar('N')
    kname(34*lenkname+3) = ichar('H')
    kname(34*lenkname+4) = ichar(' ')
    ! NO 
    kname(35*lenkname+1) = ichar('N')
    kname(35*lenkname+2) = ichar('O')
    kname(35*lenkname+3) = ichar(' ')
    ! NO2 
    kname(36*lenkname+1) = ichar('N')
    kname(36*lenkname+2) = ichar('O')
    kname(36*lenkname+3) = ichar('2')
    kname(36*lenkname+4) = ichar(' ')
    ! N2O 
    kname(37*lenkname+1) = ichar('N')
    kname(37*lenkname+2) = ichar('2')
    kname(37*lenkname+3) = ichar('O')
    kname(37*lenkname+4) = ichar(' ')
    ! HNO 
    kname(38*lenkname+1) = ichar('H')
    kname(38*lenkname+2) = ichar('N')
    kname(38*lenkname+3) = ichar('O')
    kname(38*lenkname+4) = ichar(' ')
    ! CN 
    kname(39*lenkname+1) = ichar('C')
    kname(39*lenkname+2) = ichar('N')
    kname(39*lenkname+3) = ichar(' ')
    ! HCN 
    kname(40*lenkname+1) = ichar('H')
    kname(40*lenkname+2) = ichar('C')
    kname(40*lenkname+3) = ichar('N')
    kname(40*lenkname+4) = ichar(' ')
    ! H2CN 
    kname(41*lenkname+1) = ichar('H')
    kname(41*lenkname+2) = ichar('2')
    kname(41*lenkname+3) = ichar('C')
    kname(41*lenkname+4) = ichar('N')
    kname(41*lenkname+5) = ichar(' ')
    ! HCNN 
    kname(42*lenkname+1) = ichar('H')
    kname(42*lenkname+2) = ichar('C')
    kname(42*lenkname+3) = ichar('N')
    kname(42*lenkname+4) = ichar('N')
    kname(42*lenkname+5) = ichar(' ')
    ! HCNO 
    kname(43*lenkname+1) = ichar('H')
    kname(43*lenkname+2) = ichar('C')
    kname(43*lenkname+3) = ichar('N')
    kname(43*lenkname+4) = ichar('O')
    kname(43*lenkname+5) = ichar(' ')
    ! HOCN 
    kname(44*lenkname+1) = ichar('H')
    kname(44*lenkname+2) = ichar('O')
    kname(44*lenkname+3) = ichar('C')
    kname(44*lenkname+4) = ichar('N')
    kname(44*lenkname+5) = ichar(' ')
    ! HNCO 
    kname(45*lenkname+1) = ichar('H')
    kname(45*lenkname+2) = ichar('N')
    kname(45*lenkname+3) = ichar('C')
    kname(45*lenkname+4) = ichar('O')
    kname(45*lenkname+5) = ichar(' ')
    ! NCO 
    kname(46*lenkname+1) = ichar('N')
    kname(46*lenkname+2) = ichar('C')
    kname(46*lenkname+3) = ichar('O')
    kname(46*lenkname+4) = ichar(' ')
    ! N2 
    kname(47*lenkname+1) = ichar('N')
    kname(47*lenkname+2) = ichar('2')
    kname(47*lenkname+3) = ichar(' ')
    ! AR 
    kname(48*lenkname+1) = ichar('A')
    kname(48*lenkname+2) = ichar('R')
    kname(48*lenkname+3) = ichar(' ')
    ! C3H7 
    kname(49*lenkname+1) = ichar('C')
    kname(49*lenkname+2) = ichar('3')
    kname(49*lenkname+3) = ichar('H')
    kname(49*lenkname+4) = ichar('7')
    kname(49*lenkname+5) = ichar(' ')
    ! C3H8 
    kname(50*lenkname+1) = ichar('C')
    kname(50*lenkname+2) = ichar('3')
    kname(50*lenkname+3) = ichar('H')
    kname(50*lenkname+4) = ichar('8')
    kname(50*lenkname+5) = ichar(' ')
    ! CH2CHO 
    kname(51*lenkname+1) = ichar('C')
    kname(51*lenkname+2) = ichar('H')
    kname(51*lenkname+3) = ichar('2')
    kname(51*lenkname+4) = ichar('C')
    kname(51*lenkname+5) = ichar('H')
    kname(51*lenkname+6) = ichar('O')
    kname(51*lenkname+7) = ichar(' ')
    ! CH3CHO 
    kname(52*lenkname+1) = ichar('C')
    kname(52*lenkname+2) = ichar('H')
    kname(52*lenkname+3) = ichar('3')
    kname(52*lenkname+4) = ichar('C')
    kname(52*lenkname+5) = ichar('H')
    kname(52*lenkname+6) = ichar('O')
    kname(52*lenkname+7) = ichar(' ')

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
    double precision, intent(in) :: y(53)
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
    YOW = YOW + (y(8) * imw(8)) ! H2O2
    YOW = YOW + (y(9) * imw(9)) ! C
    YOW = YOW + (y(10) * imw(10)) ! CH
    YOW = YOW + (y(11) * imw(11)) ! CH2
    YOW = YOW + (y(12) * imw(12)) ! CH2(S)
    YOW = YOW + (y(13) * imw(13)) ! CH3
    YOW = YOW + (y(14) * imw(14)) ! CH4
    YOW = YOW + (y(15) * imw(15)) ! CO
    YOW = YOW + (y(16) * imw(16)) ! CO2
    YOW = YOW + (y(17) * imw(17)) ! HCO
    YOW = YOW + (y(18) * imw(18)) ! CH2O
    YOW = YOW + (y(19) * imw(19)) ! CH2OH
    YOW = YOW + (y(20) * imw(20)) ! CH3O
    YOW = YOW + (y(21) * imw(21)) ! CH3OH
    YOW = YOW + (y(22) * imw(22)) ! C2H
    YOW = YOW + (y(23) * imw(23)) ! C2H2
    YOW = YOW + (y(24) * imw(24)) ! C2H3
    YOW = YOW + (y(25) * imw(25)) ! C2H4
    YOW = YOW + (y(26) * imw(26)) ! C2H5
    YOW = YOW + (y(27) * imw(27)) ! C2H6
    YOW = YOW + (y(28) * imw(28)) ! HCCO
    YOW = YOW + (y(29) * imw(29)) ! CH2CO
    YOW = YOW + (y(30) * imw(30)) ! HCCOH
    YOW = YOW + (y(31) * imw(31)) ! N
    YOW = YOW + (y(32) * imw(32)) ! NH
    YOW = YOW + (y(33) * imw(33)) ! NH2
    YOW = YOW + (y(34) * imw(34)) ! NH3
    YOW = YOW + (y(35) * imw(35)) ! NNH
    YOW = YOW + (y(36) * imw(36)) ! NO
    YOW = YOW + (y(37) * imw(37)) ! NO2
    YOW = YOW + (y(38) * imw(38)) ! N2O
    YOW = YOW + (y(39) * imw(39)) ! HNO
    YOW = YOW + (y(40) * imw(40)) ! CN
    YOW = YOW + (y(41) * imw(41)) ! HCN
    YOW = YOW + (y(42) * imw(42)) ! H2CN
    YOW = YOW + (y(43) * imw(43)) ! HCNN
    YOW = YOW + (y(44) * imw(44)) ! HCNO
    YOW = YOW + (y(45) * imw(45)) ! HOCN
    YOW = YOW + (y(46) * imw(46)) ! HNCO
    YOW = YOW + (y(47) * imw(47)) ! NCO
    YOW = YOW + (y(48) * imw(48)) ! N2
    YOW = YOW + (y(49) * imw(49)) ! AR
    YOW = YOW + (y(50) * imw(50)) ! C3H7
    YOW = YOW + (y(51) * imw(51)) ! C3H8
    YOW = YOW + (y(52) * imw(52)) ! CH2CHO
    YOW = YOW + (y(53) * imw(53)) ! CH3CHO

    ! YOW holds the reciprocal of the mean molecular wt
    P = rho * 8.31451000d+07 * T * YOW ! P = rho*R*T/W

end subroutine

! Compute rho = P*W(y)/RT
subroutine ckrhoy(P, T, y, iwrk, rwrk, rho)

    implicit none

    double precision, intent(in) :: P
    double precision, intent(in) :: T
    double precision, intent(in) :: y(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: rho

    double precision :: YOW, tmp(53)
    integer :: i

    YOW = 0.d0

    do i=1, 53
        tmp(i) = y(i) * imw(i)
    end do
    do i=1, 53
        YOW = YOW + tmp(i)
    end do

    rho = P / ( 8.31451000d+07 * T * YOW) ! rho = P*W/(R*T)

end subroutine

! get molecular weight for all species
subroutine ckwt(iwrk, rwrk, wt)

    implicit none

    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: wt(53)

    call molecularWeight(wt)

end subroutine

! convert y[species] (mass fracs) to x[species] (mole fracs)
subroutine ckytx(y, iwrk, rwrk, x)

    implicit none

    double precision, intent(in) :: y(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: x(53)

    double precision :: YOW, YOWINV
    integer :: i

    YOW = 0.d0

    do i=1, 53
        YOW = YOW + y(i) * imw(i)
    end do

    YOWINV = 1.d0 / YOW

    do i=1, 53
        x(i) = y(i) * imw(i) * YOWINV
    end do

end subroutine

! convert y(npoints,species) (mass fracs) to x(npoints,species) (mole fracs)
subroutine vckytx(np, y, iwrk, rwrk, x)

    implicit none

    integer, intent(in) :: np
    double precision, intent(in) :: y(np,53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: x(np,53)

    double precision :: YOW(np)
    integer :: i, n

    do i=1, np
        YOW(i) = 0.d0
    end do

    do n=1, 53
        do i=1, np
            x(i,n) = y(i,n) * imw(n)
            YOW(i) = YOW(i) + x(i,n)
        end do
    end do

    do i=1, np
        YOW(i) = 1.d0 / YOW(i)
    end do

    do n=1, 53
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
    double precision, intent(in) :: y(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: c(53)

    integer :: i

    do i=1, 53
        c(i) = rho * y(i) * imw(i)
    end do

end subroutine

! convert x[species] (mole fracs) to y[species] (mass fracs)
subroutine ckxty(x, iwrk, rwrk, y)

    implicit none

    double precision, intent(in) :: x(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: y(53)

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
    XW = XW + (x(8) * 3.40147400d+01) ! H2O2
    XW = XW + (x(9) * 1.20111500d+01) ! C
    XW = XW + (x(10) * 1.30191200d+01) ! CH
    XW = XW + (x(11) * 1.40270900d+01) ! CH2
    XW = XW + (x(12) * 1.40270900d+01) ! CH2(S)
    XW = XW + (x(13) * 1.50350600d+01) ! CH3
    XW = XW + (x(14) * 1.60430300d+01) ! CH4
    XW = XW + (x(15) * 2.80105500d+01) ! CO
    XW = XW + (x(16) * 4.40099500d+01) ! CO2
    XW = XW + (x(17) * 2.90185200d+01) ! HCO
    XW = XW + (x(18) * 3.00264900d+01) ! CH2O
    XW = XW + (x(19) * 3.10344600d+01) ! CH2OH
    XW = XW + (x(20) * 3.10344600d+01) ! CH3O
    XW = XW + (x(21) * 3.20424300d+01) ! CH3OH
    XW = XW + (x(22) * 2.50302700d+01) ! C2H
    XW = XW + (x(23) * 2.60382400d+01) ! C2H2
    XW = XW + (x(24) * 2.70462100d+01) ! C2H3
    XW = XW + (x(25) * 2.80541800d+01) ! C2H4
    XW = XW + (x(26) * 2.90621500d+01) ! C2H5
    XW = XW + (x(27) * 3.00701200d+01) ! C2H6
    XW = XW + (x(28) * 4.10296700d+01) ! HCCO
    XW = XW + (x(29) * 4.20376400d+01) ! CH2CO
    XW = XW + (x(30) * 4.20376400d+01) ! HCCOH
    XW = XW + (x(31) * 1.40067000d+01) ! N
    XW = XW + (x(32) * 1.50146700d+01) ! NH
    XW = XW + (x(33) * 1.60226400d+01) ! NH2
    XW = XW + (x(34) * 1.70306100d+01) ! NH3
    XW = XW + (x(35) * 2.90213700d+01) ! NNH
    XW = XW + (x(36) * 3.00061000d+01) ! NO
    XW = XW + (x(37) * 4.60055000d+01) ! NO2
    XW = XW + (x(38) * 4.40128000d+01) ! N2O
    XW = XW + (x(39) * 3.10140700d+01) ! HNO
    XW = XW + (x(40) * 2.60178500d+01) ! CN
    XW = XW + (x(41) * 2.70258200d+01) ! HCN
    XW = XW + (x(42) * 2.80337900d+01) ! H2CN
    XW = XW + (x(43) * 4.10325200d+01) ! HCNN
    XW = XW + (x(44) * 4.30252200d+01) ! HCNO
    XW = XW + (x(45) * 4.30252200d+01) ! HOCN
    XW = XW + (x(46) * 4.30252200d+01) ! HNCO
    XW = XW + (x(47) * 4.20172500d+01) ! NCO
    XW = XW + (x(48) * 2.80134000d+01) ! N2
    XW = XW + (x(49) * 3.99480000d+01) ! AR
    XW = XW + (x(50) * 4.30892400d+01) ! C3H7
    XW = XW + (x(51) * 4.40972100d+01) ! C3H8
    XW = XW + (x(52) * 4.30456100d+01) ! CH2CHO
    XW = XW + (x(53) * 4.40535800d+01) ! CH3CHO

    ! Now compute conversion
    XWinv = 1.d0 / XW
    y(1) = x(1) * 2.01594000d+00 * XWinv 
    y(2) = x(2) * 1.00797000d+00 * XWinv 
    y(3) = x(3) * 1.59994000d+01 * XWinv 
    y(4) = x(4) * 3.19988000d+01 * XWinv 
    y(5) = x(5) * 1.70073700d+01 * XWinv 
    y(6) = x(6) * 1.80153400d+01 * XWinv 
    y(7) = x(7) * 3.30067700d+01 * XWinv 
    y(8) = x(8) * 3.40147400d+01 * XWinv 
    y(9) = x(9) * 1.20111500d+01 * XWinv 
    y(10) = x(10) * 1.30191200d+01 * XWinv 
    y(11) = x(11) * 1.40270900d+01 * XWinv 
    y(12) = x(12) * 1.40270900d+01 * XWinv 
    y(13) = x(13) * 1.50350600d+01 * XWinv 
    y(14) = x(14) * 1.60430300d+01 * XWinv 
    y(15) = x(15) * 2.80105500d+01 * XWinv 
    y(16) = x(16) * 4.40099500d+01 * XWinv 
    y(17) = x(17) * 2.90185200d+01 * XWinv 
    y(18) = x(18) * 3.00264900d+01 * XWinv 
    y(19) = x(19) * 3.10344600d+01 * XWinv 
    y(20) = x(20) * 3.10344600d+01 * XWinv 
    y(21) = x(21) * 3.20424300d+01 * XWinv 
    y(22) = x(22) * 2.50302700d+01 * XWinv 
    y(23) = x(23) * 2.60382400d+01 * XWinv 
    y(24) = x(24) * 2.70462100d+01 * XWinv 
    y(25) = x(25) * 2.80541800d+01 * XWinv 
    y(26) = x(26) * 2.90621500d+01 * XWinv 
    y(27) = x(27) * 3.00701200d+01 * XWinv 
    y(28) = x(28) * 4.10296700d+01 * XWinv 
    y(29) = x(29) * 4.20376400d+01 * XWinv 
    y(30) = x(30) * 4.20376400d+01 * XWinv 
    y(31) = x(31) * 1.40067000d+01 * XWinv 
    y(32) = x(32) * 1.50146700d+01 * XWinv 
    y(33) = x(33) * 1.60226400d+01 * XWinv 
    y(34) = x(34) * 1.70306100d+01 * XWinv 
    y(35) = x(35) * 2.90213700d+01 * XWinv 
    y(36) = x(36) * 3.00061000d+01 * XWinv 
    y(37) = x(37) * 4.60055000d+01 * XWinv 
    y(38) = x(38) * 4.40128000d+01 * XWinv 
    y(39) = x(39) * 3.10140700d+01 * XWinv 
    y(40) = x(40) * 2.60178500d+01 * XWinv 
    y(41) = x(41) * 2.70258200d+01 * XWinv 
    y(42) = x(42) * 2.80337900d+01 * XWinv 
    y(43) = x(43) * 4.10325200d+01 * XWinv 
    y(44) = x(44) * 4.30252200d+01 * XWinv 
    y(45) = x(45) * 4.30252200d+01 * XWinv 
    y(46) = x(46) * 4.30252200d+01 * XWinv 
    y(47) = x(47) * 4.20172500d+01 * XWinv 
    y(48) = x(48) * 2.80134000d+01 * XWinv 
    y(49) = x(49) * 3.99480000d+01 * XWinv 
    y(50) = x(50) * 4.30892400d+01 * XWinv 
    y(51) = x(51) * 4.40972100d+01 * XWinv 
    y(52) = x(52) * 4.30456100d+01 * XWinv 
    y(53) = x(53) * 4.40535800d+01 * XWinv 

end subroutine

! Returns the specific heats at constant volume
! in mass units (Eq. 29)
subroutine ckcvms(T, iwrk, rwrk, cvms)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: cvms(53)

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
    cvms(8) = cvms(8) * 2.444384405113783d+06 !H2O2
    cvms(9) = cvms(9) * 6.922326338443862d+06 !C
    cvms(10) = cvms(10) * 6.386384025955671d+06 !CH
    cvms(11) = cvms(11) * 5.927466067445207d+06 !CH2
    cvms(12) = cvms(12) * 5.927466067445207d+06 !CH2(S)
    cvms(13) = cvms(13) * 5.530081023953346d+06 !CH3
    cvms(14) = cvms(14) * 5.182630712527496d+06 !CH4
    cvms(15) = cvms(15) * 2.968349425484326d+06 !CO
    cvms(16) = cvms(16) * 1.889234139098090d+06 !CO2
    cvms(17) = cvms(17) * 2.865242610581105d+06 !HCO
    cvms(18) = cvms(18) * 2.769058254894261d+06 !CH2O
    cvms(19) = cvms(19) * 2.679121853578248d+06 !CH2OH
    cvms(20) = cvms(20) * 2.679121853578248d+06 !CH3O
    cvms(21) = cvms(21) * 2.594843774332970d+06 !CH3OH
    cvms(22) = cvms(22) * 3.321781986370902d+06 !C2H
    cvms(23) = cvms(23) * 3.193192012977835d+06 !C2H2
    cvms(24) = cvms(24) * 3.074186734481467d+06 !C2H3
    cvms(25) = cvms(25) * 2.963733033722604d+06 !C2H4
    cvms(26) = cvms(26) * 2.860941121011349d+06 !C2H5
    cvms(27) = cvms(27) * 2.765040511976673d+06 !C2H6
    cvms(28) = cvms(28) * 2.026462801187531d+06 !HCCO
    cvms(29) = cvms(29) * 1.977872687429646d+06 !CH2CO
    cvms(30) = cvms(30) * 1.977872687429646d+06 !HCCOH
    cvms(31) = cvms(31) * 5.936094868884177d+06 !N
    cvms(32) = cvms(32) * 5.537590902763763d+06 !NH
    cvms(33) = cvms(33) * 5.189225995216768d+06 !NH2
    cvms(34) = cvms(34) * 4.882097587813943d+06 !NH3
    cvms(35) = cvms(35) * 2.864961233739138d+06 !NNH
    cvms(36) = cvms(36) * 2.770939908885194d+06 !NO
    cvms(37) = cvms(37) * 1.807286085359359d+06 !NO2
    cvms(38) = cvms(38) * 1.889111803838883d+06 !N2O
    cvms(39) = cvms(39) * 2.680883224936295d+06 !HNO
    cvms(40) = cvms(40) * 3.195694494356758d+06 !CN
    cvms(41) = cvms(41) * 3.076506096762281d+06 !HCN
    cvms(42) = cvms(42) * 2.965888665071686d+06 !H2CN
    cvms(43) = cvms(43) * 2.026322048950442d+06 !HCNN
    cvms(44) = cvms(44) * 1.932473558531484d+06 !HCNO
    cvms(45) = cvms(45) * 1.932473558531484d+06 !HOCN
    cvms(46) = cvms(46) * 1.932473558531484d+06 !HNCO
    cvms(47) = cvms(47) * 1.978832503317090d+06 !NCO
    cvms(48) = cvms(48) * 2.968047434442088d+06 !N2
    cvms(49) = cvms(49) * 2.081333233203164d+06 !AR
    cvms(50) = cvms(50) * 1.929602378691293d+06 !C3H7
    cvms(51) = cvms(51) * 1.885495703696447d+06 !C3H8
    cvms(52) = cvms(52) * 1.931558177477332d+06 !CH2CHO
    cvms(53) = cvms(53) * 1.887363070152301d+06 !CH3CHO

end subroutine

! Returns the specific heats at constant pressure
! in mass units (Eq. 26)
subroutine ckcpms(T, iwrk, rwrk, cpms)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: cpms(53)

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
    cpms(8) = cpms(8) * 2.444384405113783d+06 ! H2O2
    cpms(9) = cpms(9) * 6.922326338443862d+06 ! C
    cpms(10) = cpms(10) * 6.386384025955671d+06 ! CH
    cpms(11) = cpms(11) * 5.927466067445207d+06 ! CH2
    cpms(12) = cpms(12) * 5.927466067445207d+06 ! CH2(S)
    cpms(13) = cpms(13) * 5.530081023953346d+06 ! CH3
    cpms(14) = cpms(14) * 5.182630712527496d+06 ! CH4
    cpms(15) = cpms(15) * 2.968349425484326d+06 ! CO
    cpms(16) = cpms(16) * 1.889234139098090d+06 ! CO2
    cpms(17) = cpms(17) * 2.865242610581105d+06 ! HCO
    cpms(18) = cpms(18) * 2.769058254894261d+06 ! CH2O
    cpms(19) = cpms(19) * 2.679121853578248d+06 ! CH2OH
    cpms(20) = cpms(20) * 2.679121853578248d+06 ! CH3O
    cpms(21) = cpms(21) * 2.594843774332970d+06 ! CH3OH
    cpms(22) = cpms(22) * 3.321781986370902d+06 ! C2H
    cpms(23) = cpms(23) * 3.193192012977835d+06 ! C2H2
    cpms(24) = cpms(24) * 3.074186734481467d+06 ! C2H3
    cpms(25) = cpms(25) * 2.963733033722604d+06 ! C2H4
    cpms(26) = cpms(26) * 2.860941121011349d+06 ! C2H5
    cpms(27) = cpms(27) * 2.765040511976673d+06 ! C2H6
    cpms(28) = cpms(28) * 2.026462801187531d+06 ! HCCO
    cpms(29) = cpms(29) * 1.977872687429646d+06 ! CH2CO
    cpms(30) = cpms(30) * 1.977872687429646d+06 ! HCCOH
    cpms(31) = cpms(31) * 5.936094868884177d+06 ! N
    cpms(32) = cpms(32) * 5.537590902763763d+06 ! NH
    cpms(33) = cpms(33) * 5.189225995216768d+06 ! NH2
    cpms(34) = cpms(34) * 4.882097587813943d+06 ! NH3
    cpms(35) = cpms(35) * 2.864961233739138d+06 ! NNH
    cpms(36) = cpms(36) * 2.770939908885194d+06 ! NO
    cpms(37) = cpms(37) * 1.807286085359359d+06 ! NO2
    cpms(38) = cpms(38) * 1.889111803838883d+06 ! N2O
    cpms(39) = cpms(39) * 2.680883224936295d+06 ! HNO
    cpms(40) = cpms(40) * 3.195694494356758d+06 ! CN
    cpms(41) = cpms(41) * 3.076506096762281d+06 ! HCN
    cpms(42) = cpms(42) * 2.965888665071686d+06 ! H2CN
    cpms(43) = cpms(43) * 2.026322048950442d+06 ! HCNN
    cpms(44) = cpms(44) * 1.932473558531484d+06 ! HCNO
    cpms(45) = cpms(45) * 1.932473558531484d+06 ! HOCN
    cpms(46) = cpms(46) * 1.932473558531484d+06 ! HNCO
    cpms(47) = cpms(47) * 1.978832503317090d+06 ! NCO
    cpms(48) = cpms(48) * 2.968047434442088d+06 ! N2
    cpms(49) = cpms(49) * 2.081333233203164d+06 ! AR
    cpms(50) = cpms(50) * 1.929602378691293d+06 ! C3H7
    cpms(51) = cpms(51) * 1.885495703696447d+06 ! C3H8
    cpms(52) = cpms(52) * 1.931558177477332d+06 ! CH2CHO
    cpms(53) = cpms(53) * 1.887363070152301d+06 ! CH3CHO

end subroutine

! Returns internal energy in mass units (Eq 30.)
subroutine ckums(T, iwrk, rwrk, ums)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: ums(53)

    double precision :: tT, tc(5)
    double precision :: RT
    integer :: i

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesInternalEnergy(ums, tc)

    do i=1, 53
        ums(i) = ums(i) * (RT * imw(i))
    end do

end subroutine

! Returns enthalpy in mass units (Eq 27.)
subroutine ckhms(T, iwrk, rwrk, hms)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: hms(53)

    double precision :: tT, RT
    double precision :: tc(5), h(53)
    integer :: i

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesEnthalpy(hms, tc)

    do i=1, 53
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
    double precision, intent(inout) :: hms(np,53)

    double precision :: tc(5), h(53)
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
        hms(i, 22) = h(22)
        hms(i, 23) = h(23)
        hms(i, 24) = h(24)
        hms(i, 25) = h(25)
        hms(i, 26) = h(26)
        hms(i, 27) = h(27)
        hms(i, 28) = h(28)
        hms(i, 29) = h(29)
        hms(i, 30) = h(30)
        hms(i, 31) = h(31)
        hms(i, 32) = h(32)
        hms(i, 33) = h(33)
        hms(i, 34) = h(34)
        hms(i, 35) = h(35)
        hms(i, 36) = h(36)
        hms(i, 37) = h(37)
        hms(i, 38) = h(38)
        hms(i, 39) = h(39)
        hms(i, 40) = h(40)
        hms(i, 41) = h(41)
        hms(i, 42) = h(42)
        hms(i, 43) = h(43)
        hms(i, 44) = h(44)
        hms(i, 45) = h(45)
        hms(i, 46) = h(46)
        hms(i, 47) = h(47)
        hms(i, 48) = h(48)
        hms(i, 49) = h(49)
        hms(i, 50) = h(50)
        hms(i, 51) = h(51)
        hms(i, 52) = h(52)
        hms(i, 53) = h(53)
    end do

    do n=1, 53
        do i=1, np
            hms(i,n) = hms(i,n) * ( 8.31451000d+07 * T(i) * imw(n))
        end do
    end do

end subroutine

! Returns the mean specific heat at CP (Eq. 34)
subroutine ckcpbs(T, y, iwrk, rwrk, cpbs)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: cpbs

    double precision :: cpor(53)
    double precision :: tresult(53)
    double precision :: tT, tc(5)
    double precision :: res
    integer :: i

    res = 0.d0

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache

    call cp_R(cpor, tc)

    do i=1, 53
        tresult(i) = cpor(i) * y(i) * imw(i)
    end do
    do i=1, 53
        res = res + tresult(i)
    end do

    cpbs = res * 8.31451000d+07

end subroutine

! Returns the mean specific heat at CV (Eq. 36)
subroutine ckcvbs(T, y, iwrk, rwrk, cvbs)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: cvbs

    double precision :: cvor(53)
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
    res = res + (cvor(8) * y(8) * imw(8)) ! H2O2
    res = res + (cvor(9) * y(9) * imw(9)) ! C
    res = res + (cvor(10) * y(10) * imw(10)) ! CH
    res = res + (cvor(11) * y(11) * imw(11)) ! CH2
    res = res + (cvor(12) * y(12) * imw(12)) ! CH2(S)
    res = res + (cvor(13) * y(13) * imw(13)) ! CH3
    res = res + (cvor(14) * y(14) * imw(14)) ! CH4
    res = res + (cvor(15) * y(15) * imw(15)) ! CO
    res = res + (cvor(16) * y(16) * imw(16)) ! CO2
    res = res + (cvor(17) * y(17) * imw(17)) ! HCO
    res = res + (cvor(18) * y(18) * imw(18)) ! CH2O
    res = res + (cvor(19) * y(19) * imw(19)) ! CH2OH
    res = res + (cvor(20) * y(20) * imw(20)) ! CH3O
    res = res + (cvor(21) * y(21) * imw(21)) ! CH3OH
    res = res + (cvor(22) * y(22) * imw(22)) ! C2H
    res = res + (cvor(23) * y(23) * imw(23)) ! C2H2
    res = res + (cvor(24) * y(24) * imw(24)) ! C2H3
    res = res + (cvor(25) * y(25) * imw(25)) ! C2H4
    res = res + (cvor(26) * y(26) * imw(26)) ! C2H5
    res = res + (cvor(27) * y(27) * imw(27)) ! C2H6
    res = res + (cvor(28) * y(28) * imw(28)) ! HCCO
    res = res + (cvor(29) * y(29) * imw(29)) ! CH2CO
    res = res + (cvor(30) * y(30) * imw(30)) ! HCCOH
    res = res + (cvor(31) * y(31) * imw(31)) ! N
    res = res + (cvor(32) * y(32) * imw(32)) ! NH
    res = res + (cvor(33) * y(33) * imw(33)) ! NH2
    res = res + (cvor(34) * y(34) * imw(34)) ! NH3
    res = res + (cvor(35) * y(35) * imw(35)) ! NNH
    res = res + (cvor(36) * y(36) * imw(36)) ! NO
    res = res + (cvor(37) * y(37) * imw(37)) ! NO2
    res = res + (cvor(38) * y(38) * imw(38)) ! N2O
    res = res + (cvor(39) * y(39) * imw(39)) ! HNO
    res = res + (cvor(40) * y(40) * imw(40)) ! CN
    res = res + (cvor(41) * y(41) * imw(41)) ! HCN
    res = res + (cvor(42) * y(42) * imw(42)) ! H2CN
    res = res + (cvor(43) * y(43) * imw(43)) ! HCNN
    res = res + (cvor(44) * y(44) * imw(44)) ! HCNO
    res = res + (cvor(45) * y(45) * imw(45)) ! HOCN
    res = res + (cvor(46) * y(46) * imw(46)) ! HNCO
    res = res + (cvor(47) * y(47) * imw(47)) ! NCO
    res = res + (cvor(48) * y(48) * imw(48)) ! N2
    res = res + (cvor(49) * y(49) * imw(49)) ! AR
    res = res + (cvor(50) * y(50) * imw(50)) ! C3H7
    res = res + (cvor(51) * y(51) * imw(51)) ! C3H8
    res = res + (cvor(52) * y(52) * imw(52)) ! CH2CHO
    res = res + (cvor(53) * y(53) * imw(53)) ! CH3CHO

    cvbs = res * 8.31451000d+07

end subroutine

! get mean internal energy in mass units
subroutine ckubms(T, y, iwrk, rwrk, ubms)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(in) :: y(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: ubms

    double precision :: ums(53) ! temporary energy array
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
    res = res + (y(8) * ums(8) * imw(8)) ! H2O2
    res = res + (y(9) * ums(9) * imw(9)) ! C
    res = res + (y(10) * ums(10) * imw(10)) ! CH
    res = res + (y(11) * ums(11) * imw(11)) ! CH2
    res = res + (y(12) * ums(12) * imw(12)) ! CH2(S)
    res = res + (y(13) * ums(13) * imw(13)) ! CH3
    res = res + (y(14) * ums(14) * imw(14)) ! CH4
    res = res + (y(15) * ums(15) * imw(15)) ! CO
    res = res + (y(16) * ums(16) * imw(16)) ! CO2
    res = res + (y(17) * ums(17) * imw(17)) ! HCO
    res = res + (y(18) * ums(18) * imw(18)) ! CH2O
    res = res + (y(19) * ums(19) * imw(19)) ! CH2OH
    res = res + (y(20) * ums(20) * imw(20)) ! CH3O
    res = res + (y(21) * ums(21) * imw(21)) ! CH3OH
    res = res + (y(22) * ums(22) * imw(22)) ! C2H
    res = res + (y(23) * ums(23) * imw(23)) ! C2H2
    res = res + (y(24) * ums(24) * imw(24)) ! C2H3
    res = res + (y(25) * ums(25) * imw(25)) ! C2H4
    res = res + (y(26) * ums(26) * imw(26)) ! C2H5
    res = res + (y(27) * ums(27) * imw(27)) ! C2H6
    res = res + (y(28) * ums(28) * imw(28)) ! HCCO
    res = res + (y(29) * ums(29) * imw(29)) ! CH2CO
    res = res + (y(30) * ums(30) * imw(30)) ! HCCOH
    res = res + (y(31) * ums(31) * imw(31)) ! N
    res = res + (y(32) * ums(32) * imw(32)) ! NH
    res = res + (y(33) * ums(33) * imw(33)) ! NH2
    res = res + (y(34) * ums(34) * imw(34)) ! NH3
    res = res + (y(35) * ums(35) * imw(35)) ! NNH
    res = res + (y(36) * ums(36) * imw(36)) ! NO
    res = res + (y(37) * ums(37) * imw(37)) ! NO2
    res = res + (y(38) * ums(38) * imw(38)) ! N2O
    res = res + (y(39) * ums(39) * imw(39)) ! HNO
    res = res + (y(40) * ums(40) * imw(40)) ! CN
    res = res + (y(41) * ums(41) * imw(41)) ! HCN
    res = res + (y(42) * ums(42) * imw(42)) ! H2CN
    res = res + (y(43) * ums(43) * imw(43)) ! HCNN
    res = res + (y(44) * ums(44) * imw(44)) ! HCNO
    res = res + (y(45) * ums(45) * imw(45)) ! HOCN
    res = res + (y(46) * ums(46) * imw(46)) ! HNCO
    res = res + (y(47) * ums(47) * imw(47)) ! NCO
    res = res + (y(48) * ums(48) * imw(48)) ! N2
    res = res + (y(49) * ums(49) * imw(49)) ! AR
    res = res + (y(50) * ums(50) * imw(50)) ! C3H7
    res = res + (y(51) * ums(51) * imw(51)) ! C3H8
    res = res + (y(52) * ums(52) * imw(52)) ! CH2CHO
    res = res + (y(53) * ums(53) * imw(53)) ! CH3CHO

    ubms = res * RT

end subroutine

! compute the production rate for each species
subroutine ckwc(T, C, iwrk, rwrk, wdot)

    implicit none

    double precision, intent(in) :: T
    double precision, intent(inout) :: C(53)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: wdot(53)

    integer :: id ! loop counter

    ! convert to SI
    do id = 1, 53
        C(id) = C(id) * 1.0d6
    end do

    ! convert to chemkin units
    call productionRate(wdot, C, T)

    ! convert to chemkin units
    do id=1, 53
        C(id) = C(id) * 1.0d-6
        wdot(id) = wdot(id) * 1.0d-6
    end do

end subroutine

! compute the production rate for each species
subroutine productionRate(wdot, sc, T)

    implicit none

    double precision, intent(inout) :: wdot(53)
    double precision, intent(in) :: sc(53)
    double precision, intent(in) :: T

    double precision :: tc(5)
    double precision :: invT
    double precision :: qdot, q_f(325), q_r(325)
    integer :: i

    tc = (/ log(T), T, T*T, T*T*T, T*T*T*T /)
    invT = 1.d0 / tc(2)

    if (T /= T_save) then
        T_save = T
        call comp_k_f(tc,invT,k_f_save)
        call comp_Kc(tc,invT,Kc_save)
    end if

    call comp_qfqr(q_f, q_r, sc, tc, invT)

    do i=1, 53
        wdot(i) = 0.d0
    end do

    qdot = q_f(1)-q_r(1)
    wdot(2) = wdot(2) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(2)-q_r(2)
    wdot(2) = wdot(2) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot

    qdot = q_f(3)-q_r(3)
    wdot(2) = wdot(2) - qdot
    wdot(17) = wdot(17) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(4)-q_r(4)
    wdot(2) = wdot(2) - qdot
    wdot(18) = wdot(18) - qdot
    wdot(19) = wdot(19) + qdot

    qdot = q_f(5)-q_r(5)
    wdot(2) = wdot(2) - qdot
    wdot(18) = wdot(18) - qdot
    wdot(20) = wdot(20) + qdot

    qdot = q_f(6)-q_r(6)
    wdot(2) = wdot(2) - qdot
    wdot(19) = wdot(19) - qdot
    wdot(21) = wdot(21) + qdot

    qdot = q_f(7)-q_r(7)
    wdot(2) = wdot(2) - qdot
    wdot(20) = wdot(20) - qdot
    wdot(21) = wdot(21) + qdot

    qdot = q_f(8)-q_r(8)
    wdot(2) = wdot(2) - qdot
    wdot(22) = wdot(22) - qdot
    wdot(23) = wdot(23) + qdot

    qdot = q_f(9)-q_r(9)
    wdot(2) = wdot(2) - qdot
    wdot(23) = wdot(23) - qdot
    wdot(24) = wdot(24) + qdot

    qdot = q_f(10)-q_r(10)
    wdot(2) = wdot(2) - qdot
    wdot(24) = wdot(24) - qdot
    wdot(25) = wdot(25) + qdot

    qdot = q_f(11)-q_r(11)
    wdot(2) = wdot(2) - qdot
    wdot(25) = wdot(25) - qdot
    wdot(26) = wdot(26) + qdot

    qdot = q_f(12)-q_r(12)
    wdot(2) = wdot(2) - qdot
    wdot(26) = wdot(26) - qdot
    wdot(27) = wdot(27) + qdot

    qdot = q_f(13)-q_r(13)
    wdot(1) = wdot(1) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(14)-q_r(14)
    wdot(5) = wdot(5) - (2 * qdot)
    wdot(8) = wdot(8) + qdot

    qdot = q_f(15)-q_r(15)
    wdot(5) = wdot(5) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(21) = wdot(21) + qdot

    qdot = q_f(16)-q_r(16)
    wdot(10) = wdot(10) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(28) = wdot(28) + qdot

    qdot = q_f(17)-q_r(17)
    wdot(11) = wdot(11) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(29) = wdot(29) + qdot

    qdot = q_f(18)-q_r(18)
    wdot(6) = wdot(6) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(21) = wdot(21) + qdot

    qdot = q_f(19)-q_r(19)
    wdot(13) = wdot(13) - (2 * qdot)
    wdot(27) = wdot(27) + qdot

    qdot = q_f(20)-q_r(20)
    wdot(1) = wdot(1) + qdot
    wdot(23) = wdot(23) + qdot
    wdot(25) = wdot(25) - qdot

    qdot = q_f(21)-q_r(21)
    wdot(10) = wdot(10) - qdot
    wdot(43) = wdot(43) + qdot
    wdot(48) = wdot(48) - qdot

    qdot = q_f(22)-q_r(22)
    wdot(1) = wdot(1) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(23)-q_r(23)
    wdot(2) = wdot(2) - qdot
    wdot(29) = wdot(29) - qdot
    wdot(52) = wdot(52) + qdot

    qdot = q_f(24)-q_r(24)
    wdot(13) = wdot(13) - qdot
    wdot(26) = wdot(26) - qdot
    wdot(51) = wdot(51) + qdot

    qdot = q_f(25)-q_r(25)
    wdot(13) = wdot(13) - qdot
    wdot(25) = wdot(25) - qdot
    wdot(50) = wdot(50) + qdot

    qdot = q_f(26)-q_r(26)
    wdot(2) = wdot(2) - qdot
    wdot(50) = wdot(50) - qdot
    wdot(51) = wdot(51) + qdot

    qdot = q_f(27)-q_r(27)
    wdot(3) = wdot(3) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(28)-q_r(28)
    wdot(3) = wdot(3) + qdot
    wdot(38) = wdot(38) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(29)-q_r(29)
    wdot(2) = wdot(2) - qdot
    wdot(41) = wdot(41) - qdot
    wdot(42) = wdot(42) + qdot

    qdot = q_f(30)-q_r(30)
    wdot(3) = wdot(3) - (2 * qdot)
    wdot(4) = wdot(4) + qdot

    qdot = q_f(31)-q_r(31)
    wdot(2) = wdot(2) - qdot
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot

    qdot = q_f(32)-q_r(32)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot

    qdot = q_f(33)-q_r(33)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - (2 * qdot)

    qdot = q_f(34)-q_r(34)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(35)-q_r(35)
    wdot(2) = wdot(2) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(36)-q_r(36)
    wdot(3) = wdot(3) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(37) = wdot(37) + qdot

    qdot = q_f(37)-q_r(37)
    wdot(2) = wdot(2) + qdot
    wdot(35) = wdot(35) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(38)-q_r(38)
    wdot(2) = wdot(2) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(39) = wdot(39) + qdot

    qdot = q_f(39)-q_r(39)
    wdot(15) = wdot(15) + qdot
    wdot(31) = wdot(31) + qdot
    wdot(47) = wdot(47) - qdot

    qdot = q_f(40)-q_r(40)
    wdot(2) = wdot(2) + qdot
    wdot(40) = wdot(40) + qdot
    wdot(41) = wdot(41) - qdot

    qdot = q_f(41)-q_r(41)
    wdot(15) = wdot(15) + qdot
    wdot(32) = wdot(32) + qdot
    wdot(46) = wdot(46) - qdot

    qdot = q_f(42)-q_r(42)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot

    qdot = q_f(43)-q_r(43)
    wdot(3) = wdot(3) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(44)-q_r(44)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(45)-q_r(45)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(46)-q_r(46)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(17) = wdot(17) + qdot

    qdot = q_f(47)-q_r(47)
    wdot(1) = wdot(1) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(48)-q_r(48)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(17) = wdot(17) + qdot

    qdot = q_f(49)-q_r(49)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(50)-q_r(50)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(51)-q_r(51)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(52)-q_r(52)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(16) = wdot(16) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(53)-q_r(53)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) - qdot

    qdot = q_f(54)-q_r(54)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(55)-q_r(55)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(20) = wdot(20) - qdot

    qdot = q_f(56)-q_r(56)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(19) = wdot(19) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(57)-q_r(57)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(20) = wdot(20) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(58)-q_r(58)
    wdot(3) = wdot(3) - qdot
    wdot(10) = wdot(10) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(22) = wdot(22) - qdot

    qdot = q_f(59)-q_r(59)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(23) = wdot(23) - qdot
    wdot(28) = wdot(28) + qdot

    qdot = q_f(60)-q_r(60)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(22) = wdot(22) + qdot
    wdot(23) = wdot(23) - qdot

    qdot = q_f(61)-q_r(61)
    wdot(3) = wdot(3) - qdot
    wdot(11) = wdot(11) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(23) = wdot(23) - qdot

    qdot = q_f(62)-q_r(62)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(24) = wdot(24) - qdot
    wdot(29) = wdot(29) + qdot

    qdot = q_f(63)-q_r(63)
    wdot(3) = wdot(3) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(25) = wdot(25) - qdot

    qdot = q_f(64)-q_r(64)
    wdot(3) = wdot(3) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(26) = wdot(26) - qdot

    qdot = q_f(65)-q_r(65)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(27) = wdot(27) - qdot

    qdot = q_f(66)-q_r(66)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(15) = wdot(15) + (2 * qdot)
    wdot(28) = wdot(28) - qdot

    qdot = q_f(67)-q_r(67)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(28) = wdot(28) + qdot
    wdot(29) = wdot(29) - qdot

    qdot = q_f(68)-q_r(68)
    wdot(3) = wdot(3) - qdot
    wdot(11) = wdot(11) + qdot
    wdot(16) = wdot(16) + qdot
    wdot(29) = wdot(29) - qdot

    qdot = q_f(69)-q_r(69)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(70)-q_r(70)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) - qdot

    qdot = q_f(71)-q_r(71)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - (2 * qdot)
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) + qdot

    qdot = q_f(72)-q_r(72)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) + qdot

    qdot = q_f(73)-q_r(73)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(48) = wdot(48) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(74)-q_r(74)
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(49) = wdot(49) - qdot
    wdot(49) = wdot(49) + qdot

    qdot = q_f(75)-q_r(75)
    wdot(2) = wdot(2) - qdot
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot

    qdot = q_f(76)-q_r(76)
    wdot(1) = wdot(1) - qdot
    wdot(1) = wdot(1) + (2 * qdot)
    wdot(2) = wdot(2) - (2 * qdot)

    qdot = q_f(77)-q_r(77)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - (2 * qdot)
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(78)-q_r(78)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - (2 * qdot)
    wdot(16) = wdot(16) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(79)-q_r(79)
    wdot(2) = wdot(2) - qdot
    wdot(3) = wdot(3) + qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(80)-q_r(80)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(81)-q_r(81)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + (2 * qdot)
    wdot(7) = wdot(7) - qdot

    qdot = q_f(82)-q_r(82)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(83)-q_r(83)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(6) = wdot(6) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(84)-q_r(84)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(9) = wdot(9) + qdot
    wdot(10) = wdot(10) - qdot

    qdot = q_f(85)-q_r(85)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(10) = wdot(10) + qdot
    wdot(12) = wdot(12) - qdot

    qdot = q_f(86)-q_r(86)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(87)-q_r(87)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(88)-q_r(88)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) - qdot

    qdot = q_f(89)-q_r(89)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(90)-q_r(90)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(91)-q_r(91)
    wdot(2) = wdot(2) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(92)-q_r(92)
    wdot(2) = wdot(2) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(19) = wdot(19) + qdot
    wdot(20) = wdot(20) - qdot

    qdot = q_f(93)-q_r(93)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(18) = wdot(18) + qdot
    wdot(20) = wdot(20) - qdot

    qdot = q_f(94)-q_r(94)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(20) = wdot(20) - qdot

    qdot = q_f(95)-q_r(95)
    wdot(2) = wdot(2) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(20) = wdot(20) - qdot

    qdot = q_f(96)-q_r(96)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(19) = wdot(19) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(97)-q_r(97)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(20) = wdot(20) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(98)-q_r(98)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(23) = wdot(23) + qdot
    wdot(24) = wdot(24) - qdot

    qdot = q_f(99)-q_r(99)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(24) = wdot(24) + qdot
    wdot(25) = wdot(25) - qdot

    qdot = q_f(100)-q_r(100)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(25) = wdot(25) + qdot
    wdot(26) = wdot(26) - qdot

    qdot = q_f(101)-q_r(101)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(26) = wdot(26) + qdot
    wdot(27) = wdot(27) - qdot

    qdot = q_f(102)-q_r(102)
    wdot(2) = wdot(2) - qdot
    wdot(12) = wdot(12) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(28) = wdot(28) - qdot

    qdot = q_f(103)-q_r(103)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(28) = wdot(28) + qdot
    wdot(29) = wdot(29) - qdot

    qdot = q_f(104)-q_r(104)
    wdot(2) = wdot(2) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(29) = wdot(29) - qdot

    qdot = q_f(105)-q_r(105)
    wdot(2) = wdot(2) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(29) = wdot(29) + qdot
    wdot(30) = wdot(30) - qdot

    qdot = q_f(106)-q_r(106)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot

    qdot = q_f(107)-q_r(107)
    wdot(3) = wdot(3) + qdot
    wdot(5) = wdot(5) - (2 * qdot)
    wdot(6) = wdot(6) + qdot

    qdot = q_f(108)-q_r(108)
    wdot(4) = wdot(4) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(109)-q_r(109)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(110)-q_r(110)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot

    qdot = q_f(111)-q_r(111)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(9) = wdot(9) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(112)-q_r(112)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(17) = wdot(17) + qdot

    qdot = q_f(113)-q_r(113)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(114)-q_r(114)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(10) = wdot(10) + qdot
    wdot(11) = wdot(11) - qdot

    qdot = q_f(115)-q_r(115)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(116)-q_r(116)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(11) = wdot(11) + qdot
    wdot(13) = wdot(13) - qdot

    qdot = q_f(117)-q_r(117)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(12) = wdot(12) + qdot
    wdot(13) = wdot(13) - qdot

    qdot = q_f(118)-q_r(118)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(14) = wdot(14) - qdot

    qdot = q_f(119)-q_r(119)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(120)-q_r(120)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(121)-q_r(121)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) - qdot

    qdot = q_f(122)-q_r(122)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(123)-q_r(123)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(20) = wdot(20) - qdot

    qdot = q_f(124)-q_r(124)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(19) = wdot(19) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(125)-q_r(125)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(20) = wdot(20) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(126)-q_r(126)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(22) = wdot(22) - qdot
    wdot(28) = wdot(28) + qdot

    qdot = q_f(127)-q_r(127)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(23) = wdot(23) - qdot
    wdot(29) = wdot(29) + qdot

    qdot = q_f(128)-q_r(128)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(23) = wdot(23) - qdot
    wdot(30) = wdot(30) + qdot

    qdot = q_f(129)-q_r(129)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(22) = wdot(22) + qdot
    wdot(23) = wdot(23) - qdot

    qdot = q_f(130)-q_r(130)
    wdot(5) = wdot(5) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(23) = wdot(23) - qdot

    qdot = q_f(131)-q_r(131)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(23) = wdot(23) + qdot
    wdot(24) = wdot(24) - qdot

    qdot = q_f(132)-q_r(132)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(24) = wdot(24) + qdot
    wdot(25) = wdot(25) - qdot

    qdot = q_f(133)-q_r(133)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(27) = wdot(27) - qdot

    qdot = q_f(134)-q_r(134)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(28) = wdot(28) + qdot
    wdot(29) = wdot(29) - qdot

    qdot = q_f(135)-q_r(135)
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) - (2 * qdot)
    wdot(8) = wdot(8) + qdot

    qdot = q_f(136)-q_r(136)
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) - (2 * qdot)
    wdot(8) = wdot(8) + qdot

    qdot = q_f(137)-q_r(137)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(138)-q_r(138)
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot

    qdot = q_f(139)-q_r(139)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(20) = wdot(20) + qdot

    qdot = q_f(140)-q_r(140)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(141)-q_r(141)
    wdot(7) = wdot(7) - qdot
    wdot(8) = wdot(8) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) - qdot

    qdot = q_f(142)-q_r(142)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(9) = wdot(9) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(143)-q_r(143)
    wdot(2) = wdot(2) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(22) = wdot(22) + qdot

    qdot = q_f(144)-q_r(144)
    wdot(2) = wdot(2) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(23) = wdot(23) + qdot

    qdot = q_f(145)-q_r(145)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(17) = wdot(17) + qdot

    qdot = q_f(146)-q_r(146)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(11) = wdot(11) + qdot

    qdot = q_f(147)-q_r(147)
    wdot(2) = wdot(2) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(10) = wdot(10) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(148)-q_r(148)
    wdot(2) = wdot(2) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(23) = wdot(23) + qdot

    qdot = q_f(149)-q_r(149)
    wdot(2) = wdot(2) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(24) = wdot(24) + qdot

    qdot = q_f(150)-q_r(150)
    wdot(2) = wdot(2) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(14) = wdot(14) - qdot
    wdot(25) = wdot(25) + qdot

    qdot = q_f(151)-q_r(151)
    wdot(10) = wdot(10) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(16) = wdot(16) - qdot
    wdot(17) = wdot(17) + qdot

    qdot = q_f(152)-q_r(152)
    wdot(2) = wdot(2) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(18) = wdot(18) - qdot
    wdot(29) = wdot(29) + qdot

    qdot = q_f(153)-q_r(153)
    wdot(10) = wdot(10) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(23) = wdot(23) + qdot
    wdot(28) = wdot(28) - qdot

    qdot = q_f(154)-q_r(154)
    wdot(2) = wdot(2) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(11) = wdot(11) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(155)-q_r(155)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(11) = wdot(11) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(156)-q_r(156)
    wdot(1) = wdot(1) + qdot
    wdot(11) = wdot(11) - (2 * qdot)
    wdot(23) = wdot(23) + qdot

    qdot = q_f(157)-q_r(157)
    wdot(2) = wdot(2) + qdot
    wdot(11) = wdot(11) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(25) = wdot(25) + qdot

    qdot = q_f(158)-q_r(158)
    wdot(11) = wdot(11) - qdot
    wdot(13) = wdot(13) + (2 * qdot)
    wdot(14) = wdot(14) - qdot

    qdot = q_f(159)-q_r(159)
    wdot(11) = wdot(11) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(24) = wdot(24) + qdot
    wdot(28) = wdot(28) - qdot

    qdot = q_f(160)-q_r(160)
    wdot(11) = wdot(11) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(48) = wdot(48) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(161)-q_r(161)
    wdot(11) = wdot(11) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(49) = wdot(49) - qdot
    wdot(49) = wdot(49) + qdot

    qdot = q_f(162)-q_r(162)
    wdot(2) = wdot(2) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(163)-q_r(163)
    wdot(4) = wdot(4) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(164)-q_r(164)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) + qdot

    qdot = q_f(165)-q_r(165)
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(11) = wdot(11) + qdot
    wdot(12) = wdot(12) - qdot

    qdot = q_f(166)-q_r(166)
    wdot(2) = wdot(2) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(25) = wdot(25) + qdot

    qdot = q_f(167)-q_r(167)
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) + (2 * qdot)
    wdot(14) = wdot(14) - qdot

    qdot = q_f(168)-q_r(168)
    wdot(11) = wdot(11) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(15) = wdot(15) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(169)-q_r(169)
    wdot(11) = wdot(11) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(16) = wdot(16) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(170)-q_r(170)
    wdot(12) = wdot(12) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(16) = wdot(16) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(171)-q_r(171)
    wdot(12) = wdot(12) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(27) = wdot(27) - qdot

    qdot = q_f(172)-q_r(172)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(20) = wdot(20) + qdot

    qdot = q_f(173)-q_r(173)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(13) = wdot(13) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(174)-q_r(174)
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot

    qdot = q_f(175)-q_r(175)
    wdot(2) = wdot(2) + qdot
    wdot(13) = wdot(13) - (2 * qdot)
    wdot(26) = wdot(26) + qdot

    qdot = q_f(176)-q_r(176)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(177)-q_r(177)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) - qdot

    qdot = q_f(178)-q_r(178)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(19) = wdot(19) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(179)-q_r(179)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(20) = wdot(20) + qdot
    wdot(21) = wdot(21) - qdot

    qdot = q_f(180)-q_r(180)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(24) = wdot(24) + qdot
    wdot(25) = wdot(25) - qdot

    qdot = q_f(181)-q_r(181)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(27) = wdot(27) - qdot

    qdot = q_f(182)-q_r(182)
    wdot(2) = wdot(2) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(183)-q_r(183)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) - qdot

    qdot = q_f(184)-q_r(184)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(19) = wdot(19) - qdot

    qdot = q_f(185)-q_r(185)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(20) = wdot(20) - qdot

    qdot = q_f(186)-q_r(186)
    wdot(4) = wdot(4) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(22) = wdot(22) - qdot

    qdot = q_f(187)-q_r(187)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(22) = wdot(22) - qdot
    wdot(23) = wdot(23) + qdot

    qdot = q_f(188)-q_r(188)
    wdot(4) = wdot(4) - qdot
    wdot(17) = wdot(17) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(24) = wdot(24) - qdot

    qdot = q_f(189)-q_r(189)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(25) = wdot(25) + qdot
    wdot(26) = wdot(26) - qdot

    qdot = q_f(190)-q_r(190)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(15) = wdot(15) + (2 * qdot)
    wdot(28) = wdot(28) - qdot

    qdot = q_f(191)-q_r(191)
    wdot(15) = wdot(15) + (2 * qdot)
    wdot(23) = wdot(23) + qdot
    wdot(28) = wdot(28) - (2 * qdot)

    qdot = q_f(192)-q_r(192)
    wdot(3) = wdot(3) + qdot
    wdot(31) = wdot(31) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(193)-q_r(193)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(31) = wdot(31) - qdot
    wdot(36) = wdot(36) + qdot

    qdot = q_f(194)-q_r(194)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(31) = wdot(31) - qdot
    wdot(36) = wdot(36) + qdot

    qdot = q_f(195)-q_r(195)
    wdot(3) = wdot(3) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(38) = wdot(38) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(196)-q_r(196)
    wdot(3) = wdot(3) - qdot
    wdot(36) = wdot(36) + (2 * qdot)
    wdot(38) = wdot(38) - qdot

    qdot = q_f(197)-q_r(197)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(38) = wdot(38) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(198)-q_r(198)
    wdot(5) = wdot(5) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(38) = wdot(38) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(199)-q_r(199)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(37) = wdot(37) + qdot

    qdot = q_f(200)-q_r(200)
    wdot(3) = wdot(3) - qdot
    wdot(4) = wdot(4) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(37) = wdot(37) - qdot

    qdot = q_f(201)-q_r(201)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(37) = wdot(37) - qdot

    qdot = q_f(202)-q_r(202)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(32) = wdot(32) - qdot
    wdot(36) = wdot(36) + qdot

    qdot = q_f(203)-q_r(203)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(31) = wdot(31) + qdot
    wdot(32) = wdot(32) - qdot

    qdot = q_f(204)-q_r(204)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(32) = wdot(32) - qdot
    wdot(39) = wdot(39) + qdot

    qdot = q_f(205)-q_r(205)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(31) = wdot(31) + qdot
    wdot(32) = wdot(32) - qdot

    qdot = q_f(206)-q_r(206)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(32) = wdot(32) - qdot
    wdot(39) = wdot(39) + qdot

    qdot = q_f(207)-q_r(207)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(32) = wdot(32) - qdot
    wdot(36) = wdot(36) + qdot

    qdot = q_f(208)-q_r(208)
    wdot(2) = wdot(2) + qdot
    wdot(31) = wdot(31) - qdot
    wdot(32) = wdot(32) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(209)-q_r(209)
    wdot(1) = wdot(1) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(32) = wdot(32) - qdot
    wdot(39) = wdot(39) + qdot

    qdot = q_f(210)-q_r(210)
    wdot(5) = wdot(5) + qdot
    wdot(32) = wdot(32) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(211)-q_r(211)
    wdot(2) = wdot(2) + qdot
    wdot(32) = wdot(32) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(38) = wdot(38) + qdot

    qdot = q_f(212)-q_r(212)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(32) = wdot(32) + qdot
    wdot(33) = wdot(33) - qdot

    qdot = q_f(213)-q_r(213)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(33) = wdot(33) - qdot
    wdot(39) = wdot(39) + qdot

    qdot = q_f(214)-q_r(214)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(32) = wdot(32) + qdot
    wdot(33) = wdot(33) - qdot

    qdot = q_f(215)-q_r(215)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(32) = wdot(32) + qdot
    wdot(33) = wdot(33) - qdot

    qdot = q_f(216)-q_r(216)
    wdot(2) = wdot(2) + qdot
    wdot(35) = wdot(35) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(217)-q_r(217)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(35) = wdot(35) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(218)-q_r(218)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(35) = wdot(35) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(219)-q_r(219)
    wdot(3) = wdot(3) - qdot
    wdot(32) = wdot(32) + qdot
    wdot(35) = wdot(35) - qdot
    wdot(36) = wdot(36) + qdot

    qdot = q_f(220)-q_r(220)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(35) = wdot(35) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(221)-q_r(221)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(35) = wdot(35) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(222)-q_r(222)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(35) = wdot(35) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(223)-q_r(223)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(39) = wdot(39) - qdot

    qdot = q_f(224)-q_r(224)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(36) = wdot(36) + qdot
    wdot(39) = wdot(39) - qdot

    qdot = q_f(225)-q_r(225)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(39) = wdot(39) - qdot

    qdot = q_f(226)-q_r(226)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(39) = wdot(39) - qdot

    qdot = q_f(227)-q_r(227)
    wdot(3) = wdot(3) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(31) = wdot(31) + qdot
    wdot(40) = wdot(40) - qdot

    qdot = q_f(228)-q_r(228)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(40) = wdot(40) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(229)-q_r(229)
    wdot(5) = wdot(5) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(40) = wdot(40) - qdot
    wdot(41) = wdot(41) + qdot

    qdot = q_f(230)-q_r(230)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(40) = wdot(40) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(231)-q_r(231)
    wdot(1) = wdot(1) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(40) = wdot(40) - qdot
    wdot(41) = wdot(41) + qdot

    qdot = q_f(232)-q_r(232)
    wdot(3) = wdot(3) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(47) = wdot(47) - qdot

    qdot = q_f(233)-q_r(233)
    wdot(2) = wdot(2) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(32) = wdot(32) + qdot
    wdot(47) = wdot(47) - qdot

    qdot = q_f(234)-q_r(234)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(47) = wdot(47) - qdot

    qdot = q_f(235)-q_r(235)
    wdot(15) = wdot(15) + qdot
    wdot(31) = wdot(31) - qdot
    wdot(47) = wdot(47) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(236)-q_r(236)
    wdot(4) = wdot(4) - qdot
    wdot(16) = wdot(16) + qdot
    wdot(36) = wdot(36) + qdot
    wdot(47) = wdot(47) - qdot

    qdot = q_f(237)-q_r(237)
    wdot(15) = wdot(15) + qdot
    wdot(36) = wdot(36) - qdot
    wdot(38) = wdot(38) + qdot
    wdot(47) = wdot(47) - qdot

    qdot = q_f(238)-q_r(238)
    wdot(16) = wdot(16) + qdot
    wdot(36) = wdot(36) - qdot
    wdot(47) = wdot(47) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(239)-q_r(239)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(41) = wdot(41) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(240)-q_r(240)
    wdot(3) = wdot(3) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(32) = wdot(32) + qdot
    wdot(41) = wdot(41) - qdot

    qdot = q_f(241)-q_r(241)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(40) = wdot(40) + qdot
    wdot(41) = wdot(41) - qdot

    qdot = q_f(242)-q_r(242)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(41) = wdot(41) - qdot
    wdot(45) = wdot(45) + qdot

    qdot = q_f(243)-q_r(243)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(41) = wdot(41) - qdot
    wdot(46) = wdot(46) + qdot

    qdot = q_f(244)-q_r(244)
    wdot(5) = wdot(5) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(33) = wdot(33) + qdot
    wdot(41) = wdot(41) - qdot

    qdot = q_f(245)-q_r(245)
    wdot(11) = wdot(11) + qdot
    wdot(31) = wdot(31) - qdot
    wdot(42) = wdot(42) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(246)-q_r(246)
    wdot(9) = wdot(9) - qdot
    wdot(31) = wdot(31) + qdot
    wdot(40) = wdot(40) + qdot
    wdot(48) = wdot(48) - qdot

    qdot = q_f(247)-q_r(247)
    wdot(10) = wdot(10) - qdot
    wdot(31) = wdot(31) + qdot
    wdot(41) = wdot(41) + qdot
    wdot(48) = wdot(48) - qdot

    qdot = q_f(248)-q_r(248)
    wdot(11) = wdot(11) - qdot
    wdot(32) = wdot(32) + qdot
    wdot(41) = wdot(41) + qdot
    wdot(48) = wdot(48) - qdot

    qdot = q_f(249)-q_r(249)
    wdot(12) = wdot(12) - qdot
    wdot(32) = wdot(32) + qdot
    wdot(41) = wdot(41) + qdot
    wdot(48) = wdot(48) - qdot

    qdot = q_f(250)-q_r(250)
    wdot(3) = wdot(3) + qdot
    wdot(9) = wdot(9) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(40) = wdot(40) + qdot

    qdot = q_f(251)-q_r(251)
    wdot(9) = wdot(9) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(31) = wdot(31) + qdot
    wdot(36) = wdot(36) - qdot

    qdot = q_f(252)-q_r(252)
    wdot(3) = wdot(3) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(41) = wdot(41) + qdot

    qdot = q_f(253)-q_r(253)
    wdot(2) = wdot(2) + qdot
    wdot(10) = wdot(10) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(254)-q_r(254)
    wdot(10) = wdot(10) - qdot
    wdot(17) = wdot(17) + qdot
    wdot(31) = wdot(31) + qdot
    wdot(36) = wdot(36) - qdot

    qdot = q_f(255)-q_r(255)
    wdot(2) = wdot(2) + qdot
    wdot(11) = wdot(11) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(46) = wdot(46) + qdot

    qdot = q_f(256)-q_r(256)
    wdot(5) = wdot(5) + qdot
    wdot(11) = wdot(11) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(41) = wdot(41) + qdot

    qdot = q_f(257)-q_r(257)
    wdot(2) = wdot(2) + qdot
    wdot(11) = wdot(11) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(44) = wdot(44) + qdot

    qdot = q_f(258)-q_r(258)
    wdot(2) = wdot(2) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(46) = wdot(46) + qdot

    qdot = q_f(259)-q_r(259)
    wdot(5) = wdot(5) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(41) = wdot(41) + qdot

    qdot = q_f(260)-q_r(260)
    wdot(2) = wdot(2) + qdot
    wdot(12) = wdot(12) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(44) = wdot(44) + qdot

    qdot = q_f(261)-q_r(261)
    wdot(6) = wdot(6) + qdot
    wdot(13) = wdot(13) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(41) = wdot(41) + qdot

    qdot = q_f(262)-q_r(262)
    wdot(5) = wdot(5) + qdot
    wdot(13) = wdot(13) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(42) = wdot(42) + qdot

    qdot = q_f(263)-q_r(263)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(43) = wdot(43) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(264)-q_r(264)
    wdot(3) = wdot(3) - qdot
    wdot(36) = wdot(36) + qdot
    wdot(41) = wdot(41) + qdot
    wdot(43) = wdot(43) - qdot

    qdot = q_f(265)-q_r(265)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(17) = wdot(17) + qdot
    wdot(43) = wdot(43) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(266)-q_r(266)
    wdot(2) = wdot(2) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(17) = wdot(17) + qdot
    wdot(43) = wdot(43) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(267)-q_r(267)
    wdot(2) = wdot(2) - qdot
    wdot(11) = wdot(11) + qdot
    wdot(43) = wdot(43) - qdot
    wdot(48) = wdot(48) + qdot

    qdot = q_f(268)-q_r(268)
    wdot(3) = wdot(3) - qdot
    wdot(16) = wdot(16) + qdot
    wdot(32) = wdot(32) + qdot
    wdot(46) = wdot(46) - qdot

    qdot = q_f(269)-q_r(269)
    wdot(3) = wdot(3) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(39) = wdot(39) + qdot
    wdot(46) = wdot(46) - qdot

    qdot = q_f(270)-q_r(270)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(46) = wdot(46) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(271)-q_r(271)
    wdot(2) = wdot(2) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(33) = wdot(33) + qdot
    wdot(46) = wdot(46) - qdot

    qdot = q_f(272)-q_r(272)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(46) = wdot(46) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(273)-q_r(273)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(46) = wdot(46) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(274)-q_r(274)
    wdot(5) = wdot(5) - qdot
    wdot(16) = wdot(16) + qdot
    wdot(33) = wdot(33) + qdot
    wdot(46) = wdot(46) - qdot

    qdot = q_f(275)-q_r(275)
    wdot(2) = wdot(2) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(44) = wdot(44) - qdot
    wdot(46) = wdot(46) + qdot

    qdot = q_f(276)-q_r(276)
    wdot(2) = wdot(2) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(41) = wdot(41) + qdot
    wdot(44) = wdot(44) - qdot

    qdot = q_f(277)-q_r(277)
    wdot(2) = wdot(2) - qdot
    wdot(15) = wdot(15) + qdot
    wdot(33) = wdot(33) + qdot
    wdot(44) = wdot(44) - qdot

    qdot = q_f(278)-q_r(278)
    wdot(2) = wdot(2) - qdot
    wdot(2) = wdot(2) + qdot
    wdot(45) = wdot(45) - qdot
    wdot(46) = wdot(46) + qdot

    qdot = q_f(279)-q_r(279)
    wdot(15) = wdot(15) + qdot
    wdot(28) = wdot(28) - qdot
    wdot(36) = wdot(36) - qdot
    wdot(44) = wdot(44) + qdot

    qdot = q_f(280)-q_r(280)
    wdot(2) = wdot(2) + qdot
    wdot(13) = wdot(13) - qdot
    wdot(31) = wdot(31) - qdot
    wdot(42) = wdot(42) + qdot

    qdot = q_f(281)-q_r(281)
    wdot(1) = wdot(1) + qdot
    wdot(13) = wdot(13) - qdot
    wdot(31) = wdot(31) - qdot
    wdot(41) = wdot(41) + qdot

    qdot = q_f(282)-q_r(282)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(33) = wdot(33) + qdot
    wdot(34) = wdot(34) - qdot

    qdot = q_f(283)-q_r(283)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(33) = wdot(33) + qdot
    wdot(34) = wdot(34) - qdot

    qdot = q_f(284)-q_r(284)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(33) = wdot(33) + qdot
    wdot(34) = wdot(34) - qdot

    qdot = q_f(285)-q_r(285)
    wdot(15) = wdot(15) + qdot
    wdot(16) = wdot(16) - qdot
    wdot(32) = wdot(32) - qdot
    wdot(39) = wdot(39) + qdot

    qdot = q_f(286)-q_r(286)
    wdot(36) = wdot(36) + qdot
    wdot(37) = wdot(37) - qdot
    wdot(40) = wdot(40) - qdot
    wdot(47) = wdot(47) + qdot

    qdot = q_f(287)-q_r(287)
    wdot(16) = wdot(16) + qdot
    wdot(37) = wdot(37) - qdot
    wdot(38) = wdot(38) + qdot
    wdot(47) = wdot(47) - qdot

    qdot = q_f(288)-q_r(288)
    wdot(15) = wdot(15) + qdot
    wdot(16) = wdot(16) - qdot
    wdot(31) = wdot(31) - qdot
    wdot(36) = wdot(36) + qdot

    qdot = q_f(289)-q_r(289)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(15) = wdot(15) + qdot

    qdot = q_f(290)-q_r(290)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(25) = wdot(25) - qdot
    wdot(52) = wdot(52) + qdot

    qdot = q_f(291)-q_r(291)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(26) = wdot(26) - qdot
    wdot(53) = wdot(53) + qdot

    qdot = q_f(292)-q_r(292)
    wdot(4) = wdot(4) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(7) = wdot(7) - qdot

    qdot = q_f(293)-q_r(293)
    wdot(1) = wdot(1) + qdot
    wdot(5) = wdot(5) - qdot
    wdot(13) = wdot(13) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(294)-q_r(294)
    wdot(2) = wdot(2) + (2 * qdot)
    wdot(4) = wdot(4) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(16) = wdot(16) + qdot

    qdot = q_f(295)-q_r(295)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(296)-q_r(296)
    wdot(2) = wdot(2) + (2 * qdot)
    wdot(11) = wdot(11) - qdot
    wdot(11) = wdot(11) - qdot
    wdot(23) = wdot(23) + qdot

    qdot = q_f(297)-q_r(297)
    wdot(1) = wdot(1) + qdot
    wdot(6) = wdot(6) - qdot
    wdot(12) = wdot(12) - qdot
    wdot(18) = wdot(18) + qdot

    qdot = q_f(298)-q_r(298)
    wdot(3) = wdot(3) + qdot
    wdot(4) = wdot(4) - qdot
    wdot(24) = wdot(24) - qdot
    wdot(52) = wdot(52) + qdot

    qdot = q_f(299)-q_r(299)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(23) = wdot(23) + qdot
    wdot(24) = wdot(24) - qdot

    qdot = q_f(300)-q_r(300)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(52) = wdot(52) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(301)-q_r(301)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(302)-q_r(302)
    wdot(4) = wdot(4) - qdot
    wdot(7) = wdot(7) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(303)-q_r(303)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(52) = wdot(52) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(304)-q_r(304)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(305)-q_r(305)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(306)-q_r(306)
    wdot(7) = wdot(7) - qdot
    wdot(8) = wdot(8) + qdot
    wdot(13) = wdot(13) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(307)-q_r(307)
    wdot(13) = wdot(13) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(14) = wdot(14) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(53) = wdot(53) - qdot

    qdot = q_f(308)-q_r(308)
    wdot(2) = wdot(2) + qdot
    wdot(3) = wdot(3) - qdot
    wdot(11) = wdot(11) + qdot
    wdot(16) = wdot(16) + qdot
    wdot(52) = wdot(52) - qdot

    qdot = q_f(309)-q_r(309)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(15) = wdot(15) + qdot
    wdot(18) = wdot(18) + qdot
    wdot(52) = wdot(52) - qdot

    qdot = q_f(310)-q_r(310)
    wdot(4) = wdot(4) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(17) = wdot(17) + (2 * qdot)
    wdot(52) = wdot(52) - qdot

    qdot = q_f(311)-q_r(311)
    wdot(2) = wdot(2) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(17) = wdot(17) + qdot
    wdot(52) = wdot(52) - qdot

    qdot = q_f(312)-q_r(312)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(29) = wdot(29) + qdot
    wdot(52) = wdot(52) - qdot

    qdot = q_f(313)-q_r(313)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(29) = wdot(29) + qdot
    wdot(52) = wdot(52) - qdot

    qdot = q_f(314)-q_r(314)
    wdot(5) = wdot(5) - qdot
    wdot(17) = wdot(17) + qdot
    wdot(19) = wdot(19) + qdot
    wdot(52) = wdot(52) - qdot

    qdot = q_f(315)-q_r(315)
    wdot(3) = wdot(3) - qdot
    wdot(5) = wdot(5) + qdot
    wdot(50) = wdot(50) + qdot
    wdot(51) = wdot(51) - qdot

    qdot = q_f(316)-q_r(316)
    wdot(1) = wdot(1) + qdot
    wdot(2) = wdot(2) - qdot
    wdot(50) = wdot(50) + qdot
    wdot(51) = wdot(51) - qdot

    qdot = q_f(317)-q_r(317)
    wdot(5) = wdot(5) - qdot
    wdot(6) = wdot(6) + qdot
    wdot(50) = wdot(50) + qdot
    wdot(51) = wdot(51) - qdot

    qdot = q_f(318)-q_r(318)
    wdot(7) = wdot(7) + qdot
    wdot(8) = wdot(8) - qdot
    wdot(50) = wdot(50) - qdot
    wdot(51) = wdot(51) + qdot

    qdot = q_f(319)-q_r(319)
    wdot(13) = wdot(13) - qdot
    wdot(14) = wdot(14) + qdot
    wdot(50) = wdot(50) + qdot
    wdot(51) = wdot(51) - qdot

    qdot = q_f(320)-q_r(320)
    wdot(3) = wdot(3) - qdot
    wdot(18) = wdot(18) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(50) = wdot(50) - qdot

    qdot = q_f(321)-q_r(321)
    wdot(2) = wdot(2) - qdot
    wdot(13) = wdot(13) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(50) = wdot(50) - qdot

    qdot = q_f(322)-q_r(322)
    wdot(5) = wdot(5) - qdot
    wdot(19) = wdot(19) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(50) = wdot(50) - qdot

    qdot = q_f(323)-q_r(323)
    wdot(4) = wdot(4) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(50) = wdot(50) - qdot
    wdot(51) = wdot(51) + qdot

    qdot = q_f(324)-q_r(324)
    wdot(5) = wdot(5) + qdot
    wdot(7) = wdot(7) - qdot
    wdot(18) = wdot(18) + qdot
    wdot(26) = wdot(26) + qdot
    wdot(50) = wdot(50) - qdot

    qdot = q_f(325)-q_r(325)
    wdot(13) = wdot(13) - qdot
    wdot(26) = wdot(26) + (2 * qdot)
    wdot(50) = wdot(50) - qdot

end subroutine

subroutine comp_k_f(tc, invT, k_f)

    implicit none

    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT
    double precision, intent(out) :: k_f(325)

    integer :: i

    do i=1, 325
        k_f(i) = prefactor_units(i)*fwd_A(i)*exp(fwd_beta(i)*tc(1)-activation_units(i)*fwd_Ea(i)*invT)
    end do

end subroutine

subroutine comp_Kc(tc, invT, Kc)

    implicit none

    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT
    double precision, intent(inout) :: Kc(325)

    double precision :: g_RT(53)
    double precision :: refC, refCinv
    integer :: i

    ! compute the Gibbs free energy
    call gibbs(g_RT, tc)

    Kc(1) = g_RT(2) + g_RT(11) - g_RT(13)
    Kc(2) = g_RT(2) + g_RT(13) - g_RT(14)
    Kc(3) = g_RT(2) + g_RT(17) - g_RT(18)
    Kc(4) = g_RT(2) + g_RT(18) - g_RT(19)
    Kc(5) = g_RT(2) + g_RT(18) - g_RT(20)
    Kc(6) = g_RT(2) + g_RT(19) - g_RT(21)
    Kc(7) = g_RT(2) + g_RT(20) - g_RT(21)
    Kc(8) = g_RT(2) + g_RT(22) - g_RT(23)
    Kc(9) = g_RT(2) + g_RT(23) - g_RT(24)
    Kc(10) = g_RT(2) + g_RT(24) - g_RT(25)
    Kc(11) = g_RT(2) + g_RT(25) - g_RT(26)
    Kc(12) = g_RT(2) + g_RT(26) - g_RT(27)
    Kc(13) = g_RT(1) + g_RT(15) - g_RT(18)
    Kc(14) = 2*g_RT(5) - g_RT(8)
    Kc(15) = g_RT(5) + g_RT(13) - g_RT(21)
    Kc(16) = g_RT(10) + g_RT(15) - g_RT(28)
    Kc(17) = g_RT(11) + g_RT(15) - g_RT(29)
    Kc(18) = g_RT(6) + g_RT(12) - g_RT(21)
    Kc(19) = 2*g_RT(13) - g_RT(27)
    Kc(20) = -g_RT(1) - g_RT(23) + g_RT(25)
    Kc(21) = g_RT(10) - g_RT(43) + g_RT(48)
    Kc(22) = g_RT(1) + g_RT(10) - g_RT(13)
    Kc(23) = g_RT(2) + g_RT(29) - g_RT(52)
    Kc(24) = g_RT(13) + g_RT(26) - g_RT(51)
    Kc(25) = g_RT(13) + g_RT(25) - g_RT(50)
    Kc(26) = g_RT(2) + g_RT(50) - g_RT(51)
    Kc(27) = g_RT(3) + g_RT(15) - g_RT(16)
    Kc(28) = -g_RT(3) + g_RT(38) - g_RT(48)
    Kc(29) = g_RT(2) + g_RT(41) - g_RT(42)
    Kc(30) = 2*g_RT(3) - g_RT(4)
    Kc(31) = g_RT(2) + g_RT(3) - g_RT(5)
    Kc(32) = g_RT(2) + g_RT(4) - g_RT(7)
    Kc(33) = -g_RT(1) + 2*g_RT(2)
    Kc(34) = g_RT(2) + g_RT(5) - g_RT(6)
    Kc(35) = -g_RT(2) - g_RT(15) + g_RT(17)
    Kc(36) = g_RT(3) + g_RT(36) - g_RT(37)
    Kc(37) = -g_RT(2) + g_RT(35) - g_RT(48)
    Kc(38) = g_RT(2) + g_RT(36) - g_RT(39)
    Kc(39) = -g_RT(15) - g_RT(31) + g_RT(47)
    Kc(40) = -g_RT(2) - g_RT(40) + g_RT(41)
    Kc(41) = -g_RT(15) - g_RT(32) + g_RT(46)
    Kc(42) = g_RT(1) - g_RT(2) + g_RT(3) - g_RT(5)
    Kc(43) = g_RT(3) - g_RT(4) - g_RT(5) + g_RT(7)
    Kc(44) = g_RT(3) - g_RT(5) - g_RT(7) + g_RT(8)
    Kc(45) = -g_RT(2) + g_RT(3) + g_RT(10) - g_RT(15)
    Kc(46) = -g_RT(2) + g_RT(3) + g_RT(11) - g_RT(17)
    Kc(47) = -g_RT(1) + g_RT(3) + g_RT(12) - g_RT(15)
    Kc(48) = -g_RT(2) + g_RT(3) + g_RT(12) - g_RT(17)
    Kc(49) = -g_RT(2) + g_RT(3) + g_RT(13) - g_RT(18)
    Kc(50) = g_RT(3) - g_RT(5) - g_RT(13) + g_RT(14)
    Kc(51) = g_RT(3) - g_RT(5) - g_RT(15) + g_RT(17)
    Kc(52) = -g_RT(2) + g_RT(3) - g_RT(16) + g_RT(17)
    Kc(53) = g_RT(3) - g_RT(5) - g_RT(17) + g_RT(18)
    Kc(54) = g_RT(3) - g_RT(5) - g_RT(18) + g_RT(19)
    Kc(55) = g_RT(3) - g_RT(5) - g_RT(18) + g_RT(20)
    Kc(56) = g_RT(3) - g_RT(5) - g_RT(19) + g_RT(21)
    Kc(57) = g_RT(3) - g_RT(5) - g_RT(20) + g_RT(21)
    Kc(58) = g_RT(3) - g_RT(10) - g_RT(15) + g_RT(22)
    Kc(59) = -g_RT(2) + g_RT(3) + g_RT(23) - g_RT(28)
    Kc(60) = g_RT(3) - g_RT(5) - g_RT(22) + g_RT(23)
    Kc(61) = g_RT(3) - g_RT(11) - g_RT(15) + g_RT(23)
    Kc(62) = -g_RT(2) + g_RT(3) + g_RT(24) - g_RT(29)
    Kc(63) = g_RT(3) - g_RT(13) - g_RT(17) + g_RT(25)
    Kc(64) = g_RT(3) - g_RT(13) - g_RT(18) + g_RT(26)
    Kc(65) = g_RT(3) - g_RT(5) - g_RT(26) + g_RT(27)
    Kc(66) = -g_RT(2) + g_RT(3) - 2*g_RT(15) + g_RT(28)
    Kc(67) = g_RT(3) - g_RT(5) - g_RT(28) + g_RT(29)
    Kc(68) = g_RT(3) - g_RT(11) - g_RT(16) + g_RT(29)
    Kc(69) = -g_RT(3) + g_RT(4) + g_RT(15) - g_RT(16)
    Kc(70) = g_RT(4) - g_RT(7) - g_RT(17) + g_RT(18)
    Kc(71) = g_RT(2) + 2*g_RT(4) - g_RT(4) - g_RT(7)
    Kc(72) = g_RT(2) + g_RT(4) + g_RT(6) - g_RT(6) - g_RT(7)
    Kc(73) = g_RT(2) + g_RT(4) - g_RT(7) + g_RT(48) - g_RT(48)
    Kc(74) = g_RT(2) + g_RT(4) - g_RT(7) + g_RT(49) - g_RT(49)
    Kc(75) = g_RT(2) - g_RT(3) + g_RT(4) - g_RT(5)
    Kc(76) = g_RT(1) - 2*g_RT(1) + 2*g_RT(2)
    Kc(77) = -g_RT(1) + 2*g_RT(2) + g_RT(6) - g_RT(6)
    Kc(78) = -g_RT(1) + 2*g_RT(2) + g_RT(16) - g_RT(16)
    Kc(79) = g_RT(2) - g_RT(3) - g_RT(6) + g_RT(7)
    Kc(80) = -g_RT(1) + g_RT(2) - g_RT(4) + g_RT(7)
    Kc(81) = g_RT(2) - 2*g_RT(5) + g_RT(7)
    Kc(82) = -g_RT(1) + g_RT(2) - g_RT(7) + g_RT(8)
    Kc(83) = g_RT(2) - g_RT(5) - g_RT(6) + g_RT(8)
    Kc(84) = -g_RT(1) + g_RT(2) - g_RT(9) + g_RT(10)
    Kc(85) = -g_RT(1) + g_RT(2) - g_RT(10) + g_RT(12)
    Kc(86) = -g_RT(1) + g_RT(2) - g_RT(13) + g_RT(14)
    Kc(87) = -g_RT(1) + g_RT(2) - g_RT(15) + g_RT(17)
    Kc(88) = -g_RT(1) + g_RT(2) - g_RT(17) + g_RT(18)
    Kc(89) = -g_RT(1) + g_RT(2) - g_RT(18) + g_RT(19)
    Kc(90) = g_RT(2) - g_RT(5) - g_RT(13) + g_RT(19)
    Kc(91) = g_RT(2) - g_RT(6) - g_RT(12) + g_RT(19)
    Kc(92) = g_RT(2) - g_RT(2) - g_RT(19) + g_RT(20)
    Kc(93) = -g_RT(1) + g_RT(2) - g_RT(18) + g_RT(20)
    Kc(94) = g_RT(2) - g_RT(5) - g_RT(13) + g_RT(20)
    Kc(95) = g_RT(2) - g_RT(6) - g_RT(12) + g_RT(20)
    Kc(96) = -g_RT(1) + g_RT(2) - g_RT(19) + g_RT(21)
    Kc(97) = -g_RT(1) + g_RT(2) - g_RT(20) + g_RT(21)
    Kc(98) = -g_RT(1) + g_RT(2) - g_RT(23) + g_RT(24)
    Kc(99) = -g_RT(1) + g_RT(2) - g_RT(24) + g_RT(25)
    Kc(100) = -g_RT(1) + g_RT(2) - g_RT(25) + g_RT(26)
    Kc(101) = -g_RT(1) + g_RT(2) - g_RT(26) + g_RT(27)
    Kc(102) = g_RT(2) - g_RT(12) - g_RT(15) + g_RT(28)
    Kc(103) = -g_RT(1) + g_RT(2) - g_RT(28) + g_RT(29)
    Kc(104) = g_RT(2) - g_RT(13) - g_RT(15) + g_RT(29)
    Kc(105) = g_RT(2) - g_RT(2) - g_RT(29) + g_RT(30)
    Kc(106) = g_RT(1) - g_RT(2) + g_RT(5) - g_RT(6)
    Kc(107) = -g_RT(3) + 2*g_RT(5) - g_RT(6)
    Kc(108) = -g_RT(4) + g_RT(5) - g_RT(6) + g_RT(7)
    Kc(109) = g_RT(5) - g_RT(6) - g_RT(7) + g_RT(8)
    Kc(110) = g_RT(5) - g_RT(6) - g_RT(7) + g_RT(8)
    Kc(111) = -g_RT(2) + g_RT(5) + g_RT(9) - g_RT(15)
    Kc(112) = -g_RT(2) + g_RT(5) + g_RT(10) - g_RT(17)
    Kc(113) = -g_RT(2) + g_RT(5) + g_RT(11) - g_RT(18)
    Kc(114) = g_RT(5) - g_RT(6) - g_RT(10) + g_RT(11)
    Kc(115) = -g_RT(2) + g_RT(5) + g_RT(12) - g_RT(18)
    Kc(116) = g_RT(5) - g_RT(6) - g_RT(11) + g_RT(13)
    Kc(117) = g_RT(5) - g_RT(6) - g_RT(12) + g_RT(13)
    Kc(118) = g_RT(5) - g_RT(6) - g_RT(13) + g_RT(14)
    Kc(119) = -g_RT(2) + g_RT(5) + g_RT(15) - g_RT(16)
    Kc(120) = g_RT(5) - g_RT(6) - g_RT(15) + g_RT(17)
    Kc(121) = g_RT(5) - g_RT(6) - g_RT(17) + g_RT(18)
    Kc(122) = g_RT(5) - g_RT(6) - g_RT(18) + g_RT(19)
    Kc(123) = g_RT(5) - g_RT(6) - g_RT(18) + g_RT(20)
    Kc(124) = g_RT(5) - g_RT(6) - g_RT(19) + g_RT(21)
    Kc(125) = g_RT(5) - g_RT(6) - g_RT(20) + g_RT(21)
    Kc(126) = -g_RT(2) + g_RT(5) + g_RT(22) - g_RT(28)
    Kc(127) = -g_RT(2) + g_RT(5) + g_RT(23) - g_RT(29)
    Kc(128) = -g_RT(2) + g_RT(5) + g_RT(23) - g_RT(30)
    Kc(129) = g_RT(5) - g_RT(6) - g_RT(22) + g_RT(23)
    Kc(130) = g_RT(5) - g_RT(13) - g_RT(15) + g_RT(23)
    Kc(131) = g_RT(5) - g_RT(6) - g_RT(23) + g_RT(24)
    Kc(132) = g_RT(5) - g_RT(6) - g_RT(24) + g_RT(25)
    Kc(133) = g_RT(5) - g_RT(6) - g_RT(26) + g_RT(27)
    Kc(134) = g_RT(5) - g_RT(6) - g_RT(28) + g_RT(29)
    Kc(135) = -g_RT(4) + 2*g_RT(7) - g_RT(8)
    Kc(136) = -g_RT(4) + 2*g_RT(7) - g_RT(8)
    Kc(137) = -g_RT(5) + g_RT(7) + g_RT(11) - g_RT(18)
    Kc(138) = -g_RT(4) + g_RT(7) + g_RT(13) - g_RT(14)
    Kc(139) = -g_RT(5) + g_RT(7) + g_RT(13) - g_RT(20)
    Kc(140) = -g_RT(5) + g_RT(7) + g_RT(15) - g_RT(16)
    Kc(141) = g_RT(7) - g_RT(8) - g_RT(17) + g_RT(18)
    Kc(142) = -g_RT(3) + g_RT(4) + g_RT(9) - g_RT(15)
    Kc(143) = -g_RT(2) + g_RT(9) + g_RT(11) - g_RT(22)
    Kc(144) = -g_RT(2) + g_RT(9) + g_RT(13) - g_RT(23)
    Kc(145) = -g_RT(3) + g_RT(4) + g_RT(10) - g_RT(17)
    Kc(146) = g_RT(1) - g_RT(2) + g_RT(10) - g_RT(11)
    Kc(147) = -g_RT(2) + g_RT(6) + g_RT(10) - g_RT(18)
    Kc(148) = -g_RT(2) + g_RT(10) + g_RT(11) - g_RT(23)
    Kc(149) = -g_RT(2) + g_RT(10) + g_RT(13) - g_RT(24)
    Kc(150) = -g_RT(2) + g_RT(10) + g_RT(14) - g_RT(25)
    Kc(151) = g_RT(10) - g_RT(15) + g_RT(16) - g_RT(17)
    Kc(152) = -g_RT(2) + g_RT(10) + g_RT(18) - g_RT(29)
    Kc(153) = g_RT(10) - g_RT(15) - g_RT(23) + g_RT(28)
    Kc(154) = -g_RT(2) + g_RT(4) - g_RT(5) + g_RT(11) - g_RT(15)
    Kc(155) = g_RT(1) - g_RT(2) + g_RT(11) - g_RT(13)
    Kc(156) = -g_RT(1) + 2*g_RT(11) - g_RT(23)
    Kc(157) = -g_RT(2) + g_RT(11) + g_RT(13) - g_RT(25)
    Kc(158) = g_RT(11) - 2*g_RT(13) + g_RT(14)
    Kc(159) = g_RT(11) - g_RT(15) - g_RT(24) + g_RT(28)
    Kc(160) = -g_RT(11) + g_RT(12) + g_RT(48) - g_RT(48)
    Kc(161) = -g_RT(11) + g_RT(12) + g_RT(49) - g_RT(49)
    Kc(162) = -g_RT(2) + g_RT(4) - g_RT(5) + g_RT(12) - g_RT(15)
    Kc(163) = g_RT(4) - g_RT(6) + g_RT(12) - g_RT(15)
    Kc(164) = g_RT(1) - g_RT(2) + g_RT(12) - g_RT(13)
    Kc(165) = g_RT(6) - g_RT(6) - g_RT(11) + g_RT(12)
    Kc(166) = -g_RT(2) + g_RT(12) + g_RT(13) - g_RT(25)
    Kc(167) = g_RT(12) - 2*g_RT(13) + g_RT(14)
    Kc(168) = -g_RT(11) + g_RT(12) + g_RT(15) - g_RT(15)
    Kc(169) = -g_RT(11) + g_RT(12) + g_RT(16) - g_RT(16)
    Kc(170) = g_RT(12) - g_RT(15) + g_RT(16) - g_RT(18)
    Kc(171) = g_RT(12) - g_RT(13) - g_RT(26) + g_RT(27)
    Kc(172) = -g_RT(3) + g_RT(4) + g_RT(13) - g_RT(20)
    Kc(173) = g_RT(4) - g_RT(5) + g_RT(13) - g_RT(18)
    Kc(174) = -g_RT(7) + g_RT(8) + g_RT(13) - g_RT(14)
    Kc(175) = -g_RT(2) + 2*g_RT(13) - g_RT(26)
    Kc(176) = g_RT(13) - g_RT(14) - g_RT(15) + g_RT(17)
    Kc(177) = g_RT(13) - g_RT(14) - g_RT(17) + g_RT(18)
    Kc(178) = g_RT(13) - g_RT(14) - g_RT(19) + g_RT(21)
    Kc(179) = g_RT(13) - g_RT(14) - g_RT(20) + g_RT(21)
    Kc(180) = g_RT(13) - g_RT(14) - g_RT(24) + g_RT(25)
    Kc(181) = g_RT(13) - g_RT(14) - g_RT(26) + g_RT(27)
    Kc(182) = -g_RT(2) + g_RT(6) - g_RT(6) - g_RT(15) + g_RT(17)
    Kc(183) = g_RT(4) - g_RT(7) - g_RT(15) + g_RT(17)
    Kc(184) = g_RT(4) - g_RT(7) - g_RT(18) + g_RT(19)
    Kc(185) = g_RT(4) - g_RT(7) - g_RT(18) + g_RT(20)
    Kc(186) = g_RT(4) - g_RT(15) - g_RT(17) + g_RT(22)
    Kc(187) = g_RT(1) - g_RT(2) + g_RT(22) - g_RT(23)
    Kc(188) = g_RT(4) - g_RT(17) - g_RT(18) + g_RT(24)
    Kc(189) = g_RT(4) - g_RT(7) - g_RT(25) + g_RT(26)
    Kc(190) = g_RT(4) - g_RT(5) - 2*g_RT(15) + g_RT(28)
    Kc(191) = -2*g_RT(15) - g_RT(23) + 2*g_RT(28)
    Kc(192) = -g_RT(3) + g_RT(31) + g_RT(36) - g_RT(48)
    Kc(193) = -g_RT(3) + g_RT(4) + g_RT(31) - g_RT(36)
    Kc(194) = -g_RT(2) + g_RT(5) + g_RT(31) - g_RT(36)
    Kc(195) = g_RT(3) - g_RT(4) + g_RT(38) - g_RT(48)
    Kc(196) = g_RT(3) - 2*g_RT(36) + g_RT(38)
    Kc(197) = g_RT(2) - g_RT(5) + g_RT(38) - g_RT(48)
    Kc(198) = g_RT(5) - g_RT(7) + g_RT(38) - g_RT(48)
    Kc(199) = -g_RT(5) + g_RT(7) + g_RT(36) - g_RT(37)
    Kc(200) = g_RT(3) - g_RT(4) - g_RT(36) + g_RT(37)
    Kc(201) = g_RT(2) - g_RT(5) - g_RT(36) + g_RT(37)
    Kc(202) = -g_RT(2) + g_RT(3) + g_RT(32) - g_RT(36)
    Kc(203) = -g_RT(1) + g_RT(2) - g_RT(31) + g_RT(32)
    Kc(204) = -g_RT(2) + g_RT(5) + g_RT(32) - g_RT(39)
    Kc(205) = g_RT(5) - g_RT(6) - g_RT(31) + g_RT(32)
    Kc(206) = -g_RT(3) + g_RT(4) + g_RT(32) - g_RT(39)
    Kc(207) = g_RT(4) - g_RT(5) + g_RT(32) - g_RT(36)
    Kc(208) = -g_RT(2) + g_RT(31) + g_RT(32) - g_RT(48)
    Kc(209) = -g_RT(1) + g_RT(6) + g_RT(32) - g_RT(39)
    Kc(210) = -g_RT(5) + g_RT(32) + g_RT(36) - g_RT(48)
    Kc(211) = -g_RT(2) + g_RT(32) + g_RT(36) - g_RT(38)
    Kc(212) = g_RT(3) - g_RT(5) - g_RT(32) + g_RT(33)
    Kc(213) = -g_RT(2) + g_RT(3) + g_RT(33) - g_RT(39)
    Kc(214) = -g_RT(1) + g_RT(2) - g_RT(32) + g_RT(33)
    Kc(215) = g_RT(5) - g_RT(6) - g_RT(32) + g_RT(33)
    Kc(216) = -g_RT(2) + g_RT(35) - g_RT(48)
    Kc(217) = g_RT(4) - g_RT(7) + g_RT(35) - g_RT(48)
    Kc(218) = g_RT(3) - g_RT(5) + g_RT(35) - g_RT(48)
    Kc(219) = g_RT(3) - g_RT(32) + g_RT(35) - g_RT(36)
    Kc(220) = -g_RT(1) + g_RT(2) + g_RT(35) - g_RT(48)
    Kc(221) = g_RT(5) - g_RT(6) + g_RT(35) - g_RT(48)
    Kc(222) = g_RT(13) - g_RT(14) + g_RT(35) - g_RT(48)
    Kc(223) = g_RT(3) - g_RT(5) - g_RT(36) + g_RT(39)
    Kc(224) = -g_RT(1) + g_RT(2) - g_RT(36) + g_RT(39)
    Kc(225) = g_RT(5) - g_RT(6) - g_RT(36) + g_RT(39)
    Kc(226) = g_RT(4) - g_RT(7) - g_RT(36) + g_RT(39)
    Kc(227) = g_RT(3) - g_RT(15) - g_RT(31) + g_RT(40)
    Kc(228) = -g_RT(2) + g_RT(5) + g_RT(40) - g_RT(47)
    Kc(229) = -g_RT(5) + g_RT(6) + g_RT(40) - g_RT(41)
    Kc(230) = -g_RT(3) + g_RT(4) + g_RT(40) - g_RT(47)
    Kc(231) = g_RT(1) - g_RT(2) + g_RT(40) - g_RT(41)
    Kc(232) = g_RT(3) - g_RT(15) - g_RT(36) + g_RT(47)
    Kc(233) = g_RT(2) - g_RT(15) - g_RT(32) + g_RT(47)
    Kc(234) = -g_RT(2) + g_RT(5) - g_RT(15) - g_RT(36) + g_RT(47)
    Kc(235) = -g_RT(15) + g_RT(31) + g_RT(47) - g_RT(48)
    Kc(236) = g_RT(4) - g_RT(16) - g_RT(36) + g_RT(47)
    Kc(237) = -g_RT(15) + g_RT(36) - g_RT(38) + g_RT(47)
    Kc(238) = -g_RT(16) + g_RT(36) + g_RT(47) - g_RT(48)
    Kc(239) = -g_RT(2) + g_RT(3) + g_RT(41) - g_RT(47)
    Kc(240) = g_RT(3) - g_RT(15) - g_RT(32) + g_RT(41)
    Kc(241) = g_RT(3) - g_RT(5) - g_RT(40) + g_RT(41)
    Kc(242) = -g_RT(2) + g_RT(5) + g_RT(41) - g_RT(45)
    Kc(243) = -g_RT(2) + g_RT(5) + g_RT(41) - g_RT(46)
    Kc(244) = g_RT(5) - g_RT(15) - g_RT(33) + g_RT(41)
    Kc(245) = -g_RT(11) + g_RT(31) + g_RT(42) - g_RT(48)
    Kc(246) = g_RT(9) - g_RT(31) - g_RT(40) + g_RT(48)
    Kc(247) = g_RT(10) - g_RT(31) - g_RT(41) + g_RT(48)
    Kc(248) = g_RT(11) - g_RT(32) - g_RT(41) + g_RT(48)
    Kc(249) = g_RT(12) - g_RT(32) - g_RT(41) + g_RT(48)
    Kc(250) = -g_RT(3) + g_RT(9) + g_RT(36) - g_RT(40)
    Kc(251) = g_RT(9) - g_RT(15) - g_RT(31) + g_RT(36)
    Kc(252) = -g_RT(3) + g_RT(10) + g_RT(36) - g_RT(41)
    Kc(253) = -g_RT(2) + g_RT(10) + g_RT(36) - g_RT(47)
    Kc(254) = g_RT(10) - g_RT(17) - g_RT(31) + g_RT(36)
    Kc(255) = -g_RT(2) + g_RT(11) + g_RT(36) - g_RT(46)
    Kc(256) = -g_RT(5) + g_RT(11) + g_RT(36) - g_RT(41)
    Kc(257) = -g_RT(2) + g_RT(11) + g_RT(36) - g_RT(44)
    Kc(258) = -g_RT(2) + g_RT(12) + g_RT(36) - g_RT(46)
    Kc(259) = -g_RT(5) + g_RT(12) + g_RT(36) - g_RT(41)
    Kc(260) = -g_RT(2) + g_RT(12) + g_RT(36) - g_RT(44)
    Kc(261) = -g_RT(6) + g_RT(13) + g_RT(36) - g_RT(41)
    Kc(262) = -g_RT(5) + g_RT(13) + g_RT(36) - g_RT(42)
    Kc(263) = -g_RT(2) + g_RT(3) - g_RT(15) + g_RT(43) - g_RT(48)
    Kc(264) = g_RT(3) - g_RT(36) - g_RT(41) + g_RT(43)
    Kc(265) = -g_RT(3) + g_RT(4) - g_RT(17) + g_RT(43) - g_RT(48)
    Kc(266) = -g_RT(2) + g_RT(5) - g_RT(17) + g_RT(43) - g_RT(48)
    Kc(267) = g_RT(2) - g_RT(11) + g_RT(43) - g_RT(48)
    Kc(268) = g_RT(3) - g_RT(16) - g_RT(32) + g_RT(46)
    Kc(269) = g_RT(3) - g_RT(15) - g_RT(39) + g_RT(46)
    Kc(270) = g_RT(3) - g_RT(5) + g_RT(46) - g_RT(47)
    Kc(271) = g_RT(2) - g_RT(15) - g_RT(33) + g_RT(46)
    Kc(272) = -g_RT(1) + g_RT(2) + g_RT(46) - g_RT(47)
    Kc(273) = g_RT(5) - g_RT(6) + g_RT(46) - g_RT(47)
    Kc(274) = g_RT(5) - g_RT(16) - g_RT(33) + g_RT(46)
    Kc(275) = g_RT(2) - g_RT(2) + g_RT(44) - g_RT(46)
    Kc(276) = g_RT(2) - g_RT(5) - g_RT(41) + g_RT(44)
    Kc(277) = g_RT(2) - g_RT(15) - g_RT(33) + g_RT(44)
    Kc(278) = g_RT(2) - g_RT(2) + g_RT(45) - g_RT(46)
    Kc(279) = -g_RT(15) + g_RT(28) + g_RT(36) - g_RT(44)
    Kc(280) = -g_RT(2) + g_RT(13) + g_RT(31) - g_RT(42)
    Kc(281) = -g_RT(1) + g_RT(13) + g_RT(31) - g_RT(41)
    Kc(282) = -g_RT(1) + g_RT(2) - g_RT(33) + g_RT(34)
    Kc(283) = g_RT(5) - g_RT(6) - g_RT(33) + g_RT(34)
    Kc(284) = g_RT(3) - g_RT(5) - g_RT(33) + g_RT(34)
    Kc(285) = -g_RT(15) + g_RT(16) + g_RT(32) - g_RT(39)
    Kc(286) = -g_RT(36) + g_RT(37) + g_RT(40) - g_RT(47)
    Kc(287) = -g_RT(16) + g_RT(37) - g_RT(38) + g_RT(47)
    Kc(288) = -g_RT(15) + g_RT(16) + g_RT(31) - g_RT(36)
    Kc(289) = -g_RT(1) - g_RT(2) + g_RT(3) + g_RT(13) - g_RT(15)
    Kc(290) = -g_RT(2) + g_RT(3) + g_RT(25) - g_RT(52)
    Kc(291) = -g_RT(2) + g_RT(3) + g_RT(26) - g_RT(53)
    Kc(292) = -g_RT(4) + g_RT(5) - g_RT(6) + g_RT(7)
    Kc(293) = -g_RT(1) + g_RT(5) + g_RT(13) - g_RT(18)
    Kc(294) = -2*g_RT(2) + g_RT(4) + g_RT(11) - g_RT(16)
    Kc(295) = -g_RT(3) + g_RT(4) + g_RT(11) - g_RT(18)
    Kc(296) = -2*g_RT(2) + g_RT(11) + g_RT(11) - g_RT(23)
    Kc(297) = -g_RT(1) + g_RT(6) + g_RT(12) - g_RT(18)
    Kc(298) = -g_RT(3) + g_RT(4) + g_RT(24) - g_RT(52)
    Kc(299) = g_RT(4) - g_RT(7) - g_RT(23) + g_RT(24)
    Kc(300) = g_RT(3) - g_RT(5) - g_RT(52) + g_RT(53)
    Kc(301) = g_RT(3) - g_RT(5) - g_RT(13) - g_RT(15) + g_RT(53)
    Kc(302) = g_RT(4) - g_RT(7) - g_RT(13) - g_RT(15) + g_RT(53)
    Kc(303) = -g_RT(1) + g_RT(2) - g_RT(52) + g_RT(53)
    Kc(304) = -g_RT(1) + g_RT(2) - g_RT(13) - g_RT(15) + g_RT(53)
    Kc(305) = g_RT(5) - g_RT(6) - g_RT(13) - g_RT(15) + g_RT(53)
    Kc(306) = g_RT(7) - g_RT(8) - g_RT(13) - g_RT(15) + g_RT(53)
    Kc(307) = g_RT(13) - g_RT(13) - g_RT(14) - g_RT(15) + g_RT(53)
    Kc(308) = -g_RT(2) + g_RT(3) - g_RT(11) - g_RT(16) + g_RT(52)
    Kc(309) = g_RT(4) - g_RT(5) - g_RT(15) - g_RT(18) + g_RT(52)
    Kc(310) = g_RT(4) - g_RT(5) - 2*g_RT(17) + g_RT(52)
    Kc(311) = g_RT(2) - g_RT(13) - g_RT(17) + g_RT(52)
    Kc(312) = -g_RT(1) + g_RT(2) - g_RT(29) + g_RT(52)
    Kc(313) = g_RT(5) - g_RT(6) - g_RT(29) + g_RT(52)
    Kc(314) = g_RT(5) - g_RT(17) - g_RT(19) + g_RT(52)
    Kc(315) = g_RT(3) - g_RT(5) - g_RT(50) + g_RT(51)
    Kc(316) = -g_RT(1) + g_RT(2) - g_RT(50) + g_RT(51)
    Kc(317) = g_RT(5) - g_RT(6) - g_RT(50) + g_RT(51)
    Kc(318) = -g_RT(7) + g_RT(8) + g_RT(50) - g_RT(51)
    Kc(319) = g_RT(13) - g_RT(14) - g_RT(50) + g_RT(51)
    Kc(320) = g_RT(3) - g_RT(18) - g_RT(26) + g_RT(50)
    Kc(321) = g_RT(2) - g_RT(13) - g_RT(26) + g_RT(50)
    Kc(322) = g_RT(5) - g_RT(19) - g_RT(26) + g_RT(50)
    Kc(323) = -g_RT(4) + g_RT(7) + g_RT(50) - g_RT(51)
    Kc(324) = -g_RT(5) + g_RT(7) - g_RT(18) - g_RT(26) + g_RT(50)
    Kc(325) = g_RT(13) - 2*g_RT(26) + g_RT(50)

    do i=1, 325
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
    Kc(14) = Kc(14) * (refCinv)
    Kc(15) = Kc(15) * (refCinv)
    Kc(16) = Kc(16) * (refCinv)
    Kc(17) = Kc(17) * (refCinv)
    Kc(18) = Kc(18) * (refCinv)
    Kc(19) = Kc(19) * (refCinv)
    Kc(20) = Kc(20) * (refC)
    Kc(21) = Kc(21) * (refCinv)
    Kc(22) = Kc(22) * (refCinv)
    Kc(23) = Kc(23) * (refCinv)
    Kc(24) = Kc(24) * (refCinv)
    Kc(25) = Kc(25) * (refCinv)
    Kc(26) = Kc(26) * (refCinv)
    Kc(27) = Kc(27) * (refCinv)
    Kc(28) = Kc(28) * (refC)
    Kc(29) = Kc(29) * (refCinv)
    Kc(30) = Kc(30) * (refCinv)
    Kc(31) = Kc(31) * (refCinv)
    Kc(32) = Kc(32) * (refCinv)
    Kc(33) = Kc(33) * (refCinv)
    Kc(34) = Kc(34) * (refCinv)
    Kc(35) = Kc(35) * (refC)
    Kc(36) = Kc(36) * (refCinv)
    Kc(37) = Kc(37) * (refC)
    Kc(38) = Kc(38) * (refCinv)
    Kc(39) = Kc(39) * (refC)
    Kc(40) = Kc(40) * (refC)
    Kc(41) = Kc(41) * (refC)
    Kc(66) = Kc(66) * (refC)
    Kc(71) = Kc(71) * (refCinv)
    Kc(72) = Kc(72) * (refCinv)
    Kc(73) = Kc(73) * (refCinv)
    Kc(74) = Kc(74) * (refCinv)
    Kc(76) = Kc(76) * (refCinv)
    Kc(77) = Kc(77) * (refCinv)
    Kc(78) = Kc(78) * (refCinv)
    Kc(154) = Kc(154) * (refC)
    Kc(162) = Kc(162) * (refC)
    Kc(182) = Kc(182) * (refC)
    Kc(190) = Kc(190) * (refC)
    Kc(191) = Kc(191) * (refC)
    Kc(216) = Kc(216) * (refC)
    Kc(234) = Kc(234) * (refC)
    Kc(263) = Kc(263) * (refC)
    Kc(265) = Kc(265) * (refC)
    Kc(266) = Kc(266) * (refC)
    Kc(289) = Kc(289) * (refC)
    Kc(294) = Kc(294) * (refC)
    Kc(296) = Kc(296) * (refC)
    Kc(301) = Kc(301) * (refC)
    Kc(302) = Kc(302) * (refC)
    Kc(304) = Kc(304) * (refC)
    Kc(305) = Kc(305) * (refC)
    Kc(306) = Kc(306) * (refC)
    Kc(307) = Kc(307) * (refC)
    Kc(308) = Kc(308) * (refC)
    Kc(309) = Kc(309) * (refC)
    Kc(310) = Kc(310) * (refC)
    Kc(324) = Kc(324) * (refC)

end subroutine

subroutine comp_qfqr(qf, qr, sc, tc, invT)

    implicit none

    double precision, intent(out) :: qf(325)
    double precision, intent(out) :: qr(325)
    double precision, intent(in) :: sc(53)
    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT

    double precision :: T
    double precision :: mixture
    double precision :: Corr(325)
    double precision :: alpha_troe(26)
    double precision :: redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe
    double precision :: alpha_lindemann(3)
    double precision :: tmp1, tmp2, tmp3
    integer :: i

    ! reaction 1: H + CH2 (+M) <=> CH3 (+M)
    qf(1) = sc(2)*sc(11)
    qr(1) = sc(13)
    ! reaction 2: H + CH3 (+M) <=> CH4 (+M)
    qf(2) = sc(2)*sc(13)
    qr(2) = sc(14)
    ! reaction 3: H + HCO (+M) <=> CH2O (+M)
    qf(3) = sc(2)*sc(17)
    qr(3) = sc(18)
    ! reaction 4: H + CH2O (+M) <=> CH2OH (+M)
    qf(4) = sc(2)*sc(18)
    qr(4) = sc(19)
    ! reaction 5: H + CH2O (+M) <=> CH3O (+M)
    qf(5) = sc(2)*sc(18)
    qr(5) = sc(20)
    ! reaction 6: H + CH2OH (+M) <=> CH3OH (+M)
    qf(6) = sc(2)*sc(19)
    qr(6) = sc(21)
    ! reaction 7: H + CH3O (+M) <=> CH3OH (+M)
    qf(7) = sc(2)*sc(20)
    qr(7) = sc(21)
    ! reaction 8: H + C2H (+M) <=> C2H2 (+M)
    qf(8) = sc(2)*sc(22)
    qr(8) = sc(23)
    ! reaction 9: H + C2H2 (+M) <=> C2H3 (+M)
    qf(9) = sc(2)*sc(23)
    qr(9) = sc(24)
    ! reaction 10: H + C2H3 (+M) <=> C2H4 (+M)
    qf(10) = sc(2)*sc(24)
    qr(10) = sc(25)
    ! reaction 11: H + C2H4 (+M) <=> C2H5 (+M)
    qf(11) = sc(2)*sc(25)
    qr(11) = sc(26)
    ! reaction 12: H + C2H5 (+M) <=> C2H6 (+M)
    qf(12) = sc(2)*sc(26)
    qr(12) = sc(27)
    ! reaction 13: H2 + CO (+M) <=> CH2O (+M)
    qf(13) = sc(1)*sc(15)
    qr(13) = sc(18)
    ! reaction 14: 2 OH (+M) <=> H2O2 (+M)
    qf(14) = sc(5)*sc(5)
    qr(14) = sc(8)
    ! reaction 15: OH + CH3 (+M) <=> CH3OH (+M)
    qf(15) = sc(5)*sc(13)
    qr(15) = sc(21)
    ! reaction 16: CH + CO (+M) <=> HCCO (+M)
    qf(16) = sc(10)*sc(15)
    qr(16) = sc(28)
    ! reaction 17: CH2 + CO (+M) <=> CH2CO (+M)
    qf(17) = sc(11)*sc(15)
    qr(17) = sc(29)
    ! reaction 18: CH2(S) + H2O (+M) <=> CH3OH (+M)
    qf(18) = sc(6)*sc(12)
    qr(18) = sc(21)
    ! reaction 19: 2 CH3 (+M) <=> C2H6 (+M)
    qf(19) = sc(13)*sc(13)
    qr(19) = sc(27)
    ! reaction 20: C2H4 (+M) <=> H2 + C2H2 (+M)
    qf(20) = sc(25)
    qr(20) = sc(1)*sc(23)
    ! reaction 21: CH + N2 (+M) <=> HCNN (+M)
    qf(21) = sc(10)*sc(48)
    qr(21) = sc(43)
    ! reaction 22: CH + H2 (+M) <=> CH3 (+M)
    qf(22) = sc(1)*sc(10)
    qr(22) = sc(13)
    ! reaction 23: H + CH2CO (+M) <=> CH2CHO (+M)
    qf(23) = sc(2)*sc(29)
    qr(23) = sc(52)
    ! reaction 24: CH3 + C2H5 (+M) <=> C3H8 (+M)
    qf(24) = sc(13)*sc(26)
    qr(24) = sc(51)
    ! reaction 25: CH3 + C2H4 (+M) <=> C3H7 (+M)
    qf(25) = sc(13)*sc(25)
    qr(25) = sc(50)
    ! reaction 26: H + C3H7 (+M) <=> C3H8 (+M)
    qf(26) = sc(2)*sc(50)
    qr(26) = sc(51)
    ! reaction 27: O + CO (+M) <=> CO2 (+M)
    qf(27) = sc(3)*sc(15)
    qr(27) = sc(16)
    ! reaction 28: N2O (+M) <=> N2 + O (+M)
    qf(28) = sc(38)
    qr(28) = sc(3)*sc(48)
    ! reaction 29: H + HCN (+M) <=> H2CN (+M)
    qf(29) = sc(2)*sc(41)
    qr(29) = sc(42)
    ! reaction 30: 2 O + M <=> O2 + M
    qf(30) = sc(3)*sc(3)
    qr(30) = sc(4)
    ! reaction 31: O + H + M <=> OH + M
    qf(31) = sc(2)*sc(3)
    qr(31) = sc(5)
    ! reaction 32: H + O2 + M <=> HO2 + M
    qf(32) = sc(2)*sc(4)
    qr(32) = sc(7)
    ! reaction 33: 2 H + M <=> H2 + M
    qf(33) = sc(2)*sc(2)
    qr(33) = sc(1)
    ! reaction 34: H + OH + M <=> H2O + M
    qf(34) = sc(2)*sc(5)
    qr(34) = sc(6)
    ! reaction 35: HCO + M <=> H + CO + M
    qf(35) = sc(17)
    qr(35) = sc(2)*sc(15)
    ! reaction 36: NO + O + M <=> NO2 + M
    qf(36) = sc(3)*sc(36)
    qr(36) = sc(37)
    ! reaction 37: NNH + M <=> N2 + H + M
    qf(37) = sc(35)
    qr(37) = sc(2)*sc(48)
    ! reaction 38: H + NO + M <=> HNO + M
    qf(38) = sc(2)*sc(36)
    qr(38) = sc(39)
    ! reaction 39: NCO + M <=> N + CO + M
    qf(39) = sc(47)
    qr(39) = sc(15)*sc(31)
    ! reaction 40: HCN + M <=> H + CN + M
    qf(40) = sc(41)
    qr(40) = sc(2)*sc(40)
    ! reaction 41: HNCO + M <=> NH + CO + M
    qf(41) = sc(46)
    qr(41) = sc(15)*sc(32)
    ! reaction 42: O + H2 <=> H + OH
    qf(42) = sc(1)*sc(3)
    qr(42) = sc(2)*sc(5)
    ! reaction 43: O + HO2 <=> OH + O2
    qf(43) = sc(3)*sc(7)
    qr(43) = sc(4)*sc(5)
    ! reaction 44: O + H2O2 <=> OH + HO2
    qf(44) = sc(3)*sc(8)
    qr(44) = sc(5)*sc(7)
    ! reaction 45: O + CH <=> H + CO
    qf(45) = sc(3)*sc(10)
    qr(45) = sc(2)*sc(15)
    ! reaction 46: O + CH2 <=> H + HCO
    qf(46) = sc(3)*sc(11)
    qr(46) = sc(2)*sc(17)
    ! reaction 47: O + CH2(S) <=> H2 + CO
    qf(47) = sc(3)*sc(12)
    qr(47) = sc(1)*sc(15)
    ! reaction 48: O + CH2(S) <=> H + HCO
    qf(48) = sc(3)*sc(12)
    qr(48) = sc(2)*sc(17)
    ! reaction 49: O + CH3 <=> H + CH2O
    qf(49) = sc(3)*sc(13)
    qr(49) = sc(2)*sc(18)
    ! reaction 50: O + CH4 <=> OH + CH3
    qf(50) = sc(3)*sc(14)
    qr(50) = sc(5)*sc(13)
    ! reaction 51: O + HCO <=> OH + CO
    qf(51) = sc(3)*sc(17)
    qr(51) = sc(5)*sc(15)
    ! reaction 52: O + HCO <=> H + CO2
    qf(52) = sc(3)*sc(17)
    qr(52) = sc(2)*sc(16)
    ! reaction 53: O + CH2O <=> OH + HCO
    qf(53) = sc(3)*sc(18)
    qr(53) = sc(5)*sc(17)
    ! reaction 54: O + CH2OH <=> OH + CH2O
    qf(54) = sc(3)*sc(19)
    qr(54) = sc(5)*sc(18)
    ! reaction 55: O + CH3O <=> OH + CH2O
    qf(55) = sc(3)*sc(20)
    qr(55) = sc(5)*sc(18)
    ! reaction 56: O + CH3OH <=> OH + CH2OH
    qf(56) = sc(3)*sc(21)
    qr(56) = sc(5)*sc(19)
    ! reaction 57: O + CH3OH <=> OH + CH3O
    qf(57) = sc(3)*sc(21)
    qr(57) = sc(5)*sc(20)
    ! reaction 58: O + C2H <=> CH + CO
    qf(58) = sc(3)*sc(22)
    qr(58) = sc(10)*sc(15)
    ! reaction 59: O + C2H2 <=> H + HCCO
    qf(59) = sc(3)*sc(23)
    qr(59) = sc(2)*sc(28)
    ! reaction 60: O + C2H2 <=> OH + C2H
    qf(60) = sc(3)*sc(23)
    qr(60) = sc(5)*sc(22)
    ! reaction 61: O + C2H2 <=> CO + CH2
    qf(61) = sc(3)*sc(23)
    qr(61) = sc(11)*sc(15)
    ! reaction 62: O + C2H3 <=> H + CH2CO
    qf(62) = sc(3)*sc(24)
    qr(62) = sc(2)*sc(29)
    ! reaction 63: O + C2H4 <=> CH3 + HCO
    qf(63) = sc(3)*sc(25)
    qr(63) = sc(13)*sc(17)
    ! reaction 64: O + C2H5 <=> CH3 + CH2O
    qf(64) = sc(3)*sc(26)
    qr(64) = sc(13)*sc(18)
    ! reaction 65: O + C2H6 <=> OH + C2H5
    qf(65) = sc(3)*sc(27)
    qr(65) = sc(5)*sc(26)
    ! reaction 66: O + HCCO <=> H + 2 CO
    qf(66) = sc(3)*sc(28)
    qr(66) = sc(2)*sc(15)*sc(15)
    ! reaction 67: O + CH2CO <=> OH + HCCO
    qf(67) = sc(3)*sc(29)
    qr(67) = sc(5)*sc(28)
    ! reaction 68: O + CH2CO <=> CH2 + CO2
    qf(68) = sc(3)*sc(29)
    qr(68) = sc(11)*sc(16)
    ! reaction 69: O2 + CO <=> O + CO2
    qf(69) = sc(4)*sc(15)
    qr(69) = sc(3)*sc(16)
    ! reaction 70: O2 + CH2O <=> HO2 + HCO
    qf(70) = sc(4)*sc(18)
    qr(70) = sc(7)*sc(17)
    ! reaction 71: H + 2 O2 <=> HO2 + O2
    qf(71) = sc(2)*sc(4)*sc(4)
    qr(71) = sc(4)*sc(7)
    ! reaction 72: H + O2 + H2O <=> HO2 + H2O
    qf(72) = sc(2)*sc(4)*sc(6)
    qr(72) = sc(6)*sc(7)
    ! reaction 73: H + O2 + N2 <=> HO2 + N2
    qf(73) = sc(2)*sc(4)*sc(48)
    qr(73) = sc(7)*sc(48)
    ! reaction 74: H + O2 + AR <=> HO2 + AR
    qf(74) = sc(2)*sc(4)*sc(49)
    qr(74) = sc(7)*sc(49)
    ! reaction 75: H + O2 <=> O + OH
    qf(75) = sc(2)*sc(4)
    qr(75) = sc(3)*sc(5)
    ! reaction 76: 2 H + H2 <=> 2 H2
    qf(76) = sc(1)*sc(2)*sc(2)
    qr(76) = sc(1)*sc(1)
    ! reaction 77: 2 H + H2O <=> H2 + H2O
    qf(77) = sc(2)*sc(2)*sc(6)
    qr(77) = sc(1)*sc(6)
    ! reaction 78: 2 H + CO2 <=> H2 + CO2
    qf(78) = sc(2)*sc(2)*sc(16)
    qr(78) = sc(1)*sc(16)
    ! reaction 79: H + HO2 <=> O + H2O
    qf(79) = sc(2)*sc(7)
    qr(79) = sc(3)*sc(6)
    ! reaction 80: H + HO2 <=> O2 + H2
    qf(80) = sc(2)*sc(7)
    qr(80) = sc(1)*sc(4)
    ! reaction 81: H + HO2 <=> 2 OH
    qf(81) = sc(2)*sc(7)
    qr(81) = sc(5)*sc(5)
    ! reaction 82: H + H2O2 <=> HO2 + H2
    qf(82) = sc(2)*sc(8)
    qr(82) = sc(1)*sc(7)
    ! reaction 83: H + H2O2 <=> OH + H2O
    qf(83) = sc(2)*sc(8)
    qr(83) = sc(5)*sc(6)
    ! reaction 84: H + CH <=> C + H2
    qf(84) = sc(2)*sc(10)
    qr(84) = sc(1)*sc(9)
    ! reaction 85: H + CH2(S) <=> CH + H2
    qf(85) = sc(2)*sc(12)
    qr(85) = sc(1)*sc(10)
    ! reaction 86: H + CH4 <=> CH3 + H2
    qf(86) = sc(2)*sc(14)
    qr(86) = sc(1)*sc(13)
    ! reaction 87: H + HCO <=> H2 + CO
    qf(87) = sc(2)*sc(17)
    qr(87) = sc(1)*sc(15)
    ! reaction 88: H + CH2O <=> HCO + H2
    qf(88) = sc(2)*sc(18)
    qr(88) = sc(1)*sc(17)
    ! reaction 89: H + CH2OH <=> H2 + CH2O
    qf(89) = sc(2)*sc(19)
    qr(89) = sc(1)*sc(18)
    ! reaction 90: H + CH2OH <=> OH + CH3
    qf(90) = sc(2)*sc(19)
    qr(90) = sc(5)*sc(13)
    ! reaction 91: H + CH2OH <=> CH2(S) + H2O
    qf(91) = sc(2)*sc(19)
    qr(91) = sc(6)*sc(12)
    ! reaction 92: H + CH3O <=> H + CH2OH
    qf(92) = sc(2)*sc(20)
    qr(92) = sc(2)*sc(19)
    ! reaction 93: H + CH3O <=> H2 + CH2O
    qf(93) = sc(2)*sc(20)
    qr(93) = sc(1)*sc(18)
    ! reaction 94: H + CH3O <=> OH + CH3
    qf(94) = sc(2)*sc(20)
    qr(94) = sc(5)*sc(13)
    ! reaction 95: H + CH3O <=> CH2(S) + H2O
    qf(95) = sc(2)*sc(20)
    qr(95) = sc(6)*sc(12)
    ! reaction 96: H + CH3OH <=> CH2OH + H2
    qf(96) = sc(2)*sc(21)
    qr(96) = sc(1)*sc(19)
    ! reaction 97: H + CH3OH <=> CH3O + H2
    qf(97) = sc(2)*sc(21)
    qr(97) = sc(1)*sc(20)
    ! reaction 98: H + C2H3 <=> H2 + C2H2
    qf(98) = sc(2)*sc(24)
    qr(98) = sc(1)*sc(23)
    ! reaction 99: H + C2H4 <=> C2H3 + H2
    qf(99) = sc(2)*sc(25)
    qr(99) = sc(1)*sc(24)
    ! reaction 100: H + C2H5 <=> H2 + C2H4
    qf(100) = sc(2)*sc(26)
    qr(100) = sc(1)*sc(25)
    ! reaction 101: H + C2H6 <=> C2H5 + H2
    qf(101) = sc(2)*sc(27)
    qr(101) = sc(1)*sc(26)
    ! reaction 102: H + HCCO <=> CH2(S) + CO
    qf(102) = sc(2)*sc(28)
    qr(102) = sc(12)*sc(15)
    ! reaction 103: H + CH2CO <=> HCCO + H2
    qf(103) = sc(2)*sc(29)
    qr(103) = sc(1)*sc(28)
    ! reaction 104: H + CH2CO <=> CH3 + CO
    qf(104) = sc(2)*sc(29)
    qr(104) = sc(13)*sc(15)
    ! reaction 105: H + HCCOH <=> H + CH2CO
    qf(105) = sc(2)*sc(30)
    qr(105) = sc(2)*sc(29)
    ! reaction 106: OH + H2 <=> H + H2O
    qf(106) = sc(1)*sc(5)
    qr(106) = sc(2)*sc(6)
    ! reaction 107: 2 OH <=> O + H2O
    qf(107) = sc(5)*sc(5)
    qr(107) = sc(3)*sc(6)
    ! reaction 108: OH + HO2 <=> O2 + H2O
    qf(108) = sc(5)*sc(7)
    qr(108) = sc(4)*sc(6)
    ! reaction 109: OH + H2O2 <=> HO2 + H2O
    qf(109) = sc(5)*sc(8)
    qr(109) = sc(6)*sc(7)
    ! reaction 110: OH + H2O2 <=> HO2 + H2O
    qf(110) = sc(5)*sc(8)
    qr(110) = sc(6)*sc(7)
    ! reaction 111: OH + C <=> H + CO
    qf(111) = sc(5)*sc(9)
    qr(111) = sc(2)*sc(15)
    ! reaction 112: OH + CH <=> H + HCO
    qf(112) = sc(5)*sc(10)
    qr(112) = sc(2)*sc(17)
    ! reaction 113: OH + CH2 <=> H + CH2O
    qf(113) = sc(5)*sc(11)
    qr(113) = sc(2)*sc(18)
    ! reaction 114: OH + CH2 <=> CH + H2O
    qf(114) = sc(5)*sc(11)
    qr(114) = sc(6)*sc(10)
    ! reaction 115: OH + CH2(S) <=> H + CH2O
    qf(115) = sc(5)*sc(12)
    qr(115) = sc(2)*sc(18)
    ! reaction 116: OH + CH3 <=> CH2 + H2O
    qf(116) = sc(5)*sc(13)
    qr(116) = sc(6)*sc(11)
    ! reaction 117: OH + CH3 <=> CH2(S) + H2O
    qf(117) = sc(5)*sc(13)
    qr(117) = sc(6)*sc(12)
    ! reaction 118: OH + CH4 <=> CH3 + H2O
    qf(118) = sc(5)*sc(14)
    qr(118) = sc(6)*sc(13)
    ! reaction 119: OH + CO <=> H + CO2
    qf(119) = sc(5)*sc(15)
    qr(119) = sc(2)*sc(16)
    ! reaction 120: OH + HCO <=> H2O + CO
    qf(120) = sc(5)*sc(17)
    qr(120) = sc(6)*sc(15)
    ! reaction 121: OH + CH2O <=> HCO + H2O
    qf(121) = sc(5)*sc(18)
    qr(121) = sc(6)*sc(17)
    ! reaction 122: OH + CH2OH <=> H2O + CH2O
    qf(122) = sc(5)*sc(19)
    qr(122) = sc(6)*sc(18)
    ! reaction 123: OH + CH3O <=> H2O + CH2O
    qf(123) = sc(5)*sc(20)
    qr(123) = sc(6)*sc(18)
    ! reaction 124: OH + CH3OH <=> CH2OH + H2O
    qf(124) = sc(5)*sc(21)
    qr(124) = sc(6)*sc(19)
    ! reaction 125: OH + CH3OH <=> CH3O + H2O
    qf(125) = sc(5)*sc(21)
    qr(125) = sc(6)*sc(20)
    ! reaction 126: OH + C2H <=> H + HCCO
    qf(126) = sc(5)*sc(22)
    qr(126) = sc(2)*sc(28)
    ! reaction 127: OH + C2H2 <=> H + CH2CO
    qf(127) = sc(5)*sc(23)
    qr(127) = sc(2)*sc(29)
    ! reaction 128: OH + C2H2 <=> H + HCCOH
    qf(128) = sc(5)*sc(23)
    qr(128) = sc(2)*sc(30)
    ! reaction 129: OH + C2H2 <=> C2H + H2O
    qf(129) = sc(5)*sc(23)
    qr(129) = sc(6)*sc(22)
    ! reaction 130: OH + C2H2 <=> CH3 + CO
    qf(130) = sc(5)*sc(23)
    qr(130) = sc(13)*sc(15)
    ! reaction 131: OH + C2H3 <=> H2O + C2H2
    qf(131) = sc(5)*sc(24)
    qr(131) = sc(6)*sc(23)
    ! reaction 132: OH + C2H4 <=> C2H3 + H2O
    qf(132) = sc(5)*sc(25)
    qr(132) = sc(6)*sc(24)
    ! reaction 133: OH + C2H6 <=> C2H5 + H2O
    qf(133) = sc(5)*sc(27)
    qr(133) = sc(6)*sc(26)
    ! reaction 134: OH + CH2CO <=> HCCO + H2O
    qf(134) = sc(5)*sc(29)
    qr(134) = sc(6)*sc(28)
    ! reaction 135: 2 HO2 <=> O2 + H2O2
    qf(135) = sc(7)*sc(7)
    qr(135) = sc(4)*sc(8)
    ! reaction 136: 2 HO2 <=> O2 + H2O2
    qf(136) = sc(7)*sc(7)
    qr(136) = sc(4)*sc(8)
    ! reaction 137: HO2 + CH2 <=> OH + CH2O
    qf(137) = sc(7)*sc(11)
    qr(137) = sc(5)*sc(18)
    ! reaction 138: HO2 + CH3 <=> O2 + CH4
    qf(138) = sc(7)*sc(13)
    qr(138) = sc(4)*sc(14)
    ! reaction 139: HO2 + CH3 <=> OH + CH3O
    qf(139) = sc(7)*sc(13)
    qr(139) = sc(5)*sc(20)
    ! reaction 140: HO2 + CO <=> OH + CO2
    qf(140) = sc(7)*sc(15)
    qr(140) = sc(5)*sc(16)
    ! reaction 141: HO2 + CH2O <=> HCO + H2O2
    qf(141) = sc(7)*sc(18)
    qr(141) = sc(8)*sc(17)
    ! reaction 142: C + O2 <=> O + CO
    qf(142) = sc(4)*sc(9)
    qr(142) = sc(3)*sc(15)
    ! reaction 143: C + CH2 <=> H + C2H
    qf(143) = sc(9)*sc(11)
    qr(143) = sc(2)*sc(22)
    ! reaction 144: C + CH3 <=> H + C2H2
    qf(144) = sc(9)*sc(13)
    qr(144) = sc(2)*sc(23)
    ! reaction 145: CH + O2 <=> O + HCO
    qf(145) = sc(4)*sc(10)
    qr(145) = sc(3)*sc(17)
    ! reaction 146: CH + H2 <=> H + CH2
    qf(146) = sc(1)*sc(10)
    qr(146) = sc(2)*sc(11)
    ! reaction 147: CH + H2O <=> H + CH2O
    qf(147) = sc(6)*sc(10)
    qr(147) = sc(2)*sc(18)
    ! reaction 148: CH + CH2 <=> H + C2H2
    qf(148) = sc(10)*sc(11)
    qr(148) = sc(2)*sc(23)
    ! reaction 149: CH + CH3 <=> H + C2H3
    qf(149) = sc(10)*sc(13)
    qr(149) = sc(2)*sc(24)
    ! reaction 150: CH + CH4 <=> H + C2H4
    qf(150) = sc(10)*sc(14)
    qr(150) = sc(2)*sc(25)
    ! reaction 151: CH + CO2 <=> HCO + CO
    qf(151) = sc(10)*sc(16)
    qr(151) = sc(15)*sc(17)
    ! reaction 152: CH + CH2O <=> H + CH2CO
    qf(152) = sc(10)*sc(18)
    qr(152) = sc(2)*sc(29)
    ! reaction 153: CH + HCCO <=> CO + C2H2
    qf(153) = sc(10)*sc(28)
    qr(153) = sc(15)*sc(23)
    ! reaction 154: CH2 + O2 => OH + H + CO
    qf(154) = sc(4)*sc(11)
    qr(154) = 0.d0
    ! reaction 155: CH2 + H2 <=> H + CH3
    qf(155) = sc(1)*sc(11)
    qr(155) = sc(2)*sc(13)
    ! reaction 156: 2 CH2 <=> H2 + C2H2
    qf(156) = sc(11)*sc(11)
    qr(156) = sc(1)*sc(23)
    ! reaction 157: CH2 + CH3 <=> H + C2H4
    qf(157) = sc(11)*sc(13)
    qr(157) = sc(2)*sc(25)
    ! reaction 158: CH2 + CH4 <=> 2 CH3
    qf(158) = sc(11)*sc(14)
    qr(158) = sc(13)*sc(13)
    ! reaction 159: CH2 + HCCO <=> C2H3 + CO
    qf(159) = sc(11)*sc(28)
    qr(159) = sc(15)*sc(24)
    ! reaction 160: CH2(S) + N2 <=> CH2 + N2
    qf(160) = sc(12)*sc(48)
    qr(160) = sc(11)*sc(48)
    ! reaction 161: CH2(S) + AR <=> CH2 + AR
    qf(161) = sc(12)*sc(49)
    qr(161) = sc(11)*sc(49)
    ! reaction 162: CH2(S) + O2 <=> H + OH + CO
    qf(162) = sc(4)*sc(12)
    qr(162) = sc(2)*sc(5)*sc(15)
    ! reaction 163: CH2(S) + O2 <=> CO + H2O
    qf(163) = sc(4)*sc(12)
    qr(163) = sc(6)*sc(15)
    ! reaction 164: CH2(S) + H2 <=> CH3 + H
    qf(164) = sc(1)*sc(12)
    qr(164) = sc(2)*sc(13)
    ! reaction 165: CH2(S) + H2O <=> CH2 + H2O
    qf(165) = sc(6)*sc(12)
    qr(165) = sc(6)*sc(11)
    ! reaction 166: CH2(S) + CH3 <=> H + C2H4
    qf(166) = sc(12)*sc(13)
    qr(166) = sc(2)*sc(25)
    ! reaction 167: CH2(S) + CH4 <=> 2 CH3
    qf(167) = sc(12)*sc(14)
    qr(167) = sc(13)*sc(13)
    ! reaction 168: CH2(S) + CO <=> CH2 + CO
    qf(168) = sc(12)*sc(15)
    qr(168) = sc(11)*sc(15)
    ! reaction 169: CH2(S) + CO2 <=> CH2 + CO2
    qf(169) = sc(12)*sc(16)
    qr(169) = sc(11)*sc(16)
    ! reaction 170: CH2(S) + CO2 <=> CO + CH2O
    qf(170) = sc(12)*sc(16)
    qr(170) = sc(15)*sc(18)
    ! reaction 171: CH2(S) + C2H6 <=> CH3 + C2H5
    qf(171) = sc(12)*sc(27)
    qr(171) = sc(13)*sc(26)
    ! reaction 172: CH3 + O2 <=> O + CH3O
    qf(172) = sc(4)*sc(13)
    qr(172) = sc(3)*sc(20)
    ! reaction 173: CH3 + O2 <=> OH + CH2O
    qf(173) = sc(4)*sc(13)
    qr(173) = sc(5)*sc(18)
    ! reaction 174: CH3 + H2O2 <=> HO2 + CH4
    qf(174) = sc(8)*sc(13)
    qr(174) = sc(7)*sc(14)
    ! reaction 175: 2 CH3 <=> H + C2H5
    qf(175) = sc(13)*sc(13)
    qr(175) = sc(2)*sc(26)
    ! reaction 176: CH3 + HCO <=> CH4 + CO
    qf(176) = sc(13)*sc(17)
    qr(176) = sc(14)*sc(15)
    ! reaction 177: CH3 + CH2O <=> HCO + CH4
    qf(177) = sc(13)*sc(18)
    qr(177) = sc(14)*sc(17)
    ! reaction 178: CH3 + CH3OH <=> CH2OH + CH4
    qf(178) = sc(13)*sc(21)
    qr(178) = sc(14)*sc(19)
    ! reaction 179: CH3 + CH3OH <=> CH3O + CH4
    qf(179) = sc(13)*sc(21)
    qr(179) = sc(14)*sc(20)
    ! reaction 180: CH3 + C2H4 <=> C2H3 + CH4
    qf(180) = sc(13)*sc(25)
    qr(180) = sc(14)*sc(24)
    ! reaction 181: CH3 + C2H6 <=> C2H5 + CH4
    qf(181) = sc(13)*sc(27)
    qr(181) = sc(14)*sc(26)
    ! reaction 182: HCO + H2O <=> H + CO + H2O
    qf(182) = sc(6)*sc(17)
    qr(182) = sc(2)*sc(6)*sc(15)
    ! reaction 183: HCO + O2 <=> HO2 + CO
    qf(183) = sc(4)*sc(17)
    qr(183) = sc(7)*sc(15)
    ! reaction 184: CH2OH + O2 <=> HO2 + CH2O
    qf(184) = sc(4)*sc(19)
    qr(184) = sc(7)*sc(18)
    ! reaction 185: CH3O + O2 <=> HO2 + CH2O
    qf(185) = sc(4)*sc(20)
    qr(185) = sc(7)*sc(18)
    ! reaction 186: C2H + O2 <=> HCO + CO
    qf(186) = sc(4)*sc(22)
    qr(186) = sc(15)*sc(17)
    ! reaction 187: C2H + H2 <=> H + C2H2
    qf(187) = sc(1)*sc(22)
    qr(187) = sc(2)*sc(23)
    ! reaction 188: C2H3 + O2 <=> HCO + CH2O
    qf(188) = sc(4)*sc(24)
    qr(188) = sc(17)*sc(18)
    ! reaction 189: C2H5 + O2 <=> HO2 + C2H4
    qf(189) = sc(4)*sc(26)
    qr(189) = sc(7)*sc(25)
    ! reaction 190: HCCO + O2 <=> OH + 2 CO
    qf(190) = sc(4)*sc(28)
    qr(190) = sc(5)*sc(15)*sc(15)
    ! reaction 191: 2 HCCO <=> 2 CO + C2H2
    qf(191) = sc(28)*sc(28)
    qr(191) = sc(15)*sc(15)*sc(23)
    ! reaction 192: N + NO <=> N2 + O
    qf(192) = sc(31)*sc(36)
    qr(192) = sc(3)*sc(48)
    ! reaction 193: N + O2 <=> NO + O
    qf(193) = sc(4)*sc(31)
    qr(193) = sc(3)*sc(36)
    ! reaction 194: N + OH <=> NO + H
    qf(194) = sc(5)*sc(31)
    qr(194) = sc(2)*sc(36)
    ! reaction 195: N2O + O <=> N2 + O2
    qf(195) = sc(3)*sc(38)
    qr(195) = sc(4)*sc(48)
    ! reaction 196: N2O + O <=> 2 NO
    qf(196) = sc(3)*sc(38)
    qr(196) = sc(36)*sc(36)
    ! reaction 197: N2O + H <=> N2 + OH
    qf(197) = sc(2)*sc(38)
    qr(197) = sc(5)*sc(48)
    ! reaction 198: N2O + OH <=> N2 + HO2
    qf(198) = sc(5)*sc(38)
    qr(198) = sc(7)*sc(48)
    ! reaction 199: HO2 + NO <=> NO2 + OH
    qf(199) = sc(7)*sc(36)
    qr(199) = sc(5)*sc(37)
    ! reaction 200: NO2 + O <=> NO + O2
    qf(200) = sc(3)*sc(37)
    qr(200) = sc(4)*sc(36)
    ! reaction 201: NO2 + H <=> NO + OH
    qf(201) = sc(2)*sc(37)
    qr(201) = sc(5)*sc(36)
    ! reaction 202: NH + O <=> NO + H
    qf(202) = sc(3)*sc(32)
    qr(202) = sc(2)*sc(36)
    ! reaction 203: NH + H <=> N + H2
    qf(203) = sc(2)*sc(32)
    qr(203) = sc(1)*sc(31)
    ! reaction 204: NH + OH <=> HNO + H
    qf(204) = sc(5)*sc(32)
    qr(204) = sc(2)*sc(39)
    ! reaction 205: NH + OH <=> N + H2O
    qf(205) = sc(5)*sc(32)
    qr(205) = sc(6)*sc(31)
    ! reaction 206: NH + O2 <=> HNO + O
    qf(206) = sc(4)*sc(32)
    qr(206) = sc(3)*sc(39)
    ! reaction 207: NH + O2 <=> NO + OH
    qf(207) = sc(4)*sc(32)
    qr(207) = sc(5)*sc(36)
    ! reaction 208: NH + N <=> N2 + H
    qf(208) = sc(31)*sc(32)
    qr(208) = sc(2)*sc(48)
    ! reaction 209: NH + H2O <=> HNO + H2
    qf(209) = sc(6)*sc(32)
    qr(209) = sc(1)*sc(39)
    ! reaction 210: NH + NO <=> N2 + OH
    qf(210) = sc(32)*sc(36)
    qr(210) = sc(5)*sc(48)
    ! reaction 211: NH + NO <=> N2O + H
    qf(211) = sc(32)*sc(36)
    qr(211) = sc(2)*sc(38)
    ! reaction 212: NH2 + O <=> OH + NH
    qf(212) = sc(3)*sc(33)
    qr(212) = sc(5)*sc(32)
    ! reaction 213: NH2 + O <=> H + HNO
    qf(213) = sc(3)*sc(33)
    qr(213) = sc(2)*sc(39)
    ! reaction 214: NH2 + H <=> NH + H2
    qf(214) = sc(2)*sc(33)
    qr(214) = sc(1)*sc(32)
    ! reaction 215: NH2 + OH <=> NH + H2O
    qf(215) = sc(5)*sc(33)
    qr(215) = sc(6)*sc(32)
    ! reaction 216: NNH <=> N2 + H
    qf(216) = sc(35)
    qr(216) = sc(2)*sc(48)
    ! reaction 217: NNH + O2 <=> HO2 + N2
    qf(217) = sc(4)*sc(35)
    qr(217) = sc(7)*sc(48)
    ! reaction 218: NNH + O <=> OH + N2
    qf(218) = sc(3)*sc(35)
    qr(218) = sc(5)*sc(48)
    ! reaction 219: NNH + O <=> NH + NO
    qf(219) = sc(3)*sc(35)
    qr(219) = sc(32)*sc(36)
    ! reaction 220: NNH + H <=> H2 + N2
    qf(220) = sc(2)*sc(35)
    qr(220) = sc(1)*sc(48)
    ! reaction 221: NNH + OH <=> H2O + N2
    qf(221) = sc(5)*sc(35)
    qr(221) = sc(6)*sc(48)
    ! reaction 222: NNH + CH3 <=> CH4 + N2
    qf(222) = sc(13)*sc(35)
    qr(222) = sc(14)*sc(48)
    ! reaction 223: HNO + O <=> NO + OH
    qf(223) = sc(3)*sc(39)
    qr(223) = sc(5)*sc(36)
    ! reaction 224: HNO + H <=> H2 + NO
    qf(224) = sc(2)*sc(39)
    qr(224) = sc(1)*sc(36)
    ! reaction 225: HNO + OH <=> NO + H2O
    qf(225) = sc(5)*sc(39)
    qr(225) = sc(6)*sc(36)
    ! reaction 226: HNO + O2 <=> HO2 + NO
    qf(226) = sc(4)*sc(39)
    qr(226) = sc(7)*sc(36)
    ! reaction 227: CN + O <=> CO + N
    qf(227) = sc(3)*sc(40)
    qr(227) = sc(15)*sc(31)
    ! reaction 228: CN + OH <=> NCO + H
    qf(228) = sc(5)*sc(40)
    qr(228) = sc(2)*sc(47)
    ! reaction 229: CN + H2O <=> HCN + OH
    qf(229) = sc(6)*sc(40)
    qr(229) = sc(5)*sc(41)
    ! reaction 230: CN + O2 <=> NCO + O
    qf(230) = sc(4)*sc(40)
    qr(230) = sc(3)*sc(47)
    ! reaction 231: CN + H2 <=> HCN + H
    qf(231) = sc(1)*sc(40)
    qr(231) = sc(2)*sc(41)
    ! reaction 232: NCO + O <=> NO + CO
    qf(232) = sc(3)*sc(47)
    qr(232) = sc(15)*sc(36)
    ! reaction 233: NCO + H <=> NH + CO
    qf(233) = sc(2)*sc(47)
    qr(233) = sc(15)*sc(32)
    ! reaction 234: NCO + OH <=> NO + H + CO
    qf(234) = sc(5)*sc(47)
    qr(234) = sc(2)*sc(15)*sc(36)
    ! reaction 235: NCO + N <=> N2 + CO
    qf(235) = sc(31)*sc(47)
    qr(235) = sc(15)*sc(48)
    ! reaction 236: NCO + O2 <=> NO + CO2
    qf(236) = sc(4)*sc(47)
    qr(236) = sc(16)*sc(36)
    ! reaction 237: NCO + NO <=> N2O + CO
    qf(237) = sc(36)*sc(47)
    qr(237) = sc(15)*sc(38)
    ! reaction 238: NCO + NO <=> N2 + CO2
    qf(238) = sc(36)*sc(47)
    qr(238) = sc(16)*sc(48)
    ! reaction 239: HCN + O <=> NCO + H
    qf(239) = sc(3)*sc(41)
    qr(239) = sc(2)*sc(47)
    ! reaction 240: HCN + O <=> NH + CO
    qf(240) = sc(3)*sc(41)
    qr(240) = sc(15)*sc(32)
    ! reaction 241: HCN + O <=> CN + OH
    qf(241) = sc(3)*sc(41)
    qr(241) = sc(5)*sc(40)
    ! reaction 242: HCN + OH <=> HOCN + H
    qf(242) = sc(5)*sc(41)
    qr(242) = sc(2)*sc(45)
    ! reaction 243: HCN + OH <=> HNCO + H
    qf(243) = sc(5)*sc(41)
    qr(243) = sc(2)*sc(46)
    ! reaction 244: HCN + OH <=> NH2 + CO
    qf(244) = sc(5)*sc(41)
    qr(244) = sc(15)*sc(33)
    ! reaction 245: H2CN + N <=> N2 + CH2
    qf(245) = sc(31)*sc(42)
    qr(245) = sc(11)*sc(48)
    ! reaction 246: C + N2 <=> CN + N
    qf(246) = sc(9)*sc(48)
    qr(246) = sc(31)*sc(40)
    ! reaction 247: CH + N2 <=> HCN + N
    qf(247) = sc(10)*sc(48)
    qr(247) = sc(31)*sc(41)
    ! reaction 248: CH2 + N2 <=> HCN + NH
    qf(248) = sc(11)*sc(48)
    qr(248) = sc(32)*sc(41)
    ! reaction 249: CH2(S) + N2 <=> NH + HCN
    qf(249) = sc(12)*sc(48)
    qr(249) = sc(32)*sc(41)
    ! reaction 250: C + NO <=> CN + O
    qf(250) = sc(9)*sc(36)
    qr(250) = sc(3)*sc(40)
    ! reaction 251: C + NO <=> CO + N
    qf(251) = sc(9)*sc(36)
    qr(251) = sc(15)*sc(31)
    ! reaction 252: CH + NO <=> HCN + O
    qf(252) = sc(10)*sc(36)
    qr(252) = sc(3)*sc(41)
    ! reaction 253: CH + NO <=> H + NCO
    qf(253) = sc(10)*sc(36)
    qr(253) = sc(2)*sc(47)
    ! reaction 254: CH + NO <=> N + HCO
    qf(254) = sc(10)*sc(36)
    qr(254) = sc(17)*sc(31)
    ! reaction 255: CH2 + NO <=> H + HNCO
    qf(255) = sc(11)*sc(36)
    qr(255) = sc(2)*sc(46)
    ! reaction 256: CH2 + NO <=> OH + HCN
    qf(256) = sc(11)*sc(36)
    qr(256) = sc(5)*sc(41)
    ! reaction 257: CH2 + NO <=> H + HCNO
    qf(257) = sc(11)*sc(36)
    qr(257) = sc(2)*sc(44)
    ! reaction 258: CH2(S) + NO <=> H + HNCO
    qf(258) = sc(12)*sc(36)
    qr(258) = sc(2)*sc(46)
    ! reaction 259: CH2(S) + NO <=> OH + HCN
    qf(259) = sc(12)*sc(36)
    qr(259) = sc(5)*sc(41)
    ! reaction 260: CH2(S) + NO <=> H + HCNO
    qf(260) = sc(12)*sc(36)
    qr(260) = sc(2)*sc(44)
    ! reaction 261: CH3 + NO <=> HCN + H2O
    qf(261) = sc(13)*sc(36)
    qr(261) = sc(6)*sc(41)
    ! reaction 262: CH3 + NO <=> H2CN + OH
    qf(262) = sc(13)*sc(36)
    qr(262) = sc(5)*sc(42)
    ! reaction 263: HCNN + O <=> CO + H + N2
    qf(263) = sc(3)*sc(43)
    qr(263) = sc(2)*sc(15)*sc(48)
    ! reaction 264: HCNN + O <=> HCN + NO
    qf(264) = sc(3)*sc(43)
    qr(264) = sc(36)*sc(41)
    ! reaction 265: HCNN + O2 <=> O + HCO + N2
    qf(265) = sc(4)*sc(43)
    qr(265) = sc(3)*sc(17)*sc(48)
    ! reaction 266: HCNN + OH <=> H + HCO + N2
    qf(266) = sc(5)*sc(43)
    qr(266) = sc(2)*sc(17)*sc(48)
    ! reaction 267: HCNN + H <=> CH2 + N2
    qf(267) = sc(2)*sc(43)
    qr(267) = sc(11)*sc(48)
    ! reaction 268: HNCO + O <=> NH + CO2
    qf(268) = sc(3)*sc(46)
    qr(268) = sc(16)*sc(32)
    ! reaction 269: HNCO + O <=> HNO + CO
    qf(269) = sc(3)*sc(46)
    qr(269) = sc(15)*sc(39)
    ! reaction 270: HNCO + O <=> NCO + OH
    qf(270) = sc(3)*sc(46)
    qr(270) = sc(5)*sc(47)
    ! reaction 271: HNCO + H <=> NH2 + CO
    qf(271) = sc(2)*sc(46)
    qr(271) = sc(15)*sc(33)
    ! reaction 272: HNCO + H <=> H2 + NCO
    qf(272) = sc(2)*sc(46)
    qr(272) = sc(1)*sc(47)
    ! reaction 273: HNCO + OH <=> NCO + H2O
    qf(273) = sc(5)*sc(46)
    qr(273) = sc(6)*sc(47)
    ! reaction 274: HNCO + OH <=> NH2 + CO2
    qf(274) = sc(5)*sc(46)
    qr(274) = sc(16)*sc(33)
    ! reaction 275: HCNO + H <=> H + HNCO
    qf(275) = sc(2)*sc(44)
    qr(275) = sc(2)*sc(46)
    ! reaction 276: HCNO + H <=> OH + HCN
    qf(276) = sc(2)*sc(44)
    qr(276) = sc(5)*sc(41)
    ! reaction 277: HCNO + H <=> NH2 + CO
    qf(277) = sc(2)*sc(44)
    qr(277) = sc(15)*sc(33)
    ! reaction 278: HOCN + H <=> H + HNCO
    qf(278) = sc(2)*sc(45)
    qr(278) = sc(2)*sc(46)
    ! reaction 279: HCCO + NO <=> HCNO + CO
    qf(279) = sc(28)*sc(36)
    qr(279) = sc(15)*sc(44)
    ! reaction 280: CH3 + N <=> H2CN + H
    qf(280) = sc(13)*sc(31)
    qr(280) = sc(2)*sc(42)
    ! reaction 281: CH3 + N <=> HCN + H2
    qf(281) = sc(13)*sc(31)
    qr(281) = sc(1)*sc(41)
    ! reaction 282: NH3 + H <=> NH2 + H2
    qf(282) = sc(2)*sc(34)
    qr(282) = sc(1)*sc(33)
    ! reaction 283: NH3 + OH <=> NH2 + H2O
    qf(283) = sc(5)*sc(34)
    qr(283) = sc(6)*sc(33)
    ! reaction 284: NH3 + O <=> NH2 + OH
    qf(284) = sc(3)*sc(34)
    qr(284) = sc(5)*sc(33)
    ! reaction 285: NH + CO2 <=> HNO + CO
    qf(285) = sc(16)*sc(32)
    qr(285) = sc(15)*sc(39)
    ! reaction 286: CN + NO2 <=> NCO + NO
    qf(286) = sc(37)*sc(40)
    qr(286) = sc(36)*sc(47)
    ! reaction 287: NCO + NO2 <=> N2O + CO2
    qf(287) = sc(37)*sc(47)
    qr(287) = sc(16)*sc(38)
    ! reaction 288: N + CO2 <=> NO + CO
    qf(288) = sc(16)*sc(31)
    qr(288) = sc(15)*sc(36)
    ! reaction 289: O + CH3 => H + H2 + CO
    qf(289) = sc(3)*sc(13)
    qr(289) = 0.d0
    ! reaction 290: O + C2H4 <=> H + CH2CHO
    qf(290) = sc(3)*sc(25)
    qr(290) = sc(2)*sc(52)
    ! reaction 291: O + C2H5 <=> H + CH3CHO
    qf(291) = sc(3)*sc(26)
    qr(291) = sc(2)*sc(53)
    ! reaction 292: OH + HO2 <=> O2 + H2O
    qf(292) = sc(5)*sc(7)
    qr(292) = sc(4)*sc(6)
    ! reaction 293: OH + CH3 => H2 + CH2O
    qf(293) = sc(5)*sc(13)
    qr(293) = 0.d0
    ! reaction 294: CH2 + O2 => 2 H + CO2
    qf(294) = sc(4)*sc(11)
    qr(294) = 0.d0
    ! reaction 295: CH2 + O2 <=> O + CH2O
    qf(295) = sc(4)*sc(11)
    qr(295) = sc(3)*sc(18)
    ! reaction 296: CH2 + CH2 => 2 H + C2H2
    qf(296) = sc(11)*sc(11)
    qr(296) = 0.d0
    ! reaction 297: CH2(S) + H2O => H2 + CH2O
    qf(297) = sc(6)*sc(12)
    qr(297) = 0.d0
    ! reaction 298: C2H3 + O2 <=> O + CH2CHO
    qf(298) = sc(4)*sc(24)
    qr(298) = sc(3)*sc(52)
    ! reaction 299: C2H3 + O2 <=> HO2 + C2H2
    qf(299) = sc(4)*sc(24)
    qr(299) = sc(7)*sc(23)
    ! reaction 300: O + CH3CHO <=> OH + CH2CHO
    qf(300) = sc(3)*sc(53)
    qr(300) = sc(5)*sc(52)
    ! reaction 301: O + CH3CHO => OH + CH3 + CO
    qf(301) = sc(3)*sc(53)
    qr(301) = 0.d0
    ! reaction 302: O2 + CH3CHO => HO2 + CH3 + CO
    qf(302) = sc(4)*sc(53)
    qr(302) = 0.d0
    ! reaction 303: H + CH3CHO <=> CH2CHO + H2
    qf(303) = sc(2)*sc(53)
    qr(303) = sc(1)*sc(52)
    ! reaction 304: H + CH3CHO => CH3 + H2 + CO
    qf(304) = sc(2)*sc(53)
    qr(304) = 0.d0
    ! reaction 305: OH + CH3CHO => CH3 + H2O + CO
    qf(305) = sc(5)*sc(53)
    qr(305) = 0.d0
    ! reaction 306: HO2 + CH3CHO => CH3 + H2O2 + CO
    qf(306) = sc(7)*sc(53)
    qr(306) = 0.d0
    ! reaction 307: CH3 + CH3CHO => CH3 + CH4 + CO
    qf(307) = sc(13)*sc(53)
    qr(307) = 0.d0
    ! reaction 308: O + CH2CHO => H + CH2 + CO2
    qf(308) = sc(3)*sc(52)
    qr(308) = 0.d0
    ! reaction 309: O2 + CH2CHO => OH + CO + CH2O
    qf(309) = sc(4)*sc(52)
    qr(309) = 0.d0
    ! reaction 310: O2 + CH2CHO => OH + 2 HCO
    qf(310) = sc(4)*sc(52)
    qr(310) = 0.d0
    ! reaction 311: H + CH2CHO <=> CH3 + HCO
    qf(311) = sc(2)*sc(52)
    qr(311) = sc(13)*sc(17)
    ! reaction 312: H + CH2CHO <=> CH2CO + H2
    qf(312) = sc(2)*sc(52)
    qr(312) = sc(1)*sc(29)
    ! reaction 313: OH + CH2CHO <=> H2O + CH2CO
    qf(313) = sc(5)*sc(52)
    qr(313) = sc(6)*sc(29)
    ! reaction 314: OH + CH2CHO <=> HCO + CH2OH
    qf(314) = sc(5)*sc(52)
    qr(314) = sc(17)*sc(19)
    ! reaction 315: O + C3H8 <=> OH + C3H7
    qf(315) = sc(3)*sc(51)
    qr(315) = sc(5)*sc(50)
    ! reaction 316: H + C3H8 <=> C3H7 + H2
    qf(316) = sc(2)*sc(51)
    qr(316) = sc(1)*sc(50)
    ! reaction 317: OH + C3H8 <=> C3H7 + H2O
    qf(317) = sc(5)*sc(51)
    qr(317) = sc(6)*sc(50)
    ! reaction 318: C3H7 + H2O2 <=> HO2 + C3H8
    qf(318) = sc(8)*sc(50)
    qr(318) = sc(7)*sc(51)
    ! reaction 319: CH3 + C3H8 <=> C3H7 + CH4
    qf(319) = sc(13)*sc(51)
    qr(319) = sc(14)*sc(50)
    ! reaction 320: O + C3H7 <=> C2H5 + CH2O
    qf(320) = sc(3)*sc(50)
    qr(320) = sc(18)*sc(26)
    ! reaction 321: H + C3H7 <=> CH3 + C2H5
    qf(321) = sc(2)*sc(50)
    qr(321) = sc(13)*sc(26)
    ! reaction 322: OH + C3H7 <=> C2H5 + CH2OH
    qf(322) = sc(5)*sc(50)
    qr(322) = sc(19)*sc(26)
    ! reaction 323: HO2 + C3H7 <=> O2 + C3H8
    qf(323) = sc(7)*sc(50)
    qr(323) = sc(4)*sc(51)
    ! reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O
    qf(324) = sc(7)*sc(50)
    qr(324) = 0.d0
    ! reaction 325: CH3 + C3H7 <=> 2 C2H5
    qf(325) = sc(13)*sc(50)
    qr(325) = sc(26)*sc(26)

    T = tc(2)

    ! compute the mixture concentration
    mixture = 0.d0
    do i=1, 53
        mixture = mixture + sc(i)
    end do

    do i=1, 325
        Corr(i) = 1.d0
    end do

    ! troe
    alpha_troe(1) = mixture + (TB(1) % vector(1) - 1)*sc(1) + (TB(1) % vector(2) - 1)*sc(6) + (TB(1) % vector(3) - 1)*sc(14) + (TB(1) % vector(4) - 1)*sc(15) + (TB(1) % vector(5) - 1)*sc(16) + (TB(1) % vector(6) - 1)*sc(27) + (TB(1) % vector(7) - 1)*sc(49)
    alpha_troe(2) = mixture + (TB(2) % vector(1) - 1)*sc(1) + (TB(2) % vector(2) - 1)*sc(6) + (TB(2) % vector(3) - 1)*sc(14) + (TB(2) % vector(4) - 1)*sc(15) + (TB(2) % vector(5) - 1)*sc(16) + (TB(2) % vector(6) - 1)*sc(27) + (TB(2) % vector(7) - 1)*sc(49)
    alpha_troe(3) = mixture + (TB(3) % vector(1) - 1)*sc(1) + (TB(3) % vector(2) - 1)*sc(6) + (TB(3) % vector(3) - 1)*sc(14) + (TB(3) % vector(4) - 1)*sc(15) + (TB(3) % vector(5) - 1)*sc(16) + (TB(3) % vector(6) - 1)*sc(27) + (TB(3) % vector(7) - 1)*sc(49)
    alpha_troe(4) = mixture + (TB(4) % vector(1) - 1)*sc(1) + (TB(4) % vector(2) - 1)*sc(6) + (TB(4) % vector(3) - 1)*sc(14) + (TB(4) % vector(4) - 1)*sc(15) + (TB(4) % vector(5) - 1)*sc(16) + (TB(4) % vector(6) - 1)*sc(27)
    alpha_troe(5) = mixture + (TB(5) % vector(1) - 1)*sc(1) + (TB(5) % vector(2) - 1)*sc(6) + (TB(5) % vector(3) - 1)*sc(14) + (TB(5) % vector(4) - 1)*sc(15) + (TB(5) % vector(5) - 1)*sc(16) + (TB(5) % vector(6) - 1)*sc(27)
    alpha_troe(6) = mixture + (TB(6) % vector(1) - 1)*sc(1) + (TB(6) % vector(2) - 1)*sc(6) + (TB(6) % vector(3) - 1)*sc(14) + (TB(6) % vector(4) - 1)*sc(15) + (TB(6) % vector(5) - 1)*sc(16) + (TB(6) % vector(6) - 1)*sc(27)
    alpha_troe(7) = mixture + (TB(7) % vector(1) - 1)*sc(1) + (TB(7) % vector(2) - 1)*sc(6) + (TB(7) % vector(3) - 1)*sc(14) + (TB(7) % vector(4) - 1)*sc(15) + (TB(7) % vector(5) - 1)*sc(16) + (TB(7) % vector(6) - 1)*sc(27)
    alpha_troe(8) = mixture + (TB(8) % vector(1) - 1)*sc(1) + (TB(8) % vector(2) - 1)*sc(6) + (TB(8) % vector(3) - 1)*sc(14) + (TB(8) % vector(4) - 1)*sc(15) + (TB(8) % vector(5) - 1)*sc(16) + (TB(8) % vector(6) - 1)*sc(27) + (TB(8) % vector(7) - 1)*sc(49)
    alpha_troe(9) = mixture + (TB(9) % vector(1) - 1)*sc(1) + (TB(9) % vector(2) - 1)*sc(6) + (TB(9) % vector(3) - 1)*sc(14) + (TB(9) % vector(4) - 1)*sc(15) + (TB(9) % vector(5) - 1)*sc(16) + (TB(9) % vector(6) - 1)*sc(27) + (TB(9) % vector(7) - 1)*sc(49)
    alpha_troe(10) = mixture + (TB(10) % vector(1) - 1)*sc(1) + (TB(10) % vector(2) - 1)*sc(6) + (TB(10) % vector(3) - 1)*sc(14) + (TB(10) % vector(4) - 1)*sc(15) + (TB(10) % vector(5) - 1)*sc(16) + (TB(10) % vector(6) - 1)*sc(27) + (TB(10) % vector(7) - 1)*sc(49)
    alpha_troe(11) = mixture + (TB(11) % vector(1) - 1)*sc(1) + (TB(11) % vector(2) - 1)*sc(6) + (TB(11) % vector(3) - 1)*sc(14) + (TB(11) % vector(4) - 1)*sc(15) + (TB(11) % vector(5) - 1)*sc(16) + (TB(11) % vector(6) - 1)*sc(27) + (TB(11) % vector(7) - 1)*sc(49)
    alpha_troe(12) = mixture + (TB(12) % vector(1) - 1)*sc(1) + (TB(12) % vector(2) - 1)*sc(6) + (TB(12) % vector(3) - 1)*sc(14) + (TB(12) % vector(4) - 1)*sc(15) + (TB(12) % vector(5) - 1)*sc(16) + (TB(12) % vector(6) - 1)*sc(27) + (TB(12) % vector(7) - 1)*sc(49)
    alpha_troe(13) = mixture + (TB(13) % vector(1) - 1)*sc(1) + (TB(13) % vector(2) - 1)*sc(6) + (TB(13) % vector(3) - 1)*sc(14) + (TB(13) % vector(4) - 1)*sc(15) + (TB(13) % vector(5) - 1)*sc(16) + (TB(13) % vector(6) - 1)*sc(27) + (TB(13) % vector(7) - 1)*sc(49)
    alpha_troe(14) = mixture + (TB(14) % vector(1) - 1)*sc(1) + (TB(14) % vector(2) - 1)*sc(6) + (TB(14) % vector(3) - 1)*sc(14) + (TB(14) % vector(4) - 1)*sc(15) + (TB(14) % vector(5) - 1)*sc(16) + (TB(14) % vector(6) - 1)*sc(27) + (TB(14) % vector(7) - 1)*sc(49)
    alpha_troe(15) = mixture + (TB(15) % vector(1) - 1)*sc(1) + (TB(15) % vector(2) - 1)*sc(6) + (TB(15) % vector(3) - 1)*sc(14) + (TB(15) % vector(4) - 1)*sc(15) + (TB(15) % vector(5) - 1)*sc(16) + (TB(15) % vector(6) - 1)*sc(27)
    alpha_troe(16) = mixture + (TB(16) % vector(1) - 1)*sc(1) + (TB(16) % vector(2) - 1)*sc(6) + (TB(16) % vector(3) - 1)*sc(14) + (TB(16) % vector(4) - 1)*sc(15) + (TB(16) % vector(5) - 1)*sc(16) + (TB(16) % vector(6) - 1)*sc(27) + (TB(16) % vector(7) - 1)*sc(49)
    alpha_troe(17) = mixture + (TB(17) % vector(1) - 1)*sc(1) + (TB(17) % vector(2) - 1)*sc(6) + (TB(17) % vector(3) - 1)*sc(14) + (TB(17) % vector(4) - 1)*sc(15) + (TB(17) % vector(5) - 1)*sc(16) + (TB(17) % vector(6) - 1)*sc(27) + (TB(17) % vector(7) - 1)*sc(49)
    alpha_troe(18) = mixture + (TB(18) % vector(1) - 1)*sc(1) + (TB(18) % vector(2) - 1)*sc(6) + (TB(18) % vector(3) - 1)*sc(14) + (TB(18) % vector(4) - 1)*sc(15) + (TB(18) % vector(5) - 1)*sc(16) + (TB(18) % vector(6) - 1)*sc(27)
    alpha_troe(19) = mixture + (TB(19) % vector(1) - 1)*sc(1) + (TB(19) % vector(2) - 1)*sc(6) + (TB(19) % vector(3) - 1)*sc(14) + (TB(19) % vector(4) - 1)*sc(15) + (TB(19) % vector(5) - 1)*sc(16) + (TB(19) % vector(6) - 1)*sc(27) + (TB(19) % vector(7) - 1)*sc(49)
    alpha_troe(20) = mixture + (TB(20) % vector(1) - 1)*sc(1) + (TB(20) % vector(2) - 1)*sc(6) + (TB(20) % vector(3) - 1)*sc(14) + (TB(20) % vector(4) - 1)*sc(15) + (TB(20) % vector(5) - 1)*sc(16) + (TB(20) % vector(6) - 1)*sc(27) + (TB(20) % vector(7) - 1)*sc(49)
    alpha_troe(21) = mixture + (TB(21) % vector(1) - 1)*sc(1) + (TB(21) % vector(2) - 1)*sc(6) + (TB(21) % vector(3) - 1)*sc(14) + (TB(21) % vector(4) - 1)*sc(15) + (TB(21) % vector(5) - 1)*sc(16) + (TB(21) % vector(6) - 1)*sc(27) + (TB(21) % vector(7) - 1)*sc(49)
    alpha_troe(22) = mixture + (TB(22) % vector(1) - 1)*sc(1) + (TB(22) % vector(2) - 1)*sc(6) + (TB(22) % vector(3) - 1)*sc(14) + (TB(22) % vector(4) - 1)*sc(15) + (TB(22) % vector(5) - 1)*sc(16) + (TB(22) % vector(6) - 1)*sc(27) + (TB(22) % vector(7) - 1)*sc(49)
    alpha_troe(23) = mixture + (TB(23) % vector(1) - 1)*sc(1) + (TB(23) % vector(2) - 1)*sc(6) + (TB(23) % vector(3) - 1)*sc(14) + (TB(23) % vector(4) - 1)*sc(15) + (TB(23) % vector(5) - 1)*sc(16) + (TB(23) % vector(6) - 1)*sc(27) + (TB(23) % vector(7) - 1)*sc(49)
    alpha_troe(24) = mixture + (TB(24) % vector(1) - 1)*sc(1) + (TB(24) % vector(2) - 1)*sc(6) + (TB(24) % vector(3) - 1)*sc(14) + (TB(24) % vector(4) - 1)*sc(15) + (TB(24) % vector(5) - 1)*sc(16) + (TB(24) % vector(6) - 1)*sc(27) + (TB(24) % vector(7) - 1)*sc(49)
    alpha_troe(25) = mixture + (TB(25) % vector(1) - 1)*sc(1) + (TB(25) % vector(2) - 1)*sc(6) + (TB(25) % vector(3) - 1)*sc(14) + (TB(25) % vector(4) - 1)*sc(15) + (TB(25) % vector(5) - 1)*sc(16) + (TB(25) % vector(6) - 1)*sc(27) + (TB(25) % vector(7) - 1)*sc(49)
    alpha_troe(26) = mixture + (TB(26) % vector(1) - 1)*sc(1) + (TB(26) % vector(2) - 1)*sc(6) + (TB(26) % vector(3) - 1)*sc(14) + (TB(26) % vector(4) - 1)*sc(15) + (TB(26) % vector(5) - 1)*sc(16) + (TB(26) % vector(6) - 1)*sc(27) + (TB(26) % vector(7) - 1)*sc(49)

    do i=1, 26
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

    ! Lindemann
    alpha_lindemann(1) = mixture + (TB(27) % vector(1) - 1)*sc(1) + (TB(27) % vector(2) - 1)*sc(4) + (TB(27) % vector(3) - 1)*sc(6) + (TB(27) % vector(4) - 1)*sc(14) + (TB(27) % vector(5) - 1)*sc(15) + (TB(27) % vector(6) - 1)*sc(16) + (TB(27) % vector(7) - 1)*sc(27) + (TB(27) % vector(8) - 1)*sc(49)
    alpha_lindemann(2) = mixture + (TB(28) % vector(1) - 1)*sc(1) + (TB(28) % vector(2) - 1)*sc(6) + (TB(28) % vector(3) - 1)*sc(14) + (TB(28) % vector(4) - 1)*sc(15) + (TB(28) % vector(5) - 1)*sc(16) + (TB(28) % vector(6) - 1)*sc(27) + (TB(28) % vector(7) - 1)*sc(49)
    alpha_lindemann(3) = mixture + (TB(29) % vector(1) - 1)*sc(1) + (TB(29) % vector(2) - 1)*sc(6) + (TB(29) % vector(3) - 1)*sc(14) + (TB(29) % vector(4) - 1)*sc(15) + (TB(29) % vector(5) - 1)*sc(16) + (TB(29) % vector(6) - 1)*sc(27) + (TB(29) % vector(7) - 1)*sc(49)
    do i=27, 29

        redP = alpha_lindemann(i-26) / k_f_save(i) * phase_units(i) * low_A(i) * exp(low_beta(i) * tc(1) - activation_units(i) * low_Ea(i) * invT)
        Corr(i) = redP / (1.d0 + redP)
    end do

! simple three-body correction
Corr(30) = mixture + (TB(30) % vector(1) - 1)*sc(1) + (TB(30) % vector(2) - 1)*sc(6) + (TB(30) % vector(3) - 1)*sc(14) + (TB(30) % vector(4) - 1)*sc(15) + (TB(30) % vector(5) - 1)*sc(16) + (TB(30) % vector(6) - 1)*sc(27) + (TB(30) % vector(7) - 1)*sc(49)
Corr(31) = mixture + (TB(31) % vector(1) - 1)*sc(1) + (TB(31) % vector(2) - 1)*sc(6) + (TB(31) % vector(3) - 1)*sc(14) + (TB(31) % vector(4) - 1)*sc(15) + (TB(31) % vector(5) - 1)*sc(16) + (TB(31) % vector(6) - 1)*sc(27) + (TB(31) % vector(7) - 1)*sc(49)
Corr(32) = mixture + (TB(32) % vector(1) - 1)*sc(4) + (TB(32) % vector(2) - 1)*sc(6) + (TB(32) % vector(3) - 1)*sc(15) + (TB(32) % vector(4) - 1)*sc(16) + (TB(32) % vector(5) - 1)*sc(27) + (TB(32) % vector(6) - 1)*sc(48) + (TB(32) % vector(7) - 1)*sc(49)
Corr(33) = mixture + (TB(33) % vector(1) - 1)*sc(1) + (TB(33) % vector(2) - 1)*sc(6) + (TB(33) % vector(3) - 1)*sc(14) + (TB(33) % vector(4) - 1)*sc(16) + (TB(33) % vector(5) - 1)*sc(27) + (TB(33) % vector(6) - 1)*sc(49)
Corr(34) = mixture + (TB(34) % vector(1) - 1)*sc(1) + (TB(34) % vector(2) - 1)*sc(6) + (TB(34) % vector(3) - 1)*sc(14) + (TB(34) % vector(4) - 1)*sc(27) + (TB(34) % vector(5) - 1)*sc(49)
Corr(35) = mixture + (TB(35) % vector(1) - 1)*sc(1) + (TB(35) % vector(2) - 1)*sc(6) + (TB(35) % vector(3) - 1)*sc(14) + (TB(35) % vector(4) - 1)*sc(15) + (TB(35) % vector(5) - 1)*sc(16) + (TB(35) % vector(6) - 1)*sc(27)
Corr(36) = mixture + (TB(36) % vector(1) - 1)*sc(1) + (TB(36) % vector(2) - 1)*sc(6) + (TB(36) % vector(3) - 1)*sc(14) + (TB(36) % vector(4) - 1)*sc(15) + (TB(36) % vector(5) - 1)*sc(16) + (TB(36) % vector(6) - 1)*sc(27) + (TB(36) % vector(7) - 1)*sc(49)
Corr(37) = mixture + (TB(37) % vector(1) - 1)*sc(1) + (TB(37) % vector(2) - 1)*sc(6) + (TB(37) % vector(3) - 1)*sc(14) + (TB(37) % vector(4) - 1)*sc(15) + (TB(37) % vector(5) - 1)*sc(16) + (TB(37) % vector(6) - 1)*sc(27) + (TB(37) % vector(7) - 1)*sc(49)
Corr(38) = mixture + (TB(38) % vector(1) - 1)*sc(1) + (TB(38) % vector(2) - 1)*sc(6) + (TB(38) % vector(3) - 1)*sc(14) + (TB(38) % vector(4) - 1)*sc(15) + (TB(38) % vector(5) - 1)*sc(16) + (TB(38) % vector(6) - 1)*sc(27) + (TB(38) % vector(7) - 1)*sc(49)
Corr(39) = mixture + (TB(39) % vector(1) - 1)*sc(1) + (TB(39) % vector(2) - 1)*sc(6) + (TB(39) % vector(3) - 1)*sc(14) + (TB(39) % vector(4) - 1)*sc(15) + (TB(39) % vector(5) - 1)*sc(16) + (TB(39) % vector(6) - 1)*sc(27) + (TB(39) % vector(7) - 1)*sc(49)
Corr(40) = mixture + (TB(40) % vector(1) - 1)*sc(1) + (TB(40) % vector(2) - 1)*sc(6) + (TB(40) % vector(3) - 1)*sc(14) + (TB(40) % vector(4) - 1)*sc(15) + (TB(40) % vector(5) - 1)*sc(16) + (TB(40) % vector(6) - 1)*sc(27) + (TB(40) % vector(7) - 1)*sc(49)
Corr(41) = mixture + (TB(41) % vector(1) - 1)*sc(1) + (TB(41) % vector(2) - 1)*sc(6) + (TB(41) % vector(3) - 1)*sc(14) + (TB(41) % vector(4) - 1)*sc(15) + (TB(41) % vector(5) - 1)*sc(16) + (TB(41) % vector(6) - 1)*sc(27) + (TB(41) % vector(7) - 1)*sc(49)

do i=1, 325
    qf(i) = qf(i) * (Corr(i) * k_f_save(i))
    qr(i) = qr(i) * (Corr(i) * k_f_save(i) / Kc_save(i))
end do

end subroutine

! compute the g/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine gibbs(species, tc)

implicit none

double precision, intent(out) :: species(53)
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
    ! species 8: H2O2
    species(8) = &
        -1.770258210000000d+04 * invT &
        +8.410619499999998d-01 &
        -4.276112690000000d+00 * tc(1) &
        +2.714112085000000d-04 * tc(2) &
        -2.788928350000000d-06 * tc(3) &
        +1.798090108333333d-09 * tc(4) &
        -4.312271815000000d-13 * tc(5)
    ! species 9: C
    species(9) = &
        +8.544388320000000d+04 * invT &
        -1.977068930000000d+00 &
        -2.554239550000000d+00 * tc(1) &
        +1.607688620000000d-04 * tc(2) &
        -1.222987075000000d-07 * tc(3) &
        +6.101957408333333d-11 * tc(4) &
        -1.332607230000000d-14 * tc(5)
    ! species 10: CH
    species(10) = &
        +7.079729340000000d+04 * invT &
        +1.405805570000000d+00 &
        -3.489816650000000d+00 * tc(1) &
        -1.619177705000000d-04 * tc(2) &
        +2.814984416666667d-07 * tc(3) &
        -2.635144391666666d-10 * tc(4) &
        +7.030453350000001d-14 * tc(5)
    ! species 11: CH2
    species(11) = &
        +4.600404010000000d+04 * invT &
        +2.200146820000000d+00 &
        -3.762678670000000d+00 * tc(1) &
        -4.844360715000000d-04 * tc(2) &
        -4.658164016666667d-07 * tc(3) &
        +3.209092941666667d-10 * tc(4) &
        -8.437085950000000d-14 * tc(5)
    ! species 12: CH2(S)
    species(12) = &
        +5.049681630000000d+04 * invT &
        +4.967723077000000d+00 &
        -4.198604110000000d+00 * tc(1) &
        +1.183307095000000d-03 * tc(2) &
        -1.372160366666667d-06 * tc(3) &
        +5.573466508333334d-10 * tc(4) &
        -9.715736850000000d-14 * tc(5)
    ! species 13: CH3
    species(13) = &
        +1.644499880000000d+04 * invT &
        +2.069026070000000d+00 &
        -3.673590400000000d+00 * tc(1) &
        -1.005475875000000d-03 * tc(2) &
        -9.550364266666668d-07 * tc(3) &
        +5.725978541666666d-10 * tc(4) &
        -1.271928670000000d-13 * tc(5)
    ! species 14: CH4
    species(14) = &
        -1.024664760000000d+04 * invT &
        +9.791179889999999d+00 &
        -5.149876130000000d+00 * tc(1) &
        +6.835489400000000d-03 * tc(2) &
        -8.196676650000000d-06 * tc(3) &
        +4.039525216666667d-09 * tc(4) &
        -8.334697800000000d-13 * tc(5)
    ! species 15: CO
    species(15) = &
        -1.434408600000000d+04 * invT &
        +7.112418999999992d-02 &
        -3.579533470000000d+00 * tc(1) &
        +3.051768400000000d-04 * tc(2) &
        -1.694690550000000d-07 * tc(3) &
        -7.558382366666667d-11 * tc(4) &
        +4.522122495000000d-14 * tc(5)
    ! species 16: CO2
    species(16) = &
        -4.837196970000000d+04 * invT &
        -7.544278700000000d+00 &
        -2.356773520000000d+00 * tc(1) &
        -4.492298385000000d-03 * tc(2) &
        +1.187260448333333d-06 * tc(3) &
        -2.049325183333333d-10 * tc(4) &
        +7.184977399999999d-15 * tc(5)
    ! species 17: HCO
    species(17) = &
        +3.839564960000000d+03 * invT &
        +8.268134100000002d-01 &
        -4.221185840000000d+00 * tc(1) &
        +1.621962660000000d-03 * tc(2) &
        -2.296657433333333d-06 * tc(3) &
        +1.109534108333333d-09 * tc(4) &
        -2.168844325000000d-13 * tc(5)
    ! species 18: CH2O
    species(18) = &
        -1.430895670000000d+04 * invT &
        +4.190910250000000d+00 &
        -4.793723150000000d+00 * tc(1) &
        +4.954166845000000d-03 * tc(2) &
        -6.220333466666666d-06 * tc(3) &
        +3.160710508333333d-09 * tc(4) &
        -6.588632600000000d-13 * tc(5)
    ! species 19: CH2OH
    species(19) = &
        -3.193913670000000d+03 * invT &
        -1.609133250000000d+00 &
        -3.863889180000000d+00 * tc(1) &
        -2.798361520000000d-03 * tc(2) &
        -9.887863183333334d-07 * tc(3) &
        +8.711001000000001d-10 * tc(4) &
        -2.184836390000000d-13 * tc(5)
    ! species 20: CH3O
    species(20) = &
        +9.786011000000000d+02 * invT &
        -1.104597300000000d+01 &
        -2.106204000000000d+00 * tc(1) &
        -3.608297500000000d-03 * tc(2) &
        -8.897453333333333d-07 * tc(3) &
        +6.148030000000000d-10 * tc(4) &
        -1.037805000000000d-13 * tc(5)
    ! species 21: CH3OH
    species(21) = &
        -2.564276560000000d+04 * invT &
        +7.219494050000001d+00 &
        -5.715395820000000d+00 * tc(1) &
        +7.615456450000000d-03 * tc(2) &
        -1.087401925000000d-05 * tc(3) &
        +5.923390741666667d-09 * tc(4) &
        -1.306763490000000d-12 * tc(5)
    ! species 22: C2H
    species(22) = &
        +6.683939320000001d+04 * invT &
        -3.333307050000000d+00 &
        -2.889657330000000d+00 * tc(1) &
        -6.704980550000000d-03 * tc(2) &
        +4.746158350000000d-06 * tc(3) &
        -2.456592041666667d-09 * tc(4) &
        +5.466575550000000d-13 * tc(5)
    ! species 23: C2H2
    species(23) = &
        +2.642898070000000d+04 * invT &
        -1.313102400600000d+01 &
        -8.086810940000000d-01 * tc(1) &
        -1.168078145000000d-02 * tc(2) &
        +5.919530250000000d-06 * tc(3) &
        -2.334603641666667d-09 * tc(4) &
        +4.250364870000000d-13 * tc(5)
    ! species 24: C2H3
    species(24) = &
        +3.485984680000000d+04 * invT &
        -5.298073800000000d+00 &
        -3.212466450000000d+00 * tc(1) &
        -7.573958100000000d-04 * tc(2) &
        -4.320156866666666d-06 * tc(3) &
        +2.980482058333333d-09 * tc(4) &
        -7.357543650000000d-13 * tc(5)
    ! species 25: C2H4
    species(25) = &
        +5.089775930000000d+03 * invT &
        -1.381294799999999d-01 &
        -3.959201480000000d+00 * tc(1) &
        +3.785261235000000d-03 * tc(2) &
        -9.516504866666667d-06 * tc(3) &
        +5.763239608333333d-09 * tc(4) &
        -1.349421865000000d-12 * tc(5)
    ! species 26: C2H5
    species(26) = &
        +1.284162650000000d+04 * invT &
        -4.007435600000004d-01 &
        -4.306465680000000d+00 * tc(1) &
        +2.093294460000000d-03 * tc(2) &
        -8.285713450000000d-06 * tc(3) &
        +4.992721716666666d-09 * tc(4) &
        -1.152545020000000d-12 * tc(5)
    ! species 27: C2H6
    species(27) = &
        -1.152220550000000d+04 * invT &
        +1.624601760000000d+00 &
        -4.291424920000000d+00 * tc(1) &
        +2.750771350000000d-03 * tc(2) &
        -9.990638133333334d-06 * tc(3) &
        +5.903885708333334d-09 * tc(4) &
        -1.343428855000000d-12 * tc(5)
    ! species 28: HCCO
    species(28) = &
        +2.005944900000000d+04 * invT &
        -1.023869560000000d+01 &
        -2.251721400000000d+00 * tc(1) &
        -8.827510500000000d-03 * tc(2) &
        +3.954850166666666d-06 * tc(3) &
        -1.439646583333334d-09 * tc(4) &
        +2.533240550000000d-13 * tc(5)
    ! species 29: CH2CO
    species(29) = &
        -7.042918040000000d+03 * invT &
        -1.007981170000000d+01 &
        -2.135836300000000d+00 * tc(1) &
        -9.059436050000000d-03 * tc(2) &
        +2.899124566666666d-06 * tc(3) &
        -7.786646400000000d-10 * tc(4) &
        +1.007288075000000d-13 * tc(5)
    ! species 30: HCCOH
    species(30) = &
        +8.031614300000000d+03 * invT &
        -1.263194570000000d+01 &
        -1.242373300000000d+00 * tc(1) &
        -1.553610050000000d-02 * tc(2) &
        +8.477810666666667d-06 * tc(3) &
        -3.594760916666667d-09 * tc(4) &
        +7.007297000000000d-13 * tc(5)
    ! species 31: N
    species(31) = &
        +5.610463700000000d+04 * invT &
        -1.693908700000000d+00 &
        -2.500000000000000d+00 * tc(1) &
        -0.000000000000000d+00 * tc(2) &
        -0.000000000000000d+00 * tc(3) &
        -0.000000000000000d+00 * tc(4) &
        -0.000000000000000d+00 * tc(5)
    ! species 32: NH
    species(32) = &
        +4.188062900000000d+04 * invT &
        +1.644580700000000d+00 &
        -3.492908500000000d+00 * tc(1) &
        -1.558959900000000d-04 * tc(2) &
        +2.481747333333333d-07 * tc(3) &
        -2.068036833333333d-10 * tc(4) &
        +5.178483500000000d-14 * tc(5)
    ! species 33: NH2
    species(33) = &
        +2.188591000000000d+04 * invT &
        +4.345845380000000d+00 &
        -4.204002900000000d+00 * tc(1) &
        +1.053069250000000d-03 * tc(2) &
        -1.184472466666667d-06 * tc(3) &
        +4.676266416666667d-10 * tc(4) &
        -8.220358500000000d-14 * tc(5)
    ! species 34: NH3
    species(34) = &
        -6.741728500000000d+03 * invT &
        +4.911400170000000d+00 &
        -4.286027400000000d+00 * tc(1) &
        +2.330261500000000d-03 * tc(2) &
        -3.619752166666667d-06 * tc(3) &
        +1.900740583333333d-09 * tc(4) &
        -4.131902300000000d-13 * tc(5)
    ! species 35: NNH
    species(35) = &
        +2.879197300000000d+04 * invT &
        +1.366751700000000d+00 &
        -4.344692700000000d+00 * tc(1) &
        +2.424853600000000d-03 * tc(2) &
        -3.343243166666667d-06 * tc(3) &
        +1.810538666666667d-09 * tc(4) &
        -3.973476950000000d-13 * tc(5)
    ! species 36: NO
    species(36) = &
        +9.844623000000000d+03 * invT &
        +1.937629900000000d+00 &
        -4.218476300000000d+00 * tc(1) &
        +2.319488000000000d-03 * tc(2) &
        -1.840170333333333d-06 * tc(3) &
        +7.780112833333333d-10 * tc(4) &
        -1.401788500000000d-13 * tc(5)
    ! species 37: NO2
    species(37) = &
        +2.896617900000000d+03 * invT &
        -2.367960500000000d+00 &
        -3.944031200000000d+00 * tc(1) &
        +7.927145000000000d-04 * tc(2) &
        -2.776302000000000d-06 * tc(3) &
        +1.706285500000000d-09 * tc(4) &
        -3.917528200000000d-13 * tc(5)
    ! species 38: N2O
    species(38) = &
        +8.741774400000000d+03 * invT &
        -8.500841800000000d+00 &
        -2.257150200000000d+00 * tc(1) &
        -5.652364000000000d-03 * tc(2) &
        +2.278553166666666d-06 * tc(3) &
        -8.068317166666668d-10 * tc(4) &
        +1.465359100000000d-13 * tc(5)
    ! species 39: HNO
    species(39) = &
        +1.154829700000000d+04 * invT &
        +2.783649899999999d+00 &
        -4.533491600000000d+00 * tc(1) &
        +2.834808550000000d-03 * tc(2) &
        -3.078867833333333d-06 * tc(3) &
        +1.428091166666667d-09 * tc(4) &
        -2.772728650000000d-13 * tc(5)
    ! species 40: CN
    species(40) = &
        +5.170834000000000d+04 * invT &
        -3.675644000000000d-01 &
        -3.612935100000000d+00 * tc(1) &
        +4.777566350000000d-04 * tc(2) &
        -3.573829500000000d-07 * tc(3) &
        +2.626360250000000d-11 * tc(4) &
        +2.321517800000000d-14 * tc(5)
    ! species 41: HCN
    species(41) = &
        +1.471263300000000d+04 * invT &
        -6.657453300000000d+00 &
        -2.258988600000000d+00 * tc(1) &
        -5.025585000000000d-03 * tc(2) &
        +2.225293833333333d-06 * tc(3) &
        -8.410290833333334d-10 * tc(4) &
        +1.504451400000000d-13 * tc(5)
    ! species 42: H2CN
    species(42) = &
        +2.863782000000000d+04 * invT &
        -6.141090100000000d+00 &
        -2.851661000000000d+00 * tc(1) &
        -2.847616550000000d-03 * tc(2) &
        -1.785233333333333d-07 * tc(3) &
        +1.352176666666667d-10 * tc(4) &
        +1.175554050000000d-14 * tc(5)
    ! species 43: HCNN
    species(43) = &
        +5.426198400000000d+04 * invT &
        -9.151550600000000d+00 &
        -2.524319400000000d+00 * tc(1) &
        -7.980309499999999d-03 * tc(2) &
        +3.136059000000000d-06 * tc(3) &
        -1.010461666666667d-09 * tc(4) &
        +1.617868900000000d-13 * tc(5)
    ! species 47: NCO
    species(47) = &
        +1.468247700000000d+04 * invT &
        -6.723533800000000d+00 &
        -2.826930800000000d+00 * tc(1) &
        -4.402584400000000d-03 * tc(2) &
        +1.397768900000000d-06 * tc(3) &
        -4.001413666666666d-10 * tc(4) &
        +6.656797500000001d-14 * tc(5)
    ! species 48: N2
    species(48) = &
        -1.020899900000000d+03 * invT &
        -6.516950000000001d-01 &
        -3.298677000000000d+00 * tc(1) &
        -7.041202000000000d-04 * tc(2) &
        +6.605369999999999d-07 * tc(3) &
        -4.701262500000001d-10 * tc(4) &
        +1.222427000000000d-13 * tc(5)
    ! species 49: AR
    species(49) = &
        -7.453750000000000d+02 * invT &
        -1.866000000000000d+00 &
        -2.500000000000000d+00 * tc(1) &
        -0.000000000000000d+00 * tc(2) &
        -0.000000000000000d+00 * tc(3) &
        -0.000000000000000d+00 * tc(4) &
        -0.000000000000000d+00 * tc(5)
    ! species 50: C3H7
    species(50) = &
        +1.063186300000000d+04 * invT &
        -2.007100720000000d+01 &
        -1.051551800000000d+00 * tc(1) &
        -1.299599000000000d-02 * tc(2) &
        -3.966756666666666d-07 * tc(3) &
        +1.634130750000000d-09 * tc(4) &
        -4.686623500000000d-13 * tc(5)
    ! species 51: C3H8
    species(51) = &
        -1.395852000000000d+04 * invT &
        -1.826813719000000d+01 &
        -9.335538100000000d-01 * tc(1) &
        -1.321228950000000d-02 * tc(2) &
        -1.017662116666667d-06 * tc(3) &
        +1.831458250000000d-09 * tc(4) &
        -4.757462650000000d-13 * tc(5)
    ! species 52: CH2CHO
    species(52) = &
        +1.521476600000000d+03 * invT &
        -6.149227999999999d+00 &
        -3.409062000000000d+00 * tc(1) &
        -5.369287000000000d-03 * tc(2) &
        -3.152486666666667d-07 * tc(3) &
        +5.965485833333333d-10 * tc(4) &
        -1.433692500000000d-13 * tc(5)
    ! species 53: CH3CHO
    species(53) = &
        -2.157287800000000d+04 * invT &
        +6.264436000000000d-01 &
        -4.729459500000000d+00 * tc(1) &
        +1.596642900000000d-03 * tc(2) &
        -7.922486833333334d-06 * tc(3) &
        +4.788217583333333d-09 * tc(4) &
        -1.096555600000000d-12 * tc(5)
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
    !species 8: H2O2
    species(8) = &
        -1.786178770000000d+04 * invT &
        +1.248846229999999d+00 &
        -4.165002850000000d+00 * tc(1) &
        -2.454158470000000d-03 * tc(2) &
        +3.168987083333333d-07 * tc(3) &
        -3.093216550000000d-11 * tc(4) &
        +1.439541525000000d-15 * tc(5)
    !species 9: C
    species(9) = &
        +8.545129530000000d+04 * invT &
        -2.308834850000000d+00 &
        -2.492668880000000d+00 * tc(1) &
        -2.399446420000000d-05 * tc(2) &
        +1.207225033333333d-08 * tc(3) &
        -3.119091908333333d-12 * tc(4) &
        +2.436389465000000d-16 * tc(5)
    !species 10: CH
    species(10) = &
        +7.101243640000001d+04 * invT &
        -2.606515260000000d+00 &
        -2.878464730000000d+00 * tc(1) &
        -4.854568405000000d-04 * tc(2) &
        -2.407427583333333d-08 * tc(3) &
        +1.089065408333333d-11 * tc(4) &
        -8.803969149999999d-16 * tc(5)
    !species 11: CH2
    species(11) = &
        +4.626360400000000d+04 * invT &
        -3.297092110000000d+00 &
        -2.874101130000000d+00 * tc(1) &
        -1.828196460000000d-03 * tc(2) &
        +2.348243283333333d-07 * tc(3) &
        -2.168162908333333d-11 * tc(4) &
        +9.386378350000000d-16 * tc(5)
    !species 12: CH2(S)
    species(12) = &
        +5.092599970000000d+04 * invT &
        -6.334463270000000d+00 &
        -2.292038420000000d+00 * tc(1) &
        -2.327943185000000d-03 * tc(2) &
        +3.353199116666667d-07 * tc(3) &
        -3.482550000000000d-11 * tc(4) &
        +1.698581825000000d-15 * tc(5)
    !species 13: CH3
    species(13) = &
        +1.677558430000000d+04 * invT &
        -6.194354070000000d+00 &
        -2.285717720000000d+00 * tc(1) &
        -3.619950185000000d-03 * tc(2) &
        +4.978572466666667d-07 * tc(3) &
        -4.964038700000000d-11 * tc(4) &
        +2.335771970000000d-15 * tc(5)
    !species 14: CH4
    species(14) = &
        -9.468344590000001d+03 * invT &
        -1.836246650500000d+01 &
        -7.485149500000000d-02 * tc(1) &
        -6.695473350000000d-03 * tc(2) &
        +9.554763483333333d-07 * tc(3) &
        -1.019104458333333d-10 * tc(4) &
        +5.090761500000000d-15 * tc(5)
    !species 15: CO
    species(15) = &
        -1.415187240000000d+04 * invT &
        -5.103502110000000d+00 &
        -2.715185610000000d+00 * tc(1) &
        -1.031263715000000d-03 * tc(2) &
        +1.664709618333334d-07 * tc(3) &
        -1.917108400000000d-11 * tc(4) &
        +1.018238580000000d-15 * tc(5)
    !species 16: CO2
    species(16) = &
        -4.875916600000000d+04 * invT &
        +1.585822230000000d+00 &
        -3.857460290000000d+00 * tc(1) &
        -2.207185130000000d-03 * tc(2) &
        +3.691356733333334d-07 * tc(3) &
        -4.362418233333334d-11 * tc(4) &
        +2.360420820000000d-15 * tc(5)
    !species 17: HCO
    species(17) = &
        +4.011918150000000d+03 * invT &
        -7.026170540000000d+00 &
        -2.772174380000000d+00 * tc(1) &
        -2.478477630000000d-03 * tc(2) &
        +4.140760216666667d-07 * tc(3) &
        -4.909681483333334d-11 * tc(4) &
        +2.667543555000000d-15 * tc(5)
    !species 18: CH2O
    species(18) = &
        -1.399583230000000d+04 * invT &
        -1.189563292000000d+01 &
        -1.760690080000000d+00 * tc(1) &
        -4.600000410000000d-03 * tc(2) &
        +7.370980216666666d-07 * tc(3) &
        -8.386767666666666d-11 * tc(4) &
        +4.419278200000001d-15 * tc(5)
    !species 19: CH2OH
    species(19) = &
        -3.242506270000000d+03 * invT &
        -2.117766459999999d+00 &
        -3.692665690000000d+00 * tc(1) &
        -4.322883985000000d-03 * tc(2) &
        +6.251685333333334d-07 * tc(3) &
        -6.560288633333334d-11 * tc(4) &
        +3.242771005000000d-15 * tc(5)
    !species 20: CH3O
    species(20) = &
        +1.278325200000000d+02 * invT &
        +8.412240000000000d-01 &
        -3.770799000000000d+00 * tc(1) &
        -3.935748500000000d-03 * tc(2) &
        +4.427306666666667d-07 * tc(3) &
        -3.287025833333333d-11 * tc(4) &
        +1.056308000000000d-15 * tc(5)
    !species 21: CH3OH
    species(21) = &
        -2.537487470000000d+04 * invT &
        -1.271265439000000d+01 &
        -1.789707910000000d+00 * tc(1) &
        -7.046914600000000d-03 * tc(2) &
        +1.060834725000000d-06 * tc(3) &
        -1.151425708333333d-10 * tc(4) &
        +5.853011000000000d-15 * tc(5)
    !species 22: C2H
    species(22) = &
        +6.712106500000000d+04 * invT &
        -3.468088230000000d+00 &
        -3.167806520000000d+00 * tc(1) &
        -2.376109510000000d-03 * tc(2) &
        +3.063117950000000d-07 * tc(3) &
        -2.534918766666666d-11 * tc(4) &
        +8.861638500000000d-16 * tc(5)
    !species 23: C2H2
    species(23) = &
        +2.593599920000000d+04 * invT &
        +5.377850850000001d+00 &
        -4.147569640000000d+00 * tc(1) &
        -2.980833320000000d-03 * tc(2) &
        +3.954914200000000d-07 * tc(3) &
        -3.895101425000000d-11 * tc(4) &
        +1.806176065000000d-15 * tc(5)
    !species 24: C2H3
    species(24) = &
        +3.461287390000000d+04 * invT &
        -4.770599780000000d+00 &
        -3.016724000000000d+00 * tc(1) &
        -5.165114600000000d-03 * tc(2) &
        +7.801372483333333d-07 * tc(3) &
        -8.480274000000000d-11 * tc(4) &
        +4.313035205000000d-15 * tc(5)
    !species 25: C2H4
    species(25) = &
        +4.939886140000000d+03 * invT &
        -8.269258140000002d+00 &
        -2.036111160000000d+00 * tc(1) &
        -7.322707550000000d-03 * tc(2) &
        +1.118463191666667d-06 * tc(3) &
        -1.226857691666667d-10 * tc(4) &
        +6.285303050000000d-15 * tc(5)
    !species 26: C2H5
    species(26) = &
        +1.285752000000000d+04 * invT &
        -1.150777788000000d+01 &
        -1.954656420000000d+00 * tc(1) &
        -8.698636100000001d-03 * tc(2) &
        +1.330344446666667d-06 * tc(3) &
        -1.460147408333333d-10 * tc(4) &
        +7.482078800000000d-15 * tc(5)
    !species 27: C2H6
    species(27) = &
        -1.142639320000000d+04 * invT &
        -1.404372920000000d+01 &
        -1.071881500000000d+00 * tc(1) &
        -1.084263385000000d-02 * tc(2) &
        +1.670934450000000d-06 * tc(3) &
        -1.845100008333333d-10 * tc(4) &
        +9.500144500000000d-15 * tc(5)
    !species 28: HCCO
    species(28) = &
        +1.932721500000000d+04 * invT &
        +9.558465300000000d+00 &
        -5.628205800000000d+00 * tc(1) &
        -2.042670050000000d-03 * tc(2) &
        +2.655757833333333d-07 * tc(3) &
        -2.385504333333333d-11 * tc(4) &
        +9.703915999999999d-16 * tc(5)
    !species 29: CH2CO
    species(29) = &
        -7.551053110000000d+03 * invT &
        +3.879050115000000d+00 &
        -4.511297320000000d+00 * tc(1) &
        -4.501798725000000d-03 * tc(2) &
        +6.948993916666666d-07 * tc(3) &
        -7.694549016666667d-11 * tc(4) &
        +3.974191005000000d-15 * tc(5)
    !species 30: HCCOH
    species(30) = &
        +7.264626000000000d+03 * invT &
        +1.352560330000000d+01 &
        -5.923829100000000d+00 * tc(1) &
        -3.396180000000000d-03 * tc(2) &
        +4.276427333333333d-07 * tc(3) &
        -3.748986750000000d-11 * tc(4) &
        +1.497005050000000d-15 * tc(5)
    !species 31: N
    species(31) = &
        +5.613377300000000d+04 * invT &
        -2.233666700000000d+00 &
        -2.415942900000000d+00 * tc(1) &
        -8.744532500000000d-05 * tc(2) &
        +1.983728166666667d-08 * tc(3) &
        -2.518853750000000d-12 * tc(4) &
        +1.018049100000000d-16 * tc(5)
    !species 32: NH
    species(32) = &
        +4.212084800000000d+04 * invT &
        -2.957087100000000d+00 &
        -2.783692800000000d+00 * tc(1) &
        -6.649215000000000d-04 * tc(2) &
        +7.079674500000000d-08 * tc(3) &
        -6.529041750000000d-12 * tc(4) &
        +2.752223500000000d-16 * tc(5)
    !species 33: NH2
    species(33) = &
        +2.217195700000000d+04 * invT &
        -3.685674200000000d+00 &
        -2.834742100000000d+00 * tc(1) &
        -1.603654100000000d-03 * tc(2) &
        +1.556513400000000d-07 * tc(3) &
        -1.141912750000000d-11 * tc(4) &
        +3.960307200000000d-16 * tc(5)
    !species 34: NH3
    species(34) = &
        -6.544695800000000d+03 * invT &
        -3.931840700000000d+00 &
        -2.634452100000000d+00 * tc(1) &
        -2.833128000000000d-03 * tc(2) &
        +2.879779333333333d-07 * tc(3) &
        -1.988930083333333d-11 * tc(4) &
        +6.289392999999999d-16 * tc(5)
    !species 35: NNH
    species(35) = &
        +2.865069700000000d+04 * invT &
        -7.037522999999997d-01 &
        -3.766754400000000d+00 * tc(1) &
        -1.445754100000000d-03 * tc(2) &
        +1.736103333333333d-07 * tc(3) &
        -1.403549500000000d-11 * tc(4) &
        +5.045948000000000d-16 * tc(5)
    !species 36: NO
    species(36) = &
        +9.920974600000000d+03 * invT &
        -3.108697100000001d+00 &
        -3.260605600000000d+00 * tc(1) &
        -5.955521500000000d-04 * tc(2) &
        +7.152841333333333d-08 * tc(3) &
        -5.788139083333334d-12 * tc(4) &
        +2.016804950000000d-16 * tc(5)
    !species 37: NO2
    species(37) = &
        +2.316498300000000d+03 * invT &
        +5.002171150000000d+00 &
        -4.884754200000000d+00 * tc(1) &
        -1.086197800000000d-03 * tc(2) &
        +1.380115100000000d-07 * tc(3) &
        -1.312292500000000d-11 * tc(4) &
        +5.255447500000000d-16 * tc(5)
    !species 38: N2O
    species(38) = &
        +8.073404800000000d+03 * invT &
        +7.024793600000000d+00 &
        -4.823072900000000d+00 * tc(1) &
        -1.313512550000000d-03 * tc(2) &
        +1.597514566666667d-07 * tc(3) &
        -1.333392666666667d-11 * tc(4) &
        +4.887615150000000d-16 * tc(5)
    !species 39: HNO
    species(39) = &
        +1.175058200000000d+04 * invT &
        -5.627121900000001d+00 &
        -2.979250900000000d+00 * tc(1) &
        -1.747202950000000d-03 * tc(2) &
        +1.309162966666667d-07 * tc(3) &
        -4.789966166666667d-12 * tc(4) &
        +9.667958000000000d-18 * tc(5)
    !species 40: CN
    species(40) = &
        +5.153618800000000d+04 * invT &
        +9.592204000000000d-01 &
        -3.745980500000000d+00 * tc(1) &
        -2.172538750000000d-05 * tc(2) &
        -4.950997333333333d-08 * tc(3) &
        +5.720983833333333d-12 * tc(4) &
        -2.206708650000000d-16 * tc(5)
    !species 41: HCN
    species(41) = &
        +1.440729200000000d+04 * invT &
        +2.226779100000000d+00 &
        -3.802239200000000d+00 * tc(1) &
        -1.573211400000000d-03 * tc(2) &
        +1.772030833333333d-07 * tc(3) &
        -1.384979750000000d-11 * tc(4) &
        +4.899878500000000d-16 * tc(5)
    !species 42: H2CN
    species(42) = &
        +2.767710900000000d+04 * invT &
        +9.654181000000001d+00 &
        -5.209703000000000d+00 * tc(1) &
        -1.484645550000000d-03 * tc(2) &
        +4.759315166666666d-08 * tc(3) &
        +1.362958333333333d-11 * tc(4) &
        -1.521629450000000d-15 * tc(5)
    !species 43: HCNN
    species(43) = &
        +5.345294100000000d+04 * invT &
        +1.099768640000000d+01 &
        -5.894636200000000d+00 * tc(1) &
        -1.994797950000000d-03 * tc(2) &
        +2.663730000000000d-07 * tc(3) &
        -2.437449583333333d-11 * tc(4) &
        +1.004734300000000d-15 * tc(5)
    !species 47: NCO
    species(47) = &
        +1.400412300000000d+04 * invT &
        +7.696450499999999d+00 &
        -5.152184500000000d+00 * tc(1) &
        -1.152588050000000d-03 * tc(2) &
        +1.467219216666667d-07 * tc(3) &
        -1.232424833333333d-11 * tc(4) &
        +4.548899800000000d-16 * tc(5)
    !species 48: N2
    species(48) = &
        -9.227977000000000d+02 * invT &
        -3.053888000000000d+00 &
        -2.926640000000000d+00 * tc(1) &
        -7.439884000000000d-04 * tc(2) &
        +9.474600000000001d-08 * tc(3) &
        -8.414198333333333d-12 * tc(4) &
        +3.376675500000000d-16 * tc(5)
    !species 49: AR
    species(49) = &
        -7.453750000000000d+02 * invT &
        -1.866000000000000d+00 &
        -2.500000000000000d+00 * tc(1) &
        -0.000000000000000d+00 * tc(2) &
        -0.000000000000000d+00 * tc(3) &
        -0.000000000000000d+00 * tc(4) &
        -0.000000000000000d+00 * tc(5)
    !species 50: C3H7
    species(50) = &
        +8.298433600000000d+03 * invT &
        +2.318287870000000d+01 &
        -7.702698700000000d+00 * tc(1) &
        -8.022101500000000d-03 * tc(2) &
        +8.805536666666666d-07 * tc(3) &
        -6.358215833333333d-11 * tc(4) &
        +1.969614200000000d-15 * tc(5)
    !species 51: C3H8
    species(51) = &
        -1.646751600000000d+04 * invT &
        +2.542648580000000d+01 &
        -7.534136800000000d+00 * tc(1) &
        -9.436119499999999d-03 * tc(2) &
        +1.045308183333333d-06 * tc(3) &
        -7.622970750000001d-11 * tc(4) &
        +2.391903450000000d-15 * tc(5)
    !species 52: CH2CHO
    species(52) = &
        +4.903218000000000d+02 * invT &
        +1.102092100000000d+01 &
        -5.975670000000000d+00 * tc(1) &
        -4.065295500000000d-03 * tc(2) &
        +4.572706666666667d-07 * tc(3) &
        -3.391920000000000d-11 * tc(4) &
        +1.088008500000000d-15 * tc(5)
    !species 53: CH3CHO
    species(53) = &
        -2.259312200000000d+04 * invT &
        +8.884902499999999d+00 &
        -5.404110800000000d+00 * tc(1) &
        -5.861529500000000d-03 * tc(2) &
        +7.043856166666666d-07 * tc(3) &
        -5.697704250000000d-11 * tc(4) &
        +2.049243150000000d-15 * tc(5)
end if

! species with midpoint at T=1368 kelvin
if (T < 1368.000000d0) then
    ! species 45: HOCN
    species(45) = &
        -2.826984000000000d+03 * invT &
        -1.846872100000000d+00 &
        -3.786049520000000d+00 * tc(1) &
        -3.443339610000000d-03 * tc(2) &
        +5.358131066666666d-07 * tc(3) &
        -4.309964725000000d-11 * tc(4) &
        -5.968039400000000d-16 * tc(5)
else
    !species 45: HOCN
    species(45) = &
        -3.706533310000000d+03 * invT &
        +1.207952710000000d+01 &
        -5.897848850000000d+00 * tc(1) &
        -1.583946965000000d-03 * tc(2) &
        +1.863351066666667d-07 * tc(3) &
        -1.477026200000000d-11 * tc(4) &
        +5.216958849999999d-16 * tc(5)
end if

! species with midpoint at T=1478 kelvin
if (T < 1478.000000d0) then
    ! species 46: HNCO
    species(46) = &
        -1.558736360000000d+04 * invT &
        -2.563614100000000d+00 &
        -3.630963170000000d+00 * tc(1) &
        -3.651411785000000d-03 * tc(2) &
        +3.800833383333333d-07 * tc(3) &
        +5.510594150000000d-11 * tc(4) &
        -1.811178760000000d-14 * tc(5)
else
    !species 46: HNCO
    species(46) = &
        -1.665993440000000d+04 * invT &
        +1.460619875000000d+01 &
        -6.223951340000000d+00 * tc(1) &
        -1.589320020000000d-03 * tc(2) &
        +1.822979250000000d-07 * tc(3) &
        -1.422793025000000d-11 * tc(4) &
        +4.975109775000000d-16 * tc(5)
end if

! species with midpoint at T=1382 kelvin
if (T < 1382.000000d0) then
    ! species 44: HCNO
    species(44) = &
        +1.929902520000000d+04 * invT &
        -8.086017310000001d+00 &
        -2.647279890000000d+00 * tc(1) &
        -6.375267100000000d-03 * tc(2) &
        +1.746570600000000d-06 * tc(3) &
        -3.678606966666667d-10 * tc(4) &
        +3.787607330000000d-14 * tc(5)
else
    !species 44: HCNO
    species(44) = &
        +1.796613390000000d+04 * invT &
        +1.692926446000000d+01 &
        -6.598604560000000d+00 * tc(1) &
        -1.513893130000000d-03 * tc(2) &
        +1.795072433333333d-07 * tc(3) &
        -1.430554400000000d-11 * tc(4) &
        +5.071969550000000d-16 * tc(5)
end if

end subroutine


! compute Cv/R at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine cv_R(species, tc)

implicit none

double precision, intent(out) :: species(53)
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
    ! species 8: H2O2
    species(8) = &
        +3.27611269d+00 &
        -5.42822417d-04 * tc(2) &
        +1.67335701d-05 * tc(3) &
        -2.15770813d-08 * tc(4) &
        +8.62454363d-12 * tc(5)
    ! species 9: C
    species(9) = &
        +1.55423955d+00 &
        -3.21537724d-04 * tc(2) &
        +7.33792245d-07 * tc(3) &
        -7.32234889d-10 * tc(4) &
        +2.66521446d-13 * tc(5)
    ! species 10: CH
    species(10) = &
        +2.48981665d+00 &
        +3.23835541d-04 * tc(2) &
        -1.68899065d-06 * tc(3) &
        +3.16217327d-09 * tc(4) &
        -1.40609067d-12 * tc(5)
    ! species 11: CH2
    species(11) = &
        +2.76267867d+00 &
        +9.68872143d-04 * tc(2) &
        +2.79489841d-06 * tc(3) &
        -3.85091153d-09 * tc(4) &
        +1.68741719d-12 * tc(5)
    ! species 12: CH2(S)
    species(12) = &
        +3.19860411d+00 &
        -2.36661419d-03 * tc(2) &
        +8.23296220d-06 * tc(3) &
        -6.68815981d-09 * tc(4) &
        +1.94314737d-12 * tc(5)
    ! species 13: CH3
    species(13) = &
        +2.67359040d+00 &
        +2.01095175d-03 * tc(2) &
        +5.73021856d-06 * tc(3) &
        -6.87117425d-09 * tc(4) &
        +2.54385734d-12 * tc(5)
    ! species 14: CH4
    species(14) = &
        +4.14987613d+00 &
        -1.36709788d-02 * tc(2) &
        +4.91800599d-05 * tc(3) &
        -4.84743026d-08 * tc(4) &
        +1.66693956d-11 * tc(5)
    ! species 15: CO
    species(15) = &
        +2.57953347d+00 &
        -6.10353680d-04 * tc(2) &
        +1.01681433d-06 * tc(3) &
        +9.07005884d-10 * tc(4) &
        -9.04424499d-13 * tc(5)
    ! species 16: CO2
    species(16) = &
        +1.35677352d+00 &
        +8.98459677d-03 * tc(2) &
        -7.12356269d-06 * tc(3) &
        +2.45919022d-09 * tc(4) &
        -1.43699548d-13 * tc(5)
    ! species 17: HCO
    species(17) = &
        +3.22118584d+00 &
        -3.24392532d-03 * tc(2) &
        +1.37799446d-05 * tc(3) &
        -1.33144093d-08 * tc(4) &
        +4.33768865d-12 * tc(5)
    ! species 18: CH2O
    species(18) = &
        +3.79372315d+00 &
        -9.90833369d-03 * tc(2) &
        +3.73220008d-05 * tc(3) &
        -3.79285261d-08 * tc(4) &
        +1.31772652d-11 * tc(5)
    ! species 19: CH2OH
    species(19) = &
        +2.86388918d+00 &
        +5.59672304d-03 * tc(2) &
        +5.93271791d-06 * tc(3) &
        -1.04532012d-08 * tc(4) &
        +4.36967278d-12 * tc(5)
    ! species 20: CH3O
    species(20) = &
        +1.10620400d+00 &
        +7.21659500d-03 * tc(2) &
        +5.33847200d-06 * tc(3) &
        -7.37763600d-09 * tc(4) &
        +2.07561000d-12 * tc(5)
    ! species 21: CH3OH
    species(21) = &
        +4.71539582d+00 &
        -1.52309129d-02 * tc(2) &
        +6.52441155d-05 * tc(3) &
        -7.10806889d-08 * tc(4) &
        +2.61352698d-11 * tc(5)
    ! species 22: C2H
    species(22) = &
        +1.88965733d+00 &
        +1.34099611d-02 * tc(2) &
        -2.84769501d-05 * tc(3) &
        +2.94791045d-08 * tc(4) &
        -1.09331511d-11 * tc(5)
    ! species 23: C2H2
    species(23) = &
        -1.91318906d-01 &
        +2.33615629d-02 * tc(2) &
        -3.55171815d-05 * tc(3) &
        +2.80152437d-08 * tc(4) &
        -8.50072974d-12 * tc(5)
    ! species 24: C2H3
    species(24) = &
        +2.21246645d+00 &
        +1.51479162d-03 * tc(2) &
        +2.59209412d-05 * tc(3) &
        -3.57657847d-08 * tc(4) &
        +1.47150873d-11 * tc(5)
    ! species 25: C2H4
    species(25) = &
        +2.95920148d+00 &
        -7.57052247d-03 * tc(2) &
        +5.70990292d-05 * tc(3) &
        -6.91588753d-08 * tc(4) &
        +2.69884373d-11 * tc(5)
    ! species 26: C2H5
    species(26) = &
        +3.30646568d+00 &
        -4.18658892d-03 * tc(2) &
        +4.97142807d-05 * tc(3) &
        -5.99126606d-08 * tc(4) &
        +2.30509004d-11 * tc(5)
    ! species 27: C2H6
    species(27) = &
        +3.29142492d+00 &
        -5.50154270d-03 * tc(2) &
        +5.99438288d-05 * tc(3) &
        -7.08466285d-08 * tc(4) &
        +2.68685771d-11 * tc(5)
    ! species 28: HCCO
    species(28) = &
        +1.25172140d+00 &
        +1.76550210d-02 * tc(2) &
        -2.37291010d-05 * tc(3) &
        +1.72757590d-08 * tc(4) &
        -5.06648110d-12 * tc(5)
    ! species 29: CH2CO
    species(29) = &
        +1.13583630d+00 &
        +1.81188721d-02 * tc(2) &
        -1.73947474d-05 * tc(3) &
        +9.34397568d-09 * tc(4) &
        -2.01457615d-12 * tc(5)
    ! species 30: HCCOH
    species(30) = &
        +2.42373300d-01 &
        +3.10722010d-02 * tc(2) &
        -5.08668640d-05 * tc(3) &
        +4.31371310d-08 * tc(4) &
        -1.40145940d-11 * tc(5)
    ! species 31: N
    species(31) = &
        +1.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5)
    ! species 32: NH
    species(32) = &
        +2.49290850d+00 &
        +3.11791980d-04 * tc(2) &
        -1.48904840d-06 * tc(3) &
        +2.48164420d-09 * tc(4) &
        -1.03569670d-12 * tc(5)
    ! species 33: NH2
    species(33) = &
        +3.20400290d+00 &
        -2.10613850d-03 * tc(2) &
        +7.10683480d-06 * tc(3) &
        -5.61151970d-09 * tc(4) &
        +1.64407170d-12 * tc(5)
    ! species 34: NH3
    species(34) = &
        +3.28602740d+00 &
        -4.66052300d-03 * tc(2) &
        +2.17185130d-05 * tc(3) &
        -2.28088870d-08 * tc(4) &
        +8.26380460d-12 * tc(5)
    ! species 35: NNH
    species(35) = &
        +3.34469270d+00 &
        -4.84970720d-03 * tc(2) &
        +2.00594590d-05 * tc(3) &
        -2.17264640d-08 * tc(4) &
        +7.94695390d-12 * tc(5)
    ! species 36: NO
    species(36) = &
        +3.21847630d+00 &
        -4.63897600d-03 * tc(2) &
        +1.10410220d-05 * tc(3) &
        -9.33613540d-09 * tc(4) &
        +2.80357700d-12 * tc(5)
    ! species 37: NO2
    species(37) = &
        +2.94403120d+00 &
        -1.58542900d-03 * tc(2) &
        +1.66578120d-05 * tc(3) &
        -2.04754260d-08 * tc(4) &
        +7.83505640d-12 * tc(5)
    ! species 38: N2O
    species(38) = &
        +1.25715020d+00 &
        +1.13047280d-02 * tc(2) &
        -1.36713190d-05 * tc(3) &
        +9.68198060d-09 * tc(4) &
        -2.93071820d-12 * tc(5)
    ! species 39: HNO
    species(39) = &
        +3.53349160d+00 &
        -5.66961710d-03 * tc(2) &
        +1.84732070d-05 * tc(3) &
        -1.71370940d-08 * tc(4) &
        +5.54545730d-12 * tc(5)
    ! species 40: CN
    species(40) = &
        +2.61293510d+00 &
        -9.55513270d-04 * tc(2) &
        +2.14429770d-06 * tc(3) &
        -3.15163230d-10 * tc(4) &
        -4.64303560d-13 * tc(5)
    ! species 41: HCN
    species(41) = &
        +1.25898860d+00 &
        +1.00511700d-02 * tc(2) &
        -1.33517630d-05 * tc(3) &
        +1.00923490d-08 * tc(4) &
        -3.00890280d-12 * tc(5)
    ! species 42: H2CN
    species(42) = &
        +1.85166100d+00 &
        +5.69523310d-03 * tc(2) &
        +1.07114000d-06 * tc(3) &
        -1.62261200d-09 * tc(4) &
        -2.35110810d-13 * tc(5)
    ! species 43: HCNN
    species(43) = &
        +1.52431940d+00 &
        +1.59606190d-02 * tc(2) &
        -1.88163540d-05 * tc(3) &
        +1.21255400d-08 * tc(4) &
        -3.23573780d-12 * tc(5)
    ! species 47: NCO
    species(47) = &
        +1.82693080d+00 &
        +8.80516880d-03 * tc(2) &
        -8.38661340d-06 * tc(3) &
        +4.80169640d-09 * tc(4) &
        -1.33135950d-12 * tc(5)
    ! species 48: N2
    species(48) = &
        +2.29867700d+00 &
        +1.40824040d-03 * tc(2) &
        -3.96322200d-06 * tc(3) &
        +5.64151500d-09 * tc(4) &
        -2.44485400d-12 * tc(5)
    ! species 49: AR
    species(49) = &
        +1.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5)
    ! species 50: C3H7
    species(50) = &
        +5.15518000d-02 &
        +2.59919800d-02 * tc(2) &
        +2.38005400d-06 * tc(3) &
        -1.96095690d-08 * tc(4) &
        +9.37324700d-12 * tc(5)
    ! species 51: C3H8
    species(51) = &
        -6.64461900d-02 &
        +2.64245790d-02 * tc(2) &
        +6.10597270d-06 * tc(3) &
        -2.19774990d-08 * tc(4) &
        +9.51492530d-12 * tc(5)
    ! species 52: CH2CHO
    species(52) = &
        +2.40906200d+00 &
        +1.07385740d-02 * tc(2) &
        +1.89149200d-06 * tc(3) &
        -7.15858300d-09 * tc(4) &
        +2.86738500d-12 * tc(5)
    ! species 53: CH3CHO
    species(53) = &
        +3.72945950d+00 &
        -3.19328580d-03 * tc(2) &
        +4.75349210d-05 * tc(3) &
        -5.74586110d-08 * tc(4) &
        +2.19311120d-11 * tc(5)
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
    !species 8: H2O2
    species(8) = &
        +3.16500285d+00 &
        +4.90831694d-03 * tc(2) &
        -1.90139225d-06 * tc(3) &
        +3.71185986d-10 * tc(4) &
        -2.87908305d-14 * tc(5)
    !species 9: C
    species(9) = &
        +1.49266888d+00 &
        +4.79889284d-05 * tc(2) &
        -7.24335020d-08 * tc(3) &
        +3.74291029d-11 * tc(4) &
        -4.87277893d-15 * tc(5)
    !species 10: CH
    species(10) = &
        +1.87846473d+00 &
        +9.70913681d-04 * tc(2) &
        +1.44445655d-07 * tc(3) &
        -1.30687849d-10 * tc(4) &
        +1.76079383d-14 * tc(5)
    !species 11: CH2
    species(11) = &
        +1.87410113d+00 &
        +3.65639292d-03 * tc(2) &
        -1.40894597d-06 * tc(3) &
        +2.60179549d-10 * tc(4) &
        -1.87727567d-14 * tc(5)
    !species 12: CH2(S)
    species(12) = &
        +1.29203842d+00 &
        +4.65588637d-03 * tc(2) &
        -2.01191947d-06 * tc(3) &
        +4.17906000d-10 * tc(4) &
        -3.39716365d-14 * tc(5)
    !species 13: CH3
    species(13) = &
        +1.28571772d+00 &
        +7.23990037d-03 * tc(2) &
        -2.98714348d-06 * tc(3) &
        +5.95684644d-10 * tc(4) &
        -4.67154394d-14 * tc(5)
    !species 14: CH4
    species(14) = &
        -9.25148505d-01 &
        +1.33909467d-02 * tc(2) &
        -5.73285809d-06 * tc(3) &
        +1.22292535d-09 * tc(4) &
        -1.01815230d-13 * tc(5)
    !species 15: CO
    species(15) = &
        +1.71518561d+00 &
        +2.06252743d-03 * tc(2) &
        -9.98825771d-07 * tc(3) &
        +2.30053008d-10 * tc(4) &
        -2.03647716d-14 * tc(5)
    !species 16: CO2
    species(16) = &
        +2.85746029d+00 &
        +4.41437026d-03 * tc(2) &
        -2.21481404d-06 * tc(3) &
        +5.23490188d-10 * tc(4) &
        -4.72084164d-14 * tc(5)
    !species 17: HCO
    species(17) = &
        +1.77217438d+00 &
        +4.95695526d-03 * tc(2) &
        -2.48445613d-06 * tc(3) &
        +5.89161778d-10 * tc(4) &
        -5.33508711d-14 * tc(5)
    !species 18: CH2O
    species(18) = &
        +7.60690080d-01 &
        +9.20000082d-03 * tc(2) &
        -4.42258813d-06 * tc(3) &
        +1.00641212d-09 * tc(4) &
        -8.83855640d-14 * tc(5)
    !species 19: CH2OH
    species(19) = &
        +2.69266569d+00 &
        +8.64576797d-03 * tc(2) &
        -3.75101120d-06 * tc(3) &
        +7.87234636d-10 * tc(4) &
        -6.48554201d-14 * tc(5)
    !species 20: CH3O
    species(20) = &
        +2.77079900d+00 &
        +7.87149700d-03 * tc(2) &
        -2.65638400d-06 * tc(3) &
        +3.94443100d-10 * tc(4) &
        -2.11261600d-14 * tc(5)
    !species 21: CH3OH
    species(21) = &
        +7.89707910d-01 &
        +1.40938292d-02 * tc(2) &
        -6.36500835d-06 * tc(3) &
        +1.38171085d-09 * tc(4) &
        -1.17060220d-13 * tc(5)
    !species 22: C2H
    species(22) = &
        +2.16780652d+00 &
        +4.75221902d-03 * tc(2) &
        -1.83787077d-06 * tc(3) &
        +3.04190252d-10 * tc(4) &
        -1.77232770d-14 * tc(5)
    !species 23: C2H2
    species(23) = &
        +3.14756964d+00 &
        +5.96166664d-03 * tc(2) &
        -2.37294852d-06 * tc(3) &
        +4.67412171d-10 * tc(4) &
        -3.61235213d-14 * tc(5)
    !species 24: C2H3
    species(24) = &
        +2.01672400d+00 &
        +1.03302292d-02 * tc(2) &
        -4.68082349d-06 * tc(3) &
        +1.01763288d-09 * tc(4) &
        -8.62607041d-14 * tc(5)
    !species 25: C2H4
    species(25) = &
        +1.03611116d+00 &
        +1.46454151d-02 * tc(2) &
        -6.71077915d-06 * tc(3) &
        +1.47222923d-09 * tc(4) &
        -1.25706061d-13 * tc(5)
    !species 26: C2H5
    species(26) = &
        +9.54656420d-01 &
        +1.73972722d-02 * tc(2) &
        -7.98206668d-06 * tc(3) &
        +1.75217689d-09 * tc(4) &
        -1.49641576d-13 * tc(5)
    !species 27: C2H6
    species(27) = &
        +7.18815000d-02 &
        +2.16852677d-02 * tc(2) &
        -1.00256067d-05 * tc(3) &
        +2.21412001d-09 * tc(4) &
        -1.90002890d-13 * tc(5)
    !species 28: HCCO
    species(28) = &
        +4.62820580d+00 &
        +4.08534010d-03 * tc(2) &
        -1.59345470d-06 * tc(3) &
        +2.86260520d-10 * tc(4) &
        -1.94078320d-14 * tc(5)
    !species 29: CH2CO
    species(29) = &
        +3.51129732d+00 &
        +9.00359745d-03 * tc(2) &
        -4.16939635d-06 * tc(3) &
        +9.23345882d-10 * tc(4) &
        -7.94838201d-14 * tc(5)
    !species 30: HCCOH
    species(30) = &
        +4.92382910d+00 &
        +6.79236000d-03 * tc(2) &
        -2.56585640d-06 * tc(3) &
        +4.49878410d-10 * tc(4) &
        -2.99401010d-14 * tc(5)
    !species 31: N
    species(31) = &
        +1.41594290d+00 &
        +1.74890650d-04 * tc(2) &
        -1.19023690d-07 * tc(3) &
        +3.02262450d-11 * tc(4) &
        -2.03609820d-15 * tc(5)
    !species 32: NH
    species(32) = &
        +1.78369280d+00 &
        +1.32984300d-03 * tc(2) &
        -4.24780470d-07 * tc(3) &
        +7.83485010d-11 * tc(4) &
        -5.50444700d-15 * tc(5)
    !species 33: NH2
    species(33) = &
        +1.83474210d+00 &
        +3.20730820d-03 * tc(2) &
        -9.33908040d-07 * tc(3) &
        +1.37029530d-10 * tc(4) &
        -7.92061440d-15 * tc(5)
    !species 34: NH3
    species(34) = &
        +1.63445210d+00 &
        +5.66625600d-03 * tc(2) &
        -1.72786760d-06 * tc(3) &
        +2.38671610d-10 * tc(4) &
        -1.25787860d-14 * tc(5)
    !species 35: NNH
    species(35) = &
        +2.76675440d+00 &
        +2.89150820d-03 * tc(2) &
        -1.04166200d-06 * tc(3) &
        +1.68425940d-10 * tc(4) &
        -1.00918960d-14 * tc(5)
    !species 36: NO
    species(36) = &
        +2.26060560d+00 &
        +1.19110430d-03 * tc(2) &
        -4.29170480d-07 * tc(3) &
        +6.94576690d-11 * tc(4) &
        -4.03360990d-15 * tc(5)
    !species 37: NO2
    species(37) = &
        +3.88475420d+00 &
        +2.17239560d-03 * tc(2) &
        -8.28069060d-07 * tc(3) &
        +1.57475100d-10 * tc(4) &
        -1.05108950d-14 * tc(5)
    !species 38: N2O
    species(38) = &
        +3.82307290d+00 &
        +2.62702510d-03 * tc(2) &
        -9.58508740d-07 * tc(3) &
        +1.60007120d-10 * tc(4) &
        -9.77523030d-15 * tc(5)
    !species 39: HNO
    species(39) = &
        +1.97925090d+00 &
        +3.49440590d-03 * tc(2) &
        -7.85497780d-07 * tc(3) &
        +5.74795940d-11 * tc(4) &
        -1.93359160d-16 * tc(5)
    !species 40: CN
    species(40) = &
        +2.74598050d+00 &
        +4.34507750d-05 * tc(2) &
        +2.97059840d-07 * tc(3) &
        -6.86518060d-11 * tc(4) &
        +4.41341730d-15 * tc(5)
    !species 41: HCN
    species(41) = &
        +2.80223920d+00 &
        +3.14642280d-03 * tc(2) &
        -1.06321850d-06 * tc(3) &
        +1.66197570d-10 * tc(4) &
        -9.79975700d-15 * tc(5)
    !species 42: H2CN
    species(42) = &
        +4.20970300d+00 &
        +2.96929110d-03 * tc(2) &
        -2.85558910d-07 * tc(3) &
        -1.63555000d-10 * tc(4) &
        +3.04325890d-14 * tc(5)
    !species 43: HCNN
    species(43) = &
        +4.89463620d+00 &
        +3.98959590d-03 * tc(2) &
        -1.59823800d-06 * tc(3) &
        +2.92493950d-10 * tc(4) &
        -2.00946860d-14 * tc(5)
    !species 47: NCO
    species(47) = &
        +4.15218450d+00 &
        +2.30517610d-03 * tc(2) &
        -8.80331530d-07 * tc(3) &
        +1.47890980d-10 * tc(4) &
        -9.09779960d-15 * tc(5)
    !species 48: N2
    species(48) = &
        +1.92664000d+00 &
        +1.48797680d-03 * tc(2) &
        -5.68476000d-07 * tc(3) &
        +1.00970380d-10 * tc(4) &
        -6.75335100d-15 * tc(5)
    !species 49: AR
    species(49) = &
        +1.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5)
    !species 50: C3H7
    species(50) = &
        +6.70269870d+00 &
        +1.60442030d-02 * tc(2) &
        -5.28332200d-06 * tc(3) &
        +7.62985900d-10 * tc(4) &
        -3.93922840d-14 * tc(5)
    !species 51: C3H8
    species(51) = &
        +6.53413680d+00 &
        +1.88722390d-02 * tc(2) &
        -6.27184910d-06 * tc(3) &
        +9.14756490d-10 * tc(4) &
        -4.78380690d-14 * tc(5)
    !species 52: CH2CHO
    species(52) = &
        +4.97567000d+00 &
        +8.13059100d-03 * tc(2) &
        -2.74362400d-06 * tc(3) &
        +4.07030400d-10 * tc(4) &
        -2.17601700d-14 * tc(5)
    !species 53: CH3CHO
    species(53) = &
        +4.40411080d+00 &
        +1.17230590d-02 * tc(2) &
        -4.22631370d-06 * tc(3) &
        +6.83724510d-10 * tc(4) &
        -4.09848630d-14 * tc(5)
end if

! species with midpoint at T=1368 kelvin
if (T < 1368.000000d0) then
    ! species 45: HOCN
    species(45) = &
        +2.78604952d+00 &
        +6.88667922d-03 * tc(2) &
        -3.21487864d-06 * tc(3) &
        +5.17195767d-10 * tc(4) &
        +1.19360788d-14 * tc(5)
else
    !species 45: HOCN
    species(45) = &
        +4.89784885d+00 &
        +3.16789393d-03 * tc(2) &
        -1.11801064d-06 * tc(3) &
        +1.77243144d-10 * tc(4) &
        -1.04339177d-14 * tc(5)
end if

! species with midpoint at T=1478 kelvin
if (T < 1478.000000d0) then
    ! species 46: HNCO
    species(46) = &
        +2.63096317d+00 &
        +7.30282357d-03 * tc(2) &
        -2.28050003d-06 * tc(3) &
        -6.61271298d-10 * tc(4) &
        +3.62235752d-13 * tc(5)
else
    !species 46: HNCO
    species(46) = &
        +5.22395134d+00 &
        +3.17864004d-03 * tc(2) &
        -1.09378755d-06 * tc(3) &
        +1.70735163d-10 * tc(4) &
        -9.95021955d-15 * tc(5)
end if

! species with midpoint at T=1382 kelvin
if (T < 1382.000000d0) then
    ! species 44: HCNO
    species(44) = &
        +1.64727989d+00 &
        +1.27505342d-02 * tc(2) &
        -1.04794236d-05 * tc(3) &
        +4.41432836d-09 * tc(4) &
        -7.57521466d-13 * tc(5)
else
    !species 44: HCNO
    species(44) = &
        +5.59860456d+00 &
        +3.02778626d-03 * tc(2) &
        -1.07704346d-06 * tc(3) &
        +1.71666528d-10 * tc(4) &
        -1.01439391d-14 * tc(5)
end if

end subroutine


! compute Cp/R at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine cp_R(species, tc)

implicit none

double precision, intent(out) :: species(53)
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
    ! species 8: H2O2
    species(8) = &
        +4.27611269d+00 &
        -5.42822417d-04 * tc(2) &
        +1.67335701d-05 * tc(3) &
        -2.15770813d-08 * tc(4) &
        +8.62454363d-12 * tc(5)
    ! species 9: C
    species(9) = &
        +2.55423955d+00 &
        -3.21537724d-04 * tc(2) &
        +7.33792245d-07 * tc(3) &
        -7.32234889d-10 * tc(4) &
        +2.66521446d-13 * tc(5)
    ! species 10: CH
    species(10) = &
        +3.48981665d+00 &
        +3.23835541d-04 * tc(2) &
        -1.68899065d-06 * tc(3) &
        +3.16217327d-09 * tc(4) &
        -1.40609067d-12 * tc(5)
    ! species 11: CH2
    species(11) = &
        +3.76267867d+00 &
        +9.68872143d-04 * tc(2) &
        +2.79489841d-06 * tc(3) &
        -3.85091153d-09 * tc(4) &
        +1.68741719d-12 * tc(5)
    ! species 12: CH2(S)
    species(12) = &
        +4.19860411d+00 &
        -2.36661419d-03 * tc(2) &
        +8.23296220d-06 * tc(3) &
        -6.68815981d-09 * tc(4) &
        +1.94314737d-12 * tc(5)
    ! species 13: CH3
    species(13) = &
        +3.67359040d+00 &
        +2.01095175d-03 * tc(2) &
        +5.73021856d-06 * tc(3) &
        -6.87117425d-09 * tc(4) &
        +2.54385734d-12 * tc(5)
    ! species 14: CH4
    species(14) = &
        +5.14987613d+00 &
        -1.36709788d-02 * tc(2) &
        +4.91800599d-05 * tc(3) &
        -4.84743026d-08 * tc(4) &
        +1.66693956d-11 * tc(5)
    ! species 15: CO
    species(15) = &
        +3.57953347d+00 &
        -6.10353680d-04 * tc(2) &
        +1.01681433d-06 * tc(3) &
        +9.07005884d-10 * tc(4) &
        -9.04424499d-13 * tc(5)
    ! species 16: CO2
    species(16) = &
        +2.35677352d+00 &
        +8.98459677d-03 * tc(2) &
        -7.12356269d-06 * tc(3) &
        +2.45919022d-09 * tc(4) &
        -1.43699548d-13 * tc(5)
    ! species 17: HCO
    species(17) = &
        +4.22118584d+00 &
        -3.24392532d-03 * tc(2) &
        +1.37799446d-05 * tc(3) &
        -1.33144093d-08 * tc(4) &
        +4.33768865d-12 * tc(5)
    ! species 18: CH2O
    species(18) = &
        +4.79372315d+00 &
        -9.90833369d-03 * tc(2) &
        +3.73220008d-05 * tc(3) &
        -3.79285261d-08 * tc(4) &
        +1.31772652d-11 * tc(5)
    ! species 19: CH2OH
    species(19) = &
        +3.86388918d+00 &
        +5.59672304d-03 * tc(2) &
        +5.93271791d-06 * tc(3) &
        -1.04532012d-08 * tc(4) &
        +4.36967278d-12 * tc(5)
    ! species 20: CH3O
    species(20) = &
        +2.10620400d+00 &
        +7.21659500d-03 * tc(2) &
        +5.33847200d-06 * tc(3) &
        -7.37763600d-09 * tc(4) &
        +2.07561000d-12 * tc(5)
    ! species 21: CH3OH
    species(21) = &
        +5.71539582d+00 &
        -1.52309129d-02 * tc(2) &
        +6.52441155d-05 * tc(3) &
        -7.10806889d-08 * tc(4) &
        +2.61352698d-11 * tc(5)
    ! species 22: C2H
    species(22) = &
        +2.88965733d+00 &
        +1.34099611d-02 * tc(2) &
        -2.84769501d-05 * tc(3) &
        +2.94791045d-08 * tc(4) &
        -1.09331511d-11 * tc(5)
    ! species 23: C2H2
    species(23) = &
        +8.08681094d-01 &
        +2.33615629d-02 * tc(2) &
        -3.55171815d-05 * tc(3) &
        +2.80152437d-08 * tc(4) &
        -8.50072974d-12 * tc(5)
    ! species 24: C2H3
    species(24) = &
        +3.21246645d+00 &
        +1.51479162d-03 * tc(2) &
        +2.59209412d-05 * tc(3) &
        -3.57657847d-08 * tc(4) &
        +1.47150873d-11 * tc(5)
    ! species 25: C2H4
    species(25) = &
        +3.95920148d+00 &
        -7.57052247d-03 * tc(2) &
        +5.70990292d-05 * tc(3) &
        -6.91588753d-08 * tc(4) &
        +2.69884373d-11 * tc(5)
    ! species 26: C2H5
    species(26) = &
        +4.30646568d+00 &
        -4.18658892d-03 * tc(2) &
        +4.97142807d-05 * tc(3) &
        -5.99126606d-08 * tc(4) &
        +2.30509004d-11 * tc(5)
    ! species 27: C2H6
    species(27) = &
        +4.29142492d+00 &
        -5.50154270d-03 * tc(2) &
        +5.99438288d-05 * tc(3) &
        -7.08466285d-08 * tc(4) &
        +2.68685771d-11 * tc(5)
    ! species 28: HCCO
    species(28) = &
        +2.25172140d+00 &
        +1.76550210d-02 * tc(2) &
        -2.37291010d-05 * tc(3) &
        +1.72757590d-08 * tc(4) &
        -5.06648110d-12 * tc(5)
    ! species 29: CH2CO
    species(29) = &
        +2.13583630d+00 &
        +1.81188721d-02 * tc(2) &
        -1.73947474d-05 * tc(3) &
        +9.34397568d-09 * tc(4) &
        -2.01457615d-12 * tc(5)
    ! species 30: HCCOH
    species(30) = &
        +1.24237330d+00 &
        +3.10722010d-02 * tc(2) &
        -5.08668640d-05 * tc(3) &
        +4.31371310d-08 * tc(4) &
        -1.40145940d-11 * tc(5)
    ! species 31: N
    species(31) = &
        +2.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5)
    ! species 32: NH
    species(32) = &
        +3.49290850d+00 &
        +3.11791980d-04 * tc(2) &
        -1.48904840d-06 * tc(3) &
        +2.48164420d-09 * tc(4) &
        -1.03569670d-12 * tc(5)
    ! species 33: NH2
    species(33) = &
        +4.20400290d+00 &
        -2.10613850d-03 * tc(2) &
        +7.10683480d-06 * tc(3) &
        -5.61151970d-09 * tc(4) &
        +1.64407170d-12 * tc(5)
    ! species 34: NH3
    species(34) = &
        +4.28602740d+00 &
        -4.66052300d-03 * tc(2) &
        +2.17185130d-05 * tc(3) &
        -2.28088870d-08 * tc(4) &
        +8.26380460d-12 * tc(5)
    ! species 35: NNH
    species(35) = &
        +4.34469270d+00 &
        -4.84970720d-03 * tc(2) &
        +2.00594590d-05 * tc(3) &
        -2.17264640d-08 * tc(4) &
        +7.94695390d-12 * tc(5)
    ! species 36: NO
    species(36) = &
        +4.21847630d+00 &
        -4.63897600d-03 * tc(2) &
        +1.10410220d-05 * tc(3) &
        -9.33613540d-09 * tc(4) &
        +2.80357700d-12 * tc(5)
    ! species 37: NO2
    species(37) = &
        +3.94403120d+00 &
        -1.58542900d-03 * tc(2) &
        +1.66578120d-05 * tc(3) &
        -2.04754260d-08 * tc(4) &
        +7.83505640d-12 * tc(5)
    ! species 38: N2O
    species(38) = &
        +2.25715020d+00 &
        +1.13047280d-02 * tc(2) &
        -1.36713190d-05 * tc(3) &
        +9.68198060d-09 * tc(4) &
        -2.93071820d-12 * tc(5)
    ! species 39: HNO
    species(39) = &
        +4.53349160d+00 &
        -5.66961710d-03 * tc(2) &
        +1.84732070d-05 * tc(3) &
        -1.71370940d-08 * tc(4) &
        +5.54545730d-12 * tc(5)
    ! species 40: CN
    species(40) = &
        +3.61293510d+00 &
        -9.55513270d-04 * tc(2) &
        +2.14429770d-06 * tc(3) &
        -3.15163230d-10 * tc(4) &
        -4.64303560d-13 * tc(5)
    ! species 41: HCN
    species(41) = &
        +2.25898860d+00 &
        +1.00511700d-02 * tc(2) &
        -1.33517630d-05 * tc(3) &
        +1.00923490d-08 * tc(4) &
        -3.00890280d-12 * tc(5)
    ! species 42: H2CN
    species(42) = &
        +2.85166100d+00 &
        +5.69523310d-03 * tc(2) &
        +1.07114000d-06 * tc(3) &
        -1.62261200d-09 * tc(4) &
        -2.35110810d-13 * tc(5)
    ! species 43: HCNN
    species(43) = &
        +2.52431940d+00 &
        +1.59606190d-02 * tc(2) &
        -1.88163540d-05 * tc(3) &
        +1.21255400d-08 * tc(4) &
        -3.23573780d-12 * tc(5)
    ! species 47: NCO
    species(47) = &
        +2.82693080d+00 &
        +8.80516880d-03 * tc(2) &
        -8.38661340d-06 * tc(3) &
        +4.80169640d-09 * tc(4) &
        -1.33135950d-12 * tc(5)
    ! species 48: N2
    species(48) = &
        +3.29867700d+00 &
        +1.40824040d-03 * tc(2) &
        -3.96322200d-06 * tc(3) &
        +5.64151500d-09 * tc(4) &
        -2.44485400d-12 * tc(5)
    ! species 49: AR
    species(49) = &
        +2.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5)
    ! species 50: C3H7
    species(50) = &
        +1.05155180d+00 &
        +2.59919800d-02 * tc(2) &
        +2.38005400d-06 * tc(3) &
        -1.96095690d-08 * tc(4) &
        +9.37324700d-12 * tc(5)
    ! species 51: C3H8
    species(51) = &
        +9.33553810d-01 &
        +2.64245790d-02 * tc(2) &
        +6.10597270d-06 * tc(3) &
        -2.19774990d-08 * tc(4) &
        +9.51492530d-12 * tc(5)
    ! species 52: CH2CHO
    species(52) = &
        +3.40906200d+00 &
        +1.07385740d-02 * tc(2) &
        +1.89149200d-06 * tc(3) &
        -7.15858300d-09 * tc(4) &
        +2.86738500d-12 * tc(5)
    ! species 53: CH3CHO
    species(53) = &
        +4.72945950d+00 &
        -3.19328580d-03 * tc(2) &
        +4.75349210d-05 * tc(3) &
        -5.74586110d-08 * tc(4) &
        +2.19311120d-11 * tc(5)
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
    !species 8: H2O2
    species(8) = &
        +4.16500285d+00 &
        +4.90831694d-03 * tc(2) &
        -1.90139225d-06 * tc(3) &
        +3.71185986d-10 * tc(4) &
        -2.87908305d-14 * tc(5)
    !species 9: C
    species(9) = &
        +2.49266888d+00 &
        +4.79889284d-05 * tc(2) &
        -7.24335020d-08 * tc(3) &
        +3.74291029d-11 * tc(4) &
        -4.87277893d-15 * tc(5)
    !species 10: CH
    species(10) = &
        +2.87846473d+00 &
        +9.70913681d-04 * tc(2) &
        +1.44445655d-07 * tc(3) &
        -1.30687849d-10 * tc(4) &
        +1.76079383d-14 * tc(5)
    !species 11: CH2
    species(11) = &
        +2.87410113d+00 &
        +3.65639292d-03 * tc(2) &
        -1.40894597d-06 * tc(3) &
        +2.60179549d-10 * tc(4) &
        -1.87727567d-14 * tc(5)
    !species 12: CH2(S)
    species(12) = &
        +2.29203842d+00 &
        +4.65588637d-03 * tc(2) &
        -2.01191947d-06 * tc(3) &
        +4.17906000d-10 * tc(4) &
        -3.39716365d-14 * tc(5)
    !species 13: CH3
    species(13) = &
        +2.28571772d+00 &
        +7.23990037d-03 * tc(2) &
        -2.98714348d-06 * tc(3) &
        +5.95684644d-10 * tc(4) &
        -4.67154394d-14 * tc(5)
    !species 14: CH4
    species(14) = &
        +7.48514950d-02 &
        +1.33909467d-02 * tc(2) &
        -5.73285809d-06 * tc(3) &
        +1.22292535d-09 * tc(4) &
        -1.01815230d-13 * tc(5)
    !species 15: CO
    species(15) = &
        +2.71518561d+00 &
        +2.06252743d-03 * tc(2) &
        -9.98825771d-07 * tc(3) &
        +2.30053008d-10 * tc(4) &
        -2.03647716d-14 * tc(5)
    !species 16: CO2
    species(16) = &
        +3.85746029d+00 &
        +4.41437026d-03 * tc(2) &
        -2.21481404d-06 * tc(3) &
        +5.23490188d-10 * tc(4) &
        -4.72084164d-14 * tc(5)
    !species 17: HCO
    species(17) = &
        +2.77217438d+00 &
        +4.95695526d-03 * tc(2) &
        -2.48445613d-06 * tc(3) &
        +5.89161778d-10 * tc(4) &
        -5.33508711d-14 * tc(5)
    !species 18: CH2O
    species(18) = &
        +1.76069008d+00 &
        +9.20000082d-03 * tc(2) &
        -4.42258813d-06 * tc(3) &
        +1.00641212d-09 * tc(4) &
        -8.83855640d-14 * tc(5)
    !species 19: CH2OH
    species(19) = &
        +3.69266569d+00 &
        +8.64576797d-03 * tc(2) &
        -3.75101120d-06 * tc(3) &
        +7.87234636d-10 * tc(4) &
        -6.48554201d-14 * tc(5)
    !species 20: CH3O
    species(20) = &
        +3.77079900d+00 &
        +7.87149700d-03 * tc(2) &
        -2.65638400d-06 * tc(3) &
        +3.94443100d-10 * tc(4) &
        -2.11261600d-14 * tc(5)
    !species 21: CH3OH
    species(21) = &
        +1.78970791d+00 &
        +1.40938292d-02 * tc(2) &
        -6.36500835d-06 * tc(3) &
        +1.38171085d-09 * tc(4) &
        -1.17060220d-13 * tc(5)
    !species 22: C2H
    species(22) = &
        +3.16780652d+00 &
        +4.75221902d-03 * tc(2) &
        -1.83787077d-06 * tc(3) &
        +3.04190252d-10 * tc(4) &
        -1.77232770d-14 * tc(5)
    !species 23: C2H2
    species(23) = &
        +4.14756964d+00 &
        +5.96166664d-03 * tc(2) &
        -2.37294852d-06 * tc(3) &
        +4.67412171d-10 * tc(4) &
        -3.61235213d-14 * tc(5)
    !species 24: C2H3
    species(24) = &
        +3.01672400d+00 &
        +1.03302292d-02 * tc(2) &
        -4.68082349d-06 * tc(3) &
        +1.01763288d-09 * tc(4) &
        -8.62607041d-14 * tc(5)
    !species 25: C2H4
    species(25) = &
        +2.03611116d+00 &
        +1.46454151d-02 * tc(2) &
        -6.71077915d-06 * tc(3) &
        +1.47222923d-09 * tc(4) &
        -1.25706061d-13 * tc(5)
    !species 26: C2H5
    species(26) = &
        +1.95465642d+00 &
        +1.73972722d-02 * tc(2) &
        -7.98206668d-06 * tc(3) &
        +1.75217689d-09 * tc(4) &
        -1.49641576d-13 * tc(5)
    !species 27: C2H6
    species(27) = &
        +1.07188150d+00 &
        +2.16852677d-02 * tc(2) &
        -1.00256067d-05 * tc(3) &
        +2.21412001d-09 * tc(4) &
        -1.90002890d-13 * tc(5)
    !species 28: HCCO
    species(28) = &
        +5.62820580d+00 &
        +4.08534010d-03 * tc(2) &
        -1.59345470d-06 * tc(3) &
        +2.86260520d-10 * tc(4) &
        -1.94078320d-14 * tc(5)
    !species 29: CH2CO
    species(29) = &
        +4.51129732d+00 &
        +9.00359745d-03 * tc(2) &
        -4.16939635d-06 * tc(3) &
        +9.23345882d-10 * tc(4) &
        -7.94838201d-14 * tc(5)
    !species 30: HCCOH
    species(30) = &
        +5.92382910d+00 &
        +6.79236000d-03 * tc(2) &
        -2.56585640d-06 * tc(3) &
        +4.49878410d-10 * tc(4) &
        -2.99401010d-14 * tc(5)
    !species 31: N
    species(31) = &
        +2.41594290d+00 &
        +1.74890650d-04 * tc(2) &
        -1.19023690d-07 * tc(3) &
        +3.02262450d-11 * tc(4) &
        -2.03609820d-15 * tc(5)
    !species 32: NH
    species(32) = &
        +2.78369280d+00 &
        +1.32984300d-03 * tc(2) &
        -4.24780470d-07 * tc(3) &
        +7.83485010d-11 * tc(4) &
        -5.50444700d-15 * tc(5)
    !species 33: NH2
    species(33) = &
        +2.83474210d+00 &
        +3.20730820d-03 * tc(2) &
        -9.33908040d-07 * tc(3) &
        +1.37029530d-10 * tc(4) &
        -7.92061440d-15 * tc(5)
    !species 34: NH3
    species(34) = &
        +2.63445210d+00 &
        +5.66625600d-03 * tc(2) &
        -1.72786760d-06 * tc(3) &
        +2.38671610d-10 * tc(4) &
        -1.25787860d-14 * tc(5)
    !species 35: NNH
    species(35) = &
        +3.76675440d+00 &
        +2.89150820d-03 * tc(2) &
        -1.04166200d-06 * tc(3) &
        +1.68425940d-10 * tc(4) &
        -1.00918960d-14 * tc(5)
    !species 36: NO
    species(36) = &
        +3.26060560d+00 &
        +1.19110430d-03 * tc(2) &
        -4.29170480d-07 * tc(3) &
        +6.94576690d-11 * tc(4) &
        -4.03360990d-15 * tc(5)
    !species 37: NO2
    species(37) = &
        +4.88475420d+00 &
        +2.17239560d-03 * tc(2) &
        -8.28069060d-07 * tc(3) &
        +1.57475100d-10 * tc(4) &
        -1.05108950d-14 * tc(5)
    !species 38: N2O
    species(38) = &
        +4.82307290d+00 &
        +2.62702510d-03 * tc(2) &
        -9.58508740d-07 * tc(3) &
        +1.60007120d-10 * tc(4) &
        -9.77523030d-15 * tc(5)
    !species 39: HNO
    species(39) = &
        +2.97925090d+00 &
        +3.49440590d-03 * tc(2) &
        -7.85497780d-07 * tc(3) &
        +5.74795940d-11 * tc(4) &
        -1.93359160d-16 * tc(5)
    !species 40: CN
    species(40) = &
        +3.74598050d+00 &
        +4.34507750d-05 * tc(2) &
        +2.97059840d-07 * tc(3) &
        -6.86518060d-11 * tc(4) &
        +4.41341730d-15 * tc(5)
    !species 41: HCN
    species(41) = &
        +3.80223920d+00 &
        +3.14642280d-03 * tc(2) &
        -1.06321850d-06 * tc(3) &
        +1.66197570d-10 * tc(4) &
        -9.79975700d-15 * tc(5)
    !species 42: H2CN
    species(42) = &
        +5.20970300d+00 &
        +2.96929110d-03 * tc(2) &
        -2.85558910d-07 * tc(3) &
        -1.63555000d-10 * tc(4) &
        +3.04325890d-14 * tc(5)
    !species 43: HCNN
    species(43) = &
        +5.89463620d+00 &
        +3.98959590d-03 * tc(2) &
        -1.59823800d-06 * tc(3) &
        +2.92493950d-10 * tc(4) &
        -2.00946860d-14 * tc(5)
    !species 47: NCO
    species(47) = &
        +5.15218450d+00 &
        +2.30517610d-03 * tc(2) &
        -8.80331530d-07 * tc(3) &
        +1.47890980d-10 * tc(4) &
        -9.09779960d-15 * tc(5)
    !species 48: N2
    species(48) = &
        +2.92664000d+00 &
        +1.48797680d-03 * tc(2) &
        -5.68476000d-07 * tc(3) &
        +1.00970380d-10 * tc(4) &
        -6.75335100d-15 * tc(5)
    !species 49: AR
    species(49) = &
        +2.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5)
    !species 50: C3H7
    species(50) = &
        +7.70269870d+00 &
        +1.60442030d-02 * tc(2) &
        -5.28332200d-06 * tc(3) &
        +7.62985900d-10 * tc(4) &
        -3.93922840d-14 * tc(5)
    !species 51: C3H8
    species(51) = &
        +7.53413680d+00 &
        +1.88722390d-02 * tc(2) &
        -6.27184910d-06 * tc(3) &
        +9.14756490d-10 * tc(4) &
        -4.78380690d-14 * tc(5)
    !species 52: CH2CHO
    species(52) = &
        +5.97567000d+00 &
        +8.13059100d-03 * tc(2) &
        -2.74362400d-06 * tc(3) &
        +4.07030400d-10 * tc(4) &
        -2.17601700d-14 * tc(5)
    !species 53: CH3CHO
    species(53) = &
        +5.40411080d+00 &
        +1.17230590d-02 * tc(2) &
        -4.22631370d-06 * tc(3) &
        +6.83724510d-10 * tc(4) &
        -4.09848630d-14 * tc(5)
end if

! species with midpoint at T=1368 kelvin
if (T < 1368.000000d0) then
    ! species 45: HOCN
    species(45) = &
        +3.78604952d+00 &
        +6.88667922d-03 * tc(2) &
        -3.21487864d-06 * tc(3) &
        +5.17195767d-10 * tc(4) &
        +1.19360788d-14 * tc(5)
else
    !species 45: HOCN
    species(45) = &
        +5.89784885d+00 &
        +3.16789393d-03 * tc(2) &
        -1.11801064d-06 * tc(3) &
        +1.77243144d-10 * tc(4) &
        -1.04339177d-14 * tc(5)
end if

! species with midpoint at T=1478 kelvin
if (T < 1478.000000d0) then
    ! species 46: HNCO
    species(46) = &
        +3.63096317d+00 &
        +7.30282357d-03 * tc(2) &
        -2.28050003d-06 * tc(3) &
        -6.61271298d-10 * tc(4) &
        +3.62235752d-13 * tc(5)
else
    !species 46: HNCO
    species(46) = &
        +6.22395134d+00 &
        +3.17864004d-03 * tc(2) &
        -1.09378755d-06 * tc(3) &
        +1.70735163d-10 * tc(4) &
        -9.95021955d-15 * tc(5)
end if

! species with midpoint at T=1382 kelvin
if (T < 1382.000000d0) then
    ! species 44: HCNO
    species(44) = &
        +2.64727989d+00 &
        +1.27505342d-02 * tc(2) &
        -1.04794236d-05 * tc(3) &
        +4.41432836d-09 * tc(4) &
        -7.57521466d-13 * tc(5)
else
    !species 44: HCNO
    species(44) = &
        +6.59860456d+00 &
        +3.02778626d-03 * tc(2) &
        -1.07704346d-06 * tc(3) &
        +1.71666528d-10 * tc(4) &
        -1.01439391d-14 * tc(5)
end if

end subroutine

! compute the e/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine speciesInternalEnergy(species, tc)

implicit none

double precision, intent(out) :: species(53)
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
    ! species 8: H2O2
    species(8) = &
        +3.27611269d+00 &
        -2.71411208d-04 * tc(2) &
        +5.57785670d-06 * tc(3) &
        -5.39427032d-09 * tc(4) &
        +1.72490873d-12 * tc(5) &
        -1.77025821d+04 * invT
    ! species 9: C
    species(9) = &
        +1.55423955d+00 &
        -1.60768862d-04 * tc(2) &
        +2.44597415d-07 * tc(3) &
        -1.83058722d-10 * tc(4) &
        +5.33042892d-14 * tc(5) &
        +8.54438832d+04 * invT
    ! species 10: CH
    species(10) = &
        +2.48981665d+00 &
        +1.61917771d-04 * tc(2) &
        -5.62996883d-07 * tc(3) &
        +7.90543317d-10 * tc(4) &
        -2.81218134d-13 * tc(5) &
        +7.07972934d+04 * invT
    ! species 11: CH2
    species(11) = &
        +2.76267867d+00 &
        +4.84436072d-04 * tc(2) &
        +9.31632803d-07 * tc(3) &
        -9.62727883d-10 * tc(4) &
        +3.37483438d-13 * tc(5) &
        +4.60040401d+04 * invT
    ! species 12: CH2(S)
    species(12) = &
        +3.19860411d+00 &
        -1.18330710d-03 * tc(2) &
        +2.74432073d-06 * tc(3) &
        -1.67203995d-09 * tc(4) &
        +3.88629474d-13 * tc(5) &
        +5.04968163d+04 * invT
    ! species 13: CH3
    species(13) = &
        +2.67359040d+00 &
        +1.00547588d-03 * tc(2) &
        +1.91007285d-06 * tc(3) &
        -1.71779356d-09 * tc(4) &
        +5.08771468d-13 * tc(5) &
        +1.64449988d+04 * invT
    ! species 14: CH4
    species(14) = &
        +4.14987613d+00 &
        -6.83548940d-03 * tc(2) &
        +1.63933533d-05 * tc(3) &
        -1.21185757d-08 * tc(4) &
        +3.33387912d-12 * tc(5) &
        -1.02466476d+04 * invT
    ! species 15: CO
    species(15) = &
        +2.57953347d+00 &
        -3.05176840d-04 * tc(2) &
        +3.38938110d-07 * tc(3) &
        +2.26751471d-10 * tc(4) &
        -1.80884900d-13 * tc(5) &
        -1.43440860d+04 * invT
    ! species 16: CO2
    species(16) = &
        +1.35677352d+00 &
        +4.49229839d-03 * tc(2) &
        -2.37452090d-06 * tc(3) &
        +6.14797555d-10 * tc(4) &
        -2.87399096d-14 * tc(5) &
        -4.83719697d+04 * invT
    ! species 17: HCO
    species(17) = &
        +3.22118584d+00 &
        -1.62196266d-03 * tc(2) &
        +4.59331487d-06 * tc(3) &
        -3.32860233d-09 * tc(4) &
        +8.67537730d-13 * tc(5) &
        +3.83956496d+03 * invT
    ! species 18: CH2O
    species(18) = &
        +3.79372315d+00 &
        -4.95416684d-03 * tc(2) &
        +1.24406669d-05 * tc(3) &
        -9.48213152d-09 * tc(4) &
        +2.63545304d-12 * tc(5) &
        -1.43089567d+04 * invT
    ! species 19: CH2OH
    species(19) = &
        +2.86388918d+00 &
        +2.79836152d-03 * tc(2) &
        +1.97757264d-06 * tc(3) &
        -2.61330030d-09 * tc(4) &
        +8.73934556d-13 * tc(5) &
        -3.19391367d+03 * invT
    ! species 20: CH3O
    species(20) = &
        +1.10620400d+00 &
        +3.60829750d-03 * tc(2) &
        +1.77949067d-06 * tc(3) &
        -1.84440900d-09 * tc(4) &
        +4.15122000d-13 * tc(5) &
        +9.78601100d+02 * invT
    ! species 21: CH3OH
    species(21) = &
        +4.71539582d+00 &
        -7.61545645d-03 * tc(2) &
        +2.17480385d-05 * tc(3) &
        -1.77701722d-08 * tc(4) &
        +5.22705396d-12 * tc(5) &
        -2.56427656d+04 * invT
    ! species 22: C2H
    species(22) = &
        +1.88965733d+00 &
        +6.70498055d-03 * tc(2) &
        -9.49231670d-06 * tc(3) &
        +7.36977613d-09 * tc(4) &
        -2.18663022d-12 * tc(5) &
        +6.68393932d+04 * invT
    ! species 23: C2H2
    species(23) = &
        -1.91318906d-01 &
        +1.16807815d-02 * tc(2) &
        -1.18390605d-05 * tc(3) &
        +7.00381092d-09 * tc(4) &
        -1.70014595d-12 * tc(5) &
        +2.64289807d+04 * invT
    ! species 24: C2H3
    species(24) = &
        +2.21246645d+00 &
        +7.57395810d-04 * tc(2) &
        +8.64031373d-06 * tc(3) &
        -8.94144617d-09 * tc(4) &
        +2.94301746d-12 * tc(5) &
        +3.48598468d+04 * invT
    ! species 25: C2H4
    species(25) = &
        +2.95920148d+00 &
        -3.78526124d-03 * tc(2) &
        +1.90330097d-05 * tc(3) &
        -1.72897188d-08 * tc(4) &
        +5.39768746d-12 * tc(5) &
        +5.08977593d+03 * invT
    ! species 26: C2H5
    species(26) = &
        +3.30646568d+00 &
        -2.09329446d-03 * tc(2) &
        +1.65714269d-05 * tc(3) &
        -1.49781651d-08 * tc(4) &
        +4.61018008d-12 * tc(5) &
        +1.28416265d+04 * invT
    ! species 27: C2H6
    species(27) = &
        +3.29142492d+00 &
        -2.75077135d-03 * tc(2) &
        +1.99812763d-05 * tc(3) &
        -1.77116571d-08 * tc(4) &
        +5.37371542d-12 * tc(5) &
        -1.15222055d+04 * invT
    ! species 28: HCCO
    species(28) = &
        +1.25172140d+00 &
        +8.82751050d-03 * tc(2) &
        -7.90970033d-06 * tc(3) &
        +4.31893975d-09 * tc(4) &
        -1.01329622d-12 * tc(5) &
        +2.00594490d+04 * invT
    ! species 29: CH2CO
    species(29) = &
        +1.13583630d+00 &
        +9.05943605d-03 * tc(2) &
        -5.79824913d-06 * tc(3) &
        +2.33599392d-09 * tc(4) &
        -4.02915230d-13 * tc(5) &
        -7.04291804d+03 * invT
    ! species 30: HCCOH
    species(30) = &
        +2.42373300d-01 &
        +1.55361005d-02 * tc(2) &
        -1.69556213d-05 * tc(3) &
        +1.07842828d-08 * tc(4) &
        -2.80291880d-12 * tc(5) &
        +8.03161430d+03 * invT
    ! species 31: N
    species(31) = &
        +1.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5) &
        +5.61046370d+04 * invT
    ! species 32: NH
    species(32) = &
        +2.49290850d+00 &
        +1.55895990d-04 * tc(2) &
        -4.96349467d-07 * tc(3) &
        +6.20411050d-10 * tc(4) &
        -2.07139340d-13 * tc(5) &
        +4.18806290d+04 * invT
    ! species 33: NH2
    species(33) = &
        +3.20400290d+00 &
        -1.05306925d-03 * tc(2) &
        +2.36894493d-06 * tc(3) &
        -1.40287992d-09 * tc(4) &
        +3.28814340d-13 * tc(5) &
        +2.18859100d+04 * invT
    ! species 34: NH3
    species(34) = &
        +3.28602740d+00 &
        -2.33026150d-03 * tc(2) &
        +7.23950433d-06 * tc(3) &
        -5.70222175d-09 * tc(4) &
        +1.65276092d-12 * tc(5) &
        -6.74172850d+03 * invT
    ! species 35: NNH
    species(35) = &
        +3.34469270d+00 &
        -2.42485360d-03 * tc(2) &
        +6.68648633d-06 * tc(3) &
        -5.43161600d-09 * tc(4) &
        +1.58939078d-12 * tc(5) &
        +2.87919730d+04 * invT
    ! species 36: NO
    species(36) = &
        +3.21847630d+00 &
        -2.31948800d-03 * tc(2) &
        +3.68034067d-06 * tc(3) &
        -2.33403385d-09 * tc(4) &
        +5.60715400d-13 * tc(5) &
        +9.84462300d+03 * invT
    ! species 37: NO2
    species(37) = &
        +2.94403120d+00 &
        -7.92714500d-04 * tc(2) &
        +5.55260400d-06 * tc(3) &
        -5.11885650d-09 * tc(4) &
        +1.56701128d-12 * tc(5) &
        +2.89661790d+03 * invT
    ! species 38: N2O
    species(38) = &
        +1.25715020d+00 &
        +5.65236400d-03 * tc(2) &
        -4.55710633d-06 * tc(3) &
        +2.42049515d-09 * tc(4) &
        -5.86143640d-13 * tc(5) &
        +8.74177440d+03 * invT
    ! species 39: HNO
    species(39) = &
        +3.53349160d+00 &
        -2.83480855d-03 * tc(2) &
        +6.15773567d-06 * tc(3) &
        -4.28427350d-09 * tc(4) &
        +1.10909146d-12 * tc(5) &
        +1.15482970d+04 * invT
    ! species 40: CN
    species(40) = &
        +2.61293510d+00 &
        -4.77756635d-04 * tc(2) &
        +7.14765900d-07 * tc(3) &
        -7.87908075d-11 * tc(4) &
        -9.28607120d-14 * tc(5) &
        +5.17083400d+04 * invT
    ! species 41: HCN
    species(41) = &
        +1.25898860d+00 &
        +5.02558500d-03 * tc(2) &
        -4.45058767d-06 * tc(3) &
        +2.52308725d-09 * tc(4) &
        -6.01780560d-13 * tc(5) &
        +1.47126330d+04 * invT
    ! species 42: H2CN
    species(42) = &
        +1.85166100d+00 &
        +2.84761655d-03 * tc(2) &
        +3.57046667d-07 * tc(3) &
        -4.05653000d-10 * tc(4) &
        -4.70221620d-14 * tc(5) &
        +2.86378200d+04 * invT
    ! species 43: HCNN
    species(43) = &
        +1.52431940d+00 &
        +7.98030950d-03 * tc(2) &
        -6.27211800d-06 * tc(3) &
        +3.03138500d-09 * tc(4) &
        -6.47147560d-13 * tc(5) &
        +5.42619840d+04 * invT
    ! species 47: NCO
    species(47) = &
        +1.82693080d+00 &
        +4.40258440d-03 * tc(2) &
        -2.79553780d-06 * tc(3) &
        +1.20042410d-09 * tc(4) &
        -2.66271900d-13 * tc(5) &
        +1.46824770d+04 * invT
    ! species 48: N2
    species(48) = &
        +2.29867700d+00 &
        +7.04120200d-04 * tc(2) &
        -1.32107400d-06 * tc(3) &
        +1.41037875d-09 * tc(4) &
        -4.88970800d-13 * tc(5) &
        -1.02089990d+03 * invT
    ! species 49: AR
    species(49) = &
        +1.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5) &
        -7.45375000d+02 * invT
    ! species 50: C3H7
    species(50) = &
        +5.15518000d-02 &
        +1.29959900d-02 * tc(2) &
        +7.93351333d-07 * tc(3) &
        -4.90239225d-09 * tc(4) &
        +1.87464940d-12 * tc(5) &
        +1.06318630d+04 * invT
    ! species 51: C3H8
    species(51) = &
        -6.64461900d-02 &
        +1.32122895d-02 * tc(2) &
        +2.03532423d-06 * tc(3) &
        -5.49437475d-09 * tc(4) &
        +1.90298506d-12 * tc(5) &
        -1.39585200d+04 * invT
    ! species 52: CH2CHO
    species(52) = &
        +2.40906200d+00 &
        +5.36928700d-03 * tc(2) &
        +6.30497333d-07 * tc(3) &
        -1.78964575d-09 * tc(4) &
        +5.73477000d-13 * tc(5) &
        +1.52147660d+03 * invT
    ! species 53: CH3CHO
    species(53) = &
        +3.72945950d+00 &
        -1.59664290d-03 * tc(2) &
        +1.58449737d-05 * tc(3) &
        -1.43646527d-08 * tc(4) &
        +4.38622240d-12 * tc(5) &
        -2.15728780d+04 * invT
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
    !species 8: H2O2
    species(8) = &
        +3.16500285d+00 &
        +2.45415847d-03 * tc(2) &
        -6.33797417d-07 * tc(3) &
        +9.27964965d-11 * tc(4) &
        -5.75816610d-15 * tc(5) &
        -1.78617877d+04 * invT
    !species 9: C
    species(9) = &
        +1.49266888d+00 &
        +2.39944642d-05 * tc(2) &
        -2.41445007d-08 * tc(3) &
        +9.35727573d-12 * tc(4) &
        -9.74555786d-16 * tc(5) &
        +8.54512953d+04 * invT
    !species 10: CH
    species(10) = &
        +1.87846473d+00 &
        +4.85456840d-04 * tc(2) &
        +4.81485517d-08 * tc(3) &
        -3.26719623d-11 * tc(4) &
        +3.52158766d-15 * tc(5) &
        +7.10124364d+04 * invT
    !species 11: CH2
    species(11) = &
        +1.87410113d+00 &
        +1.82819646d-03 * tc(2) &
        -4.69648657d-07 * tc(3) &
        +6.50448872d-11 * tc(4) &
        -3.75455134d-15 * tc(5) &
        +4.62636040d+04 * invT
    !species 12: CH2(S)
    species(12) = &
        +1.29203842d+00 &
        +2.32794318d-03 * tc(2) &
        -6.70639823d-07 * tc(3) &
        +1.04476500d-10 * tc(4) &
        -6.79432730d-15 * tc(5) &
        +5.09259997d+04 * invT
    !species 13: CH3
    species(13) = &
        +1.28571772d+00 &
        +3.61995018d-03 * tc(2) &
        -9.95714493d-07 * tc(3) &
        +1.48921161d-10 * tc(4) &
        -9.34308788d-15 * tc(5) &
        +1.67755843d+04 * invT
    !species 14: CH4
    species(14) = &
        -9.25148505d-01 &
        +6.69547335d-03 * tc(2) &
        -1.91095270d-06 * tc(3) &
        +3.05731338d-10 * tc(4) &
        -2.03630460d-14 * tc(5) &
        -9.46834459d+03 * invT
    !species 15: CO
    species(15) = &
        +1.71518561d+00 &
        +1.03126372d-03 * tc(2) &
        -3.32941924d-07 * tc(3) &
        +5.75132520d-11 * tc(4) &
        -4.07295432d-15 * tc(5) &
        -1.41518724d+04 * invT
    !species 16: CO2
    species(16) = &
        +2.85746029d+00 &
        +2.20718513d-03 * tc(2) &
        -7.38271347d-07 * tc(3) &
        +1.30872547d-10 * tc(4) &
        -9.44168328d-15 * tc(5) &
        -4.87591660d+04 * invT
    !species 17: HCO
    species(17) = &
        +1.77217438d+00 &
        +2.47847763d-03 * tc(2) &
        -8.28152043d-07 * tc(3) &
        +1.47290445d-10 * tc(4) &
        -1.06701742d-14 * tc(5) &
        +4.01191815d+03 * invT
    !species 18: CH2O
    species(18) = &
        +7.60690080d-01 &
        +4.60000041d-03 * tc(2) &
        -1.47419604d-06 * tc(3) &
        +2.51603030d-10 * tc(4) &
        -1.76771128d-14 * tc(5) &
        -1.39958323d+04 * invT
    !species 19: CH2OH
    species(19) = &
        +2.69266569d+00 &
        +4.32288399d-03 * tc(2) &
        -1.25033707d-06 * tc(3) &
        +1.96808659d-10 * tc(4) &
        -1.29710840d-14 * tc(5) &
        -3.24250627d+03 * invT
    !species 20: CH3O
    species(20) = &
        +2.77079900d+00 &
        +3.93574850d-03 * tc(2) &
        -8.85461333d-07 * tc(3) &
        +9.86107750d-11 * tc(4) &
        -4.22523200d-15 * tc(5) &
        +1.27832520d+02 * invT
    !species 21: CH3OH
    species(21) = &
        +7.89707910d-01 &
        +7.04691460d-03 * tc(2) &
        -2.12166945d-06 * tc(3) &
        +3.45427713d-10 * tc(4) &
        -2.34120440d-14 * tc(5) &
        -2.53748747d+04 * invT
    !species 22: C2H
    species(22) = &
        +2.16780652d+00 &
        +2.37610951d-03 * tc(2) &
        -6.12623590d-07 * tc(3) &
        +7.60475630d-11 * tc(4) &
        -3.54465540d-15 * tc(5) &
        +6.71210650d+04 * invT
    !species 23: C2H2
    species(23) = &
        +3.14756964d+00 &
        +2.98083332d-03 * tc(2) &
        -7.90982840d-07 * tc(3) &
        +1.16853043d-10 * tc(4) &
        -7.22470426d-15 * tc(5) &
        +2.59359992d+04 * invT
    !species 24: C2H3
    species(24) = &
        +2.01672400d+00 &
        +5.16511460d-03 * tc(2) &
        -1.56027450d-06 * tc(3) &
        +2.54408220d-10 * tc(4) &
        -1.72521408d-14 * tc(5) &
        +3.46128739d+04 * invT
    !species 25: C2H4
    species(25) = &
        +1.03611116d+00 &
        +7.32270755d-03 * tc(2) &
        -2.23692638d-06 * tc(3) &
        +3.68057308d-10 * tc(4) &
        -2.51412122d-14 * tc(5) &
        +4.93988614d+03 * invT
    !species 26: C2H5
    species(26) = &
        +9.54656420d-01 &
        +8.69863610d-03 * tc(2) &
        -2.66068889d-06 * tc(3) &
        +4.38044223d-10 * tc(4) &
        -2.99283152d-14 * tc(5) &
        +1.28575200d+04 * invT
    !species 27: C2H6
    species(27) = &
        +7.18815000d-02 &
        +1.08426339d-02 * tc(2) &
        -3.34186890d-06 * tc(3) &
        +5.53530003d-10 * tc(4) &
        -3.80005780d-14 * tc(5) &
        -1.14263932d+04 * invT
    !species 28: HCCO
    species(28) = &
        +4.62820580d+00 &
        +2.04267005d-03 * tc(2) &
        -5.31151567d-07 * tc(3) &
        +7.15651300d-11 * tc(4) &
        -3.88156640d-15 * tc(5) &
        +1.93272150d+04 * invT
    !species 29: CH2CO
    species(29) = &
        +3.51129732d+00 &
        +4.50179872d-03 * tc(2) &
        -1.38979878d-06 * tc(3) &
        +2.30836470d-10 * tc(4) &
        -1.58967640d-14 * tc(5) &
        -7.55105311d+03 * invT
    !species 30: HCCOH
    species(30) = &
        +4.92382910d+00 &
        +3.39618000d-03 * tc(2) &
        -8.55285467d-07 * tc(3) &
        +1.12469603d-10 * tc(4) &
        -5.98802020d-15 * tc(5) &
        +7.26462600d+03 * invT
    !species 31: N
    species(31) = &
        +1.41594290d+00 &
        +8.74453250d-05 * tc(2) &
        -3.96745633d-08 * tc(3) &
        +7.55656125d-12 * tc(4) &
        -4.07219640d-16 * tc(5) &
        +5.61337730d+04 * invT
    !species 32: NH
    species(32) = &
        +1.78369280d+00 &
        +6.64921500d-04 * tc(2) &
        -1.41593490d-07 * tc(3) &
        +1.95871253d-11 * tc(4) &
        -1.10088940d-15 * tc(5) &
        +4.21208480d+04 * invT
    !species 33: NH2
    species(33) = &
        +1.83474210d+00 &
        +1.60365410d-03 * tc(2) &
        -3.11302680d-07 * tc(3) &
        +3.42573825d-11 * tc(4) &
        -1.58412288d-15 * tc(5) &
        +2.21719570d+04 * invT
    !species 34: NH3
    species(34) = &
        +1.63445210d+00 &
        +2.83312800d-03 * tc(2) &
        -5.75955867d-07 * tc(3) &
        +5.96679025d-11 * tc(4) &
        -2.51575720d-15 * tc(5) &
        -6.54469580d+03 * invT
    !species 35: NNH
    species(35) = &
        +2.76675440d+00 &
        +1.44575410d-03 * tc(2) &
        -3.47220667d-07 * tc(3) &
        +4.21064850d-11 * tc(4) &
        -2.01837920d-15 * tc(5) &
        +2.86506970d+04 * invT
    !species 36: NO
    species(36) = &
        +2.26060560d+00 &
        +5.95552150d-04 * tc(2) &
        -1.43056827d-07 * tc(3) &
        +1.73644173d-11 * tc(4) &
        -8.06721980d-16 * tc(5) &
        +9.92097460d+03 * invT
    !species 37: NO2
    species(37) = &
        +3.88475420d+00 &
        +1.08619780d-03 * tc(2) &
        -2.76023020d-07 * tc(3) &
        +3.93687750d-11 * tc(4) &
        -2.10217900d-15 * tc(5) &
        +2.31649830d+03 * invT
    !species 38: N2O
    species(38) = &
        +3.82307290d+00 &
        +1.31351255d-03 * tc(2) &
        -3.19502913d-07 * tc(3) &
        +4.00017800d-11 * tc(4) &
        -1.95504606d-15 * tc(5) &
        +8.07340480d+03 * invT
    !species 39: HNO
    species(39) = &
        +1.97925090d+00 &
        +1.74720295d-03 * tc(2) &
        -2.61832593d-07 * tc(3) &
        +1.43698985d-11 * tc(4) &
        -3.86718320d-17 * tc(5) &
        +1.17505820d+04 * invT
    !species 40: CN
    species(40) = &
        +2.74598050d+00 &
        +2.17253875d-05 * tc(2) &
        +9.90199467d-08 * tc(3) &
        -1.71629515d-11 * tc(4) &
        +8.82683460d-16 * tc(5) &
        +5.15361880d+04 * invT
    !species 41: HCN
    species(41) = &
        +2.80223920d+00 &
        +1.57321140d-03 * tc(2) &
        -3.54406167d-07 * tc(3) &
        +4.15493925d-11 * tc(4) &
        -1.95995140d-15 * tc(5) &
        +1.44072920d+04 * invT
    !species 42: H2CN
    species(42) = &
        +4.20970300d+00 &
        +1.48464555d-03 * tc(2) &
        -9.51863033d-08 * tc(3) &
        -4.08887500d-11 * tc(4) &
        +6.08651780d-15 * tc(5) &
        +2.76771090d+04 * invT
    !species 43: HCNN
    species(43) = &
        +4.89463620d+00 &
        +1.99479795d-03 * tc(2) &
        -5.32746000d-07 * tc(3) &
        +7.31234875d-11 * tc(4) &
        -4.01893720d-15 * tc(5) &
        +5.34529410d+04 * invT
    !species 47: NCO
    species(47) = &
        +4.15218450d+00 &
        +1.15258805d-03 * tc(2) &
        -2.93443843d-07 * tc(3) &
        +3.69727450d-11 * tc(4) &
        -1.81955992d-15 * tc(5) &
        +1.40041230d+04 * invT
    !species 48: N2
    species(48) = &
        +1.92664000d+00 &
        +7.43988400d-04 * tc(2) &
        -1.89492000d-07 * tc(3) &
        +2.52425950d-11 * tc(4) &
        -1.35067020d-15 * tc(5) &
        -9.22797700d+02 * invT
    !species 49: AR
    species(49) = &
        +1.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5) &
        -7.45375000d+02 * invT
    !species 50: C3H7
    species(50) = &
        +6.70269870d+00 &
        +8.02210150d-03 * tc(2) &
        -1.76110733d-06 * tc(3) &
        +1.90746475d-10 * tc(4) &
        -7.87845680d-15 * tc(5) &
        +8.29843360d+03 * invT
    !species 51: C3H8
    species(51) = &
        +6.53413680d+00 &
        +9.43611950d-03 * tc(2) &
        -2.09061637d-06 * tc(3) &
        +2.28689123d-10 * tc(4) &
        -9.56761380d-15 * tc(5) &
        -1.64675160d+04 * invT
    !species 52: CH2CHO
    species(52) = &
        +4.97567000d+00 &
        +4.06529550d-03 * tc(2) &
        -9.14541333d-07 * tc(3) &
        +1.01757600d-10 * tc(4) &
        -4.35203400d-15 * tc(5) &
        +4.90321800d+02 * invT
    !species 53: CH3CHO
    species(53) = &
        +4.40411080d+00 &
        +5.86152950d-03 * tc(2) &
        -1.40877123d-06 * tc(3) &
        +1.70931128d-10 * tc(4) &
        -8.19697260d-15 * tc(5) &
        -2.25931220d+04 * invT
end if

! species with midpoint at T=1368 kelvin
if (T < 1368.000000d0) then
    ! species 45: HOCN
    species(45) = &
        +2.78604952d+00 &
        +3.44333961d-03 * tc(2) &
        -1.07162621d-06 * tc(3) &
        +1.29298942d-10 * tc(4) &
        +2.38721576d-15 * tc(5) &
        -2.82698400d+03 * invT
else
    !species 45: HOCN
    species(45) = &
        +4.89784885d+00 &
        +1.58394696d-03 * tc(2) &
        -3.72670213d-07 * tc(3) &
        +4.43107860d-11 * tc(4) &
        -2.08678354d-15 * tc(5) &
        -3.70653331d+03 * invT
end if

! species with midpoint at T=1478 kelvin
if (T < 1478.000000d0) then
    ! species 46: HNCO
    species(46) = &
        +2.63096317d+00 &
        +3.65141179d-03 * tc(2) &
        -7.60166677d-07 * tc(3) &
        -1.65317825d-10 * tc(4) &
        +7.24471504d-14 * tc(5) &
        -1.55873636d+04 * invT
else
    !species 46: HNCO
    species(46) = &
        +5.22395134d+00 &
        +1.58932002d-03 * tc(2) &
        -3.64595850d-07 * tc(3) &
        +4.26837908d-11 * tc(4) &
        -1.99004391d-15 * tc(5) &
        -1.66599344d+04 * invT
end if

! species with midpoint at T=1382 kelvin
if (T < 1382.000000d0) then
    ! species 44: HCNO
    species(44) = &
        +1.64727989d+00 &
        +6.37526710d-03 * tc(2) &
        -3.49314120d-06 * tc(3) &
        +1.10358209d-09 * tc(4) &
        -1.51504293d-13 * tc(5) &
        +1.92990252d+04 * invT
else
    !species 44: HCNO
    species(44) = &
        +5.59860456d+00 &
        +1.51389313d-03 * tc(2) &
        -3.59014487d-07 * tc(3) &
        +4.29166320d-11 * tc(4) &
        -2.02878782d-15 * tc(5) &
        +1.79661339d+04 * invT
end if

end subroutine

! compute the h/(RT) at the given temperature (Eq 20)
! tc contains precomputed powers of T, tc(1) = log(T)
subroutine speciesEnthalpy(species, tc)

implicit none

double precision, intent(out) :: species(53)
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
    ! species 8: H2O2
    species(8) = &
        +4.27611269d+00 &
        -2.71411208d-04 * tc(2) &
        +5.57785670d-06 * tc(3) &
        -5.39427032d-09 * tc(4) &
        +1.72490873d-12 * tc(5) &
        -1.77025821d+04 * invT
    ! species 9: C
    species(9) = &
        +2.55423955d+00 &
        -1.60768862d-04 * tc(2) &
        +2.44597415d-07 * tc(3) &
        -1.83058722d-10 * tc(4) &
        +5.33042892d-14 * tc(5) &
        +8.54438832d+04 * invT
    ! species 10: CH
    species(10) = &
        +3.48981665d+00 &
        +1.61917771d-04 * tc(2) &
        -5.62996883d-07 * tc(3) &
        +7.90543317d-10 * tc(4) &
        -2.81218134d-13 * tc(5) &
        +7.07972934d+04 * invT
    ! species 11: CH2
    species(11) = &
        +3.76267867d+00 &
        +4.84436072d-04 * tc(2) &
        +9.31632803d-07 * tc(3) &
        -9.62727883d-10 * tc(4) &
        +3.37483438d-13 * tc(5) &
        +4.60040401d+04 * invT
    ! species 12: CH2(S)
    species(12) = &
        +4.19860411d+00 &
        -1.18330710d-03 * tc(2) &
        +2.74432073d-06 * tc(3) &
        -1.67203995d-09 * tc(4) &
        +3.88629474d-13 * tc(5) &
        +5.04968163d+04 * invT
    ! species 13: CH3
    species(13) = &
        +3.67359040d+00 &
        +1.00547588d-03 * tc(2) &
        +1.91007285d-06 * tc(3) &
        -1.71779356d-09 * tc(4) &
        +5.08771468d-13 * tc(5) &
        +1.64449988d+04 * invT
    ! species 14: CH4
    species(14) = &
        +5.14987613d+00 &
        -6.83548940d-03 * tc(2) &
        +1.63933533d-05 * tc(3) &
        -1.21185757d-08 * tc(4) &
        +3.33387912d-12 * tc(5) &
        -1.02466476d+04 * invT
    ! species 15: CO
    species(15) = &
        +3.57953347d+00 &
        -3.05176840d-04 * tc(2) &
        +3.38938110d-07 * tc(3) &
        +2.26751471d-10 * tc(4) &
        -1.80884900d-13 * tc(5) &
        -1.43440860d+04 * invT
    ! species 16: CO2
    species(16) = &
        +2.35677352d+00 &
        +4.49229839d-03 * tc(2) &
        -2.37452090d-06 * tc(3) &
        +6.14797555d-10 * tc(4) &
        -2.87399096d-14 * tc(5) &
        -4.83719697d+04 * invT
    ! species 17: HCO
    species(17) = &
        +4.22118584d+00 &
        -1.62196266d-03 * tc(2) &
        +4.59331487d-06 * tc(3) &
        -3.32860233d-09 * tc(4) &
        +8.67537730d-13 * tc(5) &
        +3.83956496d+03 * invT
    ! species 18: CH2O
    species(18) = &
        +4.79372315d+00 &
        -4.95416684d-03 * tc(2) &
        +1.24406669d-05 * tc(3) &
        -9.48213152d-09 * tc(4) &
        +2.63545304d-12 * tc(5) &
        -1.43089567d+04 * invT
    ! species 19: CH2OH
    species(19) = &
        +3.86388918d+00 &
        +2.79836152d-03 * tc(2) &
        +1.97757264d-06 * tc(3) &
        -2.61330030d-09 * tc(4) &
        +8.73934556d-13 * tc(5) &
        -3.19391367d+03 * invT
    ! species 20: CH3O
    species(20) = &
        +2.10620400d+00 &
        +3.60829750d-03 * tc(2) &
        +1.77949067d-06 * tc(3) &
        -1.84440900d-09 * tc(4) &
        +4.15122000d-13 * tc(5) &
        +9.78601100d+02 * invT
    ! species 21: CH3OH
    species(21) = &
        +5.71539582d+00 &
        -7.61545645d-03 * tc(2) &
        +2.17480385d-05 * tc(3) &
        -1.77701722d-08 * tc(4) &
        +5.22705396d-12 * tc(5) &
        -2.56427656d+04 * invT
    ! species 22: C2H
    species(22) = &
        +2.88965733d+00 &
        +6.70498055d-03 * tc(2) &
        -9.49231670d-06 * tc(3) &
        +7.36977613d-09 * tc(4) &
        -2.18663022d-12 * tc(5) &
        +6.68393932d+04 * invT
    ! species 23: C2H2
    species(23) = &
        +8.08681094d-01 &
        +1.16807815d-02 * tc(2) &
        -1.18390605d-05 * tc(3) &
        +7.00381092d-09 * tc(4) &
        -1.70014595d-12 * tc(5) &
        +2.64289807d+04 * invT
    ! species 24: C2H3
    species(24) = &
        +3.21246645d+00 &
        +7.57395810d-04 * tc(2) &
        +8.64031373d-06 * tc(3) &
        -8.94144617d-09 * tc(4) &
        +2.94301746d-12 * tc(5) &
        +3.48598468d+04 * invT
    ! species 25: C2H4
    species(25) = &
        +3.95920148d+00 &
        -3.78526124d-03 * tc(2) &
        +1.90330097d-05 * tc(3) &
        -1.72897188d-08 * tc(4) &
        +5.39768746d-12 * tc(5) &
        +5.08977593d+03 * invT
    ! species 26: C2H5
    species(26) = &
        +4.30646568d+00 &
        -2.09329446d-03 * tc(2) &
        +1.65714269d-05 * tc(3) &
        -1.49781651d-08 * tc(4) &
        +4.61018008d-12 * tc(5) &
        +1.28416265d+04 * invT
    ! species 27: C2H6
    species(27) = &
        +4.29142492d+00 &
        -2.75077135d-03 * tc(2) &
        +1.99812763d-05 * tc(3) &
        -1.77116571d-08 * tc(4) &
        +5.37371542d-12 * tc(5) &
        -1.15222055d+04 * invT
    ! species 28: HCCO
    species(28) = &
        +2.25172140d+00 &
        +8.82751050d-03 * tc(2) &
        -7.90970033d-06 * tc(3) &
        +4.31893975d-09 * tc(4) &
        -1.01329622d-12 * tc(5) &
        +2.00594490d+04 * invT
    ! species 29: CH2CO
    species(29) = &
        +2.13583630d+00 &
        +9.05943605d-03 * tc(2) &
        -5.79824913d-06 * tc(3) &
        +2.33599392d-09 * tc(4) &
        -4.02915230d-13 * tc(5) &
        -7.04291804d+03 * invT
    ! species 30: HCCOH
    species(30) = &
        +1.24237330d+00 &
        +1.55361005d-02 * tc(2) &
        -1.69556213d-05 * tc(3) &
        +1.07842828d-08 * tc(4) &
        -2.80291880d-12 * tc(5) &
        +8.03161430d+03 * invT
    ! species 31: N
    species(31) = &
        +2.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5) &
        +5.61046370d+04 * invT
    ! species 32: NH
    species(32) = &
        +3.49290850d+00 &
        +1.55895990d-04 * tc(2) &
        -4.96349467d-07 * tc(3) &
        +6.20411050d-10 * tc(4) &
        -2.07139340d-13 * tc(5) &
        +4.18806290d+04 * invT
    ! species 33: NH2
    species(33) = &
        +4.20400290d+00 &
        -1.05306925d-03 * tc(2) &
        +2.36894493d-06 * tc(3) &
        -1.40287992d-09 * tc(4) &
        +3.28814340d-13 * tc(5) &
        +2.18859100d+04 * invT
    ! species 34: NH3
    species(34) = &
        +4.28602740d+00 &
        -2.33026150d-03 * tc(2) &
        +7.23950433d-06 * tc(3) &
        -5.70222175d-09 * tc(4) &
        +1.65276092d-12 * tc(5) &
        -6.74172850d+03 * invT
    ! species 35: NNH
    species(35) = &
        +4.34469270d+00 &
        -2.42485360d-03 * tc(2) &
        +6.68648633d-06 * tc(3) &
        -5.43161600d-09 * tc(4) &
        +1.58939078d-12 * tc(5) &
        +2.87919730d+04 * invT
    ! species 36: NO
    species(36) = &
        +4.21847630d+00 &
        -2.31948800d-03 * tc(2) &
        +3.68034067d-06 * tc(3) &
        -2.33403385d-09 * tc(4) &
        +5.60715400d-13 * tc(5) &
        +9.84462300d+03 * invT
    ! species 37: NO2
    species(37) = &
        +3.94403120d+00 &
        -7.92714500d-04 * tc(2) &
        +5.55260400d-06 * tc(3) &
        -5.11885650d-09 * tc(4) &
        +1.56701128d-12 * tc(5) &
        +2.89661790d+03 * invT
    ! species 38: N2O
    species(38) = &
        +2.25715020d+00 &
        +5.65236400d-03 * tc(2) &
        -4.55710633d-06 * tc(3) &
        +2.42049515d-09 * tc(4) &
        -5.86143640d-13 * tc(5) &
        +8.74177440d+03 * invT
    ! species 39: HNO
    species(39) = &
        +4.53349160d+00 &
        -2.83480855d-03 * tc(2) &
        +6.15773567d-06 * tc(3) &
        -4.28427350d-09 * tc(4) &
        +1.10909146d-12 * tc(5) &
        +1.15482970d+04 * invT
    ! species 40: CN
    species(40) = &
        +3.61293510d+00 &
        -4.77756635d-04 * tc(2) &
        +7.14765900d-07 * tc(3) &
        -7.87908075d-11 * tc(4) &
        -9.28607120d-14 * tc(5) &
        +5.17083400d+04 * invT
    ! species 41: HCN
    species(41) = &
        +2.25898860d+00 &
        +5.02558500d-03 * tc(2) &
        -4.45058767d-06 * tc(3) &
        +2.52308725d-09 * tc(4) &
        -6.01780560d-13 * tc(5) &
        +1.47126330d+04 * invT
    ! species 42: H2CN
    species(42) = &
        +2.85166100d+00 &
        +2.84761655d-03 * tc(2) &
        +3.57046667d-07 * tc(3) &
        -4.05653000d-10 * tc(4) &
        -4.70221620d-14 * tc(5) &
        +2.86378200d+04 * invT
    ! species 43: HCNN
    species(43) = &
        +2.52431940d+00 &
        +7.98030950d-03 * tc(2) &
        -6.27211800d-06 * tc(3) &
        +3.03138500d-09 * tc(4) &
        -6.47147560d-13 * tc(5) &
        +5.42619840d+04 * invT
    ! species 47: NCO
    species(47) = &
        +2.82693080d+00 &
        +4.40258440d-03 * tc(2) &
        -2.79553780d-06 * tc(3) &
        +1.20042410d-09 * tc(4) &
        -2.66271900d-13 * tc(5) &
        +1.46824770d+04 * invT
    ! species 48: N2
    species(48) = &
        +3.29867700d+00 &
        +7.04120200d-04 * tc(2) &
        -1.32107400d-06 * tc(3) &
        +1.41037875d-09 * tc(4) &
        -4.88970800d-13 * tc(5) &
        -1.02089990d+03 * invT
    ! species 49: AR
    species(49) = &
        +2.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5) &
        -7.45375000d+02 * invT
    ! species 50: C3H7
    species(50) = &
        +1.05155180d+00 &
        +1.29959900d-02 * tc(2) &
        +7.93351333d-07 * tc(3) &
        -4.90239225d-09 * tc(4) &
        +1.87464940d-12 * tc(5) &
        +1.06318630d+04 * invT
    ! species 51: C3H8
    species(51) = &
        +9.33553810d-01 &
        +1.32122895d-02 * tc(2) &
        +2.03532423d-06 * tc(3) &
        -5.49437475d-09 * tc(4) &
        +1.90298506d-12 * tc(5) &
        -1.39585200d+04 * invT
    ! species 52: CH2CHO
    species(52) = &
        +3.40906200d+00 &
        +5.36928700d-03 * tc(2) &
        +6.30497333d-07 * tc(3) &
        -1.78964575d-09 * tc(4) &
        +5.73477000d-13 * tc(5) &
        +1.52147660d+03 * invT
    ! species 53: CH3CHO
    species(53) = &
        +4.72945950d+00 &
        -1.59664290d-03 * tc(2) &
        +1.58449737d-05 * tc(3) &
        -1.43646527d-08 * tc(4) &
        +4.38622240d-12 * tc(5) &
        -2.15728780d+04 * invT
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
    !species 8: H2O2
    species(8) = &
        +4.16500285d+00 &
        +2.45415847d-03 * tc(2) &
        -6.33797417d-07 * tc(3) &
        +9.27964965d-11 * tc(4) &
        -5.75816610d-15 * tc(5) &
        -1.78617877d+04 * invT
    !species 9: C
    species(9) = &
        +2.49266888d+00 &
        +2.39944642d-05 * tc(2) &
        -2.41445007d-08 * tc(3) &
        +9.35727573d-12 * tc(4) &
        -9.74555786d-16 * tc(5) &
        +8.54512953d+04 * invT
    !species 10: CH
    species(10) = &
        +2.87846473d+00 &
        +4.85456840d-04 * tc(2) &
        +4.81485517d-08 * tc(3) &
        -3.26719623d-11 * tc(4) &
        +3.52158766d-15 * tc(5) &
        +7.10124364d+04 * invT
    !species 11: CH2
    species(11) = &
        +2.87410113d+00 &
        +1.82819646d-03 * tc(2) &
        -4.69648657d-07 * tc(3) &
        +6.50448872d-11 * tc(4) &
        -3.75455134d-15 * tc(5) &
        +4.62636040d+04 * invT
    !species 12: CH2(S)
    species(12) = &
        +2.29203842d+00 &
        +2.32794318d-03 * tc(2) &
        -6.70639823d-07 * tc(3) &
        +1.04476500d-10 * tc(4) &
        -6.79432730d-15 * tc(5) &
        +5.09259997d+04 * invT
    !species 13: CH3
    species(13) = &
        +2.28571772d+00 &
        +3.61995018d-03 * tc(2) &
        -9.95714493d-07 * tc(3) &
        +1.48921161d-10 * tc(4) &
        -9.34308788d-15 * tc(5) &
        +1.67755843d+04 * invT
    !species 14: CH4
    species(14) = &
        +7.48514950d-02 &
        +6.69547335d-03 * tc(2) &
        -1.91095270d-06 * tc(3) &
        +3.05731338d-10 * tc(4) &
        -2.03630460d-14 * tc(5) &
        -9.46834459d+03 * invT
    !species 15: CO
    species(15) = &
        +2.71518561d+00 &
        +1.03126372d-03 * tc(2) &
        -3.32941924d-07 * tc(3) &
        +5.75132520d-11 * tc(4) &
        -4.07295432d-15 * tc(5) &
        -1.41518724d+04 * invT
    !species 16: CO2
    species(16) = &
        +3.85746029d+00 &
        +2.20718513d-03 * tc(2) &
        -7.38271347d-07 * tc(3) &
        +1.30872547d-10 * tc(4) &
        -9.44168328d-15 * tc(5) &
        -4.87591660d+04 * invT
    !species 17: HCO
    species(17) = &
        +2.77217438d+00 &
        +2.47847763d-03 * tc(2) &
        -8.28152043d-07 * tc(3) &
        +1.47290445d-10 * tc(4) &
        -1.06701742d-14 * tc(5) &
        +4.01191815d+03 * invT
    !species 18: CH2O
    species(18) = &
        +1.76069008d+00 &
        +4.60000041d-03 * tc(2) &
        -1.47419604d-06 * tc(3) &
        +2.51603030d-10 * tc(4) &
        -1.76771128d-14 * tc(5) &
        -1.39958323d+04 * invT
    !species 19: CH2OH
    species(19) = &
        +3.69266569d+00 &
        +4.32288399d-03 * tc(2) &
        -1.25033707d-06 * tc(3) &
        +1.96808659d-10 * tc(4) &
        -1.29710840d-14 * tc(5) &
        -3.24250627d+03 * invT
    !species 20: CH3O
    species(20) = &
        +3.77079900d+00 &
        +3.93574850d-03 * tc(2) &
        -8.85461333d-07 * tc(3) &
        +9.86107750d-11 * tc(4) &
        -4.22523200d-15 * tc(5) &
        +1.27832520d+02 * invT
    !species 21: CH3OH
    species(21) = &
        +1.78970791d+00 &
        +7.04691460d-03 * tc(2) &
        -2.12166945d-06 * tc(3) &
        +3.45427713d-10 * tc(4) &
        -2.34120440d-14 * tc(5) &
        -2.53748747d+04 * invT
    !species 22: C2H
    species(22) = &
        +3.16780652d+00 &
        +2.37610951d-03 * tc(2) &
        -6.12623590d-07 * tc(3) &
        +7.60475630d-11 * tc(4) &
        -3.54465540d-15 * tc(5) &
        +6.71210650d+04 * invT
    !species 23: C2H2
    species(23) = &
        +4.14756964d+00 &
        +2.98083332d-03 * tc(2) &
        -7.90982840d-07 * tc(3) &
        +1.16853043d-10 * tc(4) &
        -7.22470426d-15 * tc(5) &
        +2.59359992d+04 * invT
    !species 24: C2H3
    species(24) = &
        +3.01672400d+00 &
        +5.16511460d-03 * tc(2) &
        -1.56027450d-06 * tc(3) &
        +2.54408220d-10 * tc(4) &
        -1.72521408d-14 * tc(5) &
        +3.46128739d+04 * invT
    !species 25: C2H4
    species(25) = &
        +2.03611116d+00 &
        +7.32270755d-03 * tc(2) &
        -2.23692638d-06 * tc(3) &
        +3.68057308d-10 * tc(4) &
        -2.51412122d-14 * tc(5) &
        +4.93988614d+03 * invT
    !species 26: C2H5
    species(26) = &
        +1.95465642d+00 &
        +8.69863610d-03 * tc(2) &
        -2.66068889d-06 * tc(3) &
        +4.38044223d-10 * tc(4) &
        -2.99283152d-14 * tc(5) &
        +1.28575200d+04 * invT
    !species 27: C2H6
    species(27) = &
        +1.07188150d+00 &
        +1.08426339d-02 * tc(2) &
        -3.34186890d-06 * tc(3) &
        +5.53530003d-10 * tc(4) &
        -3.80005780d-14 * tc(5) &
        -1.14263932d+04 * invT
    !species 28: HCCO
    species(28) = &
        +5.62820580d+00 &
        +2.04267005d-03 * tc(2) &
        -5.31151567d-07 * tc(3) &
        +7.15651300d-11 * tc(4) &
        -3.88156640d-15 * tc(5) &
        +1.93272150d+04 * invT
    !species 29: CH2CO
    species(29) = &
        +4.51129732d+00 &
        +4.50179872d-03 * tc(2) &
        -1.38979878d-06 * tc(3) &
        +2.30836470d-10 * tc(4) &
        -1.58967640d-14 * tc(5) &
        -7.55105311d+03 * invT
    !species 30: HCCOH
    species(30) = &
        +5.92382910d+00 &
        +3.39618000d-03 * tc(2) &
        -8.55285467d-07 * tc(3) &
        +1.12469603d-10 * tc(4) &
        -5.98802020d-15 * tc(5) &
        +7.26462600d+03 * invT
    !species 31: N
    species(31) = &
        +2.41594290d+00 &
        +8.74453250d-05 * tc(2) &
        -3.96745633d-08 * tc(3) &
        +7.55656125d-12 * tc(4) &
        -4.07219640d-16 * tc(5) &
        +5.61337730d+04 * invT
    !species 32: NH
    species(32) = &
        +2.78369280d+00 &
        +6.64921500d-04 * tc(2) &
        -1.41593490d-07 * tc(3) &
        +1.95871253d-11 * tc(4) &
        -1.10088940d-15 * tc(5) &
        +4.21208480d+04 * invT
    !species 33: NH2
    species(33) = &
        +2.83474210d+00 &
        +1.60365410d-03 * tc(2) &
        -3.11302680d-07 * tc(3) &
        +3.42573825d-11 * tc(4) &
        -1.58412288d-15 * tc(5) &
        +2.21719570d+04 * invT
    !species 34: NH3
    species(34) = &
        +2.63445210d+00 &
        +2.83312800d-03 * tc(2) &
        -5.75955867d-07 * tc(3) &
        +5.96679025d-11 * tc(4) &
        -2.51575720d-15 * tc(5) &
        -6.54469580d+03 * invT
    !species 35: NNH
    species(35) = &
        +3.76675440d+00 &
        +1.44575410d-03 * tc(2) &
        -3.47220667d-07 * tc(3) &
        +4.21064850d-11 * tc(4) &
        -2.01837920d-15 * tc(5) &
        +2.86506970d+04 * invT
    !species 36: NO
    species(36) = &
        +3.26060560d+00 &
        +5.95552150d-04 * tc(2) &
        -1.43056827d-07 * tc(3) &
        +1.73644173d-11 * tc(4) &
        -8.06721980d-16 * tc(5) &
        +9.92097460d+03 * invT
    !species 37: NO2
    species(37) = &
        +4.88475420d+00 &
        +1.08619780d-03 * tc(2) &
        -2.76023020d-07 * tc(3) &
        +3.93687750d-11 * tc(4) &
        -2.10217900d-15 * tc(5) &
        +2.31649830d+03 * invT
    !species 38: N2O
    species(38) = &
        +4.82307290d+00 &
        +1.31351255d-03 * tc(2) &
        -3.19502913d-07 * tc(3) &
        +4.00017800d-11 * tc(4) &
        -1.95504606d-15 * tc(5) &
        +8.07340480d+03 * invT
    !species 39: HNO
    species(39) = &
        +2.97925090d+00 &
        +1.74720295d-03 * tc(2) &
        -2.61832593d-07 * tc(3) &
        +1.43698985d-11 * tc(4) &
        -3.86718320d-17 * tc(5) &
        +1.17505820d+04 * invT
    !species 40: CN
    species(40) = &
        +3.74598050d+00 &
        +2.17253875d-05 * tc(2) &
        +9.90199467d-08 * tc(3) &
        -1.71629515d-11 * tc(4) &
        +8.82683460d-16 * tc(5) &
        +5.15361880d+04 * invT
    !species 41: HCN
    species(41) = &
        +3.80223920d+00 &
        +1.57321140d-03 * tc(2) &
        -3.54406167d-07 * tc(3) &
        +4.15493925d-11 * tc(4) &
        -1.95995140d-15 * tc(5) &
        +1.44072920d+04 * invT
    !species 42: H2CN
    species(42) = &
        +5.20970300d+00 &
        +1.48464555d-03 * tc(2) &
        -9.51863033d-08 * tc(3) &
        -4.08887500d-11 * tc(4) &
        +6.08651780d-15 * tc(5) &
        +2.76771090d+04 * invT
    !species 43: HCNN
    species(43) = &
        +5.89463620d+00 &
        +1.99479795d-03 * tc(2) &
        -5.32746000d-07 * tc(3) &
        +7.31234875d-11 * tc(4) &
        -4.01893720d-15 * tc(5) &
        +5.34529410d+04 * invT
    !species 47: NCO
    species(47) = &
        +5.15218450d+00 &
        +1.15258805d-03 * tc(2) &
        -2.93443843d-07 * tc(3) &
        +3.69727450d-11 * tc(4) &
        -1.81955992d-15 * tc(5) &
        +1.40041230d+04 * invT
    !species 48: N2
    species(48) = &
        +2.92664000d+00 &
        +7.43988400d-04 * tc(2) &
        -1.89492000d-07 * tc(3) &
        +2.52425950d-11 * tc(4) &
        -1.35067020d-15 * tc(5) &
        -9.22797700d+02 * invT
    !species 49: AR
    species(49) = &
        +2.50000000d+00 &
        +0.00000000d+00 * tc(2) &
        +0.00000000d+00 * tc(3) &
        +0.00000000d+00 * tc(4) &
        +0.00000000d+00 * tc(5) &
        -7.45375000d+02 * invT
    !species 50: C3H7
    species(50) = &
        +7.70269870d+00 &
        +8.02210150d-03 * tc(2) &
        -1.76110733d-06 * tc(3) &
        +1.90746475d-10 * tc(4) &
        -7.87845680d-15 * tc(5) &
        +8.29843360d+03 * invT
    !species 51: C3H8
    species(51) = &
        +7.53413680d+00 &
        +9.43611950d-03 * tc(2) &
        -2.09061637d-06 * tc(3) &
        +2.28689123d-10 * tc(4) &
        -9.56761380d-15 * tc(5) &
        -1.64675160d+04 * invT
    !species 52: CH2CHO
    species(52) = &
        +5.97567000d+00 &
        +4.06529550d-03 * tc(2) &
        -9.14541333d-07 * tc(3) &
        +1.01757600d-10 * tc(4) &
        -4.35203400d-15 * tc(5) &
        +4.90321800d+02 * invT
    !species 53: CH3CHO
    species(53) = &
        +5.40411080d+00 &
        +5.86152950d-03 * tc(2) &
        -1.40877123d-06 * tc(3) &
        +1.70931128d-10 * tc(4) &
        -8.19697260d-15 * tc(5) &
        -2.25931220d+04 * invT
end if

! species with midpoint at T=1368 kelvin
if (T < 1368.000000d0) then
    ! species 45: HOCN
    species(45) = &
        +3.78604952d+00 &
        +3.44333961d-03 * tc(2) &
        -1.07162621d-06 * tc(3) &
        +1.29298942d-10 * tc(4) &
        +2.38721576d-15 * tc(5) &
        -2.82698400d+03 * invT
else
    !species 45: HOCN
    species(45) = &
        +5.89784885d+00 &
        +1.58394696d-03 * tc(2) &
        -3.72670213d-07 * tc(3) &
        +4.43107860d-11 * tc(4) &
        -2.08678354d-15 * tc(5) &
        -3.70653331d+03 * invT
end if

! species with midpoint at T=1478 kelvin
if (T < 1478.000000d0) then
    ! species 46: HNCO
    species(46) = &
        +3.63096317d+00 &
        +3.65141179d-03 * tc(2) &
        -7.60166677d-07 * tc(3) &
        -1.65317825d-10 * tc(4) &
        +7.24471504d-14 * tc(5) &
        -1.55873636d+04 * invT
else
    !species 46: HNCO
    species(46) = &
        +6.22395134d+00 &
        +1.58932002d-03 * tc(2) &
        -3.64595850d-07 * tc(3) &
        +4.26837908d-11 * tc(4) &
        -1.99004391d-15 * tc(5) &
        -1.66599344d+04 * invT
end if

! species with midpoint at T=1382 kelvin
if (T < 1382.000000d0) then
    ! species 44: HCNO
    species(44) = &
        +2.64727989d+00 &
        +6.37526710d-03 * tc(2) &
        -3.49314120d-06 * tc(3) &
        +1.10358209d-09 * tc(4) &
        -1.51504293d-13 * tc(5) &
        +1.92990252d+04 * invT
else
    !species 44: HCNO
    species(44) = &
        +6.59860456d+00 &
        +1.51389313d-03 * tc(2) &
        -3.59014487d-07 * tc(3) &
        +4.29166320d-11 * tc(4) &
        -2.02878782d-15 * tc(5) &
        +1.79661339d+04 * invT
end if

end subroutine

! save molecular weights into array
subroutine molecularWeight(wt)

implicit none

double precision, intent(out) :: wt(53)

wt(1) = 2.015940d0 ! H2
wt(2) = 1.007970d0 ! H
wt(3) = 15.999400d0 ! O
wt(4) = 31.998800d0 ! O2
wt(5) = 17.007370d0 ! OH
wt(6) = 18.015340d0 ! H2O
wt(7) = 33.006770d0 ! HO2
wt(8) = 34.014740d0 ! H2O2
wt(9) = 12.011150d0 ! C
wt(10) = 13.019120d0 ! CH
wt(11) = 14.027090d0 ! CH2
wt(12) = 14.027090d0 ! CH2(S)
wt(13) = 15.035060d0 ! CH3
wt(14) = 16.043030d0 ! CH4
wt(15) = 28.010550d0 ! CO
wt(16) = 44.009950d0 ! CO2
wt(17) = 29.018520d0 ! HCO
wt(18) = 30.026490d0 ! CH2O
wt(19) = 31.034460d0 ! CH2OH
wt(20) = 31.034460d0 ! CH3O
wt(21) = 32.042430d0 ! CH3OH
wt(22) = 25.030270d0 ! C2H
wt(23) = 26.038240d0 ! C2H2
wt(24) = 27.046210d0 ! C2H3
wt(25) = 28.054180d0 ! C2H4
wt(26) = 29.062150d0 ! C2H5
wt(27) = 30.070120d0 ! C2H6
wt(28) = 41.029670d0 ! HCCO
wt(29) = 42.037640d0 ! CH2CO
wt(30) = 42.037640d0 ! HCCOH
wt(31) = 14.006700d0 ! N
wt(32) = 15.014670d0 ! NH
wt(33) = 16.022640d0 ! NH2
wt(34) = 17.030610d0 ! NH3
wt(35) = 29.021370d0 ! NNH
wt(36) = 30.006100d0 ! NO
wt(37) = 46.005500d0 ! NO2
wt(38) = 44.012800d0 ! N2O
wt(39) = 31.014070d0 ! HNO
wt(40) = 26.017850d0 ! CN
wt(41) = 27.025820d0 ! HCN
wt(42) = 28.033790d0 ! H2CN
wt(43) = 41.032520d0 ! HCNN
wt(44) = 43.025220d0 ! HCNO
wt(45) = 43.025220d0 ! HOCN
wt(46) = 43.025220d0 ! HNCO
wt(47) = 42.017250d0 ! NCO
wt(48) = 28.013400d0 ! N2
wt(49) = 39.948000d0 ! AR
wt(50) = 43.089240d0 ! C3H7
wt(51) = 44.097210d0 ! C3H8
wt(52) = 43.045610d0 ! CH2CHO
wt(53) = 44.053580d0 ! CH3CHO

end subroutine

! get temperature given internal energy in mass units and mass fracs
subroutine get_t_given_ey(e, y, iwrk, rwrk, t, ierr)

implicit none

double precision, intent(in) :: e
double precision, intent(in) :: y(53)
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

LENIMC = 214

end subroutine

subroutine egtransetLENRMC(LENRMC)

implicit none

integer, intent(out) :: LENRMC

LENRMC = 55226

end subroutine

subroutine egtransetNO(NO)

implicit none

integer, intent(out) :: NO

NO = 4

end subroutine

subroutine egtransetKK(KK)

implicit none

integer, intent(out) :: KK

KK = 53

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

double precision, intent(out) :: WT(53)

WT(1) = 2.01594000d+00
WT(2) = 1.00797000d+00
WT(3) = 1.59994000d+01
WT(4) = 3.19988000d+01
WT(5) = 1.70073700d+01
WT(6) = 1.80153400d+01
WT(7) = 3.30067700d+01
WT(8) = 3.40147400d+01
WT(9) = 1.20111500d+01
WT(10) = 1.30191200d+01
WT(11) = 1.40270900d+01
WT(12) = 1.40270900d+01
WT(13) = 1.50350600d+01
WT(14) = 1.60430300d+01
WT(15) = 2.80105500d+01
WT(16) = 4.40099500d+01
WT(17) = 2.90185200d+01
WT(18) = 3.00264900d+01
WT(19) = 3.10344600d+01
WT(20) = 3.10344600d+01
WT(21) = 3.20424300d+01
WT(22) = 2.50302700d+01
WT(23) = 2.60382400d+01
WT(24) = 2.70462100d+01
WT(25) = 2.80541800d+01
WT(26) = 2.90621500d+01
WT(27) = 3.00701200d+01
WT(28) = 4.10296700d+01
WT(29) = 4.20376400d+01
WT(30) = 4.20376400d+01
WT(31) = 1.40067000d+01
WT(32) = 1.50146700d+01
WT(33) = 1.60226400d+01
WT(34) = 1.70306100d+01
WT(35) = 2.90213700d+01
WT(36) = 3.00061000d+01
WT(37) = 4.60055000d+01
WT(38) = 4.40128000d+01
WT(39) = 3.10140700d+01
WT(40) = 2.60178500d+01
WT(41) = 2.70258200d+01
WT(42) = 2.80337900d+01
WT(43) = 4.10325200d+01
WT(44) = 4.30252200d+01
WT(45) = 4.30252200d+01
WT(46) = 4.30252200d+01
WT(47) = 4.20172500d+01
WT(48) = 2.80134000d+01
WT(49) = 3.99480000d+01
WT(50) = 4.30892400d+01
WT(51) = 4.40972100d+01
WT(52) = 4.30456100d+01
WT(53) = 4.40535800d+01

end subroutine

! the lennard-jones potential well depth eps/kb in K
subroutine egtransetEPS(EPS)

implicit none

double precision, intent(out) :: EPS(53)

EPS(25) = 2.80800000d+02
EPS(26) = 2.52300000d+02
EPS(37) = 2.00000000d+02
EPS(27) = 2.52300000d+02
EPS(29) = 4.36000000d+02
EPS(38) = 2.32400000d+02
EPS(28) = 1.50000000d+02
EPS(3) = 8.00000000d+01
EPS(42) = 5.69000000d+02
EPS(4) = 1.07400000d+02
EPS(41) = 5.69000000d+02
EPS(44) = 2.32400000d+02
EPS(39) = 1.16700000d+02
EPS(31) = 7.14000000d+01
EPS(32) = 8.00000000d+01
EPS(2) = 1.45000000d+02
EPS(52) = 4.36000000d+02
EPS(1) = 3.80000000d+01
EPS(48) = 9.75300000d+01
EPS(5) = 8.00000000d+01
EPS(53) = 4.36000000d+02
EPS(6) = 5.72400000d+02
EPS(45) = 2.32400000d+02
EPS(7) = 1.07400000d+02
EPS(50) = 2.66800000d+02
EPS(8) = 1.07400000d+02
EPS(33) = 8.00000000d+01
EPS(9) = 7.14000000d+01
EPS(10) = 8.00000000d+01
EPS(49) = 1.36500000d+02
EPS(11) = 1.44000000d+02
EPS(12) = 1.44000000d+02
EPS(46) = 2.32400000d+02
EPS(13) = 1.44000000d+02
EPS(14) = 1.41400000d+02
EPS(34) = 4.81000000d+02
EPS(15) = 9.81000000d+01
EPS(43) = 1.50000000d+02
EPS(16) = 2.44000000d+02
EPS(51) = 2.66800000d+02
EPS(17) = 4.98000000d+02
EPS(40) = 7.50000000d+01
EPS(18) = 4.98000000d+02
EPS(47) = 2.32400000d+02
EPS(19) = 4.17000000d+02
EPS(20) = 4.17000000d+02
EPS(36) = 9.75300000d+01
EPS(21) = 4.81800000d+02
EPS(22) = 2.09000000d+02
EPS(35) = 7.14000000d+01
EPS(23) = 2.09000000d+02
EPS(24) = 2.09000000d+02
EPS(30) = 4.36000000d+02

end subroutine

! the lennard-jones collision diameter in Angstroms
subroutine egtransetSIG(SIG)

implicit none

double precision, intent(out) :: SIG(53)

SIG(25) = 3.97100000d+00
SIG(26) = 4.30200000d+00
SIG(37) = 3.50000000d+00
SIG(27) = 4.30200000d+00
SIG(29) = 3.97000000d+00
SIG(38) = 3.82800000d+00
SIG(28) = 2.50000000d+00
SIG(3) = 2.75000000d+00
SIG(42) = 3.63000000d+00
SIG(4) = 3.45800000d+00
SIG(41) = 3.63000000d+00
SIG(44) = 3.82800000d+00
SIG(39) = 3.49200000d+00
SIG(31) = 3.29800000d+00
SIG(32) = 2.65000000d+00
SIG(2) = 2.05000000d+00
SIG(52) = 3.97000000d+00
SIG(1) = 2.92000000d+00
SIG(48) = 3.62100000d+00
SIG(5) = 2.75000000d+00
SIG(53) = 3.97000000d+00
SIG(6) = 2.60500000d+00
SIG(45) = 3.82800000d+00
SIG(7) = 3.45800000d+00
SIG(50) = 4.98200000d+00
SIG(8) = 3.45800000d+00
SIG(33) = 2.65000000d+00
SIG(9) = 3.29800000d+00
SIG(10) = 2.75000000d+00
SIG(49) = 3.33000000d+00
SIG(11) = 3.80000000d+00
SIG(12) = 3.80000000d+00
SIG(46) = 3.82800000d+00
SIG(13) = 3.80000000d+00
SIG(14) = 3.74600000d+00
SIG(34) = 2.92000000d+00
SIG(15) = 3.65000000d+00
SIG(43) = 2.50000000d+00
SIG(16) = 3.76300000d+00
SIG(51) = 4.98200000d+00
SIG(17) = 3.59000000d+00
SIG(40) = 3.85600000d+00
SIG(18) = 3.59000000d+00
SIG(47) = 3.82800000d+00
SIG(19) = 3.69000000d+00
SIG(20) = 3.69000000d+00
SIG(36) = 3.62100000d+00
SIG(21) = 3.62600000d+00
SIG(22) = 4.10000000d+00
SIG(35) = 3.79800000d+00
SIG(23) = 4.10000000d+00
SIG(24) = 4.10000000d+00
SIG(30) = 3.97000000d+00

end subroutine

! the dipole moment in Debye
subroutine egtransetDIP(DIP)

implicit none

double precision, intent(out) :: DIP(53)

DIP(25) = 0.00000000d+00
DIP(26) = 0.00000000d+00
DIP(37) = 0.00000000d+00
DIP(27) = 0.00000000d+00
DIP(29) = 0.00000000d+00
DIP(38) = 0.00000000d+00
DIP(28) = 0.00000000d+00
DIP(3) = 0.00000000d+00
DIP(42) = 0.00000000d+00
DIP(4) = 0.00000000d+00
DIP(41) = 0.00000000d+00
DIP(44) = 0.00000000d+00
DIP(39) = 0.00000000d+00
DIP(31) = 0.00000000d+00
DIP(32) = 0.00000000d+00
DIP(2) = 0.00000000d+00
DIP(52) = 0.00000000d+00
DIP(1) = 0.00000000d+00
DIP(48) = 0.00000000d+00
DIP(5) = 0.00000000d+00
DIP(53) = 0.00000000d+00
DIP(6) = 1.84400000d+00
DIP(45) = 0.00000000d+00
DIP(7) = 0.00000000d+00
DIP(50) = 0.00000000d+00
DIP(8) = 0.00000000d+00
DIP(33) = 0.00000000d+00
DIP(9) = 0.00000000d+00
DIP(10) = 0.00000000d+00
DIP(49) = 0.00000000d+00
DIP(11) = 0.00000000d+00
DIP(12) = 0.00000000d+00
DIP(46) = 0.00000000d+00
DIP(13) = 0.00000000d+00
DIP(14) = 0.00000000d+00
DIP(34) = 1.47000000d+00
DIP(15) = 0.00000000d+00
DIP(43) = 0.00000000d+00
DIP(16) = 0.00000000d+00
DIP(51) = 0.00000000d+00
DIP(17) = 0.00000000d+00
DIP(40) = 0.00000000d+00
DIP(18) = 0.00000000d+00
DIP(47) = 0.00000000d+00
DIP(19) = 1.70000000d+00
DIP(20) = 1.70000000d+00
DIP(36) = 0.00000000d+00
DIP(21) = 0.00000000d+00
DIP(22) = 0.00000000d+00
DIP(35) = 0.00000000d+00
DIP(23) = 0.00000000d+00
DIP(24) = 0.00000000d+00
DIP(30) = 0.00000000d+00

end subroutine

! the polarizability in cubic Angstroms
subroutine egtransetPOL(POL)

implicit none

double precision, intent(out) :: POL(53)

POL(25) = 0.00000000d+00
POL(26) = 0.00000000d+00
POL(37) = 0.00000000d+00
POL(27) = 0.00000000d+00
POL(29) = 0.00000000d+00
POL(38) = 0.00000000d+00
POL(28) = 0.00000000d+00
POL(3) = 0.00000000d+00
POL(42) = 0.00000000d+00
POL(4) = 1.60000000d+00
POL(41) = 0.00000000d+00
POL(44) = 0.00000000d+00
POL(39) = 0.00000000d+00
POL(31) = 0.00000000d+00
POL(32) = 0.00000000d+00
POL(2) = 0.00000000d+00
POL(52) = 0.00000000d+00
POL(1) = 7.90000000d-01
POL(48) = 1.76000000d+00
POL(5) = 0.00000000d+00
POL(53) = 0.00000000d+00
POL(6) = 0.00000000d+00
POL(45) = 0.00000000d+00
POL(7) = 0.00000000d+00
POL(50) = 0.00000000d+00
POL(8) = 0.00000000d+00
POL(33) = 2.26000000d+00
POL(9) = 0.00000000d+00
POL(10) = 0.00000000d+00
POL(49) = 0.00000000d+00
POL(11) = 0.00000000d+00
POL(12) = 0.00000000d+00
POL(46) = 0.00000000d+00
POL(13) = 0.00000000d+00
POL(14) = 2.60000000d+00
POL(34) = 0.00000000d+00
POL(15) = 1.95000000d+00
POL(43) = 0.00000000d+00
POL(16) = 2.65000000d+00
POL(51) = 0.00000000d+00
POL(17) = 0.00000000d+00
POL(40) = 0.00000000d+00
POL(18) = 0.00000000d+00
POL(47) = 0.00000000d+00
POL(19) = 0.00000000d+00
POL(20) = 0.00000000d+00
POL(36) = 1.76000000d+00
POL(21) = 0.00000000d+00
POL(22) = 0.00000000d+00
POL(35) = 0.00000000d+00
POL(23) = 0.00000000d+00
POL(24) = 0.00000000d+00
POL(30) = 0.00000000d+00

end subroutine

! the rotational relaxation collision number at 298 K
subroutine egtransetZROT(ZROT)

implicit none

double precision, intent(out) :: ZROT(53)

ZROT(25) = 1.50000000d+00
ZROT(26) = 1.50000000d+00
ZROT(37) = 1.00000000d+00
ZROT(27) = 1.50000000d+00
ZROT(29) = 2.00000000d+00
ZROT(38) = 1.00000000d+00
ZROT(28) = 1.00000000d+00
ZROT(3) = 0.00000000d+00
ZROT(42) = 1.00000000d+00
ZROT(4) = 3.80000000d+00
ZROT(41) = 1.00000000d+00
ZROT(44) = 1.00000000d+00
ZROT(39) = 1.00000000d+00
ZROT(31) = 0.00000000d+00
ZROT(32) = 4.00000000d+00
ZROT(2) = 0.00000000d+00
ZROT(52) = 2.00000000d+00
ZROT(1) = 2.80000000d+02
ZROT(48) = 4.00000000d+00
ZROT(5) = 0.00000000d+00
ZROT(53) = 2.00000000d+00
ZROT(6) = 4.00000000d+00
ZROT(45) = 1.00000000d+00
ZROT(7) = 1.00000000d+00
ZROT(50) = 1.00000000d+00
ZROT(8) = 3.80000000d+00
ZROT(33) = 4.00000000d+00
ZROT(9) = 0.00000000d+00
ZROT(10) = 0.00000000d+00
ZROT(49) = 0.00000000d+00
ZROT(11) = 0.00000000d+00
ZROT(12) = 0.00000000d+00
ZROT(46) = 1.00000000d+00
ZROT(13) = 0.00000000d+00
ZROT(14) = 1.30000000d+01
ZROT(34) = 1.00000000d+01
ZROT(15) = 1.80000000d+00
ZROT(43) = 1.00000000d+00
ZROT(16) = 2.10000000d+00
ZROT(51) = 1.00000000d+00
ZROT(17) = 0.00000000d+00
ZROT(40) = 1.00000000d+00
ZROT(18) = 2.00000000d+00
ZROT(47) = 1.00000000d+00
ZROT(19) = 2.00000000d+00
ZROT(20) = 2.00000000d+00
ZROT(36) = 4.00000000d+00
ZROT(21) = 1.00000000d+00
ZROT(22) = 2.50000000d+00
ZROT(35) = 1.00000000d+00
ZROT(23) = 2.50000000d+00
ZROT(24) = 1.00000000d+00
ZROT(30) = 2.00000000d+00

end subroutine

! 0: monoatomic, 1: linear, 2: nonlinear
subroutine egtransetNLIN(NLIN)

implicit none

integer, intent(out) :: NLIN(53)

NLIN(25) = 2
NLIN(26) = 2
NLIN(37) = 2
NLIN(27) = 2
NLIN(29) = 2
NLIN(38) = 1
NLIN(28) = 2
NLIN(3) = 0
NLIN(42) = 1
NLIN(4) = 1
NLIN(41) = 1
NLIN(44) = 2
NLIN(39) = 2
NLIN(31) = 0
NLIN(32) = 1
NLIN(2) = 0
NLIN(52) = 2
NLIN(1) = 1
NLIN(48) = 1
NLIN(5) = 1
NLIN(53) = 2
NLIN(6) = 2
NLIN(45) = 2
NLIN(7) = 2
NLIN(50) = 2
NLIN(8) = 2
NLIN(33) = 2
NLIN(9) = 0
NLIN(10) = 1
NLIN(49) = 0
NLIN(11) = 1
NLIN(12) = 1
NLIN(46) = 2
NLIN(13) = 1
NLIN(14) = 2
NLIN(34) = 2
NLIN(15) = 1
NLIN(43) = 2
NLIN(16) = 1
NLIN(51) = 2
NLIN(17) = 2
NLIN(40) = 1
NLIN(18) = 2
NLIN(47) = 1
NLIN(19) = 2
NLIN(20) = 2
NLIN(36) = 1
NLIN(21) = 2
NLIN(22) = 1
NLIN(35) = 2
NLIN(23) = 1
NLIN(24) = 2
NLIN(30) = 2

end subroutine


! Poly fits for the viscosities, dim NO*KK
subroutine egtransetCOFETA(COFETA)

implicit none

double precision, intent(out) :: COFETA(212)

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
COFETA(29) = -1.71312832d+01
COFETA(30) = 2.68036374d+00
COFETA(31) = -2.72570227d-01
COFETA(32) = 1.21650964d-02
COFETA(33) = -1.50067648d+01
COFETA(34) = 1.69625877d+00
COFETA(35) = -1.42936462d-01
COFETA(36) = 6.47223426d-03
COFETA(37) = -1.51956901d+01
COFETA(38) = 1.92606504d+00
COFETA(39) = -1.73487476d-01
COFETA(40) = 7.82572931d-03
COFETA(41) = -2.02663469d+01
COFETA(42) = 3.63241793d+00
COFETA(43) = -3.95581049d-01
COFETA(44) = 1.74725495d-02
COFETA(45) = -2.02663469d+01
COFETA(46) = 3.63241793d+00
COFETA(47) = -3.95581049d-01
COFETA(48) = 1.74725495d-02
COFETA(49) = -2.02316497d+01
COFETA(50) = 3.63241793d+00
COFETA(51) = -3.95581049d-01
COFETA(52) = 1.74725495d-02
COFETA(53) = -2.00094664d+01
COFETA(54) = 3.57220167d+00
COFETA(55) = -3.87936446d-01
COFETA(56) = 1.71483254d-02
COFETA(57) = -1.66188336d+01
COFETA(58) = 2.40307799d+00
COFETA(59) = -2.36167638d-01
COFETA(60) = 1.05714061d-02
COFETA(61) = -2.40014975d+01
COFETA(62) = 5.14359547d+00
COFETA(63) = -5.74269731d-01
COFETA(64) = 2.44937679d-02
COFETA(65) = -1.98501306d+01
COFETA(66) = 2.69480162d+00
COFETA(67) = -1.65880845d-01
COFETA(68) = 3.14504769d-03
COFETA(69) = -1.98330577d+01
COFETA(70) = 2.69480162d+00
COFETA(71) = -1.65880845d-01
COFETA(72) = 3.14504769d-03
COFETA(73) = -1.99945919d+01
COFETA(74) = 2.86923313d+00
COFETA(75) = -2.03325661d-01
COFETA(76) = 5.39056989d-03
COFETA(77) = -1.99945919d+01
COFETA(78) = 2.86923313d+00
COFETA(79) = -2.03325661d-01
COFETA(80) = 5.39056989d-03
COFETA(81) = -2.05644525d+01
COFETA(82) = 3.03946431d+00
COFETA(83) = -2.16994867d-01
COFETA(84) = 5.61394012d-03
COFETA(85) = -2.33863848d+01
COFETA(86) = 4.80350223d+00
COFETA(87) = -5.38341336d-01
COFETA(88) = 2.32747213d-02
COFETA(89) = -2.33666446d+01
COFETA(90) = 4.80350223d+00
COFETA(91) = -5.38341336d-01
COFETA(92) = 2.32747213d-02
COFETA(93) = -2.33476543d+01
COFETA(94) = 4.80350223d+00
COFETA(95) = -5.38341336d-01
COFETA(96) = 2.32747213d-02
COFETA(97) = -2.50655444d+01
COFETA(98) = 5.33982977d+00
COFETA(99) = -5.89982992d-01
COFETA(100) = 2.47780650d-02
COFETA(101) = -2.46581414d+01
COFETA(102) = 5.19497183d+00
COFETA(103) = -5.78827339d-01
COFETA(104) = 2.46050921d-02
COFETA(105) = -2.46410937d+01
COFETA(106) = 5.19497183d+00
COFETA(107) = -5.78827339d-01
COFETA(108) = 2.46050921d-02
COFETA(109) = -1.92183831d+01
COFETA(110) = 3.75164499d+00
COFETA(111) = -4.10390993d-01
COFETA(112) = 1.80861665d-02
COFETA(113) = -2.23395647d+01
COFETA(114) = 3.86433912d+00
COFETA(115) = -3.41553983d-01
COFETA(116) = 1.17083447d-02
COFETA(117) = -2.23395647d+01
COFETA(118) = 3.86433912d+00
COFETA(119) = -3.41553983d-01
COFETA(120) = 1.17083447d-02
COFETA(121) = -1.49299146d+01
COFETA(122) = 1.69625877d+00
COFETA(123) = -1.42936462d-01
COFETA(124) = 6.47223426d-03
COFETA(125) = -1.50503032d+01
COFETA(126) = 1.92606504d+00
COFETA(127) = -1.73487476d-01
COFETA(128) = 7.82572931d-03
COFETA(129) = -1.50178157d+01
COFETA(130) = 1.92606504d+00
COFETA(131) = -1.73487476d-01
COFETA(132) = 7.82572931d-03
COFETA(133) = -1.62012869d+01
COFETA(134) = 1.18120643d+00
COFETA(135) = 4.63881610d-02
COFETA(136) = -6.59115575d-03
COFETA(137) = -1.48479830d+01
COFETA(138) = 1.69625877d+00
COFETA(139) = -1.42936462d-01
COFETA(140) = 6.47223426d-03
COFETA(141) = -1.65352006d+01
COFETA(142) = 2.39056562d+00
COFETA(143) = -2.34558144d-01
COFETA(144) = 1.05024037d-02
COFETA(145) = -2.24069350d+01
COFETA(146) = 4.68007883d+00
COFETA(147) = -5.23982071d-01
COFETA(148) = 2.27205629d-02
COFETA(149) = -2.38016070d+01
COFETA(150) = 5.08149651d+00
COFETA(151) = -5.69446240d-01
COFETA(152) = 2.44183705d-02
COFETA(153) = -1.77852706d+01
COFETA(154) = 2.90334768d+00
COFETA(155) = -3.01373404d-01
COFETA(156) = 1.34067714d-02
COFETA(157) = -1.51559287d+01
COFETA(158) = 1.78185997d+00
COFETA(159) = -1.54248937d-01
COFETA(160) = 6.97039632d-03
COFETA(161) = -1.67652882d+01
COFETA(162) = 1.25414938d+00
COFETA(163) = 4.57994226d-02
COFETA(164) = -6.99199537d-03
COFETA(165) = -1.67469793d+01
COFETA(166) = 1.25414938d+00
COFETA(167) = 4.57994226d-02
COFETA(168) = -6.99199537d-03
COFETA(169) = -1.92183484d+01
COFETA(170) = 3.75164499d+00
COFETA(171) = -4.10390993d-01
COFETA(172) = 1.80861665d-02
COFETA(173) = -2.38129540d+01
COFETA(174) = 5.08149651d+00
COFETA(175) = -5.69446240d-01
COFETA(176) = 2.44183705d-02
COFETA(177) = -2.38129540d+01
COFETA(178) = 5.08149651d+00
COFETA(179) = -5.69446240d-01
COFETA(180) = 2.44183705d-02
COFETA(181) = -2.38129540d+01
COFETA(182) = 5.08149651d+00
COFETA(183) = -5.69446240d-01
COFETA(184) = 2.44183705d-02
COFETA(185) = -2.38248072d+01
COFETA(186) = 5.08149651d+00
COFETA(187) = -5.69446240d-01
COFETA(188) = 2.44183705d-02
COFETA(189) = -1.65695594d+01
COFETA(190) = 2.39056562d+00
COFETA(191) = -2.34558144d-01
COFETA(192) = 1.05024037d-02
COFETA(193) = -1.90422907d+01
COFETA(194) = 3.47025711d+00
COFETA(195) = -3.75102111d-01
COFETA(196) = 1.66086076d-02
COFETA(197) = -2.51288278d+01
COFETA(198) = 5.30723075d+00
COFETA(199) = -5.89742369d-01
COFETA(200) = 2.49294033d-02
COFETA(201) = -2.51172662d+01
COFETA(202) = 5.30723075d+00
COFETA(203) = -5.89742369d-01
COFETA(204) = 2.49294033d-02
COFETA(205) = -2.23277173d+01
COFETA(206) = 3.86433912d+00
COFETA(207) = -3.41553983d-01
COFETA(208) = 1.17083447d-02
COFETA(209) = -2.23161441d+01
COFETA(210) = 3.86433912d+00
COFETA(211) = -3.41553983d-01
COFETA(212) = 1.17083447d-02

end subroutine


! Poly fits for the conductivities, dim NO*KK
subroutine egtransetCOFLAM(COFLAM)

implicit none

double precision, intent(out) :: COFLAM(212)

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
COFLAM(29) = 8.83996545d-01
COFLAM(30) = 1.31525428d+00
COFLAM(31) = 1.91774576d-02
COFLAM(32) = -4.41642722d-03
COFLAM(33) = 2.06524872d+00
COFLAM(34) = 1.69625877d+00
COFLAM(35) = -1.42936462d-01
COFLAM(36) = 6.47223426d-03
COFLAM(37) = 2.08093824d+01
COFLAM(38) = -6.24180163d+00
COFLAM(39) = 9.82387427d-01
COFLAM(40) = -4.50360664d-02
COFLAM(41) = 1.29177902d+01
COFLAM(42) = -3.73745535d+00
COFLAM(43) = 7.15831021d-01
COFLAM(44) = -3.63846910d-02
COFLAM(45) = 1.89383266d+01
COFLAM(46) = -6.51018128d+00
COFLAM(47) = 1.13292060d+00
COFLAM(48) = -5.69603226d-02
COFLAM(49) = 1.39937901d+01
COFLAM(50) = -4.64256494d+00
COFLAM(51) = 9.07728674d-01
COFLAM(52) = -4.77274469d-02
COFLAM(53) = 1.33091602d+01
COFLAM(54) = -4.96140261d+00
COFLAM(55) = 1.03295505d+00
COFLAM(56) = -5.63420090d-02
COFLAM(57) = 1.18777264d+01
COFLAM(58) = -3.15463949d+00
COFLAM(59) = 6.01973268d-01
COFLAM(60) = -3.03211261d-02
COFLAM(61) = -1.13649314d+01
COFLAM(62) = 5.88177395d+00
COFLAM(63) = -5.68651819d-01
COFLAM(64) = 2.03561485d-02
COFLAM(65) = 6.30243508d+00
COFLAM(66) = -2.22810801d+00
COFLAM(67) = 6.37340514d-01
COFLAM(68) = -3.81056018d-02
COFLAM(69) = 5.39305086d+00
COFLAM(70) = -2.39312375d+00
COFLAM(71) = 7.39585221d-01
COFLAM(72) = -4.58435589d-02
COFLAM(73) = 1.05517076d+00
COFLAM(74) = 4.29805638d-02
COFLAM(75) = 3.31460301d-01
COFLAM(76) = -2.40803935d-02
COFLAM(77) = -6.14588187d+00
COFLAM(78) = 2.47428873d+00
COFLAM(79) = 6.43999571d-02
COFLAM(80) = -1.45368336d-02
COFLAM(81) = -1.83491607d+00
COFLAM(82) = 5.55823162d-01
COFLAM(83) = 3.51309919d-01
COFLAM(84) = -2.85527489d-02
COFLAM(85) = 1.33513522d+01
COFLAM(86) = -4.32126196d+00
COFLAM(87) = 8.36087620d-01
COFLAM(88) = -4.37693086d-02
COFLAM(89) = -7.70164502d+00
COFLAM(90) = 4.56884453d+00
COFLAM(91) = -4.04747583d-01
COFLAM(92) = 1.40841060d-02
COFLAM(93) = -9.10384516d+00
COFLAM(94) = 4.54798869d+00
COFLAM(95) = -3.18114458d-01
COFLAM(96) = 6.59577576d-03
COFLAM(97) = -1.46152839d+01
COFLAM(98) = 6.36251406d+00
COFLAM(99) = -5.03832130d-01
COFLAM(100) = 1.26121050d-02
COFLAM(101) = -8.95009705d+00
COFLAM(102) = 4.02515080d+00
COFLAM(103) = -1.84063946d-01
COFLAM(104) = -1.94054752d-03
COFLAM(105) = -1.09902209d+01
COFLAM(106) = 4.70647707d+00
COFLAM(107) = -2.52272495d-01
COFLAM(108) = 1.75193258d-04
COFLAM(109) = -5.96382813d+00
COFLAM(110) = 4.39356020d+00
COFLAM(111) = -3.97598260d-01
COFLAM(112) = 1.39749624d-02
COFLAM(113) = -8.32871231d+00
COFLAM(114) = 3.97067262d+00
COFLAM(115) = -2.21252287d-01
COFLAM(116) = 1.47870329d-03
COFLAM(117) = -5.36683560d+00
COFLAM(118) = 2.98311144d+00
COFLAM(119) = -1.14643230d-01
COFLAM(120) = -2.12830349d-03
COFLAM(121) = 1.98839852d+00
COFLAM(122) = 1.69625877d+00
COFLAM(123) = -1.42936462d-01
COFLAM(124) = 6.47223426d-03
COFLAM(125) = 1.17441953d+01
COFLAM(126) = -2.21994519d+00
COFLAM(127) = 4.01105510d-01
COFLAM(128) = -1.75558984d-02
COFLAM(129) = 1.97087241d+01
COFLAM(130) = -6.06759546d+00
COFLAM(131) = 1.00806444d+00
COFLAM(132) = -4.81231868d-02
COFLAM(133) = 1.53834605d+01
COFLAM(134) = -5.95777871d+00
COFLAM(135) = 1.17737648d+00
COFLAM(136) = -6.30156178d-02
COFLAM(137) = 7.42134531d+00
COFLAM(138) = -1.55621090d+00
COFLAM(139) = 4.24683399d-01
COFLAM(140) = -2.37628990d-02
COFLAM(141) = 9.56600558d+00
COFLAM(142) = -2.13078015d+00
COFLAM(143) = 4.54319631d-01
COFLAM(144) = -2.33435085d-02
COFLAM(145) = -1.39111538d+01
COFLAM(146) = 7.26782733d+00
COFLAM(147) = -7.97501475d-01
COFLAM(148) = 3.25210655d-02
COFLAM(149) = -1.14721212d+01
COFLAM(150) = 6.04006282d+00
COFLAM(151) = -6.06920866d-01
COFLAM(152) = 2.28508533d-02
COFLAM(153) = 1.57461008d+01
COFLAM(154) = -5.14565359d+00
COFLAM(155) = 9.25368949d-01
COFLAM(156) = -4.60945622d-02
COFLAM(157) = 4.51858742d-01
COFLAM(158) = 2.19081663d+00
COFLAM(159) = -2.29436447d-01
COFLAM(160) = 1.25479265d-02
COFLAM(161) = 1.36343230d+00
COFLAM(162) = 1.72616614d-03
COFLAM(163) = 3.04578843d-01
COFLAM(164) = -2.14829136d-02
COFLAM(165) = -2.42156758d+00
COFLAM(166) = 1.16422215d+00
COFLAM(167) = 2.07536654d-01
COFLAM(168) = -1.96618407d-02
COFLAM(169) = -7.34823323d+00
COFLAM(170) = 4.96238900d+00
COFLAM(171) = -4.72909115d-01
COFLAM(172) = 1.72133122d-02
COFLAM(173) = -1.18051819d+01
COFLAM(174) = 6.05433464d+00
COFLAM(175) = -5.78395212d-01
COFLAM(176) = 2.02076911d-02
COFLAM(177) = -3.19925179d+00
COFLAM(178) = 2.43089746d+00
COFLAM(179) = -7.68461937d-02
COFLAM(180) = -2.78657323d-03
COFLAM(181) = -6.28657249d+00
COFLAM(182) = 3.68588773d+00
COFLAM(183) = -2.44665048d-01
COFLAM(184) = 4.69897992d-03
COFLAM(185) = -1.12094260d+01
COFLAM(186) = 5.94941246d+00
COFLAM(187) = -5.91419413d-01
COFLAM(188) = 2.18094767d-02
COFLAM(189) = 1.29306158d+01
COFLAM(190) = -3.52817362d+00
COFLAM(191) = 6.45499013d-01
COFLAM(192) = -3.19375299d-02
COFLAM(193) = -3.17202048d+00
COFLAM(194) = 3.47025711d+00
COFLAM(195) = -3.75102111d-01
COFLAM(196) = 1.66086076d-02
COFLAM(197) = -1.92765552d+01
COFLAM(198) = 8.39477615d+00
COFLAM(199) = -8.05413956d-01
COFLAM(200) = 2.74022595d-02
COFLAM(201) = -1.93799923d+01
COFLAM(202) = 8.30033923d+00
COFLAM(203) = -7.71586937d-01
COFLAM(204) = 2.50297465d-02
COFLAM(205) = -6.27425421d+00
COFLAM(206) = 2.90471256d+00
COFLAM(207) = -4.35144205d-02
COFLAM(208) = -7.77922561d-03
COFLAM(209) = -9.88397056d+00
COFLAM(210) = 4.13891036d+00
COFLAM(211) = -1.74918892d-01
COFLAM(212) = -3.28392643d-03

end subroutine

! Poly fits for the diffusion coefficients, dim NO*KK*KK
subroutine egtransetCOFD(COFD)

implicit none

double precision, intent(out) :: COFD(11236)

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
COFD(29) = -1.19006548d+01
COFD(30) = 2.57507000d+00
COFD(31) = -1.24033737d-01
COFD(32) = 5.56694959d-03
COFD(33) = -1.08369481d+01
COFD(34) = 2.19094415d+00
COFD(35) = -7.11992510d-02
COFD(36) = 3.14105973d-03
COFD(37) = -1.09469245d+01
COFD(38) = 2.30836460d+00
COFD(39) = -8.76339315d-02
COFD(40) = 3.90878445d-03
COFD(41) = -1.25098960d+01
COFD(42) = 2.77873601d+00
COFD(43) = -1.50637360d-01
COFD(44) = 6.72684281d-03
COFD(45) = -1.25098960d+01
COFD(46) = 2.77873601d+00
COFD(47) = -1.50637360d-01
COFD(48) = 6.72684281d-03
COFD(49) = -1.25141260d+01
COFD(50) = 2.77873601d+00
COFD(51) = -1.50637360d-01
COFD(52) = 6.72684281d-03
COFD(53) = -1.24693568d+01
COFD(54) = 2.76686648d+00
COFD(55) = -1.49120141d-01
COFD(56) = 6.66220432d-03
COFD(57) = -1.17159737d+01
COFD(58) = 2.48123210d+00
COFD(59) = -1.11322604d-01
COFD(60) = 4.99282389d-03
COFD(61) = -1.37794315d+01
COFD(62) = 3.23973858d+00
COFD(63) = -2.09989036d-01
COFD(64) = 9.27667906d-03
COFD(65) = -1.60517370d+01
COFD(66) = 4.11188603d+00
COFD(67) = -3.21540884d-01
COFD(68) = 1.40482564d-02
COFD(69) = -1.60528285d+01
COFD(70) = 4.11188603d+00
COFD(71) = -3.21540884d-01
COFD(72) = 1.40482564d-02
COFD(73) = -1.58456300d+01
COFD(74) = 4.02074783d+00
COFD(75) = -3.10018522d-01
COFD(76) = 1.35599552d-02
COFD(77) = -1.58456300d+01
COFD(78) = 4.02074783d+00
COFD(79) = -3.10018522d-01
COFD(80) = 1.35599552d-02
COFD(81) = -1.59537247d+01
COFD(82) = 4.07051484d+00
COFD(83) = -3.16303109d-01
COFD(84) = 1.38259377d-02
COFD(85) = -1.34695359d+01
COFD(86) = 3.09379603d+00
COFD(87) = -1.91268635d-01
COFD(88) = 8.47480224d-03
COFD(89) = -1.34709807d+01
COFD(90) = 3.09379603d+00
COFD(91) = -1.91268635d-01
COFD(92) = 8.47480224d-03
COFD(93) = -1.34723215d+01
COFD(94) = 3.09379603d+00
COFD(95) = -1.91268635d-01
COFD(96) = 8.47480224d-03
COFD(97) = -1.42229194d+01
COFD(98) = 3.38669384d+00
COFD(99) = -2.28784122d-01
COFD(100) = 1.00790953d-02
COFD(101) = -1.39913897d+01
COFD(102) = 3.26384506d+00
COFD(103) = -2.12947087d-01
COFD(104) = 9.39743888d-03
COFD(105) = -1.39924781d+01
COFD(106) = 3.26384506d+00
COFD(107) = -2.12947087d-01
COFD(108) = 9.39743888d-03
COFD(109) = -1.22004324d+01
COFD(110) = 2.80725489d+00
COFD(111) = -1.54291406d-01
COFD(112) = 6.88290911d-03
COFD(113) = -1.57034851d+01
COFD(114) = 3.93614244d+00
COFD(115) = -2.99111497d-01
COFD(116) = 1.30888229d-02
COFD(117) = -1.57034851d+01
COFD(118) = 3.93614244d+00
COFD(119) = -2.99111497d-01
COFD(120) = 1.30888229d-02
COFD(121) = -1.08472922d+01
COFD(122) = 2.19094415d+00
COFD(123) = -7.11992510d-02
COFD(124) = 3.14105973d-03
COFD(125) = -1.09203269d+01
COFD(126) = 2.30836460d+00
COFD(127) = -8.76339315d-02
COFD(128) = 3.90878445d-03
COFD(129) = -1.09240642d+01
COFD(130) = 2.30836460d+00
COFD(131) = -8.76339315d-02
COFD(132) = 3.90878445d-03
COFD(133) = -1.61691241d+01
COFD(134) = 4.23579857d+00
COFD(135) = -3.36814551d-01
COFD(136) = 1.46782923d-02
COFD(137) = -1.10356312d+01
COFD(138) = 2.19094415d+00
COFD(139) = -7.11992510d-02
COFD(140) = 3.14105973d-03
COFD(141) = -1.16928639d+01
COFD(142) = 2.47469981d+00
COFD(143) = -1.10436257d-01
COFD(144) = 4.95273813d-03
COFD(145) = -1.31911278d+01
COFD(146) = 3.04970299d+00
COFD(147) = -1.85555523d-01
COFD(148) = 8.22773480d-03
COFD(149) = -1.36894329d+01
COFD(150) = 3.19999740d+00
COFD(151) = -2.04999020d-01
COFD(152) = 9.06766105d-03
COFD(153) = -1.20867697d+01
COFD(154) = 2.64389960d+00
COFD(155) = -1.33241706d-01
COFD(156) = 5.97810200d-03
COFD(157) = -1.11796892d+01
COFD(158) = 2.24409765d+00
COFD(159) = -7.86374012d-02
COFD(160) = 3.48840768d-03
COFD(161) = -1.64793369d+01
COFD(162) = 4.26176920d+00
COFD(163) = -3.39986562d-01
COFD(164) = 1.48077740d-02
COFD(165) = -1.64805864d+01
COFD(166) = 4.26176920d+00
COFD(167) = -3.39986562d-01
COFD(168) = 1.48077740d-02
COFD(169) = -1.22004340d+01
COFD(170) = 2.80725489d+00
COFD(171) = -1.54291406d-01
COFD(172) = 6.88290911d-03
COFD(173) = -1.36889305d+01
COFD(174) = 3.19999740d+00
COFD(175) = -2.04999020d-01
COFD(176) = 9.06766105d-03
COFD(177) = -1.36889305d+01
COFD(178) = 3.19999740d+00
COFD(179) = -2.04999020d-01
COFD(180) = 9.06766105d-03
COFD(181) = -1.36889305d+01
COFD(182) = 3.19999740d+00
COFD(183) = -2.04999020d-01
COFD(184) = 9.06766105d-03
COFD(185) = -1.36883939d+01
COFD(186) = 3.19999740d+00
COFD(187) = -2.04999020d-01
COFD(188) = 9.06766105d-03
COFD(189) = -1.16906297d+01
COFD(190) = 2.47469981d+00
COFD(191) = -1.10436257d-01
COFD(192) = 4.95273813d-03
COFD(193) = -1.23130152d+01
COFD(194) = 2.74418790d+00
COFD(195) = -1.46230156d-01
COFD(196) = 6.53948886d-03
COFD(197) = -1.43140669d+01
COFD(198) = 3.31177824d+00
COFD(199) = -2.18945280d-01
COFD(200) = 9.64764419d-03
COFD(201) = -1.43145779d+01
COFD(202) = 3.31177824d+00
COFD(203) = -2.18945280d-01
COFD(204) = 9.64764419d-03
COFD(205) = -1.57040212d+01
COFD(206) = 3.93614244d+00
COFD(207) = -2.99111497d-01
COFD(208) = 1.30888229d-02
COFD(209) = -1.57045332d+01
COFD(210) = 3.93614244d+00
COFD(211) = -2.99111497d-01
COFD(212) = 1.30888229d-02
COFD(213) = -1.14366381d+01
COFD(214) = 2.78323501d+00
COFD(215) = -1.51214064d-01
COFD(216) = 6.75150012d-03
COFD(217) = -1.47968712d+01
COFD(218) = 4.23027636d+00
COFD(219) = -3.36139991d-01
COFD(220) = 1.46507621d-02
COFD(221) = -1.34230272d+01
COFD(222) = 3.48624238d+00
COFD(223) = -2.41554467d-01
COFD(224) = 1.06263545d-02
COFD(225) = -1.46550083d+01
COFD(226) = 3.83606243d+00
COFD(227) = -2.86076532d-01
COFD(228) = 1.25205829d-02
COFD(229) = -1.34247866d+01
COFD(230) = 3.48624238d+00
COFD(231) = -2.41554467d-01
COFD(232) = 1.06263545d-02
COFD(233) = -1.95739570d+01
COFD(234) = 5.61113230d+00
COFD(235) = -4.90190187d-01
COFD(236) = 2.03260675d-02
COFD(237) = -1.46554748d+01
COFD(238) = 3.83606243d+00
COFD(239) = -2.86076532d-01
COFD(240) = 1.25205829d-02
COFD(241) = -1.46559141d+01
COFD(242) = 3.83606243d+00
COFD(243) = -2.86076532d-01
COFD(244) = 1.25205829d-02
COFD(245) = -1.32467319d+01
COFD(246) = 3.34156587d+00
COFD(247) = -2.22853306d-01
COFD(248) = 9.81883417d-03
COFD(249) = -1.34162893d+01
COFD(250) = 3.48624238d+00
COFD(251) = -2.41554467d-01
COFD(252) = 1.06263545d-02
COFD(253) = -1.57972369d+01
COFD(254) = 4.22225052d+00
COFD(255) = -3.35156428d-01
COFD(256) = 1.46104855d-02
COFD(257) = -1.57972369d+01
COFD(258) = 4.22225052d+00
COFD(259) = -3.35156428d-01
COFD(260) = 1.46104855d-02
COFD(261) = -1.57994893d+01
COFD(262) = 4.22225052d+00
COFD(263) = -3.35156428d-01
COFD(264) = 1.46104855d-02
COFD(265) = -1.57199037d+01
COFD(266) = 4.19936335d+00
COFD(267) = -3.32311009d-01
COFD(268) = 1.44921003d-02
COFD(269) = -1.43151174d+01
COFD(270) = 3.68038508d+00
COFD(271) = -2.65779346d-01
COFD(272) = 1.16360771d-02
COFD(273) = -1.76147026d+01
COFD(274) = 4.86049500d+00
COFD(275) = -4.12200578d-01
COFD(276) = 1.77160971d-02
COFD(277) = -1.97544450d+01
COFD(278) = 5.56931926d+00
COFD(279) = -4.89105511d-01
COFD(280) = 2.04493129d-02
COFD(281) = -1.97550088d+01
COFD(282) = 5.56931926d+00
COFD(283) = -4.89105511d-01
COFD(284) = 2.04493129d-02
COFD(285) = -1.92718582d+01
COFD(286) = 5.41172124d+00
COFD(287) = -4.73213887d-01
COFD(288) = 1.99405473d-02
COFD(289) = -1.92718582d+01
COFD(290) = 5.41172124d+00
COFD(291) = -4.73213887d-01
COFD(292) = 1.99405473d-02
COFD(293) = -1.96866103d+01
COFD(294) = 5.54637286d+00
COFD(295) = -4.87070324d-01
COFD(296) = 2.03983467d-02
COFD(297) = -1.72224724d+01
COFD(298) = 4.69060745d+00
COFD(299) = -3.92369888d-01
COFD(300) = 1.69459661d-02
COFD(301) = -1.72232223d+01
COFD(302) = 4.69060745d+00
COFD(303) = -3.92369888d-01
COFD(304) = 1.69459661d-02
COFD(305) = -1.72239172d+01
COFD(306) = 4.69060745d+00
COFD(307) = -3.92369888d-01
COFD(308) = 1.69459661d-02
COFD(309) = -1.82251914d+01
COFD(310) = 5.05237312d+00
COFD(311) = -4.35182396d-01
COFD(312) = 1.86363074d-02
COFD(313) = -1.79339327d+01
COFD(314) = 4.91373893d+00
COFD(315) = -4.18747629d-01
COFD(316) = 1.79856610d-02
COFD(317) = -1.79344949d+01
COFD(318) = 4.91373893d+00
COFD(319) = -4.18747629d-01
COFD(320) = 1.79856610d-02
COFD(321) = -1.54460820d+01
COFD(322) = 4.26819983d+00
COFD(323) = -3.40766379d-01
COFD(324) = 1.48393361d-02
COFD(325) = -1.94688688d+01
COFD(326) = 5.43830787d+00
COFD(327) = -4.75472880d-01
COFD(328) = 1.99909996d-02
COFD(329) = -1.94688688d+01
COFD(330) = 5.43830787d+00
COFD(331) = -4.75472880d-01
COFD(332) = 1.99909996d-02
COFD(333) = -1.32522778d+01
COFD(334) = 3.34156587d+00
COFD(335) = -2.22853306d-01
COFD(336) = 9.81883417d-03
COFD(337) = -1.33789807d+01
COFD(338) = 3.48624238d+00
COFD(339) = -2.41554467d-01
COFD(340) = 1.06263545d-02
COFD(341) = -1.33809634d+01
COFD(342) = 3.48624238d+00
COFD(343) = -2.41554467d-01
COFD(344) = 1.06263545d-02
COFD(345) = -1.94027559d+01
COFD(346) = 5.54487588d+00
COFD(347) = -4.86919661d-01
COFD(348) = 2.03935297d-02
COFD(349) = -1.34487066d+01
COFD(350) = 3.34156587d+00
COFD(351) = -2.22853306d-01
COFD(352) = 9.81883417d-03
COFD(353) = -1.42905987d+01
COFD(354) = 3.67490723d+00
COFD(355) = -2.65114792d-01
COFD(356) = 1.16092671d-02
COFD(357) = -1.68734643d+01
COFD(358) = 4.63687143d+00
COFD(359) = -3.85900861d-01
COFD(360) = 1.66856798d-02
COFD(361) = -1.74819100d+01
COFD(362) = 4.80792005d+00
COFD(363) = -4.06126584d-01
COFD(364) = 1.74831083d-02
COFD(365) = -1.50140778d+01
COFD(366) = 3.96853403d+00
COFD(367) = -3.03320126d-01
COFD(368) = 1.32719819d-02
COFD(369) = -1.36558633d+01
COFD(370) = 3.41450296d+00
COFD(371) = -2.32431042d-01
COFD(372) = 1.02388399d-02
COFD(373) = -1.99775614d+01
COFD(374) = 5.61184673d+00
COFD(375) = -4.90499758d-01
COFD(376) = 2.03480832d-02
COFD(377) = -1.99782082d+01
COFD(378) = 5.61184673d+00
COFD(379) = -4.90499758d-01
COFD(380) = 2.03480832d-02
COFD(381) = -1.54460828d+01
COFD(382) = 4.26819983d+00
COFD(383) = -3.40766379d-01
COFD(384) = 1.48393361d-02
COFD(385) = -1.74816531d+01
COFD(386) = 4.80792005d+00
COFD(387) = -4.06126584d-01
COFD(388) = 1.74831083d-02
COFD(389) = -1.74816531d+01
COFD(390) = 4.80792005d+00
COFD(391) = -4.06126584d-01
COFD(392) = 1.74831083d-02
COFD(393) = -1.74816531d+01
COFD(394) = 4.80792005d+00
COFD(395) = -4.06126584d-01
COFD(396) = 1.74831083d-02
COFD(397) = -1.74813786d+01
COFD(398) = 4.80792005d+00
COFD(399) = -4.06126584d-01
COFD(400) = 1.74831083d-02
COFD(401) = -1.42894441d+01
COFD(402) = 3.67490723d+00
COFD(403) = -2.65114792d-01
COFD(404) = 1.16092671d-02
COFD(405) = -1.54738604d+01
COFD(406) = 4.15765300d+00
COFD(407) = -3.27126237d-01
COFD(408) = 1.42762611d-02
COFD(409) = -1.83545293d+01
COFD(410) = 4.98756925d+00
COFD(411) = -4.27526072d-01
COFD(412) = 1.83341755d-02
COFD(413) = -1.83547906d+01
COFD(414) = 4.98756925d+00
COFD(415) = -4.27526072d-01
COFD(416) = 1.83341755d-02
COFD(417) = -1.94691430d+01
COFD(418) = 5.43830787d+00
COFD(419) = -4.75472880d-01
COFD(420) = 1.99909996d-02
COFD(421) = -1.94694048d+01
COFD(422) = 5.43830787d+00
COFD(423) = -4.75472880d-01
COFD(424) = 1.99909996d-02
COFD(425) = -1.09595712d+01
COFD(426) = 2.30836460d+00
COFD(427) = -8.76339315d-02
COFD(428) = 3.90878445d-03
COFD(429) = -1.34230272d+01
COFD(430) = 3.48624238d+00
COFD(431) = -2.41554467d-01
COFD(432) = 1.06263545d-02
COFD(433) = -1.32093628d+01
COFD(434) = 2.90778936d+00
COFD(435) = -1.67388544d-01
COFD(436) = 7.45220609d-03
COFD(437) = -1.43139231d+01
COFD(438) = 3.17651319d+00
COFD(439) = -2.02028974d-01
COFD(440) = 8.94232502d-03
COFD(441) = -1.32244035d+01
COFD(442) = 2.90778936d+00
COFD(443) = -1.67388544d-01
COFD(444) = 7.45220609d-03
COFD(445) = -1.94093572d+01
COFD(446) = 5.16013126d+00
COFD(447) = -4.46824543d-01
COFD(448) = 1.90464887d-02
COFD(449) = -1.43190389d+01
COFD(450) = 3.17651319d+00
COFD(451) = -2.02028974d-01
COFD(452) = 8.94232502d-03
COFD(453) = -1.43238998d+01
COFD(454) = 3.17651319d+00
COFD(455) = -2.02028974d-01
COFD(456) = 8.94232502d-03
COFD(457) = -1.30610083d+01
COFD(458) = 2.80913567d+00
COFD(459) = -1.54536855d-01
COFD(460) = 6.89359313d-03
COFD(461) = -1.31551788d+01
COFD(462) = 2.90778936d+00
COFD(463) = -1.67388544d-01
COFD(464) = 7.45220609d-03
COFD(465) = -1.50584249d+01
COFD(466) = 3.47945612d+00
COFD(467) = -2.40703722d-01
COFD(468) = 1.05907441d-02
COFD(469) = -1.50584249d+01
COFD(470) = 3.47945612d+00
COFD(471) = -2.40703722d-01
COFD(472) = 1.05907441d-02
COFD(473) = -1.50766130d+01
COFD(474) = 3.47945612d+00
COFD(475) = -2.40703722d-01
COFD(476) = 1.05907441d-02
COFD(477) = -1.50270339d+01
COFD(478) = 3.46140064d+00
COFD(479) = -2.38440092d-01
COFD(480) = 1.04960087d-02
COFD(481) = -1.40999008d+01
COFD(482) = 3.08120012d+00
COFD(483) = -1.89629903d-01
COFD(484) = 8.40361952d-03
COFD(485) = -1.70534856d+01
COFD(486) = 4.14240922d+00
COFD(487) = -3.25239774d-01
COFD(488) = 1.41980687d-02
COFD(489) = -1.94313116d+01
COFD(490) = 5.02567894d+00
COFD(491) = -4.32045169d-01
COFD(492) = 1.85132214d-02
COFD(493) = -1.94373127d+01
COFD(494) = 5.02567894d+00
COFD(495) = -4.32045169d-01
COFD(496) = 1.85132214d-02
COFD(497) = -1.88179418d+01
COFD(498) = 4.79683898d+00
COFD(499) = -4.04829719d-01
COFD(500) = 1.74325475d-02
COFD(501) = -1.88179418d+01
COFD(502) = 4.79683898d+00
COFD(503) = -4.04829719d-01
COFD(504) = 1.74325475d-02
COFD(505) = -1.93364585d+01
COFD(506) = 4.98286777d+00
COFD(507) = -4.26970814d-01
COFD(508) = 1.83122917d-02
COFD(509) = -1.65412306d+01
COFD(510) = 3.95035840d+00
COFD(511) = -3.00959418d-01
COFD(512) = 1.31692593d-02
COFD(513) = -1.65488358d+01
COFD(514) = 3.95035840d+00
COFD(515) = -3.00959418d-01
COFD(516) = 1.31692593d-02
COFD(517) = -1.65559787d+01
COFD(518) = 3.95035840d+00
COFD(519) = -3.00959418d-01
COFD(520) = 1.31692593d-02
COFD(521) = -1.74792112d+01
COFD(522) = 4.29676909d+00
COFD(523) = -3.44085306d-01
COFD(524) = 1.49671135d-02
COFD(525) = -1.72496634d+01
COFD(526) = 4.17889917d+00
COFD(527) = -3.29752510d-01
COFD(528) = 1.43850275d-02
COFD(529) = -1.72556499d+01
COFD(530) = 4.17889917d+00
COFD(531) = -3.29752510d-01
COFD(532) = 1.43850275d-02
COFD(533) = -1.49500357d+01
COFD(534) = 3.52327209d+00
COFD(535) = -2.46286208d-01
COFD(536) = 1.08285963d-02
COFD(537) = -1.90883268d+01
COFD(538) = 4.84384483d+00
COFD(539) = -4.10265575d-01
COFD(540) = 1.76414287d-02
COFD(541) = -1.90883268d+01
COFD(542) = 4.84384483d+00
COFD(543) = -4.10265575d-01
COFD(544) = 1.76414287d-02
COFD(545) = -1.31034488d+01
COFD(546) = 2.80913567d+00
COFD(547) = -1.54536855d-01
COFD(548) = 6.89359313d-03
COFD(549) = -1.31565315d+01
COFD(550) = 2.90778936d+00
COFD(551) = -1.67388544d-01
COFD(552) = 7.45220609d-03
COFD(553) = -1.31730273d+01
COFD(554) = 2.90778936d+00
COFD(555) = -1.67388544d-01
COFD(556) = 7.45220609d-03
COFD(557) = -1.89669545d+01
COFD(558) = 4.98075560d+00
COFD(559) = -4.26721620d-01
COFD(560) = 1.83024823d-02
COFD(561) = -1.34236995d+01
COFD(562) = 2.80913567d+00
COFD(563) = -1.54536855d-01
COFD(564) = 6.89359313d-03
COFD(565) = -1.40879121d+01
COFD(566) = 3.07549274d+00
COFD(567) = -1.88889344d-01
COFD(568) = 8.37152866d-03
COFD(569) = -1.62772222d+01
COFD(570) = 3.88250968d+00
COFD(571) = -2.92155848d-01
COFD(572) = 1.27867850d-02
COFD(573) = -1.69290095d+01
COFD(574) = 4.09077642d+00
COFD(575) = -3.18894990d-01
COFD(576) = 1.39371445d-02
COFD(577) = -1.45065606d+01
COFD(578) = 3.24450689d+00
COFD(579) = -2.10570734d-01
COFD(580) = 9.30026771d-03
COFD(581) = -1.35403054d+01
COFD(582) = 2.85448298d+00
COFD(583) = -1.60491387d-01
COFD(584) = 7.15451323d-03
COFD(585) = -1.98343743d+01
COFD(586) = 5.15754807d+00
COFD(587) = -4.46654084d-01
COFD(588) = 1.90458977d-02
COFD(589) = -1.98411047d+01
COFD(590) = 5.15754807d+00
COFD(591) = -4.46654084d-01
COFD(592) = 1.90458977d-02
COFD(593) = -1.49500454d+01
COFD(594) = 3.52327209d+00
COFD(595) = -2.46286208d-01
COFD(596) = 1.08285963d-02
COFD(597) = -1.69259591d+01
COFD(598) = 4.09077642d+00
COFD(599) = -3.18894990d-01
COFD(600) = 1.39371445d-02
COFD(601) = -1.69259591d+01
COFD(602) = 4.09077642d+00
COFD(603) = -3.18894990d-01
COFD(604) = 1.39371445d-02
COFD(605) = -1.69259591d+01
COFD(606) = 4.09077642d+00
COFD(607) = -3.18894990d-01
COFD(608) = 1.39371445d-02
COFD(609) = -1.69227183d+01
COFD(610) = 4.09077642d+00
COFD(611) = -3.18894990d-01
COFD(612) = 1.39371445d-02
COFD(613) = -1.40756935d+01
COFD(614) = 3.07549274d+00
COFD(615) = -1.88889344d-01
COFD(616) = 8.37152866d-03
COFD(617) = -1.49610527d+01
COFD(618) = 3.41988961d+00
COFD(619) = -2.33128386d-01
COFD(620) = 1.02689994d-02
COFD(621) = -1.76841045d+01
COFD(622) = 4.24719726d+00
COFD(623) = -3.38206061d-01
COFD(624) = 1.47350654d-02
COFD(625) = -1.76872087d+01
COFD(626) = 4.24719726d+00
COFD(627) = -3.38206061d-01
COFD(628) = 1.47350654d-02
COFD(629) = -1.90915649d+01
COFD(630) = 4.84384483d+00
COFD(631) = -4.10265575d-01
COFD(632) = 1.76414287d-02
COFD(633) = -1.90946745d+01
COFD(634) = 4.84384483d+00
COFD(635) = -4.10265575d-01
COFD(636) = 1.76414287d-02
COFD(637) = -1.18988955d+01
COFD(638) = 2.57507000d+00
COFD(639) = -1.24033737d-01
COFD(640) = 5.56694959d-03
COFD(641) = -1.46550083d+01
COFD(642) = 3.83606243d+00
COFD(643) = -2.86076532d-01
COFD(644) = 1.25205829d-02
COFD(645) = -1.43139231d+01
COFD(646) = 3.17651319d+00
COFD(647) = -2.02028974d-01
COFD(648) = 8.94232502d-03
COFD(649) = -1.55511344d+01
COFD(650) = 3.48070094d+00
COFD(651) = -2.40859499d-01
COFD(652) = 1.05972514d-02
COFD(653) = -1.43340796d+01
COFD(654) = 3.17651319d+00
COFD(655) = -2.02028974d-01
COFD(656) = 8.94232502d-03
COFD(657) = -2.12652533d+01
COFD(658) = 5.59961818d+00
COFD(659) = -4.91624858d-01
COFD(660) = 2.05035550d-02
COFD(661) = -1.55588279d+01
COFD(662) = 3.48070094d+00
COFD(663) = -2.40859499d-01
COFD(664) = 1.05972514d-02
COFD(665) = -1.55661750d+01
COFD(666) = 3.48070094d+00
COFD(667) = -2.40859499d-01
COFD(668) = 1.05972514d-02
COFD(669) = -1.40707473d+01
COFD(670) = 3.05837263d+00
COFD(671) = -1.86672802d-01
COFD(672) = 8.27575734d-03
COFD(673) = -1.42429085d+01
COFD(674) = 3.17651319d+00
COFD(675) = -2.02028974d-01
COFD(676) = 8.94232502d-03
COFD(677) = -1.63254691d+01
COFD(678) = 3.82388595d+00
COFD(679) = -2.84480724d-01
COFD(680) = 1.24506311d-02
COFD(681) = -1.63254691d+01
COFD(682) = 3.82388595d+00
COFD(683) = -2.84480724d-01
COFD(684) = 1.24506311d-02
COFD(685) = -1.63493345d+01
COFD(686) = 3.82388595d+00
COFD(687) = -2.84480724d-01
COFD(688) = 1.24506311d-02
COFD(689) = -1.62724462d+01
COFD(690) = 3.79163564d+00
COFD(691) = -2.80257365d-01
COFD(692) = 1.22656902d-02
COFD(693) = -1.52721107d+01
COFD(694) = 3.36790500d+00
COFD(695) = -2.26321740d-01
COFD(696) = 9.97135055d-03
COFD(697) = -1.84688406d+01
COFD(698) = 4.49330851d+00
COFD(699) = -3.68208715d-01
COFD(700) = 1.59565402d-02
COFD(701) = -2.08204449d+01
COFD(702) = 5.35267674d+00
COFD(703) = -4.69010505d-01
COFD(704) = 1.98979152d-02
COFD(705) = -2.08293255d+01
COFD(706) = 5.35267674d+00
COFD(707) = -4.69010505d-01
COFD(708) = 1.98979152d-02
COFD(709) = -2.04928958d+01
COFD(710) = 5.22397933d+00
COFD(711) = -4.54138171d-01
COFD(712) = 1.93249285d-02
COFD(713) = -2.04928958d+01
COFD(714) = 5.22397933d+00
COFD(715) = -4.54138171d-01
COFD(716) = 1.93249285d-02
COFD(717) = -2.07595845d+01
COFD(718) = 5.32244593d+00
COFD(719) = -4.65829403d-01
COFD(720) = 1.97895274d-02
COFD(721) = -1.78725135d+01
COFD(722) = 4.29613154d+00
COFD(723) = -3.44012526d-01
COFD(724) = 1.49643715d-02
COFD(725) = -1.78834935d+01
COFD(726) = 4.29613154d+00
COFD(727) = -3.44012526d-01
COFD(728) = 1.49643715d-02
COFD(729) = -1.78938745d+01
COFD(730) = 4.29613154d+00
COFD(731) = -3.44012526d-01
COFD(732) = 1.49643715d-02
COFD(733) = -1.89544778d+01
COFD(734) = 4.68595732d+00
COFD(735) = -3.91842840d-01
COFD(736) = 1.69262542d-02
COFD(737) = -1.86335932d+01
COFD(738) = 4.53572533d+00
COFD(739) = -3.73386925d-01
COFD(740) = 1.61678881d-02
COFD(741) = -1.86424545d+01
COFD(742) = 4.53572533d+00
COFD(743) = -3.73386925d-01
COFD(744) = 1.61678881d-02
COFD(745) = -1.64169433d+01
COFD(746) = 3.89309916d+00
COFD(747) = -2.93528188d-01
COFD(748) = 1.28463177d-02
COFD(749) = -2.05184870d+01
COFD(750) = 5.18417470d+00
COFD(751) = -4.49491573d-01
COFD(752) = 1.91438508d-02
COFD(753) = -2.05184870d+01
COFD(754) = 5.18417470d+00
COFD(755) = -4.49491573d-01
COFD(756) = 1.91438508d-02
COFD(757) = -1.41254249d+01
COFD(758) = 3.05837263d+00
COFD(759) = -1.86672802d-01
COFD(760) = 8.27575734d-03
COFD(761) = -1.42600473d+01
COFD(762) = 3.17651319d+00
COFD(763) = -2.02028974d-01
COFD(764) = 8.94232502d-03
COFD(765) = -1.42819281d+01
COFD(766) = 3.17651319d+00
COFD(767) = -2.02028974d-01
COFD(768) = 8.94232502d-03
COFD(769) = -2.06107018d+01
COFD(770) = 5.38762344d+00
COFD(771) = -4.71526404d-01
COFD(772) = 1.99251819d-02
COFD(773) = -1.44912469d+01
COFD(774) = 3.05837263d+00
COFD(775) = -1.86672802d-01
COFD(776) = 8.27575734d-03
COFD(777) = -1.52594746d+01
COFD(778) = 3.35922578d+00
COFD(779) = -2.25181399d-01
COFD(780) = 9.92132878d-03
COFD(781) = -1.77337342d+01
COFD(782) = 4.25438185d+00
COFD(783) = -3.39084808d-01
COFD(784) = 1.47709916d-02
COFD(785) = -1.83308393d+01
COFD(786) = 4.43878381d+00
COFD(787) = -3.61697624d-01
COFD(788) = 1.56975581d-02
COFD(789) = -1.57967324d+01
COFD(790) = 3.57122377d+00
COFD(791) = -2.52409987d-01
COFD(792) = 1.10900562d-02
COFD(793) = -1.46137599d+01
COFD(794) = 3.10987301d+00
COFD(795) = -1.93373245d-01
COFD(796) = 8.56678670d-03
COFD(797) = -2.10651573d+01
COFD(798) = 5.41544337d+00
COFD(799) = -4.73380953d-01
COFD(800) = 1.99350119d-02
COFD(801) = -2.10749998d+01
COFD(802) = 5.41544337d+00
COFD(803) = -4.73380953d-01
COFD(804) = 1.99350119d-02
COFD(805) = -1.64169585d+01
COFD(806) = 3.89309916d+00
COFD(807) = -2.93528188d-01
COFD(808) = 1.28463177d-02
COFD(809) = -1.83260311d+01
COFD(810) = 4.43878381d+00
COFD(811) = -3.61697624d-01
COFD(812) = 1.56975581d-02
COFD(813) = -1.83260311d+01
COFD(814) = 4.43878381d+00
COFD(815) = -3.61697624d-01
COFD(816) = 1.56975581d-02
COFD(817) = -1.83260311d+01
COFD(818) = 4.43878381d+00
COFD(819) = -3.61697624d-01
COFD(820) = 1.56975581d-02
COFD(821) = -1.83209411d+01
COFD(822) = 4.43878381d+00
COFD(823) = -3.61697624d-01
COFD(824) = 1.56975581d-02
COFD(825) = -1.52414485d+01
COFD(826) = 3.35922578d+00
COFD(827) = -2.25181399d-01
COFD(828) = 9.92132878d-03
COFD(829) = -1.62380075d+01
COFD(830) = 3.72612300d+00
COFD(831) = -2.71663673d-01
COFD(832) = 1.18889643d-02
COFD(833) = -1.91225244d+01
COFD(834) = 4.61801405d+00
COFD(835) = -3.83535652d-01
COFD(836) = 1.65862513d-02
COFD(837) = -1.91274187d+01
COFD(838) = 4.61801405d+00
COFD(839) = -3.83535652d-01
COFD(840) = 1.65862513d-02
COFD(841) = -2.05235731d+01
COFD(842) = 5.18417470d+00
COFD(843) = -4.49491573d-01
COFD(844) = 1.91438508d-02
COFD(845) = -2.05284752d+01
COFD(846) = 5.18417470d+00
COFD(847) = -4.49491573d-01
COFD(848) = 1.91438508d-02
COFD(849) = -1.09628982d+01
COFD(850) = 2.30836460d+00
COFD(851) = -8.76339315d-02
COFD(852) = 3.90878445d-03
COFD(853) = -1.34247866d+01
COFD(854) = 3.48624238d+00
COFD(855) = -2.41554467d-01
COFD(856) = 1.06263545d-02
COFD(857) = -1.32244035d+01
COFD(858) = 2.90778936d+00
COFD(859) = -1.67388544d-01
COFD(860) = 7.45220609d-03
COFD(861) = -1.43340796d+01
COFD(862) = 3.17651319d+00
COFD(863) = -2.02028974d-01
COFD(864) = 8.94232502d-03
COFD(865) = -1.32399106d+01
COFD(866) = 2.90778936d+00
COFD(867) = -1.67388544d-01
COFD(868) = 7.45220609d-03
COFD(869) = -1.94253036d+01
COFD(870) = 5.16013126d+00
COFD(871) = -4.46824543d-01
COFD(872) = 1.90464887d-02
COFD(873) = -1.43394069d+01
COFD(874) = 3.17651319d+00
COFD(875) = -2.02028974d-01
COFD(876) = 8.94232502d-03
COFD(877) = -1.43444709d+01
COFD(878) = 3.17651319d+00
COFD(879) = -2.02028974d-01
COFD(880) = 8.94232502d-03
COFD(881) = -1.30738796d+01
COFD(882) = 2.80913567d+00
COFD(883) = -1.54536855d-01
COFD(884) = 6.89359313d-03
COFD(885) = -1.31686537d+01
COFD(886) = 2.90778936d+00
COFD(887) = -1.67388544d-01
COFD(888) = 7.45220609d-03
COFD(889) = -1.50724636d+01
COFD(890) = 3.47945612d+00
COFD(891) = -2.40703722d-01
COFD(892) = 1.05907441d-02
COFD(893) = -1.50724636d+01
COFD(894) = 3.47945612d+00
COFD(895) = -2.40703722d-01
COFD(896) = 1.05907441d-02
COFD(897) = -1.50911794d+01
COFD(898) = 3.47945612d+00
COFD(899) = -2.40703722d-01
COFD(900) = 1.05907441d-02
COFD(901) = -1.50420953d+01
COFD(902) = 3.46140064d+00
COFD(903) = -2.38440092d-01
COFD(904) = 1.04960087d-02
COFD(905) = -1.41191261d+01
COFD(906) = 3.08120012d+00
COFD(907) = -1.89629903d-01
COFD(908) = 8.40361952d-03
COFD(909) = -1.70757047d+01
COFD(910) = 4.14240922d+00
COFD(911) = -3.25239774d-01
COFD(912) = 1.41980687d-02
COFD(913) = -1.94507876d+01
COFD(914) = 5.02567894d+00
COFD(915) = -4.32045169d-01
COFD(916) = 1.85132214d-02
COFD(917) = -1.94570287d+01
COFD(918) = 5.02567894d+00
COFD(919) = -4.32045169d-01
COFD(920) = 1.85132214d-02
COFD(921) = -1.88378874d+01
COFD(922) = 4.79683898d+00
COFD(923) = -4.04829719d-01
COFD(924) = 1.74325475d-02
COFD(925) = -1.88378874d+01
COFD(926) = 4.79683898d+00
COFD(927) = -4.04829719d-01
COFD(928) = 1.74325475d-02
COFD(929) = -1.93566243d+01
COFD(930) = 4.98286777d+00
COFD(931) = -4.26970814d-01
COFD(932) = 1.83122917d-02
COFD(933) = -1.65596434d+01
COFD(934) = 3.95035840d+00
COFD(935) = -3.00959418d-01
COFD(936) = 1.31692593d-02
COFD(937) = -1.65675362d+01
COFD(938) = 3.95035840d+00
COFD(939) = -3.00959418d-01
COFD(940) = 1.31692593d-02
COFD(941) = -1.65749533d+01
COFD(942) = 3.95035840d+00
COFD(943) = -3.00959418d-01
COFD(944) = 1.31692593d-02
COFD(945) = -1.74984476d+01
COFD(946) = 4.29676909d+00
COFD(947) = -3.44085306d-01
COFD(948) = 1.49671135d-02
COFD(949) = -1.72691500d+01
COFD(950) = 4.17889917d+00
COFD(951) = -3.29752510d-01
COFD(952) = 1.43850275d-02
COFD(953) = -1.72753760d+01
COFD(954) = 4.17889917d+00
COFD(955) = -3.29752510d-01
COFD(956) = 1.43850275d-02
COFD(957) = -1.49718233d+01
COFD(958) = 3.52327209d+00
COFD(959) = -2.46286208d-01
COFD(960) = 1.08285963d-02
COFD(961) = -1.91102652d+01
COFD(962) = 4.84384483d+00
COFD(963) = -4.10265575d-01
COFD(964) = 1.76414287d-02
COFD(965) = -1.91102652d+01
COFD(966) = 4.84384483d+00
COFD(967) = -4.10265575d-01
COFD(968) = 1.76414287d-02
COFD(969) = -1.31174764d+01
COFD(970) = 2.80913567d+00
COFD(971) = -1.54536855d-01
COFD(972) = 6.89359313d-03
COFD(973) = -1.31710876d+01
COFD(974) = 2.90778936d+00
COFD(975) = -1.67388544d-01
COFD(976) = 7.45220609d-03
COFD(977) = -1.31880790d+01
COFD(978) = 2.90778936d+00
COFD(979) = -1.67388544d-01
COFD(980) = 7.45220609d-03
COFD(981) = -1.89824721d+01
COFD(982) = 4.98075560d+00
COFD(983) = -4.26721620d-01
COFD(984) = 1.83024823d-02
COFD(985) = -1.34431763d+01
COFD(986) = 2.80913567d+00
COFD(987) = -1.54536855d-01
COFD(988) = 6.89359313d-03
COFD(989) = -1.41076233d+01
COFD(990) = 3.07549274d+00
COFD(991) = -1.88889344d-01
COFD(992) = 8.37152866d-03
COFD(993) = -1.62997072d+01
COFD(994) = 3.88250968d+00
COFD(995) = -2.92155848d-01
COFD(996) = 1.27867850d-02
COFD(997) = -1.69512290d+01
COFD(998) = 4.09077642d+00
COFD(999) = -3.18894990d-01
COFD(1000) = 1.39371445d-02
COFD(1001) = -1.45265017d+01
COFD(1002) = 3.24450689d+00
COFD(1003) = -2.10570734d-01
COFD(1004) = 9.30026771d-03
COFD(1005) = -1.35590001d+01
COFD(1006) = 2.85448298d+00
COFD(1007) = -1.60491387d-01
COFD(1008) = 7.15451323d-03
COFD(1009) = -1.98533435d+01
COFD(1010) = 5.15754807d+00
COFD(1011) = -4.46654084d-01
COFD(1012) = 1.90458977d-02
COFD(1013) = -1.98603359d+01
COFD(1014) = 5.15754807d+00
COFD(1015) = -4.46654084d-01
COFD(1016) = 1.90458977d-02
COFD(1017) = -1.49718335d+01
COFD(1018) = 3.52327209d+00
COFD(1019) = -2.46286208d-01
COFD(1020) = 1.08285963d-02
COFD(1021) = -1.69480404d+01
COFD(1022) = 4.09077642d+00
COFD(1023) = -3.18894990d-01
COFD(1024) = 1.39371445d-02
COFD(1025) = -1.69480404d+01
COFD(1026) = 4.09077642d+00
COFD(1027) = -3.18894990d-01
COFD(1028) = 1.39371445d-02
COFD(1029) = -1.69480404d+01
COFD(1030) = 4.09077642d+00
COFD(1031) = -3.18894990d-01
COFD(1032) = 1.39371445d-02
COFD(1033) = -1.69446537d+01
COFD(1034) = 4.09077642d+00
COFD(1035) = -3.18894990d-01
COFD(1036) = 1.39371445d-02
COFD(1037) = -1.40949196d+01
COFD(1038) = 3.07549274d+00
COFD(1039) = -1.88889344d-01
COFD(1040) = 8.37152866d-03
COFD(1041) = -1.49826725d+01
COFD(1042) = 3.41988961d+00
COFD(1043) = -2.33128386d-01
COFD(1044) = 1.02689994d-02
COFD(1045) = -1.77061949d+01
COFD(1046) = 4.24719726d+00
COFD(1047) = -3.38206061d-01
COFD(1048) = 1.47350654d-02
COFD(1049) = -1.77094398d+01
COFD(1050) = 4.24719726d+00
COFD(1051) = -3.38206061d-01
COFD(1052) = 1.47350654d-02
COFD(1053) = -1.91136491d+01
COFD(1054) = 4.84384483d+00
COFD(1055) = -4.10265575d-01
COFD(1056) = 1.76414287d-02
COFD(1057) = -1.91168996d+01
COFD(1058) = 4.84384483d+00
COFD(1059) = -4.10265575d-01
COFD(1060) = 1.76414287d-02
COFD(1061) = -1.71982995d+01
COFD(1062) = 4.63881404d+00
COFD(1063) = -3.86139633d-01
COFD(1064) = 1.66955081d-02
COFD(1065) = -1.95739570d+01
COFD(1066) = 5.61113230d+00
COFD(1067) = -4.90190187d-01
COFD(1068) = 2.03260675d-02
COFD(1069) = -1.94093572d+01
COFD(1070) = 5.16013126d+00
COFD(1071) = -4.46824543d-01
COFD(1072) = 1.90464887d-02
COFD(1073) = -2.12652533d+01
COFD(1074) = 5.59961818d+00
COFD(1075) = -4.91624858d-01
COFD(1076) = 2.05035550d-02
COFD(1077) = -1.94253036d+01
COFD(1078) = 5.16013126d+00
COFD(1079) = -4.46824543d-01
COFD(1080) = 1.90464887d-02
COFD(1081) = -1.19157919d+01
COFD(1082) = 9.28955130d-01
COFD(1083) = 2.42107090d-01
COFD(1084) = -1.59823963d-02
COFD(1085) = -2.06463744d+01
COFD(1086) = 5.41688482d+00
COFD(1087) = -4.73387188d-01
COFD(1088) = 1.99280175d-02
COFD(1089) = -2.06516336d+01
COFD(1090) = 5.41688482d+00
COFD(1091) = -4.73387188d-01
COFD(1092) = 1.99280175d-02
COFD(1093) = -1.92003817d+01
COFD(1094) = 5.05708283d+00
COFD(1095) = -4.35739290d-01
COFD(1096) = 1.86583205d-02
COFD(1097) = -1.93521390d+01
COFD(1098) = 5.16013126d+00
COFD(1099) = -4.46824543d-01
COFD(1100) = 1.90464887d-02
COFD(1101) = -2.12639214d+01
COFD(1102) = 5.61184117d+00
COFD(1103) = -4.90532156d-01
COFD(1104) = 2.03507922d-02
COFD(1105) = -2.12639214d+01
COFD(1106) = 5.61184117d+00
COFD(1107) = -4.90532156d-01
COFD(1108) = 2.03507922d-02
COFD(1109) = -2.12831323d+01
COFD(1110) = 5.61184117d+00
COFD(1111) = -4.90532156d-01
COFD(1112) = 2.03507922d-02
COFD(1113) = -2.14087397d+01
COFD(1114) = 5.57282008d+00
COFD(1115) = -4.76690890d-01
COFD(1116) = 1.94000719d-02
COFD(1117) = -2.11388331d+01
COFD(1118) = 5.55529675d+00
COFD(1119) = -4.87942518d-01
COFD(1120) = 2.04249054d-02
COFD(1121) = -2.07653719d+01
COFD(1122) = 5.01092022d+00
COFD(1123) = -3.77985635d-01
COFD(1124) = 1.40968645d-02
COFD(1125) = -1.77498543d+01
COFD(1126) = 3.57475686d+00
COFD(1127) = -1.56396297d-01
COFD(1128) = 3.12157721d-03
COFD(1129) = -1.77563250d+01
COFD(1130) = 3.57475686d+00
COFD(1131) = -1.56396297d-01
COFD(1132) = 3.12157721d-03
COFD(1133) = -1.65295288d+01
COFD(1134) = 2.97569206d+00
COFD(1135) = -6.75652842d-02
COFD(1136) = -1.08648422d-03
COFD(1137) = -1.65295288d+01
COFD(1138) = 2.97569206d+00
COFD(1139) = -6.75652842d-02
COFD(1140) = -1.08648422d-03
COFD(1141) = -1.80253664d+01
COFD(1142) = 3.69199168d+00
COFD(1143) = -1.74005516d-01
COFD(1144) = 3.97694372d-03
COFD(1145) = -2.15014310d+01
COFD(1146) = 5.46737673d+00
COFD(1147) = -4.55696085d-01
COFD(1148) = 1.81982625d-02
COFD(1149) = -2.15095980d+01
COFD(1150) = 5.46737673d+00
COFD(1151) = -4.55696085d-01
COFD(1152) = 1.81982625d-02
COFD(1153) = -2.15172770d+01
COFD(1154) = 5.46737673d+00
COFD(1155) = -4.55696085d-01
COFD(1156) = 1.81982625d-02
COFD(1157) = -2.08812333d+01
COFD(1158) = 5.08859217d+00
COFD(1159) = -3.90525428d-01
COFD(1160) = 1.47376395d-02
COFD(1161) = -2.12597312d+01
COFD(1162) = 5.24930667d+00
COFD(1163) = -4.17435088d-01
COFD(1164) = 1.61434424d-02
COFD(1165) = -2.12661865d+01
COFD(1166) = 5.24930667d+00
COFD(1167) = -4.17435088d-01
COFD(1168) = 1.61434424d-02
COFD(1169) = -2.10440675d+01
COFD(1170) = 5.59806282d+00
COFD(1171) = -4.87109535d-01
COFD(1172) = 2.01370226d-02
COFD(1173) = -1.87383952d+01
COFD(1174) = 3.96926341d+00
COFD(1175) = -2.16412264d-01
COFD(1176) = 6.06012078d-03
COFD(1177) = -1.87383952d+01
COFD(1178) = 3.96926341d+00
COFD(1179) = -2.16412264d-01
COFD(1180) = 6.06012078d-03
COFD(1181) = -1.92450597d+01
COFD(1182) = 5.05708283d+00
COFD(1183) = -4.35739290d-01
COFD(1184) = 1.86583205d-02
COFD(1185) = -1.93545828d+01
COFD(1186) = 5.16013126d+00
COFD(1187) = -4.46824543d-01
COFD(1188) = 1.90464887d-02
COFD(1189) = -2.10227861d+01
COFD(1190) = 5.58489889d+00
COFD(1191) = -4.81572453d-01
COFD(1192) = 1.97432291d-02
COFD(1193) = -1.47127617d+01
COFD(1194) = 2.26104173d+00
COFD(1195) = 4.18075774d-02
COFD(1196) = -6.41614212d-03
COFD(1197) = -1.95796680d+01
COFD(1198) = 5.05708283d+00
COFD(1199) = -4.35739290d-01
COFD(1200) = 1.86583205d-02
COFD(1201) = -2.10774940d+01
COFD(1202) = 5.53614847d+00
COFD(1203) = -4.86046736d-01
COFD(1204) = 2.03659188d-02
COFD(1205) = -2.14542524d+01
COFD(1206) = 5.49993732d+00
COFD(1207) = -4.62042917d-01
COFD(1208) = 1.85577413d-02
COFD(1209) = -2.13656038d+01
COFD(1210) = 5.35901292d+00
COFD(1211) = -4.36172487d-01
COFD(1212) = 1.71345319d-02
COFD(1213) = -2.09199029d+01
COFD(1214) = 5.50410224d+00
COFD(1215) = -4.82760329d-01
COFD(1216) = 2.02578587d-02
COFD(1217) = -1.97463680d+01
COFD(1218) = 5.11531505d+00
COFD(1219) = -4.42537303d-01
COFD(1220) = 1.89231395d-02
COFD(1221) = -1.66996477d+01
COFD(1222) = 3.08217426d+00
COFD(1223) = -8.33350945d-02
COFD(1224) = -3.89962983d-04
COFD(1225) = -1.67068906d+01
COFD(1226) = 3.08217426d+00
COFD(1227) = -8.33350945d-02
COFD(1228) = -3.89962983d-04
COFD(1229) = -2.10440781d+01
COFD(1230) = 5.59806282d+00
COFD(1231) = -4.87109535d-01
COFD(1232) = 2.01370226d-02
COFD(1233) = -2.13622815d+01
COFD(1234) = 5.35901292d+00
COFD(1235) = -4.36172487d-01
COFD(1236) = 1.71345319d-02
COFD(1237) = -2.13622815d+01
COFD(1238) = 5.35901292d+00
COFD(1239) = -4.36172487d-01
COFD(1240) = 1.71345319d-02
COFD(1241) = -2.13622815d+01
COFD(1242) = 5.35901292d+00
COFD(1243) = -4.36172487d-01
COFD(1244) = 1.71345319d-02
COFD(1245) = -2.13587539d+01
COFD(1246) = 5.35901292d+00
COFD(1247) = -4.36172487d-01
COFD(1248) = 1.71345319d-02
COFD(1249) = -2.10643259d+01
COFD(1250) = 5.53614847d+00
COFD(1251) = -4.86046736d-01
COFD(1252) = 2.03659188d-02
COFD(1253) = -2.12755888d+01
COFD(1254) = 5.60381989d+00
COFD(1255) = -4.91225459d-01
COFD(1256) = 2.04487844d-02
COFD(1257) = -2.13919273d+01
COFD(1258) = 5.17440955d+00
COFD(1259) = -4.04678430d-01
COFD(1260) = 1.54706350d-02
COFD(1261) = -2.13953083d+01
COFD(1262) = 5.17440955d+00
COFD(1263) = -4.04678430d-01
COFD(1264) = 1.54706350d-02
COFD(1265) = -1.87419199d+01
COFD(1266) = 3.96926341d+00
COFD(1267) = -2.16412264d-01
COFD(1268) = 6.06012078d-03
COFD(1269) = -1.87453067d+01
COFD(1270) = 3.96926341d+00
COFD(1271) = -2.16412264d-01
COFD(1272) = 6.06012078d-03
COFD(1273) = -1.18998012d+01
COFD(1274) = 2.57507000d+00
COFD(1275) = -1.24033737d-01
COFD(1276) = 5.56694959d-03
COFD(1277) = -1.46554748d+01
COFD(1278) = 3.83606243d+00
COFD(1279) = -2.86076532d-01
COFD(1280) = 1.25205829d-02
COFD(1281) = -1.43190389d+01
COFD(1282) = 3.17651319d+00
COFD(1283) = -2.02028974d-01
COFD(1284) = 8.94232502d-03
COFD(1285) = -1.55588279d+01
COFD(1286) = 3.48070094d+00
COFD(1287) = -2.40859499d-01
COFD(1288) = 1.05972514d-02
COFD(1289) = -1.43394069d+01
COFD(1290) = 3.17651319d+00
COFD(1291) = -2.02028974d-01
COFD(1292) = 8.94232502d-03
COFD(1293) = -2.06463744d+01
COFD(1294) = 5.41688482d+00
COFD(1295) = -4.73387188d-01
COFD(1296) = 1.99280175d-02
COFD(1297) = -1.55666415d+01
COFD(1298) = 3.48070094d+00
COFD(1299) = -2.40859499d-01
COFD(1300) = 1.05972514d-02
COFD(1301) = -1.55741053d+01
COFD(1302) = 3.48070094d+00
COFD(1303) = -2.40859499d-01
COFD(1304) = 1.05972514d-02
COFD(1305) = -1.40749320d+01
COFD(1306) = 3.05837263d+00
COFD(1307) = -1.86672802d-01
COFD(1308) = 8.27575734d-03
COFD(1309) = -1.42473439d+01
COFD(1310) = 3.17651319d+00
COFD(1311) = -2.02028974d-01
COFD(1312) = 8.94232502d-03
COFD(1313) = -1.63301444d+01
COFD(1314) = 3.82388595d+00
COFD(1315) = -2.84480724d-01
COFD(1316) = 1.24506311d-02
COFD(1317) = -1.63301444d+01
COFD(1318) = 3.82388595d+00
COFD(1319) = -2.84480724d-01
COFD(1320) = 1.24506311d-02
COFD(1321) = -1.63542394d+01
COFD(1322) = 3.82388595d+00
COFD(1323) = -2.84480724d-01
COFD(1324) = 1.24506311d-02
COFD(1325) = -1.62775714d+01
COFD(1326) = 3.79163564d+00
COFD(1327) = -2.80257365d-01
COFD(1328) = 1.22656902d-02
COFD(1329) = -1.52792891d+01
COFD(1330) = 3.36790500d+00
COFD(1331) = -2.26321740d-01
COFD(1332) = 9.97135055d-03
COFD(1333) = -1.84777607d+01
COFD(1334) = 4.49330851d+00
COFD(1335) = -3.68208715d-01
COFD(1336) = 1.59565402d-02
COFD(1337) = -2.08277598d+01
COFD(1338) = 5.35267674d+00
COFD(1339) = -4.69010505d-01
COFD(1340) = 1.98979152d-02
COFD(1341) = -2.08367725d+01
COFD(1342) = 5.35267674d+00
COFD(1343) = -4.69010505d-01
COFD(1344) = 1.98979152d-02
COFD(1345) = -2.02637994d+01
COFD(1346) = 5.14984081d+00
COFD(1347) = -4.46093018d-01
COFD(1348) = 1.90396647d-02
COFD(1349) = -2.02637994d+01
COFD(1350) = 5.14984081d+00
COFD(1351) = -4.46093018d-01
COFD(1352) = 1.90396647d-02
COFD(1353) = -2.07672833d+01
COFD(1354) = 5.32244593d+00
COFD(1355) = -4.65829403d-01
COFD(1356) = 1.97895274d-02
COFD(1357) = -1.78792605d+01
COFD(1358) = 4.29613154d+00
COFD(1359) = -3.44012526d-01
COFD(1360) = 1.49643715d-02
COFD(1361) = -1.78903913d+01
COFD(1362) = 4.29613154d+00
COFD(1363) = -3.44012526d-01
COFD(1364) = 1.49643715d-02
COFD(1365) = -1.79009181d+01
COFD(1366) = 4.29613154d+00
COFD(1367) = -3.44012526d-01
COFD(1368) = 1.49643715d-02
COFD(1369) = -1.89616623d+01
COFD(1370) = 4.68595732d+00
COFD(1371) = -3.91842840d-01
COFD(1372) = 1.69262542d-02
COFD(1373) = -1.86409139d+01
COFD(1374) = 4.53572533d+00
COFD(1375) = -3.73386925d-01
COFD(1376) = 1.61678881d-02
COFD(1377) = -1.86499071d+01
COFD(1378) = 4.53572533d+00
COFD(1379) = -3.73386925d-01
COFD(1380) = 1.61678881d-02
COFD(1381) = -1.64255964d+01
COFD(1382) = 3.89309916d+00
COFD(1383) = -2.93528188d-01
COFD(1384) = 1.28463177d-02
COFD(1385) = -2.05272328d+01
COFD(1386) = 5.18417470d+00
COFD(1387) = -4.49491573d-01
COFD(1388) = 1.91438508d-02
COFD(1389) = -2.05272328d+01
COFD(1390) = 5.18417470d+00
COFD(1391) = -4.49491573d-01
COFD(1392) = 1.91438508d-02
COFD(1393) = -1.41300955d+01
COFD(1394) = 3.05837263d+00
COFD(1395) = -1.86672802d-01
COFD(1396) = 8.27575734d-03
COFD(1397) = -1.42649477d+01
COFD(1398) = 3.17651319d+00
COFD(1399) = -2.02028974d-01
COFD(1400) = 8.94232502d-03
COFD(1401) = -1.42870488d+01
COFD(1402) = 3.17651319d+00
COFD(1403) = -2.02028974d-01
COFD(1404) = 8.94232502d-03
COFD(1405) = -2.03661760d+01
COFD(1406) = 5.32029196d+00
COFD(1407) = -4.65578114d-01
COFD(1408) = 1.97797257d-02
COFD(1409) = -1.44985622d+01
COFD(1410) = 3.05837263d+00
COFD(1411) = -1.86672802d-01
COFD(1412) = 8.27575734d-03
COFD(1413) = -1.52669189d+01
COFD(1414) = 3.35922578d+00
COFD(1415) = -2.25181399d-01
COFD(1416) = 9.92132878d-03
COFD(1417) = -1.77428218d+01
COFD(1418) = 4.25438185d+00
COFD(1419) = -3.39084808d-01
COFD(1420) = 1.47709916d-02
COFD(1421) = -1.83397596d+01
COFD(1422) = 4.43878381d+00
COFD(1423) = -3.61697624d-01
COFD(1424) = 1.56975581d-02
COFD(1425) = -1.58043047d+01
COFD(1426) = 3.57122377d+00
COFD(1427) = -2.52409987d-01
COFD(1428) = 1.10900562d-02
COFD(1429) = -1.46206548d+01
COFD(1430) = 3.10987301d+00
COFD(1431) = -1.93373245d-01
COFD(1432) = 8.56678670d-03
COFD(1433) = -2.10721980d+01
COFD(1434) = 5.41544337d+00
COFD(1435) = -4.73380953d-01
COFD(1436) = 1.99350119d-02
COFD(1437) = -2.10821814d+01
COFD(1438) = 5.41544337d+00
COFD(1439) = -4.73380953d-01
COFD(1440) = 1.99350119d-02
COFD(1441) = -1.64256119d+01
COFD(1442) = 3.89309916d+00
COFD(1443) = -2.93528188d-01
COFD(1444) = 1.28463177d-02
COFD(1445) = -1.83348653d+01
COFD(1446) = 4.43878381d+00
COFD(1447) = -3.61697624d-01
COFD(1448) = 1.56975581d-02
COFD(1449) = -1.83348653d+01
COFD(1450) = 4.43878381d+00
COFD(1451) = -3.61697624d-01
COFD(1452) = 1.56975581d-02
COFD(1453) = -1.83348653d+01
COFD(1454) = 4.43878381d+00
COFD(1455) = -3.61697624d-01
COFD(1456) = 1.56975581d-02
COFD(1457) = -1.83296851d+01
COFD(1458) = 4.43878381d+00
COFD(1459) = -3.61697624d-01
COFD(1460) = 1.56975581d-02
COFD(1461) = -1.52486273d+01
COFD(1462) = 3.35922578d+00
COFD(1463) = -2.25181399d-01
COFD(1464) = 9.92132878d-03
COFD(1465) = -1.62465583d+01
COFD(1466) = 3.72612300d+00
COFD(1467) = -2.71663673d-01
COFD(1468) = 1.18889643d-02
COFD(1469) = -1.91313643d+01
COFD(1470) = 4.61801405d+00
COFD(1471) = -3.83535652d-01
COFD(1472) = 1.65862513d-02
COFD(1473) = -1.91363464d+01
COFD(1474) = 4.61801405d+00
COFD(1475) = -3.83535652d-01
COFD(1476) = 1.65862513d-02
COFD(1477) = -2.05324091d+01
COFD(1478) = 5.18417470d+00
COFD(1479) = -4.49491573d-01
COFD(1480) = 1.91438508d-02
COFD(1481) = -2.05373990d+01
COFD(1482) = 5.18417470d+00
COFD(1483) = -4.49491573d-01
COFD(1484) = 1.91438508d-02
COFD(1485) = -1.19006548d+01
COFD(1486) = 2.57507000d+00
COFD(1487) = -1.24033737d-01
COFD(1488) = 5.56694959d-03
COFD(1489) = -1.46559141d+01
COFD(1490) = 3.83606243d+00
COFD(1491) = -2.86076532d-01
COFD(1492) = 1.25205829d-02
COFD(1493) = -1.43238998d+01
COFD(1494) = 3.17651319d+00
COFD(1495) = -2.02028974d-01
COFD(1496) = 8.94232502d-03
COFD(1497) = -1.55661750d+01
COFD(1498) = 3.48070094d+00
COFD(1499) = -2.40859499d-01
COFD(1500) = 1.05972514d-02
COFD(1501) = -1.43444709d+01
COFD(1502) = 3.17651319d+00
COFD(1503) = -2.02028974d-01
COFD(1504) = 8.94232502d-03
COFD(1505) = -2.06516336d+01
COFD(1506) = 5.41688482d+00
COFD(1507) = -4.73387188d-01
COFD(1508) = 1.99280175d-02
COFD(1509) = -1.55741053d+01
COFD(1510) = 3.48070094d+00
COFD(1511) = -2.40859499d-01
COFD(1512) = 1.05972514d-02
COFD(1513) = -1.55816822d+01
COFD(1514) = 3.48070094d+00
COFD(1515) = -2.40859499d-01
COFD(1516) = 1.05972514d-02
COFD(1517) = -1.40789009d+01
COFD(1518) = 3.05837263d+00
COFD(1519) = -1.86672802d-01
COFD(1520) = 8.27575734d-03
COFD(1521) = -1.42515527d+01
COFD(1522) = 3.17651319d+00
COFD(1523) = -2.02028974d-01
COFD(1524) = 8.94232502d-03
COFD(1525) = -1.63345829d+01
COFD(1526) = 3.82388595d+00
COFD(1527) = -2.84480724d-01
COFD(1528) = 1.24506311d-02
COFD(1529) = -1.63345829d+01
COFD(1530) = 3.82388595d+00
COFD(1531) = -2.84480724d-01
COFD(1532) = 1.24506311d-02
COFD(1533) = -1.63588981d+01
COFD(1534) = 3.82388595d+00
COFD(1535) = -2.84480724d-01
COFD(1536) = 1.24506311d-02
COFD(1537) = -1.62824412d+01
COFD(1538) = 3.79163564d+00
COFD(1539) = -2.80257365d-01
COFD(1540) = 1.22656902d-02
COFD(1541) = -1.52861376d+01
COFD(1542) = 3.36790500d+00
COFD(1543) = -2.26321740d-01
COFD(1544) = 9.97135055d-03
COFD(1545) = -1.84863000d+01
COFD(1546) = 4.49330851d+00
COFD(1547) = -3.68208715d-01
COFD(1548) = 1.59565402d-02
COFD(1549) = -2.08347403d+01
COFD(1550) = 5.35267674d+00
COFD(1551) = -4.69010505d-01
COFD(1552) = 1.98979152d-02
COFD(1553) = -2.08438809d+01
COFD(1554) = 5.35267674d+00
COFD(1555) = -4.69010505d-01
COFD(1556) = 1.98979152d-02
COFD(1557) = -2.02710316d+01
COFD(1558) = 5.14984081d+00
COFD(1559) = -4.46093018d-01
COFD(1560) = 1.90396647d-02
COFD(1561) = -2.02710316d+01
COFD(1562) = 5.14984081d+00
COFD(1563) = -4.46093018d-01
COFD(1564) = 1.90396647d-02
COFD(1565) = -2.07746356d+01
COFD(1566) = 5.32244593d+00
COFD(1567) = -4.65829403d-01
COFD(1568) = 1.97895274d-02
COFD(1569) = -1.78856918d+01
COFD(1570) = 4.29613154d+00
COFD(1571) = -3.44012526d-01
COFD(1572) = 1.49643715d-02
COFD(1573) = -1.78969684d+01
COFD(1574) = 4.29613154d+00
COFD(1575) = -3.44012526d-01
COFD(1576) = 1.49643715d-02
COFD(1577) = -1.79076361d+01
COFD(1578) = 4.29613154d+00
COFD(1579) = -3.44012526d-01
COFD(1580) = 1.49643715d-02
COFD(1581) = -1.89685165d+01
COFD(1582) = 4.68595732d+00
COFD(1583) = -3.91842840d-01
COFD(1584) = 1.69262542d-02
COFD(1585) = -1.86479000d+01
COFD(1586) = 4.53572533d+00
COFD(1587) = -3.73386925d-01
COFD(1588) = 1.61678881d-02
COFD(1589) = -1.86570209d+01
COFD(1590) = 4.53572533d+00
COFD(1591) = -3.73386925d-01
COFD(1592) = 1.61678881d-02
COFD(1593) = -1.64338757d+01
COFD(1594) = 3.89309916d+00
COFD(1595) = -2.93528188d-01
COFD(1596) = 1.28463177d-02
COFD(1597) = -2.05356023d+01
COFD(1598) = 5.18417470d+00
COFD(1599) = -4.49491573d-01
COFD(1600) = 1.91438508d-02
COFD(1601) = -2.05356023d+01
COFD(1602) = 5.18417470d+00
COFD(1603) = -4.49491573d-01
COFD(1604) = 1.91438508d-02
COFD(1605) = -1.41345294d+01
COFD(1606) = 3.05837263d+00
COFD(1607) = -1.86672802d-01
COFD(1608) = 8.27575734d-03
COFD(1609) = -1.42696020d+01
COFD(1610) = 3.17651319d+00
COFD(1611) = -2.02028974d-01
COFD(1612) = 8.94232502d-03
COFD(1613) = -1.42919145d+01
COFD(1614) = 3.17651319d+00
COFD(1615) = -2.02028974d-01
COFD(1616) = 8.94232502d-03
COFD(1617) = -2.03712446d+01
COFD(1618) = 5.32029196d+00
COFD(1619) = -4.65578114d-01
COFD(1620) = 1.97797257d-02
COFD(1621) = -1.45055431d+01
COFD(1622) = 3.05837263d+00
COFD(1623) = -1.86672802d-01
COFD(1624) = 8.27575734d-03
COFD(1625) = -1.52740247d+01
COFD(1626) = 3.35922578d+00
COFD(1627) = -2.25181399d-01
COFD(1628) = 9.92132878d-03
COFD(1629) = -1.77515242d+01
COFD(1630) = 4.25438185d+00
COFD(1631) = -3.39084808d-01
COFD(1632) = 1.47709916d-02
COFD(1633) = -1.83482991d+01
COFD(1634) = 4.43878381d+00
COFD(1635) = -3.61697624d-01
COFD(1636) = 1.56975581d-02
COFD(1637) = -1.58115345d+01
COFD(1638) = 3.57122377d+00
COFD(1639) = -2.52409987d-01
COFD(1640) = 1.10900562d-02
COFD(1641) = -1.46272289d+01
COFD(1642) = 3.10987301d+00
COFD(1643) = -1.93373245d-01
COFD(1644) = 8.56678670d-03
COFD(1645) = -2.10789131d+01
COFD(1646) = 5.41544337d+00
COFD(1647) = -4.73380953d-01
COFD(1648) = 1.99350119d-02
COFD(1649) = -2.10890329d+01
COFD(1650) = 5.41544337d+00
COFD(1651) = -4.73380953d-01
COFD(1652) = 1.99350119d-02
COFD(1653) = -1.64338914d+01
COFD(1654) = 3.89309916d+00
COFD(1655) = -2.93528188d-01
COFD(1656) = 1.28463177d-02
COFD(1657) = -1.83433209d+01
COFD(1658) = 4.43878381d+00
COFD(1659) = -3.61697624d-01
COFD(1660) = 1.56975581d-02
COFD(1661) = -1.83433209d+01
COFD(1662) = 4.43878381d+00
COFD(1663) = -3.61697624d-01
COFD(1664) = 1.56975581d-02
COFD(1665) = -1.83433209d+01
COFD(1666) = 4.43878381d+00
COFD(1667) = -3.61697624d-01
COFD(1668) = 1.56975581d-02
COFD(1669) = -1.83380528d+01
COFD(1670) = 4.43878381d+00
COFD(1671) = -3.61697624d-01
COFD(1672) = 1.56975581d-02
COFD(1673) = -1.52554761d+01
COFD(1674) = 3.35922578d+00
COFD(1675) = -2.25181399d-01
COFD(1676) = 9.92132878d-03
COFD(1677) = -1.62547381d+01
COFD(1678) = 3.72612300d+00
COFD(1679) = -2.71663673d-01
COFD(1680) = 1.18889643d-02
COFD(1681) = -1.91398254d+01
COFD(1682) = 4.61801405d+00
COFD(1683) = -3.83535652d-01
COFD(1684) = 1.65862513d-02
COFD(1685) = -1.91448929d+01
COFD(1686) = 4.61801405d+00
COFD(1687) = -3.83535652d-01
COFD(1688) = 1.65862513d-02
COFD(1689) = -2.05408665d+01
COFD(1690) = 5.18417470d+00
COFD(1691) = -4.49491573d-01
COFD(1692) = 1.91438508d-02
COFD(1693) = -2.05459419d+01
COFD(1694) = 5.18417470d+00
COFD(1695) = -4.49491573d-01
COFD(1696) = 1.91438508d-02
COFD(1697) = -1.08369481d+01
COFD(1698) = 2.19094415d+00
COFD(1699) = -7.11992510d-02
COFD(1700) = 3.14105973d-03
COFD(1701) = -1.32467319d+01
COFD(1702) = 3.34156587d+00
COFD(1703) = -2.22853306d-01
COFD(1704) = 9.81883417d-03
COFD(1705) = -1.30610083d+01
COFD(1706) = 2.80913567d+00
COFD(1707) = -1.54536855d-01
COFD(1708) = 6.89359313d-03
COFD(1709) = -1.40707473d+01
COFD(1710) = 3.05837263d+00
COFD(1711) = -1.86672802d-01
COFD(1712) = 8.27575734d-03
COFD(1713) = -1.30738796d+01
COFD(1714) = 2.80913567d+00
COFD(1715) = -1.54536855d-01
COFD(1716) = 6.89359313d-03
COFD(1717) = -1.92003817d+01
COFD(1718) = 5.05708283d+00
COFD(1719) = -4.35739290d-01
COFD(1720) = 1.86583205d-02
COFD(1721) = -1.40749320d+01
COFD(1722) = 3.05837263d+00
COFD(1723) = -1.86672802d-01
COFD(1724) = 8.27575734d-03
COFD(1725) = -1.40789009d+01
COFD(1726) = 3.05837263d+00
COFD(1727) = -1.86672802d-01
COFD(1728) = 8.27575734d-03
COFD(1729) = -1.29573076d+01
COFD(1730) = 2.73155251d+00
COFD(1731) = -1.44594082d-01
COFD(1732) = 6.46883252d-03
COFD(1733) = -1.30141899d+01
COFD(1734) = 2.80913567d+00
COFD(1735) = -1.54536855d-01
COFD(1736) = 6.89359313d-03
COFD(1737) = -1.47558850d+01
COFD(1738) = 3.33113524d+00
COFD(1739) = -2.21479057d-01
COFD(1740) = 9.75837737d-03
COFD(1741) = -1.47558850d+01
COFD(1742) = 3.33113524d+00
COFD(1743) = -2.21479057d-01
COFD(1744) = 9.75837737d-03
COFD(1745) = -1.47715918d+01
COFD(1746) = 3.33113524d+00
COFD(1747) = -2.21479057d-01
COFD(1748) = 9.75837737d-03
COFD(1749) = -1.47048024d+01
COFD(1750) = 3.30594991d+00
COFD(1751) = -2.18182207d-01
COFD(1752) = 9.61429447d-03
COFD(1753) = -1.38862839d+01
COFD(1754) = 2.97564184d+00
COFD(1755) = -1.76025309d-01
COFD(1756) = 7.81869993d-03
COFD(1757) = -1.67390425d+01
COFD(1758) = 4.00828594d+00
COFD(1759) = -3.08414344d-01
COFD(1760) = 1.34907430d-02
COFD(1761) = -1.90526941d+01
COFD(1762) = 4.86821670d+00
COFD(1763) = -4.13144121d-01
COFD(1764) = 1.77546701d-02
COFD(1765) = -1.90576320d+01
COFD(1766) = 4.86821670d+00
COFD(1767) = -4.13144121d-01
COFD(1768) = 1.77546701d-02
COFD(1769) = -1.85157414d+01
COFD(1770) = 4.67076124d+00
COFD(1771) = -3.90022427d-01
COFD(1772) = 1.68533953d-02
COFD(1773) = -1.85157414d+01
COFD(1774) = 4.67076124d+00
COFD(1775) = -3.90022427d-01
COFD(1776) = 1.68533953d-02
COFD(1777) = -1.89671752d+01
COFD(1778) = 4.83076737d+00
COFD(1779) = -4.08802573d-01
COFD(1780) = 1.75875241d-02
COFD(1781) = -1.61038806d+01
COFD(1782) = 3.75910622d+00
COFD(1783) = -2.75986578d-01
COFD(1784) = 1.20782843d-02
COFD(1785) = -1.61101966d+01
COFD(1786) = 3.75910622d+00
COFD(1787) = -2.75986578d-01
COFD(1788) = 1.20782843d-02
COFD(1789) = -1.61161138d+01
COFD(1790) = 3.75910622d+00
COFD(1791) = -2.75986578d-01
COFD(1792) = 1.20782843d-02
COFD(1793) = -1.71884218d+01
COFD(1794) = 4.17190426d+00
COFD(1795) = -3.28894681d-01
COFD(1796) = 1.43498101d-02
COFD(1797) = -1.69491115d+01
COFD(1798) = 4.05099737d+00
COFD(1799) = -3.13841660d-01
COFD(1800) = 1.37218854d-02
COFD(1801) = -1.69540369d+01
COFD(1802) = 4.05099737d+00
COFD(1803) = -3.13841660d-01
COFD(1804) = 1.37218854d-02
COFD(1805) = -1.46907028d+01
COFD(1806) = 3.39229020d+00
COFD(1807) = -2.29520232d-01
COFD(1808) = 1.01114311d-02
COFD(1809) = -1.87688110d+01
COFD(1810) = 4.71729964d+00
COFD(1811) = -3.95432573d-01
COFD(1812) = 1.70623691d-02
COFD(1813) = -1.87688110d+01
COFD(1814) = 4.71729964d+00
COFD(1815) = -3.95432573d-01
COFD(1816) = 1.70623691d-02
COFD(1817) = -1.29942577d+01
COFD(1818) = 2.73155251d+00
COFD(1819) = -1.44594082d-01
COFD(1820) = 6.46883252d-03
COFD(1821) = -1.30137956d+01
COFD(1822) = 2.80913567d+00
COFD(1823) = -1.54536855d-01
COFD(1824) = 6.89359313d-03
COFD(1825) = -1.30279742d+01
COFD(1826) = 2.80913567d+00
COFD(1827) = -1.54536855d-01
COFD(1828) = 6.89359313d-03
COFD(1829) = -1.86394179d+01
COFD(1830) = 4.82909392d+00
COFD(1831) = -4.08610711d-01
COFD(1832) = 1.75802236d-02
COFD(1833) = -1.32768507d+01
COFD(1834) = 2.73155251d+00
COFD(1835) = -1.44594082d-01
COFD(1836) = 6.46883252d-03
COFD(1837) = -1.38762133d+01
COFD(1838) = 2.97137588d+00
COFD(1839) = -1.75491257d-01
COFD(1840) = 7.79646773d-03
COFD(1841) = -1.58184910d+01
COFD(1842) = 3.68407693d+00
COFD(1843) = -2.66228170d-01
COFD(1844) = 1.16542305d-02
COFD(1845) = -1.65721778d+01
COFD(1846) = 3.93849401d+00
COFD(1847) = -2.99416642d-01
COFD(1848) = 1.31020815d-02
COFD(1849) = -1.43049857d+01
COFD(1850) = 3.14480429d+00
COFD(1851) = -1.97906290d-01
COFD(1852) = 8.76325718d-03
COFD(1853) = -1.33664987d+01
COFD(1854) = 2.76469604d+00
COFD(1855) = -1.48843508d-01
COFD(1856) = 6.65045597d-03
COFD(1857) = -1.95705999d+01
COFD(1858) = 5.04961455d+00
COFD(1859) = -4.34856114d-01
COFD(1860) = 1.86234058d-02
COFD(1861) = -1.95761623d+01
COFD(1862) = 5.04961455d+00
COFD(1863) = -4.34856114d-01
COFD(1864) = 1.86234058d-02
COFD(1865) = -1.46907107d+01
COFD(1866) = 3.39229020d+00
COFD(1867) = -2.29520232d-01
COFD(1868) = 1.01114311d-02
COFD(1869) = -1.65697233d+01
COFD(1870) = 3.93849401d+00
COFD(1871) = -2.99416642d-01
COFD(1872) = 1.31020815d-02
COFD(1873) = -1.65697233d+01
COFD(1874) = 3.93849401d+00
COFD(1875) = -2.99416642d-01
COFD(1876) = 1.31020815d-02
COFD(1877) = -1.65697233d+01
COFD(1878) = 3.93849401d+00
COFD(1879) = -2.99416642d-01
COFD(1880) = 1.31020815d-02
COFD(1881) = -1.65671124d+01
COFD(1882) = 3.93849401d+00
COFD(1883) = -2.99416642d-01
COFD(1884) = 1.31020815d-02
COFD(1885) = -1.38661480d+01
COFD(1886) = 2.97137588d+00
COFD(1887) = -1.75491257d-01
COFD(1888) = 7.79646773d-03
COFD(1889) = -1.46461881d+01
COFD(1890) = 3.27505697d+00
COFD(1891) = -2.14306851d-01
COFD(1892) = 9.45219335d-03
COFD(1893) = -1.73541877d+01
COFD(1894) = 4.11838823d+00
COFD(1895) = -3.22329602d-01
COFD(1896) = 1.40802468d-02
COFD(1897) = -1.73566853d+01
COFD(1898) = 4.11838823d+00
COFD(1899) = -3.22329602d-01
COFD(1900) = 1.40802468d-02
COFD(1901) = -1.87714196d+01
COFD(1902) = 4.71729964d+00
COFD(1903) = -3.95432573d-01
COFD(1904) = 1.70623691d-02
COFD(1905) = -1.87739217d+01
COFD(1906) = 4.71729964d+00
COFD(1907) = -3.95432573d-01
COFD(1908) = 1.70623691d-02
COFD(1909) = -1.09469245d+01
COFD(1910) = 2.30836460d+00
COFD(1911) = -8.76339315d-02
COFD(1912) = 3.90878445d-03
COFD(1913) = -1.34162893d+01
COFD(1914) = 3.48624238d+00
COFD(1915) = -2.41554467d-01
COFD(1916) = 1.06263545d-02
COFD(1917) = -1.31551788d+01
COFD(1918) = 2.90778936d+00
COFD(1919) = -1.67388544d-01
COFD(1920) = 7.45220609d-03
COFD(1921) = -1.42429085d+01
COFD(1922) = 3.17651319d+00
COFD(1923) = -2.02028974d-01
COFD(1924) = 8.94232502d-03
COFD(1925) = -1.31686537d+01
COFD(1926) = 2.90778936d+00
COFD(1927) = -1.67388544d-01
COFD(1928) = 7.45220609d-03
COFD(1929) = -1.93521390d+01
COFD(1930) = 5.16013126d+00
COFD(1931) = -4.46824543d-01
COFD(1932) = 1.90464887d-02
COFD(1933) = -1.42473439d+01
COFD(1934) = 3.17651319d+00
COFD(1935) = -2.02028974d-01
COFD(1936) = 8.94232502d-03
COFD(1937) = -1.42515527d+01
COFD(1938) = 3.17651319d+00
COFD(1939) = -2.02028974d-01
COFD(1940) = 8.94232502d-03
COFD(1941) = -1.30141899d+01
COFD(1942) = 2.80913567d+00
COFD(1943) = -1.54536855d-01
COFD(1944) = 6.89359313d-03
COFD(1945) = -1.31062967d+01
COFD(1946) = 2.90778936d+00
COFD(1947) = -1.67388544d-01
COFD(1948) = 7.45220609d-03
COFD(1949) = -1.50076254d+01
COFD(1950) = 3.47945612d+00
COFD(1951) = -2.40703722d-01
COFD(1952) = 1.05907441d-02
COFD(1953) = -1.50076254d+01
COFD(1954) = 3.47945612d+00
COFD(1955) = -2.40703722d-01
COFD(1956) = 1.05907441d-02
COFD(1957) = -1.50240272d+01
COFD(1958) = 3.47945612d+00
COFD(1959) = -2.40703722d-01
COFD(1960) = 1.05907441d-02
COFD(1961) = -1.49727799d+01
COFD(1962) = 3.46140064d+00
COFD(1963) = -2.38440092d-01
COFD(1964) = 1.04960087d-02
COFD(1965) = -1.40318948d+01
COFD(1966) = 3.08120012d+00
COFD(1967) = -1.89629903d-01
COFD(1968) = 8.40361952d-03
COFD(1969) = -1.69758891d+01
COFD(1970) = 4.14240922d+00
COFD(1971) = -3.25239774d-01
COFD(1972) = 1.41980687d-02
COFD(1973) = -1.93624931d+01
COFD(1974) = 5.02567894d+00
COFD(1975) = -4.32045169d-01
COFD(1976) = 1.85132214d-02
COFD(1977) = -1.93677186d+01
COFD(1978) = 5.02567894d+00
COFD(1979) = -4.32045169d-01
COFD(1980) = 1.85132214d-02
COFD(1981) = -1.87476063d+01
COFD(1982) = 4.79683898d+00
COFD(1983) = -4.04829719d-01
COFD(1984) = 1.74325475d-02
COFD(1985) = -1.87476063d+01
COFD(1986) = 4.79683898d+00
COFD(1987) = -4.04829719d-01
COFD(1988) = 1.74325475d-02
COFD(1989) = -1.92654138d+01
COFD(1990) = 4.98286777d+00
COFD(1991) = -4.26970814d-01
COFD(1992) = 1.83122917d-02
COFD(1993) = -1.64758697d+01
COFD(1994) = 3.95035840d+00
COFD(1995) = -3.00959418d-01
COFD(1996) = 1.31692593d-02
COFD(1997) = -1.64825368d+01
COFD(1998) = 3.95035840d+00
COFD(1999) = -3.00959418d-01
COFD(2000) = 1.31692593d-02
COFD(2001) = -1.64887871d+01
COFD(2002) = 3.95035840d+00
COFD(2003) = -3.00959418d-01
COFD(2004) = 1.31692593d-02
COFD(2005) = -1.74111692d+01
COFD(2006) = 4.29676909d+00
COFD(2007) = -3.44085306d-01
COFD(2008) = 1.49671135d-02
COFD(2009) = -1.71808106d+01
COFD(2010) = 4.17889917d+00
COFD(2011) = -3.29752510d-01
COFD(2012) = 1.43850275d-02
COFD(2013) = -1.71860230d+01
COFD(2014) = 4.17889917d+00
COFD(2015) = -3.29752510d-01
COFD(2016) = 1.43850275d-02
COFD(2017) = -1.48738066d+01
COFD(2018) = 3.52327209d+00
COFD(2019) = -2.46286208d-01
COFD(2020) = 1.08285963d-02
COFD(2021) = -1.90116191d+01
COFD(2022) = 4.84384483d+00
COFD(2023) = -4.10265575d-01
COFD(2024) = 1.76414287d-02
COFD(2025) = -1.90116191d+01
COFD(2026) = 4.84384483d+00
COFD(2027) = -4.10265575d-01
COFD(2028) = 1.76414287d-02
COFD(2029) = -1.30526867d+01
COFD(2030) = 2.80913567d+00
COFD(2031) = -1.54536855d-01
COFD(2032) = 6.89359313d-03
COFD(2033) = -1.31039806d+01
COFD(2034) = 2.90778936d+00
COFD(2035) = -1.67388544d-01
COFD(2036) = 7.45220609d-03
COFD(2037) = -1.31188060d+01
COFD(2038) = 2.90778936d+00
COFD(2039) = -1.67388544d-01
COFD(2040) = 7.45220609d-03
COFD(2041) = -1.89111699d+01
COFD(2042) = 4.98075560d+00
COFD(2043) = -4.26721620d-01
COFD(2044) = 1.83024823d-02
COFD(2045) = -1.33548788d+01
COFD(2046) = 2.80913567d+00
COFD(2047) = -1.54536855d-01
COFD(2048) = 6.89359313d-03
COFD(2049) = -1.40183333d+01
COFD(2050) = 3.07549274d+00
COFD(2051) = -1.88889344d-01
COFD(2052) = 8.37152866d-03
COFD(2053) = -1.61987855d+01
COFD(2054) = 3.88250968d+00
COFD(2055) = -2.92155848d-01
COFD(2056) = 1.27867850d-02
COFD(2057) = -1.68514118d+01
COFD(2058) = 4.09077642d+00
COFD(2059) = -3.18894990d-01
COFD(2060) = 1.39371445d-02
COFD(2061) = -1.44362397d+01
COFD(2062) = 3.24450689d+00
COFD(2063) = -2.10570734d-01
COFD(2064) = 9.30026771d-03
COFD(2065) = -1.34740249d+01
COFD(2066) = 2.85448298d+00
COFD(2067) = -1.60491387d-01
COFD(2068) = 7.15451323d-03
COFD(2069) = -1.97672003d+01
COFD(2070) = 5.15754807d+00
COFD(2071) = -4.46654084d-01
COFD(2072) = 1.90458977d-02
COFD(2073) = -1.97730795d+01
COFD(2074) = 5.15754807d+00
COFD(2075) = -4.46654084d-01
COFD(2076) = 1.90458977d-02
COFD(2077) = -1.48738150d+01
COFD(2078) = 3.52327209d+00
COFD(2079) = -2.46286208d-01
COFD(2080) = 1.08285963d-02
COFD(2081) = -1.68487987d+01
COFD(2082) = 4.09077642d+00
COFD(2083) = -3.18894990d-01
COFD(2084) = 1.39371445d-02
COFD(2085) = -1.68487987d+01
COFD(2086) = 4.09077642d+00
COFD(2087) = -3.18894990d-01
COFD(2088) = 1.39371445d-02
COFD(2089) = -1.68487987d+01
COFD(2090) = 4.09077642d+00
COFD(2091) = -3.18894990d-01
COFD(2092) = 1.39371445d-02
COFD(2093) = -1.68460201d+01
COFD(2094) = 4.09077642d+00
COFD(2095) = -3.18894990d-01
COFD(2096) = 1.39371445d-02
COFD(2097) = -1.40076852d+01
COFD(2098) = 3.07549274d+00
COFD(2099) = -1.88889344d-01
COFD(2100) = 8.37152866d-03
COFD(2101) = -1.48853569d+01
COFD(2102) = 3.41988961d+00
COFD(2103) = -2.33128386d-01
COFD(2104) = 1.02689994d-02
COFD(2105) = -1.76069153d+01
COFD(2106) = 4.24719726d+00
COFD(2107) = -3.38206061d-01
COFD(2108) = 1.47350654d-02
COFD(2109) = -1.76095743d+01
COFD(2110) = 4.24719726d+00
COFD(2111) = -3.38206061d-01
COFD(2112) = 1.47350654d-02
COFD(2113) = -1.90143953d+01
COFD(2114) = 4.84384483d+00
COFD(2115) = -4.10265575d-01
COFD(2116) = 1.76414287d-02
COFD(2117) = -1.90170590d+01
COFD(2118) = 4.84384483d+00
COFD(2119) = -4.10265575d-01
COFD(2120) = 1.76414287d-02
COFD(2121) = -1.25098960d+01
COFD(2122) = 2.77873601d+00
COFD(2123) = -1.50637360d-01
COFD(2124) = 6.72684281d-03
COFD(2125) = -1.57972369d+01
COFD(2126) = 4.22225052d+00
COFD(2127) = -3.35156428d-01
COFD(2128) = 1.46104855d-02
COFD(2129) = -1.50584249d+01
COFD(2130) = 3.47945612d+00
COFD(2131) = -2.40703722d-01
COFD(2132) = 1.05907441d-02
COFD(2133) = -1.63254691d+01
COFD(2134) = 3.82388595d+00
COFD(2135) = -2.84480724d-01
COFD(2136) = 1.24506311d-02
COFD(2137) = -1.50724636d+01
COFD(2138) = 3.47945612d+00
COFD(2139) = -2.40703722d-01
COFD(2140) = 1.05907441d-02
COFD(2141) = -2.12639214d+01
COFD(2142) = 5.61184117d+00
COFD(2143) = -4.90532156d-01
COFD(2144) = 2.03507922d-02
COFD(2145) = -1.63301444d+01
COFD(2146) = 3.82388595d+00
COFD(2147) = -2.84480724d-01
COFD(2148) = 1.24506311d-02
COFD(2149) = -1.63345829d+01
COFD(2150) = 3.82388595d+00
COFD(2151) = -2.84480724d-01
COFD(2152) = 1.24506311d-02
COFD(2153) = -1.47558850d+01
COFD(2154) = 3.33113524d+00
COFD(2155) = -2.21479057d-01
COFD(2156) = 9.75837737d-03
COFD(2157) = -1.50076254d+01
COFD(2158) = 3.47945612d+00
COFD(2159) = -2.40703722d-01
COFD(2160) = 1.05907441d-02
COFD(2161) = -1.73027557d+01
COFD(2162) = 4.21416723d+00
COFD(2163) = -3.34163932d-01
COFD(2164) = 1.45697432d-02
COFD(2165) = -1.73027557d+01
COFD(2166) = 4.21416723d+00
COFD(2167) = -3.34163932d-01
COFD(2168) = 1.45697432d-02
COFD(2169) = -1.73198034d+01
COFD(2170) = 4.21416723d+00
COFD(2171) = -3.34163932d-01
COFD(2172) = 1.45697432d-02
COFD(2173) = -1.72556729d+01
COFD(2174) = 4.19029808d+00
COFD(2175) = -3.31177076d-01
COFD(2176) = 1.44446234d-02
COFD(2177) = -1.59634533d+01
COFD(2178) = 3.67388294d+00
COFD(2179) = -2.64990709d-01
COFD(2180) = 1.16042706d-02
COFD(2181) = -1.93015555d+01
COFD(2182) = 4.85015581d+00
COFD(2183) = -4.10945109d-01
COFD(2184) = 1.76651398d-02
COFD(2185) = -2.14160703d+01
COFD(2186) = 5.56531152d+00
COFD(2187) = -4.88789821d-01
COFD(2188) = 2.04437116d-02
COFD(2189) = -2.14215700d+01
COFD(2190) = 5.56531152d+00
COFD(2191) = -4.88789821d-01
COFD(2192) = 2.04437116d-02
COFD(2193) = -2.09376196d+01
COFD(2194) = 5.40870099d+00
COFD(2195) = -4.73017610d-01
COFD(2196) = 1.99399066d-02
COFD(2197) = -2.09376196d+01
COFD(2198) = 5.40870099d+00
COFD(2199) = -4.73017610d-01
COFD(2200) = 1.99399066d-02
COFD(2201) = -2.13538553d+01
COFD(2202) = 5.54007827d+00
COFD(2203) = -4.86434511d-01
COFD(2204) = 2.03779006d-02
COFD(2205) = -1.88171077d+01
COFD(2206) = 4.68393046d+00
COFD(2207) = -3.91610863d-01
COFD(2208) = 1.69174645d-02
COFD(2209) = -1.88241079d+01
COFD(2210) = 4.68393046d+00
COFD(2211) = -3.91610863d-01
COFD(2212) = 1.69174645d-02
COFD(2213) = -1.88306747d+01
COFD(2214) = 4.68393046d+00
COFD(2215) = -3.91610863d-01
COFD(2216) = 1.69174645d-02
COFD(2217) = -1.98418115d+01
COFD(2218) = 5.04367502d+00
COFD(2219) = -4.34153325d-01
COFD(2220) = 1.85956055d-02
COFD(2221) = -1.95263312d+01
COFD(2222) = 4.90255048d+00
COFD(2223) = -4.17368501d-01
COFD(2224) = 1.79287358d-02
COFD(2225) = -1.95318173d+01
COFD(2226) = 4.90255048d+00
COFD(2227) = -4.17368501d-01
COFD(2228) = 1.79287358d-02
COFD(2229) = -1.72572042d+01
COFD(2230) = 4.26063341d+00
COFD(2231) = -3.39848064d-01
COFD(2232) = 1.48021313d-02
COFD(2233) = -2.11349086d+01
COFD(2234) = 5.42846112d+00
COFD(2235) = -4.74321870d-01
COFD(2236) = 1.99459749d-02
COFD(2237) = -2.11349086d+01
COFD(2238) = 5.42846112d+00
COFD(2239) = -4.74321870d-01
COFD(2240) = 1.99459749d-02
COFD(2241) = -1.47958130d+01
COFD(2242) = 3.33113524d+00
COFD(2243) = -2.21479057d-01
COFD(2244) = 9.75837737d-03
COFD(2245) = -1.50125659d+01
COFD(2246) = 3.47945612d+00
COFD(2247) = -2.40703722d-01
COFD(2248) = 1.05907441d-02
COFD(2249) = -1.50279940d+01
COFD(2250) = 3.47945612d+00
COFD(2251) = -2.40703722d-01
COFD(2252) = 1.05907441d-02
COFD(2253) = -2.10302665d+01
COFD(2254) = 5.53859192d+00
COFD(2255) = -4.86285742d-01
COFD(2256) = 2.03731899d-02
COFD(2257) = -1.50817474d+01
COFD(2258) = 3.33113524d+00
COFD(2259) = -2.21479057d-01
COFD(2260) = 9.75837737d-03
COFD(2261) = -1.59516918d+01
COFD(2262) = 3.66853818d+00
COFD(2263) = -2.64346221d-01
COFD(2264) = 1.15784613d-02
COFD(2265) = -1.85866399d+01
COFD(2266) = 4.62613551d+00
COFD(2267) = -3.84556396d-01
COFD(2268) = 1.66292467d-02
COFD(2269) = -1.91704234d+01
COFD(2270) = 4.80030220d+00
COFD(2271) = -4.05235041d-01
COFD(2272) = 1.74483531d-02
COFD(2273) = -1.66814763d+01
COFD(2274) = 3.95813366d+00
COFD(2275) = -3.01970405d-01
COFD(2276) = 1.32132947d-02
COFD(2277) = -1.52696275d+01
COFD(2278) = 3.40487615d+00
COFD(2279) = -2.31175464d-01
COFD(2280) = 1.01841353d-02
COFD(2281) = -2.16351552d+01
COFD(2282) = 5.61201619d+00
COFD(2283) = -4.90763024d-01
COFD(2284) = 2.03690162d-02
COFD(2285) = -2.16413359d+01
COFD(2286) = 5.61201619d+00
COFD(2287) = -4.90763024d-01
COFD(2288) = 2.03690162d-02
COFD(2289) = -1.72572130d+01
COFD(2290) = 4.26063341d+00
COFD(2291) = -3.39848064d-01
COFD(2292) = 1.48021313d-02
COFD(2293) = -1.91676574d+01
COFD(2294) = 4.80030220d+00
COFD(2295) = -4.05235041d-01
COFD(2296) = 1.74483531d-02
COFD(2297) = -1.91676574d+01
COFD(2298) = 4.80030220d+00
COFD(2299) = -4.05235041d-01
COFD(2300) = 1.74483531d-02
COFD(2301) = -1.91676574d+01
COFD(2302) = 4.80030220d+00
COFD(2303) = -4.05235041d-01
COFD(2304) = 1.74483531d-02
COFD(2305) = -1.91647170d+01
COFD(2306) = 4.80030220d+00
COFD(2307) = -4.05235041d-01
COFD(2308) = 1.74483531d-02
COFD(2309) = -1.59404882d+01
COFD(2310) = 3.66853818d+00
COFD(2311) = -2.64346221d-01
COFD(2312) = 1.15784613d-02
COFD(2313) = -1.71942502d+01
COFD(2314) = 4.14993355d+00
COFD(2315) = -3.26168062d-01
COFD(2316) = 1.42364115d-02
COFD(2317) = -1.99607067d+01
COFD(2318) = 4.97875278d+00
COFD(2319) = -4.26485475d-01
COFD(2320) = 1.82931933d-02
COFD(2321) = -1.99635214d+01
COFD(2322) = 4.97875278d+00
COFD(2323) = -4.26485475d-01
COFD(2324) = 1.82931933d-02
COFD(2325) = -2.11378465d+01
COFD(2326) = 5.42846112d+00
COFD(2327) = -4.74321870d-01
COFD(2328) = 1.99459749d-02
COFD(2329) = -2.11406662d+01
COFD(2330) = 5.42846112d+00
COFD(2331) = -4.74321870d-01
COFD(2332) = 1.99459749d-02
COFD(2333) = -1.25098960d+01
COFD(2334) = 2.77873601d+00
COFD(2335) = -1.50637360d-01
COFD(2336) = 6.72684281d-03
COFD(2337) = -1.57972369d+01
COFD(2338) = 4.22225052d+00
COFD(2339) = -3.35156428d-01
COFD(2340) = 1.46104855d-02
COFD(2341) = -1.50584249d+01
COFD(2342) = 3.47945612d+00
COFD(2343) = -2.40703722d-01
COFD(2344) = 1.05907441d-02
COFD(2345) = -1.63254691d+01
COFD(2346) = 3.82388595d+00
COFD(2347) = -2.84480724d-01
COFD(2348) = 1.24506311d-02
COFD(2349) = -1.50724636d+01
COFD(2350) = 3.47945612d+00
COFD(2351) = -2.40703722d-01
COFD(2352) = 1.05907441d-02
COFD(2353) = -2.12639214d+01
COFD(2354) = 5.61184117d+00
COFD(2355) = -4.90532156d-01
COFD(2356) = 2.03507922d-02
COFD(2357) = -1.63301444d+01
COFD(2358) = 3.82388595d+00
COFD(2359) = -2.84480724d-01
COFD(2360) = 1.24506311d-02
COFD(2361) = -1.63345829d+01
COFD(2362) = 3.82388595d+00
COFD(2363) = -2.84480724d-01
COFD(2364) = 1.24506311d-02
COFD(2365) = -1.47558850d+01
COFD(2366) = 3.33113524d+00
COFD(2367) = -2.21479057d-01
COFD(2368) = 9.75837737d-03
COFD(2369) = -1.50076254d+01
COFD(2370) = 3.47945612d+00
COFD(2371) = -2.40703722d-01
COFD(2372) = 1.05907441d-02
COFD(2373) = -1.73027557d+01
COFD(2374) = 4.21416723d+00
COFD(2375) = -3.34163932d-01
COFD(2376) = 1.45697432d-02
COFD(2377) = -1.73027557d+01
COFD(2378) = 4.21416723d+00
COFD(2379) = -3.34163932d-01
COFD(2380) = 1.45697432d-02
COFD(2381) = -1.73198034d+01
COFD(2382) = 4.21416723d+00
COFD(2383) = -3.34163932d-01
COFD(2384) = 1.45697432d-02
COFD(2385) = -1.72556729d+01
COFD(2386) = 4.19029808d+00
COFD(2387) = -3.31177076d-01
COFD(2388) = 1.44446234d-02
COFD(2389) = -1.59634533d+01
COFD(2390) = 3.67388294d+00
COFD(2391) = -2.64990709d-01
COFD(2392) = 1.16042706d-02
COFD(2393) = -1.93015555d+01
COFD(2394) = 4.85015581d+00
COFD(2395) = -4.10945109d-01
COFD(2396) = 1.76651398d-02
COFD(2397) = -2.14160703d+01
COFD(2398) = 5.56531152d+00
COFD(2399) = -4.88789821d-01
COFD(2400) = 2.04437116d-02
COFD(2401) = -2.14215700d+01
COFD(2402) = 5.56531152d+00
COFD(2403) = -4.88789821d-01
COFD(2404) = 2.04437116d-02
COFD(2405) = -2.09376196d+01
COFD(2406) = 5.40870099d+00
COFD(2407) = -4.73017610d-01
COFD(2408) = 1.99399066d-02
COFD(2409) = -2.09376196d+01
COFD(2410) = 5.40870099d+00
COFD(2411) = -4.73017610d-01
COFD(2412) = 1.99399066d-02
COFD(2413) = -2.13538553d+01
COFD(2414) = 5.54007827d+00
COFD(2415) = -4.86434511d-01
COFD(2416) = 2.03779006d-02
COFD(2417) = -1.88171077d+01
COFD(2418) = 4.68393046d+00
COFD(2419) = -3.91610863d-01
COFD(2420) = 1.69174645d-02
COFD(2421) = -1.88241079d+01
COFD(2422) = 4.68393046d+00
COFD(2423) = -3.91610863d-01
COFD(2424) = 1.69174645d-02
COFD(2425) = -1.88306747d+01
COFD(2426) = 4.68393046d+00
COFD(2427) = -3.91610863d-01
COFD(2428) = 1.69174645d-02
COFD(2429) = -1.98418115d+01
COFD(2430) = 5.04367502d+00
COFD(2431) = -4.34153325d-01
COFD(2432) = 1.85956055d-02
COFD(2433) = -1.95263312d+01
COFD(2434) = 4.90255048d+00
COFD(2435) = -4.17368501d-01
COFD(2436) = 1.79287358d-02
COFD(2437) = -1.95318173d+01
COFD(2438) = 4.90255048d+00
COFD(2439) = -4.17368501d-01
COFD(2440) = 1.79287358d-02
COFD(2441) = -1.72572042d+01
COFD(2442) = 4.26063341d+00
COFD(2443) = -3.39848064d-01
COFD(2444) = 1.48021313d-02
COFD(2445) = -2.11349086d+01
COFD(2446) = 5.42846112d+00
COFD(2447) = -4.74321870d-01
COFD(2448) = 1.99459749d-02
COFD(2449) = -2.11349086d+01
COFD(2450) = 5.42846112d+00
COFD(2451) = -4.74321870d-01
COFD(2452) = 1.99459749d-02
COFD(2453) = -1.47958130d+01
COFD(2454) = 3.33113524d+00
COFD(2455) = -2.21479057d-01
COFD(2456) = 9.75837737d-03
COFD(2457) = -1.50125659d+01
COFD(2458) = 3.47945612d+00
COFD(2459) = -2.40703722d-01
COFD(2460) = 1.05907441d-02
COFD(2461) = -1.50279940d+01
COFD(2462) = 3.47945612d+00
COFD(2463) = -2.40703722d-01
COFD(2464) = 1.05907441d-02
COFD(2465) = -2.10302665d+01
COFD(2466) = 5.53859192d+00
COFD(2467) = -4.86285742d-01
COFD(2468) = 2.03731899d-02
COFD(2469) = -1.50817474d+01
COFD(2470) = 3.33113524d+00
COFD(2471) = -2.21479057d-01
COFD(2472) = 9.75837737d-03
COFD(2473) = -1.59516918d+01
COFD(2474) = 3.66853818d+00
COFD(2475) = -2.64346221d-01
COFD(2476) = 1.15784613d-02
COFD(2477) = -1.85866399d+01
COFD(2478) = 4.62613551d+00
COFD(2479) = -3.84556396d-01
COFD(2480) = 1.66292467d-02
COFD(2481) = -1.91704234d+01
COFD(2482) = 4.80030220d+00
COFD(2483) = -4.05235041d-01
COFD(2484) = 1.74483531d-02
COFD(2485) = -1.66814763d+01
COFD(2486) = 3.95813366d+00
COFD(2487) = -3.01970405d-01
COFD(2488) = 1.32132947d-02
COFD(2489) = -1.52696275d+01
COFD(2490) = 3.40487615d+00
COFD(2491) = -2.31175464d-01
COFD(2492) = 1.01841353d-02
COFD(2493) = -2.16351552d+01
COFD(2494) = 5.61201619d+00
COFD(2495) = -4.90763024d-01
COFD(2496) = 2.03690162d-02
COFD(2497) = -2.16413359d+01
COFD(2498) = 5.61201619d+00
COFD(2499) = -4.90763024d-01
COFD(2500) = 2.03690162d-02
COFD(2501) = -1.72572130d+01
COFD(2502) = 4.26063341d+00
COFD(2503) = -3.39848064d-01
COFD(2504) = 1.48021313d-02
COFD(2505) = -1.91676574d+01
COFD(2506) = 4.80030220d+00
COFD(2507) = -4.05235041d-01
COFD(2508) = 1.74483531d-02
COFD(2509) = -1.91676574d+01
COFD(2510) = 4.80030220d+00
COFD(2511) = -4.05235041d-01
COFD(2512) = 1.74483531d-02
COFD(2513) = -1.91676574d+01
COFD(2514) = 4.80030220d+00
COFD(2515) = -4.05235041d-01
COFD(2516) = 1.74483531d-02
COFD(2517) = -1.91647170d+01
COFD(2518) = 4.80030220d+00
COFD(2519) = -4.05235041d-01
COFD(2520) = 1.74483531d-02
COFD(2521) = -1.59404882d+01
COFD(2522) = 3.66853818d+00
COFD(2523) = -2.64346221d-01
COFD(2524) = 1.15784613d-02
COFD(2525) = -1.71942502d+01
COFD(2526) = 4.14993355d+00
COFD(2527) = -3.26168062d-01
COFD(2528) = 1.42364115d-02
COFD(2529) = -1.99607067d+01
COFD(2530) = 4.97875278d+00
COFD(2531) = -4.26485475d-01
COFD(2532) = 1.82931933d-02
COFD(2533) = -1.99635214d+01
COFD(2534) = 4.97875278d+00
COFD(2535) = -4.26485475d-01
COFD(2536) = 1.82931933d-02
COFD(2537) = -2.11378465d+01
COFD(2538) = 5.42846112d+00
COFD(2539) = -4.74321870d-01
COFD(2540) = 1.99459749d-02
COFD(2541) = -2.11406662d+01
COFD(2542) = 5.42846112d+00
COFD(2543) = -4.74321870d-01
COFD(2544) = 1.99459749d-02
COFD(2545) = -1.25141260d+01
COFD(2546) = 2.77873601d+00
COFD(2547) = -1.50637360d-01
COFD(2548) = 6.72684281d-03
COFD(2549) = -1.57994893d+01
COFD(2550) = 4.22225052d+00
COFD(2551) = -3.35156428d-01
COFD(2552) = 1.46104855d-02
COFD(2553) = -1.50766130d+01
COFD(2554) = 3.47945612d+00
COFD(2555) = -2.40703722d-01
COFD(2556) = 1.05907441d-02
COFD(2557) = -1.63493345d+01
COFD(2558) = 3.82388595d+00
COFD(2559) = -2.84480724d-01
COFD(2560) = 1.24506311d-02
COFD(2561) = -1.50911794d+01
COFD(2562) = 3.47945612d+00
COFD(2563) = -2.40703722d-01
COFD(2564) = 1.05907441d-02
COFD(2565) = -2.12831323d+01
COFD(2566) = 5.61184117d+00
COFD(2567) = -4.90532156d-01
COFD(2568) = 2.03507922d-02
COFD(2569) = -1.63542394d+01
COFD(2570) = 3.82388595d+00
COFD(2571) = -2.84480724d-01
COFD(2572) = 1.24506311d-02
COFD(2573) = -1.63588981d+01
COFD(2574) = 3.82388595d+00
COFD(2575) = -2.84480724d-01
COFD(2576) = 1.24506311d-02
COFD(2577) = -1.47715918d+01
COFD(2578) = 3.33113524d+00
COFD(2579) = -2.21479057d-01
COFD(2580) = 9.75837737d-03
COFD(2581) = -1.50240272d+01
COFD(2582) = 3.47945612d+00
COFD(2583) = -2.40703722d-01
COFD(2584) = 1.05907441d-02
COFD(2585) = -1.73198034d+01
COFD(2586) = 4.21416723d+00
COFD(2587) = -3.34163932d-01
COFD(2588) = 1.45697432d-02
COFD(2589) = -1.73198034d+01
COFD(2590) = 4.21416723d+00
COFD(2591) = -3.34163932d-01
COFD(2592) = 1.45697432d-02
COFD(2593) = -1.73374529d+01
COFD(2594) = 4.21416723d+00
COFD(2595) = -3.34163932d-01
COFD(2596) = 1.45697432d-02
COFD(2597) = -1.72738845d+01
COFD(2598) = 4.19029808d+00
COFD(2599) = -3.31177076d-01
COFD(2600) = 1.44446234d-02
COFD(2601) = -1.59863030d+01
COFD(2602) = 3.67388294d+00
COFD(2603) = -2.64990709d-01
COFD(2604) = 1.16042706d-02
COFD(2605) = -1.93276434d+01
COFD(2606) = 4.85015581d+00
COFD(2607) = -4.10945109d-01
COFD(2608) = 1.76651398d-02
COFD(2609) = -2.14391943d+01
COFD(2610) = 5.56531152d+00
COFD(2611) = -4.88789821d-01
COFD(2612) = 2.04437116d-02
COFD(2613) = -2.14449559d+01
COFD(2614) = 5.56531152d+00
COFD(2615) = -4.88789821d-01
COFD(2616) = 2.04437116d-02
COFD(2617) = -2.09612557d+01
COFD(2618) = 5.40870099d+00
COFD(2619) = -4.73017610d-01
COFD(2620) = 1.99399066d-02
COFD(2621) = -2.09612557d+01
COFD(2622) = 5.40870099d+00
COFD(2623) = -4.73017610d-01
COFD(2624) = 1.99399066d-02
COFD(2625) = -2.13777308d+01
COFD(2626) = 5.54007827d+00
COFD(2627) = -4.86434511d-01
COFD(2628) = 2.03779006d-02
COFD(2629) = -1.88390649d+01
COFD(2630) = 4.68393046d+00
COFD(2631) = -3.91610863d-01
COFD(2632) = 1.69174645d-02
COFD(2633) = -1.88463816d+01
COFD(2634) = 4.68393046d+00
COFD(2635) = -3.91610863d-01
COFD(2636) = 1.69174645d-02
COFD(2637) = -1.88532497d+01
COFD(2638) = 4.68393046d+00
COFD(2639) = -3.91610863d-01
COFD(2640) = 1.69174645d-02
COFD(2641) = -1.98646734d+01
COFD(2642) = 5.04367502d+00
COFD(2643) = -4.34153325d-01
COFD(2644) = 1.85956055d-02
COFD(2645) = -1.95494668d+01
COFD(2646) = 4.90255048d+00
COFD(2647) = -4.17368501d-01
COFD(2648) = 1.79287358d-02
COFD(2649) = -1.95552142d+01
COFD(2650) = 4.90255048d+00
COFD(2651) = -4.17368501d-01
COFD(2652) = 1.79287358d-02
COFD(2653) = -1.72828302d+01
COFD(2654) = 4.26063341d+00
COFD(2655) = -3.39848064d-01
COFD(2656) = 1.48021313d-02
COFD(2657) = -2.11606963d+01
COFD(2658) = 5.42846112d+00
COFD(2659) = -4.74321870d-01
COFD(2660) = 1.99459749d-02
COFD(2661) = -2.11606963d+01
COFD(2662) = 5.42846112d+00
COFD(2663) = -4.74321870d-01
COFD(2664) = 1.99459749d-02
COFD(2665) = -1.48128481d+01
COFD(2666) = 3.33113524d+00
COFD(2667) = -2.21479057d-01
COFD(2668) = 9.75837737d-03
COFD(2669) = -1.50302037d+01
COFD(2670) = 3.47945612d+00
COFD(2671) = -2.40703722d-01
COFD(2672) = 1.05907441d-02
COFD(2673) = -1.50461946d+01
COFD(2674) = 3.47945612d+00
COFD(2675) = -2.40703722d-01
COFD(2676) = 1.05907441d-02
COFD(2677) = -2.10489941d+01
COFD(2678) = 5.53859192d+00
COFD(2679) = -4.86285742d-01
COFD(2680) = 2.03731899d-02
COFD(2681) = -1.51048721d+01
COFD(2682) = 3.33113524d+00
COFD(2683) = -2.21479057d-01
COFD(2684) = 9.75837737d-03
COFD(2685) = -1.59750724d+01
COFD(2686) = 3.66853818d+00
COFD(2687) = -2.64346221d-01
COFD(2688) = 1.15784613d-02
COFD(2689) = -1.86130116d+01
COFD(2690) = 4.62613551d+00
COFD(2691) = -3.84556396d-01
COFD(2692) = 1.66292467d-02
COFD(2693) = -1.91965117d+01
COFD(2694) = 4.80030220d+00
COFD(2695) = -4.05235041d-01
COFD(2696) = 1.74483531d-02
COFD(2697) = -1.67051074d+01
COFD(2698) = 3.95813366d+00
COFD(2699) = -3.01970405d-01
COFD(2700) = 1.32132947d-02
COFD(2701) = -1.52918950d+01
COFD(2702) = 3.40487615d+00
COFD(2703) = -2.31175464d-01
COFD(2704) = 1.01841353d-02
COFD(2705) = -2.16577242d+01
COFD(2706) = 5.61201619d+00
COFD(2707) = -4.90763024d-01
COFD(2708) = 2.03690162d-02
COFD(2709) = -2.16641921d+01
COFD(2710) = 5.61201619d+00
COFD(2711) = -4.90763024d-01
COFD(2712) = 2.03690162d-02
COFD(2713) = -1.72828396d+01
COFD(2714) = 4.26063341d+00
COFD(2715) = -3.39848064d-01
COFD(2716) = 1.48021313d-02
COFD(2717) = -1.91935980d+01
COFD(2718) = 4.80030220d+00
COFD(2719) = -4.05235041d-01
COFD(2720) = 1.74483531d-02
COFD(2721) = -1.91935980d+01
COFD(2722) = 4.80030220d+00
COFD(2723) = -4.05235041d-01
COFD(2724) = 1.74483531d-02
COFD(2725) = -1.91935980d+01
COFD(2726) = 4.80030220d+00
COFD(2727) = -4.05235041d-01
COFD(2728) = 1.74483531d-02
COFD(2729) = -1.91905015d+01
COFD(2730) = 4.80030220d+00
COFD(2731) = -4.05235041d-01
COFD(2732) = 1.74483531d-02
COFD(2733) = -1.59633387d+01
COFD(2734) = 3.66853818d+00
COFD(2735) = -2.64346221d-01
COFD(2736) = 1.15784613d-02
COFD(2737) = -1.72196961d+01
COFD(2738) = 4.14993355d+00
COFD(2739) = -3.26168062d-01
COFD(2740) = 1.42364115d-02
COFD(2741) = -1.99866570d+01
COFD(2742) = 4.97875278d+00
COFD(2743) = -4.26485475d-01
COFD(2744) = 1.82931933d-02
COFD(2745) = -1.99896221d+01
COFD(2746) = 4.97875278d+00
COFD(2747) = -4.26485475d-01
COFD(2748) = 1.82931933d-02
COFD(2749) = -2.11637902d+01
COFD(2750) = 5.42846112d+00
COFD(2751) = -4.74321870d-01
COFD(2752) = 1.99459749d-02
COFD(2753) = -2.11667605d+01
COFD(2754) = 5.42846112d+00
COFD(2755) = -4.74321870d-01
COFD(2756) = 1.99459749d-02
COFD(2757) = -1.24693568d+01
COFD(2758) = 2.76686648d+00
COFD(2759) = -1.49120141d-01
COFD(2760) = 6.66220432d-03
COFD(2761) = -1.57199037d+01
COFD(2762) = 4.19936335d+00
COFD(2763) = -3.32311009d-01
COFD(2764) = 1.44921003d-02
COFD(2765) = -1.50270339d+01
COFD(2766) = 3.46140064d+00
COFD(2767) = -2.38440092d-01
COFD(2768) = 1.04960087d-02
COFD(2769) = -1.62724462d+01
COFD(2770) = 3.79163564d+00
COFD(2771) = -2.80257365d-01
COFD(2772) = 1.22656902d-02
COFD(2773) = -1.50420953d+01
COFD(2774) = 3.46140064d+00
COFD(2775) = -2.38440092d-01
COFD(2776) = 1.04960087d-02
COFD(2777) = -2.14087397d+01
COFD(2778) = 5.57282008d+00
COFD(2779) = -4.76690890d-01
COFD(2780) = 1.94000719d-02
COFD(2781) = -1.62775714d+01
COFD(2782) = 3.79163564d+00
COFD(2783) = -2.80257365d-01
COFD(2784) = 1.22656902d-02
COFD(2785) = -1.62824412d+01
COFD(2786) = 3.79163564d+00
COFD(2787) = -2.80257365d-01
COFD(2788) = 1.22656902d-02
COFD(2789) = -1.47048024d+01
COFD(2790) = 3.30594991d+00
COFD(2791) = -2.18182207d-01
COFD(2792) = 9.61429447d-03
COFD(2793) = -1.49727799d+01
COFD(2794) = 3.46140064d+00
COFD(2795) = -2.38440092d-01
COFD(2796) = 1.04960087d-02
COFD(2797) = -1.72556729d+01
COFD(2798) = 4.19029808d+00
COFD(2799) = -3.31177076d-01
COFD(2800) = 1.44446234d-02
COFD(2801) = -1.72556729d+01
COFD(2802) = 4.19029808d+00
COFD(2803) = -3.31177076d-01
COFD(2804) = 1.44446234d-02
COFD(2805) = -1.72738845d+01
COFD(2806) = 4.19029808d+00
COFD(2807) = -3.31177076d-01
COFD(2808) = 1.44446234d-02
COFD(2809) = -1.72167708d+01
COFD(2810) = 4.16886779d+00
COFD(2811) = -3.28518156d-01
COFD(2812) = 1.43341626d-02
COFD(2813) = -1.59525102d+01
COFD(2814) = 3.66023858d+00
COFD(2815) = -2.63401043d-01
COFD(2816) = 1.15432000d-02
COFD(2817) = -1.92867554d+01
COFD(2818) = 4.83375900d+00
COFD(2819) = -4.09146560d-01
COFD(2820) = 1.76006599d-02
COFD(2821) = -2.14022336d+01
COFD(2822) = 5.55346617d+00
COFD(2823) = -4.87783156d-01
COFD(2824) = 2.04210886d-02
COFD(2825) = -2.14082453d+01
COFD(2826) = 5.55346617d+00
COFD(2827) = -4.87783156d-01
COFD(2828) = 2.04210886d-02
COFD(2829) = -2.11381508d+01
COFD(2830) = 5.45574440d+00
COFD(2831) = -4.77436155d-01
COFD(2832) = 2.00644596d-02
COFD(2833) = -2.11381508d+01
COFD(2834) = 5.45574440d+00
COFD(2835) = -4.77436155d-01
COFD(2836) = 2.00644596d-02
COFD(2837) = -2.13319784d+01
COFD(2838) = 5.52422470d+00
COFD(2839) = -4.84872944d-01
COFD(2840) = 2.03298213d-02
COFD(2841) = -1.87821119d+01
COFD(2842) = 4.66162351d+00
COFD(2843) = -3.88920477d-01
COFD(2844) = 1.68089648d-02
COFD(2845) = -1.87897298d+01
COFD(2846) = 4.66162351d+00
COFD(2847) = -3.88920477d-01
COFD(2848) = 1.68089648d-02
COFD(2849) = -1.87968848d+01
COFD(2850) = 4.66162351d+00
COFD(2851) = -3.88920477d-01
COFD(2852) = 1.68089648d-02
COFD(2853) = -1.98075055d+01
COFD(2854) = 5.02169524d+00
COFD(2855) = -4.31582804d-01
COFD(2856) = 1.84953568d-02
COFD(2857) = -1.94763688d+01
COFD(2858) = 4.87333294d+00
COFD(2859) = -4.13769241d-01
COFD(2860) = 1.77802244d-02
COFD(2861) = -1.94823660d+01
COFD(2862) = 4.87333294d+00
COFD(2863) = -4.13769241d-01
COFD(2864) = 1.77802244d-02
COFD(2865) = -1.72316148d+01
COFD(2866) = 4.24011069d+00
COFD(2867) = -3.37339810d-01
COFD(2868) = 1.46996679d-02
COFD(2869) = -2.11309207d+01
COFD(2870) = 5.41773516d+00
COFD(2871) = -4.73414338d-01
COFD(2872) = 1.99258685d-02
COFD(2873) = -2.11309207d+01
COFD(2874) = 5.41773516d+00
COFD(2875) = -4.73414338d-01
COFD(2876) = 1.99258685d-02
COFD(2877) = -1.47472946d+01
COFD(2878) = 3.30594991d+00
COFD(2879) = -2.18182207d-01
COFD(2880) = 9.61429447d-03
COFD(2881) = -1.49798516d+01
COFD(2882) = 3.46140064d+00
COFD(2883) = -2.38440092d-01
COFD(2884) = 1.04960087d-02
COFD(2885) = -1.49963695d+01
COFD(2886) = 3.46140064d+00
COFD(2887) = -2.38440092d-01
COFD(2888) = 1.04960087d-02
COFD(2889) = -2.12890022d+01
COFD(2890) = 5.60080149d+00
COFD(2891) = -4.91551539d-01
COFD(2892) = 2.04914050d-02
COFD(2893) = -1.50460762d+01
COFD(2894) = 3.30594991d+00
COFD(2895) = -2.18182207d-01
COFD(2896) = 9.61429447d-03
COFD(2897) = -1.59449698d+01
COFD(2898) = 3.65620899d+00
COFD(2899) = -2.62933804d-01
COFD(2900) = 1.15253223d-02
COFD(2901) = -1.85418532d+01
COFD(2902) = 4.59643893d+00
COFD(2903) = -3.80823304d-01
COFD(2904) = 1.64720603d-02
COFD(2905) = -1.91494015d+01
COFD(2906) = 4.78094221d+00
COFD(2907) = -4.02985837d-01
COFD(2908) = 1.73614221d-02
COFD(2909) = -1.66377559d+01
COFD(2910) = 3.92993514d+00
COFD(2911) = -2.98305855d-01
COFD(2912) = 1.30538105d-02
COFD(2913) = -1.52271026d+01
COFD(2914) = 3.37751537d+00
COFD(2915) = -2.27579233d-01
COFD(2916) = 1.00262778d-02
COFD(2917) = -2.16493228d+01
COFD(2918) = 5.61229743d+00
COFD(2919) = -4.91428747d-01
COFD(2920) = 2.04226556d-02
COFD(2921) = -2.16560648d+01
COFD(2922) = 5.61229743d+00
COFD(2923) = -4.91428747d-01
COFD(2924) = 2.04226556d-02
COFD(2925) = -1.72316246d+01
COFD(2926) = 4.24011069d+00
COFD(2927) = -3.37339810d-01
COFD(2928) = 1.46996679d-02
COFD(2929) = -1.91463450d+01
COFD(2930) = 4.78094221d+00
COFD(2931) = -4.02985837d-01
COFD(2932) = 1.73614221d-02
COFD(2933) = -1.91463450d+01
COFD(2934) = 4.78094221d+00
COFD(2935) = -4.02985837d-01
COFD(2936) = 1.73614221d-02
COFD(2937) = -1.91463450d+01
COFD(2938) = 4.78094221d+00
COFD(2939) = -4.02985837d-01
COFD(2940) = 1.73614221d-02
COFD(2941) = -1.91430978d+01
COFD(2942) = 4.78094221d+00
COFD(2943) = -4.02985837d-01
COFD(2944) = 1.73614221d-02
COFD(2945) = -1.59327297d+01
COFD(2946) = 3.65620899d+00
COFD(2947) = -2.62933804d-01
COFD(2948) = 1.15253223d-02
COFD(2949) = -1.71754154d+01
COFD(2950) = 4.13131681d+00
COFD(2951) = -3.23897559d-01
COFD(2952) = 1.41438222d-02
COFD(2953) = -1.99301980d+01
COFD(2954) = 4.95514826d+00
COFD(2955) = -4.23691395d-01
COFD(2956) = 1.81828318d-02
COFD(2957) = -1.99333084d+01
COFD(2958) = 4.95514826d+00
COFD(2959) = -4.23691395d-01
COFD(2960) = 1.81828318d-02
COFD(2961) = -2.11341653d+01
COFD(2962) = 5.41773516d+00
COFD(2963) = -4.73414338d-01
COFD(2964) = 1.99258685d-02
COFD(2965) = -2.11372811d+01
COFD(2966) = 5.41773516d+00
COFD(2967) = -4.73414338d-01
COFD(2968) = 1.99258685d-02
COFD(2969) = -1.17159737d+01
COFD(2970) = 2.48123210d+00
COFD(2971) = -1.11322604d-01
COFD(2972) = 4.99282389d-03
COFD(2973) = -1.43151174d+01
COFD(2974) = 3.68038508d+00
COFD(2975) = -2.65779346d-01
COFD(2976) = 1.16360771d-02
COFD(2977) = -1.40999008d+01
COFD(2978) = 3.08120012d+00
COFD(2979) = -1.89629903d-01
COFD(2980) = 8.40361952d-03
COFD(2981) = -1.52721107d+01
COFD(2982) = 3.36790500d+00
COFD(2983) = -2.26321740d-01
COFD(2984) = 9.97135055d-03
COFD(2985) = -1.41191261d+01
COFD(2986) = 3.08120012d+00
COFD(2987) = -1.89629903d-01
COFD(2988) = 8.40361952d-03
COFD(2989) = -2.11388331d+01
COFD(2990) = 5.55529675d+00
COFD(2991) = -4.87942518d-01
COFD(2992) = 2.04249054d-02
COFD(2993) = -1.52792891d+01
COFD(2994) = 3.36790500d+00
COFD(2995) = -2.26321740d-01
COFD(2996) = 9.97135055d-03
COFD(2997) = -1.52861376d+01
COFD(2998) = 3.36790500d+00
COFD(2999) = -2.26321740d-01
COFD(3000) = 9.97135055d-03
COFD(3001) = -1.38862839d+01
COFD(3002) = 2.97564184d+00
COFD(3003) = -1.76025309d-01
COFD(3004) = 7.81869993d-03
COFD(3005) = -1.40318948d+01
COFD(3006) = 3.08120012d+00
COFD(3007) = -1.89629903d-01
COFD(3008) = 8.40361952d-03
COFD(3009) = -1.59634533d+01
COFD(3010) = 3.67388294d+00
COFD(3011) = -2.64990709d-01
COFD(3012) = 1.16042706d-02
COFD(3013) = -1.59634533d+01
COFD(3014) = 3.67388294d+00
COFD(3015) = -2.64990709d-01
COFD(3016) = 1.16042706d-02
COFD(3017) = -1.59863030d+01
COFD(3018) = 3.67388294d+00
COFD(3019) = -2.64990709d-01
COFD(3020) = 1.16042706d-02
COFD(3021) = -1.59525102d+01
COFD(3022) = 3.66023858d+00
COFD(3023) = -2.63401043d-01
COFD(3024) = 1.15432000d-02
COFD(3025) = -1.50233475d+01
COFD(3026) = 3.26660767d+00
COFD(3027) = -2.13287177d-01
COFD(3028) = 9.41137857d-03
COFD(3029) = -1.81735763d+01
COFD(3030) = 4.38391495d+00
COFD(3031) = -3.54941287d-01
COFD(3032) = 1.54195107d-02
COFD(3033) = -2.05045578d+01
COFD(3034) = 5.23843909d+00
COFD(3035) = -4.55815614d-01
COFD(3036) = 1.93898040d-02
COFD(3037) = -2.05128705d+01
COFD(3038) = 5.23843909d+00
COFD(3039) = -4.55815614d-01
COFD(3040) = 1.93898040d-02
COFD(3041) = -2.02642227d+01
COFD(3042) = 5.14499740d+00
COFD(3043) = -4.45694430d-01
COFD(3044) = 1.90318646d-02
COFD(3045) = -2.02642227d+01
COFD(3046) = 5.14499740d+00
COFD(3047) = -4.45694430d-01
COFD(3048) = 1.90318646d-02
COFD(3049) = -2.04144604d+01
COFD(3050) = 5.19614628d+00
COFD(3051) = -4.50889164d-01
COFD(3052) = 1.91983328d-02
COFD(3053) = -1.76182365d+01
COFD(3054) = 4.19935698d+00
COFD(3055) = -3.32310212d-01
COFD(3056) = 1.44920670d-02
COFD(3057) = -1.76285640d+01
COFD(3058) = 4.19935698d+00
COFD(3059) = -3.32310212d-01
COFD(3060) = 1.44920670d-02
COFD(3061) = -1.76383156d+01
COFD(3062) = 4.19935698d+00
COFD(3063) = -3.32310212d-01
COFD(3064) = 1.44920670d-02
COFD(3065) = -1.86157761d+01
COFD(3066) = 4.55689508d+00
COFD(3067) = -3.75937921d-01
COFD(3068) = 1.62703488d-02
COFD(3069) = -1.83455435d+01
COFD(3070) = 4.42828044d+00
COFD(3071) = -3.60417833d-01
COFD(3072) = 1.56455103d-02
COFD(3073) = -1.83538377d+01
COFD(3074) = 4.42828044d+00
COFD(3075) = -3.60417833d-01
COFD(3076) = 1.56455103d-02
COFD(3077) = -1.60261675d+01
COFD(3078) = 3.73312045d+00
COFD(3079) = -2.72579779d-01
COFD(3080) = 1.19290272d-02
COFD(3081) = -2.02922701d+01
COFD(3082) = 5.11106992d+00
COFD(3083) = -4.42047129d-01
COFD(3084) = 1.89042990d-02
COFD(3085) = -2.02922701d+01
COFD(3086) = 5.11106992d+00
COFD(3087) = -4.42047129d-01
COFD(3088) = 1.89042990d-02
COFD(3089) = -1.39388049d+01
COFD(3090) = 2.97564184d+00
COFD(3091) = -1.76025309d-01
COFD(3092) = 7.81869993d-03
COFD(3093) = -1.40479570d+01
COFD(3094) = 3.08120012d+00
COFD(3095) = -1.89629903d-01
COFD(3096) = 8.40361952d-03
COFD(3097) = -1.40688658d+01
COFD(3098) = 3.08120012d+00
COFD(3099) = -1.89629903d-01
COFD(3100) = 8.40361952d-03
COFD(3101) = -2.04839638d+01
COFD(3102) = 5.34597277d+00
COFD(3103) = -4.68432639d-01
COFD(3104) = 1.98846555d-02
COFD(3105) = -1.42892712d+01
COFD(3106) = 2.97564184d+00
COFD(3107) = -1.76025309d-01
COFD(3108) = 7.81869993d-03
COFD(3109) = -1.50200522d+01
COFD(3110) = 3.26223357d+00
COFD(3111) = -2.12746642d-01
COFD(3112) = 9.38912883d-03
COFD(3113) = -1.74527876d+01
COFD(3114) = 4.14792932d+00
COFD(3115) = -3.25920382d-01
COFD(3116) = 1.42261620d-02
COFD(3117) = -1.80065982d+01
COFD(3118) = 4.31656593d+00
COFD(3119) = -3.46539554d-01
COFD(3120) = 1.50688196d-02
COFD(3121) = -1.55536213d+01
COFD(3122) = 3.47339631d+00
COFD(3123) = -2.39946241d-01
COFD(3124) = 1.05591452d-02
COFD(3125) = -1.43882477d+01
COFD(3126) = 3.01675183d+00
COFD(3127) = -1.81281138d-01
COFD(3128) = 8.04273850d-03
COFD(3129) = -2.09183945d+01
COFD(3130) = 5.37381486d+00
COFD(3131) = -4.70526046d-01
COFD(3132) = 1.99137630d-02
COFD(3133) = -2.09276290d+01
COFD(3134) = 5.37381486d+00
COFD(3135) = -4.70526046d-01
COFD(3136) = 1.99137630d-02
COFD(3137) = -1.60261816d+01
COFD(3138) = 3.73312045d+00
COFD(3139) = -2.72579779d-01
COFD(3140) = 1.19290272d-02
COFD(3141) = -1.80021546d+01
COFD(3142) = 4.31656593d+00
COFD(3143) = -3.46539554d-01
COFD(3144) = 1.50688196d-02
COFD(3145) = -1.80021546d+01
COFD(3146) = 4.31656593d+00
COFD(3147) = -3.46539554d-01
COFD(3148) = 1.50688196d-02
COFD(3149) = -1.80021546d+01
COFD(3150) = 4.31656593d+00
COFD(3151) = -3.46539554d-01
COFD(3152) = 1.50688196d-02
COFD(3153) = -1.79974471d+01
COFD(3154) = 4.31656593d+00
COFD(3155) = -3.46539554d-01
COFD(3156) = 1.50688196d-02
COFD(3157) = -1.50031687d+01
COFD(3158) = 3.26223357d+00
COFD(3159) = -2.12746642d-01
COFD(3160) = 9.38912883d-03
COFD(3161) = -1.60074211d+01
COFD(3162) = 3.63723937d+00
COFD(3163) = -2.60754222d-01
COFD(3164) = 1.14428814d-02
COFD(3165) = -1.87780799d+01
COFD(3166) = 4.49191492d+00
COFD(3167) = -3.68041771d-01
COFD(3168) = 1.59498676d-02
COFD(3169) = -1.87826029d+01
COFD(3170) = 4.49191492d+00
COFD(3171) = -3.68041771d-01
COFD(3172) = 1.59498676d-02
COFD(3173) = -2.02969740d+01
COFD(3174) = 5.11106992d+00
COFD(3175) = -4.42047129d-01
COFD(3176) = 1.89042990d-02
COFD(3177) = -2.03015042d+01
COFD(3178) = 5.11106992d+00
COFD(3179) = -4.42047129d-01
COFD(3180) = 1.89042990d-02
COFD(3181) = -1.37794315d+01
COFD(3182) = 3.23973858d+00
COFD(3183) = -2.09989036d-01
COFD(3184) = 9.27667906d-03
COFD(3185) = -1.76147026d+01
COFD(3186) = 4.86049500d+00
COFD(3187) = -4.12200578d-01
COFD(3188) = 1.77160971d-02
COFD(3189) = -1.70534856d+01
COFD(3190) = 4.14240922d+00
COFD(3191) = -3.25239774d-01
COFD(3192) = 1.41980687d-02
COFD(3193) = -1.84688406d+01
COFD(3194) = 4.49330851d+00
COFD(3195) = -3.68208715d-01
COFD(3196) = 1.59565402d-02
COFD(3197) = -1.70757047d+01
COFD(3198) = 4.14240922d+00
COFD(3199) = -3.25239774d-01
COFD(3200) = 1.41980687d-02
COFD(3201) = -2.07653719d+01
COFD(3202) = 5.01092022d+00
COFD(3203) = -3.77985635d-01
COFD(3204) = 1.40968645d-02
COFD(3205) = -1.84777607d+01
COFD(3206) = 4.49330851d+00
COFD(3207) = -3.68208715d-01
COFD(3208) = 1.59565402d-02
COFD(3209) = -1.84863000d+01
COFD(3210) = 4.49330851d+00
COFD(3211) = -3.68208715d-01
COFD(3212) = 1.59565402d-02
COFD(3213) = -1.67390425d+01
COFD(3214) = 4.00828594d+00
COFD(3215) = -3.08414344d-01
COFD(3216) = 1.34907430d-02
COFD(3217) = -1.69758891d+01
COFD(3218) = 4.14240922d+00
COFD(3219) = -3.25239774d-01
COFD(3220) = 1.41980687d-02
COFD(3221) = -1.93015555d+01
COFD(3222) = 4.85015581d+00
COFD(3223) = -4.10945109d-01
COFD(3224) = 1.76651398d-02
COFD(3225) = -1.93015555d+01
COFD(3226) = 4.85015581d+00
COFD(3227) = -4.10945109d-01
COFD(3228) = 1.76651398d-02
COFD(3229) = -1.93276434d+01
COFD(3230) = 4.85015581d+00
COFD(3231) = -4.10945109d-01
COFD(3232) = 1.76651398d-02
COFD(3233) = -1.92867554d+01
COFD(3234) = 4.83375900d+00
COFD(3235) = -4.09146560d-01
COFD(3236) = 1.76006599d-02
COFD(3237) = -1.81735763d+01
COFD(3238) = 4.38391495d+00
COFD(3239) = -3.54941287d-01
COFD(3240) = 1.54195107d-02
COFD(3241) = -2.13425698d+01
COFD(3242) = 5.40460130d+00
COFD(3243) = -4.72718910d-01
COFD(3244) = 1.99362717d-02
COFD(3245) = -2.19215555d+01
COFD(3246) = 5.45216133d+00
COFD(3247) = -4.52916925d-01
COFD(3248) = 1.80456400d-02
COFD(3249) = -2.19317743d+01
COFD(3250) = 5.45216133d+00
COFD(3251) = -4.52916925d-01
COFD(3252) = 1.80456400d-02
COFD(3253) = -2.20421041d+01
COFD(3254) = 5.52708332d+00
COFD(3255) = -4.68000808d-01
COFD(3256) = 1.89131908d-02
COFD(3257) = -2.20421041d+01
COFD(3258) = 5.52708332d+00
COFD(3259) = -4.68000808d-01
COFD(3260) = 1.89131908d-02
COFD(3261) = -2.20063594d+01
COFD(3262) = 5.48540187d+00
COFD(3263) = -4.58962148d-01
COFD(3264) = 1.83770355d-02
COFD(3265) = -2.09066354d+01
COFD(3266) = 5.30153901d+00
COFD(3267) = -4.63335119d-01
COFD(3268) = 1.96897053d-02
COFD(3269) = -2.09191285d+01
COFD(3270) = 5.30153901d+00
COFD(3271) = -4.63335119d-01
COFD(3272) = 1.96897053d-02
COFD(3273) = -2.09309753d+01
COFD(3274) = 5.30153901d+00
COFD(3275) = -4.63335119d-01
COFD(3276) = 1.96897053d-02
COFD(3277) = -2.16802612d+01
COFD(3278) = 5.52918296d+00
COFD(3279) = -4.85360709d-01
COFD(3280) = 2.03448006d-02
COFD(3281) = -2.14224484d+01
COFD(3282) = 5.41729961d+00
COFD(3283) = -4.73400281d-01
COFD(3284) = 1.99269567d-02
COFD(3285) = -2.14326461d+01
COFD(3286) = 5.41729961d+00
COFD(3287) = -4.73400281d-01
COFD(3288) = 1.99269567d-02
COFD(3289) = -1.94485982d+01
COFD(3290) = 4.91446566d+00
COFD(3291) = -4.18837152d-01
COFD(3292) = 1.79893537d-02
COFD(3293) = -2.22116706d+01
COFD(3294) = 5.54251230d+00
COFD(3295) = -4.70946314d-01
COFD(3296) = 1.90785869d-02
COFD(3297) = -2.22116706d+01
COFD(3298) = 5.54251230d+00
COFD(3299) = -4.70946314d-01
COFD(3300) = 1.90785869d-02
COFD(3301) = -1.67983919d+01
COFD(3302) = 4.00828594d+00
COFD(3303) = -3.08414344d-01
COFD(3304) = 1.34907430d-02
COFD(3305) = -1.69990507d+01
COFD(3306) = 4.14240922d+00
COFD(3307) = -3.25239774d-01
COFD(3308) = 1.41980687d-02
COFD(3309) = -1.70230718d+01
COFD(3310) = 4.14240922d+00
COFD(3311) = -3.25239774d-01
COFD(3312) = 1.41980687d-02
COFD(3313) = -2.14900480d+01
COFD(3314) = 5.39920344d+00
COFD(3315) = -4.43232069d-01
COFD(3316) = 1.75138765d-02
COFD(3317) = -1.71843946d+01
COFD(3318) = 4.00828594d+00
COFD(3319) = -3.08414344d-01
COFD(3320) = 1.34907430d-02
COFD(3321) = -1.81639592d+01
COFD(3322) = 4.37565431d+00
COFD(3323) = -3.53906025d-01
COFD(3324) = 1.53760786d-02
COFD(3325) = -2.07349645d+01
COFD(3326) = 5.23705932d+00
COFD(3327) = -4.55655792d-01
COFD(3328) = 1.93836339d-02
COFD(3329) = -2.12671597d+01
COFD(3330) = 5.38135645d+00
COFD(3331) = -4.71058360d-01
COFD(3332) = 1.99188046d-02
COFD(3333) = -1.87873369d+01
COFD(3334) = 4.60768689d+00
COFD(3335) = -3.82235384d-01
COFD(3336) = 1.65314065d-02
COFD(3337) = -1.73328509d+01
COFD(3338) = 4.06993115d+00
COFD(3339) = -3.16228805d-01
COFD(3340) = 1.38227651d-02
COFD(3341) = -2.16858789d+01
COFD(3342) = 5.30214979d+00
COFD(3343) = -4.26479693d-01
COFD(3344) = 1.66224918d-02
COFD(3345) = -2.16971429d+01
COFD(3346) = 5.30214979d+00
COFD(3347) = -4.26479693d-01
COFD(3348) = 1.66224918d-02
COFD(3349) = -1.94486162d+01
COFD(3350) = 4.91446566d+00
COFD(3351) = -4.18837152d-01
COFD(3352) = 1.79893537d-02
COFD(3353) = -2.12614542d+01
COFD(3354) = 5.38135645d+00
COFD(3355) = -4.71058360d-01
COFD(3356) = 1.99188046d-02
COFD(3357) = -2.12614542d+01
COFD(3358) = 5.38135645d+00
COFD(3359) = -4.71058360d-01
COFD(3360) = 1.99188046d-02
COFD(3361) = -2.12614542d+01
COFD(3362) = 5.38135645d+00
COFD(3363) = -4.71058360d-01
COFD(3364) = 1.99188046d-02
COFD(3365) = -2.12554255d+01
COFD(3366) = 5.38135645d+00
COFD(3367) = -4.71058360d-01
COFD(3368) = 1.99188046d-02
COFD(3369) = -1.81432461d+01
COFD(3370) = 4.37565431d+00
COFD(3371) = -3.53906025d-01
COFD(3372) = 1.53760786d-02
COFD(3373) = -1.93483692d+01
COFD(3374) = 4.79506290d+00
COFD(3375) = -4.04621659d-01
COFD(3376) = 1.74244230d-02
COFD(3377) = -2.18713229d+01
COFD(3378) = 5.47368915d+00
COFD(3379) = -4.79424291d-01
COFD(3380) = 2.01372920d-02
COFD(3381) = -2.18771314d+01
COFD(3382) = 5.47368915d+00
COFD(3383) = -4.79424291d-01
COFD(3384) = 2.01372920d-02
COFD(3385) = -2.22176950d+01
COFD(3386) = 5.54251230d+00
COFD(3387) = -4.70946314d-01
COFD(3388) = 1.90785869d-02
COFD(3389) = -2.22235122d+01
COFD(3390) = 5.54251230d+00
COFD(3391) = -4.70946314d-01
COFD(3392) = 1.90785869d-02
COFD(3393) = -1.60517370d+01
COFD(3394) = 4.11188603d+00
COFD(3395) = -3.21540884d-01
COFD(3396) = 1.40482564d-02
COFD(3397) = -1.97544450d+01
COFD(3398) = 5.56931926d+00
COFD(3399) = -4.89105511d-01
COFD(3400) = 2.04493129d-02
COFD(3401) = -1.94313116d+01
COFD(3402) = 5.02567894d+00
COFD(3403) = -4.32045169d-01
COFD(3404) = 1.85132214d-02
COFD(3405) = -2.08204449d+01
COFD(3406) = 5.35267674d+00
COFD(3407) = -4.69010505d-01
COFD(3408) = 1.98979152d-02
COFD(3409) = -1.94507876d+01
COFD(3410) = 5.02567894d+00
COFD(3411) = -4.32045169d-01
COFD(3412) = 1.85132214d-02
COFD(3413) = -1.77498543d+01
COFD(3414) = 3.57475686d+00
COFD(3415) = -1.56396297d-01
COFD(3416) = 3.12157721d-03
COFD(3417) = -2.08277598d+01
COFD(3418) = 5.35267674d+00
COFD(3419) = -4.69010505d-01
COFD(3420) = 1.98979152d-02
COFD(3421) = -2.08347403d+01
COFD(3422) = 5.35267674d+00
COFD(3423) = -4.69010505d-01
COFD(3424) = 1.98979152d-02
COFD(3425) = -1.90526941d+01
COFD(3426) = 4.86821670d+00
COFD(3427) = -4.13144121d-01
COFD(3428) = 1.77546701d-02
COFD(3429) = -1.93624931d+01
COFD(3430) = 5.02567894d+00
COFD(3431) = -4.32045169d-01
COFD(3432) = 1.85132214d-02
COFD(3433) = -2.14160703d+01
COFD(3434) = 5.56531152d+00
COFD(3435) = -4.88789821d-01
COFD(3436) = 2.04437116d-02
COFD(3437) = -2.14160703d+01
COFD(3438) = 5.56531152d+00
COFD(3439) = -4.88789821d-01
COFD(3440) = 2.04437116d-02
COFD(3441) = -2.14391943d+01
COFD(3442) = 5.56531152d+00
COFD(3443) = -4.88789821d-01
COFD(3444) = 2.04437116d-02
COFD(3445) = -2.14022336d+01
COFD(3446) = 5.55346617d+00
COFD(3447) = -4.87783156d-01
COFD(3448) = 2.04210886d-02
COFD(3449) = -2.05045578d+01
COFD(3450) = 5.23843909d+00
COFD(3451) = -4.55815614d-01
COFD(3452) = 1.93898040d-02
COFD(3453) = -2.19215555d+01
COFD(3454) = 5.45216133d+00
COFD(3455) = -4.52916925d-01
COFD(3456) = 1.80456400d-02
COFD(3457) = -1.90328712d+01
COFD(3458) = 3.99221757d+00
COFD(3459) = -2.19854880d-01
COFD(3460) = 6.22736279d-03
COFD(3461) = -1.90413348d+01
COFD(3462) = 3.99221757d+00
COFD(3463) = -2.19854880d-01
COFD(3464) = 6.22736279d-03
COFD(3465) = -2.01801667d+01
COFD(3466) = 4.53183330d+00
COFD(3467) = -3.02186760d-01
COFD(3468) = 1.02756490d-02
COFD(3469) = -2.01801667d+01
COFD(3470) = 4.53183330d+00
COFD(3471) = -3.02186760d-01
COFD(3472) = 1.02756490d-02
COFD(3473) = -1.93125662d+01
COFD(3474) = 4.10954793d+00
COFD(3475) = -2.37523329d-01
COFD(3476) = 7.08858141d-03
COFD(3477) = -2.19851338d+01
COFD(3478) = 5.55935694d+00
COFD(3479) = -4.74154740d-01
COFD(3480) = 1.92584304d-02
COFD(3481) = -2.19956352d+01
COFD(3482) = 5.55935694d+00
COFD(3483) = -4.74154740d-01
COFD(3484) = 1.92584304d-02
COFD(3485) = -2.20055544d+01
COFD(3486) = 5.55935694d+00
COFD(3487) = -4.74154740d-01
COFD(3488) = 1.92584304d-02
COFD(3489) = -2.16296373d+01
COFD(3490) = 5.29019717d+00
COFD(3491) = -4.24502606d-01
COFD(3492) = 1.65197343d-02
COFD(3493) = -2.19229190d+01
COFD(3494) = 5.41841631d+00
COFD(3495) = -4.46818971d-01
COFD(3496) = 1.77127652d-02
COFD(3497) = -2.19313638d+01
COFD(3498) = 5.41841631d+00
COFD(3499) = -4.46818971d-01
COFD(3500) = 1.77127652d-02
COFD(3501) = -2.14204185d+01
COFD(3502) = 5.59268435d+00
COFD(3503) = -4.91232974d-01
COFD(3504) = 2.05064746d-02
COFD(3505) = -2.00915040d+01
COFD(3506) = 4.41511629d+00
COFD(3507) = -2.84086963d-01
COFD(3508) = 9.37586971d-03
COFD(3509) = -2.00915040d+01
COFD(3510) = 4.41511629d+00
COFD(3511) = -2.84086963d-01
COFD(3512) = 9.37586971d-03
COFD(3513) = -1.91057988d+01
COFD(3514) = 4.86821670d+00
COFD(3515) = -4.13144121d-01
COFD(3516) = 1.77546701d-02
COFD(3517) = -1.93788111d+01
COFD(3518) = 5.02567894d+00
COFD(3519) = -4.32045169d-01
COFD(3520) = 1.85132214d-02
COFD(3521) = -1.93999821d+01
COFD(3522) = 5.02567894d+00
COFD(3523) = -4.32045169d-01
COFD(3524) = 1.85132214d-02
COFD(3525) = -1.89427986d+01
COFD(3526) = 4.11490576d+00
COFD(3527) = -2.38333952d-01
COFD(3528) = 7.12818549d-03
COFD(3529) = -1.94605277d+01
COFD(3530) = 4.86821670d+00
COFD(3531) = -4.13144121d-01
COFD(3532) = 1.77546701d-02
COFD(3533) = -2.04922452d+01
COFD(3534) = 5.23112374d+00
COFD(3535) = -4.54967682d-01
COFD(3536) = 1.93570423d-02
COFD(3537) = -2.19815149d+01
COFD(3538) = 5.58446511d+00
COFD(3539) = -4.79399331d-01
COFD(3540) = 1.95652693d-02
COFD(3541) = -2.19960393d+01
COFD(3542) = 5.49642957d+00
COFD(3543) = -4.61132993d-01
COFD(3544) = 1.85004773d-02
COFD(3545) = -2.09853721d+01
COFD(3546) = 5.39401640d+00
COFD(3547) = -4.72026580d-01
COFD(3548) = 1.99336592d-02
COFD(3549) = -1.96553384d+01
COFD(3550) = 4.94225344d+00
COFD(3551) = -4.22163426d-01
COFD(3552) = 1.81223730d-02
COFD(3553) = -1.82195639d+01
COFD(3554) = 3.59760627d+00
COFD(3555) = -1.59809307d-01
COFD(3556) = 3.28676045d-03
COFD(3557) = -1.82289601d+01
COFD(3558) = 3.59760627d+00
COFD(3559) = -1.59809307d-01
COFD(3560) = 3.28676045d-03
COFD(3561) = -2.14204329d+01
COFD(3562) = 5.59268435d+00
COFD(3563) = -4.91232974d-01
COFD(3564) = 2.05064746d-02
COFD(3565) = -2.19914997d+01
COFD(3566) = 5.49642957d+00
COFD(3567) = -4.61132993d-01
COFD(3568) = 1.85004773d-02
COFD(3569) = -2.19914997d+01
COFD(3570) = 5.49642957d+00
COFD(3571) = -4.61132993d-01
COFD(3572) = 1.85004773d-02
COFD(3573) = -2.19914997d+01
COFD(3574) = 5.49642957d+00
COFD(3575) = -4.61132993d-01
COFD(3576) = 1.85004773d-02
COFD(3577) = -2.19866916d+01
COFD(3578) = 5.49642957d+00
COFD(3579) = -4.61132993d-01
COFD(3580) = 1.85004773d-02
COFD(3581) = -2.04750581d+01
COFD(3582) = 5.23112374d+00
COFD(3583) = -4.54967682d-01
COFD(3584) = 1.93570423d-02
COFD(3585) = -2.14255087d+01
COFD(3586) = 5.52240865d+00
COFD(3587) = -4.84699537d-01
COFD(3588) = 2.03247833d-02
COFD(3589) = -2.20946362d+01
COFD(3590) = 5.36053938d+00
COFD(3591) = -4.36434519d-01
COFD(3592) = 1.71484255d-02
COFD(3593) = -2.20992569d+01
COFD(3594) = 5.36053938d+00
COFD(3595) = -4.36434519d-01
COFD(3596) = 1.71484255d-02
COFD(3597) = -2.00963085d+01
COFD(3598) = 4.41511629d+00
COFD(3599) = -2.84086963d-01
COFD(3600) = 9.37586971d-03
COFD(3601) = -2.01009366d+01
COFD(3602) = 4.41511629d+00
COFD(3603) = -2.84086963d-01
COFD(3604) = 9.37586971d-03
COFD(3605) = -1.60528285d+01
COFD(3606) = 4.11188603d+00
COFD(3607) = -3.21540884d-01
COFD(3608) = 1.40482564d-02
COFD(3609) = -1.97550088d+01
COFD(3610) = 5.56931926d+00
COFD(3611) = -4.89105511d-01
COFD(3612) = 2.04493129d-02
COFD(3613) = -1.94373127d+01
COFD(3614) = 5.02567894d+00
COFD(3615) = -4.32045169d-01
COFD(3616) = 1.85132214d-02
COFD(3617) = -2.08293255d+01
COFD(3618) = 5.35267674d+00
COFD(3619) = -4.69010505d-01
COFD(3620) = 1.98979152d-02
COFD(3621) = -1.94570287d+01
COFD(3622) = 5.02567894d+00
COFD(3623) = -4.32045169d-01
COFD(3624) = 1.85132214d-02
COFD(3625) = -1.77563250d+01
COFD(3626) = 3.57475686d+00
COFD(3627) = -1.56396297d-01
COFD(3628) = 3.12157721d-03
COFD(3629) = -2.08367725d+01
COFD(3630) = 5.35267674d+00
COFD(3631) = -4.69010505d-01
COFD(3632) = 1.98979152d-02
COFD(3633) = -2.08438809d+01
COFD(3634) = 5.35267674d+00
COFD(3635) = -4.69010505d-01
COFD(3636) = 1.98979152d-02
COFD(3637) = -1.90576320d+01
COFD(3638) = 4.86821670d+00
COFD(3639) = -4.13144121d-01
COFD(3640) = 1.77546701d-02
COFD(3641) = -1.93677186d+01
COFD(3642) = 5.02567894d+00
COFD(3643) = -4.32045169d-01
COFD(3644) = 1.85132214d-02
COFD(3645) = -2.14215700d+01
COFD(3646) = 5.56531152d+00
COFD(3647) = -4.88789821d-01
COFD(3648) = 2.04437116d-02
COFD(3649) = -2.14215700d+01
COFD(3650) = 5.56531152d+00
COFD(3651) = -4.88789821d-01
COFD(3652) = 2.04437116d-02
COFD(3653) = -2.14449559d+01
COFD(3654) = 5.56531152d+00
COFD(3655) = -4.88789821d-01
COFD(3656) = 2.04437116d-02
COFD(3657) = -2.14082453d+01
COFD(3658) = 5.55346617d+00
COFD(3659) = -4.87783156d-01
COFD(3660) = 2.04210886d-02
COFD(3661) = -2.05128705d+01
COFD(3662) = 5.23843909d+00
COFD(3663) = -4.55815614d-01
COFD(3664) = 1.93898040d-02
COFD(3665) = -2.19317743d+01
COFD(3666) = 5.45216133d+00
COFD(3667) = -4.52916925d-01
COFD(3668) = 1.80456400d-02
COFD(3669) = -1.90413348d+01
COFD(3670) = 3.99221757d+00
COFD(3671) = -2.19854880d-01
COFD(3672) = 6.22736279d-03
COFD(3673) = -1.90499441d+01
COFD(3674) = 3.99221757d+00
COFD(3675) = -2.19854880d-01
COFD(3676) = 6.22736279d-03
COFD(3677) = -2.01889168d+01
COFD(3678) = 4.53183330d+00
COFD(3679) = -3.02186760d-01
COFD(3680) = 1.02756490d-02
COFD(3681) = -2.01889168d+01
COFD(3682) = 4.53183330d+00
COFD(3683) = -3.02186760d-01
COFD(3684) = 1.02756490d-02
COFD(3685) = -1.93214527d+01
COFD(3686) = 4.10954793d+00
COFD(3687) = -2.37523329d-01
COFD(3688) = 7.08858141d-03
COFD(3689) = -2.19929679d+01
COFD(3690) = 5.55935694d+00
COFD(3691) = -4.74154740d-01
COFD(3692) = 1.92584304d-02
COFD(3693) = -2.20036369d+01
COFD(3694) = 5.55935694d+00
COFD(3695) = -4.74154740d-01
COFD(3696) = 1.92584304d-02
COFD(3697) = -2.20137178d+01
COFD(3698) = 5.55935694d+00
COFD(3699) = -4.74154740d-01
COFD(3700) = 1.92584304d-02
COFD(3701) = -2.16379567d+01
COFD(3702) = 5.29019717d+00
COFD(3703) = -4.24502606d-01
COFD(3704) = 1.65197343d-02
COFD(3705) = -2.19313890d+01
COFD(3706) = 5.41841631d+00
COFD(3707) = -4.46818971d-01
COFD(3708) = 1.77127652d-02
COFD(3709) = -2.19399793d+01
COFD(3710) = 5.41841631d+00
COFD(3711) = -4.46818971d-01
COFD(3712) = 1.77127652d-02
COFD(3713) = -2.14303479d+01
COFD(3714) = 5.59268435d+00
COFD(3715) = -4.91232974d-01
COFD(3716) = 2.05064746d-02
COFD(3717) = -2.01015340d+01
COFD(3718) = 4.41511629d+00
COFD(3719) = -2.84086963d-01
COFD(3720) = 9.37586971d-03
COFD(3721) = -2.01015340d+01
COFD(3722) = 4.41511629d+00
COFD(3723) = -2.84086963d-01
COFD(3724) = 9.37586971d-03
COFD(3725) = -1.91112931d+01
COFD(3726) = 4.86821670d+00
COFD(3727) = -4.13144121d-01
COFD(3728) = 1.77546701d-02
COFD(3729) = -1.93845675d+01
COFD(3730) = 5.02567894d+00
COFD(3731) = -4.32045169d-01
COFD(3732) = 1.85132214d-02
COFD(3733) = -1.94059889d+01
COFD(3734) = 5.02567894d+00
COFD(3735) = -4.32045169d-01
COFD(3736) = 1.85132214d-02
COFD(3737) = -1.89490450d+01
COFD(3738) = 4.11490576d+00
COFD(3739) = -2.38333952d-01
COFD(3740) = 7.12818549d-03
COFD(3741) = -1.94689917d+01
COFD(3742) = 4.86821670d+00
COFD(3743) = -4.13144121d-01
COFD(3744) = 1.77546701d-02
COFD(3745) = -2.05008516d+01
COFD(3746) = 5.23112374d+00
COFD(3747) = -4.54967682d-01
COFD(3748) = 1.93570423d-02
COFD(3749) = -2.19919149d+01
COFD(3750) = 5.58446511d+00
COFD(3751) = -4.79399331d-01
COFD(3752) = 1.95652693d-02
COFD(3753) = -2.20062584d+01
COFD(3754) = 5.49642957d+00
COFD(3755) = -4.61132993d-01
COFD(3756) = 1.85004773d-02
COFD(3757) = -2.09941194d+01
COFD(3758) = 5.39401640d+00
COFD(3759) = -4.72026580d-01
COFD(3760) = 1.99336592d-02
COFD(3761) = -1.96633368d+01
COFD(3762) = 4.94225344d+00
COFD(3763) = -4.22163426d-01
COFD(3764) = 1.81223730d-02
COFD(3765) = -1.82277241d+01
COFD(3766) = 3.59760627d+00
COFD(3767) = -1.59809307d-01
COFD(3768) = 3.28676045d-03
COFD(3769) = -1.82372764d+01
COFD(3770) = 3.59760627d+00
COFD(3771) = -1.59809307d-01
COFD(3772) = 3.28676045d-03
COFD(3773) = -2.14303625d+01
COFD(3774) = 5.59268435d+00
COFD(3775) = -4.91232974d-01
COFD(3776) = 2.05064746d-02
COFD(3777) = -2.20016256d+01
COFD(3778) = 5.49642957d+00
COFD(3779) = -4.61132993d-01
COFD(3780) = 1.85004773d-02
COFD(3781) = -2.20016256d+01
COFD(3782) = 5.49642957d+00
COFD(3783) = -4.61132993d-01
COFD(3784) = 1.85004773d-02
COFD(3785) = -2.20016256d+01
COFD(3786) = 5.49642957d+00
COFD(3787) = -4.61132993d-01
COFD(3788) = 1.85004773d-02
COFD(3789) = -2.19967195d+01
COFD(3790) = 5.49642957d+00
COFD(3791) = -4.61132993d-01
COFD(3792) = 1.85004773d-02
COFD(3793) = -2.04833713d+01
COFD(3794) = 5.23112374d+00
COFD(3795) = -4.54967682d-01
COFD(3796) = 1.93570423d-02
COFD(3797) = -2.14353267d+01
COFD(3798) = 5.52240865d+00
COFD(3799) = -4.84699537d-01
COFD(3800) = 2.03247833d-02
COFD(3801) = -2.21047682d+01
COFD(3802) = 5.36053938d+00
COFD(3803) = -4.36434519d-01
COFD(3804) = 1.71484255d-02
COFD(3805) = -2.21094839d+01
COFD(3806) = 5.36053938d+00
COFD(3807) = -4.36434519d-01
COFD(3808) = 1.71484255d-02
COFD(3809) = -2.01064363d+01
COFD(3810) = 4.41511629d+00
COFD(3811) = -2.84086963d-01
COFD(3812) = 9.37586971d-03
COFD(3813) = -2.01111595d+01
COFD(3814) = 4.41511629d+00
COFD(3815) = -2.84086963d-01
COFD(3816) = 9.37586971d-03
COFD(3817) = -1.58456300d+01
COFD(3818) = 4.02074783d+00
COFD(3819) = -3.10018522d-01
COFD(3820) = 1.35599552d-02
COFD(3821) = -1.92718582d+01
COFD(3822) = 5.41172124d+00
COFD(3823) = -4.73213887d-01
COFD(3824) = 1.99405473d-02
COFD(3825) = -1.88179418d+01
COFD(3826) = 4.79683898d+00
COFD(3827) = -4.04829719d-01
COFD(3828) = 1.74325475d-02
COFD(3829) = -2.04928958d+01
COFD(3830) = 5.22397933d+00
COFD(3831) = -4.54138171d-01
COFD(3832) = 1.93249285d-02
COFD(3833) = -1.88378874d+01
COFD(3834) = 4.79683898d+00
COFD(3835) = -4.04829719d-01
COFD(3836) = 1.74325475d-02
COFD(3837) = -1.65295288d+01
COFD(3838) = 2.97569206d+00
COFD(3839) = -6.75652842d-02
COFD(3840) = -1.08648422d-03
COFD(3841) = -2.02637994d+01
COFD(3842) = 5.14984081d+00
COFD(3843) = -4.46093018d-01
COFD(3844) = 1.90396647d-02
COFD(3845) = -2.02710316d+01
COFD(3846) = 5.14984081d+00
COFD(3847) = -4.46093018d-01
COFD(3848) = 1.90396647d-02
COFD(3849) = -1.85157414d+01
COFD(3850) = 4.67076124d+00
COFD(3851) = -3.90022427d-01
COFD(3852) = 1.68533953d-02
COFD(3853) = -1.87476063d+01
COFD(3854) = 4.79683898d+00
COFD(3855) = -4.04829719d-01
COFD(3856) = 1.74325475d-02
COFD(3857) = -2.09376196d+01
COFD(3858) = 5.40870099d+00
COFD(3859) = -4.73017610d-01
COFD(3860) = 1.99399066d-02
COFD(3861) = -2.09376196d+01
COFD(3862) = 5.40870099d+00
COFD(3863) = -4.73017610d-01
COFD(3864) = 1.99399066d-02
COFD(3865) = -2.09612557d+01
COFD(3866) = 5.40870099d+00
COFD(3867) = -4.73017610d-01
COFD(3868) = 1.99399066d-02
COFD(3869) = -2.11381508d+01
COFD(3870) = 5.45574440d+00
COFD(3871) = -4.77436155d-01
COFD(3872) = 2.00644596d-02
COFD(3873) = -2.02642227d+01
COFD(3874) = 5.14499740d+00
COFD(3875) = -4.45694430d-01
COFD(3876) = 1.90318646d-02
COFD(3877) = -2.20421041d+01
COFD(3878) = 5.52708332d+00
COFD(3879) = -4.68000808d-01
COFD(3880) = 1.89131908d-02
COFD(3881) = -2.01801667d+01
COFD(3882) = 4.53183330d+00
COFD(3883) = -3.02186760d-01
COFD(3884) = 1.02756490d-02
COFD(3885) = -2.01889168d+01
COFD(3886) = 4.53183330d+00
COFD(3887) = -3.02186760d-01
COFD(3888) = 1.02756490d-02
COFD(3889) = -1.95877017d+01
COFD(3890) = 4.27643051d+00
COFD(3891) = -2.68040901d-01
COFD(3892) = 8.77650113d-03
COFD(3893) = -1.95877017d+01
COFD(3894) = 4.27643051d+00
COFD(3895) = -2.68040901d-01
COFD(3896) = 8.77650113d-03
COFD(3897) = -2.03599050d+01
COFD(3898) = 4.60682543d+00
COFD(3899) = -3.13971634d-01
COFD(3900) = 1.08661011d-02
COFD(3901) = -2.19383403d+01
COFD(3902) = 5.59157589d+00
COFD(3903) = -4.85617912d-01
COFD(3904) = 2.00461138d-02
COFD(3905) = -2.19491710d+01
COFD(3906) = 5.59157589d+00
COFD(3907) = -4.85617912d-01
COFD(3908) = 2.00461138d-02
COFD(3909) = -2.19594078d+01
COFD(3910) = 5.59157589d+00
COFD(3911) = -4.85617912d-01
COFD(3912) = 2.00461138d-02
COFD(3913) = -2.19670848d+01
COFD(3914) = 5.48847873d+00
COFD(3915) = -4.59558930d-01
COFD(3916) = 1.84107961d-02
COFD(3917) = -2.21070030d+01
COFD(3918) = 5.55072945d+00
COFD(3919) = -4.72525345d-01
COFD(3920) = 1.91674202d-02
COFD(3921) = -2.21157340d+01
COFD(3922) = 5.55072945d+00
COFD(3923) = -4.72525345d-01
COFD(3924) = 1.91674202d-02
COFD(3925) = -2.09241647d+01
COFD(3926) = 5.42316225d+00
COFD(3927) = -4.73702801d-01
COFD(3928) = 1.99217718d-02
COFD(3929) = -2.09222454d+01
COFD(3930) = 4.82184721d+00
COFD(3931) = -3.48128875d-01
COFD(3932) = 1.25918978d-02
COFD(3933) = -2.09222454d+01
COFD(3934) = 4.82184721d+00
COFD(3935) = -3.48128875d-01
COFD(3936) = 1.25918978d-02
COFD(3937) = -1.85699333d+01
COFD(3938) = 4.67076124d+00
COFD(3939) = -3.90022427d-01
COFD(3940) = 1.68533953d-02
COFD(3941) = -1.87654600d+01
COFD(3942) = 4.79683898d+00
COFD(3943) = -4.04829719d-01
COFD(3944) = 1.74325475d-02
COFD(3945) = -1.97755767d+01
COFD(3946) = 5.13121203d+00
COFD(3947) = -4.44310112d-01
COFD(3948) = 1.89883544d-02
COFD(3949) = -1.82953333d+01
COFD(3950) = 3.82052409d+00
COFD(3951) = -1.97327882d-01
COFD(3952) = 5.27318631d-03
COFD(3953) = -1.89285474d+01
COFD(3954) = 4.67076124d+00
COFD(3955) = -3.90022427d-01
COFD(3956) = 1.68533953d-02
COFD(3957) = -2.02446539d+01
COFD(3958) = 5.13632093d+00
COFD(3959) = -4.44839124d-01
COFD(3960) = 1.90058354d-02
COFD(3961) = -2.19270918d+01
COFD(3962) = 5.60963444d+00
COFD(3963) = -4.89805286d-01
COFD(3964) = 2.03017681d-02
COFD(3965) = -2.20870532d+01
COFD(3966) = 5.58518321d+00
COFD(3967) = -4.80534235d-01
COFD(3968) = 1.96556393d-02
COFD(3969) = -2.05164480d+01
COFD(3970) = 5.23355494d+00
COFD(3971) = -4.55249645d-01
COFD(3972) = 1.93679438d-02
COFD(3973) = -1.90658825d+01
COFD(3974) = 4.72218939d+00
COFD(3975) = -3.96003719d-01
COFD(3976) = 1.70845726d-02
COFD(3977) = -1.93898781d+01
COFD(3978) = 4.14577053d+00
COFD(3979) = -2.43003893d-01
COFD(3980) = 7.35638742d-03
COFD(3981) = -1.93995811d+01
COFD(3982) = 4.14577053d+00
COFD(3983) = -2.43003893d-01
COFD(3984) = 7.35638742d-03
COFD(3985) = -2.09241797d+01
COFD(3986) = 5.42316225d+00
COFD(3987) = -4.73702801d-01
COFD(3988) = 1.99217718d-02
COFD(3989) = -2.20823296d+01
COFD(3990) = 5.58518321d+00
COFD(3991) = -4.80534235d-01
COFD(3992) = 1.96556393d-02
COFD(3993) = -2.20823296d+01
COFD(3994) = 5.58518321d+00
COFD(3995) = -4.80534235d-01
COFD(3996) = 1.96556393d-02
COFD(3997) = -2.20823296d+01
COFD(3998) = 5.58518321d+00
COFD(3999) = -4.80534235d-01
COFD(4000) = 1.96556393d-02
COFD(4001) = -2.20773284d+01
COFD(4002) = 5.58518321d+00
COFD(4003) = -4.80534235d-01
COFD(4004) = 1.96556393d-02
COFD(4005) = -2.02268902d+01
COFD(4006) = 5.13632093d+00
COFD(4007) = -4.44839124d-01
COFD(4008) = 1.90058354d-02
COFD(4009) = -2.10026861d+01
COFD(4010) = 5.38326647d+00
COFD(4011) = -4.71201048d-01
COFD(4012) = 1.99207516d-02
COFD(4013) = -2.23155854d+01
COFD(4014) = 5.50872500d+00
COFD(4015) = -4.64427323d-01
COFD(4016) = 1.87105612d-02
COFD(4017) = -2.23203936d+01
COFD(4018) = 5.50872500d+00
COFD(4019) = -4.64427323d-01
COFD(4020) = 1.87105612d-02
COFD(4021) = -2.09272429d+01
COFD(4022) = 4.82184721d+00
COFD(4023) = -3.48128875d-01
COFD(4024) = 1.25918978d-02
COFD(4025) = -2.09320587d+01
COFD(4026) = 4.82184721d+00
COFD(4027) = -3.48128875d-01
COFD(4028) = 1.25918978d-02
COFD(4029) = -1.58456300d+01
COFD(4030) = 4.02074783d+00
COFD(4031) = -3.10018522d-01
COFD(4032) = 1.35599552d-02
COFD(4033) = -1.92718582d+01
COFD(4034) = 5.41172124d+00
COFD(4035) = -4.73213887d-01
COFD(4036) = 1.99405473d-02
COFD(4037) = -1.88179418d+01
COFD(4038) = 4.79683898d+00
COFD(4039) = -4.04829719d-01
COFD(4040) = 1.74325475d-02
COFD(4041) = -2.04928958d+01
COFD(4042) = 5.22397933d+00
COFD(4043) = -4.54138171d-01
COFD(4044) = 1.93249285d-02
COFD(4045) = -1.88378874d+01
COFD(4046) = 4.79683898d+00
COFD(4047) = -4.04829719d-01
COFD(4048) = 1.74325475d-02
COFD(4049) = -1.65295288d+01
COFD(4050) = 2.97569206d+00
COFD(4051) = -6.75652842d-02
COFD(4052) = -1.08648422d-03
COFD(4053) = -2.02637994d+01
COFD(4054) = 5.14984081d+00
COFD(4055) = -4.46093018d-01
COFD(4056) = 1.90396647d-02
COFD(4057) = -2.02710316d+01
COFD(4058) = 5.14984081d+00
COFD(4059) = -4.46093018d-01
COFD(4060) = 1.90396647d-02
COFD(4061) = -1.85157414d+01
COFD(4062) = 4.67076124d+00
COFD(4063) = -3.90022427d-01
COFD(4064) = 1.68533953d-02
COFD(4065) = -1.87476063d+01
COFD(4066) = 4.79683898d+00
COFD(4067) = -4.04829719d-01
COFD(4068) = 1.74325475d-02
COFD(4069) = -2.09376196d+01
COFD(4070) = 5.40870099d+00
COFD(4071) = -4.73017610d-01
COFD(4072) = 1.99399066d-02
COFD(4073) = -2.09376196d+01
COFD(4074) = 5.40870099d+00
COFD(4075) = -4.73017610d-01
COFD(4076) = 1.99399066d-02
COFD(4077) = -2.09612557d+01
COFD(4078) = 5.40870099d+00
COFD(4079) = -4.73017610d-01
COFD(4080) = 1.99399066d-02
COFD(4081) = -2.11381508d+01
COFD(4082) = 5.45574440d+00
COFD(4083) = -4.77436155d-01
COFD(4084) = 2.00644596d-02
COFD(4085) = -2.02642227d+01
COFD(4086) = 5.14499740d+00
COFD(4087) = -4.45694430d-01
COFD(4088) = 1.90318646d-02
COFD(4089) = -2.20421041d+01
COFD(4090) = 5.52708332d+00
COFD(4091) = -4.68000808d-01
COFD(4092) = 1.89131908d-02
COFD(4093) = -2.01801667d+01
COFD(4094) = 4.53183330d+00
COFD(4095) = -3.02186760d-01
COFD(4096) = 1.02756490d-02
COFD(4097) = -2.01889168d+01
COFD(4098) = 4.53183330d+00
COFD(4099) = -3.02186760d-01
COFD(4100) = 1.02756490d-02
COFD(4101) = -1.95877017d+01
COFD(4102) = 4.27643051d+00
COFD(4103) = -2.68040901d-01
COFD(4104) = 8.77650113d-03
COFD(4105) = -1.95877017d+01
COFD(4106) = 4.27643051d+00
COFD(4107) = -2.68040901d-01
COFD(4108) = 8.77650113d-03
COFD(4109) = -2.03599050d+01
COFD(4110) = 4.60682543d+00
COFD(4111) = -3.13971634d-01
COFD(4112) = 1.08661011d-02
COFD(4113) = -2.19383403d+01
COFD(4114) = 5.59157589d+00
COFD(4115) = -4.85617912d-01
COFD(4116) = 2.00461138d-02
COFD(4117) = -2.19491710d+01
COFD(4118) = 5.59157589d+00
COFD(4119) = -4.85617912d-01
COFD(4120) = 2.00461138d-02
COFD(4121) = -2.19594078d+01
COFD(4122) = 5.59157589d+00
COFD(4123) = -4.85617912d-01
COFD(4124) = 2.00461138d-02
COFD(4125) = -2.19670848d+01
COFD(4126) = 5.48847873d+00
COFD(4127) = -4.59558930d-01
COFD(4128) = 1.84107961d-02
COFD(4129) = -2.21070030d+01
COFD(4130) = 5.55072945d+00
COFD(4131) = -4.72525345d-01
COFD(4132) = 1.91674202d-02
COFD(4133) = -2.21157340d+01
COFD(4134) = 5.55072945d+00
COFD(4135) = -4.72525345d-01
COFD(4136) = 1.91674202d-02
COFD(4137) = -2.09241647d+01
COFD(4138) = 5.42316225d+00
COFD(4139) = -4.73702801d-01
COFD(4140) = 1.99217718d-02
COFD(4141) = -2.09222454d+01
COFD(4142) = 4.82184721d+00
COFD(4143) = -3.48128875d-01
COFD(4144) = 1.25918978d-02
COFD(4145) = -2.09222454d+01
COFD(4146) = 4.82184721d+00
COFD(4147) = -3.48128875d-01
COFD(4148) = 1.25918978d-02
COFD(4149) = -1.85699333d+01
COFD(4150) = 4.67076124d+00
COFD(4151) = -3.90022427d-01
COFD(4152) = 1.68533953d-02
COFD(4153) = -1.87654600d+01
COFD(4154) = 4.79683898d+00
COFD(4155) = -4.04829719d-01
COFD(4156) = 1.74325475d-02
COFD(4157) = -1.97755767d+01
COFD(4158) = 5.13121203d+00
COFD(4159) = -4.44310112d-01
COFD(4160) = 1.89883544d-02
COFD(4161) = -1.82953333d+01
COFD(4162) = 3.82052409d+00
COFD(4163) = -1.97327882d-01
COFD(4164) = 5.27318631d-03
COFD(4165) = -1.89285474d+01
COFD(4166) = 4.67076124d+00
COFD(4167) = -3.90022427d-01
COFD(4168) = 1.68533953d-02
COFD(4169) = -2.02446539d+01
COFD(4170) = 5.13632093d+00
COFD(4171) = -4.44839124d-01
COFD(4172) = 1.90058354d-02
COFD(4173) = -2.19270918d+01
COFD(4174) = 5.60963444d+00
COFD(4175) = -4.89805286d-01
COFD(4176) = 2.03017681d-02
COFD(4177) = -2.20870532d+01
COFD(4178) = 5.58518321d+00
COFD(4179) = -4.80534235d-01
COFD(4180) = 1.96556393d-02
COFD(4181) = -2.05164480d+01
COFD(4182) = 5.23355494d+00
COFD(4183) = -4.55249645d-01
COFD(4184) = 1.93679438d-02
COFD(4185) = -1.90658825d+01
COFD(4186) = 4.72218939d+00
COFD(4187) = -3.96003719d-01
COFD(4188) = 1.70845726d-02
COFD(4189) = -1.93898781d+01
COFD(4190) = 4.14577053d+00
COFD(4191) = -2.43003893d-01
COFD(4192) = 7.35638742d-03
COFD(4193) = -1.93995811d+01
COFD(4194) = 4.14577053d+00
COFD(4195) = -2.43003893d-01
COFD(4196) = 7.35638742d-03
COFD(4197) = -2.09241797d+01
COFD(4198) = 5.42316225d+00
COFD(4199) = -4.73702801d-01
COFD(4200) = 1.99217718d-02
COFD(4201) = -2.20823296d+01
COFD(4202) = 5.58518321d+00
COFD(4203) = -4.80534235d-01
COFD(4204) = 1.96556393d-02
COFD(4205) = -2.20823296d+01
COFD(4206) = 5.58518321d+00
COFD(4207) = -4.80534235d-01
COFD(4208) = 1.96556393d-02
COFD(4209) = -2.20823296d+01
COFD(4210) = 5.58518321d+00
COFD(4211) = -4.80534235d-01
COFD(4212) = 1.96556393d-02
COFD(4213) = -2.20773284d+01
COFD(4214) = 5.58518321d+00
COFD(4215) = -4.80534235d-01
COFD(4216) = 1.96556393d-02
COFD(4217) = -2.02268902d+01
COFD(4218) = 5.13632093d+00
COFD(4219) = -4.44839124d-01
COFD(4220) = 1.90058354d-02
COFD(4221) = -2.10026861d+01
COFD(4222) = 5.38326647d+00
COFD(4223) = -4.71201048d-01
COFD(4224) = 1.99207516d-02
COFD(4225) = -2.23155854d+01
COFD(4226) = 5.50872500d+00
COFD(4227) = -4.64427323d-01
COFD(4228) = 1.87105612d-02
COFD(4229) = -2.23203936d+01
COFD(4230) = 5.50872500d+00
COFD(4231) = -4.64427323d-01
COFD(4232) = 1.87105612d-02
COFD(4233) = -2.09272429d+01
COFD(4234) = 4.82184721d+00
COFD(4235) = -3.48128875d-01
COFD(4236) = 1.25918978d-02
COFD(4237) = -2.09320587d+01
COFD(4238) = 4.82184721d+00
COFD(4239) = -3.48128875d-01
COFD(4240) = 1.25918978d-02
COFD(4241) = -1.59537247d+01
COFD(4242) = 4.07051484d+00
COFD(4243) = -3.16303109d-01
COFD(4244) = 1.38259377d-02
COFD(4245) = -1.96866103d+01
COFD(4246) = 5.54637286d+00
COFD(4247) = -4.87070324d-01
COFD(4248) = 2.03983467d-02
COFD(4249) = -1.93364585d+01
COFD(4250) = 4.98286777d+00
COFD(4251) = -4.26970814d-01
COFD(4252) = 1.83122917d-02
COFD(4253) = -2.07595845d+01
COFD(4254) = 5.32244593d+00
COFD(4255) = -4.65829403d-01
COFD(4256) = 1.97895274d-02
COFD(4257) = -1.93566243d+01
COFD(4258) = 4.98286777d+00
COFD(4259) = -4.26970814d-01
COFD(4260) = 1.83122917d-02
COFD(4261) = -1.80253664d+01
COFD(4262) = 3.69199168d+00
COFD(4263) = -1.74005516d-01
COFD(4264) = 3.97694372d-03
COFD(4265) = -2.07672833d+01
COFD(4266) = 5.32244593d+00
COFD(4267) = -4.65829403d-01
COFD(4268) = 1.97895274d-02
COFD(4269) = -2.07746356d+01
COFD(4270) = 5.32244593d+00
COFD(4271) = -4.65829403d-01
COFD(4272) = 1.97895274d-02
COFD(4273) = -1.89671752d+01
COFD(4274) = 4.83076737d+00
COFD(4275) = -4.08802573d-01
COFD(4276) = 1.75875241d-02
COFD(4277) = -1.92654138d+01
COFD(4278) = 4.98286777d+00
COFD(4279) = -4.26970814d-01
COFD(4280) = 1.83122917d-02
COFD(4281) = -2.13538553d+01
COFD(4282) = 5.54007827d+00
COFD(4283) = -4.86434511d-01
COFD(4284) = 2.03779006d-02
COFD(4285) = -2.13538553d+01
COFD(4286) = 5.54007827d+00
COFD(4287) = -4.86434511d-01
COFD(4288) = 2.03779006d-02
COFD(4289) = -2.13777308d+01
COFD(4290) = 5.54007827d+00
COFD(4291) = -4.86434511d-01
COFD(4292) = 2.03779006d-02
COFD(4293) = -2.13319784d+01
COFD(4294) = 5.52422470d+00
COFD(4295) = -4.84872944d-01
COFD(4296) = 2.03298213d-02
COFD(4297) = -2.04144604d+01
COFD(4298) = 5.19614628d+00
COFD(4299) = -4.50889164d-01
COFD(4300) = 1.91983328d-02
COFD(4301) = -2.20063594d+01
COFD(4302) = 5.48540187d+00
COFD(4303) = -4.58962148d-01
COFD(4304) = 1.83770355d-02
COFD(4305) = -1.93125662d+01
COFD(4306) = 4.10954793d+00
COFD(4307) = -2.37523329d-01
COFD(4308) = 7.08858141d-03
COFD(4309) = -1.93214527d+01
COFD(4310) = 4.10954793d+00
COFD(4311) = -2.37523329d-01
COFD(4312) = 7.08858141d-03
COFD(4313) = -2.03599050d+01
COFD(4314) = 4.60682543d+00
COFD(4315) = -3.13971634d-01
COFD(4316) = 1.08661011d-02
COFD(4317) = -2.03599050d+01
COFD(4318) = 4.60682543d+00
COFD(4319) = -3.13971634d-01
COFD(4320) = 1.08661011d-02
COFD(4321) = -1.95785144d+01
COFD(4322) = 4.22062499d+00
COFD(4323) = -2.54326872d-01
COFD(4324) = 7.91017784d-03
COFD(4325) = -2.20378456d+01
COFD(4326) = 5.58129885d+00
COFD(4327) = -4.78532921d-01
COFD(4328) = 1.95095699d-02
COFD(4329) = -2.20488323d+01
COFD(4330) = 5.58129885d+00
COFD(4331) = -4.78532921d-01
COFD(4332) = 1.95095699d-02
COFD(4333) = -2.20592197d+01
COFD(4334) = 5.58129885d+00
COFD(4335) = -4.78532921d-01
COFD(4336) = 1.95095699d-02
COFD(4337) = -2.17419207d+01
COFD(4338) = 5.33732875d+00
COFD(4339) = -4.32453425d-01
COFD(4340) = 1.69373833d-02
COFD(4341) = -2.20028184d+01
COFD(4342) = 5.45178028d+00
COFD(4343) = -4.52847771d-01
COFD(4344) = 1.80418544d-02
COFD(4345) = -2.20116855d+01
COFD(4346) = 5.45178028d+00
COFD(4347) = -4.52847771d-01
COFD(4348) = 1.80418544d-02
COFD(4349) = -2.13796303d+01
COFD(4350) = 5.56978987d+00
COFD(4351) = -4.89141980d-01
COFD(4352) = 2.04499210d-02
COFD(4353) = -2.03036402d+01
COFD(4354) = 4.50250781d+00
COFD(4355) = -2.97622106d-01
COFD(4356) = 1.00481473d-02
COFD(4357) = -2.03036402d+01
COFD(4358) = 4.50250781d+00
COFD(4359) = -2.97622106d-01
COFD(4360) = 1.00481473d-02
COFD(4361) = -1.90218743d+01
COFD(4362) = 4.83076737d+00
COFD(4363) = -4.08802573d-01
COFD(4364) = 1.75875241d-02
COFD(4365) = -1.92834357d+01
COFD(4366) = 4.98286777d+00
COFD(4367) = -4.26970814d-01
COFD(4368) = 1.83122917d-02
COFD(4369) = -1.93053262d+01
COFD(4370) = 4.98286777d+00
COFD(4371) = -4.26970814d-01
COFD(4372) = 1.83122917d-02
COFD(4373) = -1.92028153d+01
COFD(4374) = 4.22627146d+00
COFD(4375) = -2.55181433d-01
COFD(4376) = 7.95198535d-03
COFD(4377) = -1.93844661d+01
COFD(4378) = 4.83076737d+00
COFD(4379) = -4.08802573d-01
COFD(4380) = 1.75875241d-02
COFD(4381) = -2.04024630d+01
COFD(4382) = 5.18856872d+00
COFD(4383) = -4.50001829d-01
COFD(4384) = 1.91636142d-02
COFD(4385) = -2.19932321d+01
COFD(4386) = 5.58514538d+00
COFD(4387) = -4.80745077d-01
COFD(4388) = 1.96733087d-02
COFD(4389) = -2.20260694d+01
COFD(4390) = 5.50583166d+00
COFD(4391) = -4.63753262d-01
COFD(4392) = 1.86693462d-02
COFD(4393) = -2.09549117d+01
COFD(4394) = 5.37710408d+00
COFD(4395) = -4.70743467d-01
COFD(4396) = 1.99147079d-02
COFD(4397) = -1.95512940d+01
COFD(4398) = 4.89369650d+00
COFD(4399) = -4.16273855d-01
COFD(4400) = 1.78833908d-02
COFD(4401) = -1.84848887d+01
COFD(4402) = 3.70911349d+00
COFD(4403) = -1.76610648d-01
COFD(4404) = 4.10440582d-03
COFD(4405) = -1.84947374d+01
COFD(4406) = 3.70911349d+00
COFD(4407) = -1.76610648d-01
COFD(4408) = 4.10440582d-03
COFD(4409) = -2.13796455d+01
COFD(4410) = 5.56978987d+00
COFD(4411) = -4.89141980d-01
COFD(4412) = 2.04499210d-02
COFD(4413) = -2.20212574d+01
COFD(4414) = 5.50583166d+00
COFD(4415) = -4.63753262d-01
COFD(4416) = 1.86693462d-02
COFD(4417) = -2.20212574d+01
COFD(4418) = 5.50583166d+00
COFD(4419) = -4.63753262d-01
COFD(4420) = 1.86693462d-02
COFD(4421) = -2.20212574d+01
COFD(4422) = 5.50583166d+00
COFD(4423) = -4.63753262d-01
COFD(4424) = 1.86693462d-02
COFD(4425) = -2.20161635d+01
COFD(4426) = 5.50583166d+00
COFD(4427) = -4.63753262d-01
COFD(4428) = 1.86693462d-02
COFD(4429) = -2.03844252d+01
COFD(4430) = 5.18856872d+00
COFD(4431) = -4.50001829d-01
COFD(4432) = 1.91636142d-02
COFD(4433) = -2.13502344d+01
COFD(4434) = 5.48617067d+00
COFD(4435) = -4.80816776d-01
COFD(4436) = 2.01887868d-02
COFD(4437) = -2.21897029d+01
COFD(4438) = 5.39861129d+00
COFD(4439) = -4.43119198d-01
COFD(4440) = 1.75075657d-02
COFD(4441) = -2.21946011d+01
COFD(4442) = 5.39861129d+00
COFD(4443) = -4.43119198d-01
COFD(4444) = 1.75075657d-02
COFD(4445) = -2.03087302d+01
COFD(4446) = 4.50250781d+00
COFD(4447) = -2.97622106d-01
COFD(4448) = 1.00481473d-02
COFD(4449) = -2.03136362d+01
COFD(4450) = 4.50250781d+00
COFD(4451) = -2.97622106d-01
COFD(4452) = 1.00481473d-02
COFD(4453) = -1.34695359d+01
COFD(4454) = 3.09379603d+00
COFD(4455) = -1.91268635d-01
COFD(4456) = 8.47480224d-03
COFD(4457) = -1.72224724d+01
COFD(4458) = 4.69060745d+00
COFD(4459) = -3.92369888d-01
COFD(4460) = 1.69459661d-02
COFD(4461) = -1.65412306d+01
COFD(4462) = 3.95035840d+00
COFD(4463) = -3.00959418d-01
COFD(4464) = 1.31692593d-02
COFD(4465) = -1.78725135d+01
COFD(4466) = 4.29613154d+00
COFD(4467) = -3.44012526d-01
COFD(4468) = 1.49643715d-02
COFD(4469) = -1.65596434d+01
COFD(4470) = 3.95035840d+00
COFD(4471) = -3.00959418d-01
COFD(4472) = 1.31692593d-02
COFD(4473) = -2.15014310d+01
COFD(4474) = 5.46737673d+00
COFD(4475) = -4.55696085d-01
COFD(4476) = 1.81982625d-02
COFD(4477) = -1.78792605d+01
COFD(4478) = 4.29613154d+00
COFD(4479) = -3.44012526d-01
COFD(4480) = 1.49643715d-02
COFD(4481) = -1.78856918d+01
COFD(4482) = 4.29613154d+00
COFD(4483) = -3.44012526d-01
COFD(4484) = 1.49643715d-02
COFD(4485) = -1.61038806d+01
COFD(4486) = 3.75910622d+00
COFD(4487) = -2.75986578d-01
COFD(4488) = 1.20782843d-02
COFD(4489) = -1.64758697d+01
COFD(4490) = 3.95035840d+00
COFD(4491) = -3.00959418d-01
COFD(4492) = 1.31692593d-02
COFD(4493) = -1.88171077d+01
COFD(4494) = 4.68393046d+00
COFD(4495) = -3.91610863d-01
COFD(4496) = 1.69174645d-02
COFD(4497) = -1.88171077d+01
COFD(4498) = 4.68393046d+00
COFD(4499) = -3.91610863d-01
COFD(4500) = 1.69174645d-02
COFD(4501) = -1.88390649d+01
COFD(4502) = 4.68393046d+00
COFD(4503) = -3.91610863d-01
COFD(4504) = 1.69174645d-02
COFD(4505) = -1.87821119d+01
COFD(4506) = 4.66162351d+00
COFD(4507) = -3.88920477d-01
COFD(4508) = 1.68089648d-02
COFD(4509) = -1.76182365d+01
COFD(4510) = 4.19935698d+00
COFD(4511) = -3.32310212d-01
COFD(4512) = 1.44920670d-02
COFD(4513) = -2.09066354d+01
COFD(4514) = 5.30153901d+00
COFD(4515) = -4.63335119d-01
COFD(4516) = 1.96897053d-02
COFD(4517) = -2.19851338d+01
COFD(4518) = 5.55935694d+00
COFD(4519) = -4.74154740d-01
COFD(4520) = 1.92584304d-02
COFD(4521) = -2.19929679d+01
COFD(4522) = 5.55935694d+00
COFD(4523) = -4.74154740d-01
COFD(4524) = 1.92584304d-02
COFD(4525) = -2.19383403d+01
COFD(4526) = 5.59157589d+00
COFD(4527) = -4.85617912d-01
COFD(4528) = 2.00461138d-02
COFD(4529) = -2.19383403d+01
COFD(4530) = 5.59157589d+00
COFD(4531) = -4.85617912d-01
COFD(4532) = 2.00461138d-02
COFD(4533) = -2.20378456d+01
COFD(4534) = 5.58129885d+00
COFD(4535) = -4.78532921d-01
COFD(4536) = 1.95095699d-02
COFD(4537) = -2.03569548d+01
COFD(4538) = 5.13263469d+00
COFD(4539) = -4.44457285d-01
COFD(4540) = 1.89932102d-02
COFD(4541) = -2.03667275d+01
COFD(4542) = 5.13263469d+00
COFD(4543) = -4.44457285d-01
COFD(4544) = 1.89932102d-02
COFD(4545) = -2.03759451d+01
COFD(4546) = 5.13263469d+00
COFD(4547) = -4.44457285d-01
COFD(4548) = 1.89932102d-02
COFD(4549) = -2.12018018d+01
COFD(4550) = 5.39823225d+00
COFD(4551) = -4.72294645d-01
COFD(4552) = 1.99340225d-02
COFD(4553) = -2.10759163d+01
COFD(4554) = 5.34286099d+00
COFD(4555) = -4.68100992d-01
COFD(4556) = 1.98731399d-02
COFD(4557) = -2.10837327d+01
COFD(4558) = 5.34286099d+00
COFD(4559) = -4.68100992d-01
COFD(4560) = 1.98731399d-02
COFD(4561) = -1.88524385d+01
COFD(4562) = 4.72476764d+00
COFD(4563) = -3.96306836d-01
COFD(4564) = 1.70964541d-02
COFD(4565) = -2.20942491d+01
COFD(4566) = 5.58360799d+00
COFD(4567) = -4.82701436d-01
COFD(4568) = 1.98437922d-02
COFD(4569) = -2.20942491d+01
COFD(4570) = 5.58360799d+00
COFD(4571) = -4.82701436d-01
COFD(4572) = 1.98437922d-02
COFD(4573) = -1.61544946d+01
COFD(4574) = 3.75910622d+00
COFD(4575) = -2.75986578d-01
COFD(4576) = 1.20782843d-02
COFD(4577) = -1.64922031d+01
COFD(4578) = 3.95035840d+00
COFD(4579) = -3.00959418d-01
COFD(4580) = 1.31692593d-02
COFD(4581) = -1.65122609d+01
COFD(4582) = 3.95035840d+00
COFD(4583) = -3.00959418d-01
COFD(4584) = 1.31692593d-02
COFD(4585) = -2.16823299d+01
COFD(4586) = 5.58175603d+00
COFD(4587) = -4.78660918d-01
COFD(4588) = 1.95178500d-02
COFD(4589) = -1.64868273d+01
COFD(4590) = 3.75910622d+00
COFD(4591) = -2.75986578d-01
COFD(4592) = 1.20782843d-02
COFD(4593) = -1.76057946d+01
COFD(4594) = 4.19171952d+00
COFD(4595) = -3.31354810d-01
COFD(4596) = 1.44520623d-02
COFD(4597) = -2.01897084d+01
COFD(4598) = 5.08408346d+00
COFD(4599) = -4.38907391d-01
COFD(4600) = 1.87824789d-02
COFD(4601) = -2.07211857d+01
COFD(4602) = 5.23116678d+00
COFD(4603) = -4.54972675d-01
COFD(4604) = 1.93572354d-02
COFD(4605) = -1.81843740d+01
COFD(4606) = 4.40883885d+00
COFD(4607) = -3.58021243d-01
COFD(4608) = 1.55467750d-02
COFD(4609) = -1.67047589d+01
COFD(4610) = 3.84749999d+00
COFD(4611) = -2.87575310d-01
COFD(4612) = 1.25862768d-02
COFD(4613) = -2.19022500d+01
COFD(4614) = 5.47350698d+00
COFD(4615) = -4.56806529d-01
COFD(4616) = 1.82590273d-02
COFD(4617) = -2.19109699d+01
COFD(4618) = 5.47350698d+00
COFD(4619) = -4.56806529d-01
COFD(4620) = 1.82590273d-02
COFD(4621) = -1.88524516d+01
COFD(4622) = 4.72476764d+00
COFD(4623) = -3.96306836d-01
COFD(4624) = 1.70964541d-02
COFD(4625) = -2.07170423d+01
COFD(4626) = 5.23116678d+00
COFD(4627) = -4.54972675d-01
COFD(4628) = 1.93572354d-02
COFD(4629) = -2.07170423d+01
COFD(4630) = 5.23116678d+00
COFD(4631) = -4.54972675d-01
COFD(4632) = 1.93572354d-02
COFD(4633) = -2.07170423d+01
COFD(4634) = 5.23116678d+00
COFD(4635) = -4.54972675d-01
COFD(4636) = 1.93572354d-02
COFD(4637) = -2.07126500d+01
COFD(4638) = 5.23116678d+00
COFD(4639) = -4.54972675d-01
COFD(4640) = 1.93572354d-02
COFD(4641) = -1.75898751d+01
COFD(4642) = 4.19171952d+00
COFD(4643) = -3.31354810d-01
COFD(4644) = 1.44520623d-02
COFD(4645) = -1.87596895d+01
COFD(4646) = 4.61078776d+00
COFD(4647) = -3.82625667d-01
COFD(4648) = 1.65478601d-02
COFD(4649) = -2.14290825d+01
COFD(4650) = 5.37331605d+00
COFD(4651) = -4.70491203d-01
COFD(4652) = 1.99134666d-02
COFD(4653) = -2.14332997d+01
COFD(4654) = 5.37331605d+00
COFD(4655) = -4.70491203d-01
COFD(4656) = 1.99134666d-02
COFD(4657) = -2.20986379d+01
COFD(4658) = 5.58360799d+00
COFD(4659) = -4.82701436d-01
COFD(4660) = 1.98437922d-02
COFD(4661) = -2.21028621d+01
COFD(4662) = 5.58360799d+00
COFD(4663) = -4.82701436d-01
COFD(4664) = 1.98437922d-02
COFD(4665) = -1.34709807d+01
COFD(4666) = 3.09379603d+00
COFD(4667) = -1.91268635d-01
COFD(4668) = 8.47480224d-03
COFD(4669) = -1.72232223d+01
COFD(4670) = 4.69060745d+00
COFD(4671) = -3.92369888d-01
COFD(4672) = 1.69459661d-02
COFD(4673) = -1.65488358d+01
COFD(4674) = 3.95035840d+00
COFD(4675) = -3.00959418d-01
COFD(4676) = 1.31692593d-02
COFD(4677) = -1.78834935d+01
COFD(4678) = 4.29613154d+00
COFD(4679) = -3.44012526d-01
COFD(4680) = 1.49643715d-02
COFD(4681) = -1.65675362d+01
COFD(4682) = 3.95035840d+00
COFD(4683) = -3.00959418d-01
COFD(4684) = 1.31692593d-02
COFD(4685) = -2.15095980d+01
COFD(4686) = 5.46737673d+00
COFD(4687) = -4.55696085d-01
COFD(4688) = 1.81982625d-02
COFD(4689) = -1.78903913d+01
COFD(4690) = 4.29613154d+00
COFD(4691) = -3.44012526d-01
COFD(4692) = 1.49643715d-02
COFD(4693) = -1.78969684d+01
COFD(4694) = 4.29613154d+00
COFD(4695) = -3.44012526d-01
COFD(4696) = 1.49643715d-02
COFD(4697) = -1.61101966d+01
COFD(4698) = 3.75910622d+00
COFD(4699) = -2.75986578d-01
COFD(4700) = 1.20782843d-02
COFD(4701) = -1.64825368d+01
COFD(4702) = 3.95035840d+00
COFD(4703) = -3.00959418d-01
COFD(4704) = 1.31692593d-02
COFD(4705) = -1.88241079d+01
COFD(4706) = 4.68393046d+00
COFD(4707) = -3.91610863d-01
COFD(4708) = 1.69174645d-02
COFD(4709) = -1.88241079d+01
COFD(4710) = 4.68393046d+00
COFD(4711) = -3.91610863d-01
COFD(4712) = 1.69174645d-02
COFD(4713) = -1.88463816d+01
COFD(4714) = 4.68393046d+00
COFD(4715) = -3.91610863d-01
COFD(4716) = 1.69174645d-02
COFD(4717) = -1.87897298d+01
COFD(4718) = 4.66162351d+00
COFD(4719) = -3.88920477d-01
COFD(4720) = 1.68089648d-02
COFD(4721) = -1.76285640d+01
COFD(4722) = 4.19935698d+00
COFD(4723) = -3.32310212d-01
COFD(4724) = 1.44920670d-02
COFD(4725) = -2.09191285d+01
COFD(4726) = 5.30153901d+00
COFD(4727) = -4.63335119d-01
COFD(4728) = 1.96897053d-02
COFD(4729) = -2.19956352d+01
COFD(4730) = 5.55935694d+00
COFD(4731) = -4.74154740d-01
COFD(4732) = 1.92584304d-02
COFD(4733) = -2.20036369d+01
COFD(4734) = 5.55935694d+00
COFD(4735) = -4.74154740d-01
COFD(4736) = 1.92584304d-02
COFD(4737) = -2.19491710d+01
COFD(4738) = 5.59157589d+00
COFD(4739) = -4.85617912d-01
COFD(4740) = 2.00461138d-02
COFD(4741) = -2.19491710d+01
COFD(4742) = 5.59157589d+00
COFD(4743) = -4.85617912d-01
COFD(4744) = 2.00461138d-02
COFD(4745) = -2.20488323d+01
COFD(4746) = 5.58129885d+00
COFD(4747) = -4.78532921d-01
COFD(4748) = 1.95095699d-02
COFD(4749) = -2.03667275d+01
COFD(4750) = 5.13263469d+00
COFD(4751) = -4.44457285d-01
COFD(4752) = 1.89932102d-02
COFD(4753) = -2.03766950d+01
COFD(4754) = 5.13263469d+00
COFD(4755) = -4.44457285d-01
COFD(4756) = 1.89932102d-02
COFD(4757) = -2.03861000d+01
COFD(4758) = 5.13263469d+00
COFD(4759) = -4.44457285d-01
COFD(4760) = 1.89932102d-02
COFD(4761) = -2.12121370d+01
COFD(4762) = 5.39823225d+00
COFD(4763) = -4.72294645d-01
COFD(4764) = 1.99340225d-02
COFD(4765) = -2.10864251d+01
COFD(4766) = 5.34286099d+00
COFD(4767) = -4.68100992d-01
COFD(4768) = 1.98731399d-02
COFD(4769) = -2.10944088d+01
COFD(4770) = 5.34286099d+00
COFD(4771) = -4.68100992d-01
COFD(4772) = 1.98731399d-02
COFD(4773) = -1.88646070d+01
COFD(4774) = 4.72476764d+00
COFD(4775) = -3.96306836d-01
COFD(4776) = 1.70964541d-02
COFD(4777) = -2.21065306d+01
COFD(4778) = 5.58360799d+00
COFD(4779) = -4.82701436d-01
COFD(4780) = 1.98437922d-02
COFD(4781) = -2.21065306d+01
COFD(4782) = 5.58360799d+00
COFD(4783) = -4.82701436d-01
COFD(4784) = 1.98437922d-02
COFD(4785) = -1.61614882d+01
COFD(4786) = 3.75910622d+00
COFD(4787) = -2.75986578d-01
COFD(4788) = 1.20782843d-02
COFD(4789) = -1.64995136d+01
COFD(4790) = 3.95035840d+00
COFD(4791) = -3.00959418d-01
COFD(4792) = 1.31692593d-02
COFD(4793) = -1.65198729d+01
COFD(4794) = 3.95035840d+00
COFD(4795) = -3.00959418d-01
COFD(4796) = 1.31692593d-02
COFD(4797) = -2.16902292d+01
COFD(4798) = 5.58175603d+00
COFD(4799) = -4.78660918d-01
COFD(4800) = 1.95178500d-02
COFD(4801) = -1.64973292d+01
COFD(4802) = 3.75910622d+00
COFD(4803) = -2.75986578d-01
COFD(4804) = 1.20782843d-02
COFD(4805) = -1.76164603d+01
COFD(4806) = 4.19171952d+00
COFD(4807) = -3.31354810d-01
COFD(4808) = 1.44520623d-02
COFD(4809) = -2.02024036d+01
COFD(4810) = 5.08408346d+00
COFD(4811) = -4.38907391d-01
COFD(4812) = 1.87824789d-02
COFD(4813) = -2.07336791d+01
COFD(4814) = 5.23116678d+00
COFD(4815) = -4.54972675d-01
COFD(4816) = 1.93572354d-02
COFD(4817) = -1.81952014d+01
COFD(4818) = 4.40883885d+00
COFD(4819) = -3.58021243d-01
COFD(4820) = 1.55467750d-02
COFD(4821) = -1.67147225d+01
COFD(4822) = 3.84749999d+00
COFD(4823) = -2.87575310d-01
COFD(4824) = 1.25862768d-02
COFD(4825) = -2.19124011d+01
COFD(4826) = 5.47350698d+00
COFD(4827) = -4.56806529d-01
COFD(4828) = 1.82590273d-02
COFD(4829) = -2.19213014d+01
COFD(4830) = 5.47350698d+00
COFD(4831) = -4.56806529d-01
COFD(4832) = 1.82590273d-02
COFD(4833) = -1.88646205d+01
COFD(4834) = 4.72476764d+00
COFD(4835) = -3.96306836d-01
COFD(4836) = 1.70964541d-02
COFD(4837) = -2.07294312d+01
COFD(4838) = 5.23116678d+00
COFD(4839) = -4.54972675d-01
COFD(4840) = 1.93572354d-02
COFD(4841) = -2.07294312d+01
COFD(4842) = 5.23116678d+00
COFD(4843) = -4.54972675d-01
COFD(4844) = 1.93572354d-02
COFD(4845) = -2.07294312d+01
COFD(4846) = 5.23116678d+00
COFD(4847) = -4.54972675d-01
COFD(4848) = 1.93572354d-02
COFD(4849) = -2.07249293d+01
COFD(4850) = 5.23116678d+00
COFD(4851) = -4.54972675d-01
COFD(4852) = 1.93572354d-02
COFD(4853) = -1.76002031d+01
COFD(4854) = 4.19171952d+00
COFD(4855) = -3.31354810d-01
COFD(4856) = 1.44520623d-02
COFD(4857) = -1.87717330d+01
COFD(4858) = 4.61078776d+00
COFD(4859) = -3.82625667d-01
COFD(4860) = 1.65478601d-02
COFD(4861) = -2.14414783d+01
COFD(4862) = 5.37331605d+00
COFD(4863) = -4.70491203d-01
COFD(4864) = 1.99134666d-02
COFD(4865) = -2.14458019d+01
COFD(4866) = 5.37331605d+00
COFD(4867) = -4.70491203d-01
COFD(4868) = 1.99134666d-02
COFD(4869) = -2.21110290d+01
COFD(4870) = 5.58360799d+00
COFD(4871) = -4.82701436d-01
COFD(4872) = 1.98437922d-02
COFD(4873) = -2.21153597d+01
COFD(4874) = 5.58360799d+00
COFD(4875) = -4.82701436d-01
COFD(4876) = 1.98437922d-02
COFD(4877) = -1.34723215d+01
COFD(4878) = 3.09379603d+00
COFD(4879) = -1.91268635d-01
COFD(4880) = 8.47480224d-03
COFD(4881) = -1.72239172d+01
COFD(4882) = 4.69060745d+00
COFD(4883) = -3.92369888d-01
COFD(4884) = 1.69459661d-02
COFD(4885) = -1.65559787d+01
COFD(4886) = 3.95035840d+00
COFD(4887) = -3.00959418d-01
COFD(4888) = 1.31692593d-02
COFD(4889) = -1.78938745d+01
COFD(4890) = 4.29613154d+00
COFD(4891) = -3.44012526d-01
COFD(4892) = 1.49643715d-02
COFD(4893) = -1.65749533d+01
COFD(4894) = 3.95035840d+00
COFD(4895) = -3.00959418d-01
COFD(4896) = 1.31692593d-02
COFD(4897) = -2.15172770d+01
COFD(4898) = 5.46737673d+00
COFD(4899) = -4.55696085d-01
COFD(4900) = 1.81982625d-02
COFD(4901) = -1.79009181d+01
COFD(4902) = 4.29613154d+00
COFD(4903) = -3.44012526d-01
COFD(4904) = 1.49643715d-02
COFD(4905) = -1.79076361d+01
COFD(4906) = 4.29613154d+00
COFD(4907) = -3.44012526d-01
COFD(4908) = 1.49643715d-02
COFD(4909) = -1.61161138d+01
COFD(4910) = 3.75910622d+00
COFD(4911) = -2.75986578d-01
COFD(4912) = 1.20782843d-02
COFD(4913) = -1.64887871d+01
COFD(4914) = 3.95035840d+00
COFD(4915) = -3.00959418d-01
COFD(4916) = 1.31692593d-02
COFD(4917) = -1.88306747d+01
COFD(4918) = 4.68393046d+00
COFD(4919) = -3.91610863d-01
COFD(4920) = 1.69174645d-02
COFD(4921) = -1.88306747d+01
COFD(4922) = 4.68393046d+00
COFD(4923) = -3.91610863d-01
COFD(4924) = 1.69174645d-02
COFD(4925) = -1.88532497d+01
COFD(4926) = 4.68393046d+00
COFD(4927) = -3.91610863d-01
COFD(4928) = 1.69174645d-02
COFD(4929) = -1.87968848d+01
COFD(4930) = 4.66162351d+00
COFD(4931) = -3.88920477d-01
COFD(4932) = 1.68089648d-02
COFD(4933) = -1.76383156d+01
COFD(4934) = 4.19935698d+00
COFD(4935) = -3.32310212d-01
COFD(4936) = 1.44920670d-02
COFD(4937) = -2.09309753d+01
COFD(4938) = 5.30153901d+00
COFD(4939) = -4.63335119d-01
COFD(4940) = 1.96897053d-02
COFD(4941) = -2.20055544d+01
COFD(4942) = 5.55935694d+00
COFD(4943) = -4.74154740d-01
COFD(4944) = 1.92584304d-02
COFD(4945) = -2.20137178d+01
COFD(4946) = 5.55935694d+00
COFD(4947) = -4.74154740d-01
COFD(4948) = 1.92584304d-02
COFD(4949) = -2.19594078d+01
COFD(4950) = 5.59157589d+00
COFD(4951) = -4.85617912d-01
COFD(4952) = 2.00461138d-02
COFD(4953) = -2.19594078d+01
COFD(4954) = 5.59157589d+00
COFD(4955) = -4.85617912d-01
COFD(4956) = 2.00461138d-02
COFD(4957) = -2.20592197d+01
COFD(4958) = 5.58129885d+00
COFD(4959) = -4.78532921d-01
COFD(4960) = 1.95095699d-02
COFD(4961) = -2.03759451d+01
COFD(4962) = 5.13263469d+00
COFD(4963) = -4.44457285d-01
COFD(4964) = 1.89932102d-02
COFD(4965) = -2.03861000d+01
COFD(4966) = 5.13263469d+00
COFD(4967) = -4.44457285d-01
COFD(4968) = 1.89932102d-02
COFD(4969) = -2.03956853d+01
COFD(4970) = 5.13263469d+00
COFD(4971) = -4.44457285d-01
COFD(4972) = 1.89932102d-02
COFD(4973) = -2.12218960d+01
COFD(4974) = 5.39823225d+00
COFD(4975) = -4.72294645d-01
COFD(4976) = 1.99340225d-02
COFD(4977) = -2.10963515d+01
COFD(4978) = 5.34286099d+00
COFD(4979) = -4.68100992d-01
COFD(4980) = 1.98731399d-02
COFD(4981) = -2.11044965d+01
COFD(4982) = 5.34286099d+00
COFD(4983) = -4.68100992d-01
COFD(4984) = 1.98731399d-02
COFD(4985) = -1.88761387d+01
COFD(4986) = 4.72476764d+00
COFD(4987) = -3.96306836d-01
COFD(4988) = 1.70964541d-02
COFD(4989) = -2.21181720d+01
COFD(4990) = 5.58360799d+00
COFD(4991) = -4.82701436d-01
COFD(4992) = 1.98437922d-02
COFD(4993) = -2.21181720d+01
COFD(4994) = 5.58360799d+00
COFD(4995) = -4.82701436d-01
COFD(4996) = 1.98437922d-02
COFD(4997) = -1.61680488d+01
COFD(4998) = 3.75910622d+00
COFD(4999) = -2.75986578d-01
COFD(5000) = 1.20782843d-02
COFD(5001) = -1.65063757d+01
COFD(5002) = 3.95035840d+00
COFD(5003) = -3.00959418d-01
COFD(5004) = 1.31692593d-02
COFD(5005) = -1.65270223d+01
COFD(5006) = 3.95035840d+00
COFD(5007) = -3.00959418d-01
COFD(5008) = 1.31692593d-02
COFD(5009) = -2.16976525d+01
COFD(5010) = 5.58175603d+00
COFD(5011) = -4.78660918d-01
COFD(5012) = 1.95178500d-02
COFD(5013) = -1.65072489d+01
COFD(5014) = 3.75910622d+00
COFD(5015) = -2.75986578d-01
COFD(5016) = 1.20782843d-02
COFD(5017) = -1.76265380d+01
COFD(5018) = 4.19171952d+00
COFD(5019) = -3.31354810d-01
COFD(5020) = 1.44520623d-02
COFD(5021) = -2.02144469d+01
COFD(5022) = 5.08408346d+00
COFD(5023) = -4.38907391d-01
COFD(5024) = 1.87824789d-02
COFD(5025) = -2.07455261d+01
COFD(5026) = 5.23116678d+00
COFD(5027) = -4.54972675d-01
COFD(5028) = 1.93572354d-02
COFD(5029) = -1.82054352d+01
COFD(5030) = 4.40883885d+00
COFD(5031) = -3.58021243d-01
COFD(5032) = 1.55467750d-02
COFD(5033) = -1.67241238d+01
COFD(5034) = 3.84749999d+00
COFD(5035) = -2.87575310d-01
COFD(5036) = 1.25862768d-02
COFD(5037) = -2.19219828d+01
COFD(5038) = 5.47350698d+00
COFD(5039) = -4.56806529d-01
COFD(5040) = 1.82590273d-02
COFD(5041) = -2.19310570d+01
COFD(5042) = 5.47350698d+00
COFD(5043) = -4.56806529d-01
COFD(5044) = 1.82590273d-02
COFD(5045) = -1.88761525d+01
COFD(5046) = 4.72476764d+00
COFD(5047) = -3.96306836d-01
COFD(5048) = 1.70964541d-02
COFD(5049) = -2.07411769d+01
COFD(5050) = 5.23116678d+00
COFD(5051) = -4.54972675d-01
COFD(5052) = 1.93572354d-02
COFD(5053) = -2.07411769d+01
COFD(5054) = 5.23116678d+00
COFD(5055) = -4.54972675d-01
COFD(5056) = 1.93572354d-02
COFD(5057) = -2.07411769d+01
COFD(5058) = 5.23116678d+00
COFD(5059) = -4.54972675d-01
COFD(5060) = 1.93572354d-02
COFD(5061) = -2.07365685d+01
COFD(5062) = 5.23116678d+00
COFD(5063) = -4.54972675d-01
COFD(5064) = 1.93572354d-02
COFD(5065) = -1.76099552d+01
COFD(5066) = 4.19171952d+00
COFD(5067) = -3.31354810d-01
COFD(5068) = 1.44520623d-02
COFD(5069) = -1.87831434d+01
COFD(5070) = 4.61078776d+00
COFD(5071) = -3.82625667d-01
COFD(5072) = 1.65478601d-02
COFD(5073) = -2.14532306d+01
COFD(5074) = 5.37331605d+00
COFD(5075) = -4.70491203d-01
COFD(5076) = 1.99134666d-02
COFD(5077) = -2.14576575d+01
COFD(5078) = 5.37331605d+00
COFD(5079) = -4.70491203d-01
COFD(5080) = 1.99134666d-02
COFD(5081) = -2.21227768d+01
COFD(5082) = 5.58360799d+00
COFD(5083) = -4.82701436d-01
COFD(5084) = 1.98437922d-02
COFD(5085) = -2.21272109d+01
COFD(5086) = 5.58360799d+00
COFD(5087) = -4.82701436d-01
COFD(5088) = 1.98437922d-02
COFD(5089) = -1.42229194d+01
COFD(5090) = 3.38669384d+00
COFD(5091) = -2.28784122d-01
COFD(5092) = 1.00790953d-02
COFD(5093) = -1.82251914d+01
COFD(5094) = 5.05237312d+00
COFD(5095) = -4.35182396d-01
COFD(5096) = 1.86363074d-02
COFD(5097) = -1.74792112d+01
COFD(5098) = 4.29676909d+00
COFD(5099) = -3.44085306d-01
COFD(5100) = 1.49671135d-02
COFD(5101) = -1.89544778d+01
COFD(5102) = 4.68595732d+00
COFD(5103) = -3.91842840d-01
COFD(5104) = 1.69262542d-02
COFD(5105) = -1.74984476d+01
COFD(5106) = 4.29676909d+00
COFD(5107) = -3.44085306d-01
COFD(5108) = 1.49671135d-02
COFD(5109) = -2.08812333d+01
COFD(5110) = 5.08859217d+00
COFD(5111) = -3.90525428d-01
COFD(5112) = 1.47376395d-02
COFD(5113) = -1.89616623d+01
COFD(5114) = 4.68595732d+00
COFD(5115) = -3.91842840d-01
COFD(5116) = 1.69262542d-02
COFD(5117) = -1.89685165d+01
COFD(5118) = 4.68595732d+00
COFD(5119) = -3.91842840d-01
COFD(5120) = 1.69262542d-02
COFD(5121) = -1.71884218d+01
COFD(5122) = 4.17190426d+00
COFD(5123) = -3.28894681d-01
COFD(5124) = 1.43498101d-02
COFD(5125) = -1.74111692d+01
COFD(5126) = 4.29676909d+00
COFD(5127) = -3.44085306d-01
COFD(5128) = 1.49671135d-02
COFD(5129) = -1.98418115d+01
COFD(5130) = 5.04367502d+00
COFD(5131) = -4.34153325d-01
COFD(5132) = 1.85956055d-02
COFD(5133) = -1.98418115d+01
COFD(5134) = 5.04367502d+00
COFD(5135) = -4.34153325d-01
COFD(5136) = 1.85956055d-02
COFD(5137) = -1.98646734d+01
COFD(5138) = 5.04367502d+00
COFD(5139) = -4.34153325d-01
COFD(5140) = 1.85956055d-02
COFD(5141) = -1.98075055d+01
COFD(5142) = 5.02169524d+00
COFD(5143) = -4.31582804d-01
COFD(5144) = 1.84953568d-02
COFD(5145) = -1.86157761d+01
COFD(5146) = 4.55689508d+00
COFD(5147) = -3.75937921d-01
COFD(5148) = 1.62703488d-02
COFD(5149) = -2.16802612d+01
COFD(5150) = 5.52918296d+00
COFD(5151) = -4.85360709d-01
COFD(5152) = 2.03448006d-02
COFD(5153) = -2.16296373d+01
COFD(5154) = 5.29019717d+00
COFD(5155) = -4.24502606d-01
COFD(5156) = 1.65197343d-02
COFD(5157) = -2.16379567d+01
COFD(5158) = 5.29019717d+00
COFD(5159) = -4.24502606d-01
COFD(5160) = 1.65197343d-02
COFD(5161) = -2.19670848d+01
COFD(5162) = 5.48847873d+00
COFD(5163) = -4.59558930d-01
COFD(5164) = 1.84107961d-02
COFD(5165) = -2.19670848d+01
COFD(5166) = 5.48847873d+00
COFD(5167) = -4.59558930d-01
COFD(5168) = 1.84107961d-02
COFD(5169) = -2.17419207d+01
COFD(5170) = 5.33732875d+00
COFD(5171) = -4.32453425d-01
COFD(5172) = 1.69373833d-02
COFD(5173) = -2.12018018d+01
COFD(5174) = 5.39823225d+00
COFD(5175) = -4.72294645d-01
COFD(5176) = 1.99340225d-02
COFD(5177) = -2.12121370d+01
COFD(5178) = 5.39823225d+00
COFD(5179) = -4.72294645d-01
COFD(5180) = 1.99340225d-02
COFD(5181) = -2.12218960d+01
COFD(5182) = 5.39823225d+00
COFD(5183) = -4.72294645d-01
COFD(5184) = 1.99340225d-02
COFD(5185) = -2.19327397d+01
COFD(5186) = 5.60638188d+00
COFD(5187) = -4.91272522d-01
COFD(5188) = 2.04396264d-02
COFD(5189) = -2.18190539d+01
COFD(5190) = 5.55753905d+00
COFD(5191) = -4.88136714d-01
COFD(5192) = 2.04294957d-02
COFD(5193) = -2.18273547d+01
COFD(5194) = 5.55753905d+00
COFD(5195) = -4.88136714d-01
COFD(5196) = 2.04294957d-02
COFD(5197) = -1.99081556d+01
COFD(5198) = 5.09311649d+00
COFD(5199) = -4.39965178d-01
COFD(5200) = 1.88238537d-02
COFD(5201) = -2.20453723d+01
COFD(5202) = 5.44448440d+00
COFD(5203) = -4.51529024d-01
COFD(5204) = 1.79698119d-02
COFD(5205) = -2.20453723d+01
COFD(5206) = 5.44448440d+00
COFD(5207) = -4.51529024d-01
COFD(5208) = 1.79698119d-02
COFD(5209) = -1.72409687d+01
COFD(5210) = 4.17190426d+00
COFD(5211) = -3.28894681d-01
COFD(5212) = 1.43498101d-02
COFD(5213) = -1.74287716d+01
COFD(5214) = 4.29676909d+00
COFD(5215) = -3.44085306d-01
COFD(5216) = 1.49671135d-02
COFD(5217) = -1.74496921d+01
COFD(5218) = 4.29676909d+00
COFD(5219) = -3.44085306d-01
COFD(5220) = 1.49671135d-02
COFD(5221) = -2.13779772d+01
COFD(5222) = 5.33947196d+00
COFD(5223) = -4.32820858d-01
COFD(5224) = 1.69568580d-02
COFD(5225) = -1.75856336d+01
COFD(5226) = 4.17190426d+00
COFD(5227) = -3.28894681d-01
COFD(5228) = 1.43498101d-02
COFD(5229) = -1.86033113d+01
COFD(5230) = 4.54915847d+00
COFD(5231) = -3.75000738d-01
COFD(5232) = 1.62324821d-02
COFD(5233) = -2.10998106d+01
COFD(5234) = 5.37657006d+00
COFD(5235) = -4.70707921d-01
COFD(5236) = 1.99145326d-02
COFD(5237) = -2.15315009d+01
COFD(5238) = 5.47662534d+00
COFD(5239) = -4.79750189d-01
COFD(5240) = 2.01492627d-02
COFD(5241) = -1.92234721d+01
COFD(5242) = 4.77781979d+00
COFD(5243) = -4.02620572d-01
COFD(5244) = 1.73472054d-02
COFD(5245) = -1.77393195d+01
COFD(5246) = 4.23219154d+00
COFD(5247) = -3.36374352d-01
COFD(5248) = 1.46603433d-02
COFD(5249) = -2.13056183d+01
COFD(5250) = 5.10013279d+00
COFD(5251) = -3.92393844d-01
COFD(5252) = 1.48333781d-02
COFD(5253) = -2.13148599d+01
COFD(5254) = 5.10013279d+00
COFD(5255) = -3.92393844d-01
COFD(5256) = 1.48333781d-02
COFD(5257) = -1.99081697d+01
COFD(5258) = 5.09311649d+00
COFD(5259) = -4.39965178d-01
COFD(5260) = 1.88238537d-02
COFD(5261) = -2.15270531d+01
COFD(5262) = 5.47662534d+00
COFD(5263) = -4.79750189d-01
COFD(5264) = 2.01492627d-02
COFD(5265) = -2.15270531d+01
COFD(5266) = 5.47662534d+00
COFD(5267) = -4.79750189d-01
COFD(5268) = 2.01492627d-02
COFD(5269) = -2.15270531d+01
COFD(5270) = 5.47662534d+00
COFD(5271) = -4.79750189d-01
COFD(5272) = 2.01492627d-02
COFD(5273) = -2.15223412d+01
COFD(5274) = 5.47662534d+00
COFD(5275) = -4.79750189d-01
COFD(5276) = 2.01492627d-02
COFD(5277) = -1.85864144d+01
COFD(5278) = 4.54915847d+00
COFD(5279) = -3.75000738d-01
COFD(5280) = 1.62324821d-02
COFD(5281) = -1.98040322d+01
COFD(5282) = 4.97569695d+00
COFD(5283) = -4.26123307d-01
COFD(5284) = 1.82788664d-02
COFD(5285) = -2.21980988d+01
COFD(5286) = 5.59472344d+00
COFD(5287) = -4.91421518d-01
COFD(5288) = 2.05117088d-02
COFD(5289) = -2.22026260d+01
COFD(5290) = 5.59472344d+00
COFD(5291) = -4.91421518d-01
COFD(5292) = 2.05117088d-02
COFD(5293) = -2.20500806d+01
COFD(5294) = 5.44448440d+00
COFD(5295) = -4.51529024d-01
COFD(5296) = 1.79698119d-02
COFD(5297) = -2.20546151d+01
COFD(5298) = 5.44448440d+00
COFD(5299) = -4.51529024d-01
COFD(5300) = 1.79698119d-02
COFD(5301) = -1.39913897d+01
COFD(5302) = 3.26384506d+00
COFD(5303) = -2.12947087d-01
COFD(5304) = 9.39743888d-03
COFD(5305) = -1.79339327d+01
COFD(5306) = 4.91373893d+00
COFD(5307) = -4.18747629d-01
COFD(5308) = 1.79856610d-02
COFD(5309) = -1.72496634d+01
COFD(5310) = 4.17889917d+00
COFD(5311) = -3.29752510d-01
COFD(5312) = 1.43850275d-02
COFD(5313) = -1.86335932d+01
COFD(5314) = 4.53572533d+00
COFD(5315) = -3.73386925d-01
COFD(5316) = 1.61678881d-02
COFD(5317) = -1.72691500d+01
COFD(5318) = 4.17889917d+00
COFD(5319) = -3.29752510d-01
COFD(5320) = 1.43850275d-02
COFD(5321) = -2.12597312d+01
COFD(5322) = 5.24930667d+00
COFD(5323) = -4.17435088d-01
COFD(5324) = 1.61434424d-02
COFD(5325) = -1.86409139d+01
COFD(5326) = 4.53572533d+00
COFD(5327) = -3.73386925d-01
COFD(5328) = 1.61678881d-02
COFD(5329) = -1.86479000d+01
COFD(5330) = 4.53572533d+00
COFD(5331) = -3.73386925d-01
COFD(5332) = 1.61678881d-02
COFD(5333) = -1.69491115d+01
COFD(5334) = 4.05099737d+00
COFD(5335) = -3.13841660d-01
COFD(5336) = 1.37218854d-02
COFD(5337) = -1.71808106d+01
COFD(5338) = 4.17889917d+00
COFD(5339) = -3.29752510d-01
COFD(5340) = 1.43850275d-02
COFD(5341) = -1.95263312d+01
COFD(5342) = 4.90255048d+00
COFD(5343) = -4.17368501d-01
COFD(5344) = 1.79287358d-02
COFD(5345) = -1.95263312d+01
COFD(5346) = 4.90255048d+00
COFD(5347) = -4.17368501d-01
COFD(5348) = 1.79287358d-02
COFD(5349) = -1.95494668d+01
COFD(5350) = 4.90255048d+00
COFD(5351) = -4.17368501d-01
COFD(5352) = 1.79287358d-02
COFD(5353) = -1.94763688d+01
COFD(5354) = 4.87333294d+00
COFD(5355) = -4.13769241d-01
COFD(5356) = 1.77802244d-02
COFD(5357) = -1.83455435d+01
COFD(5358) = 4.42828044d+00
COFD(5359) = -3.60417833d-01
COFD(5360) = 1.56455103d-02
COFD(5361) = -2.14224484d+01
COFD(5362) = 5.41729961d+00
COFD(5363) = -4.73400281d-01
COFD(5364) = 1.99269567d-02
COFD(5365) = -2.19229190d+01
COFD(5366) = 5.41841631d+00
COFD(5367) = -4.46818971d-01
COFD(5368) = 1.77127652d-02
COFD(5369) = -2.19313890d+01
COFD(5370) = 5.41841631d+00
COFD(5371) = -4.46818971d-01
COFD(5372) = 1.77127652d-02
COFD(5373) = -2.21070030d+01
COFD(5374) = 5.55072945d+00
COFD(5375) = -4.72525345d-01
COFD(5376) = 1.91674202d-02
COFD(5377) = -2.21070030d+01
COFD(5378) = 5.55072945d+00
COFD(5379) = -4.72525345d-01
COFD(5380) = 1.91674202d-02
COFD(5381) = -2.20028184d+01
COFD(5382) = 5.45178028d+00
COFD(5383) = -4.52847771d-01
COFD(5384) = 1.80418544d-02
COFD(5385) = -2.10759163d+01
COFD(5386) = 5.34286099d+00
COFD(5387) = -4.68100992d-01
COFD(5388) = 1.98731399d-02
COFD(5389) = -2.10864251d+01
COFD(5390) = 5.34286099d+00
COFD(5391) = -4.68100992d-01
COFD(5392) = 1.98731399d-02
COFD(5393) = -2.10963515d+01
COFD(5394) = 5.34286099d+00
COFD(5395) = -4.68100992d-01
COFD(5396) = 1.98731399d-02
COFD(5397) = -2.18190539d+01
COFD(5398) = 5.55753905d+00
COFD(5399) = -4.88136714d-01
COFD(5400) = 2.04294957d-02
COFD(5401) = -2.15575659d+01
COFD(5402) = 5.44803850d+00
COFD(5403) = -4.76610560d-01
COFD(5404) = 2.00355294d-02
COFD(5405) = -2.15660171d+01
COFD(5406) = 5.44803850d+00
COFD(5407) = -4.76610560d-01
COFD(5408) = 2.00355294d-02
COFD(5409) = -1.96309042d+01
COFD(5410) = 4.95923807d+00
COFD(5411) = -4.24176182d-01
COFD(5412) = 1.82020215d-02
COFD(5413) = -2.22059540d+01
COFD(5414) = 5.51722375d+00
COFD(5415) = -4.66081431d-01
COFD(5416) = 1.88044011d-02
COFD(5417) = -2.22059540d+01
COFD(5418) = 5.51722375d+00
COFD(5419) = -4.66081431d-01
COFD(5420) = 1.88044011d-02
COFD(5421) = -1.70022408d+01
COFD(5422) = 4.05099737d+00
COFD(5423) = -3.13841660d-01
COFD(5424) = 1.37218854d-02
COFD(5425) = -1.72003855d+01
COFD(5426) = 4.17889917d+00
COFD(5427) = -3.29752510d-01
COFD(5428) = 1.43850275d-02
COFD(5429) = -1.72215676d+01
COFD(5430) = 4.17889917d+00
COFD(5431) = -3.29752510d-01
COFD(5432) = 1.43850275d-02
COFD(5433) = -2.16435160d+01
COFD(5434) = 5.45344302d+00
COFD(5435) = -4.53149740d-01
COFD(5436) = 1.80583906d-02
COFD(5437) = -1.73443798d+01
COFD(5438) = 4.05099737d+00
COFD(5439) = -3.13841660d-01
COFD(5440) = 1.37218854d-02
COFD(5441) = -1.83338353d+01
COFD(5442) = 4.42045763d+00
COFD(5443) = -3.59451578d-01
COFD(5444) = 1.56056164d-02
COFD(5445) = -2.08992622d+01
COFD(5446) = 5.28495077d+00
COFD(5447) = -4.61326381d-01
COFD(5448) = 1.96080621d-02
COFD(5449) = -2.13547906d+01
COFD(5450) = 5.39784038d+00
COFD(5451) = -4.72269379d-01
COFD(5452) = 1.99339592d-02
COFD(5453) = -1.89704705d+01
COFD(5454) = 4.65726073d+00
COFD(5455) = -3.88398872d-01
COFD(5456) = 1.67881413d-02
COFD(5457) = -1.74955909d+01
COFD(5458) = 4.11180810d+00
COFD(5459) = -3.21531431d-01
COFD(5460) = 1.40478729d-02
COFD(5461) = -2.16618765d+01
COFD(5462) = 5.25608421d+00
COFD(5463) = -4.18625448d-01
COFD(5464) = 1.62072845d-02
COFD(5465) = -2.16712796d+01
COFD(5466) = 5.25608421d+00
COFD(5467) = -4.18625448d-01
COFD(5468) = 1.62072845d-02
COFD(5469) = -1.96309186d+01
COFD(5470) = 4.95923807d+00
COFD(5471) = -4.24176182d-01
COFD(5472) = 1.82020215d-02
COFD(5473) = -2.13502469d+01
COFD(5474) = 5.39784038d+00
COFD(5475) = -4.72269379d-01
COFD(5476) = 1.99339592d-02
COFD(5477) = -2.13502469d+01
COFD(5478) = 5.39784038d+00
COFD(5479) = -4.72269379d-01
COFD(5480) = 1.99339592d-02
COFD(5481) = -2.13502469d+01
COFD(5482) = 5.39784038d+00
COFD(5483) = -4.72269379d-01
COFD(5484) = 1.99339592d-02
COFD(5485) = -2.13454345d+01
COFD(5486) = 5.39784038d+00
COFD(5487) = -4.72269379d-01
COFD(5488) = 1.99339592d-02
COFD(5489) = -1.83166353d+01
COFD(5490) = 4.42045763d+00
COFD(5491) = -3.59451578d-01
COFD(5492) = 1.56056164d-02
COFD(5493) = -1.94928815d+01
COFD(5494) = 4.83189721d+00
COFD(5495) = -4.08932249d-01
COFD(5496) = 1.75924650d-02
COFD(5497) = -2.19981673d+01
COFD(5498) = 5.51276597d+00
COFD(5499) = -4.83701824d-01
COFD(5500) = 2.02915297d-02
COFD(5501) = -2.20027922d+01
COFD(5502) = 5.51276597d+00
COFD(5503) = -4.83701824d-01
COFD(5504) = 2.02915297d-02
COFD(5505) = -2.22107627d+01
COFD(5506) = 5.51722375d+00
COFD(5507) = -4.66081431d-01
COFD(5508) = 1.88044011d-02
COFD(5509) = -2.22153950d+01
COFD(5510) = 5.51722375d+00
COFD(5511) = -4.66081431d-01
COFD(5512) = 1.88044011d-02
COFD(5513) = -1.39924781d+01
COFD(5514) = 3.26384506d+00
COFD(5515) = -2.12947087d-01
COFD(5516) = 9.39743888d-03
COFD(5517) = -1.79344949d+01
COFD(5518) = 4.91373893d+00
COFD(5519) = -4.18747629d-01
COFD(5520) = 1.79856610d-02
COFD(5521) = -1.72556499d+01
COFD(5522) = 4.17889917d+00
COFD(5523) = -3.29752510d-01
COFD(5524) = 1.43850275d-02
COFD(5525) = -1.86424545d+01
COFD(5526) = 4.53572533d+00
COFD(5527) = -3.73386925d-01
COFD(5528) = 1.61678881d-02
COFD(5529) = -1.72753760d+01
COFD(5530) = 4.17889917d+00
COFD(5531) = -3.29752510d-01
COFD(5532) = 1.43850275d-02
COFD(5533) = -2.12661865d+01
COFD(5534) = 5.24930667d+00
COFD(5535) = -4.17435088d-01
COFD(5536) = 1.61434424d-02
COFD(5537) = -1.86499071d+01
COFD(5538) = 4.53572533d+00
COFD(5539) = -3.73386925d-01
COFD(5540) = 1.61678881d-02
COFD(5541) = -1.86570209d+01
COFD(5542) = 4.53572533d+00
COFD(5543) = -3.73386925d-01
COFD(5544) = 1.61678881d-02
COFD(5545) = -1.69540369d+01
COFD(5546) = 4.05099737d+00
COFD(5547) = -3.13841660d-01
COFD(5548) = 1.37218854d-02
COFD(5549) = -1.71860230d+01
COFD(5550) = 4.17889917d+00
COFD(5551) = -3.29752510d-01
COFD(5552) = 1.43850275d-02
COFD(5553) = -1.95318173d+01
COFD(5554) = 4.90255048d+00
COFD(5555) = -4.17368501d-01
COFD(5556) = 1.79287358d-02
COFD(5557) = -1.95318173d+01
COFD(5558) = 4.90255048d+00
COFD(5559) = -4.17368501d-01
COFD(5560) = 1.79287358d-02
COFD(5561) = -1.95552142d+01
COFD(5562) = 4.90255048d+00
COFD(5563) = -4.17368501d-01
COFD(5564) = 1.79287358d-02
COFD(5565) = -1.94823660d+01
COFD(5566) = 4.87333294d+00
COFD(5567) = -4.13769241d-01
COFD(5568) = 1.77802244d-02
COFD(5569) = -1.83538377d+01
COFD(5570) = 4.42828044d+00
COFD(5571) = -3.60417833d-01
COFD(5572) = 1.56455103d-02
COFD(5573) = -2.14326461d+01
COFD(5574) = 5.41729961d+00
COFD(5575) = -4.73400281d-01
COFD(5576) = 1.99269567d-02
COFD(5577) = -2.19313638d+01
COFD(5578) = 5.41841631d+00
COFD(5579) = -4.46818971d-01
COFD(5580) = 1.77127652d-02
COFD(5581) = -2.19399793d+01
COFD(5582) = 5.41841631d+00
COFD(5583) = -4.46818971d-01
COFD(5584) = 1.77127652d-02
COFD(5585) = -2.21157340d+01
COFD(5586) = 5.55072945d+00
COFD(5587) = -4.72525345d-01
COFD(5588) = 1.91674202d-02
COFD(5589) = -2.21157340d+01
COFD(5590) = 5.55072945d+00
COFD(5591) = -4.72525345d-01
COFD(5592) = 1.91674202d-02
COFD(5593) = -2.20116855d+01
COFD(5594) = 5.45178028d+00
COFD(5595) = -4.52847771d-01
COFD(5596) = 1.80418544d-02
COFD(5597) = -2.10837327d+01
COFD(5598) = 5.34286099d+00
COFD(5599) = -4.68100992d-01
COFD(5600) = 1.98731399d-02
COFD(5601) = -2.10944088d+01
COFD(5602) = 5.34286099d+00
COFD(5603) = -4.68100992d-01
COFD(5604) = 1.98731399d-02
COFD(5605) = -2.11044965d+01
COFD(5606) = 5.34286099d+00
COFD(5607) = -4.68100992d-01
COFD(5608) = 1.98731399d-02
COFD(5609) = -2.18273547d+01
COFD(5610) = 5.55753905d+00
COFD(5611) = -4.88136714d-01
COFD(5612) = 2.04294957d-02
COFD(5613) = -2.15660171d+01
COFD(5614) = 5.44803850d+00
COFD(5615) = -4.76610560d-01
COFD(5616) = 2.00355294d-02
COFD(5617) = -2.15746136d+01
COFD(5618) = 5.44803850d+00
COFD(5619) = -4.76610560d-01
COFD(5620) = 2.00355294d-02
COFD(5621) = -1.96408127d+01
COFD(5622) = 4.95923807d+00
COFD(5623) = -4.24176182d-01
COFD(5624) = 1.82020215d-02
COFD(5625) = -2.22159630d+01
COFD(5626) = 5.51722375d+00
COFD(5627) = -4.66081431d-01
COFD(5628) = 1.88044011d-02
COFD(5629) = -2.22159630d+01
COFD(5630) = 5.51722375d+00
COFD(5631) = -4.66081431d-01
COFD(5632) = 1.88044011d-02
COFD(5633) = -1.70077215d+01
COFD(5634) = 4.05099737d+00
COFD(5635) = -3.13841660d-01
COFD(5636) = 1.37218854d-02
COFD(5637) = -1.72061277d+01
COFD(5638) = 4.17889917d+00
COFD(5639) = -3.29752510d-01
COFD(5640) = 1.43850275d-02
COFD(5641) = -1.72275598d+01
COFD(5642) = 4.17889917d+00
COFD(5643) = -3.29752510d-01
COFD(5644) = 1.43850275d-02
COFD(5645) = -2.16497474d+01
COFD(5646) = 5.45344302d+00
COFD(5647) = -4.53149740d-01
COFD(5648) = 1.80583906d-02
COFD(5649) = -1.73528250d+01
COFD(5650) = 4.05099737d+00
COFD(5651) = -3.13841660d-01
COFD(5652) = 1.37218854d-02
COFD(5653) = -1.83424227d+01
COFD(5654) = 4.42045763d+00
COFD(5655) = -3.59451578d-01
COFD(5656) = 1.56056164d-02
COFD(5657) = -2.09096408d+01
COFD(5658) = 5.28495077d+00
COFD(5659) = -4.61326381d-01
COFD(5660) = 1.96080621d-02
COFD(5661) = -2.13649886d+01
COFD(5662) = 5.39784038d+00
COFD(5663) = -4.72269379d-01
COFD(5664) = 1.99339592d-02
COFD(5665) = -1.89791987d+01
COFD(5666) = 4.65726073d+00
COFD(5667) = -3.88398872d-01
COFD(5668) = 1.67881413d-02
COFD(5669) = -1.75035712d+01
COFD(5670) = 4.11180810d+00
COFD(5671) = -3.21531431d-01
COFD(5672) = 1.40478729d-02
COFD(5673) = -2.16700183d+01
COFD(5674) = 5.25608421d+00
COFD(5675) = -4.18625448d-01
COFD(5676) = 1.62072845d-02
COFD(5677) = -2.16795773d+01
COFD(5678) = 5.25608421d+00
COFD(5679) = -4.18625448d-01
COFD(5680) = 1.62072845d-02
COFD(5681) = -1.96408274d+01
COFD(5682) = 4.95923807d+00
COFD(5683) = -4.24176182d-01
COFD(5684) = 1.82020215d-02
COFD(5685) = -2.13603517d+01
COFD(5686) = 5.39784038d+00
COFD(5687) = -4.72269379d-01
COFD(5688) = 1.99339592d-02
COFD(5689) = -2.13603517d+01
COFD(5690) = 5.39784038d+00
COFD(5691) = -4.72269379d-01
COFD(5692) = 1.99339592d-02
COFD(5693) = -2.13603517d+01
COFD(5694) = 5.39784038d+00
COFD(5695) = -4.72269379d-01
COFD(5696) = 1.99339592d-02
COFD(5697) = -2.13554415d+01
COFD(5698) = 5.39784038d+00
COFD(5699) = -4.72269379d-01
COFD(5700) = 1.99339592d-02
COFD(5701) = -1.83249299d+01
COFD(5702) = 4.42045763d+00
COFD(5703) = -3.59451578d-01
COFD(5704) = 1.56056164d-02
COFD(5705) = -1.95026789d+01
COFD(5706) = 4.83189721d+00
COFD(5707) = -4.08932249d-01
COFD(5708) = 1.75924650d-02
COFD(5709) = -2.20082782d+01
COFD(5710) = 5.51276597d+00
COFD(5711) = -4.83701824d-01
COFD(5712) = 2.02915297d-02
COFD(5713) = -2.20129980d+01
COFD(5714) = 5.51276597d+00
COFD(5715) = -4.83701824d-01
COFD(5716) = 2.02915297d-02
COFD(5717) = -2.22208695d+01
COFD(5718) = 5.51722375d+00
COFD(5719) = -4.66081431d-01
COFD(5720) = 1.88044011d-02
COFD(5721) = -2.22255968d+01
COFD(5722) = 5.51722375d+00
COFD(5723) = -4.66081431d-01
COFD(5724) = 1.88044011d-02
COFD(5725) = -1.22004324d+01
COFD(5726) = 2.80725489d+00
COFD(5727) = -1.54291406d-01
COFD(5728) = 6.88290911d-03
COFD(5729) = -1.54460820d+01
COFD(5730) = 4.26819983d+00
COFD(5731) = -3.40766379d-01
COFD(5732) = 1.48393361d-02
COFD(5733) = -1.49500357d+01
COFD(5734) = 3.52327209d+00
COFD(5735) = -2.46286208d-01
COFD(5736) = 1.08285963d-02
COFD(5737) = -1.64169433d+01
COFD(5738) = 3.89309916d+00
COFD(5739) = -2.93528188d-01
COFD(5740) = 1.28463177d-02
COFD(5741) = -1.49718233d+01
COFD(5742) = 3.52327209d+00
COFD(5743) = -2.46286208d-01
COFD(5744) = 1.08285963d-02
COFD(5745) = -2.10440675d+01
COFD(5746) = 5.59806282d+00
COFD(5747) = -4.87109535d-01
COFD(5748) = 2.01370226d-02
COFD(5749) = -1.64255964d+01
COFD(5750) = 3.89309916d+00
COFD(5751) = -2.93528188d-01
COFD(5752) = 1.28463177d-02
COFD(5753) = -1.64338757d+01
COFD(5754) = 3.89309916d+00
COFD(5755) = -2.93528188d-01
COFD(5756) = 1.28463177d-02
COFD(5757) = -1.46907028d+01
COFD(5758) = 3.39229020d+00
COFD(5759) = -2.29520232d-01
COFD(5760) = 1.01114311d-02
COFD(5761) = -1.48738066d+01
COFD(5762) = 3.52327209d+00
COFD(5763) = -2.46286208d-01
COFD(5764) = 1.08285963d-02
COFD(5765) = -1.72572042d+01
COFD(5766) = 4.26063341d+00
COFD(5767) = -3.39848064d-01
COFD(5768) = 1.48021313d-02
COFD(5769) = -1.72572042d+01
COFD(5770) = 4.26063341d+00
COFD(5771) = -3.39848064d-01
COFD(5772) = 1.48021313d-02
COFD(5773) = -1.72828302d+01
COFD(5774) = 4.26063341d+00
COFD(5775) = -3.39848064d-01
COFD(5776) = 1.48021313d-02
COFD(5777) = -1.72316148d+01
COFD(5778) = 4.24011069d+00
COFD(5779) = -3.37339810d-01
COFD(5780) = 1.46996679d-02
COFD(5781) = -1.60261675d+01
COFD(5782) = 3.73312045d+00
COFD(5783) = -2.72579779d-01
COFD(5784) = 1.19290272d-02
COFD(5785) = -1.94485982d+01
COFD(5786) = 4.91446566d+00
COFD(5787) = -4.18837152d-01
COFD(5788) = 1.79893537d-02
COFD(5789) = -2.14204185d+01
COFD(5790) = 5.59268435d+00
COFD(5791) = -4.91232974d-01
COFD(5792) = 2.05064746d-02
COFD(5793) = -2.14303479d+01
COFD(5794) = 5.59268435d+00
COFD(5795) = -4.91232974d-01
COFD(5796) = 2.05064746d-02
COFD(5797) = -2.09241647d+01
COFD(5798) = 5.42316225d+00
COFD(5799) = -4.73702801d-01
COFD(5800) = 1.99217718d-02
COFD(5801) = -2.09241647d+01
COFD(5802) = 5.42316225d+00
COFD(5803) = -4.73702801d-01
COFD(5804) = 1.99217718d-02
COFD(5805) = -2.13796303d+01
COFD(5806) = 5.56978987d+00
COFD(5807) = -4.89141980d-01
COFD(5808) = 2.04499210d-02
COFD(5809) = -1.88524385d+01
COFD(5810) = 4.72476764d+00
COFD(5811) = -3.96306836d-01
COFD(5812) = 1.70964541d-02
COFD(5813) = -1.88646070d+01
COFD(5814) = 4.72476764d+00
COFD(5815) = -3.96306836d-01
COFD(5816) = 1.70964541d-02
COFD(5817) = -1.88761387d+01
COFD(5818) = 4.72476764d+00
COFD(5819) = -3.96306836d-01
COFD(5820) = 1.70964541d-02
COFD(5821) = -1.99081556d+01
COFD(5822) = 5.09311649d+00
COFD(5823) = -4.39965178d-01
COFD(5824) = 1.88238537d-02
COFD(5825) = -1.96309042d+01
COFD(5826) = 4.95923807d+00
COFD(5827) = -4.24176182d-01
COFD(5828) = 1.82020215d-02
COFD(5829) = -1.96408127d+01
COFD(5830) = 4.95923807d+00
COFD(5831) = -4.24176182d-01
COFD(5832) = 1.82020215d-02
COFD(5833) = -1.72414862d+01
COFD(5834) = 4.29808578d+00
COFD(5835) = -3.44235570d-01
COFD(5836) = 1.49727727d-02
COFD(5837) = -2.12621914d+01
COFD(5838) = 5.47935225d+00
COFD(5839) = -4.80056796d-01
COFD(5840) = 2.01607180d-02
COFD(5841) = -2.12621914d+01
COFD(5842) = 5.47935225d+00
COFD(5843) = -4.80056796d-01
COFD(5844) = 2.01607180d-02
COFD(5845) = -1.47490868d+01
COFD(5846) = 3.39229020d+00
COFD(5847) = -2.29520232d-01
COFD(5848) = 1.01114311d-02
COFD(5849) = -1.48885202d+01
COFD(5850) = 3.52327209d+00
COFD(5851) = -2.46286208d-01
COFD(5852) = 1.08285963d-02
COFD(5853) = -1.49120950d+01
COFD(5854) = 3.52327209d+00
COFD(5855) = -2.46286208d-01
COFD(5856) = 1.08285963d-02
COFD(5857) = -2.09300683d+01
COFD(5858) = 5.56884637d+00
COFD(5859) = -4.89068748d-01
COFD(5860) = 2.04486922d-02
COFD(5861) = -1.51581584d+01
COFD(5862) = 3.39229020d+00
COFD(5863) = -2.29520232d-01
COFD(5864) = 1.01114311d-02
COFD(5865) = -1.60085628d+01
COFD(5866) = 3.72220402d+00
COFD(5867) = -2.71150591d-01
COFD(5868) = 1.18665265d-02
COFD(5869) = -1.86961076d+01
COFD(5870) = 4.67999453d+00
COFD(5871) = -3.91135253d-01
COFD(5872) = 1.68982388d-02
COFD(5873) = -1.92693772d+01
COFD(5874) = 4.84335373d+00
COFD(5875) = -4.10212751d-01
COFD(5876) = 1.76395874d-02
COFD(5877) = -1.67419735d+01
COFD(5878) = 4.01496916d+00
COFD(5879) = -3.09274229d-01
COFD(5880) = 1.35278248d-02
COFD(5881) = -1.53109350d+01
COFD(5882) = 3.45550279d+00
COFD(5883) = -2.37693211d-01
COFD(5884) = 1.04644205d-02
COFD(5885) = -2.15423744d+01
COFD(5886) = 5.60052862d+00
COFD(5887) = -4.87676563d-01
COFD(5888) = 2.01715067d-02
COFD(5889) = -2.15533322d+01
COFD(5890) = 5.60052862d+00
COFD(5891) = -4.87676563d-01
COFD(5892) = 2.01715067d-02
COFD(5893) = -1.72415036d+01
COFD(5894) = 4.29808578d+00
COFD(5895) = -3.44235570d-01
COFD(5896) = 1.49727727d-02
COFD(5897) = -1.92638705d+01
COFD(5898) = 4.84335373d+00
COFD(5899) = -4.10212751d-01
COFD(5900) = 1.76395874d-02
COFD(5901) = -1.92638705d+01
COFD(5902) = 4.84335373d+00
COFD(5903) = -4.10212751d-01
COFD(5904) = 1.76395874d-02
COFD(5905) = -1.92638705d+01
COFD(5906) = 4.84335373d+00
COFD(5907) = -4.10212751d-01
COFD(5908) = 1.76395874d-02
COFD(5909) = -1.92580496d+01
COFD(5910) = 4.84335373d+00
COFD(5911) = -4.10212751d-01
COFD(5912) = 1.76395874d-02
COFD(5913) = -1.59884305d+01
COFD(5914) = 3.72220402d+00
COFD(5915) = -2.71150591d-01
COFD(5916) = 1.18665265d-02
COFD(5917) = -1.72570607d+01
COFD(5918) = 4.19757624d+00
COFD(5919) = -3.32087529d-01
COFD(5920) = 1.44827462d-02
COFD(5921) = -2.01341439d+01
COFD(5922) = 5.03101171d+00
COFD(5923) = -4.32665019d-01
COFD(5924) = 1.85372086d-02
COFD(5925) = -2.01397498d+01
COFD(5926) = 5.03101171d+00
COFD(5927) = -4.32665019d-01
COFD(5928) = 1.85372086d-02
COFD(5929) = -2.12680082d+01
COFD(5930) = 5.47935225d+00
COFD(5931) = -4.80056796d-01
COFD(5932) = 2.01607180d-02
COFD(5933) = -2.12736225d+01
COFD(5934) = 5.47935225d+00
COFD(5935) = -4.80056796d-01
COFD(5936) = 2.01607180d-02
COFD(5937) = -1.57034851d+01
COFD(5938) = 3.93614244d+00
COFD(5939) = -2.99111497d-01
COFD(5940) = 1.30888229d-02
COFD(5941) = -1.94688688d+01
COFD(5942) = 5.43830787d+00
COFD(5943) = -4.75472880d-01
COFD(5944) = 1.99909996d-02
COFD(5945) = -1.90883268d+01
COFD(5946) = 4.84384483d+00
COFD(5947) = -4.10265575d-01
COFD(5948) = 1.76414287d-02
COFD(5949) = -2.05184870d+01
COFD(5950) = 5.18417470d+00
COFD(5951) = -4.49491573d-01
COFD(5952) = 1.91438508d-02
COFD(5953) = -1.91102652d+01
COFD(5954) = 4.84384483d+00
COFD(5955) = -4.10265575d-01
COFD(5956) = 1.76414287d-02
COFD(5957) = -1.87383952d+01
COFD(5958) = 3.96926341d+00
COFD(5959) = -2.16412264d-01
COFD(5960) = 6.06012078d-03
COFD(5961) = -2.05272328d+01
COFD(5962) = 5.18417470d+00
COFD(5963) = -4.49491573d-01
COFD(5964) = 1.91438508d-02
COFD(5965) = -2.05356023d+01
COFD(5966) = 5.18417470d+00
COFD(5967) = -4.49491573d-01
COFD(5968) = 1.91438508d-02
COFD(5969) = -1.87688110d+01
COFD(5970) = 4.71729964d+00
COFD(5971) = -3.95432573d-01
COFD(5972) = 1.70623691d-02
COFD(5973) = -1.90116191d+01
COFD(5974) = 4.84384483d+00
COFD(5975) = -4.10265575d-01
COFD(5976) = 1.76414287d-02
COFD(5977) = -2.11349086d+01
COFD(5978) = 5.42846112d+00
COFD(5979) = -4.74321870d-01
COFD(5980) = 1.99459749d-02
COFD(5981) = -2.11349086d+01
COFD(5982) = 5.42846112d+00
COFD(5983) = -4.74321870d-01
COFD(5984) = 1.99459749d-02
COFD(5985) = -2.11606963d+01
COFD(5986) = 5.42846112d+00
COFD(5987) = -4.74321870d-01
COFD(5988) = 1.99459749d-02
COFD(5989) = -2.11309207d+01
COFD(5990) = 5.41773516d+00
COFD(5991) = -4.73414338d-01
COFD(5992) = 1.99258685d-02
COFD(5993) = -2.02922701d+01
COFD(5994) = 5.11106992d+00
COFD(5995) = -4.42047129d-01
COFD(5996) = 1.89042990d-02
COFD(5997) = -2.22116706d+01
COFD(5998) = 5.54251230d+00
COFD(5999) = -4.70946314d-01
COFD(6000) = 1.90785869d-02
COFD(6001) = -2.00915040d+01
COFD(6002) = 4.41511629d+00
COFD(6003) = -2.84086963d-01
COFD(6004) = 9.37586971d-03
COFD(6005) = -2.01015340d+01
COFD(6006) = 4.41511629d+00
COFD(6007) = -2.84086963d-01
COFD(6008) = 9.37586971d-03
COFD(6009) = -2.09222454d+01
COFD(6010) = 4.82184721d+00
COFD(6011) = -3.48128875d-01
COFD(6012) = 1.25918978d-02
COFD(6013) = -2.09222454d+01
COFD(6014) = 4.82184721d+00
COFD(6015) = -3.48128875d-01
COFD(6016) = 1.25918978d-02
COFD(6017) = -2.03036402d+01
COFD(6018) = 4.50250781d+00
COFD(6019) = -2.97622106d-01
COFD(6020) = 1.00481473d-02
COFD(6021) = -2.20942491d+01
COFD(6022) = 5.58360799d+00
COFD(6023) = -4.82701436d-01
COFD(6024) = 1.98437922d-02
COFD(6025) = -2.21065306d+01
COFD(6026) = 5.58360799d+00
COFD(6027) = -4.82701436d-01
COFD(6028) = 1.98437922d-02
COFD(6029) = -2.21181720d+01
COFD(6030) = 5.58360799d+00
COFD(6031) = -4.82701436d-01
COFD(6032) = 1.98437922d-02
COFD(6033) = -2.20453723d+01
COFD(6034) = 5.44448440d+00
COFD(6035) = -4.51529024d-01
COFD(6036) = 1.79698119d-02
COFD(6037) = -2.22059540d+01
COFD(6038) = 5.51722375d+00
COFD(6039) = -4.66081431d-01
COFD(6040) = 1.88044011d-02
COFD(6041) = -2.22159630d+01
COFD(6042) = 5.51722375d+00
COFD(6043) = -4.66081431d-01
COFD(6044) = 1.88044011d-02
COFD(6045) = -2.12621914d+01
COFD(6046) = 5.47935225d+00
COFD(6047) = -4.80056796d-01
COFD(6048) = 2.01607180d-02
COFD(6049) = -2.09002742d+01
COFD(6050) = 4.72895031d+00
COFD(6051) = -3.33332771d-01
COFD(6052) = 1.18431478d-02
COFD(6053) = -2.09002742d+01
COFD(6054) = 4.72895031d+00
COFD(6055) = -3.33332771d-01
COFD(6056) = 1.18431478d-02
COFD(6057) = -1.88275332d+01
COFD(6058) = 4.71729964d+00
COFD(6059) = -3.95432573d-01
COFD(6060) = 1.70623691d-02
COFD(6061) = -1.90351359d+01
COFD(6062) = 4.84384483d+00
COFD(6063) = -4.10265575d-01
COFD(6064) = 1.76414287d-02
COFD(6065) = -1.90588668d+01
COFD(6066) = 4.84384483d+00
COFD(6067) = -4.10265575d-01
COFD(6068) = 1.76414287d-02
COFD(6069) = -1.99140603d+01
COFD(6070) = 4.50672715d+00
COFD(6071) = -2.98278922d-01
COFD(6072) = 1.00808824d-02
COFD(6073) = -1.92061609d+01
COFD(6074) = 4.71729964d+00
COFD(6075) = -3.95432573d-01
COFD(6076) = 1.70623691d-02
COFD(6077) = -2.02849953d+01
COFD(6078) = 5.10426133d+00
COFD(6079) = -4.41256919d-01
COFD(6080) = 1.88737290d-02
COFD(6081) = -2.20848247d+01
COFD(6082) = 5.59135292d+00
COFD(6083) = -4.85565630d-01
COFD(6084) = 2.00429035d-02
COFD(6085) = -2.22703814d+01
COFD(6086) = 5.57944619d+00
COFD(6087) = -4.78032626d-01
COFD(6088) = 1.94775370d-02
COFD(6089) = -2.08481163d+01
COFD(6090) = 5.29805436d+00
COFD(6091) = -4.62913371d-01
COFD(6092) = 1.96725621d-02
COFD(6093) = -1.93551579d+01
COFD(6094) = 4.77546552d+00
COFD(6095) = -4.02345606d-01
COFD(6096) = 1.73365265d-02
COFD(6097) = -1.92068539d+01
COFD(6098) = 3.99103591d+00
COFD(6099) = -2.19677673d-01
COFD(6100) = 6.21875426d-03
COFD(6101) = -1.92179182d+01
COFD(6102) = 3.99103591d+00
COFD(6103) = -2.19677673d-01
COFD(6104) = 6.21875426d-03
COFD(6105) = -2.12622090d+01
COFD(6106) = 5.47935225d+00
COFD(6107) = -4.80056796d-01
COFD(6108) = 2.01607180d-02
COFD(6109) = -2.22648059d+01
COFD(6110) = 5.57944619d+00
COFD(6111) = -4.78032626d-01
COFD(6112) = 1.94775370d-02
COFD(6113) = -2.22648059d+01
COFD(6114) = 5.57944619d+00
COFD(6115) = -4.78032626d-01
COFD(6116) = 1.94775370d-02
COFD(6117) = -2.22648059d+01
COFD(6118) = 5.57944619d+00
COFD(6119) = -4.78032626d-01
COFD(6120) = 1.94775370d-02
COFD(6121) = -2.22589130d+01
COFD(6122) = 5.57944619d+00
COFD(6123) = -4.78032626d-01
COFD(6124) = 1.94775370d-02
COFD(6125) = -2.02646611d+01
COFD(6126) = 5.10426133d+00
COFD(6127) = -4.41256919d-01
COFD(6128) = 1.88737290d-02
COFD(6129) = -2.12450965d+01
COFD(6130) = 5.40444222d+00
COFD(6131) = -4.72708609d-01
COFD(6132) = 1.99362392d-02
COFD(6133) = -2.24674263d+01
COFD(6134) = 5.49330641d+00
COFD(6135) = -4.60498247d-01
COFD(6136) = 1.84639199d-02
COFD(6137) = -2.24731023d+01
COFD(6138) = 5.49330641d+00
COFD(6139) = -4.60498247d-01
COFD(6140) = 1.84639199d-02
COFD(6141) = -2.09061629d+01
COFD(6142) = 4.72895031d+00
COFD(6143) = -3.33332771d-01
COFD(6144) = 1.18431478d-02
COFD(6145) = -2.09118474d+01
COFD(6146) = 4.72895031d+00
COFD(6147) = -3.33332771d-01
COFD(6148) = 1.18431478d-02
COFD(6149) = -1.57034851d+01
COFD(6150) = 3.93614244d+00
COFD(6151) = -2.99111497d-01
COFD(6152) = 1.30888229d-02
COFD(6153) = -1.94688688d+01
COFD(6154) = 5.43830787d+00
COFD(6155) = -4.75472880d-01
COFD(6156) = 1.99909996d-02
COFD(6157) = -1.90883268d+01
COFD(6158) = 4.84384483d+00
COFD(6159) = -4.10265575d-01
COFD(6160) = 1.76414287d-02
COFD(6161) = -2.05184870d+01
COFD(6162) = 5.18417470d+00
COFD(6163) = -4.49491573d-01
COFD(6164) = 1.91438508d-02
COFD(6165) = -1.91102652d+01
COFD(6166) = 4.84384483d+00
COFD(6167) = -4.10265575d-01
COFD(6168) = 1.76414287d-02
COFD(6169) = -1.87383952d+01
COFD(6170) = 3.96926341d+00
COFD(6171) = -2.16412264d-01
COFD(6172) = 6.06012078d-03
COFD(6173) = -2.05272328d+01
COFD(6174) = 5.18417470d+00
COFD(6175) = -4.49491573d-01
COFD(6176) = 1.91438508d-02
COFD(6177) = -2.05356023d+01
COFD(6178) = 5.18417470d+00
COFD(6179) = -4.49491573d-01
COFD(6180) = 1.91438508d-02
COFD(6181) = -1.87688110d+01
COFD(6182) = 4.71729964d+00
COFD(6183) = -3.95432573d-01
COFD(6184) = 1.70623691d-02
COFD(6185) = -1.90116191d+01
COFD(6186) = 4.84384483d+00
COFD(6187) = -4.10265575d-01
COFD(6188) = 1.76414287d-02
COFD(6189) = -2.11349086d+01
COFD(6190) = 5.42846112d+00
COFD(6191) = -4.74321870d-01
COFD(6192) = 1.99459749d-02
COFD(6193) = -2.11349086d+01
COFD(6194) = 5.42846112d+00
COFD(6195) = -4.74321870d-01
COFD(6196) = 1.99459749d-02
COFD(6197) = -2.11606963d+01
COFD(6198) = 5.42846112d+00
COFD(6199) = -4.74321870d-01
COFD(6200) = 1.99459749d-02
COFD(6201) = -2.11309207d+01
COFD(6202) = 5.41773516d+00
COFD(6203) = -4.73414338d-01
COFD(6204) = 1.99258685d-02
COFD(6205) = -2.02922701d+01
COFD(6206) = 5.11106992d+00
COFD(6207) = -4.42047129d-01
COFD(6208) = 1.89042990d-02
COFD(6209) = -2.22116706d+01
COFD(6210) = 5.54251230d+00
COFD(6211) = -4.70946314d-01
COFD(6212) = 1.90785869d-02
COFD(6213) = -2.00915040d+01
COFD(6214) = 4.41511629d+00
COFD(6215) = -2.84086963d-01
COFD(6216) = 9.37586971d-03
COFD(6217) = -2.01015340d+01
COFD(6218) = 4.41511629d+00
COFD(6219) = -2.84086963d-01
COFD(6220) = 9.37586971d-03
COFD(6221) = -2.09222454d+01
COFD(6222) = 4.82184721d+00
COFD(6223) = -3.48128875d-01
COFD(6224) = 1.25918978d-02
COFD(6225) = -2.09222454d+01
COFD(6226) = 4.82184721d+00
COFD(6227) = -3.48128875d-01
COFD(6228) = 1.25918978d-02
COFD(6229) = -2.03036402d+01
COFD(6230) = 4.50250781d+00
COFD(6231) = -2.97622106d-01
COFD(6232) = 1.00481473d-02
COFD(6233) = -2.20942491d+01
COFD(6234) = 5.58360799d+00
COFD(6235) = -4.82701436d-01
COFD(6236) = 1.98437922d-02
COFD(6237) = -2.21065306d+01
COFD(6238) = 5.58360799d+00
COFD(6239) = -4.82701436d-01
COFD(6240) = 1.98437922d-02
COFD(6241) = -2.21181720d+01
COFD(6242) = 5.58360799d+00
COFD(6243) = -4.82701436d-01
COFD(6244) = 1.98437922d-02
COFD(6245) = -2.20453723d+01
COFD(6246) = 5.44448440d+00
COFD(6247) = -4.51529024d-01
COFD(6248) = 1.79698119d-02
COFD(6249) = -2.22059540d+01
COFD(6250) = 5.51722375d+00
COFD(6251) = -4.66081431d-01
COFD(6252) = 1.88044011d-02
COFD(6253) = -2.22159630d+01
COFD(6254) = 5.51722375d+00
COFD(6255) = -4.66081431d-01
COFD(6256) = 1.88044011d-02
COFD(6257) = -2.12621914d+01
COFD(6258) = 5.47935225d+00
COFD(6259) = -4.80056796d-01
COFD(6260) = 2.01607180d-02
COFD(6261) = -2.09002742d+01
COFD(6262) = 4.72895031d+00
COFD(6263) = -3.33332771d-01
COFD(6264) = 1.18431478d-02
COFD(6265) = -2.09002742d+01
COFD(6266) = 4.72895031d+00
COFD(6267) = -3.33332771d-01
COFD(6268) = 1.18431478d-02
COFD(6269) = -1.88275332d+01
COFD(6270) = 4.71729964d+00
COFD(6271) = -3.95432573d-01
COFD(6272) = 1.70623691d-02
COFD(6273) = -1.90351359d+01
COFD(6274) = 4.84384483d+00
COFD(6275) = -4.10265575d-01
COFD(6276) = 1.76414287d-02
COFD(6277) = -1.90588668d+01
COFD(6278) = 4.84384483d+00
COFD(6279) = -4.10265575d-01
COFD(6280) = 1.76414287d-02
COFD(6281) = -1.99140603d+01
COFD(6282) = 4.50672715d+00
COFD(6283) = -2.98278922d-01
COFD(6284) = 1.00808824d-02
COFD(6285) = -1.92061609d+01
COFD(6286) = 4.71729964d+00
COFD(6287) = -3.95432573d-01
COFD(6288) = 1.70623691d-02
COFD(6289) = -2.02849953d+01
COFD(6290) = 5.10426133d+00
COFD(6291) = -4.41256919d-01
COFD(6292) = 1.88737290d-02
COFD(6293) = -2.20848247d+01
COFD(6294) = 5.59135292d+00
COFD(6295) = -4.85565630d-01
COFD(6296) = 2.00429035d-02
COFD(6297) = -2.22703814d+01
COFD(6298) = 5.57944619d+00
COFD(6299) = -4.78032626d-01
COFD(6300) = 1.94775370d-02
COFD(6301) = -2.08481163d+01
COFD(6302) = 5.29805436d+00
COFD(6303) = -4.62913371d-01
COFD(6304) = 1.96725621d-02
COFD(6305) = -1.93551579d+01
COFD(6306) = 4.77546552d+00
COFD(6307) = -4.02345606d-01
COFD(6308) = 1.73365265d-02
COFD(6309) = -1.92068539d+01
COFD(6310) = 3.99103591d+00
COFD(6311) = -2.19677673d-01
COFD(6312) = 6.21875426d-03
COFD(6313) = -1.92179182d+01
COFD(6314) = 3.99103591d+00
COFD(6315) = -2.19677673d-01
COFD(6316) = 6.21875426d-03
COFD(6317) = -2.12622090d+01
COFD(6318) = 5.47935225d+00
COFD(6319) = -4.80056796d-01
COFD(6320) = 2.01607180d-02
COFD(6321) = -2.22648059d+01
COFD(6322) = 5.57944619d+00
COFD(6323) = -4.78032626d-01
COFD(6324) = 1.94775370d-02
COFD(6325) = -2.22648059d+01
COFD(6326) = 5.57944619d+00
COFD(6327) = -4.78032626d-01
COFD(6328) = 1.94775370d-02
COFD(6329) = -2.22648059d+01
COFD(6330) = 5.57944619d+00
COFD(6331) = -4.78032626d-01
COFD(6332) = 1.94775370d-02
COFD(6333) = -2.22589130d+01
COFD(6334) = 5.57944619d+00
COFD(6335) = -4.78032626d-01
COFD(6336) = 1.94775370d-02
COFD(6337) = -2.02646611d+01
COFD(6338) = 5.10426133d+00
COFD(6339) = -4.41256919d-01
COFD(6340) = 1.88737290d-02
COFD(6341) = -2.12450965d+01
COFD(6342) = 5.40444222d+00
COFD(6343) = -4.72708609d-01
COFD(6344) = 1.99362392d-02
COFD(6345) = -2.24674263d+01
COFD(6346) = 5.49330641d+00
COFD(6347) = -4.60498247d-01
COFD(6348) = 1.84639199d-02
COFD(6349) = -2.24731023d+01
COFD(6350) = 5.49330641d+00
COFD(6351) = -4.60498247d-01
COFD(6352) = 1.84639199d-02
COFD(6353) = -2.09061629d+01
COFD(6354) = 4.72895031d+00
COFD(6355) = -3.33332771d-01
COFD(6356) = 1.18431478d-02
COFD(6357) = -2.09118474d+01
COFD(6358) = 4.72895031d+00
COFD(6359) = -3.33332771d-01
COFD(6360) = 1.18431478d-02
COFD(6361) = -1.08472922d+01
COFD(6362) = 2.19094415d+00
COFD(6363) = -7.11992510d-02
COFD(6364) = 3.14105973d-03
COFD(6365) = -1.32522778d+01
COFD(6366) = 3.34156587d+00
COFD(6367) = -2.22853306d-01
COFD(6368) = 9.81883417d-03
COFD(6369) = -1.31034488d+01
COFD(6370) = 2.80913567d+00
COFD(6371) = -1.54536855d-01
COFD(6372) = 6.89359313d-03
COFD(6373) = -1.41254249d+01
COFD(6374) = 3.05837263d+00
COFD(6375) = -1.86672802d-01
COFD(6376) = 8.27575734d-03
COFD(6377) = -1.31174764d+01
COFD(6378) = 2.80913567d+00
COFD(6379) = -1.54536855d-01
COFD(6380) = 6.89359313d-03
COFD(6381) = -1.92450597d+01
COFD(6382) = 5.05708283d+00
COFD(6383) = -4.35739290d-01
COFD(6384) = 1.86583205d-02
COFD(6385) = -1.41300955d+01
COFD(6386) = 3.05837263d+00
COFD(6387) = -1.86672802d-01
COFD(6388) = 8.27575734d-03
COFD(6389) = -1.41345294d+01
COFD(6390) = 3.05837263d+00
COFD(6391) = -1.86672802d-01
COFD(6392) = 8.27575734d-03
COFD(6393) = -1.29942577d+01
COFD(6394) = 2.73155251d+00
COFD(6395) = -1.44594082d-01
COFD(6396) = 6.46883252d-03
COFD(6397) = -1.30526867d+01
COFD(6398) = 2.80913567d+00
COFD(6399) = -1.54536855d-01
COFD(6400) = 6.89359313d-03
COFD(6401) = -1.47958130d+01
COFD(6402) = 3.33113524d+00
COFD(6403) = -2.21479057d-01
COFD(6404) = 9.75837737d-03
COFD(6405) = -1.47958130d+01
COFD(6406) = 3.33113524d+00
COFD(6407) = -2.21479057d-01
COFD(6408) = 9.75837737d-03
COFD(6409) = -1.48128481d+01
COFD(6410) = 3.33113524d+00
COFD(6411) = -2.21479057d-01
COFD(6412) = 9.75837737d-03
COFD(6413) = -1.47472946d+01
COFD(6414) = 3.30594991d+00
COFD(6415) = -2.18182207d-01
COFD(6416) = 9.61429447d-03
COFD(6417) = -1.39388049d+01
COFD(6418) = 2.97564184d+00
COFD(6419) = -1.76025309d-01
COFD(6420) = 7.81869993d-03
COFD(6421) = -1.67983919d+01
COFD(6422) = 4.00828594d+00
COFD(6423) = -3.08414344d-01
COFD(6424) = 1.34907430d-02
COFD(6425) = -1.91057988d+01
COFD(6426) = 4.86821670d+00
COFD(6427) = -4.13144121d-01
COFD(6428) = 1.77546701d-02
COFD(6429) = -1.91112931d+01
COFD(6430) = 4.86821670d+00
COFD(6431) = -4.13144121d-01
COFD(6432) = 1.77546701d-02
COFD(6433) = -1.85699333d+01
COFD(6434) = 4.67076124d+00
COFD(6435) = -3.90022427d-01
COFD(6436) = 1.68533953d-02
COFD(6437) = -1.85699333d+01
COFD(6438) = 4.67076124d+00
COFD(6439) = -3.90022427d-01
COFD(6440) = 1.68533953d-02
COFD(6441) = -1.90218743d+01
COFD(6442) = 4.83076737d+00
COFD(6443) = -4.08802573d-01
COFD(6444) = 1.75875241d-02
COFD(6445) = -1.61544946d+01
COFD(6446) = 3.75910622d+00
COFD(6447) = -2.75986578d-01
COFD(6448) = 1.20782843d-02
COFD(6449) = -1.61614882d+01
COFD(6450) = 3.75910622d+00
COFD(6451) = -2.75986578d-01
COFD(6452) = 1.20782843d-02
COFD(6453) = -1.61680488d+01
COFD(6454) = 3.75910622d+00
COFD(6455) = -2.75986578d-01
COFD(6456) = 1.20782843d-02
COFD(6457) = -1.72409687d+01
COFD(6458) = 4.17190426d+00
COFD(6459) = -3.28894681d-01
COFD(6460) = 1.43498101d-02
COFD(6461) = -1.70022408d+01
COFD(6462) = 4.05099737d+00
COFD(6463) = -3.13841660d-01
COFD(6464) = 1.37218854d-02
COFD(6465) = -1.70077215d+01
COFD(6466) = 4.05099737d+00
COFD(6467) = -3.13841660d-01
COFD(6468) = 1.37218854d-02
COFD(6469) = -1.47490868d+01
COFD(6470) = 3.39229020d+00
COFD(6471) = -2.29520232d-01
COFD(6472) = 1.01114311d-02
COFD(6473) = -1.88275332d+01
COFD(6474) = 4.71729964d+00
COFD(6475) = -3.95432573d-01
COFD(6476) = 1.70623691d-02
COFD(6477) = -1.88275332d+01
COFD(6478) = 4.71729964d+00
COFD(6479) = -3.95432573d-01
COFD(6480) = 1.70623691d-02
COFD(6481) = -1.30341578d+01
COFD(6482) = 2.73155251d+00
COFD(6483) = -1.44594082d-01
COFD(6484) = 6.46883252d-03
COFD(6485) = -1.30550259d+01
COFD(6486) = 2.80913567d+00
COFD(6487) = -1.54536855d-01
COFD(6488) = 6.89359313d-03
COFD(6489) = -1.30704422d+01
COFD(6490) = 2.80913567d+00
COFD(6491) = -1.54536855d-01
COFD(6492) = 6.89359313d-03
COFD(6493) = -1.86830405d+01
COFD(6494) = 4.82909392d+00
COFD(6495) = -4.08610711d-01
COFD(6496) = 1.75802236d-02
COFD(6497) = -1.33299570d+01
COFD(6498) = 2.73155251d+00
COFD(6499) = -1.44594082d-01
COFD(6500) = 6.46883252d-03
COFD(6501) = -1.39298633d+01
COFD(6502) = 2.97137588d+00
COFD(6503) = -1.75491257d-01
COFD(6504) = 7.79646773d-03
COFD(6505) = -1.58784322d+01
COFD(6506) = 3.68407693d+00
COFD(6507) = -2.66228170d-01
COFD(6508) = 1.16542305d-02
COFD(6509) = -1.66315280d+01
COFD(6510) = 3.93849401d+00
COFD(6511) = -2.99416642d-01
COFD(6512) = 1.31020815d-02
COFD(6513) = -1.43591671d+01
COFD(6514) = 3.14480429d+00
COFD(6515) = -1.97906290d-01
COFD(6516) = 8.76325718d-03
COFD(6517) = -1.34177769d+01
COFD(6518) = 2.76469604d+00
COFD(6519) = -1.48843508d-01
COFD(6520) = 6.65045597d-03
COFD(6521) = -1.96225223d+01
COFD(6522) = 5.04961455d+00
COFD(6523) = -4.34856114d-01
COFD(6524) = 1.86234058d-02
COFD(6525) = -1.96286971d+01
COFD(6526) = 5.04961455d+00
COFD(6527) = -4.34856114d-01
COFD(6528) = 1.86234058d-02
COFD(6529) = -1.47490956d+01
COFD(6530) = 3.39229020d+00
COFD(6531) = -2.29520232d-01
COFD(6532) = 1.01114311d-02
COFD(6533) = -1.66287650d+01
COFD(6534) = 3.93849401d+00
COFD(6535) = -2.99416642d-01
COFD(6536) = 1.31020815d-02
COFD(6537) = -1.66287650d+01
COFD(6538) = 3.93849401d+00
COFD(6539) = -2.99416642d-01
COFD(6540) = 1.31020815d-02
COFD(6541) = -1.66287650d+01
COFD(6542) = 3.93849401d+00
COFD(6543) = -2.99416642d-01
COFD(6544) = 1.31020815d-02
COFD(6545) = -1.66258278d+01
COFD(6546) = 3.93849401d+00
COFD(6547) = -2.99416642d-01
COFD(6548) = 1.31020815d-02
COFD(6549) = -1.39186707d+01
COFD(6550) = 2.97137588d+00
COFD(6551) = -1.75491257d-01
COFD(6552) = 7.79646773d-03
COFD(6553) = -1.47041948d+01
COFD(6554) = 3.27505697d+00
COFD(6555) = -2.14306851d-01
COFD(6556) = 9.45219335d-03
COFD(6557) = -1.74132498d+01
COFD(6558) = 4.11838823d+00
COFD(6559) = -3.22329602d-01
COFD(6560) = 1.40802468d-02
COFD(6561) = -1.74160614d+01
COFD(6562) = 4.11838823d+00
COFD(6563) = -3.22329602d-01
COFD(6564) = 1.40802468d-02
COFD(6565) = -1.88304679d+01
COFD(6566) = 4.71729964d+00
COFD(6567) = -3.95432573d-01
COFD(6568) = 1.70623691d-02
COFD(6569) = -1.88332845d+01
COFD(6570) = 4.71729964d+00
COFD(6571) = -3.95432573d-01
COFD(6572) = 1.70623691d-02
COFD(6573) = -1.09203269d+01
COFD(6574) = 2.30836460d+00
COFD(6575) = -8.76339315d-02
COFD(6576) = 3.90878445d-03
COFD(6577) = -1.33789807d+01
COFD(6578) = 3.48624238d+00
COFD(6579) = -2.41554467d-01
COFD(6580) = 1.06263545d-02
COFD(6581) = -1.31565315d+01
COFD(6582) = 2.90778936d+00
COFD(6583) = -1.67388544d-01
COFD(6584) = 7.45220609d-03
COFD(6585) = -1.42600473d+01
COFD(6586) = 3.17651319d+00
COFD(6587) = -2.02028974d-01
COFD(6588) = 8.94232502d-03
COFD(6589) = -1.31710876d+01
COFD(6590) = 2.90778936d+00
COFD(6591) = -1.67388544d-01
COFD(6592) = 7.45220609d-03
COFD(6593) = -1.93545828d+01
COFD(6594) = 5.16013126d+00
COFD(6595) = -4.46824543d-01
COFD(6596) = 1.90464887d-02
COFD(6597) = -1.42649477d+01
COFD(6598) = 3.17651319d+00
COFD(6599) = -2.02028974d-01
COFD(6600) = 8.94232502d-03
COFD(6601) = -1.42696020d+01
COFD(6602) = 3.17651319d+00
COFD(6603) = -2.02028974d-01
COFD(6604) = 8.94232502d-03
COFD(6605) = -1.30137956d+01
COFD(6606) = 2.80913567d+00
COFD(6607) = -1.54536855d-01
COFD(6608) = 6.89359313d-03
COFD(6609) = -1.31039806d+01
COFD(6610) = 2.90778936d+00
COFD(6611) = -1.67388544d-01
COFD(6612) = 7.45220609d-03
COFD(6613) = -1.50125659d+01
COFD(6614) = 3.47945612d+00
COFD(6615) = -2.40703722d-01
COFD(6616) = 1.05907441d-02
COFD(6617) = -1.50125659d+01
COFD(6618) = 3.47945612d+00
COFD(6619) = -2.40703722d-01
COFD(6620) = 1.05907441d-02
COFD(6621) = -1.50302037d+01
COFD(6622) = 3.47945612d+00
COFD(6623) = -2.40703722d-01
COFD(6624) = 1.05907441d-02
COFD(6625) = -1.49798516d+01
COFD(6626) = 3.46140064d+00
COFD(6627) = -2.38440092d-01
COFD(6628) = 1.04960087d-02
COFD(6629) = -1.40479570d+01
COFD(6630) = 3.08120012d+00
COFD(6631) = -1.89629903d-01
COFD(6632) = 8.40361952d-03
COFD(6633) = -1.69990507d+01
COFD(6634) = 4.14240922d+00
COFD(6635) = -3.25239774d-01
COFD(6636) = 1.41980687d-02
COFD(6637) = -1.93788111d+01
COFD(6638) = 5.02567894d+00
COFD(6639) = -4.32045169d-01
COFD(6640) = 1.85132214d-02
COFD(6641) = -1.93845675d+01
COFD(6642) = 5.02567894d+00
COFD(6643) = -4.32045169d-01
COFD(6644) = 1.85132214d-02
COFD(6645) = -1.87654600d+01
COFD(6646) = 4.79683898d+00
COFD(6647) = -4.04829719d-01
COFD(6648) = 1.74325475d-02
COFD(6649) = -1.87654600d+01
COFD(6650) = 4.79683898d+00
COFD(6651) = -4.04829719d-01
COFD(6652) = 1.74325475d-02
COFD(6653) = -1.92834357d+01
COFD(6654) = 4.98286777d+00
COFD(6655) = -4.26970814d-01
COFD(6656) = 1.83122917d-02
COFD(6657) = -1.64922031d+01
COFD(6658) = 3.95035840d+00
COFD(6659) = -3.00959418d-01
COFD(6660) = 1.31692593d-02
COFD(6661) = -1.64995136d+01
COFD(6662) = 3.95035840d+00
COFD(6663) = -3.00959418d-01
COFD(6664) = 1.31692593d-02
COFD(6665) = -1.65063757d+01
COFD(6666) = 3.95035840d+00
COFD(6667) = -3.00959418d-01
COFD(6668) = 1.31692593d-02
COFD(6669) = -1.74287716d+01
COFD(6670) = 4.29676909d+00
COFD(6671) = -3.44085306d-01
COFD(6672) = 1.49671135d-02
COFD(6673) = -1.72003855d+01
COFD(6674) = 4.17889917d+00
COFD(6675) = -3.29752510d-01
COFD(6676) = 1.43850275d-02
COFD(6677) = -1.72061277d+01
COFD(6678) = 4.17889917d+00
COFD(6679) = -3.29752510d-01
COFD(6680) = 1.43850275d-02
COFD(6681) = -1.48885202d+01
COFD(6682) = 3.52327209d+00
COFD(6683) = -2.46286208d-01
COFD(6684) = 1.08285963d-02
COFD(6685) = -1.90351359d+01
COFD(6686) = 4.84384483d+00
COFD(6687) = -4.10265575d-01
COFD(6688) = 1.76414287d-02
COFD(6689) = -1.90351359d+01
COFD(6690) = 4.84384483d+00
COFD(6691) = -4.10265575d-01
COFD(6692) = 1.76414287d-02
COFD(6693) = -1.30550259d+01
COFD(6694) = 2.80913567d+00
COFD(6695) = -1.54536855d-01
COFD(6696) = 6.89359313d-03
COFD(6697) = -1.31035185d+01
COFD(6698) = 2.90778936d+00
COFD(6699) = -1.67388544d-01
COFD(6700) = 7.45220609d-03
COFD(6701) = -1.31194985d+01
COFD(6702) = 2.90778936d+00
COFD(6703) = -1.67388544d-01
COFD(6704) = 7.45220609d-03
COFD(6705) = -1.89147380d+01
COFD(6706) = 4.98075560d+00
COFD(6707) = -4.26721620d-01
COFD(6708) = 1.83024823d-02
COFD(6709) = -1.33722163d+01
COFD(6710) = 2.80913567d+00
COFD(6711) = -1.54536855d-01
COFD(6712) = 6.89359313d-03
COFD(6713) = -1.40353277d+01
COFD(6714) = 3.07549274d+00
COFD(6715) = -1.88889344d-01
COFD(6716) = 8.37152866d-03
COFD(6717) = -1.62212061d+01
COFD(6718) = 3.88250968d+00
COFD(6719) = -2.92155848d-01
COFD(6720) = 1.27867850d-02
COFD(6721) = -1.68748824d+01
COFD(6722) = 4.09077642d+00
COFD(6723) = -3.18894990d-01
COFD(6724) = 1.39371445d-02
COFD(6725) = -1.44530825d+01
COFD(6726) = 3.24450689d+00
COFD(6727) = -2.10570734d-01
COFD(6728) = 9.30026771d-03
COFD(6729) = -1.34898944d+01
COFD(6730) = 2.85448298d+00
COFD(6731) = -1.60491387d-01
COFD(6732) = 7.15451323d-03
COFD(6733) = -1.97825930d+01
COFD(6734) = 5.15754807d+00
COFD(6735) = -4.46654084d-01
COFD(6736) = 1.90458977d-02
COFD(6737) = -1.97890553d+01
COFD(6738) = 5.15754807d+00
COFD(6739) = -4.46654084d-01
COFD(6740) = 1.90458977d-02
COFD(6741) = -1.48885295d+01
COFD(6742) = 3.52327209d+00
COFD(6743) = -2.46286208d-01
COFD(6744) = 1.08285963d-02
COFD(6745) = -1.68719715d+01
COFD(6746) = 4.09077642d+00
COFD(6747) = -3.18894990d-01
COFD(6748) = 1.39371445d-02
COFD(6749) = -1.68719715d+01
COFD(6750) = 4.09077642d+00
COFD(6751) = -3.18894990d-01
COFD(6752) = 1.39371445d-02
COFD(6753) = -1.68719715d+01
COFD(6754) = 4.09077642d+00
COFD(6755) = -3.18894990d-01
COFD(6756) = 1.39371445d-02
COFD(6757) = -1.68688781d+01
COFD(6758) = 4.09077642d+00
COFD(6759) = -3.18894990d-01
COFD(6760) = 1.39371445d-02
COFD(6761) = -1.40236044d+01
COFD(6762) = 3.07549274d+00
COFD(6763) = -1.88889344d-01
COFD(6764) = 8.37152866d-03
COFD(6765) = -1.49050015d+01
COFD(6766) = 3.41988961d+00
COFD(6767) = -2.33128386d-01
COFD(6768) = 1.02689994d-02
COFD(6769) = -1.76347103d+01
COFD(6770) = 4.24719726d+00
COFD(6771) = -3.38206061d-01
COFD(6772) = 1.47350654d-02
COFD(6773) = -1.76376724d+01
COFD(6774) = 4.24719726d+00
COFD(6775) = -3.38206061d-01
COFD(6776) = 1.47350654d-02
COFD(6777) = -1.90382267d+01
COFD(6778) = 4.84384483d+00
COFD(6779) = -4.10265575d-01
COFD(6780) = 1.76414287d-02
COFD(6781) = -1.90411940d+01
COFD(6782) = 4.84384483d+00
COFD(6783) = -4.10265575d-01
COFD(6784) = 1.76414287d-02
COFD(6785) = -1.09240642d+01
COFD(6786) = 2.30836460d+00
COFD(6787) = -8.76339315d-02
COFD(6788) = 3.90878445d-03
COFD(6789) = -1.33809634d+01
COFD(6790) = 3.48624238d+00
COFD(6791) = -2.41554467d-01
COFD(6792) = 1.06263545d-02
COFD(6793) = -1.31730273d+01
COFD(6794) = 2.90778936d+00
COFD(6795) = -1.67388544d-01
COFD(6796) = 7.45220609d-03
COFD(6797) = -1.42819281d+01
COFD(6798) = 3.17651319d+00
COFD(6799) = -2.02028974d-01
COFD(6800) = 8.94232502d-03
COFD(6801) = -1.31880790d+01
COFD(6802) = 2.90778936d+00
COFD(6803) = -1.67388544d-01
COFD(6804) = 7.45220609d-03
COFD(6805) = -2.10227861d+01
COFD(6806) = 5.58489889d+00
COFD(6807) = -4.81572453d-01
COFD(6808) = 1.97432291d-02
COFD(6809) = -1.42870488d+01
COFD(6810) = 3.17651319d+00
COFD(6811) = -2.02028974d-01
COFD(6812) = 8.94232502d-03
COFD(6813) = -1.42919145d+01
COFD(6814) = 3.17651319d+00
COFD(6815) = -2.02028974d-01
COFD(6816) = 8.94232502d-03
COFD(6817) = -1.30279742d+01
COFD(6818) = 2.80913567d+00
COFD(6819) = -1.54536855d-01
COFD(6820) = 6.89359313d-03
COFD(6821) = -1.31188060d+01
COFD(6822) = 2.90778936d+00
COFD(6823) = -1.67388544d-01
COFD(6824) = 7.45220609d-03
COFD(6825) = -1.50279940d+01
COFD(6826) = 3.47945612d+00
COFD(6827) = -2.40703722d-01
COFD(6828) = 1.05907441d-02
COFD(6829) = -1.50279940d+01
COFD(6830) = 3.47945612d+00
COFD(6831) = -2.40703722d-01
COFD(6832) = 1.05907441d-02
COFD(6833) = -1.50461946d+01
COFD(6834) = 3.47945612d+00
COFD(6835) = -2.40703722d-01
COFD(6836) = 1.05907441d-02
COFD(6837) = -1.49963695d+01
COFD(6838) = 3.46140064d+00
COFD(6839) = -2.38440092d-01
COFD(6840) = 1.04960087d-02
COFD(6841) = -1.40688658d+01
COFD(6842) = 3.08120012d+00
COFD(6843) = -1.89629903d-01
COFD(6844) = 8.40361952d-03
COFD(6845) = -1.70230718d+01
COFD(6846) = 4.14240922d+00
COFD(6847) = -3.25239774d-01
COFD(6848) = 1.41980687d-02
COFD(6849) = -1.93999821d+01
COFD(6850) = 5.02567894d+00
COFD(6851) = -4.32045169d-01
COFD(6852) = 1.85132214d-02
COFD(6853) = -1.94059889d+01
COFD(6854) = 5.02567894d+00
COFD(6855) = -4.32045169d-01
COFD(6856) = 1.85132214d-02
COFD(6857) = -1.97755767d+01
COFD(6858) = 5.13121203d+00
COFD(6859) = -4.44310112d-01
COFD(6860) = 1.89883544d-02
COFD(6861) = -1.97755767d+01
COFD(6862) = 5.13121203d+00
COFD(6863) = -4.44310112d-01
COFD(6864) = 1.89883544d-02
COFD(6865) = -1.93053262d+01
COFD(6866) = 4.98286777d+00
COFD(6867) = -4.26970814d-01
COFD(6868) = 1.83122917d-02
COFD(6869) = -1.65122609d+01
COFD(6870) = 3.95035840d+00
COFD(6871) = -3.00959418d-01
COFD(6872) = 1.31692593d-02
COFD(6873) = -1.65198729d+01
COFD(6874) = 3.95035840d+00
COFD(6875) = -3.00959418d-01
COFD(6876) = 1.31692593d-02
COFD(6877) = -1.65270223d+01
COFD(6878) = 3.95035840d+00
COFD(6879) = -3.00959418d-01
COFD(6880) = 1.31692593d-02
COFD(6881) = -1.74496921d+01
COFD(6882) = 4.29676909d+00
COFD(6883) = -3.44085306d-01
COFD(6884) = 1.49671135d-02
COFD(6885) = -1.72215676d+01
COFD(6886) = 4.17889917d+00
COFD(6887) = -3.29752510d-01
COFD(6888) = 1.43850275d-02
COFD(6889) = -1.72275598d+01
COFD(6890) = 4.17889917d+00
COFD(6891) = -3.29752510d-01
COFD(6892) = 1.43850275d-02
COFD(6893) = -1.49120950d+01
COFD(6894) = 3.52327209d+00
COFD(6895) = -2.46286208d-01
COFD(6896) = 1.08285963d-02
COFD(6897) = -1.90588668d+01
COFD(6898) = 4.84384483d+00
COFD(6899) = -4.10265575d-01
COFD(6900) = 1.76414287d-02
COFD(6901) = -1.90588668d+01
COFD(6902) = 4.84384483d+00
COFD(6903) = -4.10265575d-01
COFD(6904) = 1.76414287d-02
COFD(6905) = -1.30704422d+01
COFD(6906) = 2.80913567d+00
COFD(6907) = -1.54536855d-01
COFD(6908) = 6.89359313d-03
COFD(6909) = -1.31194985d+01
COFD(6910) = 2.90778936d+00
COFD(6911) = -1.67388544d-01
COFD(6912) = 7.45220609d-03
COFD(6913) = -1.31360060d+01
COFD(6914) = 2.90778936d+00
COFD(6915) = -1.67388544d-01
COFD(6916) = 7.45220609d-03
COFD(6917) = -2.01575364d+01
COFD(6918) = 5.37354588d+00
COFD(6919) = -4.70508014d-01
COFD(6920) = 1.99136680d-02
COFD(6921) = -1.33933879d+01
COFD(6922) = 2.80913567d+00
COFD(6923) = -1.54536855d-01
COFD(6924) = 6.89359313d-03
COFD(6925) = -1.40567442d+01
COFD(6926) = 3.07549274d+00
COFD(6927) = -1.88889344d-01
COFD(6928) = 8.37152866d-03
COFD(6929) = -1.62455018d+01
COFD(6930) = 3.88250968d+00
COFD(6931) = -2.92155848d-01
COFD(6932) = 1.27867850d-02
COFD(6933) = -1.68989038d+01
COFD(6934) = 4.09077642d+00
COFD(6935) = -3.18894990d-01
COFD(6936) = 1.39371445d-02
COFD(6937) = -1.44747388d+01
COFD(6938) = 3.24450689d+00
COFD(6939) = -2.10570734d-01
COFD(6940) = 9.30026771d-03
COFD(6941) = -1.35102477d+01
COFD(6942) = 2.85448298d+00
COFD(6943) = -1.60491387d-01
COFD(6944) = 7.15451323d-03
COFD(6945) = -1.98032338d+01
COFD(6946) = 5.15754807d+00
COFD(6947) = -4.46654084d-01
COFD(6948) = 1.90458977d-02
COFD(6949) = -1.98099704d+01
COFD(6950) = 5.15754807d+00
COFD(6951) = -4.46654084d-01
COFD(6952) = 1.90458977d-02
COFD(6953) = -1.49121048d+01
COFD(6954) = 3.52327209d+00
COFD(6955) = -2.46286208d-01
COFD(6956) = 1.08285963d-02
COFD(6957) = -1.68958501d+01
COFD(6958) = 4.09077642d+00
COFD(6959) = -3.18894990d-01
COFD(6960) = 1.39371445d-02
COFD(6961) = -1.68958501d+01
COFD(6962) = 4.09077642d+00
COFD(6963) = -3.18894990d-01
COFD(6964) = 1.39371445d-02
COFD(6965) = -1.68958501d+01
COFD(6966) = 4.09077642d+00
COFD(6967) = -3.18894990d-01
COFD(6968) = 1.39371445d-02
COFD(6969) = -1.68926059d+01
COFD(6970) = 4.09077642d+00
COFD(6971) = -3.18894990d-01
COFD(6972) = 1.39371445d-02
COFD(6973) = -1.40445141d+01
COFD(6974) = 3.07549274d+00
COFD(6975) = -1.88889344d-01
COFD(6976) = 8.37152866d-03
COFD(6977) = -1.49284025d+01
COFD(6978) = 3.41988961d+00
COFD(6979) = -2.33128386d-01
COFD(6980) = 1.02689994d-02
COFD(6981) = -1.76585983d+01
COFD(6982) = 4.24719726d+00
COFD(6983) = -3.38206061d-01
COFD(6984) = 1.47350654d-02
COFD(6985) = -1.76617059d+01
COFD(6986) = 4.24719726d+00
COFD(6987) = -3.38206061d-01
COFD(6988) = 1.47350654d-02
COFD(6989) = -1.90621083d+01
COFD(6990) = 4.84384483d+00
COFD(6991) = -4.10265575d-01
COFD(6992) = 1.76414287d-02
COFD(6993) = -1.90652212d+01
COFD(6994) = 4.84384483d+00
COFD(6995) = -4.10265575d-01
COFD(6996) = 1.76414287d-02
COFD(6997) = -1.61691241d+01
COFD(6998) = 4.23579857d+00
COFD(6999) = -3.36814551d-01
COFD(7000) = 1.46782923d-02
COFD(7001) = -1.94027559d+01
COFD(7002) = 5.54487588d+00
COFD(7003) = -4.86919661d-01
COFD(7004) = 2.03935297d-02
COFD(7005) = -1.89669545d+01
COFD(7006) = 4.98075560d+00
COFD(7007) = -4.26721620d-01
COFD(7008) = 1.83024823d-02
COFD(7009) = -2.06107018d+01
COFD(7010) = 5.38762344d+00
COFD(7011) = -4.71526404d-01
COFD(7012) = 1.99251819d-02
COFD(7013) = -1.89824721d+01
COFD(7014) = 4.98075560d+00
COFD(7015) = -4.26721620d-01
COFD(7016) = 1.83024823d-02
COFD(7017) = -1.47127617d+01
COFD(7018) = 2.26104173d+00
COFD(7019) = 4.18075774d-02
COFD(7020) = -6.41614212d-03
COFD(7021) = -2.03661760d+01
COFD(7022) = 5.32029196d+00
COFD(7023) = -4.65578114d-01
COFD(7024) = 1.97797257d-02
COFD(7025) = -2.03712446d+01
COFD(7026) = 5.32029196d+00
COFD(7027) = -4.65578114d-01
COFD(7028) = 1.97797257d-02
COFD(7029) = -1.86394179d+01
COFD(7030) = 4.82909392d+00
COFD(7031) = -4.08610711d-01
COFD(7032) = 1.75802236d-02
COFD(7033) = -1.89111699d+01
COFD(7034) = 4.98075560d+00
COFD(7035) = -4.26721620d-01
COFD(7036) = 1.83024823d-02
COFD(7037) = -2.10302665d+01
COFD(7038) = 5.53859192d+00
COFD(7039) = -4.86285742d-01
COFD(7040) = 2.03731899d-02
COFD(7041) = -2.10302665d+01
COFD(7042) = 5.53859192d+00
COFD(7043) = -4.86285742d-01
COFD(7044) = 2.03731899d-02
COFD(7045) = -2.10489941d+01
COFD(7046) = 5.53859192d+00
COFD(7047) = -4.86285742d-01
COFD(7048) = 2.03731899d-02
COFD(7049) = -2.12890022d+01
COFD(7050) = 5.60080149d+00
COFD(7051) = -4.91551539d-01
COFD(7052) = 2.04914050d-02
COFD(7053) = -2.04839638d+01
COFD(7054) = 5.34597277d+00
COFD(7055) = -4.68432639d-01
COFD(7056) = 1.98846555d-02
COFD(7057) = -2.14900480d+01
COFD(7058) = 5.39920344d+00
COFD(7059) = -4.43232069d-01
COFD(7060) = 1.75138765d-02
COFD(7061) = -1.89427986d+01
COFD(7062) = 4.11490576d+00
COFD(7063) = -2.38333952d-01
COFD(7064) = 7.12818549d-03
COFD(7065) = -1.89490450d+01
COFD(7066) = 4.11490576d+00
COFD(7067) = -2.38333952d-01
COFD(7068) = 7.12818549d-03
COFD(7069) = -1.82953333d+01
COFD(7070) = 3.82052409d+00
COFD(7071) = -1.97327882d-01
COFD(7072) = 5.27318631d-03
COFD(7073) = -1.82953333d+01
COFD(7074) = 3.82052409d+00
COFD(7075) = -1.97327882d-01
COFD(7076) = 5.27318631d-03
COFD(7077) = -1.92028153d+01
COFD(7078) = 4.22627146d+00
COFD(7079) = -2.55181433d-01
COFD(7080) = 7.95198535d-03
COFD(7081) = -2.16823299d+01
COFD(7082) = 5.58175603d+00
COFD(7083) = -4.78660918d-01
COFD(7084) = 1.95178500d-02
COFD(7085) = -2.16902292d+01
COFD(7086) = 5.58175603d+00
COFD(7087) = -4.78660918d-01
COFD(7088) = 1.95178500d-02
COFD(7089) = -2.16976525d+01
COFD(7090) = 5.58175603d+00
COFD(7091) = -4.78660918d-01
COFD(7092) = 1.95178500d-02
COFD(7093) = -2.13779772d+01
COFD(7094) = 5.33947196d+00
COFD(7095) = -4.32820858d-01
COFD(7096) = 1.69568580d-02
COFD(7097) = -2.16435160d+01
COFD(7098) = 5.45344302d+00
COFD(7099) = -4.53149740d-01
COFD(7100) = 1.80583906d-02
COFD(7101) = -2.16497474d+01
COFD(7102) = 5.45344302d+00
COFD(7103) = -4.53149740d-01
COFD(7104) = 1.80583906d-02
COFD(7105) = -2.09300683d+01
COFD(7106) = 5.56884637d+00
COFD(7107) = -4.89068748d-01
COFD(7108) = 2.04486922d-02
COFD(7109) = -1.99140603d+01
COFD(7110) = 4.50672715d+00
COFD(7111) = -2.98278922d-01
COFD(7112) = 1.00808824d-02
COFD(7113) = -1.99140603d+01
COFD(7114) = 4.50672715d+00
COFD(7115) = -2.98278922d-01
COFD(7116) = 1.00808824d-02
COFD(7117) = -1.86830405d+01
COFD(7118) = 4.82909392d+00
COFD(7119) = -4.08610711d-01
COFD(7120) = 1.75802236d-02
COFD(7121) = -1.89147380d+01
COFD(7122) = 4.98075560d+00
COFD(7123) = -4.26721620d-01
COFD(7124) = 1.83024823d-02
COFD(7125) = -2.01575364d+01
COFD(7126) = 5.37354588d+00
COFD(7127) = -4.70508014d-01
COFD(7128) = 1.99136680d-02
COFD(7129) = -1.67852909d+01
COFD(7130) = 3.25474929d+00
COFD(7131) = -1.10627410d-01
COFD(7132) = 1.02617602d-03
COFD(7133) = -1.90046828d+01
COFD(7134) = 4.82909392d+00
COFD(7135) = -4.08610711d-01
COFD(7136) = 1.75802236d-02
COFD(7137) = -2.04444507d+01
COFD(7138) = 5.33236114d+00
COFD(7139) = -4.66969586d-01
COFD(7140) = 1.98332508d-02
COFD(7141) = -2.15739313d+01
COFD(7142) = 5.58518186d+00
COFD(7143) = -4.80813479d-01
COFD(7144) = 1.96787936d-02
COFD(7145) = -2.16211036d+01
COFD(7146) = 5.50657794d+00
COFD(7147) = -4.63927523d-01
COFD(7148) = 1.86800112d-02
COFD(7149) = -2.05628395d+01
COFD(7150) = 5.37634919d+00
COFD(7151) = -4.70693222d-01
COFD(7152) = 1.99144603d-02
COFD(7153) = -1.91791848d+01
COFD(7154) = 4.89100437d+00
COFD(7155) = -4.15941503d-01
COFD(7156) = 1.78696440d-02
COFD(7157) = -1.81204154d+01
COFD(7158) = 3.71386126d+00
COFD(7159) = -1.77333629d-01
COFD(7160) = 4.13979992d-03
COFD(7161) = -1.81274137d+01
COFD(7162) = 3.71386126d+00
COFD(7163) = -1.77333629d-01
COFD(7164) = 4.13979992d-03
COFD(7165) = -2.09300785d+01
COFD(7166) = 5.56884637d+00
COFD(7167) = -4.89068748d-01
COFD(7168) = 2.04486922d-02
COFD(7169) = -2.16179118d+01
COFD(7170) = 5.50657794d+00
COFD(7171) = -4.63927523d-01
COFD(7172) = 1.86800112d-02
COFD(7173) = -2.16179118d+01
COFD(7174) = 5.50657794d+00
COFD(7175) = -4.63927523d-01
COFD(7176) = 1.86800112d-02
COFD(7177) = -2.16179118d+01
COFD(7178) = 5.50657794d+00
COFD(7179) = -4.63927523d-01
COFD(7180) = 1.86800112d-02
COFD(7181) = -2.16145219d+01
COFD(7182) = 5.50657794d+00
COFD(7183) = -4.63927523d-01
COFD(7184) = 1.86800112d-02
COFD(7185) = -2.04317360d+01
COFD(7186) = 5.33236114d+00
COFD(7187) = -4.66969586d-01
COFD(7188) = 1.98332508d-02
COFD(7189) = -2.09309775d+01
COFD(7190) = 5.48417182d+00
COFD(7191) = -4.80595673d-01
COFD(7192) = 2.01807046d-02
COFD(7193) = -2.18154539d+01
COFD(7194) = 5.39990556d+00
COFD(7195) = -4.43366024d-01
COFD(7196) = 1.75213695d-02
COFD(7197) = -2.18187020d+01
COFD(7198) = 5.39990556d+00
COFD(7199) = -4.43366024d-01
COFD(7200) = 1.75213695d-02
COFD(7201) = -1.99174474d+01
COFD(7202) = 4.50672715d+00
COFD(7203) = -2.98278922d-01
COFD(7204) = 1.00808824d-02
COFD(7205) = -1.99207011d+01
COFD(7206) = 4.50672715d+00
COFD(7207) = -2.98278922d-01
COFD(7208) = 1.00808824d-02
COFD(7209) = -1.10356312d+01
COFD(7210) = 2.19094415d+00
COFD(7211) = -7.11992510d-02
COFD(7212) = 3.14105973d-03
COFD(7213) = -1.34487066d+01
COFD(7214) = 3.34156587d+00
COFD(7215) = -2.22853306d-01
COFD(7216) = 9.81883417d-03
COFD(7217) = -1.34236995d+01
COFD(7218) = 2.80913567d+00
COFD(7219) = -1.54536855d-01
COFD(7220) = 6.89359313d-03
COFD(7221) = -1.44912469d+01
COFD(7222) = 3.05837263d+00
COFD(7223) = -1.86672802d-01
COFD(7224) = 8.27575734d-03
COFD(7225) = -1.34431763d+01
COFD(7226) = 2.80913567d+00
COFD(7227) = -1.54536855d-01
COFD(7228) = 6.89359313d-03
COFD(7229) = -1.95796680d+01
COFD(7230) = 5.05708283d+00
COFD(7231) = -4.35739290d-01
COFD(7232) = 1.86583205d-02
COFD(7233) = -1.44985622d+01
COFD(7234) = 3.05837263d+00
COFD(7235) = -1.86672802d-01
COFD(7236) = 8.27575734d-03
COFD(7237) = -1.45055431d+01
COFD(7238) = 3.05837263d+00
COFD(7239) = -1.86672802d-01
COFD(7240) = 8.27575734d-03
COFD(7241) = -1.32768507d+01
COFD(7242) = 2.73155251d+00
COFD(7243) = -1.44594082d-01
COFD(7244) = 6.46883252d-03
COFD(7245) = -1.33548788d+01
COFD(7246) = 2.80913567d+00
COFD(7247) = -1.54536855d-01
COFD(7248) = 6.89359313d-03
COFD(7249) = -1.50817474d+01
COFD(7250) = 3.33113524d+00
COFD(7251) = -2.21479057d-01
COFD(7252) = 9.75837737d-03
COFD(7253) = -1.50817474d+01
COFD(7254) = 3.33113524d+00
COFD(7255) = -2.21479057d-01
COFD(7256) = 9.75837737d-03
COFD(7257) = -1.51048721d+01
COFD(7258) = 3.33113524d+00
COFD(7259) = -2.21479057d-01
COFD(7260) = 9.75837737d-03
COFD(7261) = -1.50460762d+01
COFD(7262) = 3.30594991d+00
COFD(7263) = -2.18182207d-01
COFD(7264) = 9.61429447d-03
COFD(7265) = -1.42892712d+01
COFD(7266) = 2.97564184d+00
COFD(7267) = -1.76025309d-01
COFD(7268) = 7.81869993d-03
COFD(7269) = -1.71843946d+01
COFD(7270) = 4.00828594d+00
COFD(7271) = -3.08414344d-01
COFD(7272) = 1.34907430d-02
COFD(7273) = -1.94605277d+01
COFD(7274) = 4.86821670d+00
COFD(7275) = -4.13144121d-01
COFD(7276) = 1.77546701d-02
COFD(7277) = -1.94689917d+01
COFD(7278) = 4.86821670d+00
COFD(7279) = -4.13144121d-01
COFD(7280) = 1.77546701d-02
COFD(7281) = -1.89285474d+01
COFD(7282) = 4.67076124d+00
COFD(7283) = -3.90022427d-01
COFD(7284) = 1.68533953d-02
COFD(7285) = -1.89285474d+01
COFD(7286) = 4.67076124d+00
COFD(7287) = -3.90022427d-01
COFD(7288) = 1.68533953d-02
COFD(7289) = -1.93844661d+01
COFD(7290) = 4.83076737d+00
COFD(7291) = -4.08802573d-01
COFD(7292) = 1.75875241d-02
COFD(7293) = -1.64868273d+01
COFD(7294) = 3.75910622d+00
COFD(7295) = -2.75986578d-01
COFD(7296) = 1.20782843d-02
COFD(7297) = -1.64973292d+01
COFD(7298) = 3.75910622d+00
COFD(7299) = -2.75986578d-01
COFD(7300) = 1.20782843d-02
COFD(7301) = -1.65072489d+01
COFD(7302) = 3.75910622d+00
COFD(7303) = -2.75986578d-01
COFD(7304) = 1.20782843d-02
COFD(7305) = -1.75856336d+01
COFD(7306) = 4.17190426d+00
COFD(7307) = -3.28894681d-01
COFD(7308) = 1.43498101d-02
COFD(7309) = -1.73443798d+01
COFD(7310) = 4.05099737d+00
COFD(7311) = -3.13841660d-01
COFD(7312) = 1.37218854d-02
COFD(7313) = -1.73528250d+01
COFD(7314) = 4.05099737d+00
COFD(7315) = -3.13841660d-01
COFD(7316) = 1.37218854d-02
COFD(7317) = -1.51581584d+01
COFD(7318) = 3.39229020d+00
COFD(7319) = -2.29520232d-01
COFD(7320) = 1.01114311d-02
COFD(7321) = -1.92061609d+01
COFD(7322) = 4.71729964d+00
COFD(7323) = -3.95432573d-01
COFD(7324) = 1.70623691d-02
COFD(7325) = -1.92061609d+01
COFD(7326) = 4.71729964d+00
COFD(7327) = -3.95432573d-01
COFD(7328) = 1.70623691d-02
COFD(7329) = -1.33299570d+01
COFD(7330) = 2.73155251d+00
COFD(7331) = -1.44594082d-01
COFD(7332) = 6.46883252d-03
COFD(7333) = -1.33722163d+01
COFD(7334) = 2.80913567d+00
COFD(7335) = -1.54536855d-01
COFD(7336) = 6.89359313d-03
COFD(7337) = -1.33933879d+01
COFD(7338) = 2.80913567d+00
COFD(7339) = -1.54536855d-01
COFD(7340) = 6.89359313d-03
COFD(7341) = -1.90046828d+01
COFD(7342) = 4.82909392d+00
COFD(7343) = -4.08610711d-01
COFD(7344) = 1.75802236d-02
COFD(7345) = -1.36807230d+01
COFD(7346) = 2.73155251d+00
COFD(7347) = -1.44594082d-01
COFD(7348) = 6.46883252d-03
COFD(7349) = -1.42868965d+01
COFD(7350) = 2.97137588d+00
COFD(7351) = -1.75491257d-01
COFD(7352) = 7.79646773d-03
COFD(7353) = -1.62729751d+01
COFD(7354) = 3.68407693d+00
COFD(7355) = -2.66228170d-01
COFD(7356) = 1.16542305d-02
COFD(7357) = -1.70163290d+01
COFD(7358) = 3.93849401d+00
COFD(7359) = -2.99416642d-01
COFD(7360) = 1.31020815d-02
COFD(7361) = -1.47216151d+01
COFD(7362) = 3.14480429d+00
COFD(7363) = -1.97906290d-01
COFD(7364) = 8.76325718d-03
COFD(7365) = -1.37578623d+01
COFD(7366) = 2.76469604d+00
COFD(7367) = -1.48843508d-01
COFD(7368) = 6.65045597d-03
COFD(7369) = -1.99702270d+01
COFD(7370) = 5.04961455d+00
COFD(7371) = -4.34856114d-01
COFD(7372) = 1.86234058d-02
COFD(7373) = -1.99796237d+01
COFD(7374) = 5.04961455d+00
COFD(7375) = -4.34856114d-01
COFD(7376) = 1.86234058d-02
COFD(7377) = -1.51581727d+01
COFD(7378) = 3.39229020d+00
COFD(7379) = -2.29520232d-01
COFD(7380) = 1.01114311d-02
COFD(7381) = -1.70117892d+01
COFD(7382) = 3.93849401d+00
COFD(7383) = -2.99416642d-01
COFD(7384) = 1.31020815d-02
COFD(7385) = -1.70117892d+01
COFD(7386) = 3.93849401d+00
COFD(7387) = -2.99416642d-01
COFD(7388) = 1.31020815d-02
COFD(7389) = -1.70117892d+01
COFD(7390) = 3.93849401d+00
COFD(7391) = -2.99416642d-01
COFD(7392) = 1.31020815d-02
COFD(7393) = -1.70069808d+01
COFD(7394) = 3.93849401d+00
COFD(7395) = -2.99416642d-01
COFD(7396) = 1.31020815d-02
COFD(7397) = -1.42697086d+01
COFD(7398) = 2.97137588d+00
COFD(7399) = -1.75491257d-01
COFD(7400) = 7.79646773d-03
COFD(7401) = -1.50911395d+01
COFD(7402) = 3.27505697d+00
COFD(7403) = -2.14306851d-01
COFD(7404) = 9.45219335d-03
COFD(7405) = -1.77780309d+01
COFD(7406) = 4.11838823d+00
COFD(7407) = -3.22329602d-01
COFD(7408) = 1.40802468d-02
COFD(7409) = -1.77826519d+01
COFD(7410) = 4.11838823d+00
COFD(7411) = -3.22329602d-01
COFD(7412) = 1.40802468d-02
COFD(7413) = -1.92109657d+01
COFD(7414) = 4.71729964d+00
COFD(7415) = -3.95432573d-01
COFD(7416) = 1.70623691d-02
COFD(7417) = -1.92155940d+01
COFD(7418) = 4.71729964d+00
COFD(7419) = -3.95432573d-01
COFD(7420) = 1.70623691d-02
COFD(7421) = -1.16928639d+01
COFD(7422) = 2.47469981d+00
COFD(7423) = -1.10436257d-01
COFD(7424) = 4.95273813d-03
COFD(7425) = -1.42905987d+01
COFD(7426) = 3.67490723d+00
COFD(7427) = -2.65114792d-01
COFD(7428) = 1.16092671d-02
COFD(7429) = -1.40879121d+01
COFD(7430) = 3.07549274d+00
COFD(7431) = -1.88889344d-01
COFD(7432) = 8.37152866d-03
COFD(7433) = -1.52594746d+01
COFD(7434) = 3.35922578d+00
COFD(7435) = -2.25181399d-01
COFD(7436) = 9.92132878d-03
COFD(7437) = -1.41076233d+01
COFD(7438) = 3.07549274d+00
COFD(7439) = -1.88889344d-01
COFD(7440) = 8.37152866d-03
COFD(7441) = -2.10774940d+01
COFD(7442) = 5.53614847d+00
COFD(7443) = -4.86046736d-01
COFD(7444) = 2.03659188d-02
COFD(7445) = -1.52669189d+01
COFD(7446) = 3.35922578d+00
COFD(7447) = -2.25181399d-01
COFD(7448) = 9.92132878d-03
COFD(7449) = -1.52740247d+01
COFD(7450) = 3.35922578d+00
COFD(7451) = -2.25181399d-01
COFD(7452) = 9.92132878d-03
COFD(7453) = -1.38762133d+01
COFD(7454) = 2.97137588d+00
COFD(7455) = -1.75491257d-01
COFD(7456) = 7.79646773d-03
COFD(7457) = -1.40183333d+01
COFD(7458) = 3.07549274d+00
COFD(7459) = -1.88889344d-01
COFD(7460) = 8.37152866d-03
COFD(7461) = -1.59516918d+01
COFD(7462) = 3.66853818d+00
COFD(7463) = -2.64346221d-01
COFD(7464) = 1.15784613d-02
COFD(7465) = -1.59516918d+01
COFD(7466) = 3.66853818d+00
COFD(7467) = -2.64346221d-01
COFD(7468) = 1.15784613d-02
COFD(7469) = -1.59750724d+01
COFD(7470) = 3.66853818d+00
COFD(7471) = -2.64346221d-01
COFD(7472) = 1.15784613d-02
COFD(7473) = -1.59449698d+01
COFD(7474) = 3.65620899d+00
COFD(7475) = -2.62933804d-01
COFD(7476) = 1.15253223d-02
COFD(7477) = -1.50200522d+01
COFD(7478) = 3.26223357d+00
COFD(7479) = -2.12746642d-01
COFD(7480) = 9.38912883d-03
COFD(7481) = -1.81639592d+01
COFD(7482) = 4.37565431d+00
COFD(7483) = -3.53906025d-01
COFD(7484) = 1.53760786d-02
COFD(7485) = -2.04922452d+01
COFD(7486) = 5.23112374d+00
COFD(7487) = -4.54967682d-01
COFD(7488) = 1.93570423d-02
COFD(7489) = -2.05008516d+01
COFD(7490) = 5.23112374d+00
COFD(7491) = -4.54967682d-01
COFD(7492) = 1.93570423d-02
COFD(7493) = -2.02446539d+01
COFD(7494) = 5.13632093d+00
COFD(7495) = -4.44839124d-01
COFD(7496) = 1.90058354d-02
COFD(7497) = -2.02446539d+01
COFD(7498) = 5.13632093d+00
COFD(7499) = -4.44839124d-01
COFD(7500) = 1.90058354d-02
COFD(7501) = -2.04024630d+01
COFD(7502) = 5.18856872d+00
COFD(7503) = -4.50001829d-01
COFD(7504) = 1.91636142d-02
COFD(7505) = -1.76057946d+01
COFD(7506) = 4.19171952d+00
COFD(7507) = -3.31354810d-01
COFD(7508) = 1.44520623d-02
COFD(7509) = -1.76164603d+01
COFD(7510) = 4.19171952d+00
COFD(7511) = -3.31354810d-01
COFD(7512) = 1.44520623d-02
COFD(7513) = -1.76265380d+01
COFD(7514) = 4.19171952d+00
COFD(7515) = -3.31354810d-01
COFD(7516) = 1.44520623d-02
COFD(7517) = -1.86033113d+01
COFD(7518) = 4.54915847d+00
COFD(7519) = -3.75000738d-01
COFD(7520) = 1.62324821d-02
COFD(7521) = -1.83338353d+01
COFD(7522) = 4.42045763d+00
COFD(7523) = -3.59451578d-01
COFD(7524) = 1.56056164d-02
COFD(7525) = -1.83424227d+01
COFD(7526) = 4.42045763d+00
COFD(7527) = -3.59451578d-01
COFD(7528) = 1.56056164d-02
COFD(7529) = -1.60085628d+01
COFD(7530) = 3.72220402d+00
COFD(7531) = -2.71150591d-01
COFD(7532) = 1.18665265d-02
COFD(7533) = -2.02849953d+01
COFD(7534) = 5.10426133d+00
COFD(7535) = -4.41256919d-01
COFD(7536) = 1.88737290d-02
COFD(7537) = -2.02849953d+01
COFD(7538) = 5.10426133d+00
COFD(7539) = -4.41256919d-01
COFD(7540) = 1.88737290d-02
COFD(7541) = -1.39298633d+01
COFD(7542) = 2.97137588d+00
COFD(7543) = -1.75491257d-01
COFD(7544) = 7.79646773d-03
COFD(7545) = -1.40353277d+01
COFD(7546) = 3.07549274d+00
COFD(7547) = -1.88889344d-01
COFD(7548) = 8.37152866d-03
COFD(7549) = -1.40567442d+01
COFD(7550) = 3.07549274d+00
COFD(7551) = -1.88889344d-01
COFD(7552) = 8.37152866d-03
COFD(7553) = -2.04444507d+01
COFD(7554) = 5.33236114d+00
COFD(7555) = -4.66969586d-01
COFD(7556) = 1.98332508d-02
COFD(7557) = -1.42868965d+01
COFD(7558) = 2.97137588d+00
COFD(7559) = -1.75491257d-01
COFD(7560) = 7.79646773d-03
COFD(7561) = -1.50172019d+01
COFD(7562) = 3.25781069d+00
COFD(7563) = -2.12199367d-01
COFD(7564) = 9.36657283d-03
COFD(7565) = -1.74483060d+01
COFD(7566) = 4.14166966d+00
COFD(7567) = -3.25149462d-01
COFD(7568) = 1.41943811d-02
COFD(7569) = -1.79972553d+01
COFD(7570) = 4.30841971d+00
COFD(7571) = -3.45524579d-01
COFD(7572) = 1.50265381d-02
COFD(7573) = -1.55475346d+01
COFD(7574) = 3.46766436d+00
COFD(7575) = -2.39228775d-01
COFD(7576) = 1.05291747d-02
COFD(7577) = -1.43825950d+01
COFD(7578) = 3.01144603d+00
COFD(7579) = -1.80597248d-01
COFD(7580) = 8.01333257d-03
COFD(7581) = -2.09155238d+01
COFD(7582) = 5.37089858d+00
COFD(7583) = -4.70314009d-01
COFD(7584) = 1.99113180d-02
COFD(7585) = -2.09250730d+01
COFD(7586) = 5.37089858d+00
COFD(7587) = -4.70314009d-01
COFD(7588) = 1.99113180d-02
COFD(7589) = -1.60085775d+01
COFD(7590) = 3.72220402d+00
COFD(7591) = -2.71150591d-01
COFD(7592) = 1.18665265d-02
COFD(7593) = -1.79926243d+01
COFD(7594) = 4.30841971d+00
COFD(7595) = -3.45524579d-01
COFD(7596) = 1.50265381d-02
COFD(7597) = -1.79926243d+01
COFD(7598) = 4.30841971d+00
COFD(7599) = -3.45524579d-01
COFD(7600) = 1.50265381d-02
COFD(7601) = -1.79926243d+01
COFD(7602) = 4.30841971d+00
COFD(7603) = -3.45524579d-01
COFD(7604) = 1.50265381d-02
COFD(7605) = -1.79877202d+01
COFD(7606) = 4.30841971d+00
COFD(7607) = -3.45524579d-01
COFD(7608) = 1.50265381d-02
COFD(7609) = -1.49997274d+01
COFD(7610) = 3.25781069d+00
COFD(7611) = -2.12199367d-01
COFD(7612) = 9.36657283d-03
COFD(7613) = -1.60076874d+01
COFD(7614) = 3.63340763d+00
COFD(7615) = -2.60307961d-01
COFD(7616) = 1.14256954d-02
COFD(7617) = -1.87735510d+01
COFD(7618) = 4.48550694d+00
COFD(7619) = -3.67277454d-01
COFD(7620) = 1.59194755d-02
COFD(7621) = -1.87782648d+01
COFD(7622) = 4.48550694d+00
COFD(7623) = -3.67277454d-01
COFD(7624) = 1.59194755d-02
COFD(7625) = -2.02898957d+01
COFD(7626) = 5.10426133d+00
COFD(7627) = -4.41256919d-01
COFD(7628) = 1.88737290d-02
COFD(7629) = -2.02946170d+01
COFD(7630) = 5.10426133d+00
COFD(7631) = -4.41256919d-01
COFD(7632) = 1.88737290d-02
COFD(7633) = -1.31911278d+01
COFD(7634) = 3.04970299d+00
COFD(7635) = -1.85555523d-01
COFD(7636) = 8.22773480d-03
COFD(7637) = -1.68734643d+01
COFD(7638) = 4.63687143d+00
COFD(7639) = -3.85900861d-01
COFD(7640) = 1.66856798d-02
COFD(7641) = -1.62772222d+01
COFD(7642) = 3.88250968d+00
COFD(7643) = -2.92155848d-01
COFD(7644) = 1.27867850d-02
COFD(7645) = -1.77337342d+01
COFD(7646) = 4.25438185d+00
COFD(7647) = -3.39084808d-01
COFD(7648) = 1.47709916d-02
COFD(7649) = -1.62997072d+01
COFD(7650) = 3.88250968d+00
COFD(7651) = -2.92155848d-01
COFD(7652) = 1.27867850d-02
COFD(7653) = -2.14542524d+01
COFD(7654) = 5.49993732d+00
COFD(7655) = -4.62042917d-01
COFD(7656) = 1.85577413d-02
COFD(7657) = -1.77428218d+01
COFD(7658) = 4.25438185d+00
COFD(7659) = -3.39084808d-01
COFD(7660) = 1.47709916d-02
COFD(7661) = -1.77515242d+01
COFD(7662) = 4.25438185d+00
COFD(7663) = -3.39084808d-01
COFD(7664) = 1.47709916d-02
COFD(7665) = -1.58184910d+01
COFD(7666) = 3.68407693d+00
COFD(7667) = -2.66228170d-01
COFD(7668) = 1.16542305d-02
COFD(7669) = -1.61987855d+01
COFD(7670) = 3.88250968d+00
COFD(7671) = -2.92155848d-01
COFD(7672) = 1.27867850d-02
COFD(7673) = -1.85866399d+01
COFD(7674) = 4.62613551d+00
COFD(7675) = -3.84556396d-01
COFD(7676) = 1.66292467d-02
COFD(7677) = -1.85866399d+01
COFD(7678) = 4.62613551d+00
COFD(7679) = -3.84556396d-01
COFD(7680) = 1.66292467d-02
COFD(7681) = -1.86130116d+01
COFD(7682) = 4.62613551d+00
COFD(7683) = -3.84556396d-01
COFD(7684) = 1.66292467d-02
COFD(7685) = -1.85418532d+01
COFD(7686) = 4.59643893d+00
COFD(7687) = -3.80823304d-01
COFD(7688) = 1.64720603d-02
COFD(7689) = -1.74527876d+01
COFD(7690) = 4.14792932d+00
COFD(7691) = -3.25920382d-01
COFD(7692) = 1.42261620d-02
COFD(7693) = -2.07349645d+01
COFD(7694) = 5.23705932d+00
COFD(7695) = -4.55655792d-01
COFD(7696) = 1.93836339d-02
COFD(7697) = -2.19815149d+01
COFD(7698) = 5.58446511d+00
COFD(7699) = -4.79399331d-01
COFD(7700) = 1.95652693d-02
COFD(7701) = -2.19919149d+01
COFD(7702) = 5.58446511d+00
COFD(7703) = -4.79399331d-01
COFD(7704) = 1.95652693d-02
COFD(7705) = -2.19270918d+01
COFD(7706) = 5.60963444d+00
COFD(7707) = -4.89805286d-01
COFD(7708) = 2.03017681d-02
COFD(7709) = -2.19270918d+01
COFD(7710) = 5.60963444d+00
COFD(7711) = -4.89805286d-01
COFD(7712) = 2.03017681d-02
COFD(7713) = -2.19932321d+01
COFD(7714) = 5.58514538d+00
COFD(7715) = -4.80745077d-01
COFD(7716) = 1.96733087d-02
COFD(7717) = -2.01897084d+01
COFD(7718) = 5.08408346d+00
COFD(7719) = -4.38907391d-01
COFD(7720) = 1.87824789d-02
COFD(7721) = -2.02024036d+01
COFD(7722) = 5.08408346d+00
COFD(7723) = -4.38907391d-01
COFD(7724) = 1.87824789d-02
COFD(7725) = -2.02144469d+01
COFD(7726) = 5.08408346d+00
COFD(7727) = -4.38907391d-01
COFD(7728) = 1.87824789d-02
COFD(7729) = -2.10998106d+01
COFD(7730) = 5.37657006d+00
COFD(7731) = -4.70707921d-01
COFD(7732) = 1.99145326d-02
COFD(7733) = -2.08992622d+01
COFD(7734) = 5.28495077d+00
COFD(7735) = -4.61326381d-01
COFD(7736) = 1.96080621d-02
COFD(7737) = -2.09096408d+01
COFD(7738) = 5.28495077d+00
COFD(7739) = -4.61326381d-01
COFD(7740) = 1.96080621d-02
COFD(7741) = -1.86961076d+01
COFD(7742) = 4.67999453d+00
COFD(7743) = -3.91135253d-01
COFD(7744) = 1.68982388d-02
COFD(7745) = -2.20848247d+01
COFD(7746) = 5.59135292d+00
COFD(7747) = -4.85565630d-01
COFD(7748) = 2.00429035d-02
COFD(7749) = -2.20848247d+01
COFD(7750) = 5.59135292d+00
COFD(7751) = -4.85565630d-01
COFD(7752) = 2.00429035d-02
COFD(7753) = -1.58784322d+01
COFD(7754) = 3.68407693d+00
COFD(7755) = -2.66228170d-01
COFD(7756) = 1.16542305d-02
COFD(7757) = -1.62212061d+01
COFD(7758) = 3.88250968d+00
COFD(7759) = -2.92155848d-01
COFD(7760) = 1.27867850d-02
COFD(7761) = -1.62455018d+01
COFD(7762) = 3.88250968d+00
COFD(7763) = -2.92155848d-01
COFD(7764) = 1.27867850d-02
COFD(7765) = -2.15739313d+01
COFD(7766) = 5.58518186d+00
COFD(7767) = -4.80813479d-01
COFD(7768) = 1.96787936d-02
COFD(7769) = -1.62729751d+01
COFD(7770) = 3.68407693d+00
COFD(7771) = -2.66228170d-01
COFD(7772) = 1.16542305d-02
COFD(7773) = -1.74483060d+01
COFD(7774) = 4.14166966d+00
COFD(7775) = -3.25149462d-01
COFD(7776) = 1.41943811d-02
COFD(7777) = -2.00443698d+01
COFD(7778) = 5.03042083d+00
COFD(7779) = -4.32596342d-01
COFD(7780) = 1.85345510d-02
COFD(7781) = -2.05716765d+01
COFD(7782) = 5.17526774d+00
COFD(7783) = -4.48472252d-01
COFD(7784) = 1.91050891d-02
COFD(7785) = -1.79999640d+01
COFD(7786) = 4.34871946d+00
COFD(7787) = -3.50542624d-01
COFD(7788) = 1.52355319d-02
COFD(7789) = -1.64757540d+01
COFD(7790) = 3.76873738d+00
COFD(7791) = -2.77251292d-01
COFD(7792) = 1.21337933d-02
COFD(7793) = -2.18991553d+01
COFD(7794) = 5.50125258d+00
COFD(7795) = -4.62458530d-01
COFD(7796) = 1.85853825d-02
COFD(7797) = -2.19106105d+01
COFD(7798) = 5.50125258d+00
COFD(7799) = -4.62458530d-01
COFD(7800) = 1.85853825d-02
COFD(7801) = -1.86961260d+01
COFD(7802) = 4.67999453d+00
COFD(7803) = -3.91135253d-01
COFD(7804) = 1.68982388d-02
COFD(7805) = -2.05658452d+01
COFD(7806) = 5.17526774d+00
COFD(7807) = -4.48472252d-01
COFD(7808) = 1.91050891d-02
COFD(7809) = -2.05658452d+01
COFD(7810) = 5.17526774d+00
COFD(7811) = -4.48472252d-01
COFD(7812) = 1.91050891d-02
COFD(7813) = -2.05658452d+01
COFD(7814) = 5.17526774d+00
COFD(7815) = -4.48472252d-01
COFD(7816) = 1.91050891d-02
COFD(7817) = -2.05596852d+01
COFD(7818) = 5.17526774d+00
COFD(7819) = -4.48472252d-01
COFD(7820) = 1.91050891d-02
COFD(7821) = -1.74272299d+01
COFD(7822) = 4.14166966d+00
COFD(7823) = -3.25149462d-01
COFD(7824) = 1.41943811d-02
COFD(7825) = -1.85744661d+01
COFD(7826) = 4.54507313d+00
COFD(7827) = -3.74508073d-01
COFD(7828) = 1.62126770d-02
COFD(7829) = -2.13758293d+01
COFD(7830) = 5.35135690d+00
COFD(7831) = -4.68903699d-01
COFD(7832) = 1.98958921d-02
COFD(7833) = -2.13817659d+01
COFD(7834) = 5.35135690d+00
COFD(7835) = -4.68903699d-01
COFD(7836) = 1.98958921d-02
COFD(7837) = -2.20909803d+01
COFD(7838) = 5.59135292d+00
COFD(7839) = -4.85565630d-01
COFD(7840) = 2.00429035d-02
COFD(7841) = -2.20969258d+01
COFD(7842) = 5.59135292d+00
COFD(7843) = -4.85565630d-01
COFD(7844) = 2.00429035d-02
COFD(7845) = -1.36894329d+01
COFD(7846) = 3.19999740d+00
COFD(7847) = -2.04999020d-01
COFD(7848) = 9.06766105d-03
COFD(7849) = -1.74819100d+01
COFD(7850) = 4.80792005d+00
COFD(7851) = -4.06126584d-01
COFD(7852) = 1.74831083d-02
COFD(7853) = -1.69290095d+01
COFD(7854) = 4.09077642d+00
COFD(7855) = -3.18894990d-01
COFD(7856) = 1.39371445d-02
COFD(7857) = -1.83308393d+01
COFD(7858) = 4.43878381d+00
COFD(7859) = -3.61697624d-01
COFD(7860) = 1.56975581d-02
COFD(7861) = -1.69512290d+01
COFD(7862) = 4.09077642d+00
COFD(7863) = -3.18894990d-01
COFD(7864) = 1.39371445d-02
COFD(7865) = -2.13656038d+01
COFD(7866) = 5.35901292d+00
COFD(7867) = -4.36172487d-01
COFD(7868) = 1.71345319d-02
COFD(7869) = -1.83397596d+01
COFD(7870) = 4.43878381d+00
COFD(7871) = -3.61697624d-01
COFD(7872) = 1.56975581d-02
COFD(7873) = -1.83482991d+01
COFD(7874) = 4.43878381d+00
COFD(7875) = -3.61697624d-01
COFD(7876) = 1.56975581d-02
COFD(7877) = -1.65721778d+01
COFD(7878) = 3.93849401d+00
COFD(7879) = -2.99416642d-01
COFD(7880) = 1.31020815d-02
COFD(7881) = -1.68514118d+01
COFD(7882) = 4.09077642d+00
COFD(7883) = -3.18894990d-01
COFD(7884) = 1.39371445d-02
COFD(7885) = -1.91704234d+01
COFD(7886) = 4.80030220d+00
COFD(7887) = -4.05235041d-01
COFD(7888) = 1.74483531d-02
COFD(7889) = -1.91704234d+01
COFD(7890) = 4.80030220d+00
COFD(7891) = -4.05235041d-01
COFD(7892) = 1.74483531d-02
COFD(7893) = -1.91965117d+01
COFD(7894) = 4.80030220d+00
COFD(7895) = -4.05235041d-01
COFD(7896) = 1.74483531d-02
COFD(7897) = -1.91494015d+01
COFD(7898) = 4.78094221d+00
COFD(7899) = -4.02985837d-01
COFD(7900) = 1.73614221d-02
COFD(7901) = -1.80065982d+01
COFD(7902) = 4.31656593d+00
COFD(7903) = -3.46539554d-01
COFD(7904) = 1.50688196d-02
COFD(7905) = -2.12671597d+01
COFD(7906) = 5.38135645d+00
COFD(7907) = -4.71058360d-01
COFD(7908) = 1.99188046d-02
COFD(7909) = -2.19960393d+01
COFD(7910) = 5.49642957d+00
COFD(7911) = -4.61132993d-01
COFD(7912) = 1.85004773d-02
COFD(7913) = -2.20062584d+01
COFD(7914) = 5.49642957d+00
COFD(7915) = -4.61132993d-01
COFD(7916) = 1.85004773d-02
COFD(7917) = -2.20870532d+01
COFD(7918) = 5.58518321d+00
COFD(7919) = -4.80534235d-01
COFD(7920) = 1.96556393d-02
COFD(7921) = -2.20870532d+01
COFD(7922) = 5.58518321d+00
COFD(7923) = -4.80534235d-01
COFD(7924) = 1.96556393d-02
COFD(7925) = -2.20260694d+01
COFD(7926) = 5.50583166d+00
COFD(7927) = -4.63753262d-01
COFD(7928) = 1.86693462d-02
COFD(7929) = -2.07211857d+01
COFD(7930) = 5.23116678d+00
COFD(7931) = -4.54972675d-01
COFD(7932) = 1.93572354d-02
COFD(7933) = -2.07336791d+01
COFD(7934) = 5.23116678d+00
COFD(7935) = -4.54972675d-01
COFD(7936) = 1.93572354d-02
COFD(7937) = -2.07455261d+01
COFD(7938) = 5.23116678d+00
COFD(7939) = -4.54972675d-01
COFD(7940) = 1.93572354d-02
COFD(7941) = -2.15315009d+01
COFD(7942) = 5.47662534d+00
COFD(7943) = -4.79750189d-01
COFD(7944) = 2.01492627d-02
COFD(7945) = -2.13547906d+01
COFD(7946) = 5.39784038d+00
COFD(7947) = -4.72269379d-01
COFD(7948) = 1.99339592d-02
COFD(7949) = -2.13649886d+01
COFD(7950) = 5.39784038d+00
COFD(7951) = -4.72269379d-01
COFD(7952) = 1.99339592d-02
COFD(7953) = -1.92693772d+01
COFD(7954) = 4.84335373d+00
COFD(7955) = -4.10212751d-01
COFD(7956) = 1.76395874d-02
COFD(7957) = -2.22703814d+01
COFD(7958) = 5.57944619d+00
COFD(7959) = -4.78032626d-01
COFD(7960) = 1.94775370d-02
COFD(7961) = -2.22703814d+01
COFD(7962) = 5.57944619d+00
COFD(7963) = -4.78032626d-01
COFD(7964) = 1.94775370d-02
COFD(7965) = -1.66315280d+01
COFD(7966) = 3.93849401d+00
COFD(7967) = -2.99416642d-01
COFD(7968) = 1.31020815d-02
COFD(7969) = -1.68748824d+01
COFD(7970) = 4.09077642d+00
COFD(7971) = -3.18894990d-01
COFD(7972) = 1.39371445d-02
COFD(7973) = -1.68989038d+01
COFD(7974) = 4.09077642d+00
COFD(7975) = -3.18894990d-01
COFD(7976) = 1.39371445d-02
COFD(7977) = -2.16211036d+01
COFD(7978) = 5.50657794d+00
COFD(7979) = -4.63927523d-01
COFD(7980) = 1.86800112d-02
COFD(7981) = -1.70163290d+01
COFD(7982) = 3.93849401d+00
COFD(7983) = -2.99416642d-01
COFD(7984) = 1.31020815d-02
COFD(7985) = -1.79972553d+01
COFD(7986) = 4.30841971d+00
COFD(7987) = -3.45524579d-01
COFD(7988) = 1.50265381d-02
COFD(7989) = -2.05716765d+01
COFD(7990) = 5.17526774d+00
COFD(7991) = -4.48472252d-01
COFD(7992) = 1.91050891d-02
COFD(7993) = -2.11909744d+01
COFD(7994) = 5.35817304d+00
COFD(7995) = -4.69455306d-01
COFD(7996) = 1.99063292d-02
COFD(7997) = -1.86088628d+01
COFD(7998) = 4.53681554d+00
COFD(7999) = -3.73517224d-01
COFD(8000) = 1.61730724d-02
COFD(8001) = -1.71848523d+01
COFD(8002) = 4.00896225d+00
COFD(8003) = -3.08501394d-01
COFD(8004) = 1.34944985d-02
COFD(8005) = -2.18091268d+01
COFD(8006) = 5.36691865d+00
COFD(8007) = -4.37522912d-01
COFD(8008) = 1.72059807d-02
COFD(8009) = -2.18203911d+01
COFD(8010) = 5.36691865d+00
COFD(8011) = -4.37522912d-01
COFD(8012) = 1.72059807d-02
COFD(8013) = -1.92693952d+01
COFD(8014) = 4.84335373d+00
COFD(8015) = -4.10212751d-01
COFD(8016) = 1.76395874d-02
COFD(8017) = -2.11852687d+01
COFD(8018) = 5.35817304d+00
COFD(8019) = -4.69455306d-01
COFD(8020) = 1.99063292d-02
COFD(8021) = -2.11852687d+01
COFD(8022) = 5.35817304d+00
COFD(8023) = -4.69455306d-01
COFD(8024) = 1.99063292d-02
COFD(8025) = -2.11852687d+01
COFD(8026) = 5.35817304d+00
COFD(8027) = -4.69455306d-01
COFD(8028) = 1.99063292d-02
COFD(8029) = -2.11792397d+01
COFD(8030) = 5.35817304d+00
COFD(8031) = -4.69455306d-01
COFD(8032) = 1.99063292d-02
COFD(8033) = -1.79765417d+01
COFD(8034) = 4.30841971d+00
COFD(8035) = -3.45524579d-01
COFD(8036) = 1.50265381d-02
COFD(8037) = -1.92029732d+01
COFD(8038) = 4.73837374d+00
COFD(8039) = -3.97916523d-01
COFD(8040) = 1.71599378d-02
COFD(8041) = -2.17190071d+01
COFD(8042) = 5.41948945d+00
COFD(8043) = -4.73474138d-01
COFD(8044) = 1.99217500d-02
COFD(8045) = -2.17248157d+01
COFD(8046) = 5.41948945d+00
COFD(8047) = -4.73474138d-01
COFD(8048) = 1.99217500d-02
COFD(8049) = -2.22764060d+01
COFD(8050) = 5.57944619d+00
COFD(8051) = -4.78032626d-01
COFD(8052) = 1.94775370d-02
COFD(8053) = -2.22822234d+01
COFD(8054) = 5.57944619d+00
COFD(8055) = -4.78032626d-01
COFD(8056) = 1.94775370d-02
COFD(8057) = -1.20867697d+01
COFD(8058) = 2.64389960d+00
COFD(8059) = -1.33241706d-01
COFD(8060) = 5.97810200d-03
COFD(8061) = -1.50140778d+01
COFD(8062) = 3.96853403d+00
COFD(8063) = -3.03320126d-01
COFD(8064) = 1.32719819d-02
COFD(8065) = -1.45065606d+01
COFD(8066) = 3.24450689d+00
COFD(8067) = -2.10570734d-01
COFD(8068) = 9.30026771d-03
COFD(8069) = -1.57967324d+01
COFD(8070) = 3.57122377d+00
COFD(8071) = -2.52409987d-01
COFD(8072) = 1.10900562d-02
COFD(8073) = -1.45265017d+01
COFD(8074) = 3.24450689d+00
COFD(8075) = -2.10570734d-01
COFD(8076) = 9.30026771d-03
COFD(8077) = -2.09199029d+01
COFD(8078) = 5.50410224d+00
COFD(8079) = -4.82760329d-01
COFD(8080) = 2.02578587d-02
COFD(8081) = -1.58043047d+01
COFD(8082) = 3.57122377d+00
COFD(8083) = -2.52409987d-01
COFD(8084) = 1.10900562d-02
COFD(8085) = -1.58115345d+01
COFD(8086) = 3.57122377d+00
COFD(8087) = -2.52409987d-01
COFD(8088) = 1.10900562d-02
COFD(8089) = -1.43049857d+01
COFD(8090) = 3.14480429d+00
COFD(8091) = -1.97906290d-01
COFD(8092) = 8.76325718d-03
COFD(8093) = -1.44362397d+01
COFD(8094) = 3.24450689d+00
COFD(8095) = -2.10570734d-01
COFD(8096) = 9.30026771d-03
COFD(8097) = -1.66814763d+01
COFD(8098) = 3.95813366d+00
COFD(8099) = -3.01970405d-01
COFD(8100) = 1.32132947d-02
COFD(8101) = -1.66814763d+01
COFD(8102) = 3.95813366d+00
COFD(8103) = -3.01970405d-01
COFD(8104) = 1.32132947d-02
COFD(8105) = -1.67051074d+01
COFD(8106) = 3.95813366d+00
COFD(8107) = -3.01970405d-01
COFD(8108) = 1.32132947d-02
COFD(8109) = -1.66377559d+01
COFD(8110) = 3.92993514d+00
COFD(8111) = -2.98305855d-01
COFD(8112) = 1.30538105d-02
COFD(8113) = -1.55536213d+01
COFD(8114) = 3.47339631d+00
COFD(8115) = -2.39946241d-01
COFD(8116) = 1.05591452d-02
COFD(8117) = -1.87873369d+01
COFD(8118) = 4.60768689d+00
COFD(8119) = -3.82235384d-01
COFD(8120) = 1.65314065d-02
COFD(8121) = -2.09853721d+01
COFD(8122) = 5.39401640d+00
COFD(8123) = -4.72026580d-01
COFD(8124) = 1.99336592d-02
COFD(8125) = -2.09941194d+01
COFD(8126) = 5.39401640d+00
COFD(8127) = -4.72026580d-01
COFD(8128) = 1.99336592d-02
COFD(8129) = -2.05164480d+01
COFD(8130) = 5.23355494d+00
COFD(8131) = -4.55249645d-01
COFD(8132) = 1.93679438d-02
COFD(8133) = -2.05164480d+01
COFD(8134) = 5.23355494d+00
COFD(8135) = -4.55249645d-01
COFD(8136) = 1.93679438d-02
COFD(8137) = -2.09549117d+01
COFD(8138) = 5.37710408d+00
COFD(8139) = -4.70743467d-01
COFD(8140) = 1.99147079d-02
COFD(8141) = -1.81843740d+01
COFD(8142) = 4.40883885d+00
COFD(8143) = -3.58021243d-01
COFD(8144) = 1.55467750d-02
COFD(8145) = -1.81952014d+01
COFD(8146) = 4.40883885d+00
COFD(8147) = -3.58021243d-01
COFD(8148) = 1.55467750d-02
COFD(8149) = -1.82054352d+01
COFD(8150) = 4.40883885d+00
COFD(8151) = -3.58021243d-01
COFD(8152) = 1.55467750d-02
COFD(8153) = -1.92234721d+01
COFD(8154) = 4.77781979d+00
COFD(8155) = -4.02620572d-01
COFD(8156) = 1.73472054d-02
COFD(8157) = -1.89704705d+01
COFD(8158) = 4.65726073d+00
COFD(8159) = -3.88398872d-01
COFD(8160) = 1.67881413d-02
COFD(8161) = -1.89791987d+01
COFD(8162) = 4.65726073d+00
COFD(8163) = -3.88398872d-01
COFD(8164) = 1.67881413d-02
COFD(8165) = -1.67419735d+01
COFD(8166) = 4.01496916d+00
COFD(8167) = -3.09274229d-01
COFD(8168) = 1.35278248d-02
COFD(8169) = -2.08481163d+01
COFD(8170) = 5.29805436d+00
COFD(8171) = -4.62913371d-01
COFD(8172) = 1.96725621d-02
COFD(8173) = -2.08481163d+01
COFD(8174) = 5.29805436d+00
COFD(8175) = -4.62913371d-01
COFD(8176) = 1.96725621d-02
COFD(8177) = -1.43591671d+01
COFD(8178) = 3.14480429d+00
COFD(8179) = -1.97906290d-01
COFD(8180) = 8.76325718d-03
COFD(8181) = -1.44530825d+01
COFD(8182) = 3.24450689d+00
COFD(8183) = -2.10570734d-01
COFD(8184) = 9.30026771d-03
COFD(8185) = -1.44747388d+01
COFD(8186) = 3.24450689d+00
COFD(8187) = -2.10570734d-01
COFD(8188) = 9.30026771d-03
COFD(8189) = -2.05628395d+01
COFD(8190) = 5.37634919d+00
COFD(8191) = -4.70693222d-01
COFD(8192) = 1.99144603d-02
COFD(8193) = -1.47216151d+01
COFD(8194) = 3.14480429d+00
COFD(8195) = -1.97906290d-01
COFD(8196) = 8.76325718d-03
COFD(8197) = -1.55475346d+01
COFD(8198) = 3.46766436d+00
COFD(8199) = -2.39228775d-01
COFD(8200) = 1.05291747d-02
COFD(8201) = -1.79999640d+01
COFD(8202) = 4.34871946d+00
COFD(8203) = -3.50542624d-01
COFD(8204) = 1.52355319d-02
COFD(8205) = -1.86088628d+01
COFD(8206) = 4.53681554d+00
COFD(8207) = -3.73517224d-01
COFD(8208) = 1.61730724d-02
COFD(8209) = -1.60107192d+01
COFD(8210) = 3.64795923d+00
COFD(8211) = -2.61984156d-01
COFD(8212) = 1.14893124d-02
COFD(8213) = -1.48360101d+01
COFD(8214) = 3.19255389d+00
COFD(8215) = -2.04063463d-01
COFD(8216) = 9.02844862d-03
COFD(8217) = -2.13278787d+01
COFD(8218) = 5.49726595d+00
COFD(8219) = -4.82018917d-01
COFD(8220) = 2.02314792d-02
COFD(8221) = -2.13375787d+01
COFD(8222) = 5.49726595d+00
COFD(8223) = -4.82018917d-01
COFD(8224) = 2.02314792d-02
COFD(8225) = -1.67419884d+01
COFD(8226) = 4.01496916d+00
COFD(8227) = -3.09274229d-01
COFD(8228) = 1.35278248d-02
COFD(8229) = -1.86041410d+01
COFD(8230) = 4.53681554d+00
COFD(8231) = -3.73517224d-01
COFD(8232) = 1.61730724d-02
COFD(8233) = -1.86041410d+01
COFD(8234) = 4.53681554d+00
COFD(8235) = -3.73517224d-01
COFD(8236) = 1.61730724d-02
COFD(8237) = -1.86041410d+01
COFD(8238) = 4.53681554d+00
COFD(8239) = -3.73517224d-01
COFD(8240) = 1.61730724d-02
COFD(8241) = -1.85991417d+01
COFD(8242) = 4.53681554d+00
COFD(8243) = -3.73517224d-01
COFD(8244) = 1.61730724d-02
COFD(8245) = -1.55297765d+01
COFD(8246) = 3.46766436d+00
COFD(8247) = -2.39228775d-01
COFD(8248) = 1.05291747d-02
COFD(8249) = -1.66272428d+01
COFD(8250) = 3.87567029d+00
COFD(8251) = -2.91269312d-01
COFD(8252) = 1.27483193d-02
COFD(8253) = -1.94052185d+01
COFD(8254) = 4.71746683d+00
COFD(8255) = -3.95451856d-01
COFD(8256) = 1.70631070d-02
COFD(8257) = -1.94100248d+01
COFD(8258) = 4.71746683d+00
COFD(8259) = -3.95451856d-01
COFD(8260) = 1.70631070d-02
COFD(8261) = -2.08531119d+01
COFD(8262) = 5.29805436d+00
COFD(8263) = -4.62913371d-01
COFD(8264) = 1.96725621d-02
COFD(8265) = -2.08579258d+01
COFD(8266) = 5.29805436d+00
COFD(8267) = -4.62913371d-01
COFD(8268) = 1.96725621d-02
COFD(8269) = -1.11796892d+01
COFD(8270) = 2.24409765d+00
COFD(8271) = -7.86374012d-02
COFD(8272) = 3.48840768d-03
COFD(8273) = -1.36558633d+01
COFD(8274) = 3.41450296d+00
COFD(8275) = -2.32431042d-01
COFD(8276) = 1.02388399d-02
COFD(8277) = -1.35403054d+01
COFD(8278) = 2.85448298d+00
COFD(8279) = -1.60491387d-01
COFD(8280) = 7.15451323d-03
COFD(8281) = -1.46137599d+01
COFD(8282) = 3.10987301d+00
COFD(8283) = -1.93373245d-01
COFD(8284) = 8.56678670d-03
COFD(8285) = -1.35590001d+01
COFD(8286) = 2.85448298d+00
COFD(8287) = -1.60491387d-01
COFD(8288) = 7.15451323d-03
COFD(8289) = -1.97463680d+01
COFD(8290) = 5.11531505d+00
COFD(8291) = -4.42537303d-01
COFD(8292) = 1.89231395d-02
COFD(8293) = -1.46206548d+01
COFD(8294) = 3.10987301d+00
COFD(8295) = -1.93373245d-01
COFD(8296) = 8.56678670d-03
COFD(8297) = -1.46272289d+01
COFD(8298) = 3.10987301d+00
COFD(8299) = -1.93373245d-01
COFD(8300) = 8.56678670d-03
COFD(8301) = -1.33664987d+01
COFD(8302) = 2.76469604d+00
COFD(8303) = -1.48843508d-01
COFD(8304) = 6.65045597d-03
COFD(8305) = -1.34740249d+01
COFD(8306) = 2.85448298d+00
COFD(8307) = -1.60491387d-01
COFD(8308) = 7.15451323d-03
COFD(8309) = -1.52696275d+01
COFD(8310) = 3.40487615d+00
COFD(8311) = -2.31175464d-01
COFD(8312) = 1.01841353d-02
COFD(8313) = -1.52696275d+01
COFD(8314) = 3.40487615d+00
COFD(8315) = -2.31175464d-01
COFD(8316) = 1.01841353d-02
COFD(8317) = -1.52918950d+01
COFD(8318) = 3.40487615d+00
COFD(8319) = -2.31175464d-01
COFD(8320) = 1.01841353d-02
COFD(8321) = -1.52271026d+01
COFD(8322) = 3.37751537d+00
COFD(8323) = -2.27579233d-01
COFD(8324) = 1.00262778d-02
COFD(8325) = -1.43882477d+01
COFD(8326) = 3.01675183d+00
COFD(8327) = -1.81281138d-01
COFD(8328) = 8.04273850d-03
COFD(8329) = -1.73328509d+01
COFD(8330) = 4.06993115d+00
COFD(8331) = -3.16228805d-01
COFD(8332) = 1.38227651d-02
COFD(8333) = -1.96553384d+01
COFD(8334) = 4.94225344d+00
COFD(8335) = -4.22163426d-01
COFD(8336) = 1.81223730d-02
COFD(8337) = -1.96633368d+01
COFD(8338) = 4.94225344d+00
COFD(8339) = -4.22163426d-01
COFD(8340) = 1.81223730d-02
COFD(8341) = -1.90658825d+01
COFD(8342) = 4.72218939d+00
COFD(8343) = -3.96003719d-01
COFD(8344) = 1.70845726d-02
COFD(8345) = -1.90658825d+01
COFD(8346) = 4.72218939d+00
COFD(8347) = -3.96003719d-01
COFD(8348) = 1.70845726d-02
COFD(8349) = -1.95512940d+01
COFD(8350) = 4.89369650d+00
COFD(8351) = -4.16273855d-01
COFD(8352) = 1.78833908d-02
COFD(8353) = -1.67047589d+01
COFD(8354) = 3.84749999d+00
COFD(8355) = -2.87575310d-01
COFD(8356) = 1.25862768d-02
COFD(8357) = -1.67147225d+01
COFD(8358) = 3.84749999d+00
COFD(8359) = -2.87575310d-01
COFD(8360) = 1.25862768d-02
COFD(8361) = -1.67241238d+01
COFD(8362) = 3.84749999d+00
COFD(8363) = -2.87575310d-01
COFD(8364) = 1.25862768d-02
COFD(8365) = -1.77393195d+01
COFD(8366) = 4.23219154d+00
COFD(8367) = -3.36374352d-01
COFD(8368) = 1.46603433d-02
COFD(8369) = -1.74955909d+01
COFD(8370) = 4.11180810d+00
COFD(8371) = -3.21531431d-01
COFD(8372) = 1.40478729d-02
COFD(8373) = -1.75035712d+01
COFD(8374) = 4.11180810d+00
COFD(8375) = -3.21531431d-01
COFD(8376) = 1.40478729d-02
COFD(8377) = -1.53109350d+01
COFD(8378) = 3.45550279d+00
COFD(8379) = -2.37693211d-01
COFD(8380) = 1.04644205d-02
COFD(8381) = -1.93551579d+01
COFD(8382) = 4.77546552d+00
COFD(8383) = -4.02345606d-01
COFD(8384) = 1.73365265d-02
COFD(8385) = -1.93551579d+01
COFD(8386) = 4.77546552d+00
COFD(8387) = -4.02345606d-01
COFD(8388) = 1.73365265d-02
COFD(8389) = -1.34177769d+01
COFD(8390) = 2.76469604d+00
COFD(8391) = -1.48843508d-01
COFD(8392) = 6.65045597d-03
COFD(8393) = -1.34898944d+01
COFD(8394) = 2.85448298d+00
COFD(8395) = -1.60491387d-01
COFD(8396) = 7.15451323d-03
COFD(8397) = -1.35102477d+01
COFD(8398) = 2.85448298d+00
COFD(8399) = -1.60491387d-01
COFD(8400) = 7.15451323d-03
COFD(8401) = -1.91791848d+01
COFD(8402) = 4.89100437d+00
COFD(8403) = -4.15941503d-01
COFD(8404) = 1.78696440d-02
COFD(8405) = -1.37578623d+01
COFD(8406) = 2.76469604d+00
COFD(8407) = -1.48843508d-01
COFD(8408) = 6.65045597d-03
COFD(8409) = -1.43825950d+01
COFD(8410) = 3.01144603d+00
COFD(8411) = -1.80597248d-01
COFD(8412) = 8.01333257d-03
COFD(8413) = -1.64757540d+01
COFD(8414) = 3.76873738d+00
COFD(8415) = -2.77251292d-01
COFD(8416) = 1.21337933d-02
COFD(8417) = -1.71848523d+01
COFD(8418) = 4.00896225d+00
COFD(8419) = -3.08501394d-01
COFD(8420) = 1.34944985d-02
COFD(8421) = -1.48360101d+01
COFD(8422) = 3.19255389d+00
COFD(8423) = -2.04063463d-01
COFD(8424) = 9.02844862d-03
COFD(8425) = -1.38288862d+01
COFD(8426) = 2.79461504d+00
COFD(8427) = -1.52633727d-01
COFD(8428) = 6.81041543d-03
COFD(8429) = -2.01303629d+01
COFD(8430) = 5.10847039d+00
COFD(8431) = -4.41746871d-01
COFD(8432) = 1.88927540d-02
COFD(8433) = -2.01392596d+01
COFD(8434) = 5.10847039d+00
COFD(8435) = -4.41746871d-01
COFD(8436) = 1.88927540d-02
COFD(8437) = -1.53109485d+01
COFD(8438) = 3.45550279d+00
COFD(8439) = -2.37693211d-01
COFD(8440) = 1.04644205d-02
COFD(8441) = -1.71806065d+01
COFD(8442) = 4.00896225d+00
COFD(8443) = -3.08501394d-01
COFD(8444) = 1.34944985d-02
COFD(8445) = -1.71806065d+01
COFD(8446) = 4.00896225d+00
COFD(8447) = -3.08501394d-01
COFD(8448) = 1.34944985d-02
COFD(8449) = -1.71806065d+01
COFD(8450) = 4.00896225d+00
COFD(8451) = -3.08501394d-01
COFD(8452) = 1.34944985d-02
COFD(8453) = -1.71761068d+01
COFD(8454) = 4.00896225d+00
COFD(8455) = -3.08501394d-01
COFD(8456) = 1.34944985d-02
COFD(8457) = -1.43663445d+01
COFD(8458) = 3.01144603d+00
COFD(8459) = -1.80597248d-01
COFD(8460) = 8.01333257d-03
COFD(8461) = -1.52099661d+01
COFD(8462) = 3.32517045d+00
COFD(8463) = -2.20697887d-01
COFD(8464) = 9.72421963d-03
COFD(8465) = -1.79018076d+01
COFD(8466) = 4.16975222d+00
COFD(8467) = -3.28627847d-01
COFD(8468) = 1.43387221d-02
COFD(8469) = -1.79061291d+01
COFD(8470) = 4.16975222d+00
COFD(8471) = -3.28627847d-01
COFD(8472) = 1.43387221d-02
COFD(8473) = -1.93596541d+01
COFD(8474) = 4.77546552d+00
COFD(8475) = -4.02345606d-01
COFD(8476) = 1.73365265d-02
COFD(8477) = -1.93639826d+01
COFD(8478) = 4.77546552d+00
COFD(8479) = -4.02345606d-01
COFD(8480) = 1.73365265d-02
COFD(8481) = -1.64793369d+01
COFD(8482) = 4.26176920d+00
COFD(8483) = -3.39986562d-01
COFD(8484) = 1.48077740d-02
COFD(8485) = -1.99775614d+01
COFD(8486) = 5.61184673d+00
COFD(8487) = -4.90499758d-01
COFD(8488) = 2.03480832d-02
COFD(8489) = -1.98343743d+01
COFD(8490) = 5.15754807d+00
COFD(8491) = -4.46654084d-01
COFD(8492) = 1.90458977d-02
COFD(8493) = -2.10651573d+01
COFD(8494) = 5.41544337d+00
COFD(8495) = -4.73380953d-01
COFD(8496) = 1.99350119d-02
COFD(8497) = -1.98533435d+01
COFD(8498) = 5.15754807d+00
COFD(8499) = -4.46654084d-01
COFD(8500) = 1.90458977d-02
COFD(8501) = -1.66996477d+01
COFD(8502) = 3.08217426d+00
COFD(8503) = -8.33350945d-02
COFD(8504) = -3.89962983d-04
COFD(8505) = -2.10721980d+01
COFD(8506) = 5.41544337d+00
COFD(8507) = -4.73380953d-01
COFD(8508) = 1.99350119d-02
COFD(8509) = -2.10789131d+01
COFD(8510) = 5.41544337d+00
COFD(8511) = -4.73380953d-01
COFD(8512) = 1.99350119d-02
COFD(8513) = -1.95705999d+01
COFD(8514) = 5.04961455d+00
COFD(8515) = -4.34856114d-01
COFD(8516) = 1.86234058d-02
COFD(8517) = -1.97672003d+01
COFD(8518) = 5.15754807d+00
COFD(8519) = -4.46654084d-01
COFD(8520) = 1.90458977d-02
COFD(8521) = -2.16351552d+01
COFD(8522) = 5.61201619d+00
COFD(8523) = -4.90763024d-01
COFD(8524) = 2.03690162d-02
COFD(8525) = -2.16351552d+01
COFD(8526) = 5.61201619d+00
COFD(8527) = -4.90763024d-01
COFD(8528) = 2.03690162d-02
COFD(8529) = -2.16577242d+01
COFD(8530) = 5.61201619d+00
COFD(8531) = -4.90763024d-01
COFD(8532) = 2.03690162d-02
COFD(8533) = -2.16493228d+01
COFD(8534) = 5.61229743d+00
COFD(8535) = -4.91428747d-01
COFD(8536) = 2.04226556d-02
COFD(8537) = -2.09183945d+01
COFD(8538) = 5.37381486d+00
COFD(8539) = -4.70526046d-01
COFD(8540) = 1.99137630d-02
COFD(8541) = -2.16858789d+01
COFD(8542) = 5.30214979d+00
COFD(8543) = -4.26479693d-01
COFD(8544) = 1.66224918d-02
COFD(8545) = -1.82195639d+01
COFD(8546) = 3.59760627d+00
COFD(8547) = -1.59809307d-01
COFD(8548) = 3.28676045d-03
COFD(8549) = -1.82277241d+01
COFD(8550) = 3.59760627d+00
COFD(8551) = -1.59809307d-01
COFD(8552) = 3.28676045d-03
COFD(8553) = -1.93898781d+01
COFD(8554) = 4.14577053d+00
COFD(8555) = -2.43003893d-01
COFD(8556) = 7.35638742d-03
COFD(8557) = -1.93898781d+01
COFD(8558) = 4.14577053d+00
COFD(8559) = -2.43003893d-01
COFD(8560) = 7.35638742d-03
COFD(8561) = -1.84848887d+01
COFD(8562) = 3.70911349d+00
COFD(8563) = -1.76610648d-01
COFD(8564) = 4.10440582d-03
COFD(8565) = -2.19022500d+01
COFD(8566) = 5.47350698d+00
COFD(8567) = -4.56806529d-01
COFD(8568) = 1.82590273d-02
COFD(8569) = -2.19124011d+01
COFD(8570) = 5.47350698d+00
COFD(8571) = -4.56806529d-01
COFD(8572) = 1.82590273d-02
COFD(8573) = -2.19219828d+01
COFD(8574) = 5.47350698d+00
COFD(8575) = -4.56806529d-01
COFD(8576) = 1.82590273d-02
COFD(8577) = -2.13056183d+01
COFD(8578) = 5.10013279d+00
COFD(8579) = -3.92393844d-01
COFD(8580) = 1.48333781d-02
COFD(8581) = -2.16618765d+01
COFD(8582) = 5.25608421d+00
COFD(8583) = -4.18625448d-01
COFD(8584) = 1.62072845d-02
COFD(8585) = -2.16700183d+01
COFD(8586) = 5.25608421d+00
COFD(8587) = -4.18625448d-01
COFD(8588) = 1.62072845d-02
COFD(8589) = -2.15423744d+01
COFD(8590) = 5.60052862d+00
COFD(8591) = -4.87676563d-01
COFD(8592) = 2.01715067d-02
COFD(8593) = -1.92068539d+01
COFD(8594) = 3.99103591d+00
COFD(8595) = -2.19677673d-01
COFD(8596) = 6.21875426d-03
COFD(8597) = -1.92068539d+01
COFD(8598) = 3.99103591d+00
COFD(8599) = -2.19677673d-01
COFD(8600) = 6.21875426d-03
COFD(8601) = -1.96225223d+01
COFD(8602) = 5.04961455d+00
COFD(8603) = -4.34856114d-01
COFD(8604) = 1.86234058d-02
COFD(8605) = -1.97825930d+01
COFD(8606) = 5.15754807d+00
COFD(8607) = -4.46654084d-01
COFD(8608) = 1.90458977d-02
COFD(8609) = -1.98032338d+01
COFD(8610) = 5.15754807d+00
COFD(8611) = -4.46654084d-01
COFD(8612) = 1.90458977d-02
COFD(8613) = -1.81204154d+01
COFD(8614) = 3.71386126d+00
COFD(8615) = -1.77333629d-01
COFD(8616) = 4.13979992d-03
COFD(8617) = -1.99702270d+01
COFD(8618) = 5.04961455d+00
COFD(8619) = -4.34856114d-01
COFD(8620) = 1.86234058d-02
COFD(8621) = -2.09155238d+01
COFD(8622) = 5.37089858d+00
COFD(8623) = -4.70314009d-01
COFD(8624) = 1.99113180d-02
COFD(8625) = -2.18991553d+01
COFD(8626) = 5.50125258d+00
COFD(8627) = -4.62458530d-01
COFD(8628) = 1.85853825d-02
COFD(8629) = -2.18091268d+01
COFD(8630) = 5.36691865d+00
COFD(8631) = -4.37522912d-01
COFD(8632) = 1.72059807d-02
COFD(8633) = -2.13278787d+01
COFD(8634) = 5.49726595d+00
COFD(8635) = -4.82018917d-01
COFD(8636) = 2.02314792d-02
COFD(8637) = -2.01303629d+01
COFD(8638) = 5.10847039d+00
COFD(8639) = -4.41746871d-01
COFD(8640) = 1.88927540d-02
COFD(8641) = -1.71692288d+01
COFD(8642) = 3.10693405d+00
COFD(8643) = -8.69724915d-02
COFD(8644) = -2.16324201d-04
COFD(8645) = -1.71782995d+01
COFD(8646) = 3.10693405d+00
COFD(8647) = -8.69724915d-02
COFD(8648) = -2.16324201d-04
COFD(8649) = -2.15423882d+01
COFD(8650) = 5.60052862d+00
COFD(8651) = -4.87676563d-01
COFD(8652) = 2.01715067d-02
COFD(8653) = -2.18047796d+01
COFD(8654) = 5.36691865d+00
COFD(8655) = -4.37522912d-01
COFD(8656) = 1.72059807d-02
COFD(8657) = -2.18047796d+01
COFD(8658) = 5.36691865d+00
COFD(8659) = -4.37522912d-01
COFD(8660) = 1.72059807d-02
COFD(8661) = -2.18047796d+01
COFD(8662) = 5.36691865d+00
COFD(8663) = -4.37522912d-01
COFD(8664) = 1.72059807d-02
COFD(8665) = -2.18001733d+01
COFD(8666) = 5.36691865d+00
COFD(8667) = -4.37522912d-01
COFD(8668) = 1.72059807d-02
COFD(8669) = -2.08989475d+01
COFD(8670) = 5.37089858d+00
COFD(8671) = -4.70314009d-01
COFD(8672) = 1.99113180d-02
COFD(8673) = -2.17175057d+01
COFD(8674) = 5.60287101d+00
COFD(8675) = -4.91302853d-01
COFD(8676) = 2.04601038d-02
COFD(8677) = -2.17968118d+01
COFD(8678) = 5.18435436d+00
COFD(8679) = -4.06325299d-01
COFD(8680) = 1.55562327d-02
COFD(8681) = -2.18012366d+01
COFD(8682) = 5.18435436d+00
COFD(8683) = -4.06325299d-01
COFD(8684) = 1.55562327d-02
COFD(8685) = -1.92114566d+01
COFD(8686) = 3.99103591d+00
COFD(8687) = -2.19677673d-01
COFD(8688) = 6.21875426d-03
COFD(8689) = -1.92158886d+01
COFD(8690) = 3.99103591d+00
COFD(8691) = -2.19677673d-01
COFD(8692) = 6.21875426d-03
COFD(8693) = -1.64805864d+01
COFD(8694) = 4.26176920d+00
COFD(8695) = -3.39986562d-01
COFD(8696) = 1.48077740d-02
COFD(8697) = -1.99782082d+01
COFD(8698) = 5.61184673d+00
COFD(8699) = -4.90499758d-01
COFD(8700) = 2.03480832d-02
COFD(8701) = -1.98411047d+01
COFD(8702) = 5.15754807d+00
COFD(8703) = -4.46654084d-01
COFD(8704) = 1.90458977d-02
COFD(8705) = -2.10749998d+01
COFD(8706) = 5.41544337d+00
COFD(8707) = -4.73380953d-01
COFD(8708) = 1.99350119d-02
COFD(8709) = -1.98603359d+01
COFD(8710) = 5.15754807d+00
COFD(8711) = -4.46654084d-01
COFD(8712) = 1.90458977d-02
COFD(8713) = -1.67068906d+01
COFD(8714) = 3.08217426d+00
COFD(8715) = -8.33350945d-02
COFD(8716) = -3.89962983d-04
COFD(8717) = -2.10821814d+01
COFD(8718) = 5.41544337d+00
COFD(8719) = -4.73380953d-01
COFD(8720) = 1.99350119d-02
COFD(8721) = -2.10890329d+01
COFD(8722) = 5.41544337d+00
COFD(8723) = -4.73380953d-01
COFD(8724) = 1.99350119d-02
COFD(8725) = -1.95761623d+01
COFD(8726) = 5.04961455d+00
COFD(8727) = -4.34856114d-01
COFD(8728) = 1.86234058d-02
COFD(8729) = -1.97730795d+01
COFD(8730) = 5.15754807d+00
COFD(8731) = -4.46654084d-01
COFD(8732) = 1.90458977d-02
COFD(8733) = -2.16413359d+01
COFD(8734) = 5.61201619d+00
COFD(8735) = -4.90763024d-01
COFD(8736) = 2.03690162d-02
COFD(8737) = -2.16413359d+01
COFD(8738) = 5.61201619d+00
COFD(8739) = -4.90763024d-01
COFD(8740) = 2.03690162d-02
COFD(8741) = -2.16641921d+01
COFD(8742) = 5.61201619d+00
COFD(8743) = -4.90763024d-01
COFD(8744) = 2.03690162d-02
COFD(8745) = -2.16560648d+01
COFD(8746) = 5.61229743d+00
COFD(8747) = -4.91428747d-01
COFD(8748) = 2.04226556d-02
COFD(8749) = -2.09276290d+01
COFD(8750) = 5.37381486d+00
COFD(8751) = -4.70526046d-01
COFD(8752) = 1.99137630d-02
COFD(8753) = -2.16971429d+01
COFD(8754) = 5.30214979d+00
COFD(8755) = -4.26479693d-01
COFD(8756) = 1.66224918d-02
COFD(8757) = -1.82289601d+01
COFD(8758) = 3.59760627d+00
COFD(8759) = -1.59809307d-01
COFD(8760) = 3.28676045d-03
COFD(8761) = -1.82372764d+01
COFD(8762) = 3.59760627d+00
COFD(8763) = -1.59809307d-01
COFD(8764) = 3.28676045d-03
COFD(8765) = -1.93995811d+01
COFD(8766) = 4.14577053d+00
COFD(8767) = -2.43003893d-01
COFD(8768) = 7.35638742d-03
COFD(8769) = -1.93995811d+01
COFD(8770) = 4.14577053d+00
COFD(8771) = -2.43003893d-01
COFD(8772) = 7.35638742d-03
COFD(8773) = -1.84947374d+01
COFD(8774) = 3.70911349d+00
COFD(8775) = -1.76610648d-01
COFD(8776) = 4.10440582d-03
COFD(8777) = -2.19109699d+01
COFD(8778) = 5.47350698d+00
COFD(8779) = -4.56806529d-01
COFD(8780) = 1.82590273d-02
COFD(8781) = -2.19213014d+01
COFD(8782) = 5.47350698d+00
COFD(8783) = -4.56806529d-01
COFD(8784) = 1.82590273d-02
COFD(8785) = -2.19310570d+01
COFD(8786) = 5.47350698d+00
COFD(8787) = -4.56806529d-01
COFD(8788) = 1.82590273d-02
COFD(8789) = -2.13148599d+01
COFD(8790) = 5.10013279d+00
COFD(8791) = -3.92393844d-01
COFD(8792) = 1.48333781d-02
COFD(8793) = -2.16712796d+01
COFD(8794) = 5.25608421d+00
COFD(8795) = -4.18625448d-01
COFD(8796) = 1.62072845d-02
COFD(8797) = -2.16795773d+01
COFD(8798) = 5.25608421d+00
COFD(8799) = -4.18625448d-01
COFD(8800) = 1.62072845d-02
COFD(8801) = -2.15533322d+01
COFD(8802) = 5.60052862d+00
COFD(8803) = -4.87676563d-01
COFD(8804) = 2.01715067d-02
COFD(8805) = -1.92179182d+01
COFD(8806) = 3.99103591d+00
COFD(8807) = -2.19677673d-01
COFD(8808) = 6.21875426d-03
COFD(8809) = -1.92179182d+01
COFD(8810) = 3.99103591d+00
COFD(8811) = -2.19677673d-01
COFD(8812) = 6.21875426d-03
COFD(8813) = -1.96286971d+01
COFD(8814) = 5.04961455d+00
COFD(8815) = -4.34856114d-01
COFD(8816) = 1.86234058d-02
COFD(8817) = -1.97890553d+01
COFD(8818) = 5.15754807d+00
COFD(8819) = -4.46654084d-01
COFD(8820) = 1.90458977d-02
COFD(8821) = -1.98099704d+01
COFD(8822) = 5.15754807d+00
COFD(8823) = -4.46654084d-01
COFD(8824) = 1.90458977d-02
COFD(8825) = -1.81274137d+01
COFD(8826) = 3.71386126d+00
COFD(8827) = -1.77333629d-01
COFD(8828) = 4.13979992d-03
COFD(8829) = -1.99796237d+01
COFD(8830) = 5.04961455d+00
COFD(8831) = -4.34856114d-01
COFD(8832) = 1.86234058d-02
COFD(8833) = -2.09250730d+01
COFD(8834) = 5.37089858d+00
COFD(8835) = -4.70314009d-01
COFD(8836) = 1.99113180d-02
COFD(8837) = -2.19106105d+01
COFD(8838) = 5.50125258d+00
COFD(8839) = -4.62458530d-01
COFD(8840) = 1.85853825d-02
COFD(8841) = -2.18203911d+01
COFD(8842) = 5.36691865d+00
COFD(8843) = -4.37522912d-01
COFD(8844) = 1.72059807d-02
COFD(8845) = -2.13375787d+01
COFD(8846) = 5.49726595d+00
COFD(8847) = -4.82018917d-01
COFD(8848) = 2.02314792d-02
COFD(8849) = -2.01392596d+01
COFD(8850) = 5.10847039d+00
COFD(8851) = -4.41746871d-01
COFD(8852) = 1.88927540d-02
COFD(8853) = -1.71782995d+01
COFD(8854) = 3.10693405d+00
COFD(8855) = -8.69724915d-02
COFD(8856) = -2.16324201d-04
COFD(8857) = -1.71875378d+01
COFD(8858) = 3.10693405d+00
COFD(8859) = -8.69724915d-02
COFD(8860) = -2.16324201d-04
COFD(8861) = -2.15533463d+01
COFD(8862) = 5.60052862d+00
COFD(8863) = -4.87676563d-01
COFD(8864) = 2.01715067d-02
COFD(8865) = -2.18159452d+01
COFD(8866) = 5.36691865d+00
COFD(8867) = -4.37522912d-01
COFD(8868) = 1.72059807d-02
COFD(8869) = -2.18159452d+01
COFD(8870) = 5.36691865d+00
COFD(8871) = -4.37522912d-01
COFD(8872) = 1.72059807d-02
COFD(8873) = -2.18159452d+01
COFD(8874) = 5.36691865d+00
COFD(8875) = -4.37522912d-01
COFD(8876) = 1.72059807d-02
COFD(8877) = -2.18112354d+01
COFD(8878) = 5.36691865d+00
COFD(8879) = -4.37522912d-01
COFD(8880) = 1.72059807d-02
COFD(8881) = -2.09081824d+01
COFD(8882) = 5.37089858d+00
COFD(8883) = -4.70314009d-01
COFD(8884) = 1.99113180d-02
COFD(8885) = -2.17283455d+01
COFD(8886) = 5.60287101d+00
COFD(8887) = -4.91302853d-01
COFD(8888) = 2.04601038d-02
COFD(8889) = -2.18079839d+01
COFD(8890) = 5.18435436d+00
COFD(8891) = -4.06325299d-01
COFD(8892) = 1.55562327d-02
COFD(8893) = -2.18125092d+01
COFD(8894) = 5.18435436d+00
COFD(8895) = -4.06325299d-01
COFD(8896) = 1.55562327d-02
COFD(8897) = -1.92226244d+01
COFD(8898) = 3.99103591d+00
COFD(8899) = -2.19677673d-01
COFD(8900) = 6.21875426d-03
COFD(8901) = -1.92271569d+01
COFD(8902) = 3.99103591d+00
COFD(8903) = -2.19677673d-01
COFD(8904) = 6.21875426d-03
COFD(8905) = -1.22004340d+01
COFD(8906) = 2.80725489d+00
COFD(8907) = -1.54291406d-01
COFD(8908) = 6.88290911d-03
COFD(8909) = -1.54460828d+01
COFD(8910) = 4.26819983d+00
COFD(8911) = -3.40766379d-01
COFD(8912) = 1.48393361d-02
COFD(8913) = -1.49500454d+01
COFD(8914) = 3.52327209d+00
COFD(8915) = -2.46286208d-01
COFD(8916) = 1.08285963d-02
COFD(8917) = -1.64169585d+01
COFD(8918) = 3.89309916d+00
COFD(8919) = -2.93528188d-01
COFD(8920) = 1.28463177d-02
COFD(8921) = -1.49718335d+01
COFD(8922) = 3.52327209d+00
COFD(8923) = -2.46286208d-01
COFD(8924) = 1.08285963d-02
COFD(8925) = -2.10440781d+01
COFD(8926) = 5.59806282d+00
COFD(8927) = -4.87109535d-01
COFD(8928) = 2.01370226d-02
COFD(8929) = -1.64256119d+01
COFD(8930) = 3.89309916d+00
COFD(8931) = -2.93528188d-01
COFD(8932) = 1.28463177d-02
COFD(8933) = -1.64338914d+01
COFD(8934) = 3.89309916d+00
COFD(8935) = -2.93528188d-01
COFD(8936) = 1.28463177d-02
COFD(8937) = -1.46907107d+01
COFD(8938) = 3.39229020d+00
COFD(8939) = -2.29520232d-01
COFD(8940) = 1.01114311d-02
COFD(8941) = -1.48738150d+01
COFD(8942) = 3.52327209d+00
COFD(8943) = -2.46286208d-01
COFD(8944) = 1.08285963d-02
COFD(8945) = -1.72572130d+01
COFD(8946) = 4.26063341d+00
COFD(8947) = -3.39848064d-01
COFD(8948) = 1.48021313d-02
COFD(8949) = -1.72572130d+01
COFD(8950) = 4.26063341d+00
COFD(8951) = -3.39848064d-01
COFD(8952) = 1.48021313d-02
COFD(8953) = -1.72828396d+01
COFD(8954) = 4.26063341d+00
COFD(8955) = -3.39848064d-01
COFD(8956) = 1.48021313d-02
COFD(8957) = -1.72316246d+01
COFD(8958) = 4.24011069d+00
COFD(8959) = -3.37339810d-01
COFD(8960) = 1.46996679d-02
COFD(8961) = -1.60261816d+01
COFD(8962) = 3.73312045d+00
COFD(8963) = -2.72579779d-01
COFD(8964) = 1.19290272d-02
COFD(8965) = -1.94486162d+01
COFD(8966) = 4.91446566d+00
COFD(8967) = -4.18837152d-01
COFD(8968) = 1.79893537d-02
COFD(8969) = -2.14204329d+01
COFD(8970) = 5.59268435d+00
COFD(8971) = -4.91232974d-01
COFD(8972) = 2.05064746d-02
COFD(8973) = -2.14303625d+01
COFD(8974) = 5.59268435d+00
COFD(8975) = -4.91232974d-01
COFD(8976) = 2.05064746d-02
COFD(8977) = -2.09241797d+01
COFD(8978) = 5.42316225d+00
COFD(8979) = -4.73702801d-01
COFD(8980) = 1.99217718d-02
COFD(8981) = -2.09241797d+01
COFD(8982) = 5.42316225d+00
COFD(8983) = -4.73702801d-01
COFD(8984) = 1.99217718d-02
COFD(8985) = -2.13796455d+01
COFD(8986) = 5.56978987d+00
COFD(8987) = -4.89141980d-01
COFD(8988) = 2.04499210d-02
COFD(8989) = -1.88524516d+01
COFD(8990) = 4.72476764d+00
COFD(8991) = -3.96306836d-01
COFD(8992) = 1.70964541d-02
COFD(8993) = -1.88646205d+01
COFD(8994) = 4.72476764d+00
COFD(8995) = -3.96306836d-01
COFD(8996) = 1.70964541d-02
COFD(8997) = -1.88761525d+01
COFD(8998) = 4.72476764d+00
COFD(8999) = -3.96306836d-01
COFD(9000) = 1.70964541d-02
COFD(9001) = -1.99081697d+01
COFD(9002) = 5.09311649d+00
COFD(9003) = -4.39965178d-01
COFD(9004) = 1.88238537d-02
COFD(9005) = -1.96309186d+01
COFD(9006) = 4.95923807d+00
COFD(9007) = -4.24176182d-01
COFD(9008) = 1.82020215d-02
COFD(9009) = -1.96408274d+01
COFD(9010) = 4.95923807d+00
COFD(9011) = -4.24176182d-01
COFD(9012) = 1.82020215d-02
COFD(9013) = -1.72415036d+01
COFD(9014) = 4.29808578d+00
COFD(9015) = -3.44235570d-01
COFD(9016) = 1.49727727d-02
COFD(9017) = -2.12622090d+01
COFD(9018) = 5.47935225d+00
COFD(9019) = -4.80056796d-01
COFD(9020) = 2.01607180d-02
COFD(9021) = -2.12622090d+01
COFD(9022) = 5.47935225d+00
COFD(9023) = -4.80056796d-01
COFD(9024) = 2.01607180d-02
COFD(9025) = -1.47490956d+01
COFD(9026) = 3.39229020d+00
COFD(9027) = -2.29520232d-01
COFD(9028) = 1.01114311d-02
COFD(9029) = -1.48885295d+01
COFD(9030) = 3.52327209d+00
COFD(9031) = -2.46286208d-01
COFD(9032) = 1.08285963d-02
COFD(9033) = -1.49121048d+01
COFD(9034) = 3.52327209d+00
COFD(9035) = -2.46286208d-01
COFD(9036) = 1.08285963d-02
COFD(9037) = -2.09300785d+01
COFD(9038) = 5.56884637d+00
COFD(9039) = -4.89068748d-01
COFD(9040) = 2.04486922d-02
COFD(9041) = -1.51581727d+01
COFD(9042) = 3.39229020d+00
COFD(9043) = -2.29520232d-01
COFD(9044) = 1.01114311d-02
COFD(9045) = -1.60085775d+01
COFD(9046) = 3.72220402d+00
COFD(9047) = -2.71150591d-01
COFD(9048) = 1.18665265d-02
COFD(9049) = -1.86961260d+01
COFD(9050) = 4.67999453d+00
COFD(9051) = -3.91135253d-01
COFD(9052) = 1.68982388d-02
COFD(9053) = -1.92693952d+01
COFD(9054) = 4.84335373d+00
COFD(9055) = -4.10212751d-01
COFD(9056) = 1.76395874d-02
COFD(9057) = -1.67419884d+01
COFD(9058) = 4.01496916d+00
COFD(9059) = -3.09274229d-01
COFD(9060) = 1.35278248d-02
COFD(9061) = -1.53109485d+01
COFD(9062) = 3.45550279d+00
COFD(9063) = -2.37693211d-01
COFD(9064) = 1.04644205d-02
COFD(9065) = -2.15423882d+01
COFD(9066) = 5.60052862d+00
COFD(9067) = -4.87676563d-01
COFD(9068) = 2.01715067d-02
COFD(9069) = -2.15533463d+01
COFD(9070) = 5.60052862d+00
COFD(9071) = -4.87676563d-01
COFD(9072) = 2.01715067d-02
COFD(9073) = -1.72415209d+01
COFD(9074) = 4.29808578d+00
COFD(9075) = -3.44235570d-01
COFD(9076) = 1.49727727d-02
COFD(9077) = -1.92638883d+01
COFD(9078) = 4.84335373d+00
COFD(9079) = -4.10212751d-01
COFD(9080) = 1.76395874d-02
COFD(9081) = -1.92638883d+01
COFD(9082) = 4.84335373d+00
COFD(9083) = -4.10212751d-01
COFD(9084) = 1.76395874d-02
COFD(9085) = -1.92638883d+01
COFD(9086) = 4.84335373d+00
COFD(9087) = -4.10212751d-01
COFD(9088) = 1.76395874d-02
COFD(9089) = -1.92580672d+01
COFD(9090) = 4.84335373d+00
COFD(9091) = -4.10212751d-01
COFD(9092) = 1.76395874d-02
COFD(9093) = -1.59884445d+01
COFD(9094) = 3.72220402d+00
COFD(9095) = -2.71150591d-01
COFD(9096) = 1.18665265d-02
COFD(9097) = -1.72570778d+01
COFD(9098) = 4.19757624d+00
COFD(9099) = -3.32087529d-01
COFD(9100) = 1.44827462d-02
COFD(9101) = -2.01341617d+01
COFD(9102) = 5.03101171d+00
COFD(9103) = -4.32665019d-01
COFD(9104) = 1.85372086d-02
COFD(9105) = -2.01397678d+01
COFD(9106) = 5.03101171d+00
COFD(9107) = -4.32665019d-01
COFD(9108) = 1.85372086d-02
COFD(9109) = -2.12680259d+01
COFD(9110) = 5.47935225d+00
COFD(9111) = -4.80056796d-01
COFD(9112) = 2.01607180d-02
COFD(9113) = -2.12736405d+01
COFD(9114) = 5.47935225d+00
COFD(9115) = -4.80056796d-01
COFD(9116) = 2.01607180d-02
COFD(9117) = -1.36889305d+01
COFD(9118) = 3.19999740d+00
COFD(9119) = -2.04999020d-01
COFD(9120) = 9.06766105d-03
COFD(9121) = -1.74816531d+01
COFD(9122) = 4.80792005d+00
COFD(9123) = -4.06126584d-01
COFD(9124) = 1.74831083d-02
COFD(9125) = -1.69259591d+01
COFD(9126) = 4.09077642d+00
COFD(9127) = -3.18894990d-01
COFD(9128) = 1.39371445d-02
COFD(9129) = -1.83260311d+01
COFD(9130) = 4.43878381d+00
COFD(9131) = -3.61697624d-01
COFD(9132) = 1.56975581d-02
COFD(9133) = -1.69480404d+01
COFD(9134) = 4.09077642d+00
COFD(9135) = -3.18894990d-01
COFD(9136) = 1.39371445d-02
COFD(9137) = -2.13622815d+01
COFD(9138) = 5.35901292d+00
COFD(9139) = -4.36172487d-01
COFD(9140) = 1.71345319d-02
COFD(9141) = -1.83348653d+01
COFD(9142) = 4.43878381d+00
COFD(9143) = -3.61697624d-01
COFD(9144) = 1.56975581d-02
COFD(9145) = -1.83433209d+01
COFD(9146) = 4.43878381d+00
COFD(9147) = -3.61697624d-01
COFD(9148) = 1.56975581d-02
COFD(9149) = -1.65697233d+01
COFD(9150) = 3.93849401d+00
COFD(9151) = -2.99416642d-01
COFD(9152) = 1.31020815d-02
COFD(9153) = -1.68487987d+01
COFD(9154) = 4.09077642d+00
COFD(9155) = -3.18894990d-01
COFD(9156) = 1.39371445d-02
COFD(9157) = -1.91676574d+01
COFD(9158) = 4.80030220d+00
COFD(9159) = -4.05235041d-01
COFD(9160) = 1.74483531d-02
COFD(9161) = -1.91676574d+01
COFD(9162) = 4.80030220d+00
COFD(9163) = -4.05235041d-01
COFD(9164) = 1.74483531d-02
COFD(9165) = -1.91935980d+01
COFD(9166) = 4.80030220d+00
COFD(9167) = -4.05235041d-01
COFD(9168) = 1.74483531d-02
COFD(9169) = -1.91463450d+01
COFD(9170) = 4.78094221d+00
COFD(9171) = -4.02985837d-01
COFD(9172) = 1.73614221d-02
COFD(9173) = -1.80021546d+01
COFD(9174) = 4.31656593d+00
COFD(9175) = -3.46539554d-01
COFD(9176) = 1.50688196d-02
COFD(9177) = -2.12614542d+01
COFD(9178) = 5.38135645d+00
COFD(9179) = -4.71058360d-01
COFD(9180) = 1.99188046d-02
COFD(9181) = -2.19914997d+01
COFD(9182) = 5.49642957d+00
COFD(9183) = -4.61132993d-01
COFD(9184) = 1.85004773d-02
COFD(9185) = -2.20016256d+01
COFD(9186) = 5.49642957d+00
COFD(9187) = -4.61132993d-01
COFD(9188) = 1.85004773d-02
COFD(9189) = -2.20823296d+01
COFD(9190) = 5.58518321d+00
COFD(9191) = -4.80534235d-01
COFD(9192) = 1.96556393d-02
COFD(9193) = -2.20823296d+01
COFD(9194) = 5.58518321d+00
COFD(9195) = -4.80534235d-01
COFD(9196) = 1.96556393d-02
COFD(9197) = -2.20212574d+01
COFD(9198) = 5.50583166d+00
COFD(9199) = -4.63753262d-01
COFD(9200) = 1.86693462d-02
COFD(9201) = -2.07170423d+01
COFD(9202) = 5.23116678d+00
COFD(9203) = -4.54972675d-01
COFD(9204) = 1.93572354d-02
COFD(9205) = -2.07294312d+01
COFD(9206) = 5.23116678d+00
COFD(9207) = -4.54972675d-01
COFD(9208) = 1.93572354d-02
COFD(9209) = -2.07411769d+01
COFD(9210) = 5.23116678d+00
COFD(9211) = -4.54972675d-01
COFD(9212) = 1.93572354d-02
COFD(9213) = -2.15270531d+01
COFD(9214) = 5.47662534d+00
COFD(9215) = -4.79750189d-01
COFD(9216) = 2.01492627d-02
COFD(9217) = -2.13502469d+01
COFD(9218) = 5.39784038d+00
COFD(9219) = -4.72269379d-01
COFD(9220) = 1.99339592d-02
COFD(9221) = -2.13603517d+01
COFD(9222) = 5.39784038d+00
COFD(9223) = -4.72269379d-01
COFD(9224) = 1.99339592d-02
COFD(9225) = -1.92638705d+01
COFD(9226) = 4.84335373d+00
COFD(9227) = -4.10212751d-01
COFD(9228) = 1.76395874d-02
COFD(9229) = -2.22648059d+01
COFD(9230) = 5.57944619d+00
COFD(9231) = -4.78032626d-01
COFD(9232) = 1.94775370d-02
COFD(9233) = -2.22648059d+01
COFD(9234) = 5.57944619d+00
COFD(9235) = -4.78032626d-01
COFD(9236) = 1.94775370d-02
COFD(9237) = -1.66287650d+01
COFD(9238) = 3.93849401d+00
COFD(9239) = -2.99416642d-01
COFD(9240) = 1.31020815d-02
COFD(9241) = -1.68719715d+01
COFD(9242) = 4.09077642d+00
COFD(9243) = -3.18894990d-01
COFD(9244) = 1.39371445d-02
COFD(9245) = -1.68958501d+01
COFD(9246) = 4.09077642d+00
COFD(9247) = -3.18894990d-01
COFD(9248) = 1.39371445d-02
COFD(9249) = -2.16179118d+01
COFD(9250) = 5.50657794d+00
COFD(9251) = -4.63927523d-01
COFD(9252) = 1.86800112d-02
COFD(9253) = -1.70117892d+01
COFD(9254) = 3.93849401d+00
COFD(9255) = -2.99416642d-01
COFD(9256) = 1.31020815d-02
COFD(9257) = -1.79926243d+01
COFD(9258) = 4.30841971d+00
COFD(9259) = -3.45524579d-01
COFD(9260) = 1.50265381d-02
COFD(9261) = -2.05658452d+01
COFD(9262) = 5.17526774d+00
COFD(9263) = -4.48472252d-01
COFD(9264) = 1.91050891d-02
COFD(9265) = -2.11852687d+01
COFD(9266) = 5.35817304d+00
COFD(9267) = -4.69455306d-01
COFD(9268) = 1.99063292d-02
COFD(9269) = -1.86041410d+01
COFD(9270) = 4.53681554d+00
COFD(9271) = -3.73517224d-01
COFD(9272) = 1.61730724d-02
COFD(9273) = -1.71806065d+01
COFD(9274) = 4.00896225d+00
COFD(9275) = -3.08501394d-01
COFD(9276) = 1.34944985d-02
COFD(9277) = -2.18047796d+01
COFD(9278) = 5.36691865d+00
COFD(9279) = -4.37522912d-01
COFD(9280) = 1.72059807d-02
COFD(9281) = -2.18159452d+01
COFD(9282) = 5.36691865d+00
COFD(9283) = -4.37522912d-01
COFD(9284) = 1.72059807d-02
COFD(9285) = -1.92638883d+01
COFD(9286) = 4.84335373d+00
COFD(9287) = -4.10212751d-01
COFD(9288) = 1.76395874d-02
COFD(9289) = -2.11796273d+01
COFD(9290) = 5.35817304d+00
COFD(9291) = -4.69455306d-01
COFD(9292) = 1.99063292d-02
COFD(9293) = -2.11796273d+01
COFD(9294) = 5.35817304d+00
COFD(9295) = -4.69455306d-01
COFD(9296) = 1.99063292d-02
COFD(9297) = -2.11796273d+01
COFD(9298) = 5.35817304d+00
COFD(9299) = -4.69455306d-01
COFD(9300) = 1.99063292d-02
COFD(9301) = -2.11736657d+01
COFD(9302) = 5.35817304d+00
COFD(9303) = -4.69455306d-01
COFD(9304) = 1.99063292d-02
COFD(9305) = -1.79720978d+01
COFD(9306) = 4.30841971d+00
COFD(9307) = -3.45524579d-01
COFD(9308) = 1.50265381d-02
COFD(9309) = -1.91975423d+01
COFD(9310) = 4.73837374d+00
COFD(9311) = -3.97916523d-01
COFD(9312) = 1.71599378d-02
COFD(9313) = -2.17133615d+01
COFD(9314) = 5.41948945d+00
COFD(9315) = -4.73474138d-01
COFD(9316) = 1.99217500d-02
COFD(9317) = -2.17191046d+01
COFD(9318) = 5.41948945d+00
COFD(9319) = -4.73474138d-01
COFD(9320) = 1.99217500d-02
COFD(9321) = -2.22707633d+01
COFD(9322) = 5.57944619d+00
COFD(9323) = -4.78032626d-01
COFD(9324) = 1.94775370d-02
COFD(9325) = -2.22765150d+01
COFD(9326) = 5.57944619d+00
COFD(9327) = -4.78032626d-01
COFD(9328) = 1.94775370d-02
COFD(9329) = -1.36889305d+01
COFD(9330) = 3.19999740d+00
COFD(9331) = -2.04999020d-01
COFD(9332) = 9.06766105d-03
COFD(9333) = -1.74816531d+01
COFD(9334) = 4.80792005d+00
COFD(9335) = -4.06126584d-01
COFD(9336) = 1.74831083d-02
COFD(9337) = -1.69259591d+01
COFD(9338) = 4.09077642d+00
COFD(9339) = -3.18894990d-01
COFD(9340) = 1.39371445d-02
COFD(9341) = -1.83260311d+01
COFD(9342) = 4.43878381d+00
COFD(9343) = -3.61697624d-01
COFD(9344) = 1.56975581d-02
COFD(9345) = -1.69480404d+01
COFD(9346) = 4.09077642d+00
COFD(9347) = -3.18894990d-01
COFD(9348) = 1.39371445d-02
COFD(9349) = -2.13622815d+01
COFD(9350) = 5.35901292d+00
COFD(9351) = -4.36172487d-01
COFD(9352) = 1.71345319d-02
COFD(9353) = -1.83348653d+01
COFD(9354) = 4.43878381d+00
COFD(9355) = -3.61697624d-01
COFD(9356) = 1.56975581d-02
COFD(9357) = -1.83433209d+01
COFD(9358) = 4.43878381d+00
COFD(9359) = -3.61697624d-01
COFD(9360) = 1.56975581d-02
COFD(9361) = -1.65697233d+01
COFD(9362) = 3.93849401d+00
COFD(9363) = -2.99416642d-01
COFD(9364) = 1.31020815d-02
COFD(9365) = -1.68487987d+01
COFD(9366) = 4.09077642d+00
COFD(9367) = -3.18894990d-01
COFD(9368) = 1.39371445d-02
COFD(9369) = -1.91676574d+01
COFD(9370) = 4.80030220d+00
COFD(9371) = -4.05235041d-01
COFD(9372) = 1.74483531d-02
COFD(9373) = -1.91676574d+01
COFD(9374) = 4.80030220d+00
COFD(9375) = -4.05235041d-01
COFD(9376) = 1.74483531d-02
COFD(9377) = -1.91935980d+01
COFD(9378) = 4.80030220d+00
COFD(9379) = -4.05235041d-01
COFD(9380) = 1.74483531d-02
COFD(9381) = -1.91463450d+01
COFD(9382) = 4.78094221d+00
COFD(9383) = -4.02985837d-01
COFD(9384) = 1.73614221d-02
COFD(9385) = -1.80021546d+01
COFD(9386) = 4.31656593d+00
COFD(9387) = -3.46539554d-01
COFD(9388) = 1.50688196d-02
COFD(9389) = -2.12614542d+01
COFD(9390) = 5.38135645d+00
COFD(9391) = -4.71058360d-01
COFD(9392) = 1.99188046d-02
COFD(9393) = -2.19914997d+01
COFD(9394) = 5.49642957d+00
COFD(9395) = -4.61132993d-01
COFD(9396) = 1.85004773d-02
COFD(9397) = -2.20016256d+01
COFD(9398) = 5.49642957d+00
COFD(9399) = -4.61132993d-01
COFD(9400) = 1.85004773d-02
COFD(9401) = -2.20823296d+01
COFD(9402) = 5.58518321d+00
COFD(9403) = -4.80534235d-01
COFD(9404) = 1.96556393d-02
COFD(9405) = -2.20823296d+01
COFD(9406) = 5.58518321d+00
COFD(9407) = -4.80534235d-01
COFD(9408) = 1.96556393d-02
COFD(9409) = -2.20212574d+01
COFD(9410) = 5.50583166d+00
COFD(9411) = -4.63753262d-01
COFD(9412) = 1.86693462d-02
COFD(9413) = -2.07170423d+01
COFD(9414) = 5.23116678d+00
COFD(9415) = -4.54972675d-01
COFD(9416) = 1.93572354d-02
COFD(9417) = -2.07294312d+01
COFD(9418) = 5.23116678d+00
COFD(9419) = -4.54972675d-01
COFD(9420) = 1.93572354d-02
COFD(9421) = -2.07411769d+01
COFD(9422) = 5.23116678d+00
COFD(9423) = -4.54972675d-01
COFD(9424) = 1.93572354d-02
COFD(9425) = -2.15270531d+01
COFD(9426) = 5.47662534d+00
COFD(9427) = -4.79750189d-01
COFD(9428) = 2.01492627d-02
COFD(9429) = -2.13502469d+01
COFD(9430) = 5.39784038d+00
COFD(9431) = -4.72269379d-01
COFD(9432) = 1.99339592d-02
COFD(9433) = -2.13603517d+01
COFD(9434) = 5.39784038d+00
COFD(9435) = -4.72269379d-01
COFD(9436) = 1.99339592d-02
COFD(9437) = -1.92638705d+01
COFD(9438) = 4.84335373d+00
COFD(9439) = -4.10212751d-01
COFD(9440) = 1.76395874d-02
COFD(9441) = -2.22648059d+01
COFD(9442) = 5.57944619d+00
COFD(9443) = -4.78032626d-01
COFD(9444) = 1.94775370d-02
COFD(9445) = -2.22648059d+01
COFD(9446) = 5.57944619d+00
COFD(9447) = -4.78032626d-01
COFD(9448) = 1.94775370d-02
COFD(9449) = -1.66287650d+01
COFD(9450) = 3.93849401d+00
COFD(9451) = -2.99416642d-01
COFD(9452) = 1.31020815d-02
COFD(9453) = -1.68719715d+01
COFD(9454) = 4.09077642d+00
COFD(9455) = -3.18894990d-01
COFD(9456) = 1.39371445d-02
COFD(9457) = -1.68958501d+01
COFD(9458) = 4.09077642d+00
COFD(9459) = -3.18894990d-01
COFD(9460) = 1.39371445d-02
COFD(9461) = -2.16179118d+01
COFD(9462) = 5.50657794d+00
COFD(9463) = -4.63927523d-01
COFD(9464) = 1.86800112d-02
COFD(9465) = -1.70117892d+01
COFD(9466) = 3.93849401d+00
COFD(9467) = -2.99416642d-01
COFD(9468) = 1.31020815d-02
COFD(9469) = -1.79926243d+01
COFD(9470) = 4.30841971d+00
COFD(9471) = -3.45524579d-01
COFD(9472) = 1.50265381d-02
COFD(9473) = -2.05658452d+01
COFD(9474) = 5.17526774d+00
COFD(9475) = -4.48472252d-01
COFD(9476) = 1.91050891d-02
COFD(9477) = -2.11852687d+01
COFD(9478) = 5.35817304d+00
COFD(9479) = -4.69455306d-01
COFD(9480) = 1.99063292d-02
COFD(9481) = -1.86041410d+01
COFD(9482) = 4.53681554d+00
COFD(9483) = -3.73517224d-01
COFD(9484) = 1.61730724d-02
COFD(9485) = -1.71806065d+01
COFD(9486) = 4.00896225d+00
COFD(9487) = -3.08501394d-01
COFD(9488) = 1.34944985d-02
COFD(9489) = -2.18047796d+01
COFD(9490) = 5.36691865d+00
COFD(9491) = -4.37522912d-01
COFD(9492) = 1.72059807d-02
COFD(9493) = -2.18159452d+01
COFD(9494) = 5.36691865d+00
COFD(9495) = -4.37522912d-01
COFD(9496) = 1.72059807d-02
COFD(9497) = -1.92638883d+01
COFD(9498) = 4.84335373d+00
COFD(9499) = -4.10212751d-01
COFD(9500) = 1.76395874d-02
COFD(9501) = -2.11796273d+01
COFD(9502) = 5.35817304d+00
COFD(9503) = -4.69455306d-01
COFD(9504) = 1.99063292d-02
COFD(9505) = -2.11796273d+01
COFD(9506) = 5.35817304d+00
COFD(9507) = -4.69455306d-01
COFD(9508) = 1.99063292d-02
COFD(9509) = -2.11796273d+01
COFD(9510) = 5.35817304d+00
COFD(9511) = -4.69455306d-01
COFD(9512) = 1.99063292d-02
COFD(9513) = -2.11736657d+01
COFD(9514) = 5.35817304d+00
COFD(9515) = -4.69455306d-01
COFD(9516) = 1.99063292d-02
COFD(9517) = -1.79720978d+01
COFD(9518) = 4.30841971d+00
COFD(9519) = -3.45524579d-01
COFD(9520) = 1.50265381d-02
COFD(9521) = -1.91975423d+01
COFD(9522) = 4.73837374d+00
COFD(9523) = -3.97916523d-01
COFD(9524) = 1.71599378d-02
COFD(9525) = -2.17133615d+01
COFD(9526) = 5.41948945d+00
COFD(9527) = -4.73474138d-01
COFD(9528) = 1.99217500d-02
COFD(9529) = -2.17191046d+01
COFD(9530) = 5.41948945d+00
COFD(9531) = -4.73474138d-01
COFD(9532) = 1.99217500d-02
COFD(9533) = -2.22707633d+01
COFD(9534) = 5.57944619d+00
COFD(9535) = -4.78032626d-01
COFD(9536) = 1.94775370d-02
COFD(9537) = -2.22765150d+01
COFD(9538) = 5.57944619d+00
COFD(9539) = -4.78032626d-01
COFD(9540) = 1.94775370d-02
COFD(9541) = -1.36889305d+01
COFD(9542) = 3.19999740d+00
COFD(9543) = -2.04999020d-01
COFD(9544) = 9.06766105d-03
COFD(9545) = -1.74816531d+01
COFD(9546) = 4.80792005d+00
COFD(9547) = -4.06126584d-01
COFD(9548) = 1.74831083d-02
COFD(9549) = -1.69259591d+01
COFD(9550) = 4.09077642d+00
COFD(9551) = -3.18894990d-01
COFD(9552) = 1.39371445d-02
COFD(9553) = -1.83260311d+01
COFD(9554) = 4.43878381d+00
COFD(9555) = -3.61697624d-01
COFD(9556) = 1.56975581d-02
COFD(9557) = -1.69480404d+01
COFD(9558) = 4.09077642d+00
COFD(9559) = -3.18894990d-01
COFD(9560) = 1.39371445d-02
COFD(9561) = -2.13622815d+01
COFD(9562) = 5.35901292d+00
COFD(9563) = -4.36172487d-01
COFD(9564) = 1.71345319d-02
COFD(9565) = -1.83348653d+01
COFD(9566) = 4.43878381d+00
COFD(9567) = -3.61697624d-01
COFD(9568) = 1.56975581d-02
COFD(9569) = -1.83433209d+01
COFD(9570) = 4.43878381d+00
COFD(9571) = -3.61697624d-01
COFD(9572) = 1.56975581d-02
COFD(9573) = -1.65697233d+01
COFD(9574) = 3.93849401d+00
COFD(9575) = -2.99416642d-01
COFD(9576) = 1.31020815d-02
COFD(9577) = -1.68487987d+01
COFD(9578) = 4.09077642d+00
COFD(9579) = -3.18894990d-01
COFD(9580) = 1.39371445d-02
COFD(9581) = -1.91676574d+01
COFD(9582) = 4.80030220d+00
COFD(9583) = -4.05235041d-01
COFD(9584) = 1.74483531d-02
COFD(9585) = -1.91676574d+01
COFD(9586) = 4.80030220d+00
COFD(9587) = -4.05235041d-01
COFD(9588) = 1.74483531d-02
COFD(9589) = -1.91935980d+01
COFD(9590) = 4.80030220d+00
COFD(9591) = -4.05235041d-01
COFD(9592) = 1.74483531d-02
COFD(9593) = -1.91463450d+01
COFD(9594) = 4.78094221d+00
COFD(9595) = -4.02985837d-01
COFD(9596) = 1.73614221d-02
COFD(9597) = -1.80021546d+01
COFD(9598) = 4.31656593d+00
COFD(9599) = -3.46539554d-01
COFD(9600) = 1.50688196d-02
COFD(9601) = -2.12614542d+01
COFD(9602) = 5.38135645d+00
COFD(9603) = -4.71058360d-01
COFD(9604) = 1.99188046d-02
COFD(9605) = -2.19914997d+01
COFD(9606) = 5.49642957d+00
COFD(9607) = -4.61132993d-01
COFD(9608) = 1.85004773d-02
COFD(9609) = -2.20016256d+01
COFD(9610) = 5.49642957d+00
COFD(9611) = -4.61132993d-01
COFD(9612) = 1.85004773d-02
COFD(9613) = -2.20823296d+01
COFD(9614) = 5.58518321d+00
COFD(9615) = -4.80534235d-01
COFD(9616) = 1.96556393d-02
COFD(9617) = -2.20823296d+01
COFD(9618) = 5.58518321d+00
COFD(9619) = -4.80534235d-01
COFD(9620) = 1.96556393d-02
COFD(9621) = -2.20212574d+01
COFD(9622) = 5.50583166d+00
COFD(9623) = -4.63753262d-01
COFD(9624) = 1.86693462d-02
COFD(9625) = -2.07170423d+01
COFD(9626) = 5.23116678d+00
COFD(9627) = -4.54972675d-01
COFD(9628) = 1.93572354d-02
COFD(9629) = -2.07294312d+01
COFD(9630) = 5.23116678d+00
COFD(9631) = -4.54972675d-01
COFD(9632) = 1.93572354d-02
COFD(9633) = -2.07411769d+01
COFD(9634) = 5.23116678d+00
COFD(9635) = -4.54972675d-01
COFD(9636) = 1.93572354d-02
COFD(9637) = -2.15270531d+01
COFD(9638) = 5.47662534d+00
COFD(9639) = -4.79750189d-01
COFD(9640) = 2.01492627d-02
COFD(9641) = -2.13502469d+01
COFD(9642) = 5.39784038d+00
COFD(9643) = -4.72269379d-01
COFD(9644) = 1.99339592d-02
COFD(9645) = -2.13603517d+01
COFD(9646) = 5.39784038d+00
COFD(9647) = -4.72269379d-01
COFD(9648) = 1.99339592d-02
COFD(9649) = -1.92638705d+01
COFD(9650) = 4.84335373d+00
COFD(9651) = -4.10212751d-01
COFD(9652) = 1.76395874d-02
COFD(9653) = -2.22648059d+01
COFD(9654) = 5.57944619d+00
COFD(9655) = -4.78032626d-01
COFD(9656) = 1.94775370d-02
COFD(9657) = -2.22648059d+01
COFD(9658) = 5.57944619d+00
COFD(9659) = -4.78032626d-01
COFD(9660) = 1.94775370d-02
COFD(9661) = -1.66287650d+01
COFD(9662) = 3.93849401d+00
COFD(9663) = -2.99416642d-01
COFD(9664) = 1.31020815d-02
COFD(9665) = -1.68719715d+01
COFD(9666) = 4.09077642d+00
COFD(9667) = -3.18894990d-01
COFD(9668) = 1.39371445d-02
COFD(9669) = -1.68958501d+01
COFD(9670) = 4.09077642d+00
COFD(9671) = -3.18894990d-01
COFD(9672) = 1.39371445d-02
COFD(9673) = -2.16179118d+01
COFD(9674) = 5.50657794d+00
COFD(9675) = -4.63927523d-01
COFD(9676) = 1.86800112d-02
COFD(9677) = -1.70117892d+01
COFD(9678) = 3.93849401d+00
COFD(9679) = -2.99416642d-01
COFD(9680) = 1.31020815d-02
COFD(9681) = -1.79926243d+01
COFD(9682) = 4.30841971d+00
COFD(9683) = -3.45524579d-01
COFD(9684) = 1.50265381d-02
COFD(9685) = -2.05658452d+01
COFD(9686) = 5.17526774d+00
COFD(9687) = -4.48472252d-01
COFD(9688) = 1.91050891d-02
COFD(9689) = -2.11852687d+01
COFD(9690) = 5.35817304d+00
COFD(9691) = -4.69455306d-01
COFD(9692) = 1.99063292d-02
COFD(9693) = -1.86041410d+01
COFD(9694) = 4.53681554d+00
COFD(9695) = -3.73517224d-01
COFD(9696) = 1.61730724d-02
COFD(9697) = -1.71806065d+01
COFD(9698) = 4.00896225d+00
COFD(9699) = -3.08501394d-01
COFD(9700) = 1.34944985d-02
COFD(9701) = -2.18047796d+01
COFD(9702) = 5.36691865d+00
COFD(9703) = -4.37522912d-01
COFD(9704) = 1.72059807d-02
COFD(9705) = -2.18159452d+01
COFD(9706) = 5.36691865d+00
COFD(9707) = -4.37522912d-01
COFD(9708) = 1.72059807d-02
COFD(9709) = -1.92638883d+01
COFD(9710) = 4.84335373d+00
COFD(9711) = -4.10212751d-01
COFD(9712) = 1.76395874d-02
COFD(9713) = -2.11796273d+01
COFD(9714) = 5.35817304d+00
COFD(9715) = -4.69455306d-01
COFD(9716) = 1.99063292d-02
COFD(9717) = -2.11796273d+01
COFD(9718) = 5.35817304d+00
COFD(9719) = -4.69455306d-01
COFD(9720) = 1.99063292d-02
COFD(9721) = -2.11796273d+01
COFD(9722) = 5.35817304d+00
COFD(9723) = -4.69455306d-01
COFD(9724) = 1.99063292d-02
COFD(9725) = -2.11736657d+01
COFD(9726) = 5.35817304d+00
COFD(9727) = -4.69455306d-01
COFD(9728) = 1.99063292d-02
COFD(9729) = -1.79720978d+01
COFD(9730) = 4.30841971d+00
COFD(9731) = -3.45524579d-01
COFD(9732) = 1.50265381d-02
COFD(9733) = -1.91975423d+01
COFD(9734) = 4.73837374d+00
COFD(9735) = -3.97916523d-01
COFD(9736) = 1.71599378d-02
COFD(9737) = -2.17133615d+01
COFD(9738) = 5.41948945d+00
COFD(9739) = -4.73474138d-01
COFD(9740) = 1.99217500d-02
COFD(9741) = -2.17191046d+01
COFD(9742) = 5.41948945d+00
COFD(9743) = -4.73474138d-01
COFD(9744) = 1.99217500d-02
COFD(9745) = -2.22707633d+01
COFD(9746) = 5.57944619d+00
COFD(9747) = -4.78032626d-01
COFD(9748) = 1.94775370d-02
COFD(9749) = -2.22765150d+01
COFD(9750) = 5.57944619d+00
COFD(9751) = -4.78032626d-01
COFD(9752) = 1.94775370d-02
COFD(9753) = -1.36883939d+01
COFD(9754) = 3.19999740d+00
COFD(9755) = -2.04999020d-01
COFD(9756) = 9.06766105d-03
COFD(9757) = -1.74813786d+01
COFD(9758) = 4.80792005d+00
COFD(9759) = -4.06126584d-01
COFD(9760) = 1.74831083d-02
COFD(9761) = -1.69227183d+01
COFD(9762) = 4.09077642d+00
COFD(9763) = -3.18894990d-01
COFD(9764) = 1.39371445d-02
COFD(9765) = -1.83209411d+01
COFD(9766) = 4.43878381d+00
COFD(9767) = -3.61697624d-01
COFD(9768) = 1.56975581d-02
COFD(9769) = -1.69446537d+01
COFD(9770) = 4.09077642d+00
COFD(9771) = -3.18894990d-01
COFD(9772) = 1.39371445d-02
COFD(9773) = -2.13587539d+01
COFD(9774) = 5.35901292d+00
COFD(9775) = -4.36172487d-01
COFD(9776) = 1.71345319d-02
COFD(9777) = -1.83296851d+01
COFD(9778) = 4.43878381d+00
COFD(9779) = -3.61697624d-01
COFD(9780) = 1.56975581d-02
COFD(9781) = -1.83380528d+01
COFD(9782) = 4.43878381d+00
COFD(9783) = -3.61697624d-01
COFD(9784) = 1.56975581d-02
COFD(9785) = -1.65671124d+01
COFD(9786) = 3.93849401d+00
COFD(9787) = -2.99416642d-01
COFD(9788) = 1.31020815d-02
COFD(9789) = -1.68460201d+01
COFD(9790) = 4.09077642d+00
COFD(9791) = -3.18894990d-01
COFD(9792) = 1.39371445d-02
COFD(9793) = -1.91647170d+01
COFD(9794) = 4.80030220d+00
COFD(9795) = -4.05235041d-01
COFD(9796) = 1.74483531d-02
COFD(9797) = -1.91647170d+01
COFD(9798) = 4.80030220d+00
COFD(9799) = -4.05235041d-01
COFD(9800) = 1.74483531d-02
COFD(9801) = -1.91905015d+01
COFD(9802) = 4.80030220d+00
COFD(9803) = -4.05235041d-01
COFD(9804) = 1.74483531d-02
COFD(9805) = -1.91430978d+01
COFD(9806) = 4.78094221d+00
COFD(9807) = -4.02985837d-01
COFD(9808) = 1.73614221d-02
COFD(9809) = -1.79974471d+01
COFD(9810) = 4.31656593d+00
COFD(9811) = -3.46539554d-01
COFD(9812) = 1.50688196d-02
COFD(9813) = -2.12554255d+01
COFD(9814) = 5.38135645d+00
COFD(9815) = -4.71058360d-01
COFD(9816) = 1.99188046d-02
COFD(9817) = -2.19866916d+01
COFD(9818) = 5.49642957d+00
COFD(9819) = -4.61132993d-01
COFD(9820) = 1.85004773d-02
COFD(9821) = -2.19967195d+01
COFD(9822) = 5.49642957d+00
COFD(9823) = -4.61132993d-01
COFD(9824) = 1.85004773d-02
COFD(9825) = -2.20773284d+01
COFD(9826) = 5.58518321d+00
COFD(9827) = -4.80534235d-01
COFD(9828) = 1.96556393d-02
COFD(9829) = -2.20773284d+01
COFD(9830) = 5.58518321d+00
COFD(9831) = -4.80534235d-01
COFD(9832) = 1.96556393d-02
COFD(9833) = -2.20161635d+01
COFD(9834) = 5.50583166d+00
COFD(9835) = -4.63753262d-01
COFD(9836) = 1.86693462d-02
COFD(9837) = -2.07126500d+01
COFD(9838) = 5.23116678d+00
COFD(9839) = -4.54972675d-01
COFD(9840) = 1.93572354d-02
COFD(9841) = -2.07249293d+01
COFD(9842) = 5.23116678d+00
COFD(9843) = -4.54972675d-01
COFD(9844) = 1.93572354d-02
COFD(9845) = -2.07365685d+01
COFD(9846) = 5.23116678d+00
COFD(9847) = -4.54972675d-01
COFD(9848) = 1.93572354d-02
COFD(9849) = -2.15223412d+01
COFD(9850) = 5.47662534d+00
COFD(9851) = -4.79750189d-01
COFD(9852) = 2.01492627d-02
COFD(9853) = -2.13454345d+01
COFD(9854) = 5.39784038d+00
COFD(9855) = -4.72269379d-01
COFD(9856) = 1.99339592d-02
COFD(9857) = -2.13554415d+01
COFD(9858) = 5.39784038d+00
COFD(9859) = -4.72269379d-01
COFD(9860) = 1.99339592d-02
COFD(9861) = -1.92580496d+01
COFD(9862) = 4.84335373d+00
COFD(9863) = -4.10212751d-01
COFD(9864) = 1.76395874d-02
COFD(9865) = -2.22589130d+01
COFD(9866) = 5.57944619d+00
COFD(9867) = -4.78032626d-01
COFD(9868) = 1.94775370d-02
COFD(9869) = -2.22589130d+01
COFD(9870) = 5.57944619d+00
COFD(9871) = -4.78032626d-01
COFD(9872) = 1.94775370d-02
COFD(9873) = -1.66258278d+01
COFD(9874) = 3.93849401d+00
COFD(9875) = -2.99416642d-01
COFD(9876) = 1.31020815d-02
COFD(9877) = -1.68688781d+01
COFD(9878) = 4.09077642d+00
COFD(9879) = -3.18894990d-01
COFD(9880) = 1.39371445d-02
COFD(9881) = -1.68926059d+01
COFD(9882) = 4.09077642d+00
COFD(9883) = -3.18894990d-01
COFD(9884) = 1.39371445d-02
COFD(9885) = -2.16145219d+01
COFD(9886) = 5.50657794d+00
COFD(9887) = -4.63927523d-01
COFD(9888) = 1.86800112d-02
COFD(9889) = -1.70069808d+01
COFD(9890) = 3.93849401d+00
COFD(9891) = -2.99416642d-01
COFD(9892) = 1.31020815d-02
COFD(9893) = -1.79877202d+01
COFD(9894) = 4.30841971d+00
COFD(9895) = -3.45524579d-01
COFD(9896) = 1.50265381d-02
COFD(9897) = -2.05596852d+01
COFD(9898) = 5.17526774d+00
COFD(9899) = -4.48472252d-01
COFD(9900) = 1.91050891d-02
COFD(9901) = -2.11792397d+01
COFD(9902) = 5.35817304d+00
COFD(9903) = -4.69455306d-01
COFD(9904) = 1.99063292d-02
COFD(9905) = -1.85991417d+01
COFD(9906) = 4.53681554d+00
COFD(9907) = -3.73517224d-01
COFD(9908) = 1.61730724d-02
COFD(9909) = -1.71761068d+01
COFD(9910) = 4.00896225d+00
COFD(9911) = -3.08501394d-01
COFD(9912) = 1.34944985d-02
COFD(9913) = -2.18001733d+01
COFD(9914) = 5.36691865d+00
COFD(9915) = -4.37522912d-01
COFD(9916) = 1.72059807d-02
COFD(9917) = -2.18112354d+01
COFD(9918) = 5.36691865d+00
COFD(9919) = -4.37522912d-01
COFD(9920) = 1.72059807d-02
COFD(9921) = -1.92580672d+01
COFD(9922) = 4.84335373d+00
COFD(9923) = -4.10212751d-01
COFD(9924) = 1.76395874d-02
COFD(9925) = -2.11736657d+01
COFD(9926) = 5.35817304d+00
COFD(9927) = -4.69455306d-01
COFD(9928) = 1.99063292d-02
COFD(9929) = -2.11736657d+01
COFD(9930) = 5.35817304d+00
COFD(9931) = -4.69455306d-01
COFD(9932) = 1.99063292d-02
COFD(9933) = -2.11736657d+01
COFD(9934) = 5.35817304d+00
COFD(9935) = -4.69455306d-01
COFD(9936) = 1.99063292d-02
COFD(9937) = -2.11677742d+01
COFD(9938) = 5.35817304d+00
COFD(9939) = -4.69455306d-01
COFD(9940) = 1.99063292d-02
COFD(9941) = -1.79673900d+01
COFD(9942) = 4.30841971d+00
COFD(9943) = -3.45524579d-01
COFD(9944) = 1.50265381d-02
COFD(9945) = -1.91918004d+01
COFD(9946) = 4.73837374d+00
COFD(9947) = -3.97916523d-01
COFD(9948) = 1.71599378d-02
COFD(9949) = -2.17073954d+01
COFD(9950) = 5.41948945d+00
COFD(9951) = -4.73474138d-01
COFD(9952) = 1.99217500d-02
COFD(9953) = -2.17130700d+01
COFD(9954) = 5.41948945d+00
COFD(9955) = -4.73474138d-01
COFD(9956) = 1.99217500d-02
COFD(9957) = -2.22648002d+01
COFD(9958) = 5.57944619d+00
COFD(9959) = -4.78032626d-01
COFD(9960) = 1.94775370d-02
COFD(9961) = -2.22704834d+01
COFD(9962) = 5.57944619d+00
COFD(9963) = -4.78032626d-01
COFD(9964) = 1.94775370d-02
COFD(9965) = -1.16906297d+01
COFD(9966) = 2.47469981d+00
COFD(9967) = -1.10436257d-01
COFD(9968) = 4.95273813d-03
COFD(9969) = -1.42894441d+01
COFD(9970) = 3.67490723d+00
COFD(9971) = -2.65114792d-01
COFD(9972) = 1.16092671d-02
COFD(9973) = -1.40756935d+01
COFD(9974) = 3.07549274d+00
COFD(9975) = -1.88889344d-01
COFD(9976) = 8.37152866d-03
COFD(9977) = -1.52414485d+01
COFD(9978) = 3.35922578d+00
COFD(9979) = -2.25181399d-01
COFD(9980) = 9.92132878d-03
COFD(9981) = -1.40949196d+01
COFD(9982) = 3.07549274d+00
COFD(9983) = -1.88889344d-01
COFD(9984) = 8.37152866d-03
COFD(9985) = -2.10643259d+01
COFD(9986) = 5.53614847d+00
COFD(9987) = -4.86046736d-01
COFD(9988) = 2.03659188d-02
COFD(9989) = -1.52486273d+01
COFD(9990) = 3.35922578d+00
COFD(9991) = -2.25181399d-01
COFD(9992) = 9.92132878d-03
COFD(9993) = -1.52554761d+01
COFD(9994) = 3.35922578d+00
COFD(9995) = -2.25181399d-01
COFD(9996) = 9.92132878d-03
COFD(9997) = -1.38661480d+01
COFD(9998) = 2.97137588d+00
COFD(9999) = -1.75491257d-01
COFD(10000) = 7.79646773d-03
COFD(10001) = -1.40076852d+01
COFD(10002) = 3.07549274d+00
COFD(10003) = -1.88889344d-01
COFD(10004) = 8.37152866d-03
COFD(10005) = -1.59404882d+01
COFD(10006) = 3.66853818d+00
COFD(10007) = -2.64346221d-01
COFD(10008) = 1.15784613d-02
COFD(10009) = -1.59404882d+01
COFD(10010) = 3.66853818d+00
COFD(10011) = -2.64346221d-01
COFD(10012) = 1.15784613d-02
COFD(10013) = -1.59633387d+01
COFD(10014) = 3.66853818d+00
COFD(10015) = -2.64346221d-01
COFD(10016) = 1.15784613d-02
COFD(10017) = -1.59327297d+01
COFD(10018) = 3.65620899d+00
COFD(10019) = -2.62933804d-01
COFD(10020) = 1.15253223d-02
COFD(10021) = -1.50031687d+01
COFD(10022) = 3.26223357d+00
COFD(10023) = -2.12746642d-01
COFD(10024) = 9.38912883d-03
COFD(10025) = -1.81432461d+01
COFD(10026) = 4.37565431d+00
COFD(10027) = -3.53906025d-01
COFD(10028) = 1.53760786d-02
COFD(10029) = -2.04750581d+01
COFD(10030) = 5.23112374d+00
COFD(10031) = -4.54967682d-01
COFD(10032) = 1.93570423d-02
COFD(10033) = -2.04833713d+01
COFD(10034) = 5.23112374d+00
COFD(10035) = -4.54967682d-01
COFD(10036) = 1.93570423d-02
COFD(10037) = -2.02268902d+01
COFD(10038) = 5.13632093d+00
COFD(10039) = -4.44839124d-01
COFD(10040) = 1.90058354d-02
COFD(10041) = -2.02268902d+01
COFD(10042) = 5.13632093d+00
COFD(10043) = -4.44839124d-01
COFD(10044) = 1.90058354d-02
COFD(10045) = -2.03844252d+01
COFD(10046) = 5.18856872d+00
COFD(10047) = -4.50001829d-01
COFD(10048) = 1.91636142d-02
COFD(10049) = -1.75898751d+01
COFD(10050) = 4.19171952d+00
COFD(10051) = -3.31354810d-01
COFD(10052) = 1.44520623d-02
COFD(10053) = -1.76002031d+01
COFD(10054) = 4.19171952d+00
COFD(10055) = -3.31354810d-01
COFD(10056) = 1.44520623d-02
COFD(10057) = -1.76099552d+01
COFD(10058) = 4.19171952d+00
COFD(10059) = -3.31354810d-01
COFD(10060) = 1.44520623d-02
COFD(10061) = -1.85864144d+01
COFD(10062) = 4.54915847d+00
COFD(10063) = -3.75000738d-01
COFD(10064) = 1.62324821d-02
COFD(10065) = -1.83166353d+01
COFD(10066) = 4.42045763d+00
COFD(10067) = -3.59451578d-01
COFD(10068) = 1.56056164d-02
COFD(10069) = -1.83249299d+01
COFD(10070) = 4.42045763d+00
COFD(10071) = -3.59451578d-01
COFD(10072) = 1.56056164d-02
COFD(10073) = -1.59884305d+01
COFD(10074) = 3.72220402d+00
COFD(10075) = -2.71150591d-01
COFD(10076) = 1.18665265d-02
COFD(10077) = -2.02646611d+01
COFD(10078) = 5.10426133d+00
COFD(10079) = -4.41256919d-01
COFD(10080) = 1.88737290d-02
COFD(10081) = -2.02646611d+01
COFD(10082) = 5.10426133d+00
COFD(10083) = -4.41256919d-01
COFD(10084) = 1.88737290d-02
COFD(10085) = -1.39186707d+01
COFD(10086) = 2.97137588d+00
COFD(10087) = -1.75491257d-01
COFD(10088) = 7.79646773d-03
COFD(10089) = -1.40236044d+01
COFD(10090) = 3.07549274d+00
COFD(10091) = -1.88889344d-01
COFD(10092) = 8.37152866d-03
COFD(10093) = -1.40445141d+01
COFD(10094) = 3.07549274d+00
COFD(10095) = -1.88889344d-01
COFD(10096) = 8.37152866d-03
COFD(10097) = -2.04317360d+01
COFD(10098) = 5.33236114d+00
COFD(10099) = -4.66969586d-01
COFD(10100) = 1.98332508d-02
COFD(10101) = -1.42697086d+01
COFD(10102) = 2.97137588d+00
COFD(10103) = -1.75491257d-01
COFD(10104) = 7.79646773d-03
COFD(10105) = -1.49997274d+01
COFD(10106) = 3.25781069d+00
COFD(10107) = -2.12199367d-01
COFD(10108) = 9.36657283d-03
COFD(10109) = -1.74272299d+01
COFD(10110) = 4.14166966d+00
COFD(10111) = -3.25149462d-01
COFD(10112) = 1.41943811d-02
COFD(10113) = -1.79765417d+01
COFD(10114) = 4.30841971d+00
COFD(10115) = -3.45524579d-01
COFD(10116) = 1.50265381d-02
COFD(10117) = -1.55297765d+01
COFD(10118) = 3.46766436d+00
COFD(10119) = -2.39228775d-01
COFD(10120) = 1.05291747d-02
COFD(10121) = -1.43663445d+01
COFD(10122) = 3.01144603d+00
COFD(10123) = -1.80597248d-01
COFD(10124) = 8.01333257d-03
COFD(10125) = -2.08989475d+01
COFD(10126) = 5.37089858d+00
COFD(10127) = -4.70314009d-01
COFD(10128) = 1.99113180d-02
COFD(10129) = -2.09081824d+01
COFD(10130) = 5.37089858d+00
COFD(10131) = -4.70314009d-01
COFD(10132) = 1.99113180d-02
COFD(10133) = -1.59884445d+01
COFD(10134) = 3.72220402d+00
COFD(10135) = -2.71150591d-01
COFD(10136) = 1.18665265d-02
COFD(10137) = -1.79720978d+01
COFD(10138) = 4.30841971d+00
COFD(10139) = -3.45524579d-01
COFD(10140) = 1.50265381d-02
COFD(10141) = -1.79720978d+01
COFD(10142) = 4.30841971d+00
COFD(10143) = -3.45524579d-01
COFD(10144) = 1.50265381d-02
COFD(10145) = -1.79720978d+01
COFD(10146) = 4.30841971d+00
COFD(10147) = -3.45524579d-01
COFD(10148) = 1.50265381d-02
COFD(10149) = -1.79673900d+01
COFD(10150) = 4.30841971d+00
COFD(10151) = -3.45524579d-01
COFD(10152) = 1.50265381d-02
COFD(10153) = -1.49828430d+01
COFD(10154) = 3.25781069d+00
COFD(10155) = -2.12199367d-01
COFD(10156) = 9.36657283d-03
COFD(10157) = -1.59877782d+01
COFD(10158) = 3.63340763d+00
COFD(10159) = -2.60307961d-01
COFD(10160) = 1.14256954d-02
COFD(10161) = -1.87530122d+01
COFD(10162) = 4.48550694d+00
COFD(10163) = -3.67277454d-01
COFD(10164) = 1.59194755d-02
COFD(10165) = -1.87575355d+01
COFD(10166) = 4.48550694d+00
COFD(10167) = -3.67277454d-01
COFD(10168) = 1.59194755d-02
COFD(10169) = -2.02693653d+01
COFD(10170) = 5.10426133d+00
COFD(10171) = -4.41256919d-01
COFD(10172) = 1.88737290d-02
COFD(10173) = -2.02738958d+01
COFD(10174) = 5.10426133d+00
COFD(10175) = -4.41256919d-01
COFD(10176) = 1.88737290d-02
COFD(10177) = -1.23130152d+01
COFD(10178) = 2.74418790d+00
COFD(10179) = -1.46230156d-01
COFD(10180) = 6.53948886d-03
COFD(10181) = -1.54738604d+01
COFD(10182) = 4.15765300d+00
COFD(10183) = -3.27126237d-01
COFD(10184) = 1.42762611d-02
COFD(10185) = -1.49610527d+01
COFD(10186) = 3.41988961d+00
COFD(10187) = -2.33128386d-01
COFD(10188) = 1.02689994d-02
COFD(10189) = -1.62380075d+01
COFD(10190) = 3.72612300d+00
COFD(10191) = -2.71663673d-01
COFD(10192) = 1.18889643d-02
COFD(10193) = -1.49826725d+01
COFD(10194) = 3.41988961d+00
COFD(10195) = -2.33128386d-01
COFD(10196) = 1.02689994d-02
COFD(10197) = -2.12755888d+01
COFD(10198) = 5.60381989d+00
COFD(10199) = -4.91225459d-01
COFD(10200) = 2.04487844d-02
COFD(10201) = -1.62465583d+01
COFD(10202) = 3.72612300d+00
COFD(10203) = -2.71663673d-01
COFD(10204) = 1.18889643d-02
COFD(10205) = -1.62547381d+01
COFD(10206) = 3.72612300d+00
COFD(10207) = -2.71663673d-01
COFD(10208) = 1.18889643d-02
COFD(10209) = -1.46461881d+01
COFD(10210) = 3.27505697d+00
COFD(10211) = -2.14306851d-01
COFD(10212) = 9.45219335d-03
COFD(10213) = -1.48853569d+01
COFD(10214) = 3.41988961d+00
COFD(10215) = -2.33128386d-01
COFD(10216) = 1.02689994d-02
COFD(10217) = -1.71942502d+01
COFD(10218) = 4.14993355d+00
COFD(10219) = -3.26168062d-01
COFD(10220) = 1.42364115d-02
COFD(10221) = -1.71942502d+01
COFD(10222) = 4.14993355d+00
COFD(10223) = -3.26168062d-01
COFD(10224) = 1.42364115d-02
COFD(10225) = -1.72196961d+01
COFD(10226) = 4.14993355d+00
COFD(10227) = -3.26168062d-01
COFD(10228) = 1.42364115d-02
COFD(10229) = -1.71754154d+01
COFD(10230) = 4.13131681d+00
COFD(10231) = -3.23897559d-01
COFD(10232) = 1.41438222d-02
COFD(10233) = -1.60074211d+01
COFD(10234) = 3.63723937d+00
COFD(10235) = -2.60754222d-01
COFD(10236) = 1.14428814d-02
COFD(10237) = -1.93483692d+01
COFD(10238) = 4.79506290d+00
COFD(10239) = -4.04621659d-01
COFD(10240) = 1.74244230d-02
COFD(10241) = -2.14255087d+01
COFD(10242) = 5.52240865d+00
COFD(10243) = -4.84699537d-01
COFD(10244) = 2.03247833d-02
COFD(10245) = -2.14353267d+01
COFD(10246) = 5.52240865d+00
COFD(10247) = -4.84699537d-01
COFD(10248) = 2.03247833d-02
COFD(10249) = -2.10026861d+01
COFD(10250) = 5.38326647d+00
COFD(10251) = -4.71201048d-01
COFD(10252) = 1.99207516d-02
COFD(10253) = -2.10026861d+01
COFD(10254) = 5.38326647d+00
COFD(10255) = -4.71201048d-01
COFD(10256) = 1.99207516d-02
COFD(10257) = -2.13502344d+01
COFD(10258) = 5.48617067d+00
COFD(10259) = -4.80816776d-01
COFD(10260) = 2.01887868d-02
COFD(10261) = -1.87596895d+01
COFD(10262) = 4.61078776d+00
COFD(10263) = -3.82625667d-01
COFD(10264) = 1.65478601d-02
COFD(10265) = -1.87717330d+01
COFD(10266) = 4.61078776d+00
COFD(10267) = -3.82625667d-01
COFD(10268) = 1.65478601d-02
COFD(10269) = -1.87831434d+01
COFD(10270) = 4.61078776d+00
COFD(10271) = -3.82625667d-01
COFD(10272) = 1.65478601d-02
COFD(10273) = -1.98040322d+01
COFD(10274) = 4.97569695d+00
COFD(10275) = -4.26123307d-01
COFD(10276) = 1.82788664d-02
COFD(10277) = -1.94928815d+01
COFD(10278) = 4.83189721d+00
COFD(10279) = -4.08932249d-01
COFD(10280) = 1.75924650d-02
COFD(10281) = -1.95026789d+01
COFD(10282) = 4.83189721d+00
COFD(10283) = -4.08932249d-01
COFD(10284) = 1.75924650d-02
COFD(10285) = -1.72570607d+01
COFD(10286) = 4.19757624d+00
COFD(10287) = -3.32087529d-01
COFD(10288) = 1.44827462d-02
COFD(10289) = -2.12450965d+01
COFD(10290) = 5.40444222d+00
COFD(10291) = -4.72708609d-01
COFD(10292) = 1.99362392d-02
COFD(10293) = -2.12450965d+01
COFD(10294) = 5.40444222d+00
COFD(10295) = -4.72708609d-01
COFD(10296) = 1.99362392d-02
COFD(10297) = -1.47041948d+01
COFD(10298) = 3.27505697d+00
COFD(10299) = -2.14306851d-01
COFD(10300) = 9.45219335d-03
COFD(10301) = -1.49050015d+01
COFD(10302) = 3.41988961d+00
COFD(10303) = -2.33128386d-01
COFD(10304) = 1.02689994d-02
COFD(10305) = -1.49284025d+01
COFD(10306) = 3.41988961d+00
COFD(10307) = -2.33128386d-01
COFD(10308) = 1.02689994d-02
COFD(10309) = -2.09309775d+01
COFD(10310) = 5.48417182d+00
COFD(10311) = -4.80595673d-01
COFD(10312) = 2.01807046d-02
COFD(10313) = -1.50911395d+01
COFD(10314) = 3.27505697d+00
COFD(10315) = -2.14306851d-01
COFD(10316) = 9.45219335d-03
COFD(10317) = -1.60076874d+01
COFD(10318) = 3.63340763d+00
COFD(10319) = -2.60307961d-01
COFD(10320) = 1.14256954d-02
COFD(10321) = -1.85744661d+01
COFD(10322) = 4.54507313d+00
COFD(10323) = -3.74508073d-01
COFD(10324) = 1.62126770d-02
COFD(10325) = -1.92029732d+01
COFD(10326) = 4.73837374d+00
COFD(10327) = -3.97916523d-01
COFD(10328) = 1.71599378d-02
COFD(10329) = -1.66272428d+01
COFD(10330) = 3.87567029d+00
COFD(10331) = -2.91269312d-01
COFD(10332) = 1.27483193d-02
COFD(10333) = -1.52099661d+01
COFD(10334) = 3.32517045d+00
COFD(10335) = -2.20697887d-01
COFD(10336) = 9.72421963d-03
COFD(10337) = -2.17175057d+01
COFD(10338) = 5.60287101d+00
COFD(10339) = -4.91302853d-01
COFD(10340) = 2.04601038d-02
COFD(10341) = -2.17283455d+01
COFD(10342) = 5.60287101d+00
COFD(10343) = -4.91302853d-01
COFD(10344) = 2.04601038d-02
COFD(10345) = -1.72570778d+01
COFD(10346) = 4.19757624d+00
COFD(10347) = -3.32087529d-01
COFD(10348) = 1.44827462d-02
COFD(10349) = -1.91975423d+01
COFD(10350) = 4.73837374d+00
COFD(10351) = -3.97916523d-01
COFD(10352) = 1.71599378d-02
COFD(10353) = -1.91975423d+01
COFD(10354) = 4.73837374d+00
COFD(10355) = -3.97916523d-01
COFD(10356) = 1.71599378d-02
COFD(10357) = -1.91975423d+01
COFD(10358) = 4.73837374d+00
COFD(10359) = -3.97916523d-01
COFD(10360) = 1.71599378d-02
COFD(10361) = -1.91918004d+01
COFD(10362) = 4.73837374d+00
COFD(10363) = -3.97916523d-01
COFD(10364) = 1.71599378d-02
COFD(10365) = -1.59877782d+01
COFD(10366) = 3.63340763d+00
COFD(10367) = -2.60307961d-01
COFD(10368) = 1.14256954d-02
COFD(10369) = -1.72273911d+01
COFD(10370) = 4.09361913d+00
COFD(10371) = -3.19258125d-01
COFD(10372) = 1.39526981d-02
COFD(10373) = -1.99804429d+01
COFD(10374) = 4.90643348d+00
COFD(10375) = -4.17847566d-01
COFD(10376) = 1.79485301d-02
COFD(10377) = -1.99859716d+01
COFD(10378) = 4.90643348d+00
COFD(10379) = -4.17847566d-01
COFD(10380) = 1.79485301d-02
COFD(10381) = -2.12508342d+01
COFD(10382) = 5.40444222d+00
COFD(10383) = -4.72708609d-01
COFD(10384) = 1.99362392d-02
COFD(10385) = -2.12563714d+01
COFD(10386) = 5.40444222d+00
COFD(10387) = -4.72708609d-01
COFD(10388) = 1.99362392d-02
COFD(10389) = -1.43140669d+01
COFD(10390) = 3.31177824d+00
COFD(10391) = -2.18945280d-01
COFD(10392) = 9.64764419d-03
COFD(10393) = -1.83545293d+01
COFD(10394) = 4.98756925d+00
COFD(10395) = -4.27526072d-01
COFD(10396) = 1.83341755d-02
COFD(10397) = -1.76841045d+01
COFD(10398) = 4.24719726d+00
COFD(10399) = -3.38206061d-01
COFD(10400) = 1.47350654d-02
COFD(10401) = -1.91225244d+01
COFD(10402) = 4.61801405d+00
COFD(10403) = -3.83535652d-01
COFD(10404) = 1.65862513d-02
COFD(10405) = -1.77061949d+01
COFD(10406) = 4.24719726d+00
COFD(10407) = -3.38206061d-01
COFD(10408) = 1.47350654d-02
COFD(10409) = -2.13919273d+01
COFD(10410) = 5.17440955d+00
COFD(10411) = -4.04678430d-01
COFD(10412) = 1.54706350d-02
COFD(10413) = -1.91313643d+01
COFD(10414) = 4.61801405d+00
COFD(10415) = -3.83535652d-01
COFD(10416) = 1.65862513d-02
COFD(10417) = -1.91398254d+01
COFD(10418) = 4.61801405d+00
COFD(10419) = -3.83535652d-01
COFD(10420) = 1.65862513d-02
COFD(10421) = -1.73541877d+01
COFD(10422) = 4.11838823d+00
COFD(10423) = -3.22329602d-01
COFD(10424) = 1.40802468d-02
COFD(10425) = -1.76069153d+01
COFD(10426) = 4.24719726d+00
COFD(10427) = -3.38206061d-01
COFD(10428) = 1.47350654d-02
COFD(10429) = -1.99607067d+01
COFD(10430) = 4.97875278d+00
COFD(10431) = -4.26485475d-01
COFD(10432) = 1.82931933d-02
COFD(10433) = -1.99607067d+01
COFD(10434) = 4.97875278d+00
COFD(10435) = -4.26485475d-01
COFD(10436) = 1.82931933d-02
COFD(10437) = -1.99866570d+01
COFD(10438) = 4.97875278d+00
COFD(10439) = -4.26485475d-01
COFD(10440) = 1.82931933d-02
COFD(10441) = -1.99301980d+01
COFD(10442) = 4.95514826d+00
COFD(10443) = -4.23691395d-01
COFD(10444) = 1.81828318d-02
COFD(10445) = -1.87780799d+01
COFD(10446) = 4.49191492d+00
COFD(10447) = -3.68041771d-01
COFD(10448) = 1.59498676d-02
COFD(10449) = -2.18713229d+01
COFD(10450) = 5.47368915d+00
COFD(10451) = -4.79424291d-01
COFD(10452) = 2.01372920d-02
COFD(10453) = -2.20946362d+01
COFD(10454) = 5.36053938d+00
COFD(10455) = -4.36434519d-01
COFD(10456) = 1.71484255d-02
COFD(10457) = -2.21047682d+01
COFD(10458) = 5.36053938d+00
COFD(10459) = -4.36434519d-01
COFD(10460) = 1.71484255d-02
COFD(10461) = -2.23155854d+01
COFD(10462) = 5.50872500d+00
COFD(10463) = -4.64427323d-01
COFD(10464) = 1.87105612d-02
COFD(10465) = -2.23155854d+01
COFD(10466) = 5.50872500d+00
COFD(10467) = -4.64427323d-01
COFD(10468) = 1.87105612d-02
COFD(10469) = -2.21897029d+01
COFD(10470) = 5.39861129d+00
COFD(10471) = -4.43119198d-01
COFD(10472) = 1.75075657d-02
COFD(10473) = -2.14290825d+01
COFD(10474) = 5.37331605d+00
COFD(10475) = -4.70491203d-01
COFD(10476) = 1.99134666d-02
COFD(10477) = -2.14414783d+01
COFD(10478) = 5.37331605d+00
COFD(10479) = -4.70491203d-01
COFD(10480) = 1.99134666d-02
COFD(10481) = -2.14532306d+01
COFD(10482) = 5.37331605d+00
COFD(10483) = -4.70491203d-01
COFD(10484) = 1.99134666d-02
COFD(10485) = -2.21980988d+01
COFD(10486) = 5.59472344d+00
COFD(10487) = -4.91421518d-01
COFD(10488) = 2.05117088d-02
COFD(10489) = -2.19981673d+01
COFD(10490) = 5.51276597d+00
COFD(10491) = -4.83701824d-01
COFD(10492) = 2.02915297d-02
COFD(10493) = -2.20082782d+01
COFD(10494) = 5.51276597d+00
COFD(10495) = -4.83701824d-01
COFD(10496) = 2.02915297d-02
COFD(10497) = -2.01341439d+01
COFD(10498) = 5.03101171d+00
COFD(10499) = -4.32665019d-01
COFD(10500) = 1.85372086d-02
COFD(10501) = -2.24674263d+01
COFD(10502) = 5.49330641d+00
COFD(10503) = -4.60498247d-01
COFD(10504) = 1.84639199d-02
COFD(10505) = -2.24674263d+01
COFD(10506) = 5.49330641d+00
COFD(10507) = -4.60498247d-01
COFD(10508) = 1.84639199d-02
COFD(10509) = -1.74132498d+01
COFD(10510) = 4.11838823d+00
COFD(10511) = -3.22329602d-01
COFD(10512) = 1.40802468d-02
COFD(10513) = -1.76347103d+01
COFD(10514) = 4.24719726d+00
COFD(10515) = -3.38206061d-01
COFD(10516) = 1.47350654d-02
COFD(10517) = -1.76585983d+01
COFD(10518) = 4.24719726d+00
COFD(10519) = -3.38206061d-01
COFD(10520) = 1.47350654d-02
COFD(10521) = -2.18154539d+01
COFD(10522) = 5.39990556d+00
COFD(10523) = -4.43366024d-01
COFD(10524) = 1.75213695d-02
COFD(10525) = -1.77780309d+01
COFD(10526) = 4.11838823d+00
COFD(10527) = -3.22329602d-01
COFD(10528) = 1.40802468d-02
COFD(10529) = -1.87735510d+01
COFD(10530) = 4.48550694d+00
COFD(10531) = -3.67277454d-01
COFD(10532) = 1.59194755d-02
COFD(10533) = -2.13758293d+01
COFD(10534) = 5.35135690d+00
COFD(10535) = -4.68903699d-01
COFD(10536) = 1.98958921d-02
COFD(10537) = -2.17190071d+01
COFD(10538) = 5.41948945d+00
COFD(10539) = -4.73474138d-01
COFD(10540) = 1.99217500d-02
COFD(10541) = -1.94052185d+01
COFD(10542) = 4.71746683d+00
COFD(10543) = -3.95451856d-01
COFD(10544) = 1.70631070d-02
COFD(10545) = -1.79018076d+01
COFD(10546) = 4.16975222d+00
COFD(10547) = -3.28627847d-01
COFD(10548) = 1.43387221d-02
COFD(10549) = -2.17968118d+01
COFD(10550) = 5.18435436d+00
COFD(10551) = -4.06325299d-01
COFD(10552) = 1.55562327d-02
COFD(10553) = -2.18079839d+01
COFD(10554) = 5.18435436d+00
COFD(10555) = -4.06325299d-01
COFD(10556) = 1.55562327d-02
COFD(10557) = -2.01341617d+01
COFD(10558) = 5.03101171d+00
COFD(10559) = -4.32665019d-01
COFD(10560) = 1.85372086d-02
COFD(10561) = -2.17133615d+01
COFD(10562) = 5.41948945d+00
COFD(10563) = -4.73474138d-01
COFD(10564) = 1.99217500d-02
COFD(10565) = -2.17133615d+01
COFD(10566) = 5.41948945d+00
COFD(10567) = -4.73474138d-01
COFD(10568) = 1.99217500d-02
COFD(10569) = -2.17133615d+01
COFD(10570) = 5.41948945d+00
COFD(10571) = -4.73474138d-01
COFD(10572) = 1.99217500d-02
COFD(10573) = -2.17073954d+01
COFD(10574) = 5.41948945d+00
COFD(10575) = -4.73474138d-01
COFD(10576) = 1.99217500d-02
COFD(10577) = -1.87530122d+01
COFD(10578) = 4.48550694d+00
COFD(10579) = -3.67277454d-01
COFD(10580) = 1.59194755d-02
COFD(10581) = -1.99804429d+01
COFD(10582) = 4.90643348d+00
COFD(10583) = -4.17847566d-01
COFD(10584) = 1.79485301d-02
COFD(10585) = -2.24082390d+01
COFD(10586) = 5.56066804d+00
COFD(10587) = -4.88405706d-01
COFD(10588) = 2.04357330d-02
COFD(10589) = -2.24139864d+01
COFD(10590) = 5.56066804d+00
COFD(10591) = -4.88405706d-01
COFD(10592) = 2.04357330d-02
COFD(10593) = -2.24733881d+01
COFD(10594) = 5.49330641d+00
COFD(10595) = -4.60498247d-01
COFD(10596) = 1.84639199d-02
COFD(10597) = -2.24791442d+01
COFD(10598) = 5.49330641d+00
COFD(10599) = -4.60498247d-01
COFD(10600) = 1.84639199d-02
COFD(10601) = -1.43145779d+01
COFD(10602) = 3.31177824d+00
COFD(10603) = -2.18945280d-01
COFD(10604) = 9.64764419d-03
COFD(10605) = -1.83547906d+01
COFD(10606) = 4.98756925d+00
COFD(10607) = -4.27526072d-01
COFD(10608) = 1.83341755d-02
COFD(10609) = -1.76872087d+01
COFD(10610) = 4.24719726d+00
COFD(10611) = -3.38206061d-01
COFD(10612) = 1.47350654d-02
COFD(10613) = -1.91274187d+01
COFD(10614) = 4.61801405d+00
COFD(10615) = -3.83535652d-01
COFD(10616) = 1.65862513d-02
COFD(10617) = -1.77094398d+01
COFD(10618) = 4.24719726d+00
COFD(10619) = -3.38206061d-01
COFD(10620) = 1.47350654d-02
COFD(10621) = -2.13953083d+01
COFD(10622) = 5.17440955d+00
COFD(10623) = -4.04678430d-01
COFD(10624) = 1.54706350d-02
COFD(10625) = -1.91363464d+01
COFD(10626) = 4.61801405d+00
COFD(10627) = -3.83535652d-01
COFD(10628) = 1.65862513d-02
COFD(10629) = -1.91448929d+01
COFD(10630) = 4.61801405d+00
COFD(10631) = -3.83535652d-01
COFD(10632) = 1.65862513d-02
COFD(10633) = -1.73566853d+01
COFD(10634) = 4.11838823d+00
COFD(10635) = -3.22329602d-01
COFD(10636) = 1.40802468d-02
COFD(10637) = -1.76095743d+01
COFD(10638) = 4.24719726d+00
COFD(10639) = -3.38206061d-01
COFD(10640) = 1.47350654d-02
COFD(10641) = -1.99635214d+01
COFD(10642) = 4.97875278d+00
COFD(10643) = -4.26485475d-01
COFD(10644) = 1.82931933d-02
COFD(10645) = -1.99635214d+01
COFD(10646) = 4.97875278d+00
COFD(10647) = -4.26485475d-01
COFD(10648) = 1.82931933d-02
COFD(10649) = -1.99896221d+01
COFD(10650) = 4.97875278d+00
COFD(10651) = -4.26485475d-01
COFD(10652) = 1.82931933d-02
COFD(10653) = -1.99333084d+01
COFD(10654) = 4.95514826d+00
COFD(10655) = -4.23691395d-01
COFD(10656) = 1.81828318d-02
COFD(10657) = -1.87826029d+01
COFD(10658) = 4.49191492d+00
COFD(10659) = -3.68041771d-01
COFD(10660) = 1.59498676d-02
COFD(10661) = -2.18771314d+01
COFD(10662) = 5.47368915d+00
COFD(10663) = -4.79424291d-01
COFD(10664) = 2.01372920d-02
COFD(10665) = -2.20992569d+01
COFD(10666) = 5.36053938d+00
COFD(10667) = -4.36434519d-01
COFD(10668) = 1.71484255d-02
COFD(10669) = -2.21094839d+01
COFD(10670) = 5.36053938d+00
COFD(10671) = -4.36434519d-01
COFD(10672) = 1.71484255d-02
COFD(10673) = -2.23203936d+01
COFD(10674) = 5.50872500d+00
COFD(10675) = -4.64427323d-01
COFD(10676) = 1.87105612d-02
COFD(10677) = -2.23203936d+01
COFD(10678) = 5.50872500d+00
COFD(10679) = -4.64427323d-01
COFD(10680) = 1.87105612d-02
COFD(10681) = -2.21946011d+01
COFD(10682) = 5.39861129d+00
COFD(10683) = -4.43119198d-01
COFD(10684) = 1.75075657d-02
COFD(10685) = -2.14332997d+01
COFD(10686) = 5.37331605d+00
COFD(10687) = -4.70491203d-01
COFD(10688) = 1.99134666d-02
COFD(10689) = -2.14458019d+01
COFD(10690) = 5.37331605d+00
COFD(10691) = -4.70491203d-01
COFD(10692) = 1.99134666d-02
COFD(10693) = -2.14576575d+01
COFD(10694) = 5.37331605d+00
COFD(10695) = -4.70491203d-01
COFD(10696) = 1.99134666d-02
COFD(10697) = -2.22026260d+01
COFD(10698) = 5.59472344d+00
COFD(10699) = -4.91421518d-01
COFD(10700) = 2.05117088d-02
COFD(10701) = -2.20027922d+01
COFD(10702) = 5.51276597d+00
COFD(10703) = -4.83701824d-01
COFD(10704) = 2.02915297d-02
COFD(10705) = -2.20129980d+01
COFD(10706) = 5.51276597d+00
COFD(10707) = -4.83701824d-01
COFD(10708) = 2.02915297d-02
COFD(10709) = -2.01397498d+01
COFD(10710) = 5.03101171d+00
COFD(10711) = -4.32665019d-01
COFD(10712) = 1.85372086d-02
COFD(10713) = -2.24731023d+01
COFD(10714) = 5.49330641d+00
COFD(10715) = -4.60498247d-01
COFD(10716) = 1.84639199d-02
COFD(10717) = -2.24731023d+01
COFD(10718) = 5.49330641d+00
COFD(10719) = -4.60498247d-01
COFD(10720) = 1.84639199d-02
COFD(10721) = -1.74160614d+01
COFD(10722) = 4.11838823d+00
COFD(10723) = -3.22329602d-01
COFD(10724) = 1.40802468d-02
COFD(10725) = -1.76376724d+01
COFD(10726) = 4.24719726d+00
COFD(10727) = -3.38206061d-01
COFD(10728) = 1.47350654d-02
COFD(10729) = -1.76617059d+01
COFD(10730) = 4.24719726d+00
COFD(10731) = -3.38206061d-01
COFD(10732) = 1.47350654d-02
COFD(10733) = -2.18187020d+01
COFD(10734) = 5.39990556d+00
COFD(10735) = -4.43366024d-01
COFD(10736) = 1.75213695d-02
COFD(10737) = -1.77826519d+01
COFD(10738) = 4.11838823d+00
COFD(10739) = -3.22329602d-01
COFD(10740) = 1.40802468d-02
COFD(10741) = -1.87782648d+01
COFD(10742) = 4.48550694d+00
COFD(10743) = -3.67277454d-01
COFD(10744) = 1.59194755d-02
COFD(10745) = -2.13817659d+01
COFD(10746) = 5.35135690d+00
COFD(10747) = -4.68903699d-01
COFD(10748) = 1.98958921d-02
COFD(10749) = -2.17248157d+01
COFD(10750) = 5.41948945d+00
COFD(10751) = -4.73474138d-01
COFD(10752) = 1.99217500d-02
COFD(10753) = -1.94100248d+01
COFD(10754) = 4.71746683d+00
COFD(10755) = -3.95451856d-01
COFD(10756) = 1.70631070d-02
COFD(10757) = -1.79061291d+01
COFD(10758) = 4.16975222d+00
COFD(10759) = -3.28627847d-01
COFD(10760) = 1.43387221d-02
COFD(10761) = -2.18012366d+01
COFD(10762) = 5.18435436d+00
COFD(10763) = -4.06325299d-01
COFD(10764) = 1.55562327d-02
COFD(10765) = -2.18125092d+01
COFD(10766) = 5.18435436d+00
COFD(10767) = -4.06325299d-01
COFD(10768) = 1.55562327d-02
COFD(10769) = -2.01397678d+01
COFD(10770) = 5.03101171d+00
COFD(10771) = -4.32665019d-01
COFD(10772) = 1.85372086d-02
COFD(10773) = -2.17191046d+01
COFD(10774) = 5.41948945d+00
COFD(10775) = -4.73474138d-01
COFD(10776) = 1.99217500d-02
COFD(10777) = -2.17191046d+01
COFD(10778) = 5.41948945d+00
COFD(10779) = -4.73474138d-01
COFD(10780) = 1.99217500d-02
COFD(10781) = -2.17191046d+01
COFD(10782) = 5.41948945d+00
COFD(10783) = -4.73474138d-01
COFD(10784) = 1.99217500d-02
COFD(10785) = -2.17130700d+01
COFD(10786) = 5.41948945d+00
COFD(10787) = -4.73474138d-01
COFD(10788) = 1.99217500d-02
COFD(10789) = -1.87575355d+01
COFD(10790) = 4.48550694d+00
COFD(10791) = -3.67277454d-01
COFD(10792) = 1.59194755d-02
COFD(10793) = -1.99859716d+01
COFD(10794) = 4.90643348d+00
COFD(10795) = -4.17847566d-01
COFD(10796) = 1.79485301d-02
COFD(10797) = -2.24139864d+01
COFD(10798) = 5.56066804d+00
COFD(10799) = -4.88405706d-01
COFD(10800) = 2.04357330d-02
COFD(10801) = -2.24198006d+01
COFD(10802) = 5.56066804d+00
COFD(10803) = -4.88405706d-01
COFD(10804) = 2.04357330d-02
COFD(10805) = -2.24791326d+01
COFD(10806) = 5.49330641d+00
COFD(10807) = -4.60498247d-01
COFD(10808) = 1.84639199d-02
COFD(10809) = -2.24849555d+01
COFD(10810) = 5.49330641d+00
COFD(10811) = -4.60498247d-01
COFD(10812) = 1.84639199d-02
COFD(10813) = -1.57040212d+01
COFD(10814) = 3.93614244d+00
COFD(10815) = -2.99111497d-01
COFD(10816) = 1.30888229d-02
COFD(10817) = -1.94691430d+01
COFD(10818) = 5.43830787d+00
COFD(10819) = -4.75472880d-01
COFD(10820) = 1.99909996d-02
COFD(10821) = -1.90915649d+01
COFD(10822) = 4.84384483d+00
COFD(10823) = -4.10265575d-01
COFD(10824) = 1.76414287d-02
COFD(10825) = -2.05235731d+01
COFD(10826) = 5.18417470d+00
COFD(10827) = -4.49491573d-01
COFD(10828) = 1.91438508d-02
COFD(10829) = -1.91136491d+01
COFD(10830) = 4.84384483d+00
COFD(10831) = -4.10265575d-01
COFD(10832) = 1.76414287d-02
COFD(10833) = -1.87419199d+01
COFD(10834) = 3.96926341d+00
COFD(10835) = -2.16412264d-01
COFD(10836) = 6.06012078d-03
COFD(10837) = -2.05324091d+01
COFD(10838) = 5.18417470d+00
COFD(10839) = -4.49491573d-01
COFD(10840) = 1.91438508d-02
COFD(10841) = -2.05408665d+01
COFD(10842) = 5.18417470d+00
COFD(10843) = -4.49491573d-01
COFD(10844) = 1.91438508d-02
COFD(10845) = -1.87714196d+01
COFD(10846) = 4.71729964d+00
COFD(10847) = -3.95432573d-01
COFD(10848) = 1.70623691d-02
COFD(10849) = -1.90143953d+01
COFD(10850) = 4.84384483d+00
COFD(10851) = -4.10265575d-01
COFD(10852) = 1.76414287d-02
COFD(10853) = -2.11378465d+01
COFD(10854) = 5.42846112d+00
COFD(10855) = -4.74321870d-01
COFD(10856) = 1.99459749d-02
COFD(10857) = -2.11378465d+01
COFD(10858) = 5.42846112d+00
COFD(10859) = -4.74321870d-01
COFD(10860) = 1.99459749d-02
COFD(10861) = -2.11637902d+01
COFD(10862) = 5.42846112d+00
COFD(10863) = -4.74321870d-01
COFD(10864) = 1.99459749d-02
COFD(10865) = -2.11341653d+01
COFD(10866) = 5.41773516d+00
COFD(10867) = -4.73414338d-01
COFD(10868) = 1.99258685d-02
COFD(10869) = -2.02969740d+01
COFD(10870) = 5.11106992d+00
COFD(10871) = -4.42047129d-01
COFD(10872) = 1.89042990d-02
COFD(10873) = -2.22176950d+01
COFD(10874) = 5.54251230d+00
COFD(10875) = -4.70946314d-01
COFD(10876) = 1.90785869d-02
COFD(10877) = -2.00963085d+01
COFD(10878) = 4.41511629d+00
COFD(10879) = -2.84086963d-01
COFD(10880) = 9.37586971d-03
COFD(10881) = -2.01064363d+01
COFD(10882) = 4.41511629d+00
COFD(10883) = -2.84086963d-01
COFD(10884) = 9.37586971d-03
COFD(10885) = -2.09272429d+01
COFD(10886) = 4.82184721d+00
COFD(10887) = -3.48128875d-01
COFD(10888) = 1.25918978d-02
COFD(10889) = -2.09272429d+01
COFD(10890) = 4.82184721d+00
COFD(10891) = -3.48128875d-01
COFD(10892) = 1.25918978d-02
COFD(10893) = -2.03087302d+01
COFD(10894) = 4.50250781d+00
COFD(10895) = -2.97622106d-01
COFD(10896) = 1.00481473d-02
COFD(10897) = -2.20986379d+01
COFD(10898) = 5.58360799d+00
COFD(10899) = -4.82701436d-01
COFD(10900) = 1.98437922d-02
COFD(10901) = -2.21110290d+01
COFD(10902) = 5.58360799d+00
COFD(10903) = -4.82701436d-01
COFD(10904) = 1.98437922d-02
COFD(10905) = -2.21227768d+01
COFD(10906) = 5.58360799d+00
COFD(10907) = -4.82701436d-01
COFD(10908) = 1.98437922d-02
COFD(10909) = -2.20500806d+01
COFD(10910) = 5.44448440d+00
COFD(10911) = -4.51529024d-01
COFD(10912) = 1.79698119d-02
COFD(10913) = -2.22107627d+01
COFD(10914) = 5.51722375d+00
COFD(10915) = -4.66081431d-01
COFD(10916) = 1.88044011d-02
COFD(10917) = -2.22208695d+01
COFD(10918) = 5.51722375d+00
COFD(10919) = -4.66081431d-01
COFD(10920) = 1.88044011d-02
COFD(10921) = -2.12680082d+01
COFD(10922) = 5.47935225d+00
COFD(10923) = -4.80056796d-01
COFD(10924) = 2.01607180d-02
COFD(10925) = -2.09061629d+01
COFD(10926) = 4.72895031d+00
COFD(10927) = -3.33332771d-01
COFD(10928) = 1.18431478d-02
COFD(10929) = -2.09061629d+01
COFD(10930) = 4.72895031d+00
COFD(10931) = -3.33332771d-01
COFD(10932) = 1.18431478d-02
COFD(10933) = -1.88304679d+01
COFD(10934) = 4.71729964d+00
COFD(10935) = -3.95432573d-01
COFD(10936) = 1.70623691d-02
COFD(10937) = -1.90382267d+01
COFD(10938) = 4.84384483d+00
COFD(10939) = -4.10265575d-01
COFD(10940) = 1.76414287d-02
COFD(10941) = -1.90621083d+01
COFD(10942) = 4.84384483d+00
COFD(10943) = -4.10265575d-01
COFD(10944) = 1.76414287d-02
COFD(10945) = -1.99174474d+01
COFD(10946) = 4.50672715d+00
COFD(10947) = -2.98278922d-01
COFD(10948) = 1.00808824d-02
COFD(10949) = -1.92109657d+01
COFD(10950) = 4.71729964d+00
COFD(10951) = -3.95432573d-01
COFD(10952) = 1.70623691d-02
COFD(10953) = -2.02898957d+01
COFD(10954) = 5.10426133d+00
COFD(10955) = -4.41256919d-01
COFD(10956) = 1.88737290d-02
COFD(10957) = -2.20909803d+01
COFD(10958) = 5.59135292d+00
COFD(10959) = -4.85565630d-01
COFD(10960) = 2.00429035d-02
COFD(10961) = -2.22764060d+01
COFD(10962) = 5.57944619d+00
COFD(10963) = -4.78032626d-01
COFD(10964) = 1.94775370d-02
COFD(10965) = -2.08531119d+01
COFD(10966) = 5.29805436d+00
COFD(10967) = -4.62913371d-01
COFD(10968) = 1.96725621d-02
COFD(10969) = -1.93596541d+01
COFD(10970) = 4.77546552d+00
COFD(10971) = -4.02345606d-01
COFD(10972) = 1.73365265d-02
COFD(10973) = -1.92114566d+01
COFD(10974) = 3.99103591d+00
COFD(10975) = -2.19677673d-01
COFD(10976) = 6.21875426d-03
COFD(10977) = -1.92226244d+01
COFD(10978) = 3.99103591d+00
COFD(10979) = -2.19677673d-01
COFD(10980) = 6.21875426d-03
COFD(10981) = -2.12680259d+01
COFD(10982) = 5.47935225d+00
COFD(10983) = -4.80056796d-01
COFD(10984) = 2.01607180d-02
COFD(10985) = -2.22707633d+01
COFD(10986) = 5.57944619d+00
COFD(10987) = -4.78032626d-01
COFD(10988) = 1.94775370d-02
COFD(10989) = -2.22707633d+01
COFD(10990) = 5.57944619d+00
COFD(10991) = -4.78032626d-01
COFD(10992) = 1.94775370d-02
COFD(10993) = -2.22707633d+01
COFD(10994) = 5.57944619d+00
COFD(10995) = -4.78032626d-01
COFD(10996) = 1.94775370d-02
COFD(10997) = -2.22648002d+01
COFD(10998) = 5.57944619d+00
COFD(10999) = -4.78032626d-01
COFD(11000) = 1.94775370d-02
COFD(11001) = -2.02693653d+01
COFD(11002) = 5.10426133d+00
COFD(11003) = -4.41256919d-01
COFD(11004) = 1.88737290d-02
COFD(11005) = -2.12508342d+01
COFD(11006) = 5.40444222d+00
COFD(11007) = -4.72708609d-01
COFD(11008) = 1.99362392d-02
COFD(11009) = -2.24733881d+01
COFD(11010) = 5.49330641d+00
COFD(11011) = -4.60498247d-01
COFD(11012) = 1.84639199d-02
COFD(11013) = -2.24791326d+01
COFD(11014) = 5.49330641d+00
COFD(11015) = -4.60498247d-01
COFD(11016) = 1.84639199d-02
COFD(11017) = -2.09121217d+01
COFD(11018) = 4.72895031d+00
COFD(11019) = -3.33332771d-01
COFD(11020) = 1.18431478d-02
COFD(11021) = -2.09178748d+01
COFD(11022) = 4.72895031d+00
COFD(11023) = -3.33332771d-01
COFD(11024) = 1.18431478d-02
COFD(11025) = -1.57045332d+01
COFD(11026) = 3.93614244d+00
COFD(11027) = -2.99111497d-01
COFD(11028) = 1.30888229d-02
COFD(11029) = -1.94694048d+01
COFD(11030) = 5.43830787d+00
COFD(11031) = -4.75472880d-01
COFD(11032) = 1.99909996d-02
COFD(11033) = -1.90946745d+01
COFD(11034) = 4.84384483d+00
COFD(11035) = -4.10265575d-01
COFD(11036) = 1.76414287d-02
COFD(11037) = -2.05284752d+01
COFD(11038) = 5.18417470d+00
COFD(11039) = -4.49491573d-01
COFD(11040) = 1.91438508d-02
COFD(11041) = -1.91168996d+01
COFD(11042) = 4.84384483d+00
COFD(11043) = -4.10265575d-01
COFD(11044) = 1.76414287d-02
COFD(11045) = -1.87453067d+01
COFD(11046) = 3.96926341d+00
COFD(11047) = -2.16412264d-01
COFD(11048) = 6.06012078d-03
COFD(11049) = -2.05373990d+01
COFD(11050) = 5.18417470d+00
COFD(11051) = -4.49491573d-01
COFD(11052) = 1.91438508d-02
COFD(11053) = -2.05459419d+01
COFD(11054) = 5.18417470d+00
COFD(11055) = -4.49491573d-01
COFD(11056) = 1.91438508d-02
COFD(11057) = -1.87739217d+01
COFD(11058) = 4.71729964d+00
COFD(11059) = -3.95432573d-01
COFD(11060) = 1.70623691d-02
COFD(11061) = -1.90170590d+01
COFD(11062) = 4.84384483d+00
COFD(11063) = -4.10265575d-01
COFD(11064) = 1.76414287d-02
COFD(11065) = -2.11406662d+01
COFD(11066) = 5.42846112d+00
COFD(11067) = -4.74321870d-01
COFD(11068) = 1.99459749d-02
COFD(11069) = -2.11406662d+01
COFD(11070) = 5.42846112d+00
COFD(11071) = -4.74321870d-01
COFD(11072) = 1.99459749d-02
COFD(11073) = -2.11667605d+01
COFD(11074) = 5.42846112d+00
COFD(11075) = -4.74321870d-01
COFD(11076) = 1.99459749d-02
COFD(11077) = -2.11372811d+01
COFD(11078) = 5.41773516d+00
COFD(11079) = -4.73414338d-01
COFD(11080) = 1.99258685d-02
COFD(11081) = -2.03015042d+01
COFD(11082) = 5.11106992d+00
COFD(11083) = -4.42047129d-01
COFD(11084) = 1.89042990d-02
COFD(11085) = -2.22235122d+01
COFD(11086) = 5.54251230d+00
COFD(11087) = -4.70946314d-01
COFD(11088) = 1.90785869d-02
COFD(11089) = -2.01009366d+01
COFD(11090) = 4.41511629d+00
COFD(11091) = -2.84086963d-01
COFD(11092) = 9.37586971d-03
COFD(11093) = -2.01111595d+01
COFD(11094) = 4.41511629d+00
COFD(11095) = -2.84086963d-01
COFD(11096) = 9.37586971d-03
COFD(11097) = -2.09320587d+01
COFD(11098) = 4.82184721d+00
COFD(11099) = -3.48128875d-01
COFD(11100) = 1.25918978d-02
COFD(11101) = -2.09320587d+01
COFD(11102) = 4.82184721d+00
COFD(11103) = -3.48128875d-01
COFD(11104) = 1.25918978d-02
COFD(11105) = -2.03136362d+01
COFD(11106) = 4.50250781d+00
COFD(11107) = -2.97622106d-01
COFD(11108) = 1.00481473d-02
COFD(11109) = -2.21028621d+01
COFD(11110) = 5.58360799d+00
COFD(11111) = -4.82701436d-01
COFD(11112) = 1.98437922d-02
COFD(11113) = -2.21153597d+01
COFD(11114) = 5.58360799d+00
COFD(11115) = -4.82701436d-01
COFD(11116) = 1.98437922d-02
COFD(11117) = -2.21272109d+01
COFD(11118) = 5.58360799d+00
COFD(11119) = -4.82701436d-01
COFD(11120) = 1.98437922d-02
COFD(11121) = -2.20546151d+01
COFD(11122) = 5.44448440d+00
COFD(11123) = -4.51529024d-01
COFD(11124) = 1.79698119d-02
COFD(11125) = -2.22153950d+01
COFD(11126) = 5.51722375d+00
COFD(11127) = -4.66081431d-01
COFD(11128) = 1.88044011d-02
COFD(11129) = -2.22255968d+01
COFD(11130) = 5.51722375d+00
COFD(11131) = -4.66081431d-01
COFD(11132) = 1.88044011d-02
COFD(11133) = -2.12736225d+01
COFD(11134) = 5.47935225d+00
COFD(11135) = -4.80056796d-01
COFD(11136) = 2.01607180d-02
COFD(11137) = -2.09118474d+01
COFD(11138) = 4.72895031d+00
COFD(11139) = -3.33332771d-01
COFD(11140) = 1.18431478d-02
COFD(11141) = -2.09118474d+01
COFD(11142) = 4.72895031d+00
COFD(11143) = -3.33332771d-01
COFD(11144) = 1.18431478d-02
COFD(11145) = -1.88332845d+01
COFD(11146) = 4.71729964d+00
COFD(11147) = -3.95432573d-01
COFD(11148) = 1.70623691d-02
COFD(11149) = -1.90411940d+01
COFD(11150) = 4.84384483d+00
COFD(11151) = -4.10265575d-01
COFD(11152) = 1.76414287d-02
COFD(11153) = -1.90652212d+01
COFD(11154) = 4.84384483d+00
COFD(11155) = -4.10265575d-01
COFD(11156) = 1.76414287d-02
COFD(11157) = -1.99207011d+01
COFD(11158) = 4.50672715d+00
COFD(11159) = -2.98278922d-01
COFD(11160) = 1.00808824d-02
COFD(11161) = -1.92155940d+01
COFD(11162) = 4.71729964d+00
COFD(11163) = -3.95432573d-01
COFD(11164) = 1.70623691d-02
COFD(11165) = -2.02946170d+01
COFD(11166) = 5.10426133d+00
COFD(11167) = -4.41256919d-01
COFD(11168) = 1.88737290d-02
COFD(11169) = -2.20969258d+01
COFD(11170) = 5.59135292d+00
COFD(11171) = -4.85565630d-01
COFD(11172) = 2.00429035d-02
COFD(11173) = -2.22822234d+01
COFD(11174) = 5.57944619d+00
COFD(11175) = -4.78032626d-01
COFD(11176) = 1.94775370d-02
COFD(11177) = -2.08579258d+01
COFD(11178) = 5.29805436d+00
COFD(11179) = -4.62913371d-01
COFD(11180) = 1.96725621d-02
COFD(11181) = -1.93639826d+01
COFD(11182) = 4.77546552d+00
COFD(11183) = -4.02345606d-01
COFD(11184) = 1.73365265d-02
COFD(11185) = -1.92158886d+01
COFD(11186) = 3.99103591d+00
COFD(11187) = -2.19677673d-01
COFD(11188) = 6.21875426d-03
COFD(11189) = -1.92271569d+01
COFD(11190) = 3.99103591d+00
COFD(11191) = -2.19677673d-01
COFD(11192) = 6.21875426d-03
COFD(11193) = -2.12736405d+01
COFD(11194) = 5.47935225d+00
COFD(11195) = -4.80056796d-01
COFD(11196) = 2.01607180d-02
COFD(11197) = -2.22765150d+01
COFD(11198) = 5.57944619d+00
COFD(11199) = -4.78032626d-01
COFD(11200) = 1.94775370d-02
COFD(11201) = -2.22765150d+01
COFD(11202) = 5.57944619d+00
COFD(11203) = -4.78032626d-01
COFD(11204) = 1.94775370d-02
COFD(11205) = -2.22765150d+01
COFD(11206) = 5.57944619d+00
COFD(11207) = -4.78032626d-01
COFD(11208) = 1.94775370d-02
COFD(11209) = -2.22704834d+01
COFD(11210) = 5.57944619d+00
COFD(11211) = -4.78032626d-01
COFD(11212) = 1.94775370d-02
COFD(11213) = -2.02738958d+01
COFD(11214) = 5.10426133d+00
COFD(11215) = -4.41256919d-01
COFD(11216) = 1.88737290d-02
COFD(11217) = -2.12563714d+01
COFD(11218) = 5.40444222d+00
COFD(11219) = -4.72708609d-01
COFD(11220) = 1.99362392d-02
COFD(11221) = -2.24791442d+01
COFD(11222) = 5.49330641d+00
COFD(11223) = -4.60498247d-01
COFD(11224) = 1.84639199d-02
COFD(11225) = -2.24849555d+01
COFD(11226) = 5.49330641d+00
COFD(11227) = -4.60498247d-01
COFD(11228) = 1.84639199d-02
COFD(11229) = -2.09178748d+01
COFD(11230) = 4.72895031d+00
COFD(11231) = -3.33332771d-01
COFD(11232) = 1.18431478d-02
COFD(11233) = -2.09236948d+01
COFD(11234) = 4.72895031d+00
COFD(11235) = -3.33332771d-01
COFD(11236) = 1.18431478d-02

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

double precision, intent(out) :: COFTD(424)

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
COFTD(29) = 4.29789463d-01
COFTD(30) = 1.21313199d-04
COFTD(31) = -7.73071792d-08
COFTD(32) = 1.53234639d-11
COFTD(33) = 3.82245475d-01
COFTD(34) = 1.47167068d-05
COFTD(35) = -9.75995257d-09
COFTD(36) = 2.83217152d-12
COFTD(37) = 3.83439056d-01
COFTD(38) = 3.62717894d-05
COFTD(39) = -2.40281409d-08
COFTD(40) = 5.52792966d-12
COFTD(41) = 3.24747031d-01
COFTD(42) = 1.77798548d-04
COFTD(43) = -1.08934732d-07
COFTD(44) = 2.03595881d-11
COFTD(45) = 3.24747031d-01
COFTD(46) = 1.77798548d-04
COFTD(47) = -1.08934732d-07
COFTD(48) = 2.03595881d-11
COFTD(49) = 3.31191185d-01
COFTD(50) = 1.81326714d-04
COFTD(51) = -1.11096391d-07
COFTD(52) = 2.07635959d-11
COFTD(53) = 3.39557243d-01
COFTD(54) = 1.79335036d-04
COFTD(55) = -1.10135705d-07
COFTD(56) = 2.06427239d-11
COFTD(57) = 4.30605547d-01
COFTD(58) = 9.35961902d-05
COFTD(59) = -6.03983623d-08
COFTD(60) = 1.23115170d-11
COFTD(61) = 2.93191523d-01
COFTD(62) = 4.01430006d-04
COFTD(63) = -2.30705763d-07
COFTD(64) = 4.05176586d-11
COFTD(65) = 1.22119780d-01
COFTD(66) = 6.18373616d-04
COFTD(67) = -3.28422593d-07
COFTD(68) = 5.44603522d-11
COFTD(69) = 1.22693382d-01
COFTD(70) = 6.21278143d-04
COFTD(71) = -3.29965208d-07
COFTD(72) = 5.47161548d-11
COFTD(73) = 1.40314191d-01
COFTD(74) = 6.01266129d-04
COFTD(75) = -3.21915137d-07
COFTD(76) = 5.36679068d-11
COFTD(77) = 1.40314191d-01
COFTD(78) = 6.01266129d-04
COFTD(79) = -3.21915137d-07
COFTD(80) = 5.36679068d-11
COFTD(81) = 1.31424053d-01
COFTD(82) = 6.16429134d-04
COFTD(83) = -3.28571348d-07
COFTD(84) = 5.46153434d-11
COFTD(85) = 3.03701584d-01
COFTD(86) = 3.22476070d-04
COFTD(87) = -1.88701794d-07
COFTD(88) = 3.36545091d-11
COFTD(89) = 3.05613225d-01
COFTD(90) = 3.24505886d-04
COFTD(91) = -1.89889572d-07
COFTD(92) = 3.38663465d-11
COFTD(93) = 3.07392263d-01
COFTD(94) = 3.26394901d-04
COFTD(95) = -1.90994958d-07
COFTD(96) = 3.40634894d-11
COFTD(97) = 2.49017478d-01
COFTD(98) = 4.29036573d-04
COFTD(99) = -2.42668617d-07
COFTD(100) = 4.20801371d-11
COFTD(101) = 2.72759599d-01
COFTD(102) = 3.94402719d-04
COFTD(103) = -2.25800520d-07
COFTD(104) = 3.95325634d-11
COFTD(105) = 2.74036956d-01
COFTD(106) = 3.96249742d-04
COFTD(107) = -2.26857964d-07
COFTD(108) = 3.97176979d-11
COFTD(109) = 3.86107464d-01
COFTD(110) = 2.28760446d-04
COFTD(111) = -1.39425040d-07
COFTD(112) = 2.58989754d-11
COFTD(113) = 1.59288984d-01
COFTD(114) = 6.02833801d-04
COFTD(115) = -3.24837576d-07
COFTD(116) = 5.43909010d-11
COFTD(117) = 1.59288984d-01
COFTD(118) = 6.02833801d-04
COFTD(119) = -3.24837576d-07
COFTD(120) = 5.43909010d-11
COFTD(121) = 4.01449248d-01
COFTD(122) = 1.54560650d-05
COFTD(123) = -1.02502865d-08
COFTD(124) = 2.97445804d-12
COFTD(125) = 3.99902403d-01
COFTD(126) = 3.78291557d-05
COFTD(127) = -2.50598137d-08
COFTD(128) = 5.76527696d-12
COFTD(129) = 4.06833564d-01
COFTD(130) = 3.84848155d-05
COFTD(131) = -2.54941536d-08
COFTD(132) = 5.86520149d-12
COFTD(133) = 8.77465301d-02
COFTD(134) = 5.89528909d-04
COFTD(135) = -3.09359591d-07
COFTD(136) = 5.08970415d-11
COFTD(137) = 4.66750762d-01
COFTD(138) = 1.79702170d-05
COFTD(139) = -1.19176435d-08
COFTD(140) = 3.45829657d-12
COFTD(141) = 4.35494021d-01
COFTD(142) = 9.29420846d-05
COFTD(143) = -6.00247192d-08
COFTD(144) = 1.22609983d-11
COFTD(145) = 3.35878046d-01
COFTD(146) = 3.31392140d-04
COFTD(147) = -1.94935193d-07
COFTD(148) = 3.49301020d-11
COFTD(149) = 3.03583862d-01
COFTD(150) = 3.83874304d-04
COFTD(151) = -2.21858899d-07
COFTD(152) = 3.91461168d-11
COFTD(153) = 4.13278096d-01
COFTD(154) = 1.43904794d-04
COFTD(155) = -9.06639714d-08
COFTD(156) = 1.76105782d-11
COFTD(157) = 4.54833899d-01
COFTD(158) = 2.80897652d-05
COFTD(159) = -1.86905534d-08
COFTD(160) = 4.71583220d-12
COFTD(161) = 9.06108625d-02
COFTD(162) = 6.50510405d-04
COFTD(163) = -3.40496245d-07
COFTD(164) = 5.59280162d-11
COFTD(165) = 9.11008480d-02
COFTD(166) = 6.54028093d-04
COFTD(167) = -3.42337506d-07
COFTD(168) = 5.62304515d-11
COFTD(169) = 3.86110105d-01
COFTD(170) = 2.28762012d-04
COFTD(171) = -1.39425994d-07
COFTD(172) = 2.58991526d-11
COFTD(173) = 3.02944815d-01
COFTD(174) = 3.83066246d-04
COFTD(175) = -2.21391884d-07
COFTD(176) = 3.90637139d-11
COFTD(177) = 3.02944815d-01
COFTD(178) = 3.83066246d-04
COFTD(179) = -2.21391884d-07
COFTD(180) = 3.90637139d-11
COFTD(181) = 3.02944815d-01
COFTD(182) = 3.83066246d-04
COFTD(183) = -2.21391884d-07
COFTD(184) = 3.90637139d-11
COFTD(185) = 3.02263016d-01
COFTD(186) = 3.82204128d-04
COFTD(187) = -2.20893626d-07
COFTD(188) = 3.89757982d-11
COFTD(189) = 4.31331269d-01
COFTD(190) = 9.20536800d-05
COFTD(191) = -5.94509616d-08
COFTD(192) = 1.21437993d-11
COFTD(193) = 4.01012808d-01
COFTD(194) = 1.97252826d-04
COFTD(195) = -1.21698146d-07
COFTD(196) = 2.29408847d-11
COFTD(197) = 2.73200750d-01
COFTD(198) = 4.32801981d-04
COFTD(199) = -2.46215081d-07
COFTD(200) = 4.28882665d-11
COFTD(201) = 2.73786959d-01
COFTD(202) = 4.33730648d-04
COFTD(203) = -2.46743387d-07
COFTD(204) = 4.29802922d-11
COFTD(205) = 1.59647939d-01
COFTD(206) = 6.04192274d-04
COFTD(207) = -3.25569591d-07
COFTD(208) = 5.45134698d-11
COFTD(209) = 1.59991186d-01
COFTD(210) = 6.05491303d-04
COFTD(211) = -3.26269573d-07
COFTD(212) = 5.46306751d-11
COFTD(213) = 1.44152190d-01
COFTD(214) = 7.99993584d-05
COFTD(215) = -4.89707442d-08
COFTD(216) = 9.14277269d-12
COFTD(217) = 0.00000000d+00
COFTD(218) = 0.00000000d+00
COFTD(219) = 0.00000000d+00
COFTD(220) = 0.00000000d+00
COFTD(221) = 2.35283119d-01
COFTD(222) = 4.65670599d-04
COFTD(223) = -2.60939824d-07
COFTD(224) = 4.49271822d-11
COFTD(225) = 1.79840299d-01
COFTD(226) = 6.01722902d-04
COFTD(227) = -3.26433894d-07
COFTD(228) = 5.49112302d-11
COFTD(229) = 2.37053352d-01
COFTD(230) = 4.69174231d-04
COFTD(231) = -2.62903094d-07
COFTD(232) = 4.52652072d-11
COFTD(233) = -1.74352698d-01
COFTD(234) = 8.62246873d-04
COFTD(235) = -3.79545489d-07
COFTD(236) = 5.60262093d-11
COFTD(237) = 1.80186965d-01
COFTD(238) = 6.02882805d-04
COFTD(239) = -3.27063140d-07
COFTD(240) = 5.50170790d-11
COFTD(241) = 1.80513677d-01
COFTD(242) = 6.03975942d-04
COFTD(243) = -3.27656165d-07
COFTD(244) = 5.51168351d-11
COFTD(245) = 2.49272491d-01
COFTD(246) = 4.08682510d-04
COFTD(247) = -2.31943878d-07
COFTD(248) = 4.03271405d-11
COFTD(249) = 2.28560867d-01
COFTD(250) = 4.52365967d-04
COFTD(251) = -2.53484536d-07
COFTD(252) = 4.36435719d-11
COFTD(253) = 9.90752318d-02
COFTD(254) = 6.44201384d-04
COFTD(255) = -3.38485953d-07
COFTD(256) = 5.57356746d-11
COFTD(257) = 9.90752318d-02
COFTD(258) = 6.44201384d-04
COFTD(259) = -3.38485953d-07
COFTD(260) = 5.57356746d-11
COFTD(261) = 1.00039110d-01
COFTD(262) = 6.50468660d-04
COFTD(263) = -3.41778999d-07
COFTD(264) = 5.62779132d-11
COFTD(265) = 1.05124122d-01
COFTD(266) = 6.50665957d-04
COFTD(267) = -3.42564538d-07
COFTD(268) = 5.64804120d-11
COFTD(269) = 2.00119897d-01
COFTD(270) = 5.64793704d-04
COFTD(271) = -3.09445484d-07
COFTD(272) = 5.24139335d-11
COFTD(273) = -2.00309448d-02
COFTD(274) = 8.50440115d-04
COFTD(275) = -4.21064468d-07
COFTD(276) = 6.67959710d-11
COFTD(277) = -1.60981264d-01
COFTD(278) = 9.03807572d-04
COFTD(279) = -4.06927941d-07
COFTD(280) = 6.09202254d-11
COFTD(281) = -1.61357564d-01
COFTD(282) = 9.05920260d-04
COFTD(283) = -4.07879153d-07
COFTD(284) = 6.10626290d-11
COFTD(285) = -1.31244519d-01
COFTD(286) = 9.03901384d-04
COFTD(287) = -4.17831507d-07
COFTD(288) = 6.35725667d-11
COFTD(289) = -1.31244519d-01
COFTD(290) = 9.03901384d-04
COFTD(291) = -4.17831507d-07
COFTD(292) = 6.35725667d-11
COFTD(293) = -1.56651581d-01
COFTD(294) = 9.09789751d-04
COFTD(295) = -4.11714242d-07
COFTD(296) = 6.18310893d-11
COFTD(297) = 1.62736132d-02
COFTD(298) = 7.87669911d-04
COFTD(299) = -3.97050662d-07
COFTD(300) = 6.36859622d-11
COFTD(301) = 1.63245097d-02
COFTD(302) = 7.90133388d-04
COFTD(303) = -3.98292458d-07
COFTD(304) = 6.38851432d-11
COFTD(305) = 1.63717489d-02
COFTD(306) = 7.92419842d-04
COFTD(307) = -3.99445020d-07
COFTD(308) = 6.40700113d-11
COFTD(309) = -5.08744745d-02
COFTD(310) = 8.54342586d-04
COFTD(311) = -4.15926453d-07
COFTD(312) = 6.53063261d-11
COFTD(313) = -2.71690558d-02
COFTD(314) = 8.37233133d-04
COFTD(315) = -4.12887636d-07
COFTD(316) = 6.53405197d-11
COFTD(317) = -2.72323768d-02
COFTD(318) = 8.39184413d-04
COFTD(319) = -4.13849924d-07
COFTD(320) = 6.54928043d-11
COFTD(321) = 9.86934401d-02
COFTD(322) = 7.20974863d-04
COFTD(323) = -3.77135221d-07
COFTD(324) = 6.19202579d-11
COFTD(325) = -1.41640506d-01
COFTD(326) = 9.21404324d-04
COFTD(327) = -4.23210110d-07
COFTD(328) = 6.41400322d-11
COFTD(329) = -1.41640506d-01
COFTD(330) = 9.21404324d-04
COFTD(331) = -4.23210110d-07
COFTD(332) = 6.41400322d-11
COFTD(333) = 2.55342377d-01
COFTD(334) = 4.18634095d-04
COFTD(335) = -2.37591806d-07
COFTD(336) = 4.13091228d-11
COFTD(337) = 2.33338617d-01
COFTD(338) = 4.61822055d-04
COFTD(339) = -2.58783281d-07
COFTD(340) = 4.45558806d-11
COFTD(341) = 2.35326294d-01
COFTD(342) = 4.65756050d-04
COFTD(343) = -2.60987707d-07
COFTD(344) = 4.49354264d-11
COFTD(345) = -1.47923508d-01
COFTD(346) = 8.60600299d-04
COFTD(347) = -3.89552943d-07
COFTD(348) = 5.85120641d-11
COFTD(349) = 2.75142577d-01
COFTD(350) = 4.51096544d-04
COFTD(351) = -2.56015560d-07
COFTD(352) = 4.45123862d-11
COFTD(353) = 2.02488208d-01
COFTD(354) = 5.65443199d-04
COFTD(355) = -3.09999001d-07
COFTD(356) = 5.25313540d-11
COFTD(357) = 2.76047996d-02
COFTD(358) = 8.06201373d-04
COFTD(359) = -4.08432273d-07
COFTD(360) = 6.57153466d-11
COFTD(361) = -8.57036272d-03
COFTD(362) = 8.40084342d-04
COFTD(363) = -4.18319087d-07
COFTD(364) = 6.65930617d-11
COFTD(365) = 1.59106721d-01
COFTD(366) = 6.28946469d-04
COFTD(367) = -3.38130427d-07
COFTD(368) = 5.65284977d-11
COFTD(369) = 2.61811629d-01
COFTD(370) = 4.65435180d-04
COFTD(371) = -2.62695240d-07
COFTD(372) = 4.54779557d-11
COFTD(373) = -1.80158237d-01
COFTD(374) = 8.95404521d-04
COFTD(375) = -3.94538141d-07
COFTD(376) = 5.82762309d-11
COFTD(377) = -1.80642727d-01
COFTD(378) = 8.97812486d-04
COFTD(379) = -3.95599152d-07
COFTD(380) = 5.84329501d-11
COFTD(381) = 9.86937771d-02
COFTD(382) = 7.20977325d-04
COFTD(383) = -3.77136509d-07
COFTD(384) = 6.19204694d-11
COFTD(385) = -8.56135217d-03
COFTD(386) = 8.39201108d-04
COFTD(387) = -4.17879282d-07
COFTD(388) = 6.65230483d-11
COFTD(389) = -8.56135217d-03
COFTD(390) = 8.39201108d-04
COFTD(391) = -4.17879282d-07
COFTD(392) = 6.65230483d-11
COFTD(393) = -8.56135217d-03
COFTD(394) = 8.39201108d-04
COFTD(395) = -4.17879282d-07
COFTD(396) = 6.65230483d-11
COFTD(397) = -8.55172903d-03
COFTD(398) = 8.38257829d-04
COFTD(399) = -4.17409577d-07
COFTD(400) = 6.64482750d-11
COFTD(401) = 2.01521643d-01
COFTD(402) = 5.62744089d-04
COFTD(403) = -3.08519239d-07
COFTD(404) = 5.22805986d-11
COFTD(405) = 1.22193921d-01
COFTD(406) = 6.90321128d-04
COFTD(407) = -3.64844875d-07
COFTD(408) = 6.03054876d-11
COFTD(409) = -4.06197305d-02
COFTD(410) = 8.67025869d-04
COFTD(411) = -4.24730138d-07
COFTD(412) = 6.69410527d-11
COFTD(413) = -4.06632162d-02
COFTD(414) = 8.67954070d-04
COFTD(415) = -4.25184836d-07
COFTD(416) = 6.70127170d-11
COFTD(417) = -1.41799739d-01
COFTD(418) = 9.22440172d-04
COFTD(419) = -4.23685885d-07
COFTD(420) = 6.42121388d-11
COFTD(421) = -1.41951848d-01
COFTD(422) = 9.23429679d-04
COFTD(423) = -4.24140376d-07
COFTD(424) = 6.42810196d-11

end subroutine

end module fuego_module



