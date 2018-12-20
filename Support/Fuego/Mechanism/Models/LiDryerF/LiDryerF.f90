module fuego_module

  implicit none
  private
  public :: ckcpms
  public :: ckcvms
  public :: ckxty
  public :: ckytx
  public :: ckhms
  public :: vckytx
  public :: vckhms

! Inverse molecular weights
double precision, parameter :: imw(9) = (/ &
    1.d0 / 2.015940d0,  & ! H2
    1.d0 / 31.998800d0,  & ! O2
    1.d0 / 18.015340d0,  & ! H2O
    1.d0 / 1.007970d0,  & ! H
    1.d0 / 15.999400d0,  & ! O
    1.d0 / 17.007370d0,  & ! OH
    1.d0 / 33.006770d0,  & ! HO2
    1.d0 / 34.014740d0,  & ! H2O2
    1.d0 / 28.013400d0/)  ! N2

contains

! convert y[species] (mass fracs) to x[species] (mole fracs)
subroutine ckytx(y, iwrk, rwrk, x)

    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: x(9)

    double precision :: YOW, YOWINV
    double precision :: tmp(9)
    integer :: i

    do i=1, 9
        tmp(i) = y(i)*imw(i)
    end do
    do i=1, 9
        YOW = YOW + tmp(i)
    end do

    YOWINV = 1.d0/YOW

    do i=1, 9
        x(i) = y(i)*imw(i)*YOWINV
    end do

end subroutine

! convert y(npoints,species) (mass fracs) to x(npoints,species) (mole fracs)
subroutine vckytx(np, y, iwrk, rwrk, x)

    integer, intent(in) :: np
    double precision, intent(in) :: y(np,9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: x(np,9)

    double precision :: YOW(np)
    integer :: i, n

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
        YOW(i) = 1.d0/YOW(i)
    end do

    do n=1, 9
        do i=1, np
            x(i,n) = x(i,n) * YOW(i)
        end do
    end do

end subroutine

! convert x[species] (mole fracs) to y[species] (mass fracs)
subroutine ckxty(x, iwrk, rwrk, y)

    double precision, intent(in) :: x(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: y(9)

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
    XWinv = 1.d0/XW
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
subroutine ckcvms(T, iwrk, rwrk, cvms)

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: cvms(9)

    double precision :: tc(5)
    double precision :: tT

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
subroutine ckcpms(T, iwrk, rwrk, cpms)

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: cpms(9)

    double precision :: tc(5)
    double precision :: tT

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

! Returns enthalpy in mass units (Eq 27.)
subroutine ckhms(T, iwrk, rwrk, hms)

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: hms(9)

    double precision :: tT, RT
    double precision :: tc(5), h(9)
    integer :: i

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451d+07*tT ! R*T

    call speciesEnthalpy(hms, tc)

    do i=1, 9
        hms(i) = hms(i) * (RT*imw(i))
    end do

end subroutine


! Returns enthalpy in mass units (Eq 27.)
subroutine vckhms(np, T, iwrk, rwrk, hms)
    integer, intent(in) :: np
    double precision, intent(in) :: T(np)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: hms(np,9)

    double precision :: tc(5), h(9)
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
    end do

    do n=1, 9
        do i=1, np
            hms(i,n) = hms(i,n) * (8.31451d+07 * T(i) * imw(n))
        end do
    end do

end subroutine


! compute Cv/R at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine cv_R(species, tc)

    double precision, intent(inout) :: species(9)
    double precision, intent(in) :: tc(5)
    ! temperature
    double precision :: T

    T = tc(2)

    ! species with midpoint at T=1000 kelvin
    if (T <  1.00000000d+03) then
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

    double precision, intent(inout) :: species(9)
    double precision, intent(in) :: tc(5)
    ! temperature
    double precision :: T

    T = tc(2)

    ! species with midpoint at T=1000 kelvin
    if (T <  1.00000000d+03) then
        ! species 1: H2
        species(1) = &
            +3.29812431d+00 &
            +8.24944174d-04 * tc(1) &
            -8.14301529d-07 * tc(2) &
            -9.47543433d-11 * tc(3) &
            +4.13487224d-13 * tc(4)
        ! species 2: O2
        species(2) = &
            +3.21293640d+00 &
            +1.12748635d-03 * tc(1) &
            -5.75615047d-07 * tc(2) &
            +1.31387723d-09 * tc(3) &
            -8.76855392d-13 * tc(4)
        ! species 3: H2O
        species(3) = &
            +3.38684249d+00 &
            +3.47498246d-03 * tc(1) &
            -6.35469633d-06 * tc(2) &
            +6.96858127d-09 * tc(3) &
            -2.50658847d-12 * tc(4)
        ! species 4: H
        species(4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(1) &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4)
        ! species 5: O
        species(5) = &
            +2.94642878d+00 &
            -1.63816649d-03 * tc(1) &
            +2.42103170d-06 * tc(2) &
            -1.60284319d-09 * tc(3) &
            +3.89069636d-13 * tc(4)
        ! species 6: OH
        species(6) = &
            +4.12530561d+00 &
            -3.22544939d-03 * tc(1) &
            +6.52764691d-06 * tc(2) &
            -5.79853643d-09 * tc(3) &
            +2.06237379d-12 * tc(4)
        ! species 7: HO2
        species(7) = &
            +4.30179801d+00 &
            -4.74912051d-03 * tc(1) &
            +2.11582891d-05 * tc(2) &
            -2.42763894d-08 * tc(3) &
            +9.29225124d-12 * tc(4)
        ! species 8: H2O2
        species(8) = &
            +3.38875365d+00 &
            +6.56922581d-03 * tc(1) &
            -1.48501258d-07 * tc(2) &
            -4.62580552d-09 * tc(3) &
            +2.47151475d-12 * tc(4)
        ! species 9: N2
        species(9) = &
            +3.29867700d+00 &
            +1.40824000d-03 * tc(1) &
            -3.96322200d-06 * tc(2) &
            +5.64151500d-09 * tc(3) &
            -2.44485500d-12 * tc(4)
    else
        !species 1: H2
        species(1) = &
            +2.99142337d+00 &
            +7.00064411d-04 * tc(1) &
            -5.63382869d-08 * tc(2) &
            -9.23157818d-12 * tc(3) &
            +1.58275179d-15 * tc(4)
        !species 2: O2
        species(2) = &
            +3.69757819d+00 &
            +6.13519689d-04 * tc(1) &
            -1.25884199d-07 * tc(2) &
            +1.77528148d-11 * tc(3) &
            -1.13643531d-15 * tc(4)
        !species 3: H2O
        species(3) = &
            +2.67214561d+00 &
            +3.05629289d-03 * tc(1) &
            -8.73026011d-07 * tc(2) &
            +1.20099639d-10 * tc(3) &
            -6.39161787d-15 * tc(4)
        !species 4: H
        species(4) = &
            +2.50000000d+00 &
            +0.00000000d+00 * tc(1) &
            +0.00000000d+00 * tc(2) &
            +0.00000000d+00 * tc(3) &
            +0.00000000d+00 * tc(4)
        !species 5: O
        species(5) = &
            +2.54205966d+00 &
            -2.75506191d-05 * tc(1) &
            -3.10280335d-09 * tc(2) &
            +4.55106742d-12 * tc(3) &
            -4.36805150d-16 * tc(4)
        !species 6: OH
        species(6) = &
            +2.86472886d+00 &
            +1.05650448d-03 * tc(1) &
            -2.59082758d-07 * tc(2) &
            +3.05218674d-11 * tc(3) &
            -1.33195876d-15 * tc(4)
        !species 7: HO2
        species(7) = &
            +4.01721090d+00 &
            +2.23982013d-03 * tc(1) &
            -6.33658150d-07 * tc(2) &
            +1.14246370d-10 * tc(3) &
            -1.07908535d-14 * tc(4)
        !species 8: H2O2
        species(8) = &
            +4.57316685d+00 &
            +4.33613639d-03 * tc(1) &
            -1.47468882d-06 * tc(2) &
            +2.34890357d-10 * tc(3) &
            -1.43165356d-14 * tc(4)
        !species 9: N2
        species(9) = &
            +2.92664000d+00 &
            +1.48797700d-03 * tc(1) &
            -5.68476100d-07 * tc(2) &
            +1.00970400d-10 * tc(3) &
            -6.75335100d-15 * tc(4)
    end if

end subroutine

! compute the h/(RT) at the given temperature (Eq 20)
! tc contains precomputed powers of T, tc(1) = log(T)
subroutine speciesEnthalpy(species, tc)

    double precision, intent(inout) :: species(9)
    double precision, intent(in) :: tc(5)
    ! temperature
    double precision :: T
    double precision :: invT

    T = tc(2)
    invT = 1.d0 / T

    ! species with midpoint at T=1000 kelvin
    if (T <  1.00000000d+03) then
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

end module fuego_module



