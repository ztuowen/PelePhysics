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



double precision :: fwd_A(21), fwd_beta(21), fwd_Ea(21)
double precision :: low_A(21), low_beta(21), low_Ea(21)
double precision :: troe_a(21),troe_Ts(21), troe_Tss(21), troe_Tsss(21)
double precision :: activation_units(21), prefactor_units(21), phase_units(21)
integer :: is_PD(21), troe_len(21), sri_len(21), nTB(21), TBid(21)
double precision :: TB(21,21)
contains

! A few mechanism parameters
subroutine ckindx(iwrk, rwrk, mm, kk, ii, nfit)

    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
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

    integer, intent(out) :: kname(plenkname*3)
    integer, intent(in) :: plenkname

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

    integer, intent(out) :: kname(plenkname*9)
    integer, intent(in) :: plenkname

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
subroutine ckrp(ickwrk, rckwrk, ru, ruc, pa)

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

    double precision, intent(in) :: rho
    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: P

    double precision :: YOW ! for computing mean MW

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
subroutine ckrhoy(P, T, y, iwrk, rwrk, rho)

    double precision, intent(in) :: P
    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: rho

    double precision :: YOW, tmp(9)
    integer :: i

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
subroutine ckwt(iwrk, rwrk, wt)

    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: wt(9)

    call molecularWeight(wt)

end subroutine

! convert y[species] (mass fracs) to x[species] (mole fracs)
subroutine ckytx(y, iwrk, rwrk, x)

    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: x(9)

    double precision :: YOW, YOWINV
    double precision :: tmp(9)
    integer :: i

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
        YOW(i) = 1.d0 / YOW(i)
    end do

    do n=1, 9
        do i=1, np
            x(i,n) = x(i,n) * YOW(i)
        end do
    end do

end subroutine

! convert y[species] (mass fracs) to c[species] (molar conc)
subroutine ckytcr(rho, T, y, iwrk, rwrk, c)

    double precision, intent(in) :: rho
    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: c(9)

    integer :: i

    do i=1, 9
        c(i) = rho * y(i) * imw(i)
    end do

end subroutine

! convert x[species] (mole fracs) to y[species] (mass fracs)
subroutine ckxty(x, iwrk, rwrk, y)

    double precision, intent(in) :: x(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
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
subroutine ckcvms(T, iwrk, rwrk, cvms)

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
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
subroutine ckcpms(T, iwrk, rwrk, cpms)

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
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
subroutine ckums(T, iwrk, rwrk, ums)

    double precision, intent(in) :: T
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(inout) :: ums(9)

    double precision :: tT, tc(5)
    double precision :: RT
    integer :: i

    tT = T ! temporary temperature
    tc = (/ 0.d0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT /) ! temperature cache
    RT = 8.31451000d+07 * tT ! R*T

    call speciesInternalEnergy(ums, tc)

    do i=1, 9
        ums(i) = ums(i) * (RT * imw(i))
    end do

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
    RT = 8.31451000d+07 * tT ! R*T

    call speciesEnthalpy(hms, tc)

    do i=1, 9
        hms(i) = hms(i) * (RT * imw(i))
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
            hms(i,n) = hms(i,n) * ( 8.31451000d+07 * T(i) * imw(n))
        end do
    end do

end subroutine

! Returns the mean specific heat at CP (Eq. 34)
subroutine ckcpbs(T, y, iwrk, rwrk, cpbs)

    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: cpbs

    double precision :: cpor(9)
    double precision :: tresult(9)
    double precision :: tT, tc(5)
    double precision :: res
    integer :: i

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
subroutine ckcvbs(T, y, iwrk, rwrk, cvbs)

    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: cvbs

    double precision :: cvor(9)
    double precision :: tT, tc(5)
    double precision :: res

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
subroutine ckubms(T, y, iwrk, rwrk, ubms)

    double precision, intent(in) :: T
    double precision, intent(in) :: y(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
    double precision, intent(out) :: ubms

    double precision :: ums(9) ! temporary energy array
    double precision :: res
    double precision :: RT, tT, tc(5)

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
subroutine ckwc(T, C, iwrk, rwrk, wdot)

    double precision, intent(in) :: T
    double precision, intent(inout) :: C(9)
    integer, intent(in) :: iwrk
    double precision, intent(in) :: rwrk
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

    double precision, intent(inout) :: wdot(9)
    double precision, intent(in) :: sc(9)
    double precision, intent(in) :: T

    double precision :: tc(5)
    double precision :: invT
    double precision :: T_save
    double precision :: k_f_save(21)
    double precision :: Kc_save(21)
    double precision :: qdot, q_f(21), q_r(21)
    integer :: i

    tc = (/ log(T), T, T*T, T*T*T, T*T*T*T /) ! temperature cache
    invT = 1.d0 / tc(2)
    T_save = -1.d0

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

    double precision, intent(in) :: tc(5)
    double precision, intent(in) :: invT
    double precision, intent(out) :: k_f(21)

    integer :: i

    do i=1, 21
        k_f(i) = prefactor_units(i)*fwd_A(i)*exp(fwd_beta(i)*tc(1)-activation_units(i)*fwd_Ea(i)*invT)
    end do

end subroutine

subroutine comp_Kc(tc, invT, Kc)

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
    double precision :: X, F_sri
    double precision :: k_f_save(21)
    double precision :: tmp1, tmp2, tmp3
    integer :: i

    ! reaction 1: H + O2 (+M) <=> HO2 (+M)
    qf(1) = sc(1+1)*sc(3+1)
    qr(1) = sc(6+1)
    ! reaction 2: H2O2 (+M) <=> OH + OH (+M)
    qf(2) = sc(7+1)
    qr(2) = sc(5+1)*sc(5+1)
    ! reaction 3: H2 + M <=> H + H + M
    qf(3) = sc(0+1)
    qr(3) = sc(3+1)*sc(3+1)
    ! reaction 4: O + O + M <=> O2 + M
    qf(4) = sc(4+1)*sc(4+1)
    qr(4) = sc(1+1)
    ! reaction 5: O + H + M <=> OH + M
    qf(5) = sc(3+1)*sc(4+1)
    qr(5) = sc(5+1)
    ! reaction 6: H + OH + M <=> H2O + M
    qf(6) = sc(3+1)*sc(5+1)
    qr(6) = sc(2+1)
    ! reaction 7: H + O2 <=> O + OH
    qf(7) = sc(1+1)*sc(3+1)
    qr(7) = sc(4+1)*sc(5+1)
    ! reaction 8: O + H2 <=> H + OH
    qf(8) = sc(0+1)*sc(4+1)
    qr(8) = sc(3+1)*sc(5+1)
    ! reaction 9: H2 + OH <=> H2O + H
    qf(9) = sc(0+1)*sc(5+1)
    qr(9) = sc(2+1)*sc(3+1)
    ! reaction 10: O + H2O <=> OH + OH
    qf(10) = sc(2+1)*sc(4+1)
    qr(10) = sc(5+1)*sc(5+1)
    ! reaction 11: HO2 + H <=> H2 + O2
    qf(11) = sc(3+1)*sc(6+1)
    qr(11) = sc(0+1)*sc(1+1)
    ! reaction 12: HO2 + H <=> OH + OH
    qf(12) = sc(3+1)*sc(6+1)
    qr(12) = sc(5+1)*sc(5+1)
    ! reaction 13: HO2 + O <=> O2 + OH
    qf(13) = sc(4+1)*sc(6+1)
    qr(13) = sc(1+1)*sc(5+1)
    ! reaction 14: HO2 + OH <=> H2O + O2
    qf(14) = sc(5+1)*sc(6+1)
    qr(14) = sc(1+1)*sc(2+1)
    ! reaction 15: HO2 + HO2 <=> H2O2 + O2
    qf(15) = sc(6+1)*sc(6+1)
    qr(15) = sc(1+1)*sc(7+1)
    ! reaction 16: HO2 + HO2 <=> H2O2 + O2
    qf(16) = sc(6+1)*sc(6+1)
    qr(16) = sc(1+1)*sc(7+1)
    ! reaction 17: H2O2 + H <=> H2O + OH
    qf(17) = sc(3+1)*sc(7+1)
    qr(17) = sc(2+1)*sc(5+1)
    ! reaction 18: H2O2 + H <=> HO2 + H2
    qf(18) = sc(3+1)*sc(7+1)
    qr(18) = sc(0+1)*sc(6+1)
    ! reaction 19: H2O2 + O <=> OH + HO2
    qf(19) = sc(4+1)*sc(7+1)
    qr(19) = sc(5+1)*sc(6+1)
    ! reaction 20: H2O2 + OH <=> HO2 + H2O
    qf(20) = sc(5+1)*sc(7+1)
    qr(20) = sc(2+1)*sc(6+1)
    ! reaction 21: H2O2 + OH <=> HO2 + H2O
    qf(21) = sc(5+1)*sc(7+1)
    qr(21) = sc(2+1)*sc(6+1)

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
    alpha_troe(1) = mixture + (TB(1,1) - 1)*sc(0+1) + (TB(1,2) - 1)*sc(2+1) + (TB(1,3) - 1)*sc(1+1)
    alpha_troe(2) = mixture + (TB(2,1) - 1)*sc(0+1) + (TB(2,2) - 1)*sc(2+1)
    do i=1, 2

        redP = alpha_troe(i-0) / k_f_save(i) * phase_units(i) * low_A(i) * exp(low_beta(i) * tc(1) - activation_units(i) * low_Ea(i) *invT)
        F = redP / (1.d0 + redP)
        logPred = log10(redP)

        if (troe_Tsss(i) > 1.d-100) then
           tmp1 = abs((1.d0-troe_a(i))*exp(-T/troe_Tsss(i)))
        else
           tmp1 = 0.d0
        end if

        if (troe_Ts(i) > 1.d-100) then
           tmp2 = abs(troe_a(i) * exp(-T/troe_Ts(i)))
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
        F_troe = (logFcent / (1.d0 + troe*troe)) ** 10.d0
        Corr(i) = F * F_troe
    end do

end subroutine

! compute the g/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine gibbs(species, tc)

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

    double precision, intent(out) :: species(9)
    double precision, intent(in) :: tc(5)

    double precision :: T

    T = tc(2)

    ! species with midpoint at T=1000 kelvin
    if (T < 1000.000000d0) then
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

! compute the e/(RT) at the given temperature
! tc contains precomputed powers of T, tc[0] = log(T)
subroutine speciesInternalEnergy(species, tc)

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

! save molecular weights into array
subroutine molecularWeight(wt)

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
subroutine get_t_given_ey(e, y, iwrk, rwrk, t, ierr)

    double precision, intent(in) :: e
    double precision, intent(in) :: y(9)
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
        t1 = t1 + dt
        end if
    end do

    t = t1
    ierr = 0

end subroutine

end module fuego_module



