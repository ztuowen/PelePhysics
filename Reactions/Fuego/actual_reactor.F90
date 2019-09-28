module reactor_module

  use amrex_fort_module, only : amrex_real
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use network, only: nspecies, spec_names
  use react_type_module
  use eos_type_module
  use amrex_error_module

  implicit none

  integer :: reac_verbose = 1
  
  integer, parameter :: eint_mass_ID = 1
  integer, parameter :: enth_pres_ID = 5
  integer, parameter :: enth_masT_ID = 100
  integer, parameter :: enth_masH_ID = 200

  character (len=*), parameter :: eint_mass_label = "Constant Internal Energy and Mass via (Y,e)"
  character (len=*), parameter :: enth_pres_label = "Constant Enthalpy and Pressure via (Y,e)"
  character (len=*), parameter :: enth_masT_label = "Constant Enthalpy and Mass via (rhoY,T)"
  character (len=*), parameter :: enth_masH_label = "Constant Enthalpy and Mass via (rhoY,rhoH)"

  integer,          parameter, private :: neq = nspecies + 1
  integer,          parameter, private :: itol = 1
  integer,          parameter, private :: order = 0
  integer,          parameter, private :: maxstep = 10000
  logical,          parameter, private :: use_ajac = .false.
  logical,          parameter, private :: save_ajac = .false.
  logical,          parameter, private :: stiff = .true.
  real(amrex_real), parameter, private :: rtol = 1.d-10
  real(amrex_real), parameter, private :: atol = 1.d-10
  real(amrex_real), parameter, private :: low_Y_thresh = -1.d-3;
  real(amrex_real), parameter, private :: low_T_thresh = 250.d0;

  type (eos_t) :: eos_state
  
  logical, save, private    :: reactor_initialized = .false.
  integer, private          :: iE, low_Y_test, low_T_test, strang_fix
  real(amrex_real), private :: vodeVec(neq),cdot(nspecies),rhoydot_ext(nspecies),ydot_ext(nspecies),Tcell
  real(amrex_real), private :: rhoedot_ext, rhoe_init, time_init, rhohdot_ext, rhoh_init, hdot_ext, h_init, pressureInit
  integer, private          :: viwsave(12)
  
  !$omp threadprivate(reactor_initialized,vodeVec,cdot,rhoydot_ext,ydot_ext,rhoedot_ext,rhoe_init,time_init,rhohdot_ext,Tcell)
  !$omp threadprivate(rhoh_init,hdot_ext,h_init,pressureInit,iE,low_Y_test,low_T_test,strang_fix,eos_state,viwsave)

contains

  subroutine reactor_init(iE_in, Ncells) bind(C, name="reactor_init")

    use, intrinsic :: iso_c_binding
    use vode_module, only : vode_init
    use extern_probin_module, only : new_Jacobian_each_cell
    use amrex_omp_module

    implicit none
    integer(c_int),  intent(in   ) :: iE_in
    integer(c_int),  intent(in   ), optional :: Ncells
    logical :: always_new_j_loc, isio

    if (new_Jacobian_each_cell .ne. 0) then
       always_new_j_loc = .true.
    else
       always_new_j_loc = .false.
    endif

    call vode_init(neq,reac_verbose,itol,rtol,atol,order,&
         maxstep,use_ajac,save_ajac,always_new_j_loc,stiff)

    iE = iE_in
    
    isio = parallel_IOProcessor() .and. omp_get_thread_num().eq.0
    if (isio .and. reac_verbose.ge.1) then
       print *,"Using good ol' dvode", Ncells
       print *,"--> DENSE solver without Analytical J"
       print *,"--> Always new J ? ",always_new_j_loc

       if (iE .eq. eint_mass_ID) then
          print *,"--> reactor is: ",eint_mass_label
       else if (iE .eq. enth_pres_ID) then
          print *,"--> reactor is: ",enth_pres_label
       else if (iE .eq. enth_masT_ID) then
          print *,"--> reactor is: ",enth_masT_label
       else if (iE .eq. enth_masH_ID) then
          print *,"--> reactor is: ",enth_masH_label
       else
          print *,"--> reactor specified:",iE
          call amrex_error('unknown reactor type')
       endif
    end if 

    call build(eos_state)

    call XSETF(MAX(reac_verbose-1,0)) ! suppress DVODE messages unless reac_verbose >= 2

    reactor_initialized = .true.

  end subroutine reactor_init


  function react(rY_in,rY_src_in,rX_in,rX_src_in,P_in,dt_react,time,i,j,k) bind(C, name="react") result(cost_value)
    
    use vode_module, only : itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    use eos_module

    implicit none
    real(amrex_real),   intent(inout) :: rY_in(nspecies+1),rY_src_in(nspecies)
    real(amrex_real),   intent(inout) :: rX_in,rX_src_in,P_in
    real(amrex_real),   intent(inout) :: dt_react, time
    integer,            intent(inout) :: i, j, k

    integer :: cost_value
    
    ! For compatibility to remove later
    type(react_t) :: react_state_in

    external dvode

    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail, n
    real(amrex_real) :: vodeTime, vodeEndTime, rho, rhoInv, rhooldInv, sumRY_temp, rY_temp(nspecies)

    ! For compatibility to remove later
    call build(react_state_in)

    react_state_in %              T = rY_in(nspecies+1)
    react_state_in %        rhoY(:) = rY_in(1:nspecies)
    react_state_in %            rho = sum(react_state_in % rhoY(:))
    react_state_in % rhoYdot_ext(:) = rY_src_in(1:nspecies)

    if (iE .eq. eint_mass_ID) then
        react_state_in %              e = rX_in
        react_state_in %    rhoedot_ext = rX_src_in
    else if (iE .eq. enth_pres_ID) then
        react_state_in %              p = P_in
        react_state_in %              h = rX_in
        react_state_in %    rhohdot_ext = rX_src_in
    else ! enth_masT_ID, enth_masH_ID
        react_state_in %              h = rX_in / react_state_in % rho
        react_state_in %    rhohdot_ext = rX_src_in
    end if
    ! END For compatibility to remove later

    if (.not. reactor_initialized) then
       call amrex_error('reactor::react called before initialized')
    endif

    if ( .not. ok_to_react(react_state_in) ) then
       call amrex_error('reactor::react Not Ok To React')
       return
    end if

    eos_state % rho               = sum(react_state_in % rhoY(:))
    eos_state % T                 = react_state_in % T
    rhoInv                         = 1.d0 / eos_state % rho
    eos_state % massfrac(1:nspecies) = react_state_in % rhoY(1:nspecies) * rhoInv

    if (iE .eq. eint_mass_ID) then
        eos_state % e = react_state_in % e
        call eos_re(eos_state)
    else if (iE .eq. enth_pres_ID) then
        pressureInit  = react_state_in % p
        eos_state % p = react_state_in % p
        eos_state % h = react_state_in % h
        call eos_ph(eos_state)
    else ! enth_masT_ID, enth_masH_ID
        eos_state % h = react_state_in % h
        call eos_rh(eos_state)
    end if

    if (always_new_j) call setfirst(.true.)

    ! Fill initial state and sources
    if (iE .eq. eint_mass_ID) then
        rhoe_init               = eos_state % e  *  eos_state % rho
        rhoedot_ext             = react_state_in % rhoedot_ext
        rhoydot_ext(1:nspecies) = react_state_in % rhoydot_ext(1:nspecies)
        vodeVec(1:nspecies)     = react_state_in % rhoY(:)
        vodeVec(neq)            = eos_state % T
    else if (iE .eq. enth_pres_ID) then
        h_init                  = eos_state % h  
        hdot_ext                = react_state_in % rhohdot_ext / eos_state % rho
        ydot_ext(1:nspecies)    = react_state_in % rhoydot_ext(1:nspecies) / eos_state % rho
        vodeVec(1:nspecies)     = react_state_in % rhoY(:) / eos_state % rho
        vodeVec(neq)            = eos_state % T
    else ! enth_masT_ID, enth_masH_ID
        rhoh_init               = eos_state % h  *  eos_state % rho
        rhohdot_ext             = react_state_in % rhohdot_ext 
        rhoydot_ext(1:nspecies) = react_state_in % rhoydot_ext(1:nspecies)
        vodeVec(1:nspecies)     = react_state_in % rhoY(:)
        if (iE .eq. enth_masH_ID) then
           vodeVec(neq)         = rhoh_init
           Tcell                = react_state_in % T
        else
           vodeVec(neq) = eos_state % T
        end if
    end if

    time_init   = time ! used inside rhs to integrate rhoh
    MF          = vode_MF
    vodeTime    = time
    vodeEndTime = time + dt_react

    ! Vode: istate
    ! in:  1: init, 2: continue (no change), 3: continue (w/changes)
    ! out: 1: nothing was done, 2: success
    !      -1: excessive work, -2: too much accuracy, -3: bad input
    !      -4: repeated step failure, -5: repeated conv failure
    !      -6: EWT became 0
    !      -7: state out of bounds
    istate = 1
    low_Y_test = 0
    low_T_test = 0
    strang_fix = 0

    call dvode(f_rhs, neq, vodeVec(:), vodeTime, vodeEndTime,&
         itol, rtol, atol, itask, istate, iopt, voderwork, lvoderwork, &
         vodeiwork, lvodeiwork, f_jac, MF, voderpar, vodeipar)

#ifdef MOD_REACTOR
    time = vodeTime
#endif

    if (istate .gt. 0 .and. reac_verbose .ge. 3) then
       write(6,*) '......dvode done at:', i, j, k
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '    number of fs(RHS EVAL) = ',vodeiwork(12)
       write(6,*) '              number of Js = ',vodeiwork(13)
       write(6,*) '    method order last used = ',vodeiwork(14)
       write(6,*) '   method order to be used = ',vodeiwork(15)
       write(6,*) '            number of LUDs = ',vodeiwork(19)
       write(6,*) ' number of Newton iterations ',vodeiwork(20)
       write(6,*) ' number of Newton failures = ',vodeiwork(21)
       if (istate.eq.-4 .or. istate.eq.-5) then
          ifail = vodeiwork(16)
          if (ifail .eq. nspecies+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
       end if
    end if

    cost_value = DBLE(vodeiwork(12)) ! number of f evaluations
    
    if (istate.le.0 .or. low_Y_test.gt.0 .or. low_T_test.gt.0) then

       do n=1,12
          viwsave(n) = vodeiwork(10+n)
       enddo

       strang_fix = 1

       ! Add explit update from F into in states
       do n=1,nspecies
          rY_in(n) = rY_in(n) + dt_react * rY_src_in(n)
       enddo
       rX_in = rX_in + dt_react * rX_src_in

       rho = sum(rY_in(1:nspecies))
       rhoInv = 1.d0 / rho

       ! Make temporary version of above, but guaranteed non-negative
       sumRY_temp = 0.d0
       do n=1,nspecies
          rY_temp(n) = MAX(rY_in(n),0.d0)
          sumRY_temp = sumRY_temp + rY_in(n) * rhoInv
       enddo
       
       ! Reload initial state and (zeroed) sources
       if (iE .eq. eint_mass_ID) then
          rhoedot_ext             = 0.d0
          rhoydot_ext(1:nspecies) = 0.d0
          vodeVec(1:nspecies)     = rY_temp(1:nspecies)
          vodeVec(neq) = eos_state % T
       else if (iE .eq. enth_pres_ID) then
          hdot_ext                = 0.d0
          ydot_ext(1:nspecies)    = 0.d0
          vodeVec(1:nspecies)     = rY_temp(1:nspecies) / sumRY_temp ! Here, Y>=0 and Sum(Y)=1
          vodeVec(neq) = eos_state % T
       else ! enth_masT_ID, enth_masH_ID
          rhohdot_ext             = 0.d0
          rhoydot_ext(1:nspecies) = 0.d0
          vodeVec(1:nspecies)     = rY_temp(1:nspecies)
           if (iE .eq. enth_masH_ID) then
              vodeVec(neq)        = rhoh_init
           else
              vodeVec(neq) = eos_state % T
           end if
       end if

       vodeTime    = time
       vodeEndTime = time + dt_react
       time_init = time

       istate = 1
       low_Y_test = 0
       low_T_test = 0
       strang_fix = 0
       call dvode(f_rhs, neq, vodeVec(:), vodeTime, vodeEndTime,&
            itol, rtol, atol, itask, istate, iopt, voderwork, lvoderwork, &
            vodeiwork, lvodeiwork, f_jac, MF, voderpar, vodeipar)

       if (istate.le.0) then
          print *,'vode failed again at ',i,j,k
          print *,'input state:'
          print *,'T',react_state_in%T
          if (iE .eq. 1) then
             print *,'e',react_state_in%e
          else
             print *,'h',react_state_in%h
          end if
          print *,'rho',eos_state%rho
          print *,'rhoY',react_state_in%rhoY
          if (iE .eq. 1) then
             print *,'rhoe forcing',react_state_in%rhoedot_ext
          else
             print *,'rhoh forcing',react_state_in%rhohdot_ext
          end if
          print *,'rhoY forcing',react_state_in%rhoydot_ext(1:nspecies)

          write(6,*) '......dvode data:'
          write(6,*) ' last successful step size = ',voderwork(11)
          write(6,*) '          next step to try = ',voderwork(12)
          write(6,*) '   integrated time reached = ',voderwork(13)
          write(6,*) '      number of time steps = ',vodeiwork(11)+viwsave(1)
          write(6,*) '              number of fs = ',vodeiwork(12)+viwsave(2)
          write(6,*) '              number of Js = ',vodeiwork(13)+viwsave(3)
          write(6,*) '    method order last used = ',vodeiwork(14)
          write(6,*) '   method order to be used = ',vodeiwork(15)
          write(6,*) '            number of LUDs = ',vodeiwork(19)+viwsave(9)
          write(6,*) ' number of Newton iterations ',vodeiwork(20)+viwsave(10)
          write(6,*) ' number of Newton failures = ',vodeiwork(21)+viwsave(11)
          if (istate.eq.-4 .or. istate.eq.-5) then
             ifail = vodeiwork(16)
             if (ifail .eq. nspecies+1) then
                write(6,*) '   T has the largest error'
             else
                write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
             end if
          end if

          print *,'Final T',vodeVec(neq)
          print *,'Final rhoY',vodeVec(1:nspecies)

          call amrex_error('vode failed')
       else
          if (reac_verbose .ge. 2) then
             write(*,*) '--> strang_fix recovery succesful at ',i,j,k

             write(6,*) '......dvode done at:', i, j, k
             write(6,*) ' last successful step size = ',voderwork(11)
             write(6,*) '          next step to try = ',voderwork(12)
             write(6,*) '   integrated time reached = ',voderwork(13)
             write(6,*) '          total time steps = ',vodeiwork(11)+viwsave(1) 
             write(6,*) '             total f evals = ',vodeiwork(12)+viwsave(2) 
             write(6,*) '             total J evals = ',vodeiwork(13)+viwsave(3) 
             write(6,*) '    method order last used = ',vodeiwork(14)
             write(6,*) '   method order to be used = ',vodeiwork(15)
             write(6,*) '                total LUDs = ',vodeiwork(19)+viwsave(9) 
             write(6,*) '        total Newton iters = ',vodeiwork(20)+viwsave(10)
             write(6,*) '     total Newton failures = ',vodeiwork(21)+viwsave(11)

          endif

       endif

       cost_value = cost_value + DBLE(vodeiwork(12)) ! number of f evaluations

    endif
    
    if (iE .eq. eint_mass_ID) then
       eos_state % rho               = sum(vodeVec(1:nspecies))
       rhoInv                        = 1.d0 / eos_state % rho
       eos_state % massfrac(1:nspecies) = vodeVec(1:nspecies) * rhoInv
       eos_state % T                 = vodeVec(neq)
       eos_state % e                 = (rhoe_init  +  dt_react*rhoedot_ext) /eos_state % rho
       call eos_re(eos_state) ! computes T, e, p
       rY_in(1:nspecies)             = vodeVec(1:nspecies)
       rX_in                         = eos_state % e
    else if (iE .eq. enth_pres_ID) then
       eos_state % p                 = pressureInit  
       eos_state % massfrac(1:nspecies) = vodeVec(1:nspecies)
       eos_state % T                 = vodeVec(neq)
       eos_state % h                 = (h_init  +  dt_react*hdot_ext)
       call eos_ph(eos_state) ! computes T, rho
       rY_in(1:nspecies)             = vodeVec(1:nspecies) * eos_state % rho
       rX_in                         = eos_state % h
    else ! enth_masT_ID, enth_masH_ID
       eos_state % rho               = sum(vodeVec(1:nspecies))
       rhoInv                        = 1.d0 / eos_state % rho
       if (strang_fix .eq. 1) then
          vodeVec(1:nspecies) = vodeVec(1:nspecies) + rY_in(n) - rY_temp(n)
       endif
       eos_state % massfrac(1:nspecies) = vodeVec(1:nspecies) * rhoInv
       eos_state % h                 = (rhoh_init  +  dt_react*rhohdot_ext) * rhoInv
       if (iE .eq. enth_masH_ID) then
          eos_state % T              = Tcell
       else
          eos_state % T              = vodeVec(neq)
       end if
       call eos_rh(eos_state) ! computes T, p
       rY_in(1:nspecies)             = vodeVec(1:nspecies)
       rX_in                         = eos_state % h * eos_state % rho
    end if

    rY_in(1+nspecies)      = eos_state % T

  end function react


  subroutine f_rhs(neq_in, time, y, ydot, rpar, ipar, ierr)

    use chemistry_module, only : molecular_weight
    use eos_module

    implicit none
    integer,         intent(in)   :: neq_in, ipar(*)
    real(amrex_real), intent(in)  :: y(neq), time, rpar(*)
    real(amrex_real), intent(out) :: ydot(neq)
    integer,         intent(out)  :: ierr
    integer          :: n
    real(amrex_real) :: rhoInv

    if (neq .ne. neq_in) call amrex_error('f_rhs: neq trashed somehow')
    ierr = 0
    if (iE .eq. eint_mass_ID) then
        eos_state % rho = sum(y(1:nspecies))
        rhoInv = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspecies) = y(1:nspecies) * rhoInv
        do n=1,nspecies
           if (eos_state%massfrac(n) .lt. low_Y_thresh) then
              low_Y_test = 1
              ierr = 1
              return
           endif
        enddo
        eos_state % T = y(neq) ! guess
        eos_state % e = (rhoe_init + (time - time_init) * rhoedot_ext) * rhoInv
        call eos_re(eos_state)
        call eos_get_activity(eos_state)
    else if (iE .eq. enth_pres_ID) then
        eos_state % massfrac(1:nspecies) = y(1:nspecies) 
        do n=1,nspecies
           if (eos_state%massfrac(n) .lt. low_Y_thresh) then
              low_Y_test = 1
              ierr = 1
              return
           endif
        enddo
        eos_state % T = y(neq) ! guess
        eos_state % h = h_init + (time - time_init) * hdot_ext
        eos_state % p = pressureInit  
        call eos_ph(eos_state)
        call eos_get_activity_h(eos_state)
    else ! enth_masT_ID, enth_masH_ID
        eos_state % rho = sum(y(1:nspecies))
        rhoInv = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspecies) = y(1:nspecies) * rhoInv 
        do n=1,nspecies
           if (eos_state%massfrac(n) .lt. low_Y_thresh) then
              low_Y_test = 1
              ierr = 1
              return
           endif
        enddo
        if (iE .eq. enth_masH_ID) then
           eos_state % T = Tcell ! guess
        else
           eos_state % T = y(neq) ! guess
        endif
        eos_state % h = (rhoh_init + (time - time_init) * rhohdot_ext) * rhoInv
        call eos_rh(eos_state)
        call eos_get_activity_h(eos_state)
    end if
    if (eos_state%T .lt. low_T_thresh) then
       low_T_test = 1
       ierr = 2
       return
    endif

    call ckwc(eos_state % T, eos_state % Acti, cdot)

    if (iE .eq. eint_mass_ID) then
        ydot(neq)    = rhoedot_ext 
        do n=1,nspecies
           ydot(n)   = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%ei(n)*ydot(n)
        end do
        ydot(neq)    = ydot(neq)/(eos_state%rho * eos_state%cv)
    else if (iE .eq. enth_pres_ID) then
        ydot(neq)    = hdot_ext
        do n=1,nspecies
           ydot(n)   = cdot(n) * molecular_weight(n) / eos_state%rho + ydot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%hi(n)*ydot(n) 
        end do
        ydot(neq)    = ydot(neq)/eos_state%cp
    else ! enth_masT_ID, enth_masH_ID
        ydot(neq)    = rhohdot_ext
        do n=1,nspecies
           ydot(n)   = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%hi(n)*ydot(n) 
        end do
        if (iE .eq. enth_masH_ID) then
           ydot(neq) = rhohdot_ext
        else
           ydot(neq) = ydot(neq)/(eos_state%rho * eos_state%cp)
        end if
    end if

  end subroutine f_rhs


  subroutine f_jac(neq_in, npt, y, t, pd)

    implicit none
    integer,        intent(in)  :: neq_in, npt
    real(amrex_real),intent(in)  :: y(neq,npt), t
    real(amrex_real),intent(inout) :: pd(neq,neq)

    if (neq .ne. neq_in) call amrex_error('f_jac: neq trashed somehow')
    
    call amrex_error('DVODE version: Analytic Jacobian not yet implemented')

  end subroutine f_jac


  subroutine reactor_close() bind(C, name="reactor_close")

    implicit none
    call destroy(eos_state)

    reactor_initialized = .false.
   
  end subroutine reactor_close


  ! TODO: Remove this silly function
  function ok_to_react(state)

    use extern_probin_module, only: react_T_min, react_T_max, react_rho_min, react_rho_max

    implicit none

    type (react_t),intent(in) :: state
    logical                   :: ok_to_react
    real(amrex_real)           :: rho

    ok_to_react = .true.

    rho = sum(state % rhoY)
    if (state % T   .lt. react_T_min   .or. state % T   .gt. react_T_max .or. &
        rho         .lt. react_rho_min .or. rho         .gt. react_rho_max) then

       ok_to_react = .false.

    endif

  end function ok_to_react

end module reactor_module
