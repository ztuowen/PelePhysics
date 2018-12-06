module reactor_module

  use amrex_fort_module, only : amrex_real
  use chemistry_module, only: nspecies, spec_names
  use react_type_module
  use eos_type_module

  implicit none
  real(amrex_real), private :: vodeVec(nspecies+1),cdot(nspecies),rhoydot_ext(nspecies)
  real(amrex_real), private:: rhoedot_ext, rhoe_init, time_init
  integer,private :: iloc, jloc, kloc
  type (eos_t) :: eos_state
  !$omp threadprivate(vodeVec,cdot,rhoydot_ext,rhoedot_ext,rhoe_init,time_init,iloc,jloc,kloc,eos_state)

contains

  subroutine reactor_init()

    use vode_module, only : vode_init
    integer :: neq, verbose, itol, order, maxstep
    real(amrex_real) :: rtol, atol
    logical :: use_ajac, save_ajac, always_new_j, stiff

    neq = nspecies + 1
    verbose = 0
    itol = 1
    order = 2
    maxstep = 10000
    use_ajac = .false.
    save_ajac = .false.
    always_new_J = .false.
    stiff = .true.
    rtol = 1.d-10
    atol = 1.d-10
    call vode_init(neq,verbose,itol,rtol,atol,order,&
         maxstep,use_ajac,save_ajac,always_new_j,stiff)

    call build(eos_state)

  end subroutine reactor_init


  subroutine reactor_close()

    call destroy(eos_state)
   
  end subroutine reactor_close


  function react_null(react_state_in, react_state_out, dt_react, time)
    
    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: react_null

    react_state_out = react_state_in
    react_null % cost_value = 0.d0
    react_null % reactions_succesful = .true.

  end function react_null


  function react(react_state_in, react_state_out, dt_react, time)
    
    use amrex_error_module
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    use eos_module

    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: react

    external dvode

    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail, neq
    real(amrex_real) :: vodeTime, vodeEndTime, rhoInv

    eos_state % rho = sum(react_state_in % rhoY(:))
    eos_state % T = react_state_in % T
    rhoInv = 1.d0 / eos_state % rho
    eos_state % massfrac(1:nspecies) = react_state_in % rhoY(1:nspecies) * rhoInv
    eos_state % e = react_state_in % e

    call eos_re(eos_state)

    if (always_new_j) call setfirst(.true.)

    MF = vode_MF
    vodeTime = time
    vodeEndTime = time + dt_react
    neq = nspecies + 1

    vodeVec(1:nspecies) = react_state_in % rhoY(:)
    vodeVec(neq) = eos_state % T

    rhoe_init = eos_state % e  *  eos_state % rho
    rhoedot_ext = react_state_in % rhoedot_ext
    rhoydot_ext(1:nspecies) = react_state_in % rhoydot_ext(1:nspecies)
    time_init = time
    iloc = react_state_in % i
    jloc = react_state_in % j
    kloc = react_state_in % k

    ! Vode: istate
    ! in:  1: init, 2: continue (no change), 3: continue (w/changes)
    ! out: 1: nothing was done, 2: success
    !      -1: excessive work, -2: too much accuracy, -3: bad input
    !      -4: repeated step failure, -5: repeated conv failure
    !      -6: EWT became 0
    istate = 1

    call dvode(f_rhs, neq, vodeVec(:), vodeTime, vodeEndTime,&
         itol, rtol, atol, itask, istate, iopt, voderwork, lvoderwork, &
         vodeiwork, lvodeiwork, f_jac, MF, voderpar, vodeipar)

    if (verbose .ge. 1) then
       write(6,*) '......dvode done:'
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '              number of fs = ',vodeiwork(12)
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

    if (istate > 0) then

       react % reactions_succesful = .true.
       react % cost_value = DBLE(vodeiwork(12)) ! number of f evaluations

       eos_state % rho = sum(vodeVec(1:nspecies))
       eos_state % T = vodeVec(neq)
       rhoInv = 1.d0 / eos_state % rho
       eos_state % massfrac(1:nspecies) = vodeVec(1:nspecies) * rhoInv
       eos_state % e = (rhoe_init  +  dt_react*rhoedot_ext)/eos_state % rho
       call eos_re(eos_state)

       react_state_out % rhoY(:) = vodeVec(1:nspecies)
       react_state_out % rho = eos_state % rho
       react_state_out % T = eos_state % T
       react_state_out % e = eos_state % e

    else

       react % reactions_succesful = .false.

       print *,'vode failed at',react_state_in % i,react_state_in % j,react_state_in % k
       print *,'input state:'
       print *,'T',react_state_in%T
       print *,'e',react_state_in%e
       print *,'rho',react_state_in%rho
       print *,'rhoY',react_state_in%rhoY
       print *,'rhoe forcing',react_state_in%rhoedot_ext
       print *,'rhoY forcing',react_state_in%rhoydot_ext(1:nspecies)

       write(6,*) '......dvode data:'
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '              number of fs = ',vodeiwork(12)
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

       print *,'Final T',vodeVec(neq)
       print *,'Final rhoY',vodeVec(1:nspecies)

       call amrex_error('vode failed')

    end if

  end function react


  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)
    use chemistry_module, only : molecular_weight
    use eos_module
    integer,         intent(in)  :: neq, ipar(*)
    real(amrex_real), intent(in)  :: y(neq), time, rpar(*)
    real(amrex_real), intent(out) :: ydot(neq)
    integer         :: n
    real(amrex_real) :: rhoInv

    eos_state % rho = sum(y(1:nspecies))
    rhoInv = 1.d0 / eos_state % rho
    eos_state % massfrac(1:nspecies) = y(1:nspecies) * rhoInv
    eos_state % T = y(neq) ! guess
    eos_state % e = (rhoe_init + (time - time_init) * rhoedot_ext) * rhoInv
    call eos_re(eos_state)

    call eos_get_activity(eos_state)

    call ckwc(eos_state % T, eos_state % Acti, iwrk, rwrk, cdot)

    ydot(neq) = rhoedot_ext
    do n=1,nspecies
       ydot(n) = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
       ydot(neq) = ydot(neq) - eos_state%ei(n)*ydot(n)
    end do
    ydot(neq) = ydot(neq)/(eos_state%rho * eos_state%cv)
  end subroutine f_rhs


  subroutine f_jac(neq, npt, y, t, pd)
    use amrex_error_module

    integer,        intent(in)  :: neq, npt
    real(amrex_real),intent(in)  :: y(neq,npt), t
    real(amrex_real),intent(out) :: pd(neq,neq)

    call amrex_error('Analytic Jacobian not yet implemented')

  end subroutine f_jac

end module reactor_module
