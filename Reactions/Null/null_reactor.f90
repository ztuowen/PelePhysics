module reactor_module

  use amrex_fort_module, only : amrex_real
  use react_type_module

  implicit none

contains

  subroutine reactor_init()

    !nothing needed here

  end subroutine reactor_init


  subroutine reactor_close()

    ! nothing needed here

  end subroutine reactor_close


  function ok_to_react(state)

    implicit none

    type (react_t),intent(in) :: state
    logical                   :: ok_to_react

    ok_to_react = .true.

  end function ok_to_react


  function react_null(react_state_in, react_state_out, dt_react, time) result(stat)
    
    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: stat

    react_state_out = react_state_in
    stat % cost_value = 0.d0
    stat % reactions_succesful = .true.

  end function react_null


  function react(react_state_in, react_state_out, dt_react, time) result(stat)
    
    use eos_module

    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: stat

    stat = react_null(react_state_in, react_state_out, dt_react, time)

  end function react


end module reactor_module
