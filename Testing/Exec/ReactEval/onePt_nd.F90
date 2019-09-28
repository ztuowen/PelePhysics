module main_module

  use amrex_fort_module, only : amrex_real
  use reactor_module, only : enth_masT_ID, enth_masH_ID

  implicit none

  integer,parameter :: iE = enth_masH_ID
  integer,parameter :: ncells = 1

contains

    subroutine extern_init(name,namlen) bind(C, name="extern_init")

    use network
    use eos_module
    use transport_module

    implicit none
    integer :: namlen
    integer :: name(namlen)

    real (kind=amrex_real) :: small_temp = 1.d-200
    real (kind=amrex_real) :: small_dens = 1.d-200

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

    call network_init()

    eos_verbose = 0
    call eos_init(small_temp, small_dens)

    call transport_init()

  end subroutine extern_init

  subroutine extern_init_reactor() bind(C, name="extern_init_reactor")

#ifdef USE_SUNDIALS_PP
    use cvode_module, only : reactor_init
#else
    use reactor_module, only : reactor_init, reac_verbose
    reac_verbose = 1
#endif

    call reactor_init(iE,ncells)

  end subroutine extern_init_reactor


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module
    implicit none

    call transport_close()

  end subroutine extern_close


end module main_module
