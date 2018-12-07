module chemistry_module

  use amrex_fort_module, only : amrex_real
  implicit none

  integer, parameter :: naux = 0       ! number of auxiliary components
  integer, parameter :: nelements = 1  ! number of elements
  integer, parameter :: nspecies = 1   ! number of species
  integer, parameter :: nreactions = 0 ! number of reactions

  logical, save :: chemistry_initialized = .false.

  integer, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), save :: elem_names(nelements)

  integer, parameter :: L_spec_name = 8 ! Each species name has at most 8 characters
  character*(L_spec_name), save :: spec_names(nspecies)

  real(amrex_real), save :: molecular_weight(nspecies), inv_mwt(nspecies)

  real(amrex_real), parameter :: Ru=8.31451d+07, Ruc=1.98721558317399615845, Patm=1.01325d+06

contains

  subroutine chemistry_init()
    spec_names(1) = 'X'
    molecular_weight(1) = 28.d0
    inv_mwt(1) = 1.d0 / molecular_weight(1)
    chemistry_initialized = .true.
  end subroutine chemistry_init

  subroutine chemistry_close()
  end subroutine chemistry_close
  
end module chemistry_module
