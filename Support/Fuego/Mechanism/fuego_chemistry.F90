module chemistry_module

#include "mechanism.h"

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, parameter :: naux = 0   ! number of auxiliary components
  integer, parameter :: nspecies = NUM_SPECIES
  integer, parameter :: nreactions = NUM_REACTIONS
  integer, parameter :: nelements = NUM_ELEMENTS

  logical, save :: chemistry_initialized = .false.

  integer, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), allocatable, save :: elem_names(:)

  integer, parameter :: L_spec_name = 16 ! Each species name has at most 16 characters
  character*(L_spec_name), allocatable, save :: spec_names(:)

  integer, parameter :: L_aux_name = 16 ! Each aux name has at most 16 characters
  character*(L_aux_name), allocatable, save :: aux_names(:)

#ifdef AMREX_USE_CUDA
  real(amrex_real), allocatable, managed , save :: molecular_weight(:), inv_mwt(:)
#else
  real(amrex_real), allocatable, save :: molecular_weight(:), inv_mwt(:)
#endif


  real(amrex_real), save :: Ru, Ruc, Patm, rwrk
  integer, save          :: iwrk

  integer, private :: names(nspecies * L_spec_name)
  
contains

  subroutine chemistry_init()
    integer :: nfit, i, ic, ii
    real(amrex_real) :: T0
    integer, allocatable :: names(:)

    call ckinit()
    call ckindx(iwrk, rwrk, nelements, nspecies, nreactions, nfit)

    allocate(aux_names(naux))
    allocate(elem_names(nelements))
    allocate(spec_names(nspecies))
    allocate(molecular_weight(nspecies))
    allocate(inv_mwt(nspecies))

    allocate(names(nspecies*L_spec_name))  

    call cksyme(names, L_elem_name) 

    ic = 1
    do i = 1, nelements
       do ii=1, L_elem_name
          elem_names(i)(ii:ii) = char(names(ic))
          ic = ic + 1
       end do
    end do

    call cksyms(names, L_spec_name) 

    ic = 1
    do i = 1, nspecies
       do ii=1, L_spec_name
          spec_names(i)(ii:ii) = char(names(ic))
          ic = ic+1
       end do
    end do

    deallocate(names)

    call ckwt(iwrk, rwrk, molecular_weight)
    inv_mwt = 1.d0 / molecular_weight

    call ckrp(iwrk, rwrk, Ru, Ruc, Patm)

    chemistry_initialized = .true.

  end subroutine chemistry_init


  subroutine chemistry_close()
    call ckfinalize()
    deallocate(elem_names,spec_names,aux_names,molecular_weight,inv_mwt)
  end subroutine chemistry_close


  function get_species_index(name) result (iname)
    character(len=*), intent(in) :: name
    integer :: iname
    integer :: i
    iname = -1
    do i = 1, nspecies
       if (trim(spec_names(i)) .eq. trim(name)) then
          iname = i
          exit
       end if
    end do
  end function get_species_index


end module chemistry_module
