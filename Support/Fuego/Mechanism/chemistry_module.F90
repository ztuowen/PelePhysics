module chemistry_module

#include "mechanism.h"
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, parameter :: naux = 0   ! number of auxiliary components
  integer, parameter :: nspecies = NSP   ! number of species
  integer, parameter :: nreactions = FWD_RATES ! number of reactions

  logical, save :: chemistry_initialized = .false.

  integer, parameter :: L_spec_name = 16 ! Each species name has at most 16 characters
  character*(L_spec_name), save :: spec_names(nspecies)

  integer, parameter :: L_aux_name = 16 ! Each aux name has at most 16 characters
  character*(L_aux_name), save :: aux_names(naux)

  real(amrex_real), save :: molecular_weight(nspecies), inv_mwt(nspecies)

  real(amrex_real), save :: Ru, Ruc, Patm, rwrk
  integer, save          :: iwrk

contains

  subroutine chemistry_init()
    integer :: nfit, i, ic, ii
    real(amrex_real) :: T0
    integer, allocatable :: names(:)

    call ckinit()

    allocate(names(nspecies*L_spec_name))  
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
