module settings
  use vector
  implicit none
  integer, parameter :: STDIN  = 5
  integer, parameter :: STDOUT = 6
  integer :: num_loop
  integer :: num_molecule
  real(kind=8) :: dt                            ! pico seconds
  real(kind=8) :: mass                          ! atomic unit
  real(kind=8) :: argon_sigma                   ! Angstrom
  real(kind=8) :: argon_eps4                    ! kJ/mol x 4

contains

  subroutine settings_read(FILE)
    use property
    integer, intent(IN) :: FILE
    !local variables
    character(len=5) :: tag
    integer          :: i
    do
       read(FILE,*, END=99) tag
       if(tag == "@DTPS")then
          read(FILE,*) dt
       endif
       if(tag == "@MDLP")then
          read(FILE,*) num_loop
       endif
       if(tag == "@AR3A")then
          read(FILE,*) num_molecule
          call property_prepare(num_molecule)
          call property_read_ar3a(FILE, num_molecule)
       endif
    enddo
99  continue
    mass = 39.95                         !atomic unit
    argon_sigma = 3.41d0                 !Angstrom
    argon_eps4 = 120d0 * 0.008314 * 4d0  !kJ/mol x 4
  end subroutine settings_read
  
  subroutine settings_write
  end subroutine settings_write
end module settings
