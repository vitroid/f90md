module settings
  use vector
  implicit none
  integer, parameter :: STDIN  = 5
  integer, parameter :: STDOUT = 6
  integer :: num_loop
  integer :: logging_interval
  real(kind=8) :: dt                            ! pico seconds

contains

  subroutine settings_read(FILE)
    use property
    use berendsen
    integer, intent(IN) :: FILE
    !local variables
    character(len=5) :: tag
    integer          :: i
    logging_interval = 10
    do
       read(FILE,*, END=99) tag
       if(tag == "@DTPS")then
          read(FILE,*) dt
       endif
       if(tag == "@MDLP")then
          read(FILE,*) num_loop
       endif
       if(tag == "@NLOG")then
          read(FILE,*) logging_interval
       endif
       if(tag == "@AR3A")then
          call property_read_ar3a(FILE)
       endif
       if(tag == "@BERE")then
          call berendsen_read_bere(FILE)
       endif
    enddo
99  continue
  end subroutine settings_read
  
  subroutine settings_write
  end subroutine settings_write
end module settings
