program clustermd
  use settings
  use force
  use proceed
  use property
  use berendsen
  implicit none
  !local variables
  integer :: i
  call settings_read(STDIN)
  do i=1,num_loop
     ! velocity verlet
     call proceed_position(dt/2)
     call force_calculate()
     call force_to_accel()
     call proceed_velocity(dt)
     call property_kineticenergy()
     call berendsen_proceed(dt)
     call proceed_position(dt/2)
     if (logging_interval > 0 .and. mod(i,logging_interval) == 0) then
        write(STDOUT,*) i,ep,ek,ep+ek, temperature
     endif
  enddo
  call settings_write
end program clustermd
