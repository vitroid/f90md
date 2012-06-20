program clustermd
  use settings
  use force
  use proceed
  use property
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
     write(STDOUT,*) i,ep,ek,ep+ek
     call proceed_position(dt/2)
  enddo
  call settings_write
end program clustermd
