program clustermd
  use settings
  use force
  use proceed
  use property
  use berendsen
  use mpi_replica
  implicit none
  !local variables
  integer :: i
  !reading from command line does not work with MPI
  open(10,FILE="test.input")
  call settings_read(10)
  !All the subprocesses share the settings. Then fork.
  call mpi_replica_new
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
        if ( MYRANK == NPROCS-1 ) then
           write(STDOUT,fmt="(i12,1x,5(f12.7,1x),i4)") i,ep,ek,ep+ek, temp0, temperature, MYRANK
        endif
     endif
     call mpi_replica_exchange(i)
  enddo
  call settings_write
  call mpi_replica_done
end program clustermd
