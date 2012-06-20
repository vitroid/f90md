module mpi_replica
  use mpi
  implicit none
  !MPI variables
  integer :: NPROCS, MYRANK
  integer, parameter :: ROOT=0
  !Parameters
  integer :: exchange_interval

contains

  subroutine mpi_replica_new
    !local variables
    integer :: ierr
    ! start parallel processes.
    call mpif_new( NPROCS, ierr )
    ! Obtain who am i (starting from 0)
    call mpif_comm_rank( MYRANK, ierr )
    ! Obtain number of processes
    call mpif_comm_size( NPROCS, ierr )
    ! interval is fixed at now
    exchange_interval = 1000   !md steps
  end subroutine mpi_replica_new
  
  subroutine mpi_replica_exchange(step)
    use berendsen
    use property
    integer, intent(IN) :: step
    !local variable
    integer :: ierr
    integer :: parity, myparity, theother
    logical :: imleader, havearest
    real(kind=8) :: mpi_ep(NPROCS)
    real(kind=8) :: mpi_temp0(NPROCS)
    if ( mod(step, exchange_interval) == 0 ) then
       !Let's try exchanges
       !To simplify the problem, trials are done on the root (0th) node.
       call mpi_gather(ep,1, MPI_REAL8, mpi_ep, 1, MPI_REAL8, ROOT, COMM, ierr)
       call mpi_gather(temp0,1, MPI_REAL8, mpi_temp0, 1, MPI_REAL8, ROOT, COMM, ierr)
       if ( MYRANK == ROOT ) then
          parity = step / exhange_interval
          parity = mod(parity, 2)
          do i=parity, NPROCS-2, 2
             !note that array in fortran starts from 1, not 0
             j = i + 1
             k = j + 1
             betaj = 1.0 / (kB * mpi_temp0(j))
             betak = 1.0 / (kB * mpi_temp0(k))
             before = exp(-betaj*mpi_ep(j))*exp(-betak*mpi_ep(k))  ! likelihoods
             after  = exp(-betaj*mpi_ep(k))*exp(-betak*mpi_ep(j))  ! likelihoods
             accept = .FALSE.
             if ( after > before ) then
                accept = .TRUE.
             else if ( rand < after / before ) then
                accept = .TRUE.
             endif
             if ( accept ) then
                tempj = mpi_temp0(j)
                tempk = mpi_temp0(k)
                mpi_temp0(j) = tempk
                mpi_temp0(k) = tempj
             endif
          enddo
       endif
       call mpi_scatter(mpi_temp0,1, MPI_REAL8, temp0, 1, MPI_REAL8, ROOT, COMM, ierr)
    endif
  end subroutine mpi_replica_exchange
end module mpi_replica
       
                
