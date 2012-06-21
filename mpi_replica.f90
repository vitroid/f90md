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
    !call mpi_new( NPROCS, ierr )
    call mpi_init(ierr)
    ! Obtain who am i (starting from 0)
    call mpi_comm_rank( MPI_COMM_WORLD, MYRANK, ierr )
    ! Obtain number of processes
    call mpi_comm_size( MPI_COMM_WORLD, NPROCS, ierr )
    ! interval is fixed at now
    exchange_interval = 1000   !md steps
  end subroutine mpi_replica_new
  
  subroutine mpi_replica_exchange(step)
    use berendsen
    use property
    integer, intent(IN) :: step
    !local variable
    integer i,j,k
    integer :: ierr
    integer :: parity, myparity, theother
    real(kind=8) :: betaj, betak, before, after, tempj, tempk, r
    logical :: imleader, havearest, accept
    real(kind=8) :: mpi_ep(NPROCS)
    real(kind=8) :: mpi_temp0(NPROCS)
    if ( mod(step, exchange_interval) == 0 ) then
       !Let's try exchanges
       !To simplify the problem, trials are done on the root (0th) node.
       call mpi_gather(ep,1, MPI_REAL8, mpi_ep, 1, MPI_REAL8, ROOT, MPI_COMM_WORLD, ierr)
       call mpi_gather(temp0,1, MPI_REAL8, mpi_temp0, 1, MPI_REAL8, ROOT, MPI_COMM_WORLD, ierr)
       if ( MYRANK == ROOT ) then
          parity = step / exchange_interval
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
             call random_number(r)
             if ( after > before ) then
                accept = .TRUE.
             else if ( r < after / before ) then
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
       call mpi_scatter(mpi_temp0,1, MPI_REAL8, temp0, 1, MPI_REAL8, ROOT, MPI_COMM_WORLD, ierr)
    endif
  end subroutine mpi_replica_exchange
end module mpi_replica
       
                
