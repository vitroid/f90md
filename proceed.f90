module proceed
  implicit none
  
contains

  subroutine proceed_position(deltat)
    use settings
    use property
    real(kind=8), intent(IN) :: deltat    ! picoseconds
    !local variables
    integer                  :: i,k
    do i=1,num_molecule
       do k=1,3
          position(i)%vec(k) = position(i)%vec(k) + velocity(i)%vec(k) * deltat
       enddo
    enddo
  end subroutine proceed_position

  subroutine proceed_velocity(deltat)
    use settings
    use property
    real(kind=8), intent(IN) :: deltat    ! picoseconds
    !local variables
    integer                  :: i,k
    do i=1,num_molecule
       do k=1,3
          velocity(i)%vec(k) = velocity(i)%vec(k) + force_acc(i)%vec(k) * deltat
       enddo
    enddo
  end subroutine proceed_velocity

end module proceed
