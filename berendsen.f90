module berendsen
  implicit none
  real(kind=8) :: tau
  real(kind=8) :: temp0

contains

  subroutine berendsen_read_bere(FILE)
    integer, intent(IN) :: FILE
    read(FILE,*) temp0, tau
  end subroutine berendsen_read_bere

  ! must be called after property_kineticenergy
  subroutine berendsen_proceed(deltat)
    use property
    real(kind=8), intent(IN) :: deltat
    !local variables
    real(kind=8) :: diff, newtemp, ratio
    diff = temp0 - temperature  !must be calculated in advance
    diff = diff / tau
    newtemp = temperature + deltat * diff
    ratio = sqrt(newtemp / temperature)
    call property_scalevelocity(ratio)
  end subroutine berendsen_proceed
end module berendsen
  
