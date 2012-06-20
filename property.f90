module property
  use vector
  implicit none
  real(kind=8) :: ep,ek                         ! kJ/mol
  type(vector3), allocatable :: position(:)      ! Angstrom
  type(vector3), allocatable :: velocity(:)      ! Angstrom/ps
  type(vector3), allocatable :: force_acc(:)     ! Force or Acceleration

contains

  subroutine property_prepare(num_molecule)
    integer, intent(IN) :: num_molecule
    allocate(position(num_molecule))
    allocate(velocity(num_molecule))
    allocate(force_acc(num_molecule))
  end subroutine property_prepare

  subroutine property_read_ar3a(FILE, num_molecule)
    integer, intent(IN) :: FILE, num_molecule
    !local variables
    integer :: i
    do i=1,num_molecule
       call vector3_read(FILE,position(i))
    enddo
  end subroutine property_read_ar3a

  subroutine property_kineticenergy
    use settings
    !local variables
    integer :: i
    ek = 0d0
    do i=1,num_molecule
       ek = ek + 0.5 * mass * vector3_inner_product(velocity(i), velocity(i))
    enddo
    ! ek in g/mol * (A/ps)^2 == 10 J/mol
    ek = ek * 0.01 ! kJ/mol
  end subroutine property_kineticenergy
end module property
