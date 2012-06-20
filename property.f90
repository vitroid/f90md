module property
  use vector
  use physconst
  implicit none
  real(kind=8) :: ep,ek                         ! kJ/mol
  real(kind=8) :: temperature                   ! K
  integer :: num_molecule
  type(vector3), allocatable :: position(:)      ! Angstrom
  type(vector3), allocatable :: velocity(:)      ! Angstrom/ps
  type(vector3), allocatable :: force_acc(:)     ! Force or Acceleration
  real(kind=8) :: mass                          ! atomic unit
  real(kind=8) :: argon_sigma                   ! Angstrom
  real(kind=8) :: argon_eps4                    ! kJ/mol x 4

contains

  subroutine property_prepare
    allocate(position(num_molecule))
    allocate(velocity(num_molecule))
    allocate(force_acc(num_molecule))
    mass = 39.95                         !atomic unit
    argon_sigma = 3.41d0                 !Angstrom
    argon_eps4 = 120d0 * 1d-3 * kB * 4d0  !kJ/mol x 4
  end subroutine property_prepare

  subroutine property_read_ar3a(FILE)
    integer, intent(IN) :: FILE
    !local variables
    integer :: i
    read(FILE,*) num_molecule
    call property_prepare
    do i=1,num_molecule
       call vector3_read(FILE,position(i))
    enddo
  end subroutine property_read_ar3a

  subroutine property_kineticenergy
    !local variables
    integer :: i
    ek = 0d0
    do i=1,num_molecule
       ek = ek + 0.5 * mass * vector3_inner_product(velocity(i), velocity(i))
    enddo
    ! ek in g/mol * (A/ps)^2 == 10 J/mol
    ek = ek * 0.01 ! kJ/mol 3NkbT = 2Ek
    temperature = ek * 1d-3 * 2d0 /(3d0 * num_molecule * kB)
  end subroutine property_kineticenergy
end module property
