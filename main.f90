program clustermd
  implicit none
  !molecular properties
  !Make arrays allocatable
  integer      :: num_molecule
  real(kind=8), allocatable :: position(:,:)  ! Angstrom
  real(kind=8), allocatable :: velocity(:,:)  ! Angstrom / ps
  real(kind=8), allocatable :: accel(:,:)     ! Angstrom / ps**2
  real(kind=8), allocatable :: force(:,:)     ! N
  real(kind=8) :: argon_mass                  ! atomic mass
  real(kind=8) :: argon_eps,argon_sig         ! kJ/mol, Angstrom
  real(kind=8) :: dt                          ! picoseconds
  real(kind=8) :: k
  real(kind=8) :: temperature
  !local variables
  integer :: num_loop
  integer :: i,j
  integer :: ix,iy,iz
  integer :: j1,j2
  real(kind=8) :: delta(3),dd
  real(kind=8) :: ep,ek
  !read number of molecules from file
  read(5,*) num_molecule
  allocate(position(3,num_molecule))
  allocate(velocity(3,num_molecule))
  allocate(accel(3,num_molecule))
  allocate(force(3,num_molecule))
  do i=1,num_molecule
     read(5,*) (position(j,i),j=1,3)   ! coordinates in Angstrom
  enddo
  velocity(:,1:num_molecule) = 0.0d0
  accel(:,1:num_molecule) = 0.0d0
  dt = 0.001

  !realistic interaction parameters for Argon.
  argon_mass = 39.95d0
  argon_eps = 120d0 * 0.008314d0
  argon_sig = 3.41d0
  num_loop = 1000000
  do i=1,num_loop
     !calculate force
     !force will be reset each step.
     force(:,1:num_molecule) = 0.0d0
     ep = 0.0
     !interaction between all pairs of molecules
     do j1 = 1,num_molecule
        do j2 = j1+1, num_molecule
           delta(:) = position(:,j2) - position(:,j1)
           dd = delta(1)**2 + delta(2)**2 + delta(3)**2
           ep = ep + 4*argon_eps*(argon_sig**12/dd**6 - argon_sig**6/dd**3)      ! kJ/mol
           k  =    - 4*argon_eps*(12*argon_sig**12/dd**7 - 6*argon_sig**6/dd**4) ! kJ/mol/Ang**2
           force(:,j1) = force(:,j1) + k * delta(:)    ! kJ/mol/Ang
           force(:,j2) = force(:,j2) - k * delta(:)
        enddo
     enddo
     !calculate accel
     !force is in kJ/mol/Angstrom
     !acc is in Angstrom/ps/ps
     accel(:,1:num_molecule) = force(:,1:num_molecule) / argon_mass * 100.0
                                               !  * 1d3      & ! J/mol/A
                                               !             & ! (J == kg m^2/s^2)
                                               !             & ! kg m^2/s^2/mol/A
                                               !  / 1d24     & ! kg m^2/ps^2/mol/A
                                               !  * 1d20     & ! kG A^2/ps^2/mol/A
                                               !  / mass     & ! k A/ps^2
                                               !  * 1d3        ! A/ps^2
 

     !calculate velocity
     velocity(:,1:num_molecule) = velocity(:,1:num_molecule) + accel(:,1:num_molecule)*dt
     !calculate position
     position(:,1:num_molecule) = position(:,1:num_molecule) + velocity(:,1:num_molecule)*dt
     !calculate energies
     ek = 0.0
     do j=1,num_molecule
        ek = ek + argon_mass * (velocity(1,j)**2 + velocity(2,j)**2 + velocity(3,j)**2) / 2.0
     enddo
     ! (g/mol) * (A/ps)**2
     ! = (g A**2 / ps**2) / mol
     ! = (1e-23 kg m**2/ (1e-24 s**2) / mol
     ! = 10 J/mol
     ek = ek * 0.01 ! kJ/mol
     ! 3NkbT = 2Ek
     temperature = ek * 2d0 /(3d0 * num_molecule * 0.008314)
     ! Write not so frequently
     if ( mod(i,1000) == 0 ) then
        ! time step in ps, total energies in kJ/mol, temperature in K
        write(6,*) i,i*dt, ep,ek,ep+ek, temperature
     endif
  enddo
  write(6,*) num_molecule
  do i=1,num_molecule
     write(6,*) (position(j,i),j=1,3),(velocity(j,i),j=1,3)
  enddo
  deallocate(position)
  deallocate(velocity)
  deallocate(accel)
  deallocate(force)
end program clustermd
