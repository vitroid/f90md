program clustermd
  implicit none
  !molecular properties
  !Make arrays allocatable
  integer      :: num_molecule
  real(kind=8), allocatable :: position(:,:)  ! Angstrom
  real(kind=8), allocatable :: velocity(:,:)  ! Angstrom / ps
  real(kind=8), allocatable :: accel(:,:)     ! Angstrom / ps**2
  real(kind=8), allocatable :: force(:,:)     ! N
  real(kind=8) :: mass                        ! atomic mass
  real(kind=8) :: eps,sig                     ! kJ/mol, Angstrom
  real(kind=8) :: dt                          ! picoseconds
  real(kind=8) :: k
  real(kind=8) :: temperature
  ! volume in Ang**3, ideal and excess virial in kJ/mol, and pressure in Pa.
  real(kind=8) :: volume, vir_id, vir_ex, pressure
  !local variables
  integer :: num_loop
  integer :: i,j
  integer :: ix,iy,iz
  integer :: j1,j2
  integer :: log_interval
  real(kind=8) :: delta(3),dd
  real(kind=8) :: ep,ek
  real(kind=8) :: lasttime
  real(kind=8) :: box(3)
  character(len=1000) :: tag
  !default values
  dt = 0.001
  num_loop = 1000000
  mass = 39.95d0
  eps = 0.99768d0
  sig = 3.41d0
  lasttime = 0.0d0
  log_interval = 0
  box(1) = -1d0 ! negative value == no box
  !read various data according to the tags.
  !tags should not be very descriptive. Tags should rather be a symbolic name.
  !You may want to change the data format later.
  !But if you change the format, you must also change the tag because the tag represents
  !the following data format.
  !For example, you first decide to use the tag [Argon position] for the coordinate of Argon 
  !atoms, and later you found the file format was incomplete and want to change format.
  !The new data format should not be tagged [Argon position] because the format is different.
  !If you use the same tag for different input data,
  !the old programs that want the [Argon position] type data will confuse.
  do
     read(5,*,end=999) tag
     if ( tag == "[ATOMPOS]" ) then
        read(5,*) num_molecule
        allocate(position(3,num_molecule))
        allocate(velocity(3,num_molecule))
        do i=1,num_molecule
           read(5,*) (position(j,i),j=1,3)   ! coordinates in Angstrom
        enddo
        velocity(:,1:num_molecule) = 0.0d0
     else if ( tag == "[ATOMPOSVEL]" ) then
        read(5,*) num_molecule
        allocate(position(3,num_molecule))
        allocate(velocity(3,num_molecule))
        do i=1,num_molecule
           read(5,*) (position(j,i),j=1,3),(velocity(j,i),j=1,3)
        enddo
     else if ( tag == "[INTVPS]" ) then
        read(5,*) dt
     else if ( tag == "[ATOMMAS]" ) then
        read(5,*) mass
     else if ( tag == "[LJPARAM]" ) then
        read(5,*) eps,sig
     else if ( tag == "[STEPS]" ) then
        read(5,*) num_loop
     else if ( tag == "[LASTTIME]" ) then
        read(5,*) lasttime
     else if ( tag == "[LOGINTV]" ) then
        read(5,*) log_interval
     else if ( tag == "[CUBOIDBOX]" ) then
        read(5,*) (box(i),i=1,3)
        volume = box(1)*box(2)*box(3)
     end if
  end do
999 continue
  allocate(accel(3,num_molecule))
  allocate(force(3,num_molecule))
  accel(:,1:num_molecule) = 0.0d0

  !realistic interaction parameters for Argon.
  do i=1,num_loop
     !calculate force
     !force will be reset each step.
     force(:,1:num_molecule) = 0.0d0
     ep = 0.0
     !excess virial by interactions
     vir_ex = 0.0 !kJ/mol
     !interaction between all pairs of molecules
     do j1 = 1,num_molecule
        do j2 = j1+1, num_molecule
           delta(:) = position(:,j2) - position(:,j1)
           !minimum image
           if ( box(1) > 0.0d0 ) then
              delta(:) = delta(:) - dnint(delta(:) / box(:)) * box(:)
           endif
           dd = delta(1)**2 + delta(2)**2 + delta(3)**2
           ep = ep + 4*eps*(sig**12/dd**6 - sig**6/dd**3)      ! kJ/mol
           k  =    - 4*eps*(12*sig**12/dd**7 - 6*sig**6/dd**4) ! kJ/mol/Ang**2
           force(:,j1) = force(:,j1) + k * delta(:)    ! kJ/mol/Ang
           force(:,j2) = force(:,j2) - k * delta(:)
           do j = 1,3
              vir_ex = vir_ex + delta(j) * ( -k * delta(j) )
           enddo
        enddo
     enddo
     vir_ex = vir_ex / 3.0d0  ! W defined in Allen&Tildesley p.48
     !calculate accel
     !force is in kJ/mol/Angstrom
     !acc is in Angstrom/ps/ps
     accel(:,1:num_molecule) = force(:,1:num_molecule) / mass * 100.0
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
        ek = ek + mass * (velocity(1,j)**2 + velocity(2,j)**2 + velocity(3,j)**2) / 2.0
     enddo
     ! (g/mol) * (A/ps)**2
     ! = (g A**2 / ps**2) / mol
     ! = (1e-23 kg m**2/ (1e-24 s**2) / mol
     ! = 10 J/mol
     ek = ek * 0.01 ! kJ/mol
     ! 3NkbT = 2Ek
     temperature = ek * 2d0 /(3d0 * num_molecule * 0.008314)
     if ( box(1) > 0.0 ) then
        vir_id = 2.0 * ek / 3.0
        pressure = (vir_id + vir_ex) / volume *1e33 / 6.022e23
        ! *1.0     ! kJ/mol / A**3
        ! *1e3     ! J/mol / A**3
        ! *1e30    ! J/mol / m**3
        !          ! J / m / mol / m**2
        !          ! Pa / mol
        ! /NA      ! Pa
     endif

     lasttime = lasttime + dt
     ! Write not so frequently
     if ( log_interval > 0 .and. mod(i,log_interval) == 0 ) then
        ! time step in ps, total energies in kJ/mol, temperature in K, pressure in Pa
        if ( box(1) > 0.0 ) then
           write(6,'("[LOGV2]")')
           write(6,*) i,lasttime, ep,ek,ep+ek, temperature, pressure, volume
        else
           write(6,'("[LOGV1]")')
           write(6,*) i,lasttime, ep,ek,ep+ek, temperature
        endif
     endif
  enddo
  write(6,'("[INTVPS]")')
  write(6,*) dt
  write(6,'("[LJPARAM]")')
  write(6,*) eps,sig
  write(6,'("[ATOMMAS]")')
  write(6,*) mass
  write(6,'("[LASTTIME]")')
  write(6,*) lasttime
  write(6,'("[STEPS]")')
  write(6,*) num_loop
  write(6,'("[LOGINTV]")')
  write(6,*) log_interval
  if ( box(1) > 0.0d0 ) then
     write(6,'("[CUBOIDBOX]")')
     write(6,*) (box(i),i=1,3)
  endif
  write(6,'("[ATOMPOSVEL]")')
  write(6,*) num_molecule
  do i=1,num_molecule
     write(6,*) (position(j,i),j=1,3),(velocity(j,i),j=1,3)
  enddo
  deallocate(position)
  deallocate(velocity)
  deallocate(accel)
  deallocate(force)
end program clustermd
