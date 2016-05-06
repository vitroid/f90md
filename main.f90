module box_module
  implicit none
  real(kind=8) :: box(3), volume

contains

  subroutine box_init()
    box(1) = -1d0 ! negative value == no box
  end subroutine box_init

  
  subroutine box_read(tag, filehandle)
    character(len=*), intent(IN) :: tag
    integer, intent(IN)          :: filehandle
    integer                      :: i
    if ( tag == "[CUBOIDBOX]" ) then
       read(filehandle,*) (box(i),i=1,3)
       volume = box(1)*box(2)*box(3)
    end if
  end subroutine box_read


  subroutine box_write(filehandle)
    integer, intent(IN)          :: filehandle
    integer                      :: i
    if ( box(1) > 0.0d0 ) then
       write(filehandle,'("[CUBOIDBOX]")')
       write(filehandle,*) (box(i),i=1,3)
    endif
  end subroutine box_write


  function box_is_available() result(logic)
    logical :: logic
    logic = box(1) > 0d0
  end function box_is_available

end module box_module



module physconst_module
  implicit none
  real(kind=8), parameter :: kB = 0.008314 ! kJ/mol/K
  real(kind=8), parameter :: NA = 6.022d23
contains
  subroutine physconst_init
  end subroutine physconst_init

end module physconst_module



module properties_module
  implicit none
  integer      :: num_molecule
  real(kind=8), allocatable :: position(:,:)  ! Angstrom
  real(kind=8), allocatable :: velocity(:,:)  ! Angstrom / ps
  real(kind=8), allocatable :: accel(:,:)     ! Angstrom / ps**2
  real(kind=8), allocatable :: force(:,:)     ! N
  real(kind=8) :: mass                        ! atomic mass
  real(kind=8) :: eps,sig                     ! kJ/mol, Angstrom

contains

  subroutine properties_init
    mass = 39.95d0
    eps = 0.99768d0
    sig = 3.41d0
  end subroutine properties_init


  subroutine properties_read(tag, filehandle)
    character(len=*), intent(IN) :: tag
    integer, intent(IN)          :: filehandle
    integer                      :: i,j
    if ( tag == "[ATOMPOS]" ) then
       read(filehandle,*) num_molecule
       allocate(position(3,num_molecule))
       allocate(velocity(3,num_molecule))
       allocate(accel(3,num_molecule))
       allocate(force(3,num_molecule))
       do i=1,num_molecule
          read(filehandle,*) (position(j,i),j=1,3)   ! coordinates in Angstrom
       enddo
       velocity(:,1:num_molecule) = 0.0d0
    else if ( tag == "[ATOMPOSVEL]" ) then
       read(filehandle,*) num_molecule
       allocate(position(3,num_molecule))
       allocate(velocity(3,num_molecule))
       allocate(accel(3,num_molecule))
       allocate(force(3,num_molecule))
       do i=1,num_molecule
          read(filehandle,*) (position(j,i),j=1,3),(velocity(j,i),j=1,3)
       enddo
    else if ( tag == "[ATOMMAS]" ) then
       read(filehandle,*) mass
    else if ( tag == "[LJPARAM]" ) then
       read(filehandle,*) eps,sig
    end if
  end subroutine properties_read


  subroutine properties_write(filehandle)
    integer, intent(IN) :: filehandle
    integer :: i,j
    write(filehandle,'("[LJPARAM]")')
    write(filehandle,*) eps,sig
    write(filehandle,'("[ATOMMAS]")')
    write(filehandle,*) mass
    write(filehandle,'("[ATOMPOSVEL]")')
    write(filehandle,*) num_molecule
    do i=1,num_molecule
       write(filehandle,*) (position(j,i),j=1,3),(velocity(j,i),j=1,3)
    enddo
  end subroutine properties_write


  subroutine properties_accel_from_force()
    accel(:,1:num_molecule) = force(:,1:num_molecule) / mass * 100.0
  end subroutine properties_accel_from_force


  function properties_kinetic_energy() result(ek)
    real(kind=8) :: ek
    integer      :: j
    ek = 0.0
    do j=1,num_molecule
       ek = ek + mass * (velocity(1,j)**2 + velocity(2,j)**2 + velocity(3,j)**2) / 2.0
    enddo
    ! (g/mol) * (A/ps)**2
    ! = (g A**2 / ps**2) / mol
    ! = (1e-23 kg m**2/ (1e-24 s**2) / mol
    ! = 10 J/mol
    ek = ek * 0.01 ! kJ/mol
  end function properties_kinetic_energy


  function properties_temperature(ek) result(temperature)
    use physconst_module
    real(kind=8), intent(IN) :: ek
    real(kind=8)             :: temperature
    temperature = ek * 2d0 /(3d0 * num_molecule * kB)
  end function properties_temperature


  function properties_pressure(vir_ex, ek, volume) result(pressure)
    use physconst_module
    real(kind=8), intent(IN) :: vir_ex, ek, volume
    real(kind=8) :: pressure, vir_id
    vir_id = 2.0 * ek / 3.0
    pressure = (vir_id + vir_ex) / volume *1e33 / NA
    ! *1.0     ! kJ/mol / A**3
    ! *1e3     ! J/mol / A**3
    ! *1e30    ! J/mol / m**3
    !          ! J / m / mol / m**2
    !          ! Pa / mol
    ! /NA      ! Pa
   end function properties_pressure


   subroutine properties_done
     deallocate(position)
     deallocate(velocity)
     deallocate(accel)
     deallocate(force)
   end subroutine properties_done

end module properties_module



module force_module
  implicit none

contains

  subroutine force_init
  end subroutine force_init


  subroutine force_calculate(ep, vir_ex)
    use properties_module
    use box_module
    real(kind=8), intent(OUT):: ep,vir_ex
    real(kind=8)             :: delta(3),dd
    real(kind=8)             :: k
    integer :: j1,j2,j
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
          if ( box_is_available() ) then
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
  end subroutine force_calculate

end module force_module



module integrator_module
  implicit none

contains

  subroutine integrator_init
  end subroutine integrator_init


  subroutine integrator_proceed_velocity(deltat)
    use properties_module
    real(kind=8), intent(IN) :: deltat
    velocity(:,1:num_molecule) = velocity(:,1:num_molecule) + accel(:,1:num_molecule)*deltat
  end subroutine integrator_proceed_velocity


  subroutine integrator_proceed_position(deltat)
    use properties_module
    real(kind=8), intent(IN) :: deltat
    position(:,1:num_molecule) = position(:,1:num_molecule) + velocity(:,1:num_molecule)*deltat
  end subroutine integrator_proceed_position

end module integrator_module



module settings_module
  implicit none
  real(kind=8) :: dt                          ! picoseconds
  integer :: num_loop
  integer :: log_interval
  real(kind=8) :: lasttime

contains

  subroutine settings_init()
    dt = 0.001
    num_loop = 1000000
    lasttime = 0.0d0
    log_interval = 0
  end subroutine settings_init


  subroutine settings_read(tag, filehandle)
    character(len=*), intent(IN) :: tag
    integer, intent(IN)          :: filehandle
    if ( tag == "[INTVPS]" ) then
       read(filehandle,*) dt
    else if ( tag == "[STEPS]" ) then
       read(filehandle,*) num_loop
    else if ( tag == "[LASTTIME]" ) then
       read(filehandle,*) lasttime
    else if ( tag == "[LOGINTV]" ) then
       read(filehandle,*) log_interval
    endif
  end subroutine settings_read


  subroutine settings_write(filehandle)
    integer, intent(IN)          :: filehandle
    write(filehandle,'("[INTVPS]")')
    write(filehandle,*) dt
    write(filehandle,'("[LASTTIME]")')
    write(filehandle,*) lasttime
    write(filehandle,'("[STEPS]")')
    write(filehandle,*) num_loop
    write(filehandle,'("[LOGINTV]")')
    write(filehandle,*) log_interval
  end subroutine settings_write


  subroutine settings_write_log(i, ep, vir_ex, filehandle)
    use properties_module
    use box_module
    integer, intent(IN)          :: i, filehandle
    real(kind=8), intent(IN)     :: ep, vir_ex
    real(kind=8)                 :: ek
    real(kind=8)                 :: temperature, pressure
    ! Write not so frequently
    if ( log_interval > 0 .and. mod(i,log_interval) == 0 ) then
       !calculate energies
       ek = properties_kinetic_energy()
       temperature = properties_temperature(ek)
       ! 3NkbT = 2Ek
       if ( box(1) > 0.0 ) then
          pressure = properties_pressure(vir_ex, ek, volume)
       endif
       
       ! time step in ps, total energies in kJ/mol, temperature in K, pressure in Pa
       if ( box(1) > 0.0 ) then
          write(filehandle,'("[LOGV2]")')
          write(filehandle,*) i,lasttime, ep,ek,ep+ek, temperature, pressure, volume
       else
          write(filehandle,'("[LOGV1]")')
          write(filehandle,*) i,lasttime, ep,ek,ep+ek, temperature
       endif
    endif
  end subroutine settings_write_log

end module settings_module



program main
  use physconst_module
  use properties_module
  use integrator_module
  use force_module
  use settings_module
  use box_module
  implicit none
  !system constants
  integer, parameter :: STDIN = 5, STDOUT = 6
  !local variables
  real(kind=8) :: ep, vir_ex
  integer      :: i
  character(len=1000) :: tag
  !default values
  call box_init
  call physconst_init
  call properties_init
  call integrator_init
  call force_init
  call settings_init
  do
     read(STDIN,*,end=999) tag
     call properties_read(tag, STDIN)
     call settings_read(tag, STDIN)
     call box_read(tag, STDIN)
  end do
999 continue

  !realistic interaction parameters for Argon.
  do i=1,num_loop
     lasttime = lasttime + dt
     !Position Verlet
     !See https://en.wikipedia.org/wiki/Verlet_integration
     !calculate position at the future by dt/2
     call integrator_proceed_position(dt/2)
     !calculate force at the future by dt/2
     call force_calculate(ep, vir_ex)
     !calculate accel at the future by dt/2
     call properties_accel_from_force
     !calculate velocity at the future by dt
     call integrator_proceed_velocity(dt)
     !calculate position at the future by dt (=dt/2+dt/2)
     call integrator_proceed_position(dt/2)
     !Now both the position and the velocity are proceeded.
     ! logging
     call settings_write_log(i, ep, vir_ex, STDOUT)
  enddo
  call box_write(STDOUT)
  call settings_write(STDOUT)
  call properties_write(STDOUT)
  !deallocate arrays
  call properties_done
end program main
