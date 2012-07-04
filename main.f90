module system_module
  implicit none
  !system constants
  integer, parameter :: STDIN = 5, STDOUT = 6, STDERR = 7
end module system_module


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
    !reserved for future extension
  end subroutine physconst_init

end module physconst_module



module monatom_module
  implicit none
  type monatom_type
     integer      :: num_molecule
     real(kind=8), pointer, dimension(:,:) :: position  ! Angstrom
     real(kind=8), pointer, dimension(:,:) :: velocity  ! Angstrom / ps
     real(kind=8), pointer, dimension(:,:) :: accel     ! Angstrom / ps**2
     real(kind=8), pointer, dimension(:,:) :: force     ! N
     real(kind=8) :: mass                        ! atomic mass
  end type monatom_type

contains

  subroutine monatom_allocate(m, nmol)
    type(monatom_type), intent(OUT) :: m
    integer, intent(IN) :: nmol
    m%num_molecule = nmol
    allocate(m%position(3,nmol))
    allocate(m%velocity(3,nmol))
    allocate(m%accel(3,nmol))
    allocate(m%force(3,nmol))
  end subroutine monatom_allocate


  subroutine monatom_done(m)
    type(monatom_type), intent(INOUT) :: m
    deallocate(m%position)
    deallocate(m%velocity)
    deallocate(m%accel)
    deallocate(m%force)
  end subroutine monatom_done


  subroutine monatom_resetforce(m)
    type(monatom_type), intent(INOUT) :: m
    m%force(:,:) = 0d0
  end subroutine monatom_resetforce


  subroutine monatom_proceed_velocity(m, deltat)
    type(monatom_type), intent(INOUT) :: m
    real(kind=8), intent(IN) :: deltat
    m%velocity(:,1:m%num_molecule) = m%velocity(:,1:m%num_molecule) + m%accel(:,1:m%num_molecule)*deltat
  end subroutine monatom_proceed_velocity


  subroutine monatom_proceed_position(m, deltat)
    type(monatom_type), intent(INOUT) :: m
    real(kind=8), intent(IN) :: deltat
    m%position(:,1:m%num_molecule) = m%position(:,1:m%num_molecule) + m%velocity(:,1:m%num_molecule)*deltat
  end subroutine monatom_proceed_position


  subroutine monatom_write(m,filehandle)
    type(monatom_type), intent(IN) :: m
    integer, intent(IN) :: filehandle
    integer :: i,j
    write(filehandle,'("[ATOMMAS]")')
    write(filehandle,*) m%mass
    write(filehandle,'("[ATOMPOSVEL]")')
    write(filehandle,*) m%num_molecule
    do i=1,m%num_molecule
       write(filehandle,*) (m%position(j,i),j=1,3),(m%velocity(j,i),j=1,3)
    enddo
  end subroutine monatom_write


  subroutine monatom_accel_from_force(m)
    type(monatom_type), intent(INOUT) :: m
    m%accel(:,1:m%num_molecule) = m%force(:,1:m%num_molecule) / m%mass * 100.0
  end subroutine monatom_accel_from_force


  function monatom_kinetic_energy(m) result(ek)
    type(monatom_type), intent(IN) :: m
    real(kind=8) :: ek
    integer      :: j
    ek = 0.0
    do j=1,m%num_molecule
       ek = ek + m%mass * (m%velocity(1,j)**2 + m%velocity(2,j)**2 + m%velocity(3,j)**2) / 2.0
    enddo
    ! (g/mol) * (A/ps)**2
    ! = (g A**2 / ps**2) / mol
    ! = (1e-23 kg m**2/ (1e-24 s**2) / mol
    ! = 10 J/mol
    ek = ek * 0.01 ! kJ/mol
  end function monatom_kinetic_energy


  function monatom_degree_of_freedom(m) result(dof)
    type(monatom_type), intent(INOUT) :: m
    integer :: dof
    dof = m%num_molecule * 3
  end function monatom_degree_of_freedom


  subroutine monatom_read_atomposvel(m, filehandle)
    type(monatom_type), intent(INOUT) :: m
    integer, intent(IN) :: filehandle
    integer :: i,j
    do i=1,m%num_molecule
       read(filehandle,*) (m%position(j,i),j=1,3), (m%velocity(j,i),j=1,3)
    enddo
  end subroutine monatom_read_atomposvel


  subroutine monatom_read_atompos(m, filehandle)
    type(monatom_type), intent(INOUT) :: m
    integer, intent(IN) :: filehandle
    integer :: i,j
    do i=1,m%num_molecule
       read(filehandle,*) (m%position(j,i),j=1,3)
    enddo
    m%velocity(:,1:m%num_molecule) = 0.0d0
  end subroutine monatom_read_atompos

end module monatom_module


module properties_module
  use monatom_module
  implicit none
  integer :: num_monatom
  type(monatom_type) :: monatoms(100)
  character(len=256) :: label(100)
contains

  subroutine properties_init
    num_monatom = 0
  end subroutine properties_init

  subroutine properties_read(tag, filehandle)
    character(len=*), intent(IN) :: tag
    integer, intent(IN)          :: filehandle
    integer                      :: nmol
    if ( tag == "[COMPONENT]" ) then
       num_monatom = num_monatom + 1
       read(filehandle,*) label(num_monatom)
    else if ( tag == "[ATOMPOS]" ) then
       read(filehandle,*) nmol
       call monatom_allocate(monatoms(num_monatom), nmol)
       call monatom_read_atompos(monatoms(num_monatom), filehandle)
    else if ( tag == "[ATOMPOSVEL]" ) then
       read(filehandle,*) nmol
       call monatom_allocate(monatoms(num_monatom), nmol)
       call monatom_read_atomposvel(monatoms(num_monatom), filehandle)
    else if ( tag == "[ATOMMAS]" ) then
       read(filehandle,*) monatoms(num_monatom)%mass
    end if
  end subroutine properties_read


  subroutine properties_write(filehandle)
    integer, intent(IN) :: filehandle
    integer :: k
    do k=1,num_monatom
       write(filehandle,'("[COMPONENT]")')
       write(filehandle,'(a)') label(k)
       call monatom_write(monatoms(k), filehandle)
    end do
  end subroutine properties_write


  subroutine properties_accel_from_force
    integer :: k
    do k=1,num_monatom
       call monatom_accel_from_force(monatoms(k))
    end do
  end subroutine properties_accel_from_force


  function properties_kinetic_energy() result(ek)
    integer :: k
    real(kind=8) :: ek
    ek = 0d0
    do k=1,num_monatom
       ek = ek + monatom_kinetic_energy(monatoms(k))
    end do
  end function properties_kinetic_energy


  function properties_degree_of_freedom() result(dof)
    integer :: k
    integer :: dof
    dof = 0d0
    do k=1,num_monatom
       dof = dof + monatom_degree_of_freedom(monatoms(k))
    end do
  end function properties_degree_of_freedom


  function properties_temperature(ek) result(temperature)
    use physconst_module
    real(kind=8), intent(IN) :: ek
    real(kind=8)             :: temperature
    temperature = ek * 2d0 /(properties_degree_of_freedom() * kB)
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
    integer :: k
    do k=1,num_monatom
       call monatom_done(monatoms(k))
    end do
  end subroutine properties_done


  subroutine properties_preforce
    integer :: k
    do k=1,num_monatom
       call monatom_resetforce(monatoms(k))
    end do
  end subroutine properties_preforce


  subroutine properties_postforce
  end subroutine properties_postforce

end module properties_module



module lj_module
  implicit none
  type lj_type
     real(kind=8) :: eps,sig                     ! kJ/mol, Angstrom
  end type lj_type
  
contains

  subroutine lj_init(x)
    type(lj_type), intent(OUT) :: x
    x%eps = 0.0
    x%sig = 0.0
  end subroutine lj_init


  function lj_is_available(x) result(logic)
    type(lj_type), intent(OUT) :: x
    logical :: logic
    logic = x%eps /= 0.0
  end function lj_is_available


  subroutine lj_write(x, filehandle)
    type(lj_type), intent(IN) :: x
    integer, intent(IN) :: filehandle
    write(filehandle,'("[LJPARAM]")')
    write(filehandle,*) x%eps,x%sig
  end subroutine lj_write


  subroutine lj_calculate_homo(x, m, ep, vir_ex)
    use properties_module
    use box_module
    type(lj_type), intent(IN) :: x
    type(monatom_type), intent(INOUT) :: m
    real(kind=8), intent(OUT):: ep,vir_ex
    real(kind=8)             :: delta(3),dd
    real(kind=8)             :: k
    integer :: j1,j2,j
    ep = 0.0
    !excess virial by interactions
    vir_ex = 0.0 !kJ/mol
    !interaction between all pairs of molecules
    do j1 = 1,m%num_molecule
       do j2 = j1+1, m%num_molecule
          delta(:) = m%position(:,j2) - m%position(:,j1)
          !minimum image
          if ( box_is_available() ) then
             delta(:) = delta(:) - dnint(delta(:) / box(:)) * box(:)
          endif
          dd = delta(1)**2 + delta(2)**2 + delta(3)**2
          ep = ep + 4*x%eps*(x%sig**12/dd**6 - x%sig**6/dd**3)      ! kJ/mol
          k  =    - 4*x%eps*(12*x%sig**12/dd**7 - 6*x%sig**6/dd**4) ! kJ/mol/Ang**2
          m%force(:,j1) = m%force(:,j1) + k * delta(:)    ! kJ/mol/Ang
          m%force(:,j2) = m%force(:,j2) - k * delta(:)
          do j = 1,3
             vir_ex = vir_ex + delta(j) * ( -k * delta(j) )
          enddo
       enddo
    enddo
    vir_ex = vir_ex / 3.0d0  ! W defined in Allen&Tildesley p.48
  end subroutine lj_calculate_homo


  subroutine lj_calculate_hetero(x, m1, m2, ep, vir_ex)
    use properties_module
    use box_module
    type(lj_type), intent(IN) :: x
    type(monatom_type), intent(INOUT) :: m1,m2
    real(kind=8), intent(OUT):: ep,vir_ex
    real(kind=8)             :: delta(3),dd
    real(kind=8)             :: k
    integer :: j1,j2,j
    ep = 0.0
    !excess virial by interactions
    vir_ex = 0.0 !kJ/mol
    !interaction between all pairs of molecules
    do j1 = 1,m1%num_molecule
       do j2 = 1, m2%num_molecule
          delta(:) = m2%position(:,j2) - m1%position(:,j1)
          !minimum image
          if ( box_is_available() ) then
             delta(:) = delta(:) - dnint(delta(:) / box(:)) * box(:)
          endif
          dd = delta(1)**2 + delta(2)**2 + delta(3)**2
          ep = ep + 4*x%eps*(x%sig**12/dd**6 - x%sig**6/dd**3)      ! kJ/mol
          k  =    - 4*x%eps*(12*x%sig**12/dd**7 - 6*x%sig**6/dd**4) ! kJ/mol/Ang**2
          m1%force(:,j1) = m1%force(:,j1) + k * delta(:)    ! kJ/mol/Ang
          m2%force(:,j2) = m2%force(:,j2) - k * delta(:)
          do j = 1,3
             vir_ex = vir_ex + delta(j) * ( -k * delta(j) )
          enddo
       enddo
    enddo
    vir_ex = vir_ex / 3.0d0  ! W defined in Allen&Tildesley p.48
  end subroutine lj_calculate_hetero


  subroutine lj_lorentz_berthelot_rule(lj1,lj2,eps,sig)
    type(lj_type), intent(IN) :: lj1,lj2
    real(kind=8), intent(OUT) :: eps, sig
    eps = sqrt(lj1%eps * lj2%eps)
    sig = (lj1%sig + lj2%sig)/2d0
  end subroutine lj_lorentz_berthelot_rule

end module lj_module



module interaction_module
  use lj_module
  implicit none
  type(lj_type) :: lj_pair(100,100)

contains

  subroutine interaction_init
    integer :: i,j
    do i=1,100
       do j=1,100
          call lj_init(lj_pair(i,j))
       enddo
    enddo
  end subroutine interaction_init


  subroutine interaction_read(tag, filehandle)
    use system_module
    character(len=*), intent(IN) :: tag
    integer, intent(IN)          :: filehandle
    integer :: i,j
    character(len=256) :: intrtype
    real(kind=8) :: eps, sig
    if ( tag == "[INTRPAIR]" ) then
       read(filehandle,*) i,j,intrtype
       if ( intrtype == "LJ" ) then
          read(filehandle,*) eps,sig
          lj_pair(i,j)%eps = eps
          lj_pair(i,j)%sig = sig
       else if ( intrtype == "LB" ) then
          call lj_lorentz_berthelot_rule(lj_pair(i,i),lj_pair(j,j),eps,sig)
          lj_pair(i,j)%eps = eps
          lj_pair(i,j)%sig = sig
          lj_pair(j,i)%eps = eps
          lj_pair(j,i)%sig = sig
       else
          write(STDERR,*) "UNKNOWN INTERACTION TYPE:", intrtype
       endif
    endif
  end subroutine interaction_read


  subroutine interaction_write(filehandle)
    integer, intent(IN)          :: filehandle
    integer :: i,j
    do i=1,100
       do j=i,100
          if ( lj_is_available(lj_pair(i,j)) ) then
             write(filehandle,'("[INTRPAIR]")')
             write(filehandle,*) i,j,'LJ'
             write(filehandle,*) lj_pair(i,j)%eps, lj_pair(i,j)%sig
          endif
       enddo
    enddo
  end subroutine interaction_write


  
  subroutine interaction_calculate(ep, vir_ex)
    use properties_module
    real(kind=8), intent(OUT):: ep,vir_ex
    real(kind=8)             :: ep1,vir_ex1
    integer :: k, k1, k2
    ep = 0.0
    !excess virial by interactions
    vir_ex = 0.0 !kJ/mol
    do k=1,num_monatom
       if ( lj_is_available(lj_pair(k,k)) ) then
          call lj_calculate_homo(lj_pair(k,k), monatoms(k), ep1, vir_ex1)
          ep = ep + ep1
          vir_ex = vir_ex + vir_ex1
       endif
    enddo
    do k1=1,num_monatom
       do k2=k1+1,num_monatom
          if ( lj_is_available(lj_pair(k1,k2)) ) then
             call lj_calculate_hetero(lj_pair(k1,k2), monatoms(k1), monatoms(k2), ep1, vir_ex1)
             ep = ep + ep1
             vir_ex = vir_ex + vir_ex1
          endif
       enddo
    enddo
  end subroutine interaction_calculate
       
end module interaction_module



module integrator_module
  implicit none

contains

  subroutine integrator_init
  end subroutine integrator_init


  subroutine integrator_proceed_velocity(deltat)
    use properties_module
    real(kind=8), intent(IN) :: deltat
    integer :: k
    do k=1,num_monatom
       call monatom_proceed_velocity(monatoms(k), deltat)
    end do
  end subroutine integrator_proceed_velocity


  subroutine integrator_proceed_position(deltat)
    use properties_module
    real(kind=8), intent(IN) :: deltat
    integer :: k
    do k=1,num_monatom
       call monatom_proceed_position(monatoms(k), deltat)
    end do
  end subroutine integrator_proceed_position

end module integrator_module



module settings_module
  implicit none
  real(kind=8) :: dt                          ! picoseconds
  integer :: num_loop
  integer :: log_interval
  real(kind=8) :: lasttime

contains

  subroutine settings_init
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
  use interaction_module
  use settings_module
  use box_module
  use system_module
  implicit none
  !local variables
  real(kind=8) :: ep, vir_ex
  integer      :: i
  character(len=1000) :: tag
  !default values
  call box_init
  call physconst_init
  call properties_init
  call integrator_init
  call interaction_init
  call settings_init
  do
     read(STDIN,*,end=999) tag
     call properties_read(tag, STDIN)
     call settings_read(tag, STDIN)
     call box_read(tag, STDIN)
     call interaction_read(tag, STDIN)
  end do
999 continue

  !realistic interaction parameters for Argon.
  do i=1,num_loop
     lasttime = lasttime + dt
     !calculate position
     call integrator_proceed_position(dt/2)
     !calculate force
     call properties_preforce
     call interaction_calculate(ep, vir_ex)
     call properties_postforce
     !calculate accel
     call properties_accel_from_force
     !calculate velocity
     call integrator_proceed_velocity(dt)
     !calculate position
     call integrator_proceed_position(dt/2)
     ! logging
     call settings_write_log(i, ep, vir_ex, STDOUT)
  enddo
  call box_write(STDOUT)
  call settings_write(STDOUT)
  call properties_write(STDOUT)
  call interaction_write(STDOUT)
  !deallocate arrays
  call properties_done
end program main
