program clustermd
  implicit none
  !molecular properties
  !Increase number of molecules
  integer      :: num_molecule
  real(kind=8) :: x(8),y(8),z(8)
  real(kind=8) :: vx(8),vy(8),vz(8)
  real(kind=8) :: ax(8),ay(8),az(8)
  real(kind=8) :: fx(8),fy(8),fz(8)
  real(kind=8) :: mass
  real(kind=8) :: dt
  real(kind=8) :: k

  !local variables
  integer :: num_loop
  integer :: i,j
  integer :: ix,iy,iz
  integer :: j1,j2
  real(kind=8) :: dx,dy,dz,dd
  real(kind=8) :: ep,ek
  real(kind=8) :: eps,sig
  i = 0
  !Place molecules on the lattice.
  do ix = 1,2
     do iy = 1,2
        do iz = 1,2
           i = i + 1
           x(i) = ix * 4.0d0
           y(i) = iy * 4.0d0
           z(i) = iz * 4.0d0
        enddo
     enddo
  enddo
  num_molecule = i
  do i=1,num_molecule
     vx(i) = 0.0d0
     vy(i) = 0.0d0
     vz(i) = 0.0d0
     ax(i) = 0.0d0
     ay(i) = 0.0d0
     az(i) = 0.0d0
     !force will be reset each step.
  enddo
  dt = 0.001
  !realistic interaction parameters for Argon.
  mass = 39.95d0
  eps = 120d0 * 0.008314d0
  sig = 3.41d0
  num_loop = 100000
  do i=1,num_loop
     !calculate force
     !force will be reset each step.
     do j=1,num_molecule
        fx(j) = 0.0d0
        fy(j) = 0.0d0
        fz(j) = 0.0d0
     enddo
     ep = 0.0
     !interaction between all pairs of molecules
     do j1 = 1,num_molecule
        do j2 = j1+1, num_molecule
           dx = x(j2) - x(j1)
           dy = y(j2) - y(j1)
           dz = z(j2) - z(j1)
           dd = dx**2 + dy**2 + dz**2
           ep = ep + 4*eps*(sig**12/dd**6 - sig**6/dd**3)
           k  =    - 4*eps*(12*sig**12/dd**7 - 6*sig**6/dd**4)
           fx(j1) = fx(j1) + k * dx
           fy(j1) = fy(j1) + k * dy
           fz(j1) = fz(j1) + k * dz
           fx(j2) = fx(j2) - k * dx
           fy(j2) = fy(j2) - k * dy
           fz(j2) = fz(j2) - k * dz
        enddo
     enddo
     !calculate accel
     do j=1,num_molecule
        ax(j) = fx(j) / mass
        ay(j) = fy(j) / mass
        az(j) = fz(j) / mass
     enddo
     !calculate velocity
     do j=1,num_molecule
        vx(j) = vx(j) + ax(j)*dt
        vy(j) = vy(j) + ay(j)*dt
        vz(j) = vz(j) + az(j)*dt
     enddo
     !calculate position
     do j=1,num_molecule
        x(j) = x(j) + vx(j)*dt
        y(j) = y(j) + vy(j)*dt
        z(j) = z(j) + vz(j)*dt
     enddo
     !calculate energies
     ek = 0.0
     do j=1,num_molecule
        ek = ek + mass * (vx(j)**2 + vy(j)**2 + vz(j)**2) / 2.0
     enddo
     write(6,*) i,ep,ek,ep+ek
  enddo
end program clustermd
