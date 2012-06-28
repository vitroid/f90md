program clustermd
  implicit none
  !molecular properties
  real(kind=8) :: x(2),y(2),z(2)
  real(kind=8) :: vx(2),vy(2),vz(2)
  real(kind=8) :: ax(2),ay(2),az(2)
  real(kind=8) :: fx(2),fy(2),fz(2)
  real(kind=8) :: mass
  real(kind=8) :: dt
  real(kind=8) :: k

  !local variables
  integer :: num_loop
  integer :: i,j
  real(kind=8) :: dx,dy,dz,dd
  real(kind=8) :: ep,ek
  real(kind=8) :: eps,sig
  x(1) = 0.0d0
  y(1) = 0.0d0
  z(1) = 0.0d0
  x(2) = 3.0d0
  y(2) = 3.0d0
  z(2) = 3.0d0
  do i=1,2
     vx(i) = 0.0d0
     vy(i) = 0.0d0
     vz(i) = 0.0d0
     ax(i) = 0.0d0
     ay(i) = 0.0d0
     az(i) = 0.0d0
     fx(i) = 0.0d0
     fy(i) = 0.0d0
     fz(i) = 0.0d0
  enddo
  dt = 0.001
  mass = 1.0
  eps = 3.0
  sig = 3.0
  num_loop = 100000
  do i=1,num_loop
     !calculate force
     ep = 0.0
     dx = x(2) - x(1)
     dy = y(2) - y(1)
     dz = z(2) - z(1)
     dd = dx**2 + dy**2 + dz**2
     ep = 4*eps*(sig**12/dd**6 - sig**6/dd**3)
     k  =-4*eps*(12*sig**12/dd**7 - 6*sig**6/dd**4)
     fx(1) = +k * dx
     fy(1) = +k * dy
     fz(1) = +k * dz
     fx(2) = -k * dx
     fy(2) = -k * dy
     fz(2) = -k * dz
     !calculate accel
     do j=1,2
        ax(j) = fx(j) / mass
        ay(j) = fy(j) / mass
        az(j) = fz(j) / mass
     enddo
     !calculate velocity
     do j=1,2
        vx(j) = vx(j) + ax(j)*dt
        vy(j) = vy(j) + ay(j)*dt
        vz(j) = vz(j) + az(j)*dt
     enddo
     !calculate position
     do j=1,2
        x(j) = x(j) + vx(j)*dt
        y(j) = y(j) + vy(j)*dt
        z(j) = z(j) + vz(j)*dt
     enddo
     !calculate energies
     ek = 0.0
     do j=1,2
        ek = ek + mass * (vx(j)**2 + vy(j)**2 + vz(j)**2) / 2.0
     enddo
     write(6,*) i,ep,ek,ep+ek,sqrt(dd)
  enddo
end program clustermd
