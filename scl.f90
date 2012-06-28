program simplecubiclattice
  implicit none
  integer      :: num_edge
  integer      :: ix,iy,iz
  real(kind=8) :: latticeconst
  num_edge = 3
  latticeconst = 4.0
  write(6,fmt='("[ATOMPOS]")')
  write(6,*) num_edge**3
  !Place molecules on the lattice.
  do ix = 1,num_edge
     do iy = 1,num_edge
        do iz = 1,num_edge
           write(6,*) ix*latticeconst, iy*latticeconst, iz*latticeconst
        enddo
     enddo
  enddo
end program simplecubiclattice
