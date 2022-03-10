module force
  implicit none

contains

  subroutine force_calculate
    use settings
    use property
    !local variables
    type(vector3) :: delta,f
    integer       :: i,j,k
    real(kind=8)  :: dd,ddi,sddi,d6,d12
    
    do i=1,num_molecule
       force_acc(i)%vec(:) = 0d0
    enddo
    ep = 0d0
    do i=1,num_molecule
       do j=i+1,num_molecule
          delta%vec(:) = position(j)%vec(:) - position(i)%vec(:)
          dd = vector3_inner_product(delta,delta)
          ddi = 1d0 / dd
          sddi = argon_sigma**2 * ddi
          d6 = sddi**3
          d12 = d6**2
          ep = ep + argon_eps4 * (d12 - d6)
          f%vec(:) = argon_eps4 * (12d0*d12 - 6d0*d6)*delta%vec(:) * ddi
          force_acc(j)%vec(:) = force_acc(j)%vec(:) + f%vec(:)
          force_acc(i)%vec(:) = force_acc(i)%vec(:) - f%vec(:)
       enddo
    enddo
  end subroutine force_calculate

  subroutine force_to_accel
    use settings
    use property
    !local variables
    integer :: i
    do i=1,num_molecule
       !force is in kJ/mol/Angstrom
       !acc is in Angstrom/ps/ps
       force_acc(i)%vec(:) = force_acc(i)%vec(:) *1d2 / mass 
                                               !  * 1d3      & ! J/mol/A
                                               !             & ! (J == kg m^2/s^2)
                                               !             & ! kg m^2/s^2/mol/A
                                               !  / 1d24     & ! kg m^2/ps^2/mol/A
                                               !  * 1d20     & ! kG A^2/ps^2/mol/A
                                               !  / mass     & ! k A/ps^2
                                               !  * 1d3        ! A/ps^2
    enddo
  end subroutine force_to_accel
end module force
