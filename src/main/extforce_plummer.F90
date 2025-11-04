
module extforce_plummer 

 implicit none 
 public :: plummer_potential
 public :: read_infile_plummer,write_infile_plummer

 real, public :: Mclust = 1e3 
 real, public :: Rcore  = 0.2 

 private 

 namelist /plummer_params/ Mclust,Rcore 

contains 

!---------------------------------------------------------------
! Total force per unit mass and potential 
!---------------------------------------------------------------
subroutine plummer_potential(xi,yi,zi,extfxi,extfyi,extfzi,phi)
 real, intent(in)  :: xi,yi,zi
 real, intent(out) :: extfxi,extfyi,extfzi,phi
 real :: r2,Rcore2

 r2 = xi**2 + yi**2 + zi**2
 Rcore2 = Rcore**2 

 phi = -1.d0*Mclust/sqrt(r2 + Rcore2)

 extfxi = -1.d0*Mclust*(r2 + Rcore2)**(-3.d0/2.d0)*xi
 extfyi = -1.d0*Mclust*(r2 + Rcore2)**(-3.d0/2.d0)*yi
 extfzi = -1.d0*Mclust*(r2 + Rcore2)**(-3.d0/2.d0)*zi
 
end subroutine plummer_potential


!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_plummer(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=plummer_params,iostat=rc)
 if (rc /= 0) stop 'cannot read plummer options'

end subroutine read_infile_plummer


subroutine write_infile_plummer(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=plummer_params,iostat=rc)
 if (rc /= 0) stop 'cannot write plummer options'

end subroutine write_infile_plummer


end module extforce_plummer 