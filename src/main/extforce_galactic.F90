
module extforce_galactic 

 implicit none 
 public :: galactic_potential
 public :: read_infile_galactic,write_infile_galactic

 real, public :: Rc = 2.d-1
 real, public :: vc = 1.d0
 real, public :: q1 = 1.d0
 real, public :: q2 = 1.d0

 private 

 namelist /galactic_params/ Rc,vc,q1,q2

contains 

!---------------------------------------------------------------
! Total force per unit mass and potential 
!---------------------------------------------------------------
subroutine galactic_potential(xi,yi,zi,extfxi,extfyi,extfzi,phi)
 real, intent(in)  :: xi,yi,zi
 real, intent(out) :: extfxi,extfyi,extfzi,phi
 real :: fac,dphidx,dphidy,dphidz

 phi = 5.d-1*vc**2 * log(Rc**2 + xi**2 + yi**2/q1**2 + zi**2/q2**2)

 fac = vc**2/2.d0 * (Rc**2 + xi**2 + yi**2/q1**2 + zi**2/q2**2)**(-1)
 dphidx = fac*2.d0 * xi
 dphidy = fac*2.d0/q1**2 * yi
 dphidz = fac*2.d0/q2**2 * zi

 extfxi = -dphidx
 extfyi = -dphidy
 extfzi = -dphidz

end subroutine galactic_potential


!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_galactic(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=galactic_params,iostat=rc)
 if (rc /= 0) stop 'cannot read galactic options'

end subroutine read_infile_galactic


subroutine write_infile_galactic(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=galactic_params,iostat=rc)
 if (rc /= 0) stop 'cannot write galactic options'

end subroutine write_infile_galactic


end module extforce_galactic 