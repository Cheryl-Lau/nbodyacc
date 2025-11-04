
module centreofmass

 implicit none  
 public :: transform_com_frame

 private 

 contains 

subroutine transform_com_frame(nptmass,xyzhm_ptmass,vxyz_ptmass)
 integer, intent(in)    :: nptmass
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 integer :: i
 real    :: mass_tot,xyz_com(3),vxyz_com(3),mass_i,xyz_i(3),vxyz_i(3)

 !--Total mass of system
 mass_tot = 0.d0
 do i = 1,nptmass
    mass_tot = mass_tot + xyzhm_ptmass(5,i)
 enddo

 !--Compute position and velocity of COM
 xyz_com = (/ 0.d0, 0.d0, 0.d0 /)
 vxyz_com = (/ 0.d0, 0.d0, 0.d0 /)
 over_parts: do i = 1,nptmass
    mass_i = xyzhm_ptmass(5,i)
    xyz_i  = xyzhm_ptmass(1:3,i)
    vxyz_i = vxyz_ptmass(1:3,i)
    xyz_com  = xyz_com  + mass_i/mass_tot * xyz_i
    vxyz_com = vxyz_com + mass_i/mass_tot * vxyz_i
 enddo over_parts

 !--Transform all x and v
 each_body: do i = 1,nptmass
    xyzhm_ptmass(1:3,i) = xyzhm_ptmass(1:3,i) - xyz_com
    vxyz_ptmass(1:3,i)  = vxyz_ptmass(1:3,i) - vxyz_com
 enddo each_body

end subroutine transform_com_frame

end module centreofmass






































