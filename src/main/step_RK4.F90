
module step_RK4

 implicit none 
 public :: step 

 private 

contains 

subroutine step(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass)
 integer, intent(in)    :: nptmass
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 

end subroutine step

end module step_RK4