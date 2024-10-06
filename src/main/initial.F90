
module initial

 implicit none 
 public :: init

 private

contains 

subroutine init(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,massq_ptmass)
 integer, intent(inout) :: nptmass 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 real,    intent(inout), optional :: massq_ptmass(:)


end subroutine init


end module initial