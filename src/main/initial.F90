
module initial

 implicit none 
 public :: init

 private

contains 

subroutine init(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,sq_ptmass)
 use force,         only:iextforce,iext_king
 use extforce_king, only:cluster_profile
 integer, intent(inout) :: nptmass 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 real,    intent(inout), optional :: sq_ptmass(:,:)

 if (iextforce == iext_king) call cluster_profile

end subroutine init


end module initial