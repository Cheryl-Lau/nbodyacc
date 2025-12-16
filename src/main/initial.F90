
module initial

 implicit none 
 public :: init

 private

contains 

subroutine init(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,sq_ptmass)
 use force,         only:compute_forces,iextforce,iext_king
 use extforce_king, only:cluster_profile
 use accrete,       only:print_r_acc
 integer, intent(inout) :: nptmass 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 real,    intent(inout), optional :: sq_ptmass(:,:)

 if (iextforce == iext_king) call cluster_profile 

 call compute_forces(nptmass,xyzhm_ptmass,fxyz_ptmass)

 if (print_r_acc) then 
    open(2090,file='accretion_radii.ev',status='replace')
    write(2090,'(5A20)') 'time','sink ID','r_hill','r_bondihoyle','r_tidalneigh'
    close(2090)
 endif 

end subroutine init


end module initial