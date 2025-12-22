
module initial

 implicit none 
 public :: init

 private

contains 

subroutine init(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,sq_ptmass)
 use force,         only:compute_forces,iextforce,iext_king
 use extforce_king, only:cluster_profile
 use accrete,       only:print_r_acc,print_j_range,print_v_rel
 integer, intent(inout) :: nptmass 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 real,    intent(inout), optional :: sq_ptmass(:,:)


 !-- Compute cluster potential 
 if (iextforce == iext_king) call cluster_profile 


 !--Init for dt constraint / leapfrogs
 call compute_forces(nptmass,xyzhm_ptmass,fxyz_ptmass)


 !--Write accretion info to file 
 if (print_r_acc) then 
    open(2090,file='accretion_radii.ev',status='replace')
    write(2090,'(5A20)') 'time','sink ID','r_hill','r_bondihoyle','r_tidalneigh'
    close(2090)
 endif 
 if (print_j_range) then 
    open(2100,file='accretable_j_range.ev',status='replace')
    write(2100,'(9A20)') 'time','sink ID','jx_sink','jy_sink','jz_sink','jacc_min', &
                       & 'jacc_max','jgas_min','jgas_max'
    close(2100)
 endif 
 if (print_v_rel) then 
    open(2200,file='vrel_sinkgas.ev',status='replace')
    write(2200,'(5A20)') 'time','sink ID','vtan_sink','vtan_gas','vrel_sinkgas'
    close(2200)
 endif 

end subroutine init


end module initial