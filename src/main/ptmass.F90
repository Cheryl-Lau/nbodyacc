
module ptmass


 implicit none 
 public :: allocate_ptmass,deallocate_ptmass 

 integer, public :: maxptmass = 1e3
 integer, public :: nptmass 
 real,    public, allocatable :: xyzhm_ptmass(:,:) ! position, accretion radius, mass 
 real,    public, allocatable :: vxyz_ptmass(:,:)  ! velocity 
 real,    public, allocatable :: fxyz_ptmass(:,:)  ! forces   
 real,    public, allocatable :: Lxyz_ptmass(:,:)  ! total angular momentum in sim frame
#ifdef BINARY
 real,    public, allocatable :: sep_ptmass(:)
 real,    public, allocatable :: massq_ptmass(:)
 real,    public, allocatable :: Lspin_ptmass(:,:) ! angular momentum around COM of binary pair 
#endif 

 private

contains 



subroutine allocate_ptmass

 allocate(xyzhm_ptmass(5,maxptmass))
 allocate(vxyz_ptmass(3,maxptmass))
 allocate(fxyz_ptmass(3,maxptmass))
 allocate(Lxyz_ptmass(3,maxptmass))
#ifdef BINARY
 allocate(sep_ptmass(maxptmass))
 allocate(massq_ptmass(maxptmass))
 allocate(Lspin_ptmass(maxptmass))
#endif 

end subroutine allocate_ptmass
 


subroutine deallocate_ptmass

 deallocate(xyzhm_ptmass)
 deallocate(vxyz_ptmass)
 deallocate(fxyz_ptmass)
 deallocate(Lxyz_ptmass)
#ifdef BINARY
 deallocate(sep_ptmass)
 deallocate(massq_ptmass)
 deallocate(Lspin_ptmass)
#endif 

end subroutine deallocate_ptmass
 
end module ptmass