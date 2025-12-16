
module ptmass


 implicit none 
 public :: compute_Lxyz
#ifdef BINARY
 public :: update_sq,compute_Lspin
#endif 
 public :: allocate_ptmass,deallocate_ptmass 

 integer, public :: maxptmass = 1e2
 integer, public :: nptmass 
 real,    public, allocatable :: xyzhm_ptmass(:,:) ! position, accretion radius, mass 
 real,    public, allocatable :: vxyz_ptmass(:,:)  ! velocity 
 real,    public, allocatable :: fxyz_ptmass(:,:)  ! forces   
 real,    public, allocatable :: poten_ptmass(:)   ! potentials 
 real,    public, allocatable :: Lxyz_ptmass(:,:)  ! angular momentum in sim frame
 real,    public, allocatable :: jxyz_ptmass(:,:)  ! specific angular momentum in sim frame 
#ifdef BINARY
 real,    public, allocatable :: sq_ptmass(:,:)    ! separation and mass-ratio of binary pair 
 real,    public, allocatable :: Lspin_ptmass(:,:) ! angular momentum around COM of binary 
 real,    public, allocatable :: jspin_ptmass(:,:) ! specific angular momentum around COM of binary
#endif 
 real,    public, allocatable :: racc_ptmass(:,:)  ! individual r_acc terms 

 private

contains 

!
! Compute (specific) orbital angular momentum for a given particle 
!
subroutine compute_Lxyz(ip,xi,yi,zi,mi,vxi,vyi,vzi,jxyz,Lxyz)
 integer, intent(in)  :: ip
 real,    intent(in)  :: xi,yi,zi,mi,vxi,vyi,vzi
 real,    intent(out) :: jxyz(3),Lxyz(3)
 real     :: ri(3),vi(3)

 ri = (/ xi,  yi,  zi  /)
 vi = (/ vxi, vyi, vzi /)
 jxyz = cross_product(ri,vi)
 Lxyz = mi*cross_product(ri,vi)

 jxyz_ptmass(1:3,ip) = jxyz 
 Lxyz_ptmass(1:3,ip) = Lxyz 

end subroutine compute_Lxyz 

!
! Compute (specific) spin angular momentum for a given particle 
! assuming that it only rotates about the z-axis 
!
#ifdef BINARY
subroutine compute_Lspin(ip,mi,si,qi,jspin,Lspin)
 integer, intent(in)  :: ip
 real,    intent(in)  :: mi,si,qi
 real,    intent(out) :: jspin(3),Lspin(3)
 real     :: jspinz,Lspinz 

 Lspinz = sqrt(mi**3*si) * qi/(1.d0+qi)**2
 jspinz = Lspinz/mi
 jspin  = (/ 0.d0, 0.d0, jspinz /)
 Lspin  = (/ 0.d0, 0.d0, Lspinz /)

 jspin_ptmass(1:3,ip) = jspin 
 Lspin_ptmass(1:3,ip) = Lspin 

end subroutine compute_Lspin 
#endif 

!
! Compute binary separation with given Lspin of a particle 
!
#ifdef BINARY
subroutine update_sq(mi,qi,Lspinz,si)
 real, intent(in)  :: mi,qi,Lspinz 
 real, intent(out) :: si

 if (qi < tiny(qi)) return      ! single star

 si = Lspinz**2 * qi**(-2) * (1.d0+qi)**4 * mi**(-3)

end subroutine update_sq
#endif 


!
! Math tools 
!
function cross_product(vec1,vec2)
 real, intent(in)    :: vec1(3),vec2(3)
 real, dimension(3)  :: cross_product
 cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
 cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
 cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end function cross_product



!
! Allocate memory 
!
subroutine allocate_ptmass

 allocate(xyzhm_ptmass(5,maxptmass))
 allocate(vxyz_ptmass(3,maxptmass))
 allocate(fxyz_ptmass(3,maxptmass))
 allocate(poten_ptmass(maxptmass))
 allocate(Lxyz_ptmass(3,maxptmass))
 allocate(jxyz_ptmass(3,maxptmass))
#ifdef BINARY
 allocate(sq_ptmass(2,maxptmass))
 allocate(Lspin_ptmass(3,maxptmass))
 allocate(jspin_ptmass(3,maxptmass))
#endif 
 allocate(racc_ptmass(3,maxptmass))

end subroutine allocate_ptmass


subroutine deallocate_ptmass

 deallocate(xyzhm_ptmass)
 deallocate(vxyz_ptmass)
 deallocate(fxyz_ptmass)
 deallocate(poten_ptmass)
 deallocate(Lxyz_ptmass)
 deallocate(jxyz_ptmass)
#ifdef BINARY
 deallocate(sq_ptmass)
 deallocate(Lspin_ptmass)
 deallocate(jspin_ptmass)
#endif 
 deallocate(racc_ptmass)

end subroutine deallocate_ptmass
 
end module ptmass