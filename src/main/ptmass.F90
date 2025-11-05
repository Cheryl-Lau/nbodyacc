
module ptmass


 implicit none 
 public :: update_sep,accrete_gas,get_accretion_rad
 public :: compute_Lxyz,compute_Lspin 
 public :: allocate_ptmass,deallocate_ptmass 

 integer, public :: maxptmass = 1e3
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

 private

contains 


subroutine compute_Lxyz(nptmass,xyzhm_ptmass,vxyz_ptmass,Lxyz_tot,jxyz_tot)
 integer, intent(in)  :: nptmass 
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(in)  :: vxyz_ptmass(:,:)
 real,    intent(out) :: Lxyz_tot(3),jxyz_tot(3)   ! summed over all ptmass 
 integer :: i
 real    :: m,r(3),v(3),r_cross_v(3)

 Lxyz_tot = 0.d0 
 jxyz_tot = 0.d0 

 do i = 1,nptmass
    r = xyzhm_ptmass(1:3,i)
    v = vxyz_ptmass(1:3,i)
    m = xyzhm_ptmass(5,i)
    r_cross_v = cross_product(r,v)

    Lxyz_ptmass(1:3,i) = m*r_cross_v
    Lxyz_tot = Lxyz_tot + Lxyz_ptmass(1:3,i)

    jxyz_ptmass(1:3,i) = r_cross_v
    jxyz_tot = jxyz_tot + jxyz_ptmass(1:3,i)
 enddo 

end subroutine compute_Lxyz


subroutine compute_Lspin(nptmass,xyzhm_ptmass,sq_ptmass,Lspin_tot,jspin_tot)
 integer, intent(in)  :: nptmass 
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(in)  :: sq_ptmass(:,:)
 real,    intent(out) :: Lspin_tot(3),jspin_tot(3)
 integer :: i
 real    :: mass,sep,q,Lspinz,jspinz 

 Lspin_tot = 0.d0 
 jspin_tot = 0.d0 

 do i = 1,nptmass
     mass = xyzhm_ptmass(5,i)
     sep = sq_ptmass(1,i)
     q = sq_ptmass(2,i)

     Lspinz = 2.5d-1*sqrt(mass**3*sep) * q/(1.d0+q)**2
     Lspin_ptmass(1:3,i) = (/ 0.d0, 0.d0, Lspinz /)
     Lspin_tot = Lspin_tot + Lspin_ptmass(1:3,i)

     jspinz = Lspinz/mass 
     jspin_ptmass(1:3,i) = (/ 0.d0, 0.d0, jspinz /)
     jspin_tot = jspin_tot + jspin_ptmass(1:3,i)
 enddo 

end subroutine compute_Lspin


function cross_product(vec1,vec2)
 real, intent(in)    :: vec1(3),vec2(3)
 real, dimension(3)  :: cross_product
 cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
 cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
 cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end function cross_product



subroutine accrete_gas()

end subroutine accrete_gas


subroutine get_accretion_rad()


end subroutine get_accretion_rad


subroutine get_tidal_rad()


end subroutine get_tidal_rad


subroutine get_bondi_hyole()

end subroutine get_bondi_hyole


subroutine update_sep(nptmass,sq_ptmass)
 integer, intent(in)    :: nptmass 
 real,    intent(inout) :: sq_ptmass(:,:)

 

end subroutine update_sep




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

end subroutine deallocate_ptmass
 
end module ptmass