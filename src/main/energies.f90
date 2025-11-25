
module energy 

 implicit none 
 public :: get_energies,get_angmom

 private

contains

!
! Compute total energy of all particles 
!
subroutine get_energies(nptmass,xyzhm_ptmass,vxyz_ptmass,ekin,epot,etot)
 use ptmass, only:poten_ptmass  ! filled after calling forces 
 integer, intent(in)  :: nptmass 
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(in)  :: vxyz_ptmass(:,:)
 real,    intent(out) :: ekin,epot,etot 
 integer :: i
 real    :: mi,vx,vy,vz,v2,phi

 !--Kinetic energies 
 ekin = 0.d0
 do i = 1,nptmass 
    mi = xyzhm_ptmass(5,i)
    vx = vxyz_ptmass(1,i)
    vy = vxyz_ptmass(2,i)
    vz = vxyz_ptmass(3,i)
    v2 = vx**2 + vy**2 + vz**2 
    ekin = ekin + 5.d-1*mi*v2 
 enddo 

 !--Potential energies 
 epot = 0.d0 
 do i = 1,nptmass 
    mi  = xyzhm_ptmass(5,i)
    phi = poten_ptmass(i)
    epot = epot + mi*phi 
 enddo 

 !--Total energy 
 etot = ekin + epot 

end subroutine get_energies

!
! Compute total angular momentum of all particles 
!
subroutine get_angmom(nptmass,xyzhm_ptmass,vxyz_ptmass,jxyz,jtot,sq_ptmass)
 use ptmass, only:compute_Lxyz
#ifdef BINARY
 use ptmass, only:compute_Lspin
#endif 
 integer, intent(in)  :: nptmass 
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(in)  :: vxyz_ptmass(:,:)
 real,    intent(out) :: jxyz(3),jtot(3)
 real,    intent(in), optional :: sq_ptmass(:,:)
 integer :: i
 real    :: xi,yi,zi,mi,vxi,vyi,vzi,si,qi
 real    :: jxyzi(3),Lxyzi(3),jspini(3),Lspini(3),jspin(3)

 jxyz  = 0.d0 
 jspin = 0.d0 

 do i = 1,nptmass 
    xi  = xyzhm_ptmass(1,i)
    yi  = xyzhm_ptmass(2,i)
    zi  = xyzhm_ptmass(3,i) 
    mi  = xyzhm_ptmass(5,i)
    vxi = vxyz_ptmass(1,i)
    vyi = vxyz_ptmass(2,i)
    vzi = vxyz_ptmass(3,i)

    !--Orbital angular momentum 
    call compute_Lxyz(i,xi,yi,zi,mi,vxi,vyi,vzi,jxyzi,Lxyzi)
    jxyz = jxyz + jxyzi

    !--Spin angular momentum 
#ifdef BINARY
     si = sq_ptmass(1,i)
     qi = sq_ptmass(2,i)
     call compute_Lspin(i,mi,si,qi,jspini,Lspini)
     jspin = jspin + jspini 
#endif 
 enddo 

 jtot = jxyz + jspin 

end subroutine get_angmom

end module energy 