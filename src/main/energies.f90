
module energy 

 implicit none 
 public :: get_energies,get_angmom

 private

contains

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


subroutine get_angmom(nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass,jxyz,jspin,jtot)
 use ptmass, only:compute_Lxyz,compute_Lspin 
 integer, intent(in)  :: nptmass 
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(in)  :: vxyz_ptmass(:,:)
 real,    intent(in),  optional :: sq_ptmass(:,:)
 real,    intent(out), optional :: jspin(3)
 real,    intent(out) :: jxyz(3),jtot(3)
 real    :: Lxyz_tot(3),Lspin_tot(3)

 call compute_Lxyz(nptmass,xyzhm_ptmass,vxyz_ptmass,Lxyz_tot,jxyz)

 binary: if (present(sq_ptmass)) then 
    call compute_Lspin(nptmass,xyzhm_ptmass,sq_ptmass,Lspin_tot,jspin)
    jtot = jxyz + jspin
 endif binary 

end subroutine get_angmom

end module energy 