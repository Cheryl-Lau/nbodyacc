
module step_RK4

 implicit none 
 public :: step 

 private 

contains 

!-------------------------------------------------------------------
! Fourth-order Runge-Kutta integrator 
!-------------------------------------------------------------------
subroutine step(nptmass,xyzhm,vxyz,dt)
 use force,    only:compute_forces
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: dt
 real,    intent(inout) :: xyzhm(:,:) 
 real,    intent(inout) :: vxyz(:,:)  
 integer :: i
 real    :: hdt 
 real    :: xyzhm1(5,nptmass),xyzhm2(5,nptmass),xyzhm3(5,nptmass)   ! dummy vars for intermediate steps 
 real    :: vxyz1(3,nptmass),vxyz2(3,nptmass),vxyz3(3,nptmass)
 real    :: fxyz0(3,nptmass),fxyz1(3,nptmass),fxyz2(3,nptmass),fxyz3(3,nptmass)

 hdt = dt/2.d0  ! half-step

 call compute_forces(nptmass,xyzhm,fxyz0)
 do i = 1,nptmass 
    vxyz1(:,i) = vxyz(:,i) + fxyz0(:,i)*hdt
    xyzhm1(1:3,i) = xyzhm(1:3,i) + vxyz1(:,i)*hdt
 enddo 

 call compute_forces(nptmass,xyzhm1,fxyz1)
 do i = 1,nptmass 
    vxyz2(:,i) = vxyz(:,i) + fxyz1(:,i)*hdt 
    xyzhm2(1:3,i) = xyzhm(1:3,i) + vxyz2(:,i)*hdt 
 enddo 

 call compute_forces(nptmass,xyzhm2,fxyz2)
 do i = 1,nptmass 
    vxyz3(:,i) = vxyz(:,i) + fxyz2(:,i)*hdt 
    xyzhm3(1:3,i) = xyzhm(1:3,i) + vxyz3(:,i)*hdt 
 enddo 

 call compute_forces(nptmass,xyzhm3,fxyz3)

 !- Actual update 
 do i = 1,nptmass
    vxyz(:,i) = vxyz(:,i) + 1.d0/6.d0 * (fxyz0(:,i) + 2.d0*fxyz1(:,i) + 2.d0*fxyz2(:,i) + fxyz3(:,i)) * dt 
    xyzhm(1:3,i) = xyzhm(1:3,i) + 1.d0/6.d0 * (vxyz(:,i) + 2.d0*vxyz1(:,i) + 2.d0*vxyz2(:,i) + vxyz3(:,i)) * dt
 enddo 

end subroutine step

end module step_RK4