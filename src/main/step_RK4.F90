
module step_RK4

 implicit none 
 public :: step 

 private 

contains 

!-------------------------------------------------------------------
! Fourth-order Runge-Kutta integrator 
!-------------------------------------------------------------------
subroutine step(nptmass,xyzhm,vxyz,fxyz)
 use timestep, only:dt
 use force,    only:compute_forces
 integer, intent(in)    :: nptmass
 real,    intent(inout) :: xyzhm(:,:) ! pos 
 real,    intent(inout) :: vxyz(:,:)  ! vel
 real,    intent(inout) :: fxyz(:,:)  ! accel
 integer :: ip
 real    :: hdt 
 !- dummy vars for RK4 intermediate steps 
 real    :: xyzhm1(5,nptmass),xyzhm2(5,nptmass),xyzhm3(5,nptmass)
 real    :: vxyz1(3,nptmass),vxyz2(3,nptmass),vxyz3(3,nptmass)
 real    :: fxyz0(3,nptmass),fxyz1(3,nptmass),fxyz2(3,nptmass),fxyz3(3,nptmass)

 hdt = dt/2.d0  ! half-step

 call compute_forces(nptmass,xyzhm,fxyz0)
 do ip = 1,nptmass 
    vxyz1(:,ip) = vxyz(:,ip) + fxyz0(:,ip)*hdt
    xyzhm1(1:3,ip) = xyzhm(1:3,ip) + vxyz1(:,ip)*hdt
 enddo 

 call compute_forces(nptmass,xyzhm1,fxyz1)
 do ip = 1,nptmass 
    vxyz2(:,ip) = vxyz(:,ip) + fxyz1(:,ip)*hdt 
    xyzhm2(1:3,ip) = xyzhm(1:3,ip) + vxyz2(:,ip)*hdt 
 enddo 

 call compute_forces(nptmass,xyzhm2,fxyz2)
 do ip = 1,nptmass 
    vxyz3(:,ip) = vxyz(:,ip) + fxyz2(:,ip)*hdt 
    xyzhm3(1:3,ip) = xyzhm(1:3,ip) + vxyz3(:,ip)*hdt 
 enddo 

 call compute_forces(nptmass,xyzhm3,fxyz3)

 !- Actual update 
 do ip = 1,nptmass
    vxyz(:,ip) = vxyz(:,ip) + 1.d0/6.d0 * (fxyz0(:,ip) + 2.d0*fxyz1(:,ip) + 2.d0*fxyz2(:,ip) + fxyz3(:,ip)) * dt 
    xyzhm(1:3,ip) = xyzhm(1:3,ip) + 1.d0/6.d0 * (vxyz(:,ip) + 2.d0*vxyz1(:,ip) + 2.d0*vxyz2(:,ip) + vxyz3(:,ip)) * dt
    fxyz(:,ip) = 1.d0/6.d0 * (fxyz0(:,ip) + 2.d0*fxyz1(:,ip) + 2.d0*fxyz2(:,ip) + fxyz3(:,ip))  ! store
 enddo 

end subroutine step

end module step_RK4