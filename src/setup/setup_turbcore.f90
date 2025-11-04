
module setup 

 implicit none 
 public :: set_ptmass
 private

 integer :: iseed = -4321
 real    :: umass,udist 

contains 

!-------------------------------------------------------------------
! Set up the initial positions and velocities of point masses
!-------------------------------------------------------------------
subroutine set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass)
 use physcon,      only:pc,solarm
 use units,        only:set_units,utime,unit_velocity
 use centreofmass, only:transform_com_frame
 use turb_grid,    only:map_turbvel
 integer, intent(out)   :: nptmass
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 integer :: i
 real    :: Rcore_pc,Rcore,Mclust_solarm,Mclust,angvel_cgs,angvel,r_acc_pc,r_acc
 real    :: rms_mach,boxsize,cs_cgs,cs,x,y,z,mass 

 !--Define unit system in sim
 udist = pc
 umass = solarm
 call set_units(dist=udist,mass=umass,G=1.d0)

 !--Inputs 
 nptmass       = 30        ! Number of stars 
 Rcore_pc      = 0.2       ! Core radius [pc]
 Mclust_solarm = 1.d3      ! Cluster total mass [msun]
 angvel_cgs    = 3.d-14    ! Rotation angular vxyzocity [rad/s]
 r_acc_pc      = 1.d-5     ! Accretion radius [pc] (to be re-evaluated during runtime) 
 rms_mach      = 1.        ! Turbulence Mach number 
 cs_cgs        = 2.19d4    ! Sound speed in cm/s

 !--Convert to code units 
 Rcore  = Rcore_pc*pc/udist 
 Mclust = Mclust_solarm*solarm/umass
 angvel = angvel_cgs/unit_velocity
 r_acc  = r_acc_pc*pc/udist 
 cs     = cs_cgs/unit_velocity 

 !--Set particle position, accretion radius, and mass
 do i = 1,nptmass
    call gen_random_pos(Rcore,x,y,z)
    xyzhm_ptmass(1:3,i) = (/ x, y, z /)
    xyzhm_ptmass(4,i) = r_acc
    call gen_random_mass(2.d0,5.d0,mass)
    xyzhm_ptmass(5,i) = mass
 enddo 

 !--Adding turbulence
 boxsize = Rcore 
 call map_turbvel(boxsize,nptmass,xyzhm_ptmass,vxyz_ptmass,cs,rms_mach)

 !--Adding rotation 
 do i = 1,nptmass
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) - angvel*y
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + angvel*x
 enddo 

 !--Reset COM 
 call transform_com_frame(nptmass,xyzhm_ptmass,vxyz_ptmass)

end subroutine set_ptmass

!
! Routine to randomly generarte xyzitions 
!
subroutine gen_random_pos(rmax,x,y,z)
 use random, only:ran2
 real,   intent(in)  :: rmax 
 real,   intent(out) :: x,y,z
 real    :: r,rmin

 rmin = 0.2*rmax  ! not too close to origin

 r = rmax + 1  ! dummy 
 do while (r > rmax .or. r < rmin)
    x = 2.*(ran2(iseed)-0.5) *rmax 
    y = 2.*(ran2(iseed)-0.5) *rmax
    z = 2.*(ran2(iseed)-0.5) *rmax
    r = sqrt(x**2 + y**2 + z**2) 
 enddo 

end subroutine gen_random_pos 

!
! Routine to randomly generarte masses
!
subroutine gen_random_mass(lowbound,uppbound,mass)
 use random, only:ran2
 real,   intent(in)  :: lowbound,uppbound 
 real,   intent(out) :: mass 

 mass = lowbound + ran2(iseed)*(uppbound-lowbound)

end subroutine gen_random_mass 


end module setup 