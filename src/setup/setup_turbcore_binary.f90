
module setup 

 implicit none 
 public :: set_ptmass

 private

 real :: umass,udist 

contains 

!-------------------------------------------------------------------
! Set up the initial positions and velocities of point masses
!-------------------------------------------------------------------
subroutine set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)
 use physcon,   only:pc,solarm
 use units,     only:set_units,utime,unit_velocity
 use turb_grid, only:map_turbvel,rms_mach,c_sound
 integer, intent(out)   :: nptmass
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: massq_ptmass(:)
 integer :: ip,idim
 real    :: box_size,mass,r_acc,ran_pos

 !- Define unit system in sim
 udist = pc
 umass = solarm
 call set_units(dist=udist,mass=umass,G=1.d0)

 !- Set positions and masses 
 nptmass  = 10 
 box_size = 80.d0 
 r_acc    = 2.d-1  ! to be re-evaluated during runtime

 do ip = 1,nptmass
    do idim = 1,3
       call random_number(ran_pos)
       xyzhm_ptmass(idim,ip) = box_size/2.d0 * (2.d0*ran_pos-1.d0)
    enddo 
    xyzhm_ptmass(4,ip) = r_acc
    xyzhm_ptmass(5,ip) = 5.d0
    massq_ptmass(ip)   = 1.     ! equal mass binaries 
 enddo 

 !- Adding turbulence
 rms_mach = 1. 
 c_sound = 2.
 call map_turbvel(box_size,nptmass,xyzhm_ptmass,vxyz_ptmass)

end subroutine set_ptmass

end module setup 