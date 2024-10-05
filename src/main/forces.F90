
module force 

 implicit none 
 public :: compute_forces
 public :: read_infile_force,write_infile_force

 logical, public :: use_plumerpot = .true. 

 private 
 namelist /force_params/ use_plumerpot

contains 

!
! Total force per unit mass (actually acceleration)
!
subroutine compute_forces(nptmass,xyzhm,fxyz)
 integer, intent(in)  :: nptmass
 real,    intent(in)  :: xyzhm(5,nptmass)
 real,    intent(out) :: fxyz(3,nptmass)

 call plumer_potential()
#ifdef GRAVITY
 call self_grav()
#endif 

end subroutine compute_forces


subroutine self_grav()

end subroutine self_grav


subroutine plumer_potential()

end subroutine plumer_potential


subroutine galdisc_potential()

end subroutine galdisc_potential





subroutine read_infile_force(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=force_params,iostat=rc)
 if (rc /= 0) stop 'cannot read force options'

end subroutine read_infile_force


subroutine write_infile_force(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=force_params,iostat=rc)
 if (rc /= 0) stop 'cannot write timestep options'

end subroutine write_infile_force


end module force 