
module timestep 

 implicit none 
 public :: constrain_dt 
 public :: read_infile_timestep,write_infile_timestep
 
 integer, public :: nout = 100    ! write dump every <nout> dt
 real,    public :: dtmax  = 1d-3
 real,    public :: t_end  = 1d2
 real,    public :: t_init = 0.d0
 
 private
 namelist /timestep_params/ dtmax,t_end,nout

contains 

subroutine constrain_dt(nptmass,xyzhm_ptmass,vxyz_ptmass,dt)
 integer, intent(in) :: nptmass
 real,    intent(in) :: xyzhm_ptmass(:,:)
 real,    intent(in) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: dt

end subroutine constrain_dt


!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_timestep(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=timestep_params,iostat=rc)
 if (rc /= 0) stop 'cannot read timestep options'

end subroutine read_infile_timestep


subroutine write_infile_timestep(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=timestep_params,iostat=rc)
 if (rc /= 0) stop 'cannot write timestep options'

end subroutine write_infile_timestep

end module timestep 