
module timestep 

 implicit none 
 public :: constrain_dt 
 public :: read_infile_timestep,write_infile_timestep
 
 real,  public :: dtmax = 1e-3
 real,  public :: t_end = 1e2
 real,  public :: dt
 
 private
 namelist /timestep_params/ dtmax,t_end

contains 

subroutine constrain_dt()

end subroutine constrain_dt




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