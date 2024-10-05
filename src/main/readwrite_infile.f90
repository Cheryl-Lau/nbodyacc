
module readwrite_infile 

 implicit none 
 public :: write_infile,read_infile

 private

 integer :: unit_infile = 2025

contains 

subroutine write_infile()
 use timestep, only:write_infile_timestep
 use force,    only:write_infile_force
 integer :: rc 

 open(unit_infile,file='input_params.in',iostat=rc) 

 ! Write options from each physics module 
 call write_infile_timestep(unit_infile)
 call write_infile_force(unit_infile)

 close(unit_infile)

end subroutine write_infile


subroutine read_infile()
 use timestep, only:read_infile_timestep
 use force,    only:read_infile_force
 integer :: rc 

 inquire (file='runtime_params.in',iostat=rc)
 if (rc /= 0) stop 'no input files found'
 open(unit_infile,file='input_params.in',iostat=rc) 

 ! Read options from each physics module 
 call read_infile_timestep(unit_infile)
 call read_infile_force(unit_infile)

 close(unit_infile)

end subroutine read_infile


end module readwrite_infile 