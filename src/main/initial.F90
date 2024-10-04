
module initial

 use ptmass, only:nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,allocate_ptmass
 implicit none 
 public :: init

 private

contains 

subroutine init()

 call read_infile
 call read_startdump
 call allocate_ptmass



end subroutine init

subroutine read_infile()
 ! write with defaults if does not exist 

end subroutine read_infile

subroutine read_startdump()


end subroutine read_startdump

end module initial