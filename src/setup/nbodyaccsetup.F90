
program nbodyaccsetup 

 use setup,  only:set_ptmass,write_dump
 use ptmass, only:nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass
 use ptmass, only:allocate_ptmass,deallocate_ptmass
 use readwrite_infile, only:write_infile 
 implicit none

 call allocate_ptmass

#ifdef BINARY
 call set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)
 call write_dump(nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)    ! the starting dump 
#else 
 call set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass)
 call write_dump(nptmass,xyzhm_ptmass,vxyz_ptmass)
#endif 

 call deallocate_ptmass

 call write_infile ! with default settings 

end program
