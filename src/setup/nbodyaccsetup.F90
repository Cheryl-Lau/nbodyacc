
program nbodyaccsetup 

 use setup,  only:set_ptmass
 use units,  only:utime,udist,umass,write_units
 use ptmass, only:nptmass,xyzhm_ptmass,vxyz_ptmass
#ifdef BINARY
 use ptmass, only:sq_ptmass
#endif 
 use ptmass, only:allocate_ptmass,deallocate_ptmass
 use readwrite_infile, only:write_infile 
 use readwrite_dump,   only:write_first_dump
 implicit none


 call allocate_ptmass

#ifdef BINARY
 print*,'binary ON'
 call set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)
 call write_first_dump(0.d0,nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)    ! the starting dump 
#else 
 print*,'binary OFF'
 call set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass)
 call write_first_dump(0.d0,nptmass,xyzhm_ptmass,vxyz_ptmass)
#endif 

 call deallocate_ptmass

 call write_infile ! with default settings 

 call write_units(utime,udist,umass)  ! from setup 

end program
