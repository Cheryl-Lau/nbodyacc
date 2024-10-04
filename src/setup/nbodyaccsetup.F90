
program nbodyaccsetup 

 use setup,  only:set_ptmass,read_setup_infile,write_dump,write_infile
 use ptmass, only:maxptmass,nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass
 use ptmass, only:allocate_ptmass,deallocate_ptmass
 implicit none

 call read_setup_infile()

 call allocate_ptmass()

 call set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass)
 call write_dump()    ! the starting dump 

 call deallocate_ptmass()

 call write_infile()  ! with default settings 

end 
