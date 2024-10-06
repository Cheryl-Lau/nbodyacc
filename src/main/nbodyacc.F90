
program nbodyacc 

 use initial,  only:init
 use evolve,   only:evol
 use ptmass,   only:nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass
#ifdef BINARY
 use ptmass,   only:massq_ptmass
#endif 
 use ptmass,   only:allocate_ptmass,deallocate_ptmass
 use timestep, only:t_init
 use readwrite_infile, only:read_infile
 use readwrite_dump,   only:get_first_dump,read_dump,restart_evfile

 implicit none
 
 character(len=16) :: starting_dump

#ifdef BINARY   
 print*,'binary on from main'
#else 
 print*,'binary off from main'
#endif 

 !- Initial checks
 call allocate_ptmass
 call restart_evfile

 call read_infile

 call get_first_dump(starting_dump)

#ifdef BINARY
 call read_dump(starting_dump,t_init,nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)
 call init(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,massq_ptmass)
 !- Evolve the simulation 
 call evol(t_init,nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,massq_ptmass)
#else 

 call read_dump(starting_dump,t_init,nptmass,xyzhm_ptmass,vxyz_ptmass)
 call init(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass)
 !- Evolve the simulation 
 call evol(t_init,nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass)
#endif 

 !- End 
 call deallocate_ptmass


end program 