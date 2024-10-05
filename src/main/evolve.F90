
module evolve 

 use ptmass,   only:nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass
 use accrete,  only:accrete_gas
 use step_RK4, only:step
 use timestep, only:t_end
 implicit none 
 public :: evol

 private

contains

subroutine evol()
 real :: t

 evol_loop: do while (t <= t_end)

    call step(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass) 

    call accrete_gas()

    call write_dump()

 enddo evol_loop

end subroutine evol 


subroutine write_dump()

end subroutine write_dump

end module evolve 
