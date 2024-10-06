
module evolve 

 use readwrite_dump, only:write_dump,write_evfile
 implicit none 
 public :: evol

 private

contains

subroutine evol(t_init,nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,massq_ptmass)
 use accrete,  only:get_accretion_rad,accrete_gas
 use step_RK4, only:step
 use timestep, only:dtmax,t_end,nout,constrain_dt
 use ptmass,   only:Lxyz_ptmass
 use energy,   only:get_tot_angmomen,get_tot_energy
#ifdef BINARY
 use ptmass,   only:sep_ptmass,Lspin_ptmass,update_sep
#endif 
 real,    intent(in)    :: t_init
 integer, intent(inout) :: nptmass 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 real,    intent(inout), optional :: massq_ptmass(:)
 integer :: istep,iout
 real    :: t,dt,ekin,epot 

 if (t_init > t_end) stop 't_end needs to be greater than t_init'
 t  = t_init
 dt = dtmax
 iout = 0
 istep = 0

 evol_loop: do while (t <= t_end)

    call step(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,dt) 

    call get_tot_energy(nptmass,xyzhm_ptmass,vxyz_ptmass,ekin,epot)
    call get_tot_angmomen(nptmass,xyzhm_ptmass,vxyz_ptmass,Lxyz_ptmass)

    call get_accretion_rad()
    call accrete_gas()

#ifdef BINARY
    call update_sep(nptmass,massq_ptmass)
#endif 

    !- Write outputs 
    if (iout == nout) then 
#ifdef BINARY
       call write_dump(t,nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)
#else 
       call write_dump(t,nptmass,xyzhm_ptmass,vxyz_ptmass)
#endif 
       iout = 0
    endif 
    call write_evfile(t,ekin,epot)

    !- Control timestep for next iteration 
    call constrain_dt(nptmass,xyzhm_ptmass,vxyz_ptmass,dt)
    t = t + dt


    iout = iout + 1
    istep = istep + 1 
    if (mod(istep,nint(t_end/dtmax)/10) == 0) print*,'We are ',nint(t/t_end*100.),'% done'

 enddo evol_loop

end subroutine evol 


end module evolve 
