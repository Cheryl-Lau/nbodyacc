
module evolve 

 use readwrite_dump, only:write_dump,write_evfile
 implicit none 
 public :: evol

 private

contains

subroutine evol(t_init,nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)
 use step_RK4, only:step
 use timestep, only:dtmax,t_end,nout,constrain_dt
 use energy,   only:get_energies,get_angmom
 use ptmass,   only:get_accretion_rad,accrete_gas
#ifdef BINARY
 use ptmass,   only:Lspin_ptmass,update_sep
#endif 
 real,    intent(in)    :: t_init
 integer, intent(inout) :: nptmass 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout), optional :: sq_ptmass(:,:)
 integer :: istep,iout
 real    :: t,dt,ekin,epot,etot,jspin(3),jxyz(3),jtot(3)


 if (t_init > t_end) stop 't_end needs to be greater than t_init'
 t  = t_init
 dt = dtmax
 iout = 0
 istep = 0

 evol_loop: do while (t <= t_end)

    !--RK4 integrator 
    call step(nptmass,xyzhm_ptmass,vxyz_ptmass,dt) 

    !--Accrete and update particles 
    call get_accretion_rad()
    call accrete_gas()
#ifdef BINARY
    call update_sep(nptmass,sq_ptmass)
#endif 

    !--Compute energies and specific angular momentum 
    call get_energies(nptmass,xyzhm_ptmass,vxyz_ptmass,ekin,epot,etot)
#ifdef BINARY
    call get_angmom(nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass,jxyz,jspin,jtot)
#else 
    call get_angmom(nptmass,xyzhm_ptmass,vxyz_ptmass,jxyz,jtot)
#endif 

    !--Write dumps
    if (iout == nout) then 
#ifdef BINARY
       call write_dump(t,nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)
#else 
       call write_dump(t,nptmass,xyzhm_ptmass,vxyz_ptmass)
#endif 
       iout = 0
    endif 

    !--Write evol of energy/angmom 
#ifdef BINARY
    call write_evfile(t,ekin,epot,etot,jxyz,jspin,jtot)
#else 
    call write_evfile(t,ekin,epot,etot,jxyz)
#endif 

    !--Control timestep for next iteration 
    call constrain_dt(nptmass,xyzhm_ptmass,vxyz_ptmass,dt)
    t = t + dt

    iout = iout + 1
    istep = istep + 1 
    if (mod(istep,nint(t_end/dtmax)/10) == 0) print*,nint(t/t_end*100.),'% done'

 enddo evol_loop

end subroutine evol 


end module evolve 
