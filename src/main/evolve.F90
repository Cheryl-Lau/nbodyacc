
module evolve 

 use readwrite_dump, only:write_dump,write_evfile
 implicit none 
 public :: evol

 private

contains

subroutine evol(t_init,nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,sq_ptmass)
 use step_RK4, only:step
 use timestep, only:dtmax,t_end,nout,maxdump,constrain_dt
 use energy,   only:get_energies,get_angmom
 use ptmass,   only:get_accretion_rad,accrete_gas
#ifdef BINARY
 use ptmass,   only:Lspin_ptmass,update_sep
#endif 
 real,    intent(in)    :: t_init
 integer, intent(inout) :: nptmass 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 real,    intent(inout), optional :: sq_ptmass(:,:)
 integer :: iout,ndump,nbin
 real    :: t,dt,t_substep,ekin,epot,etot,jspin(3),jxyz(3),jtot(3)


 if (t_init > t_end) stop 't_end needs to be greater than t_init'
 t     = t_init
 dt    = dtmax
 iout  = 0
 ndump = 0 

 evol_loop: do while (t <= t_end .and. ndump < maxdump)

    !--Control timestep 
    call constrain_dt(nptmass,vxyz_ptmass,fxyz_ptmass,nbin,dt)
    print*,'time = ',t,'; dt = ',dt,' ; nbin = ',nbin

    t_substep = 0.d0 
    substep: do while (t_substep <= dtmax+tiny(dtmax))

       !--RK4 integrator 
       call step(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,dt) 

       !--Accrete and update particles 
       call get_accretion_rad()
       call accrete_gas()
#ifdef BINARY
       call update_sep(nptmass,sq_ptmass)
#endif 

       t = t + dt
       t_substep = t_substep + dt 

    enddo substep 

    !--Compute energies and specific angular momentum 
    call get_energies(nptmass,xyzhm_ptmass,vxyz_ptmass,ekin,epot,etot)
#ifdef BINARY
    call get_angmom(nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass,jxyz,jspin,jtot)
#else 
    call get_angmom(nptmass,xyzhm_ptmass,vxyz_ptmass,jxyz,jtot)
#endif 

    !--Write total energy/angmom 
#ifdef BINARY
    call write_evfile(t,ekin,epot,etot,jxyz,jspin,jtot)
#else 
    call write_evfile(t,ekin,epot,etot,jxyz)
#endif 

    !--Write dump every <nout> dtmax 
    if (iout == nout) then 
#ifdef BINARY
       call write_dump(t,nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)
#else 
       call write_dump(t,nptmass,xyzhm_ptmass,vxyz_ptmass)
#endif 
       iout  = 0
       ndump = ndump + 1 
    endif
    iout = iout + 1

 enddo evol_loop

 print*,'Run completed'
 if (t < t_end)       print*,'Number of dumps reached ',ndump,'/',maxdump 
 if (ndump < maxdump) print*,'Time reached specified t_end: ',t,'/',t_end 

end subroutine evol 

end module evolve 
