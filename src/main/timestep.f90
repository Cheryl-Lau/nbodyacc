
module timestep 

 implicit none 
 public :: constrain_dt 
 public :: read_infile_timestep,write_infile_timestep
 
 integer, public :: nout    = 10      ! write dump every <nout> dtmax 
 integer, public :: maxdump = 300     ! max number of dumpfiles 
 real,    public :: dtmax   = 1.d-5   ! sim max timestep 
 real,    public :: t_init  = 0.d0    ! sim start time 
 real,    public :: t_end   = 2.d-2   ! sim end time 

 integer, public :: nbinmax  = 20
 real,    public :: alpha_dt = 1.d-4    ! timestep constraint param 
 

 private

 namelist /timestep_params/ dtmax,t_end,nout,maxdump,alpha_dt

contains 

subroutine constrain_dt(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,nbin,dt)
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: xyzhm_ptmass(:,:)
 real,    intent(in)    :: vxyz_ptmass(:,:)
 real,    intent(in)    :: fxyz_ptmass(:,:)
 real,    intent(inout) :: dt
 integer, intent(out)   :: nbin
 integer :: i,j
 real    :: ri(3),rj(3),vi(3),vj(3),fi(3),fj(3),r_over_v,v_over_f
 real    :: dtfrac
 
 dt = huge(dt)
 do i = 1,nptmass
    ri = xyzhm_ptmass(1:3,i)
    vi = vxyz_ptmass(1:3,i)
    fi = fxyz_ptmass(1:3,i)
    do j = 1,nptmass
        rj = xyzhm_ptmass(1:3,j)
        vj = vxyz_ptmass(1:3,j)
        fj = fxyz_ptmass(1:3,j)
        r_over_v = mag(ri-rj)/mag(vi-vj)
        v_over_f = mag(vi-vj)/mag(fi-fj)
        dt = min(dt,r_over_v,v_over_f)
    enddo 
 enddo 
 dt = alpha_dt * dt

 !--Assign timestep bin 
 dtfrac = dtmax/dt
 nbin = int(log(dtfrac)/log(2.d0))
 nbin = min(nbin,nbinmax)
 nbin = max(nbin,0)

 dt  = dtmax/(2**nbin)

end subroutine constrain_dt


real function mag(vec)
 real, intent(in) :: vec(3)
 mag = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
end function mag


!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_timestep(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=timestep_params,iostat=rc)
 if (rc /= 0) stop 'cannot read timestep options'

end subroutine read_infile_timestep


subroutine write_infile_timestep(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=timestep_params,iostat=rc)
 if (rc /= 0) stop 'cannot write timestep options'

end subroutine write_infile_timestep

end module timestep 