
module turb_grid

 implicit none  
 public :: map_turbvel

 real,  public :: rms_mach = 3.0   ! Turbulent mach number
 real,  public :: c_sound  = 2.   ! Sound speed in molecular cloud

 private 

 contains 

!-------------------------------------------------------------------
! Reads the turbulence cubes generated with velfield.f and
! maps the velocities onto the given point mass locations.
! ***NOTE: all grids must be centred at the origin***
!-------------------------------------------------------------------
subroutine map_turbvel(box_size,nptmass,xyzhm_ptmass,vxyz_ptmass)
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: box_size
 real,    intent(in)    :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 integer :: nspace,nbytes_per_vel,nbytes_fheader
 integer :: ifilesizex,ifilesizey,ifilesizez
 integer :: iunit,ierr,ierrx,ierry,ierrz,iposx,iposy,iposz
 real(kind=4), allocatable :: velx(:,:,:),vely(:,:,:),velz(:,:,:)  ! turbulence cube 
 real,         allocatable :: xyz_grid(:,:,:)                      ! input grid cell positions 
 integer :: ip
 real    :: xi,yi,zi,deli,delx,dely,delz,radnorm,radius,factor
 real    :: turbboxsize,rmsmach,turbfac,v2i
 logical :: iexistx,iexisty,iexistz
 character(len=40), parameter :: filevx = 'cube_v1.dat'
 character(len=40), parameter :: filevy = 'cube_v2.dat'
 character(len=40), parameter :: filevz = 'cube_v3.dat'

 ierr = 0
 if (rms_mach <= 0.) stop 'Mach number must be a positive value'

 !- compute box size 
 turbboxsize = box_size/2.  !rmax(nptmass,xyzhm_ptmass)

 ! Check if the turb cubes exist 
 inquire(file=trim(filevx),exist=iexistx,size=ifilesizex)
 inquire(file=trim(filevy),exist=iexisty,size=ifilesizey)
 inquire(file=trim(filevz),exist=iexistz,size=ifilesizez)
 if (.not.iexistx .or. .not.iexisty .or. .not.iexistz) stop 'Please copy &
    &nbodyacc/data/turb_cubes/*.dat to your work-directory '

 !- work out the size of the velocity files from the cubes' filesize
 nbytes_per_vel = kind(velx)
 nbytes_fheader = 8   
 nspace = nint(((ifilesizex - nbytes_fheader)/nbytes_per_vel)**(1./3.))


 allocate(velx(nspace,nspace,nspace),stat=ierrx)
 allocate(vely(nspace,nspace,nspace),stat=ierry)
 allocate(velz(nspace,nspace,nspace),stat=ierrz)
 if (ierrx /= 0 .or. ierry /= 0 .or. ierrz /= 0) stop 'error allocating memory'

 !- read data cubes
 iunit = 2010
 call read_cube(trim(filevx),velx,nspace,iunit,ierrx)
 call read_cube(trim(filevy),vely,nspace,iunit,ierry)
 call read_cube(trim(filevz),velz,nspace,iunit,ierrz)
 if (ierrx /= 0 .or. ierry /= 0 .or. ierrz /= 0) stop

 !- map the velocity from cubes to ptmass locations 
 radnorm = turbboxsize
 deli = radnorm/real(nspace/2)
 do ip = 1,nptmass
    xi = xyzhm_ptmass(1,ip)
    yi = xyzhm_ptmass(2,ip)
    zi = xyzhm_ptmass(3,ip)

    iposx = int(xi/radnorm*(nspace/2)+(nspace/2)+0.5)
    iposx = min(max(iposx, 1),nspace-1)
    iposy = int(yi/radnorm*(nspace/2)+(nspace/2)+0.5)
    iposy = min(max(iposy, 1),nspace-1)
    iposz = int(zi/radnorm*(nspace/2)+(nspace/2)+0.5)
    iposz = min(max(iposz, 1),nspace-1)
    delx = xi - (iposx-(nspace/2)-0.5)/real(nspace/2)*radnorm
    dely = yi - (iposy-(nspace/2)-0.5)/real(nspace/2)*radnorm
    delz = zi - (iposz-(nspace/2)-0.5)/real(nspace/2)*radnorm
    radius = sqrt(xi**2 + yi**2 + zi**2)
    !- Find interpolated velocities
    vxyz_ptmass(1,ip) = vfield_interp(velx,iposx,iposy,iposz,delx,dely,delz,deli)
    vxyz_ptmass(2,ip) = vfield_interp(vely,iposx,iposy,iposz,delx,dely,delz,deli)
    vxyz_ptmass(3,ip) = vfield_interp(velz,iposx,iposy,iposz,delx,dely,delz,deli)
 enddo

 !- scale all velocities to the desired mach number 
 print*,'Setting turbulence to Mach number',rms_mach
 rmsmach = 0.0
 do ip = 1,nptmass
    v2i     = dot_product(vxyz_ptmass(1:3,ip),vxyz_ptmass(1:3,ip))
    rmsmach = rmsmach + v2i/c_sound**2
 enddo
 rmsmach = sqrt(rmsmach/nptmass)
 if (rmsmach > 0.) then
    turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
    else
    turbfac = 0.
 endif
 do ip = 1,nptmass
    vxyz_ptmass(1:3,ip) = turbfac*vxyz_ptmass(1:3,ip)
 enddo

end subroutine map_turbvel

!-------------------------------------------------------------------
! Routine for reading one cube i.e. one component of the velocities
!-------------------------------------------------------------------
subroutine read_cube(filename,vgrid,nspace,iunit,ierr)
 implicit none 
 character(len=*), intent(in)  :: filename
 real(kind=4),     intent(out) :: vgrid(:,:,:)
 integer,          intent(in)  :: iunit,nspace
 integer,          intent(out) :: ierr
 integer :: i,j,k

 open(unit=iunit,file=filename,form='unformatted',iostat=ierr)
 read(iunit,iostat=ierr) (((vgrid(i,j,k),i=1,nspace),j=1,nspace),k=1,nspace)
 if (ierr /= 0) stop 'error reading turb cubes'
 close(iunit)

end subroutine read_cube

!-------------------------------------------------------------------
! Routine for estimating the size of turb box required 
!-------------------------------------------------------------------
real function rmax(nptmass,xyzhm_ptmass)
 implicit none 
 integer, intent(in)  :: nptmass
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 integer :: ip
 real    :: r2,r2max

 r2max = tiny(r2max)
 do ip = 1,nptmass
    r2 = dot_product(xyzhm_ptmass(1:3,ip),xyzhm_ptmass(1:3,ip))
    r2max = max(r2,r2max)
 enddo 
 rmax = sqrt(r2max)
 rmax = 1.2 * rmax   ! adding some distance to the egde 

 end function rmax

!-------------------------------------------------------------------
! Trilinear interpolation from the grid
!-------------------------------------------------------------------
pure real function vfield_interp(vgrid,iposx,iposy,iposz,delx,dely,delz,deli)
 implicit none 
 integer,      intent(in) :: iposx,iposy,iposz
 real(kind=4), intent(in) :: vgrid(:,:,:)
 real,         intent(in) :: delx,dely,delz,deli
 real :: velx1,velx2,vely1,vely2

 velx1 = vgrid(iposx,iposy,iposz)   + delx/deli*(vgrid(iposx+1,iposy,iposz)-vgrid(iposx,iposy,iposz))
 velx2 = vgrid(iposx,iposy+1,iposz) + delx/deli*(vgrid(iposx+1,iposy+1,iposz)-vgrid(iposx,iposy+1,iposz))
 vely1 = velx1 + dely/deli*(velx2-velx1)

 velx1 = vgrid(iposx,iposy,iposz+1)   + delx/deli*(vgrid(iposx+1,iposy,iposz+1)-vgrid(iposx,iposy,iposz+1))
 velx2 = vgrid(iposx,iposy+1,iposz+1) + delx/deli*(vgrid(iposx+1,iposy+1,iposz+1)-vgrid(iposx,iposy+1,iposz+1))
 vely2 = velx1 + dely/deli*(velx2-velx1)

 vfield_interp = vely1 + delz/deli*(vely2-vely1)

end function vfield_interp

end module turb_grid 






































