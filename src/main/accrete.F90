
module accrete

 implicit none 
 public :: get_accretion_radius
 public :: accrete_gas
 public :: read_infile_accrete,write_infile_accrete,write_accfile
 
 integer, public :: iaccrete      = 1         ! Option to let sinks accrete
 real,    public :: Mcloud_solarm = 1.d+3     ! Cloud mass [msun]
 real,    public :: Rcloud_pc     = 0.2       ! Cloud radius [pc]
 real,    public :: rho_cgs       = 4.d-20    ! Cloud density [g/cm3]
 real,    public :: angvel_cgs    = 3.d-14    ! Cloud angular velocity [rad/s]
 real,    public :: cs_cgs        = 2.19d+4   ! Sound speed [cm/s]

 integer, public :: isink_print   = 1         ! ID of sink of concern 
 logical, public :: print_r_acc   = .true.    ! Option to write individual r_acc terms 
 logical, public :: print_j_range = .true.    ! Option to write accretable j range of sinks
 logical, public :: print_v_rel   = .true.    ! Option to write relative sink-gas velocity 

 private

 namelist /accrete_params/ iaccrete,Mcloud_solarm,Rcloud_pc,rho_cgs,angvel_cgs,cs_cgs 

contains 
!+--------------------------------------------------------------
!
! Gas motion within cloud rotating about the z-axis
!
!+--------------------------------------------------------------
subroutine get_cloud_rotation(x,y,z,r,vx_gas,vy_gas,vtan_gas,jz_gas)
 use units, only:utime 
 real, intent(in)  :: x,y,z,r
 real, intent(out) :: vx_gas,vy_gas,vtan_gas,jz_gas
 real  :: angvel 

 angvel = angvel_cgs*utime
 vtan_gas = r * angvel        ! tangential velocity 
 jz_gas   = r * vtan_gas      ! Jz 

 vx_gas = vtan_gas*(-y/r)     ! x-component of gas velocity
 vy_gas = vtan_gas*(x/r)      ! y-component of gas velocity 

end subroutine get_cloud_rotation

!
! Routine to compute relative sink-gas velocity 
! in tangential direction on xy-plane
!
subroutine get_vtan_sinkgas(xi,yi,vxi,vyi,vzi,vtan_gas,vtan_sink,dv_sinkgas)
 real, intent(in)  :: xi,yi,vxi,vyi,vzi,vtan_gas
 real, intent(out) :: vtan_sink,dv_sinkgas
 real  :: rtan_mag
 real  :: v_vec(3),rtan_vec(3),rtan_uvec(3)

 v_vec = (/ vxi, vyi, vzi /)

 rtan_vec  = (/ -yi, xi, 0.d0 /)  ! tangential vec anticlockwise
 rtan_mag  = sqrt(xi*xi + yi*yi)
 rtan_uvec = (/ -yi/rtan_mag, xi/rtan_mag, 0.d0 /) 
 vtan_sink = dot_product(v_vec,rtan_uvec)

 dv_sinkgas = abs(vtan_sink - vtan_gas)

end subroutine get_vtan_sinkgas



!+--------------------------------------------------------------
!
! Estimate the accretion radii and store as h
!
!+--------------------------------------------------------------
subroutine get_accretion_radius(nptmass,xyzhm_ptmass,vxyz_ptmass)
 use ptmass, only:racc_ptmass,vrel_ptmass
 integer, intent(in)    :: nptmass 
 real,    intent(in)    :: vxyz_ptmass(:,:)
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 integer :: i,j
 real    :: r_acc,r_hill,r_bondihoyle,r_tidalneigh,rmin_tidalneigh
 real    :: xi,yi,zi,ri,mi,xj,yj,zj,rj,mj,rij
 real    :: vxi,vyi,vzi,vi,vtan_sink,vtan_gas,vx_gas,vy_gas,dv_sinkgas,jz_gas
 
 do i = 1,nptmass 
    r_acc = huge(r_acc)

    !--Roche/Hill sphere radius 
    xi = xyzhm_ptmass(1,i)
    yi = xyzhm_ptmass(2,i)
    zi = xyzhm_ptmass(3,i)
    ri = sqrt(xi*xi + yi*yi + zi*zi)
    mi = xyzhm_ptmass(5,i)
    r_hill = hillsphere_radius(ri,mi,nptmass,xyzhm_ptmass)
    r_acc = min(r_acc,r_hill)

    !--Bondi-Hyole radius 
    vxi = vxyz_ptmass(1,i)
    vyi = vxyz_ptmass(2,i)
    vzi = vxyz_ptmass(3,i)
    vi  = sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
    call get_cloud_rotation(xi,yi,zi,ri,vx_gas,vy_gas,vtan_gas,jz_gas)
    call get_vtan_sinkgas(xi,yi,vxi,vyi,vzi,vtan_gas,vtan_sink,dv_sinkgas)
    r_bondihoyle = bondihoyle_radius(mi,dv_sinkgas)
    r_acc = min(r_acc,r_bondihoyle)

    !--Tidal radius wrt neighbours 
    rmin_tidalneigh = huge(rmin_tidalneigh)
    do j = 1,nptmass 
       if (i /= j) then 
          xj  = xyzhm_ptmass(1,j)
          yj  = xyzhm_ptmass(2,j)
          zj  = xyzhm_ptmass(3,j)
          rj  = sqrt(xj*xj + yj*yj + zj*zj)
          rij = abs(ri - rj)
          mj  = xyzhm_ptmass(5,j)
          r_tidalneigh = tidalneigh_radius(mi,mj,rij)
          rmin_tidalneigh = min(rmin_tidalneigh,r_tidalneigh)
       endif 
    enddo 
    r_acc = min(r_acc,rmin_tidalneigh)

    !--Store results 
    xyzhm_ptmass(4,i) = r_acc
    racc_ptmass(1,i)  = r_hill
    racc_ptmass(2,i)  = r_bondihoyle
    racc_ptmass(3,i)  = rmin_tidalneigh
    vrel_ptmass(1,i)  = vtan_sink
    vrel_ptmass(2,i)  = vtan_gas
    vrel_ptmass(3,i)  = dv_sinkgas
 enddo 

end subroutine get_accretion_radius

!
! Tidal radius with respect to the whole cluster 
! aka Jacobi radius, or Hill sphere, or Roche sphere 
!
real function hillsphere_radius(ri,mi,nptmass,xyzhm_ptmass)
 use physcon, only:solarm,pc 
 use units,   only:umass,udist 
 integer, intent(in) :: nptmass 
 real,    intent(in) :: ri,mi
 real,    intent(in) :: xyzhm_ptmass(:,:)
 integer :: j 
 real    :: xj,yj,zj,rj,mj
 real    :: Mcloud,Rcloud,M_enclosed,r_enclosed

 Mcloud = Mcloud_solarm*solarm/umass 
 Rcloud = Rcloud_pc*pc/udist 
 r_enclosed = min(ri,Rcloud)
 M_enclosed = Mcloud * (r_enclosed/Rcloud)**3

 !--include sink masses 
 do j = 1,nptmass 
    xj = xyzhm_ptmass(1,j)
    yj = xyzhm_ptmass(2,j)
    zj = xyzhm_ptmass(3,j)
    rj = sqrt(xj*xj + yj*yj * zj*zj)
    if (rj < r_enclosed-tiny(r_enclosed)) then 
       mj = xyzhm_ptmass(5,j)
       M_enclosed = M_enclosed + mj 
    endif 
 enddo 

 hillsphere_radius = (mi/(3.d0*M_enclosed))**(1.d0/3.d0) * ri

end function hillsphere_radius

!
! Tidal radius with respect to its neighbours j
!
real function tidalneigh_radius(mi,mj,rij)
 real, intent(in)  :: mi,mj,rij

 tidalneigh_radius = (mi/(2.d0*mj))**(1.d0/3.d0) * rij 
 
end function tidalneigh_radius

!
! Bondi-Hyole radius 
!
real function bondihoyle_radius(mi,dv_sinkgas)
 use units, only:unit_velocity 
 real, intent(in)  :: mi,dv_sinkgas
 real  :: cs 

 cs = cs_cgs/unit_velocity
 bondihoyle_radius = 2.d0*mi/(dv_sinkgas**2+cs**2)

end function bondihoyle_radius


!+--------------------------------------------------------------
!
! Accrete and update ptmass properties 
!
!+--------------------------------------------------------------
subroutine accrete_gas(dt,nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass,sq_ptmass)
 use units,   only:unit_density 
 use physcon, only:pi 
 use ptmass,  only:jrange_ptmass 
 use ptmass,  only:compute_Lxyz
#ifdef BINARY
 use ptmass,  only:compute_Lspin,update_sq
#endif 
 integer, intent(in)    :: nptmass 
 real,    intent(in)    :: dt 
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 real,    intent(inout), optional :: sq_ptmass(:,:)
 integer :: i
 real    :: xi,yi,zi,mi,ri,vxi,vyi,vzi,fxi,fyi,fzi,si,qi
 real    :: jspin(3),Lspin(3),jxyz(3),Lxyz(3)
 real    :: r_acc,jx_sink,jy_sink,jz_sink,jz_min,jz_max,jz_gas
 real    :: jz_gas_inner,jz_gas_outer
 real    :: vtan_gas,vtan_sink
 real    :: dv_sinkgas,vx_gas,vy_gas,rho,jgas_min,jgas_max
 real    :: dm,dxm,dym,dzm,dvxm,dvym,dvzm,dfxm,dfym,dfzm,dLx,dLy,dLz
 real    :: mnew,mnew1 
 logical :: accretable

 do i = 1,nptmass 

    !--Unpack variables
    xi = xyzhm_ptmass(1,i)
    yi = xyzhm_ptmass(2,i)
    zi = xyzhm_ptmass(3,i)
    ri = sqrt(xi*xi + yi*yi + zi*zi)
    mi = xyzhm_ptmass(5,i)
    r_acc = xyzhm_ptmass(4,i)
    vxi = vxyz_ptmass(1,i)
    vyi = vxyz_ptmass(2,i)
    vzi = vxyz_ptmass(3,i)
    fxi = fxyz_ptmass(1,i)
    fyi = fxyz_ptmass(2,i)
    fzi = fxyz_ptmass(3,i)
#ifdef BINARY
    si = sq_ptmass(1,i)
    qi = sq_ptmass(2,i)
#endif 

    !--Compute the accretable range of j 
    jxyz  = 0.d0 
    Lxyz  = 0.d0 
    jspin = 0.d0 
    Lspin = 0.d0 
    call compute_Lxyz(i,xi,yi,zi,mi,vxi,vyi,vzi,jxyz,Lxyz)
#ifdef BINARY
    call compute_Lspin(i,mi,si,qi,jspin,Lspin)
#endif 
    jx_sink = jxyz(1) + jspin(1)
    jy_sink = jxyz(2) + jspin(2)
    jz_sink = jxyz(3) + jspin(3)
    call compute_accretable_jrange(r_acc,mi,jz_sink,jz_min,jz_max)

    !--Check if the cloud gas around ri falls within j-range 
    call check_gas_accretable(ri,r_acc,jz_min,jz_max,jgas_min,jgas_max,accretable)

    !--Store accretable j-range info 
    jrange_ptmass(1,i) = jx_sink 
    jrange_ptmass(2,i) = jy_sink 
    jrange_ptmass(3,i) = jz_sink 
    jrange_ptmass(4,i) = jz_min 
    jrange_ptmass(5,i) = jz_max 
    jrange_ptmass(6,i) = jgas_min 
    jrange_ptmass(7,i) = jgas_max 


    if (iaccrete > 0 .and. accretable) then 

       call get_cloud_rotation(xi,yi,zi,ri,vx_gas,vy_gas,vtan_gas,jz_gas)

       !--Accrete mass 
       call get_vtan_sinkgas(xi,yi,vxi,vyi,vzi,vtan_gas,vtan_sink,dv_sinkgas)
       rho = rho_cgs/unit_density 
       dm  = pi*r_acc**2 * rho * dv_sinkgas * dt 

       !--Estimate the accreted dr,dv,df,dL 
       ! [Note: In sph, this would be summed over the accreted gas particles]
       dxm  = xi*dm
       dym  = yi*dm 
       dzm  = zi*dm 
       dvxm = vx_gas*dm 
       dvym = vy_gas*dm 
       dvzm = 0.d0 
       dfxm = fxi*dm       ! assuming the gas experiences similar accel as the sink
       dfym = fyi*dm 
       dfzm = fzi*dm 
       dLx  = 0.d0         ! assuming it only accretes in the xy-plane
       dLy  = 0.d0 
       dLz  = jz_gas*dm    

       !--Update sink properties 
       ! [Note: dLspin = Lxyz_before + L_accreted - Lxyz_after]
       Lspin(1) = Lspin(1) + Lxyz(1) + dLx 
       Lspin(2) = Lspin(2) + Lxyz(2) + dLy
       Lspin(3) = Lspin(3) + Lxyz(3) + dLz 
       mnew  = mi + dm 
       mnew1 = 1.d0/mnew 
       xi  = mnew1*(dxm + xi*mi)       ! Reset COM 
       yi  = mnew1*(dym + yi*mi)
       zi  = mnew1*(dzm + zi*mi)
       vxi = mnew1*(dvxm + vxi*mi)     ! Conserve linear momentum 
       vyi = mnew1*(dvym + vyi*mi)
       vzi = mnew1*(dvzm + vzi*mi)
       fxi = mnew1*(dfxm + fxi*mi)     ! Update net force 
       fyi = mnew1*(dfym + fyi*mi)
       fzi = mnew1*(dfzm + fzi*mi)
       mi  = mnew 
       call compute_Lxyz(i,xi,yi,zi,mi,vxi,vyi,vzi,jxyz,Lxyz)     ! Lxyz_after
       Lspin(1) = Lspin(1) - Lxyz(1)
       Lspin(2) = Lspin(2) - Lxyz(2)
       Lspin(3) = Lspin(3) - Lxyz(3)
#ifdef BINARY
       call update_sq(mi,qi,Lspin(3),si)
       call compute_Lspin(i,mi,si,qi,jspin,Lspin)
#endif 
       
       !--Store updated properties 
       xyzhm_ptmass(1,i) = xi 
       xyzhm_ptmass(2,i) = yi
       xyzhm_ptmass(3,i) = zi 
       xyzhm_ptmass(5,i) = mi
       vxyz_ptmass(1,i)  = vxi 
       vxyz_ptmass(2,i)  = vyi 
       vxyz_ptmass(3,i)  = vzi 
       fxyz_ptmass(1,i)  = fxi 
       fxyz_ptmass(2,i)  = fyi 
       fxyz_ptmass(3,i)  = fzi 
#ifdef BINARY
       sq_ptmass(1,i) = si 
       sq_ptmass(2,i) = qi
#endif 
    endif 
 enddo 

end subroutine accrete_gas 


!
! Window of j_gas within which the gas parcel would be bound to the sink 
!
subroutine compute_accretable_jrange(r_acc,mi,jz_sink,jz_min,jz_max)
 real, intent(in)  :: mi,r_acc,jz_sink
 real, intent(out) :: jz_min,jz_max 
 real  :: deltajz

 !--Equate centripetal to gravity, or 2KE to GPE 
 deltajz = sqrt(mi*r_acc)

 jz_min = jz_sink - deltajz
 jz_max = jz_sink + deltajz 

end subroutine compute_accretable_jrange


!
! Check if the j of gas around ri falls within [jz_min,jz_max]
!
subroutine check_gas_accretable(ri,r_acc,jz_min,jz_max,jz_gas_inner,jz_gas_outer,accretable)
 real,    intent(in)  :: ri,r_acc
 real,    intent(in)  :: jz_min,jz_max
 real,    intent(out) :: jz_gas_inner,jz_gas_outer
 logical, intent(out) :: accretable
 real    :: xi,yi,zi,vx_gas,vy_gas,vtan_gas

 xi = -1  ! dummy 
 yi = -1 
 zi = -1 
 call get_cloud_rotation(xi,yi,zi,ri-r_acc,vx_gas,vy_gas,vtan_gas,jz_gas_inner)
 call get_cloud_rotation(xi,yi,zi,ri+r_acc,vx_gas,vy_gas,vtan_gas,jz_gas_outer) 

 accretable = .false. 
 if (jz_gas_inner > jz_min .and. jz_gas_inner < jz_max .and. &
   & jz_gas_outer > jz_min .and. jz_gas_outer < jz_max) then 
    accretable = .true. 
 endif 

end subroutine check_gas_accretable



!--------------------------------------------------------------
! Write computed accretion properties to file
!--------------------------------------------------------------
subroutine write_accfile(time)
 use ptmass, only:racc_ptmass,jrange_ptmass,vrel_ptmass
 real, intent(in) :: time 

 if (print_r_acc) then 
    open(2090,file='accretion_radii.ev',status='old',position='append')
    write(2090,'(1E20.10,I20,3E20.10)') time,isink_print,racc_ptmass(1:3,isink_print)
    close(2090)
 endif 

 if (print_j_range) then 
    open(2100,file='accretable_j_range.ev',status='old',position='append')
    write(2100,'(1E20.10,I20,7E20.10)') time,isink_print,jrange_ptmass(1:7,isink_print)
    close(2100)
 endif 

 if (print_v_rel) then 
    open(2200,file='vrel_sinkgas.ev',status='old',position='append')
    write(2200,'(1E20.10,I20,7E20.10)') time,isink_print,vrel_ptmass(1:3,isink_print)
    close(2200)
 endif 

end subroutine write_accfile



!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_accrete(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=accrete_params,iostat=rc)
 if (rc /= 0) stop 'cannot read accrete options'

end subroutine read_infile_accrete


subroutine write_infile_accrete(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=accrete_params,iostat=rc)
 if (rc /= 0) stop 'cannot write accrete options'

end subroutine write_infile_accrete


end module accrete 