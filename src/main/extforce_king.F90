
module extforce_king 

 implicit none 
 public :: king_potential,cluster_profile
 public :: read_infile_king,write_infile_king

 real, public :: Mclust_msun = 1d3 
 real, public :: Rcore_pc    = 0.2 
 real, public :: sigma_cgs   = 1.8d5


 private 

 !--Storage for cluster profiles 
 integer, parameter :: nRmax = 1e6
 integer :: nR 
 real    :: r_profile(nRmax),rho_profile(nRmax)
 real    :: phi_profile(nRmax),force_profile(nRmax)
 real    :: Rcore,Rclust,Mclust,sigma 

 logical :: print_profile = .false. 
 logical :: print_forces  = .false. 

 namelist /extforce_kingmodel_params/ Mclust_msun,Rcore_pc,sigma_cgs 

contains 

!---------------------------------------------------------------
! Total force per unit mass and potential 
!---------------------------------------------------------------
subroutine king_potential(xi,yi,zi,extfxi,extfyi,extfzi,phi)
 real, intent(in)  :: xi,yi,zi
 real, intent(out) :: extfxi,extfyi,extfzi,phi
 real    :: r2i,ri,rho,fr,theta_angle,phi_angle 

 r2i = xi**2 + yi**2 + zi**2 
 ri  = sqrt(r2i)

 if (ri > Rclust) then 
    extfxi = 0.d0 
    extfyi = 0.d0 
    extfzi = 0.d0 
    phi = 0.d0
 else 
    phi = interp_from_profile(ri,nR,r_profile,phi_profile)
    fr  = interp_from_profile(ri,nR,r_profile,force_profile)

    if (ri < tiny(ri)) then 
       extfxi = 0.d0 
       extfyi = 0.d0 
       extfzi = 0.d0 
    else 
       theta_angle = acos(zi/ri)
       phi_angle   = atan2(yi,xi)
       extfxi = fr*sin(theta_angle)*cos(phi_angle)
       extfyi = fr*sin(theta_angle)*sin(phi_angle)
       extfzi = fr*cos(theta_angle)
    endif 
 endif 

 if (print_forces) then 
    open(2060,file='extpot_forces.dat',position='append')
    write(2060,'(9E20.10)') xi, yi, zi, ri, fr, extfxi, extfyi, extfzi, phi 
    close(2060)
 endif 

end subroutine king_potential


!
! Function to extract/interpolate from a tabulated cluster profile
!
real function interp_from_profile(ri,nR,rad_profile,A_profile)
 integer, intent(in) :: nR 
 real,    intent(in) :: ri
 real,    intent(in) :: rad_profile(nR),A_profile(nR) 
 integer :: i,j

 !--Locate the closest entries and bracket ri 
 ! ([i] closest index; [j] lower bound; [j+1] upper bound around ri)
 i = minloc(abs(rad_profile(:)-ri),1)
 j = 0
 if (i == 1) then
    j = 1
 elseif (i == nR) then
    j = i-1
 elseif (rad_profile(i) >= ri .and. rad_profile(i-1) < ri) then
    j = i-1
 elseif (rad_profile(i) < ri .and. rad_profile(i+1) >= ri) then
    j = i
 endif
 
 !--Interpolate 
 interp_from_profile = A_profile(j) + (ri-rad_profile(j))*(A_profile(j+1)  &
                     & - A_profile(j))/(rad_profile(j+1)-rad_profile(j))

end function interp_from_profile


!-----------------------------------------------------------------
!+
! Called during init; Computes cluster profile with King models 
! Produces tabulated density, potential and force as functions of R
!+
!-----------------------------------------------------------------
subroutine cluster_profile
 use units,   only:unit_density,unit_velocity,utime,udist,umass 
 use physcon, only:pc,solarm
 integer :: iR,io_clusterfile
 real    :: j,j2,k,ve2,R,dWdR,W,W0,dR,rhomin,phi,phi0,rho,dphidr,force
 real    :: dRmax_dW,dRmax_dR
 real    :: r_pc,rho_cgs,phi_cgs,force_cgs 

 print*,'Computing cluster profile with Mclust = ',Mclust_msun,' solarm;'

 !--Convert to code units 
 Mclust = Mclust_msun*solarm/umass 
 Rcore  = Rcore_pc*pc/udist 
 sigma  = sigma_cgs/unit_velocity 

 !--Use current Mclust to compute W0 with Plummer model
 phi0 = -Mclust/Rcore 
 j2   = 1.d0/(2*sigma**2)
 j    = sqrt(j2)
 W0   = -2.d0*j2*phi0

 !--Estimate k 
 ve2  = -2.d0*phi0 
 k    = (1.d0 - exp(-1.d0*j2*ve2))**(-1)

 !--Initialize 
 R    = 1.d-3
 dWdR = tiny(dWdR) 
 W    = W0 
 dR   = 1.d-2

 rhomin = 6.77d-23/unit_density
 phi = -1   ! dummy 
 rho = 10.  ! dummy 
 iR  = 0

 scan_over_R: do while (rho > rhomin .and. phi < 0.)

    phi = W/(-2.d0*j2)
    phi = min(phi,0.d0)

    dphidr = dWdR / (-2.d0*j2*Rcore) 
    force  = -1.d0*dphidr 

    !--Integrate the 2nd-order ODE to update W & dWdR, also computes rho(W)
    call Euler(W,dWdR,W0,k,j,j2,dR,R,rho)

    !--Store results in code units 
    if (iR > nRmax) stop 'number of R entries exceeded limit'
    if (iR > 0) then  ! to skip first entry 
       r_profile(iR)       = R*Rcore
       rho_profile(iR)     = rho
       phi_profile(iR)     = phi
       force_profile(iR)   = force 
    endif 
    iR = iR + 1 

    !--Constrain dR and update 
    dRmax_dR = 1.d-1 * R
    dRmax_dW = 1.d-3 * abs(W/dWdR)
    dR = min(dRmax_dR,dRmax_dW)
    R = R + dR

 enddo scan_over_R
 
 !--Current number of entries in profile 
 nR = iR - 1 

 !--Write profile to file for checking 
 if (print_profile) then 
    open(2050,file='cluster_profile.dat',status='replace',iostat=io_clusterfile)
    if (io_clusterfile /= 0) stop 'error opening cluster profile file'
    write(2050,'(4A20)') 'r [pc]','rho [g/cm3]','potential [cm2/s2]','force [cm/s2]'
    do iR = 1,nR
       r_pc       = r_profile(iR)
       rho_cgs    = rho_profile(iR)*unit_density
       phi_cgs    = phi_profile(iR)*udist**2/utime**2
       force_cgs  = force_profile(iR)*umass*udist/utime**2
       write(2050,'(4E20.10)') r_pc, rho_cgs, phi_cgs, force_cgs 
    enddo 
    close(2050)
 endif 

 !--Cluster truncation radius 
 Rclust = R*Rcore

end subroutine cluster_profile

!
! Integrate W and dWdR with Euler method 
!
subroutine Euler(W,dWdR,W0,k,j,j2,dR,R,rho)
 real, intent(in)    :: W0,k,j,j2,dR,R
 real, intent(inout) :: W,dWdR
 real, intent(out)   :: rho
 real :: d2WdR2

 rho = rho_as_func_of_W(W,W0,k,j)
 d2WdR2 = d2WdR2_poisson(dWdR,R,rho,j2)

 dWdR = dWdR + d2WdR2 * dR 
 W = W + dWdR * dR

end subroutine Euler

!
! Expression for rho(W) derived from distribution function 
!
real function rho_as_func_of_W(W,W0,k,j)
 use physcon, only:pi
 real, intent(in) :: W,W0,k,j
 real :: integral 

 !print*,'W,W0,k,j in func',W,W0,k,j
 integral = simpson_rule_equi(rhoW_integrand,0.d0,W,1000)
 !print*,'integral',integral 
 rho_as_func_of_W = 4.d0/3.d0*pi*k*j**(-3)*exp(W-W0)*integral

end function rho_as_func_of_W

!
! Integrand in the expression of rho(W)
!
real function rhoW_integrand(eta)
 real, intent(in) :: eta

 rhoW_integrand = exp(-1.d0*eta) * eta**(3.d0/2.d0)

end function rhoW_integrand

!
! Simpson-rule with equidistant spacing
!
real function simpson_rule_equi(func,a,b,ninterval)
 integer, intent(in) :: ninterval    ! number of sub-intervals
 real,    intent(in) :: a,b          ! boundary values
 real,    external   :: func         ! the function to be integrated
 integer :: i
 real    :: dx,x1,x2,xm,f1,f2,fm,intsum

 dx  = (b-a)/dble(ninterval)
 x1  = a                      ! left
 f1  = func(a)
 intsum = 0.d0
 do i = 1,ninterval
    x2  = a + dble(i)*dx      ! right
    xm  = 0.5d0*(x1+x2)       ! midpoint
    f2  = func(x2)
    fm  = func(xm)
    intsum = intsum + (f1+4.d0*fm+f2)/6.d0*(x2-x1)  ! Simpson rule
    x1  = x2
    f1  = f2                  ! save for next subinterval
 enddo
 simpson_rule_equi = intsum

end function simpson_rule_equi

!
! Expression for d^2W/dR^2 derived from Poisson equation 
!
real function d2WdR2_poisson(dWdR,R,rhoW,j2)
 use physcon, only:pi
 real, intent(in) :: dWdR,R,rhoW,j2

 d2WdR2_poisson = -8.d0*pi*j2*Rcore**2*rhoW - 2.d0/R*dWdR 

end function d2WdR2_poisson


!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_king(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=extforce_kingmodel_params,iostat=rc)
 if (rc /= 0) stop 'cannot read king options'

end subroutine read_infile_king


subroutine write_infile_king(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=extforce_kingmodel_params,iostat=rc)
 if (rc /= 0) stop 'cannot write king options'

end subroutine write_infile_king


end module extforce_king 