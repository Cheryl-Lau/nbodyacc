
module force 

 implicit none 
 public :: compute_forces
 public :: read_infile_force,write_infile_force

 logical, public :: use_plumerpot   = .true. 
 logical, public :: use_galdisc_pot = .false. 

 private 
 namelist /force_params/ use_plumerpot,use_galdisc_pot

contains 

!---------------------------------------------------------------
! Total force per unit mass (acceleration)
!---------------------------------------------------------------
subroutine compute_forces(nptmass,xyzhm,fxyz)
 integer, intent(in)  :: nptmass
 real,    intent(in)  :: xyzhm(:,:)
 real,    intent(out) :: fxyz(:,:)

 fxyz = 0.  ! init 

#ifdef GRAVITY
 call self_grav(nptmass,xyzhm,fxyz)
#endif 

 if (use_plumerpot) then 
    call plumer_potential(nptmass,xyzhm,fxyz)
 endif 
 if (use_galdisc_pot) then 
    call galdisc_potential(nptmass,xyzhm,fxyz)
 endif

end subroutine compute_forces

!---------------------------------------------------------------
! Self-gravity - consider doing a tree in the future 
!---------------------------------------------------------------
subroutine self_grav(nptmass,xyzhm,fxyz)
 use physcon, only:gg
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: xyzhm(:,:)
 real,    intent(inout) :: fxyz(:,:)
 integer :: i,j
 real    :: r_ij(3),fsum(3),f_ij(3),absr

 do i = 1,nptmass
    fsum = 0.d0 
    over_neigh: do j = 1,nptmass
        if (i /= j) then 
            r_ij = xyzhm(1:3,i) - xyzhm(1:3,j)
            absr = sqrt(dot_product(r_ij,r_ij))
            f_ij = -xyzhm(5,j)*r_ij/absr**3    ! per mass(i); G=1 in code units 
            fsum = fsum + f_ij 
        endif 
    enddo over_neigh 
    fxyz(1:3,i) = fxyz(1:3,i) + fsum
 enddo 

end subroutine self_grav

!---------------------------------------------------------------
! Plummer potential for modelling a cluster 
!---------------------------------------------------------------
subroutine plumer_potential(nptmass,xyzhm,fxyz)
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: xyzhm(:,:)
 real,    intent(inout) :: fxyz(:,:)

end subroutine plumer_potential

!---------------------------------------------------------------
! Galactic disc potential 
!---------------------------------------------------------------
subroutine galdisc_potential(nptmass,xyzhm,fxyz)
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: xyzhm(:,:)
 real,    intent(inout) :: fxyz(:,:)
 integer :: i,j
 real    :: r_ij(3),fsum(3),f_ij(3)
 real    :: x,y,z,fac,dphidx,dphidy,dphidz
 real    :: Rc,vc,q1,q2

 Rc = 2.d-1
 vc = 1.d0
 q1 = 1.d0
 q2 = 1.d0

 do i = 1,nptmass
    fsum = 0.d0 
    x = xyzhm(1,i)
    y = xyzhm(2,i)
    z = xyzhm(3,i)
    fac = vc**2/2.d0 * (Rc**2 + x**2 + y**2/q1**2 + z**2/q2**2)**(-1)
    dphidx = fac*2.d0 * x
    dphidy = fac*2.d0/q1**2 * y
    dphidz = fac*2.d0/q2**2 * z
    f_ij = (/-dphidx,-dphidy,-dphidz/)
    fxyz(1:3,i) = fxyz(1:3,i) + f_ij 
 enddo 

end subroutine galdisc_potential


!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_force(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=force_params,iostat=rc)
 if (rc /= 0) stop 'cannot read force options'

end subroutine read_infile_force

subroutine write_infile_force(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=force_params,iostat=rc)
 if (rc /= 0) stop 'cannot write timestep options'

end subroutine write_infile_force


end module force 