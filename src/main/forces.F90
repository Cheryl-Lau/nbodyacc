
module force 

 implicit none 
 public :: compute_forces
 public :: read_infile_force,write_infile_force

 integer, public :: iextforce = 3

 integer, parameter, public :: &
   iext_plummer  = 1, &
   iext_galactic = 2, &
   iext_king     = 3

 private 

 namelist /force_params/ iextforce

contains 

!---------------------------------------------------------------
! Total force per unit mass (acceleration)
!---------------------------------------------------------------
subroutine compute_forces(nptmass,xyzhm_ptmass,fxyz_ptmass)
 use ptmass, only:poten_ptmass 
 integer, intent(in)  :: nptmass
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(out) :: fxyz_ptmass(:,:)
 integer :: i
 real    :: xi,yi,zi,extfxi,extfyi,extfzi,phi

 fxyz_ptmass  = 0.d0  ! init 
 poten_ptmass = 0.d0

#ifdef GRAVITY
 call self_grav(nptmass,xyzhm_ptmass,fxyz_ptmass,poten_ptmass)
#endif 

 if (iextforce > 0) then 
    do i = 1,nptmass
       xi = xyzhm_ptmass(1,i)
       yi = xyzhm_ptmass(2,i)
       zi = xyzhm_ptmass(3,i)
       call external_force(xi,yi,zi,extfxi,extfyi,extfzi,phi)
       fxyz_ptmass(1,i) = fxyz_ptmass(1,i) + extfxi
       fxyz_ptmass(2,i) = fxyz_ptmass(2,i) + extfyi
       fxyz_ptmass(3,i) = fxyz_ptmass(3,i) + extfzi
       poten_ptmass(i)  = poten_ptmass(i)  + phi 
    enddo
 endif  

end subroutine compute_forces

!
! Self-gravity 
!
subroutine self_grav(nptmass,xyzhm_ptmass,fxyz_ptmass,poten_ptmass)
 integer, intent(in)    :: nptmass
 real,    intent(in)    :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:),poten_ptmass(:)
 integer :: i,j
 real    :: r_ij(3),fsum(3),f_ij(3),absr,mi,mj,phi,phisum,hi,hj,hsoft 

 do i = 1,nptmass
    fsum   = 0.d0 
    phisum = 0.d0 
    mi     = xyzhm_ptmass(5,i)
    hi     = xyzhm_ptmass(4,i)
    over_neigh: do j = 1,nptmass
        if (i /= j) then 
            mj    = xyzhm_ptmass(5,j)
            r_ij  = xyzhm_ptmass(1:3,i) - xyzhm_ptmass(1:3,j)
            absr  = sqrt(dot_product(r_ij,r_ij))
            hj    = xyzhm_ptmass(4,j)
            hsoft = 2.d0*max(hi,hj)

            softening: if (absr < hsoft) then 
               f_ij = -mj*r_ij/(absr+hsoft)**3
               fsum = fsum + f_ij 
               phi  = -mi/(absr+hsoft)
               phisum = phisum + phi 
            else 
               f_ij = -mj*r_ij/absr**3    ! per mass(i); G=1 in code units 
               fsum = fsum + f_ij 
               phi  = -mi/absr
               phisum = phisum + phi 
            endif softening 
        endif 
    enddo over_neigh 
    fxyz_ptmass(1:3,i) = fxyz_ptmass(1:3,i) + fsum
    poten_ptmass(i)    = poten_ptmass(i)    + phisum 
 enddo 

end subroutine self_grav

!
! External potential acting on a given ptmass
!
subroutine external_force(xi,yi,zi,extfxi,extfyi,extfzi,phi) 
 use extforce_plummer,  only:plummer_potential
 use extforce_galactic, only:galactic_potential
 use extforce_king,     only:king_potential
 real, intent(in)  :: xi,yi,zi
 real, intent(out) :: extfxi,extfyi,extfzi,phi

 select case(iextforce)

 case(iext_plummer) 
    call plummer_potential(xi,yi,zi,extfxi,extfyi,extfzi,phi)

 case(iext_galactic)
    call galactic_potential(xi,yi,zi,extfxi,extfyi,extfzi,phi)

 case(iext_king)
    call king_potential(xi,yi,zi,extfxi,extfyi,extfzi,phi)

 case default
    stop 'invalid iextforce'

 end select 

end subroutine external_force 


!--------------------------------------------------------------
! Write module options to input file
!--------------------------------------------------------------
subroutine read_infile_force(unit_infile)
 use extforce_plummer,  only:read_infile_plummer
 use extforce_galactic, only:read_infile_galactic
 use extforce_king,     only:read_infile_king
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=force_params,iostat=rc)
 if (rc /= 0) stop 'cannot read force options'

 select case(iextforce)
 case(iext_plummer) 
    call read_infile_plummer(unit_infile)
 case(iext_galactic)
    call read_infile_galactic(unit_infile)
 case(iext_king)
    call read_infile_king(unit_infile)
 case default
    stop 'invalid iextforce'
 end select 

end subroutine read_infile_force


subroutine write_infile_force(unit_infile)
 use extforce_plummer,  only:write_infile_plummer
 use extforce_galactic, only:write_infile_galactic
 use extforce_king,     only:write_infile_king 
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=force_params,iostat=rc)
 if (rc /= 0) stop 'cannot write force options'

 select case(iextforce)
 case(iext_plummer) 
    call write_infile_plummer(unit_infile)
 case(iext_galactic)
    call write_infile_galactic(unit_infile)
 case(iext_king)
    call write_infile_king(unit_infile)
 case default
    stop 'invalid iextforce'
 end select 

end subroutine write_infile_force


end module force 