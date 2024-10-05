
module units 

 implicit none 
 public :: set_units
 real, public :: udist = 1.d0, umass = 1.d0, utime = 1.d0
 real, public :: unit_velocity,unit_density,unit_pressure,unit_ergg,unit_energ
 private

contains 

!-------------------------------------------------------------------
! Converting between code units and cgs units (stolen from phantom)
!-------------------------------------------------------------------
subroutine set_units(dist,mass,time,G,c)
 use physcon, only:gg,clight=>c
 real, intent(in), optional :: dist,mass,time,G,c

 if (present(dist)) then
    udist = dist
 else
    udist = 1.d0
 endif
 if (present(mass)) then
    umass = mass
 else
    umass = 1.d0
 endif
 if (present(time)) then
    utime = time
 else
    utime = 1.d0
 endif

 if (present(c)) then
    if (present(mass)) then
       udist = gg*umass/clight**2
       utime = udist/clight
       if (present(dist)) print "(a)",' WARNING: over-riding length unit with c=1 assumption'
    elseif (present(dist)) then
       utime = udist/clight
       umass = clight*clight*udist/gg
       if (present(time)) print "(a)",' WARNING: over-riding time unit with c=1 assumption'
    elseif (present(time)) then
       udist = utime*clight
       umass = clight*clight*udist/gg
    else
       udist = gg*umass/clight**2   ! umass is 1
       utime = udist/clight
    endif
 elseif (present(G)) then
    if (present(mass) .and. present(dist)) then
       utime = sqrt(udist**3/(gg*umass))
       if (present(time)) print "(a)",' WARNING: over-riding time unit with G=1 assumption'
    elseif (present(dist) .and. present(time)) then
       umass = udist**2/(gg*utime**2)
       if (present(mass)) print "(a)",' WARNING: over-riding mass unit with G=1 assumption'
    elseif (present(mass) .and. present(time)) then
       udist = (utime**2*(gg*umass))**(1.d0/3.d0)
       if (present(dist)) print "(a)",' WARNING: over-riding length unit with G=1 assumption'
    elseif (present(time)) then
       umass = udist**2/(gg*utime**2)     ! udist is 1
    else
       utime = sqrt(udist**3/(gg*umass))  ! udist and umass are 1
    endif
 endif

 ! derived units 
 unit_velocity = udist/utime
 unit_density  = umass/udist**3
 unit_pressure = umass/(udist*utime**2)
 unit_ergg     = unit_velocity**2
 unit_energ    = umass*unit_ergg

end subroutine set_units

end module units 