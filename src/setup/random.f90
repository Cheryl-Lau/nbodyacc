
module random
 !
 ! Taken from phantom 
 !
 implicit none  
 public :: ran2

 private 

 contains 

 real function ran2(s1)
  integer, intent(inout) :: s1
  integer, save :: s2 = 123456789
  if (s1 < 0) s2 = 123456789
  ran2 = get_random(s1,s2)
 end function ran2

 real function get_random(s1, s2)
  integer  :: k,s1,s2,z
  k = s1 / 53668
  s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
  if ( s1 < 0 ) then
     s1 = s1 + 2147483563
  endif
  k = s2 / 52774
  s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
  if ( s2 < 0 ) then
     s2 = s2 + 2147483399
  endif
  z = s1 - s2
  if ( z < 1 ) then
     z = z + 2147483562
  endif
  get_random = real ( z ) / 2147483563.0D+00
  return
 end function get_random

end module random






































