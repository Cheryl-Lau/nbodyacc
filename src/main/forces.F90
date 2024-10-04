
module force 

 implicit none 
 public :: compute_forces

 private 

contains 

subroutine compute_forces()

 call plumer_potential()
#ifdef GRAVITY
 call self_grav()
#endif 

end subroutine compute_forces


subroutine self_grav()

end subroutine self_grav


subroutine plumer_potential()

end subroutine plumer_potential


subroutine galdisc_potential()

end subroutine galdisc_potential

end module force 