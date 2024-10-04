
module timestep 

 implicit none 
 public :: constrain_dt 
 
 real,  public :: dtmax = 1e-3
 real,  public :: dt 
 
 private

contains 

subroutine constrain_dt()

end subroutine constrain_dt

end module timestep 