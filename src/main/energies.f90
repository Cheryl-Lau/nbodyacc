
module energy 

 implicit none 
 public :: get_tot_energy,get_tot_angmomen

 private

contains

subroutine get_tot_energy(nptmass,xyzhm_ptmass,vxyz_ptmass,ekin,epot)
 integer, intent(in)  :: nptmass 
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(in)  :: vxyz_ptmass(:,:)
 real,    intent(out) :: ekin,epot 


end subroutine get_tot_energy

subroutine get_tot_angmomen(nptmass,xyzhm_ptmass,vxyz_ptmass,Lxyz_ptmass)
 integer, intent(in)  :: nptmass 
 real,    intent(in)  :: xyzhm_ptmass(:,:)
 real,    intent(in)  :: vxyz_ptmass(:,:)
 real,    intent(out) :: Lxyz_ptmass(:,:)


end subroutine get_tot_angmomen

end module energy 