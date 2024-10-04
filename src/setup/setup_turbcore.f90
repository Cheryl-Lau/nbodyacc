
module setup 

 implicit none 
 public :: set_ptmass,read_setup_infile,write_dump,write_infile

 private

contains 

subroutine set_ptmass(nptmass,xyzhm_ptmass,vxyz_ptmass,fxyz_ptmass)
 integer, intent(out)   :: nptmass
 real,    intent(inout) :: xyzhm_ptmass(:,:)
 real,    intent(inout) :: vxyz_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:)
 integer :: ip 

 !nptmass = 10 
 !do ip = 1,nptmass
 !  xyzhm_ptmass(:,ip) = (/1,1,1,1,1/)
 !   vxyz_ptmass(:,ip) = (/1,1,1/)
 !   fxyz_ptmass(:,ip) = (/1,1,1/)
 !enddo 

end subroutine set_ptmass


subroutine read_setup_infile()


end subroutine read_setup_infile

subroutine write_dump()

end subroutine write_dump


subroutine write_infile()

end subroutine write_infile

end module setup 