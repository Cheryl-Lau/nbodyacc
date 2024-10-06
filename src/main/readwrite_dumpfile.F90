
module readwrite_dump

 implicit none 
 public :: write_first_dump,get_first_dump,read_dump,write_dump
 public :: restart_evfile,write_evfile
 public :: read_infile_startdump,write_infile_startdump
 character(len=16), public :: start_dump = 'ptmass_00000.tmp'  ! default 

 private
 namelist /startdump_params/ start_dump

 integer :: iunit = 3000
 integer :: idump = 0

contains 

!-------------------------------------------------------------------
! Write the starting dumpfile 
!-------------------------------------------------------------------
subroutine write_first_dump(time,nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)
 integer, intent(in) :: nptmass 
 real,    intent(in) :: time 
 real,    intent(in) :: xyzhm_ptmass(:,:)
 real,    intent(in) :: vxyz_ptmass(:,:)
 real,    intent(in), optional :: massq_ptmass(:)
 integer :: rc,ip

 open(2010,file='ptmass_00000.tmp',iostat=rc,status='replace')
 if (rc /= 0) stop 'error writing first dump'

 write(2010,*) time

 binary: if (present(massq_ptmass)) then 
    do ip = 1,nptmass 
        write(2010,*) xyzhm_ptmass(:,ip), massq_ptmass(ip), vxyz_ptmass(:,ip)
    enddo 
 else 
    do ip = 1,nptmass 
        write(2010,*) xyzhm_ptmass(:,ip), vxyz_ptmass(:,ip)
    enddo 
 endif binary 

 close(2010)

 print*,'First dump written to ptmass_00000.tmp'

end subroutine write_first_dump

!-------------------------------------------------------------------
! Read a specified dumpfile 
!-------------------------------------------------------------------
subroutine read_dump(dumpfile,time,nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)
 character(len=*), intent(in) :: dumpfile
 integer, intent(out) :: nptmass 
 real,    intent(out) :: time 
 real,    intent(out) :: xyzhm_ptmass(:,:)
 real,    intent(out) :: vxyz_ptmass(:,:)
 real,    intent(out), optional :: massq_ptmass(:)
 integer :: ip,rc,nentry
 logical :: iex,binary

 inquire(file=trim(adjustl(dumpfile)),exist=iex)
 if (.not.iex) stop 'error finding dump'

 open(2011,file=trim(adjustl(dumpfile)),iostat=rc,status='old')
 if (rc /= 0) stop 'error reading dump'

 !- Count number of lines to find out nptmass 
 nentry = 0
 do 
    read(2011,*,iostat=rc)
    if (rc < 0) exit 
    nentry = nentry + 1 
 enddo 
 nptmass = nentry - 1  ! exclude the line for time

 rewind(2011)

 !- Read file 
 read(2011,*) time 
 binary = .false. 
 if (present(massq_ptmass)) then 
    do ip = 1,nptmass 
        read(2011,*) xyzhm_ptmass(:,ip), massq_ptmass(ip), vxyz_ptmass(:,ip)
    enddo 
    binary = .true. 
 else 
    do ip = 1,nptmass 
        read(2011,*) xyzhm_ptmass(:,ip), vxyz_ptmass(:,ip)
    enddo 
 endif

 close(2011)

#ifdef BINARY
 print*,'BINARY on'
 if (.not.binary) stop 'Please re-compile with BINARY=no'
#else 
 print*,'BINARY off'
 if (binary) stop 'Please re-compile with BINARY=yes'
#endif 

end subroutine read_dump

!-------------------------------------------------------------------
! Set the first dumpfile to begin
!-------------------------------------------------------------------
subroutine get_first_dump(starting_dump)
 character(len=16), intent(out) :: starting_dump
 integer :: ifile_search
 logical :: lastfile_found,iexist
 character(len=16) :: filename_search,lastfile

 !- Check the last snapshot file saved
 ifile_search = 100
 lastfile_found = .false.
 do while (.not.lastfile_found)
    call gen_filename(ifile_search,filename_search)
    inquire(file=trim(adjustl(filename_search)),exist=iexist)
    if (iexist .and. ifile_search > 0) then
        lastfile_found = .true.
        exit 
    else
        ifile_search = ifile_search - 1
        if (ifile_search < 0) exit 
    endif
 enddo
 lastfile = trim(adjustl(filename_search))
 print*,'Newest dumpfile found ',lastfile
 
 if (lastfile_found .and. start_dump /= lastfile) then 
    print*,'ARE YOU SURE YOU DO NOT WANT TO START FROM ',lastfile,' ?'
    call sleep(3)
 endif 

 if (start_dump == 'ptmass_00000.tmp') then
    call rename('ptmass_00000.tmp','ptmass_00000.dat')
    starting_dump = 'ptmass_00000.dat'
    idump = 0 
 else 
    starting_dump = start_dump  ! from input file 
    read(starting_dump(8:12),*) idump 
 endif 

 print*,'Beginning simulation from:', starting_dump,'; label: ', idump 

end subroutine get_first_dump


!-------------------------------------------------------------------
! Wrapper to write dumpfiles
!-------------------------------------------------------------------

subroutine write_dump(time,nptmass,xyzhm_ptmass,vxyz_ptmass,massq_ptmass)
 integer, intent(in) :: nptmass 
 real,    intent(in) :: time 
 real,    intent(in) :: xyzhm_ptmass(:,:)
 real,    intent(in) :: vxyz_ptmass(:,:)
 real,    intent(in), optional :: massq_ptmass(:)
 integer :: ip,rc
 character(len=16) :: dumpfilename

 iunit = iunit + 1 
 idump = idump + 1 

 call gen_filename(idump,dumpfilename)
 open(iunit,file=dumpfilename,iostat=rc,status='replace')
 if (rc /= 0) stop 'error writing dump'

 write(iunit,*) time 
 binary: if (present(massq_ptmass)) then 
    do ip = 1,nptmass 
        write(iunit,*) xyzhm_ptmass(:,ip), massq_ptmass(ip), vxyz_ptmass(:,ip)
    enddo 
 else 
    do ip = 1,nptmass 
        write(iunit,*) xyzhm_ptmass(:,ip), vxyz_ptmass(:,ip)
    enddo 
 endif binary 

 close(iunit)

end subroutine write_dump


!-------------------------------------------------------------------
! Generate filenames 
!-------------------------------------------------------------------
subroutine gen_filename(ifile,filename)
 integer,           intent(in)  :: ifile
 character(len=16), intent(out) :: filename
 character(len=5)   :: ifile_char

 write(ifile_char,'(i5.5)') ifile  !- convert to str
 filename = 'ptmass_'//trim(ifile_char)//'.dat'

end subroutine gen_filename


!-------------------------------------------------------------------
! Write ev file 
!-------------------------------------------------------------------
subroutine restart_evfile()

 open(2040,file='ptmass.ev',status='replace')
 write(2040,*) 'time','ekin','epot'
 close(2040)

end subroutine restart_evfile


subroutine write_evfile(time,ekin,epot)
 real, intent(in) :: time,ekin,epot

 open(2040, file='ptmass.ev',status='old',position='append')
 write(2040,*) time,ekin,epot
 close(2040)

end subroutine write_evfile


!-------------------------------------------------------------------
! Write into / Read from the input file 
!-------------------------------------------------------------------
subroutine read_infile_startdump(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 read(unit_infile,nml=startdump_params,iostat=rc)
 if (rc /= 0) stop 'cannot read startdump options'

end subroutine read_infile_startdump


subroutine write_infile_startdump(unit_infile)
 integer, intent(in) :: unit_infile
 integer :: rc

 write(unit_infile,nml=startdump_params,iostat=rc)
 if (rc /= 0) stop 'cannot write startdump options'

end subroutine write_infile_startdump

end module readwrite_dump