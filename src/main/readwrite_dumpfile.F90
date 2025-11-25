
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
subroutine write_first_dump(time,nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)
 integer, intent(in) :: nptmass 
 real,    intent(in) :: time 
 real,    intent(in) :: xyzhm_ptmass(:,:)
 real,    intent(in) :: vxyz_ptmass(:,:)
 real,    intent(in), optional :: sq_ptmass(:,:)
 integer :: rc,i

 open(2010,file='ptmass_00000.tmp',iostat=rc,status='replace')
 if (rc /= 0) stop 'error writing first dump'

 write(2010,'(A20)') 'time'
 write(2010,'(E20.10)') time

 binary: if (present(sq_ptmass)) then 
    write(2010,'(11A20)') 'ID','x','y','z','r_acc','mass','sep','q','vx','vy','vz'
    do i = 1,nptmass 
        write(2010,'(I20,10E20.10)') i, xyzhm_ptmass(:,i), sq_ptmass(:,i), vxyz_ptmass(:,i)
    enddo 
 else 
    write(2010,'(9A20)') 'ID','x','y','z','r_acc','mass','vx','vy','vz'
    do i = 1,nptmass 
        write(2010,'(I20,8E20.10)') i, xyzhm_ptmass(:,i), vxyz_ptmass(:,i)
    enddo 
 endif binary 

 close(2010)

 print*,'First dump written to ptmass_00000.tmp'

end subroutine write_first_dump

!-------------------------------------------------------------------
! Read a specified dumpfile 
!-------------------------------------------------------------------
subroutine read_dump(dumpfile,time,nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)
 character(len=*), intent(in) :: dumpfile
 integer, intent(out) :: nptmass 
 real,    intent(out) :: time 
 real,    intent(out) :: xyzhm_ptmass(:,:)
 real,    intent(out) :: vxyz_ptmass(:,:)
 real,    intent(out), optional :: sq_ptmass(:,:)
 integer :: i,id,rc,nentry
 logical :: iex,binary

 inquire(file=trim(adjustl(dumpfile)),exist=iex)
 if (.not.iex) stop 'error finding dump'

 open(2011,file=trim(adjustl(dumpfile)),iostat=rc,status='old')
 if (rc /= 0) stop 'error reading dump'

 !--Count number of lines to find out nptmass 
 nentry = 0
 do 
    read(2011,*,iostat=rc)
    if (rc < 0) exit 
    nentry = nentry + 1 
 enddo 
 nptmass = nentry - 3   ! exclude the line for time and headings 

 rewind(2011)

 !--Read file 
 read(2011,*) 
 read(2011,'(10E20.10)') time 
 read(2011,*) 

 binary = .false. 
 if (present(sq_ptmass)) then 
    do i = 1,nptmass 
        read(2011,'(I20,10E20.10)') id, xyzhm_ptmass(:,i), sq_ptmass(:,i), vxyz_ptmass(:,i)
    enddo 
    binary = .true. 
 else 
    do i = 1,nptmass 
        read(2011,'(I20,8E20.10)') id, xyzhm_ptmass(:,i), vxyz_ptmass(:,i)
    enddo 
 endif

 close(2011)


#ifdef BINARY
 print*,'BINARY ON'
 if (.not.binary) stop 'Please re-compile with BINARY=no'
#else 
 print*,'BINARY OFF'
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

 !--Check the last snapshot file saved
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
 print*,'Newest dumpfile found: ',lastfile
 
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

subroutine write_dump(time,nptmass,xyzhm_ptmass,vxyz_ptmass,sq_ptmass)
 integer, intent(in) :: nptmass 
 real,    intent(in) :: time 
 real,    intent(in) :: xyzhm_ptmass(:,:)
 real,    intent(in) :: vxyz_ptmass(:,:)
 real,    intent(in), optional :: sq_ptmass(:,:)
 integer :: i,rc
 character(len=16) :: dumpfilename

 iunit = iunit + 1 
 idump = idump + 1 

 call gen_filename(idump,dumpfilename)
 open(iunit,file=dumpfilename,iostat=rc,status='replace')
 if (rc /= 0) stop 'error writing dump'
 print*,'writing dump to ',trim(dumpfilename)

 write(iunit,'(A20)') 'time'
 write(iunit,'(E20.10)') time 

 binary: if (present(sq_ptmass)) then 
    write(iunit,'(11A20)') 'ID','x','y','z','r_acc','mass','sep','q','vx','vy','vz'
    do i = 1,nptmass 
        write(iunit,'(I20,10E20.10)') i, xyzhm_ptmass(:,i), sq_ptmass(:,i), vxyz_ptmass(:,i)
    enddo 
 else 
    write(iunit,'(9A20)') 'ID','x','y','z','r_acc','mass','vx','vy','vz'
    do i = 1,nptmass 
        write(iunit,'(I20,8E20.10)') i, xyzhm_ptmass(:,i), vxyz_ptmass(:,i)
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

 open(2040,file='energ_angmom.ev',status='replace')
 write(2040,'(13A20)') 'time','ekin','epot','etot','jx','jy','jz','jbx','jby','jbz','jtotx','jtoty','jtotz'
 close(2040)

end subroutine restart_evfile



subroutine write_evfile(time,ekin,epot,etot,jxyz,jspin,jtot)
 real, intent(in) :: time,ekin,epot,etot 
 real, intent(in) :: jxyz(3)
 real, intent(in), optional :: jspin(3)
 real, intent(in), optional :: jtot(3)

 open(2040, file='energ_angmom.ev',status='old',position='append')
 binary: if (present(jspin)) then 
    write(2040,'(13E20.10)') time,ekin,epot,etot,jxyz(1:3),jspin(1:3),jtot(1:3)
 else 
    write(2040,'(7E20.10)')  time,ekin,epot,etot,jxyz(1:3)
 endif binary 
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