program yspec_compare
  !================================================================!
  ! Compares a yspec output file with a reference file. Each       !
  ! seismogram column passes if the largest absolute difference    !
  ! is no more than tol times the peak absolute value of the       !
  ! reference column. The time columns must agree.                 !
  !                                                                !
  ! usage: yspec_compare output_file reference_file tol            !
  !================================================================!
  use iso_fortran_env, only : dp => real64
  implicit none
  integer, parameter :: io1=7,io2=8
  character(len=256) :: file_out,file_ref,string
  character(len=1), dimension(3), parameter :: comp = (/'Z','N','E'/)
  integer :: ios1,ios2,n,i
  real(dp) :: tol
  real(dp), dimension(4) :: a,b
  real(dp), dimension(3) :: maxdiff,peak,err
  logical :: pass

  if(command_argument_count() /= 3) then
     error stop 'usage: yspec_compare output_file reference_file tol'
  end if
  call get_command_argument(1,file_out)
  call get_command_argument(2,file_ref)
  call get_command_argument(3,string)
  read(string,*) tol

  open(io1,file=trim(file_out),status='old',action='read',iostat=ios1)
  if(ios1 /= 0) then
     write(6,'(a)') 'FAIL '//trim(file_out)//': cannot open file'
     stop 1
  end if
  open(io2,file=trim(file_ref),status='old',action='read',iostat=ios2)
  if(ios2 /= 0) then
     write(6,'(a)') 'FAIL '//trim(file_ref)//': cannot open file'
     stop 1
  end if

  n = 0
  maxdiff = 0.0_dp
  peak = 0.0_dp
  do
     read(io1,*,iostat=ios1) a
     read(io2,*,iostat=ios2) b
     if(ios1 /= 0 .or. ios2 /= 0) exit
     n = n+1
     if(abs(a(1)-b(1)) > 1.0e-6_dp*max(1.0_dp,abs(b(1)))) then
        write(6,'(a,i0)') 'FAIL '//trim(file_out)//': time differs at line ',n
        stop 1
     end if
     maxdiff = max(maxdiff,abs(a(2:4)-b(2:4)))
     peak = max(peak,abs(b(2:4)))
  end do
  close(io1)
  close(io2)

  if(ios1 /= ios2 .or. n == 0) then
     write(6,'(a)') 'FAIL '//trim(file_out)//': files have different lengths'
     stop 1
  end if

  do i = 1,3
     if(peak(i) > 0.0_dp) then
        err(i) = maxdiff(i)/peak(i)
     else
        err(i) = maxdiff(i)
     end if
  end do
  pass = all(err <= tol)

  if(pass) then
     string = 'ok   '
  else
     string = 'FAIL '
  end if
  write(6,'(a,3(2x,a,"=",es9.2))') trim(string)//' '//trim(file_out)//':', &
       (comp(i),err(i),i=1,3)
  if(.not.pass) stop 1

end program yspec_compare
