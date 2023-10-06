module module_util

  use nrtype
  
  interface my_reallocate
     module procedure my_reallocate_iv,my_reallocate_ia,my_reallocate_rv, &
          my_reallocate_dv,my_reallocate_cv,my_reallocate_zv
  end interface

  interface poly_interp
     module procedure poly_interp_r, poly_interp_z
  end interface

  interface get_float
     module procedure get_float_sp, get_float_dp
  end interface get_float

  interface count_columns
     procedure :: count_columns,count_columns_opened
  end interface count_columns


  interface found_command_argument
     procedure :: found_command_argument_character, &
                  found_command_argument_integer,   &
                  found_command_argument_real,      &
                  found_command_argument_logical,   &
                  found_command_argument_flag 
  end interface found_command_argument

    ! default string length
  integer(i4b), parameter :: max_string_length = 4096
  
  contains


    subroutine my_reallocate_iv(a,n)
      use nrtype
      implicit none
      integer(i4b), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      integer(i4b), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_iv


    subroutine my_reallocate_ia(a,n,m)
      use nrtype
      implicit none
      integer(i4b), dimension(:,:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n,m
      integer(i4b), dimension(:,:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store,m_tmp,m_store
      n_tmp=size(a,1); m_tmp=size(a,2)
      n_store=min(n,n_tmp); m_store=min(m,m_tmp)      
      allocate(a_tmp(n_store,m_store))
      a_tmp(1:n_store,1:m_store)=a(1:n_store,1:m_store)
      deallocate(a); allocate(a(n,m))
      a(1:n_store,1:m_store)=a_tmp(1:n_store,1:m_store)
      deallocate(a_tmp)
      return
    end subroutine my_reallocate_ia



    subroutine my_reallocate_rv(a,n)
      use nrtype
      implicit none
      real(sp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      real(sp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_rv


    subroutine my_reallocate_dv(a,n)
      use nrtype
      implicit none
      real(dp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      real(dp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_dv
      

    subroutine my_reallocate_cv(a,n)
      use nrtype
      implicit none
      complex(sp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      complex(sp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_cv


    subroutine my_reallocate_zv(a,n)
      use nrtype
      implicit none
      complex(dpc), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      complex(dpc), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_zv


    subroutine poly_interp_r(x1,f1,fp1,x2,f2,fp2,x,f,fp)
      use nrtype
      implicit none
      real(dp), intent(in) :: x,x1,x2
      real(dp), intent(in) :: f1,fp1,f2,fp2
      real(dp), intent(out) :: f,fp
      real(dp) :: c2,c3,delta,deltai,deltai2

      if(x == x1) then
        f=f1
        fp=fp1
        return
      else if(x == x2) then
         f=f2
         fp=fp2
         return
      end if

      delta=x2-x1

      if(delta == 0) stop 'poly_interp: delta = 0'
      
      deltai=1.0_dp/delta
      deltai2=deltai*deltai

      c2=3.0_dp*deltai2*(f2-f1)-deltai*(fp2+2.0_dp*fp1)
      c3=deltai2*(fp2+fp1)-2.0_dp*deltai*deltai2*(f2-f1)

      delta=x-x1
      deltai=delta*delta
      deltai2=deltai*delta

      f=f1+fp1*delta+c2*deltai+c3*deltai2
      fp=fp1+2.0_dp*c2*delta+3.0_dp*c3*deltai


      return
    end subroutine poly_interp_r


    subroutine poly_interp_z(x1,f1,fp1,x2,f2,fp2,x,f,fp)
      use nrtype
      implicit none
      real(dp), intent(in) :: x,x1,x2
      complex(dpc), intent(in) :: f1,fp1,f2,fp2
      complex(dpc), intent(out) :: f,fp
      real(dp) :: delta,deltai,deltai2
      complex(dpc) :: c2,c3

      if(x == x1) then
        f=f1
        fp=fp1
        return
      else if(x == x2) then
         f=f2
         fp=fp2
         return
      end if

      delta=x2-x1

      if(delta == 0) stop 'poly_interp: delta = 0'
      
      deltai=1.0_dp/delta
      deltai2=deltai*deltai

      c2=3.0_dp*deltai2*(f2-f1)-deltai*(fp2+2.0_dp*fp1)
      c3=deltai2*(fp2+fp1)-2.0_dp*deltai*deltai2*(f2-f1)

      delta=x-x1
      deltai=delta*delta
      deltai2=deltai*delta

      f=f1+fp1*delta+c2*deltai+c3*deltai2
      fp=fp1+2.0_dp*c2*delta+3.0_dp*c3*deltai


      return
    end subroutine poly_interp_z


    subroutine get_string(tag,string)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      character(len=*), intent(out) :: string

      write(6,'("'//tag//'")',advance='no')
      read(5,*) string
      string=adjustl(string)
      return
    end subroutine get_string

    
    subroutine get_integer(tag,int)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      integer(i4b), intent(out) :: int
      character(len=256) :: string,form,sls
      integer(i4b) :: sl,ios
      
      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         write(sls,*) sl
         sls=adjustl(sls)
         sl=len_trim(sls)
         form='(i'//sls(1:sl)//')'
         read(string,form,iostat=ios) int
         if(ios == 0) exit
         print *, ' input must be an integer'
      end do
            
      return
    end subroutine get_integer


    subroutine get_float_dp(tag,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      real(dp), intent(out) :: fpn
      character(len=256) :: string,form,sls,ips
      integer(i4b) :: sl,ios,ip

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         ip=index(string,'.')
         if(ip == 0) ip=sl
         write(sls,*) sl
         write(ips,*) sl-ip
         sls=adjustl(sls)
         sl=len_trim(sls)
         ips=adjustl(ips)
         ip=len_trim(ips)
         form='(f'//sls(1:sl)//'.'//ips(1:ip)//')'
         read(string,form,iostat=ios) fpn
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_float_dp

    subroutine get_float_sp(tag,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      real(sp), intent(out) :: fpn
      character(len=256) :: string,form,sls,ips
      integer(i4b) :: sl,ios,ip

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         ip=index(string,'.')
         if(ip == 0) ip=sl
         write(sls,*) sl
         write(ips,*) sl-ip
         sls=adjustl(sls)
         sl=len_trim(sls)
         ips=adjustl(ips)
         ip=len_trim(ips)
         form='(f'//sls(1:sl)//'.'//ips(1:ip)//')'
         read(string,form,iostat=ios) fpn
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_float_sp


    subroutine read_string(io,string,skip)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: io
      character(len=*), intent(out) :: string
      integer(i4b), intent(in), optional :: skip
      integer(i4b) :: i

      if(present(skip)) then
         do i = 1,skip
            read(io,*) 
         end do
      end if
      read(io,*) string
      return
    end subroutine read_string


    subroutine read_int(io,j,skip)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: io
      integer(i4b), intent(out) :: j
      integer(i4b), intent(in), optional :: skip
      integer(i4b) :: i

      if(present(skip)) then
         do i = 1,skip
            read(io,*) 
         end do
      end if
      read(io,*) j
      return
    end subroutine read_int


 subroutine read_float(io,x,skip)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: io
      real(dp), intent(out) :: x
      integer(i4b), intent(in), optional :: skip
      integer(i4b) :: i

      if(present(skip)) then
         do i = 1,skip
            read(io,*) 
         end do
      end if
      read(io,*) x
      return
    end subroutine read_float




    function string_length(string,sign)
      use nrtype
      implicit none
      integer(i4b) :: string_length
      integer(i4b), intent(in) :: sign
      character (len=*), intent(in) :: string
      select case(sign)
      case(1)
         string_length=0
         do 
            if(string(string_length+sign:string_length+sign) == ' ') exit
            string_length=string_length+sign
         end do
      case(-1)
         string_length=len(string)
         do 
            if(string(string_length:string_length) == ' ') exit
            string_length=string_length+sign
         end do
     case default
        stop 'bad input of sign to string_length'
     end select
     return
   end function string_length


    subroutine bound(xx,x,i1,i2,err)
      ! Given a monotonically increasing array, xx, and a value
      ! x, returns the indices i1 and i2 such that
      ! xx(i1) <= x <= xx(i2). The routine checks that x is in 
      ! range, and if not returns err=.true.
      use nrtype
      implicit none
      real(dp), dimension(:), intent(in) :: xx
      real(dp), intent(in) :: x
      integer(i4b), intent(out) :: i1,i2
      logical(lgt), intent(out) :: err
      integer(i4b) :: i,n
      real(dp), parameter :: tol=1.0e-5_dp
      real(dp) :: rtol
      
      err=.false.
      n=size(xx)
      rtol=tol*xx(n)
      if(x > xx(n) .and. x < xx(1)) then
         err=.true.
         return
      end if

      if(abs(x-xx(n)) < rtol) then
         i2=n
         i1=n-1
         return
      end if
      if(abs(x-xx(1)) < rtol) then
         i2=2
         i1=1
         return
      end if
      do i=n-1,1,-1
         if(x <= xx(i+1) .and. x >= xx(i)) then
            i2=i+1
            i1=i
            return
         end if
      end do
      err=.true.
      return
    end subroutine bound


    subroutine string_cat_int(stri,i,stro)
      use nrtype
      implicit none
      
      character(len=*), intent(in) :: stri
      integer(i4b), intent(in) :: i
      character(len=*), intent(out) :: stro

      character(len=10) :: stmp
      integer(i4b) :: j

      write(stmp,'(i10)') i
      j = floor(log10(real(i)))
      stro = trim(stri)//stmp(10-j:10)
      
      return
    end subroutine string_cat_int


    !==============================================!
    !                  IO routines                 !
    !==============================================!
  
  function count_columns(file) result(ncol)
    implicit none
    character(len = *), intent(in) :: file
    integer(i4b) :: ncol

    character(len=:), allocatable :: line
    integer(i4b) :: io,ios
    real(dp), dimension(:), allocatable :: tmp
    open(newunit = io,file = trim(file))
    line = readline(io)
    allocate(tmp(len(line)))
    close(io)
    ncol = 0
    do
       read(line,*,iostat = ios) tmp(1:ncol)
       if(ios /= 0) exit
       ncol = ncol + 1
    end do 
    ncol = ncol-1    
    return
  end function count_columns

  function count_columns_opened(io) result(ncol)
    implicit none
    integer(i4b), intent(in) :: io
    integer(i4b) :: ncol

    logical :: ok
    character(len=:), allocatable :: line
    integer(i4b) :: ios
    real(dp), dimension(:), allocatable :: tmp
    line = readline(io)
    allocate(tmp(len(line)))
    backspace(io)
    ncol = 0
    do
       read(line,*,iostat = ios) tmp(1:ncol)
       if(ios /= 0) exit
       ncol = ncol + 1
    end do 
    ncol = ncol-1
    
    return
  end function count_columns_opened

  
  function readline(io) result(line)
    implicit none
    integer, intent(IN) :: io
    character(len=:), allocatable :: line
    integer, parameter :: buf_len= max_string_length
    character(len=buf_len) :: buf
    logical :: okay, set
    integer status, size    
    okay = .false.
    set = .true.
    do
       read (io,'(a)',advance='NO',iostat=status, size=size) buf
       okay = .not. IS_IOSTAT_END(status)
       if (.not. okay) return
       if (set) then
          line = buf(1:size)
          set=.false.
       else
          line = line // buf(1:size)
       end if
       if (IS_IOSTAT_EOR(status)) exit
    end do
  end function readline


  !==============================================!
  !             command line routines            !
  !==============================================!

  subroutine print_argument_info(nreq,ireq,dreq,nopt,iopt,dopt) 
    integer(i4b), intent(in) :: nreq
    character(len=*), dimension(nreq), intent(in) :: ireq
    logical, dimension(nreq), intent(in) :: dreq
    integer(i4b), intent(in), optional :: nopt
    character(len=*), dimension(:), intent(in), optional :: iopt
    logical, dimension(:), intent(in), optional :: dopt

    integer(i4b) :: i,nc
    
    ! work out the minimum number of command line arguments
    nc = 0
    do i = 1,nreq
       if(dreq(i)) then
          nc = nc+2
       else
          nc = nc+1
       end if
    end do

    if(command_argument_count() < nc) then
       do i = 1,nreq
          print *, trim(ireq(i))
       end do
       if(present(nopt) .and. present(iopt) .and. present(dopt)) then
          do i = 1,nopt
             print *, trim(iopt(i))
          end do
       end if
       stop
    end if
  end subroutine print_argument_info

  logical function found_command_argument_flag(tag) result(found)
    character(len=*), intent(in) :: tag
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          found = .true.          
          return
       end if
    end do
    return
  end function found_command_argument_flag

  
  logical function found_command_argument_character(tag,val) result(found)
    character(len=*), intent(in) :: tag
    character(len=:), allocatable, intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          found = .true.
          val = trim(string)         
       end if
    end do
    return
  end function found_command_argument_character

  logical function found_command_argument_integer(tag,val) result(found)
    character(len=*), intent(in) :: tag
    integer(i4b), intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          read(string,*,iostat = ios) val
          if(ios /= 0)  return
          found = .true.
          return
       end if
    end do
    return
  end function found_command_argument_integer


  logical function found_command_argument_real(tag,val) result(found)
    character(len=*), intent(in) :: tag
    real(dp), intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          read(string,*,iostat = ios) val
          if(ios /= 0)  return
          found = .true.
          return
       end if
    end do
    return
  end function found_command_argument_real


  logical function found_command_argument_logical(tag,val) result(found)
    character(len=*), intent(in) :: tag
    logical, intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          read(string,*,iostat = ios) val
          if(ios /= 0)  return
          found = .true.
          return
       end if
    end do
    return
  end function found_command_argument_logical
  



end module module_util
