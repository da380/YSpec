module module_util

  ! routines for reading the parameter file and naming output files

  contains


    subroutine read_string(io,string,skip)
      use module_kinds
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
      use module_kinds
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
      use module_kinds
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


    subroutine string_cat_int(stri,i,stro)
      use module_kinds
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



end module module_util
