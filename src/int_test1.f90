program int_test1

  use nrtype
  use module_util
  use module_model
  implicit none

  character(len=256) :: model_file

  integer(i4b), parameter :: io1 = 7
  integer(i4b), parameter :: io2 = 8

  integer(i4b) :: ios,l

  real(dp) :: wr
  complex(dpc) :: w

  ! set normalization paramerters
  call set_parameters

  
  ! get the model file
  call get_string(' model file = ',model_file)
  open(io2,file=model_file,form='formatted', & 
       action='read',iostat=ios)
  if(ios /= 0) stop ' problem opening model file'
  call  read_model(io2)
  close(io2)


  ! get value of l
  call get_integer(' l = ',l)
  
  ! get omega
  call get_float(' frequency (mHz) = ',wr)
  wr = wr/fre_norm*pi_d/500.0_dp

  w = wr

  open(io1,file='int_test1.out')

  call  sph_int(io1,l,w)
  
  close(io1)


  contains



    subroutine sph_int(io,ll_in,w)
      use nrtype
      use module_model
      use module_int
      implicit none

      integer(i4b), intent(in) :: io
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w


      integer(i4b) :: i,i1,i2,n,j,k

      real(dp) :: r1,rr1,r2,rr2,drr,rr

      complex(dpc) :: m20
      complex(dpc), dimension(6,3) :: y
      complex(dpc), dimension(6,3) :: dydr
      complex(dpc), dimension(14) :: m
      complex(dpc), dimension(14) :: dmdr
      

      call set_l(ll_in)
      call set_ats(0.0_dp)
      call set_omega(w)
      call sph_steps





      ! set initial conditions at the cmb
      y = 0.0_dp
      do i = 1,3
         y(i,i) = 1.0_dp
      end do
      m = 0.0_dp
      m(1) = 1.0_dp

      i1 = noc+1
      i2 = nsl

      r1 = r(i1)
      r2 = r(i2)


      do i=i1,i2
         if(i == i1) then
            rr1=r1
            else
               rr1=r(i)
            end if
            if(i == i2) then
               rr2=r2
            else
               rr2=r(i+1)
            end if
            if(rr1 == rr2) cycle
            call set_layer(i)
            rr=rr1
            drr = steps(i)
            n = 2+floor((rr2-rr1)/drr)
            drr = (rr2-rr1)/(n-1)
            do j = 1,n-1
               rr = rr1+(j-1)*drr

               ! update the non-minor system
               do k = 1,3
                  call sph_solid_derivs(rr,y(:,k),dydr(:,k))
                  call rkck(y(:,k),dydr(:,k),rr,drr, & 
                       y(:,k),sph_solid_derivs)        
               end do

               ! update the minor system
               call sph_minors_solid_derivs(rr,m,dmdr)
               call rkck(m,dmdr,rr,drr,m,sph_minors_solid_derivs)   

               ! compute the m_{20} minor component from Y
               m20 = det3(y((/4,5,6/),1),y((/4,5,6/),2),y((/4,5,6/),3))

               ! write out the two versions of m_{20}
               write(io,*) rr*r_norm/1000.0_dp,real(m20),real(m(14))

            end do
         end do
        


      return
    end subroutine sph_int



end program int_test1
