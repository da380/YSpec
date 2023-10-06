module module_fourier



  contains

    subroutine fourrow_dp(data,isign)
      use nrtype; use nrutil, only : assert,swap
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      integer(i4b), intent(in) :: isign
      integer(i4b) :: n,i,istep,j,m,mmax,n2
      real(dp) :: theta
      complex(dpc), dimension(size(data,1)) :: temp
      complex(dpc) :: w,wp
      complex(dpc) :: ws
      n=size(data,2)
      call assert(iand(n,n-1) == 0, 'n must b a power of 2 in fourrow_dp')
      n2=n/2
      j=n2
      do i=1,n-2
         if(j > i) call swap(data(:,j+1),data(:,i+1))
         m=n2
         do
            if(m < 2 .or. j < m) exit
            j=j-m
            m=m/2
         end do
         j=j+m
      end do
      mmax=1
      do 
         if(n <= mmax) exit
         istep=2*mmax
         theta=pi_d/(isign*mmax)
         wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
         w=cmplx(1.0_dp,0.0_dp,kind=dpc)
         do m=1,mmax
            ws=w
            do i=m,n,istep
               j=i+mmax
               temp=ws*data(:,j)
               data(:,j)=data(:,i)-temp
               data(:,i)=data(:,i)+temp
            end do
            w=w*wp+w
         end do
         mmax=istep
      end do
    end subroutine fourrow_dp




    subroutine ifft(nt,i1,i2,df,w1,w2,w3,w4,dt,ur,ut,up)
      use nrtype
      implicit none
      ! input/output
      integer(i4b), intent(in) :: nt
      integer(i4b), intent(in) :: i1,i2
      real(dp), intent(in) :: df
      real(dp), intent(in) :: w1,w2,w3,w4
      real(dp), intent(in) :: dt
      complex(dpc), dimension(:,:), allocatable, intent(inout) :: ur,ut,up
      logical(lgt) :: ltmp
      ! local variables
      integer(i4b) :: n,ne,i,j,nr,k
      real(dp) :: w,filt
      complex(dp), dimension(:,:), allocatable :: data

      nr = size(ur,2)

      n = nt




      ! copy spectra into temporary array
      allocate(data(3*nr,n))

      do k = 1,nr
         do i = 1,i1-1
            data(:,i) = 0.0_dp
         end do
         do i = i1,i2
            data(1+(k-1)*3,i) = ur(i-i1+1,k)
            data(2+(k-1)*3,i) = ut(i-i1+1,k)
            data(3+(k-1)*3,i) = up(i-i1+1,k)
         end do
         do i = i2+1,n/2+1
            data(:,i) = 0.0_dp
         end do
      end do

      ! filter the spectra
      ltmp = (w1 == 0.0_dp .and. w2 == 0.0_dp & 
           .and. w3 == 0.0_dp .and. w4 == 0.0_dp)
      if(.not.ltmp) then
         do i=1,n/2+1
            w = (i-1)*df
            if(w < w1) then
               filt = 0.0_dp
            else if(w >= w1 .and. w < w2) then
               filt = pi_d*(w-w1)/(w2-w1)
               filt = 0.5_dp*(1.0_dp-cos(filt))
            else if(w >= w2 .and. w < w3) then
               filt = 1.0_dp
            else if(w >= w3 .and. w < w4) then
               filt = pi_d*(w4-w)/(w4-w3)
               filt = 0.5_dp*(1.0_dp-cos(filt))
            else if(w >= w4) then
               filt = 0.0_dp
            end if
            data(:,i) = filt*data(:,i)
         end do
      end if


      ! do the negative frequencies
      j = 0
      do i = n/2+2,n
         j = j+1
         data(:,i) = conjg(data(:,n/2-j+1))
      end do
      
               

      ! do the inverse Fourier transform
      call fourrow_dp(data,1)

      data = data/(dt*n)

      ! put seismograms back into the displacement
      ! vector arrays
      deallocate(ur,ut,up)
      allocate(ur(n,nr),ut(n,nr),up(n,nr))
      do k = 1,nr
         ur(:,k) = data(1+(k-1)*3,:)
         ut(:,k) = data(2+(k-1)*3,:)
         up(:,k) = data(3+(k-1)*3,:)
      end do
      deallocate(data)

      return
    end subroutine ifft


    subroutine ifft4(nt,i1,i2,df,w1,w2,w3,w4,dt,ur,ut,up,phi)
      use nrtype
      implicit none
      ! input/output
      integer(i4b), intent(in) :: nt
      integer(i4b), intent(in) :: i1,i2
      real(dp), intent(in) :: df
      real(dp), intent(in) :: w1,w2,w3,w4
      real(dp), intent(in) :: dt
      complex(dpc), dimension(:,:), allocatable, intent(inout) :: ur
      complex(dpc), dimension(:,:), allocatable, intent(inout) :: ut
      complex(dpc), dimension(:,:), allocatable, intent(inout) :: up
      complex(dpc), dimension(:,:), allocatable, intent(inout) :: phi
      logical(lgt) :: ltmp
      ! local variables
      integer(i4b) :: n,ne,i,j,nr,k
      real(dp) :: w,filt
      complex(dp), dimension(:,:), allocatable :: data

      nr = size(ur,2)
      n = nt

      ! copy spectra into temporary array
      allocate(data(4*nr,n))

      do k = 1,nr
         do i = 1,i1-1
            data(:,i) = 0.0_dp
         end do
         do i = i1,i2
            data(1+(k-1)*4,i) = ur(i-i1+1,k)
            data(2+(k-1)*4,i) = ut(i-i1+1,k)
            data(3+(k-1)*4,i) = up(i-i1+1,k)
            data(4+(k-1)*4,i) = phi(i-i1+1,k)
         end do
         do i = i2+1,n/2+1
            data(:,i) = 0.0_dp
         end do
      end do

      ! filter the spectra
      ltmp = (w1 == 0.0_dp .and. w2 == 0.0_dp & 
           .and. w3 == 0.0_dp .and. w4 == 0.0_dp)
      if(.not.ltmp) then
         do i=1,n/2+1
            w = (i-1)*df
            if(w < w1) then
               filt = 0.0_dp
            else if(w >= w1 .and. w < w2) then
               filt = pi_d*(w-w1)/(w2-w1)
               filt = 0.5_dp*(1.0_dp-cos(filt))
            else if(w >= w2 .and. w < w3) then
               filt = 1.0_dp
            else if(w >= w3 .and. w < w4) then
               filt = pi_d*(w4-w)/(w4-w3)
               filt = 0.5_dp*(1.0_dp-cos(filt))
            else if(w >= w4) then
               filt = 0.0_dp
            end if
            data(:,i) = filt*data(:,i)
         end do
      end if


      ! do the negative frequencies
      j = 0
      do i = n/2+2,n
         j = j+1
         data(:,i) = conjg(data(:,n/2-j+1))
      end do
      
               

      ! do the inverse Fourier transform
      call fourrow_dp(data,1)

      data = data/(dt*n)

      ! put seismograms back into the displacement
      ! vector arrays
      deallocate(ur,ut,up,phi)
      allocate(ur(n,nr),ut(n,nr),up(n,nr),phi(n,nr))
      do k = 1,nr
         ur(:,k)  = data(1+(k-1)*4,:)
         ut(:,k)  = data(2+(k-1)*4,:)
         up(:,k)  = data(3+(k-1)*4,:)
         phi(:,k) = data(4+(k-1)*4,:)
      end do
      deallocate(data)

      return
    end subroutine ifft4





    subroutine ifft_phi(nt,i1,i2,df,w1,w2,w3,w4,dt,phi_coef,phi_coef_time)
      use nrtype
      implicit none
      ! input/output
      integer(i4b), intent(in) :: nt
      integer(i4b), intent(in) :: i1,i2
      real(dp), intent(in) :: df
      real(dp), intent(in) :: w1,w2,w3,w4
      real(dp), intent(in) :: dt
      complex(dpc), dimension(:,:), intent(in) :: phi_coef
      complex(dpc), dimension(:,:), intent(out) :: phi_coef_time
      logical(lgt) :: ltmp
      ! local variables
      integer(i4b) :: n,ne,i,j,k
      real(dp) :: w,filt
      complex(dp), dimension(:,:), allocatable :: data



      n = nt




      ! copy spectra into temporary array
      allocate(data(5,n))

      do i = 1,i1-1
         data(:,i) = 0.0_dp
      end do
      do i = i1,i2
            data(:,i) = phi_coef(i-i1+1,:)
      end do
      do i = i2+1,n/2+1
         data(:,i) = 0.0_dp
      end do


      ! filter the spectra
      ltmp = (w1 == 0.0_dp .and. w2 == 0.0_dp & 
           .and. w3 == 0.0_dp .and. w4 == 0.0_dp)
      if(.not.ltmp) then
         do i=1,n/2+1
            w = (i-1)*df
            if(w < w1) then
               filt = 0.0_dp
            else if(w >= w1 .and. w < w2) then
               filt = pi_d*(w-w1)/(w2-w1)
               filt = 0.5_dp*(1.0_dp-cos(filt))
            else if(w >= w2 .and. w < w3) then
               filt = 1.0_dp
            else if(w >= w3 .and. w < w4) then
               filt = pi_d*(w4-w)/(w4-w3)
               filt = 0.5_dp*(1.0_dp-cos(filt))
            else if(w >= w4) then
               filt = 0.0_dp
            end if
            data(:,i) = filt*data(:,i)
         end do
      end if


      ! do the negative frequencies
      j = 0
      do i = n/2+2,n
         j = j+1
         data(:,i) = conjg(data(:,n/2-j+1))
      end do
      
               

      ! do the inverse Fourier transform
      call fourrow_dp(data,1)

      data = data/(dt*n)

      
      do k = 1,5
         phi_coef_time(:,k) = data(k,:)
      end do
      deallocate(data)

      return
    end subroutine ifft_phi
    



  function hann(t,t11,t12,t21,t22)
    use nrtype
    implicit none
    logical(lgt) :: ltmp
    real(dp) :: hann
    real(dp), intent(in) :: t,t11,t12,t21,t22

    ltmp = (t11 == 0.0_dp .and. t12 == 0.0_dp & 
         .and. t21 == 0.0_dp .and. t22 == 0.0_dp)
    if(.not.ltmp) then
       if(t < t11) then
          hann = 0.0_dp
       else if(t >= t11 .and. t < t12) then
          hann = pi_d*(t-t11)/(t12-t11)
          hann = 0.5_dp*(1.0_dp-cos(hann))
       else if(t >= t12 .and. t < t21) then
          hann = 1.0_dp
       else if(t >= t21 .and. t < t22) then
          hann = pi_d*(t22-t)/(t22-t21)
          hann = 0.5_dp*(1.0_dp-cos(hann))
       else if(t >= t22) then
          hann = 0.0_dp
       end if       
    end if
    return
  end function hann



  subroutine fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
    use nrtype
    implicit none
    real(dp), intent(inout) :: f1,f2
    real(dp), intent(in) :: dt,tout
    integer(i4b), intent(in) :: mex,qex
    real(dp), intent(out) :: df,ep
    integer(i4b), intent(out) :: nt,i1,i2
    
    integer(i4b) :: ne
    real(dp) :: fn
    
    ! check that the Nyquist frequency for the given
    ! value of dt lies above f2
    fn = 0.5_sp/dt
    if(fn < f2) stop
        
    ep = mex/tout
    
    df = ep/(twopi_d*qex)
    
    nt = 1.0_sp/(df*dt)
    
    ne = log(real(nt))/log(2.0_sp)+1
    nt = 2**ne
    
    df = 1.0_sp/(nt*dt)
    
    i1 = max(floor(f1/df),2)
    f1 = (i1-1)*df          
    
    i2 = floor(f2/df)+2
    f2 = (i2-1)*df

    return
  end subroutine fcal



end module module_fourier
