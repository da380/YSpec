module module_fourier
  !================================================================!
  ! Routines for setting up the frequency discretisation, and for  !
  ! filtering the spectra and transforming them to the time        !
  ! domain using FFTW.                                             !
  !================================================================!
  use, intrinsic :: iso_c_binding
  use module_kinds
  implicit none
  include 'fftw3.f03'
  private
  public :: ifft,fcal

contains


  subroutine ifft(nt,i1,i2,df,w1,w2,w3,w4,dt,ur,ut,up,seis)
    !------------------------------------------------------------!
    ! The spectra ur, ut and up hold, for each receiver, the     !
    ! frequencies with indices i1 to i2 of the nt/2+1            !
    ! non-negative frequencies (i-1)*df. A cosine high-pass      !
    ! filter between w1 and w2 and a cosine low-pass filter      !
    ! between w3 and w4 are applied unless all four are zero,    !
    ! and the spectra are then transformed to time series with   !
    ! nt samples spaced by dt. On output, seis(:,1:3,k) holds    !
    ! the three components of the time series at receiver k.     !
    !------------------------------------------------------------!
    integer(i4b), intent(in) :: nt,i1,i2
    real(dp), intent(in) :: df,w1,w2,w3,w4,dt
    complex(dpc), dimension(:,:), intent(in) :: ur,ut,up
    real(dp), dimension(:,:,:), allocatable, intent(out) :: seis
    integer(i4b) :: nr,nh,i,j,k
    integer(c_int), dimension(1) :: nn,inembed,onembed
    real(dp) :: w,filt
    complex(c_double_complex), dimension(:,:), allocatable :: spec
    real(c_double), dimension(:,:), allocatable :: series
    type(c_ptr) :: plan

    nr = size(ur,2)
    nh = nt/2+1
    allocate(spec(nh,3*nr),series(nt,3*nr))

    ! copy the spectra into the array to be transformed; any
    ! frequency above the Nyquist frequency is dropped
    spec = 0.0_dp
    do k = 1,nr
       do i = i1,min(i2,nh)
          spec(i,3*k-2) = ur(i-i1+1,k)
          spec(i,3*k-1) = ut(i-i1+1,k)
          spec(i,3*k)   = up(i-i1+1,k)
       end do
    end do

    ! filter the spectra
    if(.not.(w1 == 0.0_dp .and. w2 == 0.0_dp .and. &
             w3 == 0.0_dp .and. w4 == 0.0_dp)) then
       do i = 1,nh
          w = (i-1)*df
          if(w < w1) then
             filt = 0.0_dp
          else if(w < w2) then
             filt = 0.5_dp*(1.0_dp-cos(pi_d*(w-w1)/(w2-w1)))
          else if(w < w3) then
             filt = 1.0_dp
          else if(w < w4) then
             filt = 0.5_dp*(1.0_dp-cos(pi_d*(w4-w)/(w4-w3)))
          else
             filt = 0.0_dp
          end if
          spec(i,:) = filt*spec(i,:)
       end do
    end if

    ! transform all the components at once; the complex-to-real
    ! transform makes use of the conjugate symmetry of the spectra
    nn = int(nt,c_int)
    inembed = int(nh,c_int)
    onembed = int(nt,c_int)
    plan = fftw_plan_many_dft_c2r(1_c_int,nn,int(3*nr,c_int),         &
                                  spec,inembed,1_c_int,int(nh,c_int), &
                                  series,onembed,1_c_int,int(nt,c_int), &
                                  FFTW_ESTIMATE)
    call fftw_execute_dft_c2r(plan,spec,series)
    call fftw_destroy_plan(plan)

    ! normalise and store the time series
    allocate(seis(nt,3,nr))
    do k = 1,nr
       do j = 1,3
          seis(:,j,k) = series(:,3*(k-1)+j)/(dt*nt)
       end do
    end do
    deallocate(spec,series)

    return
  end subroutine ifft


  subroutine fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
    !------------------------------------------------------------!
    ! Sets up the frequency discretisation. The imaginary part   !
    ! ep of the complex frequency is chosen so that the response !
    ! is damped by the factor exp(-mex) over the length tout of  !
    ! the time series. The number of samples nt is the smallest  !
    ! power of two larger than 1/(df*dt) with df = ep/(2*pi*qex),!
    ! and df is then reset to 1/(nt*dt). On output f1 and f2 are !
    ! rounded to the frequency grid, and i1 and i2 are their     !
    ! indices, with zero frequency at index 1.                   !
    !------------------------------------------------------------!
    real(dp), intent(inout) :: f1,f2
    real(dp), intent(in) :: dt,tout
    integer(i4b), intent(in) :: mex,qex
    real(dp), intent(out) :: df,ep
    integer(i4b), intent(out) :: nt,i1,i2
    integer(i4b) :: nmin

    ! check that the Nyquist frequency lies above f2
    if(0.5_dp/dt < f2) then
       error stop 'the Nyquist frequency is below fmax: reduce the time step'
    end if

    ep = mex/tout
    df = ep/(twopi_d*qex)
    nmin = int(1.0_dp/(df*dt))
    nt = 1
    do while(nt <= nmin)
       nt = 2*nt
    end do
    df = 1.0_dp/(nt*dt)

    i1 = max(floor(f1/df),2)
    f1 = (i1-1)*df
    i2 = floor(f2/df)+2
    f2 = (i2-1)*df

    return
  end subroutine fcal


end module module_fourier
