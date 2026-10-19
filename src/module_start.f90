module module_start
  !================================================================!
  ! Routines that choose the radius at which to start the radial   !
  ! integrations. Below the turning point of a given (l,w) the     !
  ! solutions are evanescent, so the integration can start at the  !
  ! depth where a WKB estimate of the solution has decayed by the  !
  ! factor exp(-tol_start) from its value at the turning point.    !
  !================================================================!
  use module_kinds
  implicit none

  ! With tol_start = 15 the seismograms differ from those obtained by
  ! always starting at the base of the model by about 1e-8 of their
  ! peak values (1e-6 for the earlier value of 11), at a cost of a
  ! few percent in run time.
  real(dp), parameter :: tol_start = 15.0_dp

contains


  subroutine tor_start_level(l,w,it,up,is)
    !------------------------------------------------------------!
    ! Returns the index is of the knot at which to start the     !
    ! integration of the toroidal equations, given the index it  !
    ! of the knot below the source or receiver. The search is    !
    ! upwards from the core-mantle boundary if up = 1, and       !
    ! downwards to it if up = -1.                                !
    !------------------------------------------------------------!
    use module_model
    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: w
    integer(i4b), intent(in) :: it,up
    integer(i4b), intent(out) :: is

    if(up == 1) then
       call start_level(l,w,noc+1,it,1,is)
    else if(up == -1) then
       call start_level(l,w,it,noc+1,1,is)
    else
       error stop 'tor_start_level: up should equal 1 or -1'
    end if

    return
  end subroutine tor_start_level


  subroutine sph_start_level(l,w,it,up,is)
    !------------------------------------------------------------!
    ! Returns the index is of the knot at which to start the     !
    ! integration of the spheroidal equations, given the index   !
    ! it of the knot below the source or receiver. The mantle is !
    ! tried first, then the outer core and then the inner core,  !
    ! stopping in the first region where the solution decays     !
    ! sufficiently. The argument up is not used.                 !
    !------------------------------------------------------------!
    use module_model
    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: w
    integer(i4b), intent(in) :: it,up
    integer(i4b), intent(out) :: is

    ! start within the mantle
    call start_level(l,w,noc+1,it,1,is)
    if(is == nsl) is = is-1
    if(is > noc+1) return

    ! start within the outer core
    call start_level(l,w,nic+1,noc,2,is)
    if(is == noc) is = is-1
    if(is > nic+1) return

    ! start within the inner core
    if(fic) then
       call start_level(l,w,2,nic,2,is)
    else
       call start_level(l,w,2,nic,1,is)
    end if
    if(is == nic) is = is-1

    return
  end subroutine sph_start_level


  subroutine start_level(l,w,i1,i2,wt,is)
    !------------------------------------------------------------!
    ! Searches the knots i1 to i2 for the turning point of the   !
    ! S waves (wt = 1) or the P waves (wt = 2) at degree l and   !
    ! frequency w, and then moves down from it until the WKB     !
    ! estimate of the solution has decayed by exp(-tol_start).   !
    ! Returns the index of that knot, or of the bottom of the    !
    ! interval if the decay is not reached within it.            !
    !------------------------------------------------------------!
    use module_model
    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: w
    integer(i4b), intent(in) :: i1,i2,wt
    integer(i4b), intent(out) :: is
    integer(i4b) :: i,j,isign,j1,j2
    real(dp) :: pp,q1,q2,sum

    ! search direction up or down
    if(i2 >= i1) then
       isign =  1
    else
       isign = -1
    end if

    ! compute squared ray parameter
    pp = real(l*(l+1),dp)/(w*w)

    ! look for a turning point
    j = i2
    do i = i1,i2,isign
       if(wt == 1) then
          q1 = 1.0_dp/(beta(i)*beta(i))-pp/(r(i)*r(i))
       else
          q1 = 1.0_dp/(alpha(i)*alpha(i))-pp/(r(i)*r(i))
       end if
       if(isign == 1) then
          if(q1 > 0.0_dp) then
             j = max(i1,i-1)
             exit
          end if
       else
          if(q1 < 0.0_dp) then
             j = i
             exit
          end if
       end if
    end do
    is = j

    ! return if turning point is at the
    ! base of the interval
    if(is == min(i1,i2)) return

    ! move down from the turning point  so that
    ! the asymptotic solution  will decay by the
    ! factor exp(-tol_start) from its value at the
    ! turning point
    j1 = is
    j2 = min(i1,i2)
    is = j2
    sum = 0.0_dp
    do i = j1,j2+1,-1
       if(wt == 1) then
          q1 = 1.0_dp/(beta(i)*beta(i))-pp/(r(i)*r(i))
          q2 = 1.0_dp/(beta(i-1)*beta(i-1))-pp/(r(i-1)*r(i-1))
       else
          q1 = 1.0_dp/(alpha(i)*alpha(i))-pp/(r(i)*r(i))
          q2 = 1.0_dp/(alpha(i-1)*alpha(i-1))-pp/(r(i-1)*r(i-1))
       end if
       ! if solution becomes oscillatory again
       ! re-initialise the sum, else add the
       ! term to the sum
       if(q1 > 0.0_dp .or. q2 > 0.0_dp) then
          sum = 0.0_dp
       else
          q1 = sqrt(-q1)
          q2 = sqrt(-q2)
          sum = sum+0.5_dp*w*(q1+q2)*(r(i)-r(i-1))
       end if
       ! if the sum has become sufficiently large
       ! the starting point has been found
       if(sum > tol_start) then
          is = i-1
          exit
       end if
    end do

    return
  end subroutine start_level


end module module_start
