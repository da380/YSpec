module module_start

  use nrtype
  implicit none
  
  real(dp), parameter :: tol_start = 11

  
contains




  subroutine tor_start_level(l,w,it,up,is)
    use nrtype
    use module_model
    implicit none
    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: w
    integer(i4b), intent(in) :: it,up
    integer(i4b), intent(out) :: is
    
    if(up == 1) then
       call start_level(l,w,noc+1,it,1,is)
    else if(up == -1) then
       call start_level(l,w,it,noc+1,1,is)
    else
       stop 'tor_start_level: up should equal 1 or -1'
    end if
    
    return
  end subroutine tor_start_level


  subroutine sph_start_level(l,w,it,up,is)
     use nrtype
    use module_model
    implicit none
    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: w
    integer(i4b), intent(in) :: it,up
    integer(i4b), intent(out) :: is
    integer(i4b) :: isp,iss
    real(dp), parameter :: aw = -2.00e-3_dp
    real(dp), parameter :: bw =  2.25e-3_dp
    real(dp), parameter :: dw =  1.28e-3_dp
    real(dp) :: wsoc,wsic,wdim,sum


 
    wdim = w*fre_norm
    wsoc = aw+dw*l
    if(wdim <= wsoc) then
       call start_level(l,w,noc+1,it,1,is)
       if(is == nsl) is = is-1
       if(is > noc+1) return
    end if
    wsic = aw+bw*l
    if(wdim <= wsic) then
       call start_level(l,w,nic+1,noc,2,is)
       if(is == noc) is = is-1
       if(is > nic+1) return
    end if
    if(fic) then
       call start_level(l,w,2,nic,2,is)
    else
       call start_level(l,w,2,nic,1,is)
    end if
    if(is == nic) is = is-1
 



    return
  end subroutine sph_start_level




  subroutine start_level(l,w,i1,i2,wt,is,sum_in)
    use nrtype
    use module_model
    implicit none

    ! input/output
    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: w
    integer(i4b), intent(in) :: i1,i2,wt
    integer(i4b), intent(out) :: is
    real(dp), intent(inout), optional  :: sum_in

    ! local variables
    integer(i4b) :: i,j,isign,j1,j2
    real(dp) :: pp,q1,q2,sum

    ! search direction up or down
    if(i2 >= i1) then
       isign =  1
    else
       isign = -1 
    end if

    ! compute squared ray parameter
    pp = real(l*(l+1))/(w*w)


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

    if(present(sum_in)) then
       sum = sum_in
    else
       sum = 0.0_dp
    end if
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
    
    if(present(sum_in)) sum_in = sum

    return
  end subroutine start_level





  

end module module_start
