module module_spline
  !================================================================!
  ! Natural cubic splines of real-valued functions. A spline is    !
  ! stored as its values y and second derivatives y2 at the knots  !
  ! x. Discontinuous functions are represented by building a       !
  ! separate spline on each layer of the earth model, with the     !
  ! knots at a discontinuity repeated.                             !
  !================================================================!
  use module_kinds
  implicit none

contains


  subroutine natural_spline(x,y,y2)
    !------------------------------------------------------------!
    ! Returns the second derivatives y2 at the knots x of the    !
    ! natural cubic spline through the values y, for which y2    !
    ! vanishes at both ends. The knots must be increasing.       !
    !                                                            !
    ! Continuity of the first derivative at the interior knots   !
    ! gives a symmetric, diagonally dominant tridiagonal system  !
    ! for y2(2:n-1), which is solved by elimination without      !
    ! pivoting.                                                  !
    !------------------------------------------------------------!
    real(dp), dimension(:), intent(in) :: x,y
    real(dp), dimension(:), intent(out) :: y2
    integer(i4b) :: i,n
    real(dp) :: h1,h2,w
    real(dp), dimension(size(x)) :: diag,rhs

    n = size(x)
    y2 = 0.0_dp
    if(n < 3) return

    ! forward elimination: row i couples y2(i-1), y2(i) and y2(i+1)
    ! with coefficients h1, 2*(h1+h2) and h2
    do i = 2,n-1
       h1 = x(i)-x(i-1)
       h2 = x(i+1)-x(i)
       diag(i) = 2.0_dp*(h1+h2)
       rhs(i)  = 6.0_dp*((y(i+1)-y(i))/h2-(y(i)-y(i-1))/h1)
       if(i > 2) then
          w = h1/diag(i-1)
          diag(i) = diag(i)-w*h1
          rhs(i)  = rhs(i)-w*rhs(i-1)
       end if
    end do

    ! back substitution
    y2(n-1) = rhs(n-1)/diag(n-1)
    do i = n-2,2,-1
       y2(i) = (rhs(i)-(x(i+1)-x(i))*y2(i+1))/diag(i)
    end do

    return
  end subroutine natural_spline


  function spline_value(x,y,y2,xx,k)
    !------------------------------------------------------------!
    ! Evaluates at xx the cubic spline with knots x, values y    !
    ! and second derivatives y2, given the index k of the knot   !
    ! at the bottom of the interval containing xx. If this       !
    ! interval has zero length (a discontinuity in the model)    !
    ! the value at knot k is returned.                           !
    !------------------------------------------------------------!
    real(dp), dimension(:), intent(in) :: x,y,y2
    real(dp), intent(in) :: xx
    integer(i4b), intent(in) :: k
    real(dp) :: spline_value
    real(dp) :: h,a,b

    h = x(k+1)-x(k)
    if(h == 0.0_dp) then
       spline_value = y(k)
       return
    end if
    a = (x(k+1)-xx)/h
    b = (xx-x(k))/h
    spline_value = a*y(k)+b*y(k+1) &
         +((a*a-1.0_dp)*a*y2(k)+(b*b-1.0_dp)*b*y2(k+1))*h*h/6.0_dp

    return
  end function spline_value


end module module_spline
