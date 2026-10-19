module module_function


  contains

    subroutine legendre(theta,l,m,x,xp,xcosec)
      
      ! This routine computes the associated legendre 
      ! polynomials P_{lm}(cos(theta)) along with
      ! their angular derivative and  products with
      ! cosec(theta). This routine is a f90 translation
      ! of the routine legndr by John Woodhouse
      
      ! the output arrays should have dimension at least as 
      ! big as m+1
      
      ! inputs: 
      !          theta = angular argument in radians
      !          l     = angular degree of the Legendre functions
      !          m     = maximum angular order to calculate
      !
      ! outputs:
      !          x  = array of the Legendre functions; the array
      !               order is such that x(1) = P_{l0}, x(2) = P_{l1},
      !               etc upto x(m+1) = P_{lm}. 
      !          xp = theta derivatives of the Legendre functions
      !          xcosec = legendre functions times by cosec(theta).
      
      ! the legendre functions returned are such that the fully normalized 
      ! spherical harmonic functions Y_{lm} are given by
      ! 
      ! Y_{lm}(theta,phi) = P_{lm}(cos(theta))*exp(i*m*phi),
      !
      ! The P_{lm}(cos(theta)) calculated are equal to the
      ! X_{lm}(theta) described in appendix B of Dahlen & Tromp (1998).
      
      
      use module_kinds
      implicit none
      
      ! inputs:
      real(dp), intent(in) :: theta
      integer(i4b), intent(in) :: l,m
      ! outputs:
      real(dp), dimension(:), intent(out) :: x,xp,xcosec
      
      ! local variables
      integer(i4b) :: lp1,mp1,i,k
      real(dp) :: sum,th,ct,st,fct,sfl3,compar,dsfl3, &
           cot,cosec,x1,x2,x3,f1,f2,xm,small
      
      sum = 0.0_dp
      lp1 = l+1
      mp1 = m+1
      
      th = theta
      ct = cos(th)
      st = sin(th)
      
      fct    = sqrt(real(2*l+1)/fourpi_d)
      sfl3   = sqrt(real(l*(l+1),dp))
      compar = real(2*l+1)/fourpi_d
      dsfl3  = sfl3
      small  = 1.0e-16_dp*compar
      
      x      = 0.0_dp
      xp     = 0.0_dp
      xcosec = 0.0_dp
      
      if(l <= 1 .or. abs(theta) <= 1.0e-5_dp) then
         x(1) = fct
         if(l == 0) return     
         x(1)  = ct*fct
         x(2)  = -st*fct/dsfl3
         xp(1) = -st*fct
         xp(2) = -0.5_dp*ct*fct*dsfl3
         if(abs(theta) <  1.0e-5_dp) then
            xcosec(2) = xp(2)
         else
            xcosec(2) = x(2)/st
         end if
         return
      end if
      
      x1 = 1.0_dp
      x2 = ct
      
      do i = 2,l
         x3 = (real(2*i-1)*ct*x2-real(i-1)*x1)/real(i)
         x1 = x2
         x2 = x3
      end do
      
      cot   = ct/st;
      cosec = 1.0_dp/st
      
      x3 = x2*fct
      x2 = real(l)*(x1-ct*x2)*fct/st

      x(1) = x3
      x(2) = x2
      sum  = x3*x3
      
      xp(1) = -x2
      xp(2) = real(l*(l+1))*x3-cot*x2
      
      x(2)      = -x(2)/sfl3
      xcosec(2) = x(2)*cosec
      xp(2)     = -xp(2)/sfl3
      
      sum = sum+2.0_dp*x(2)*x(2)
      if(sum-compar > small) return

      x1 =  x3
      x2 = -x2/sqrt(real(l*(l+1),dp))
      do i = 3,mp1
         k   = i-1
         f1  = sqrt(real((l+i-1)*(l-i+2),dp))
         f2  = sqrt(real((l+i-2)*(l-i+3),dp))
         xm  = k
         x3  = -(2.0_dp*cot*(xm-1.0_dp)*x2+f2*x1)/f1
         sum = sum+2.0_dp*x3*x3
         if(sum-compar > small .and. i /= lp1) return
         x(i)      = x3
         xcosec(i) = x(i)*cosec
         x1        = x2
         xp(i)     = -(f1*x2+xm*cot*x3)
         x2        = x3
      end do
      
      return
    end subroutine legendre



    
    subroutine delaz(eplati,eplong,stlati,stlong,delta,azep,azst)
      use module_kinds
      implicit none
      ! in/out
      real(dp), intent(in) :: eplati,eplong,stlati,stlong
      real(dp), intent(out) :: delta,azep,azst
      
      ! parameters
      real(dp), parameter :: geoco  = 0.993277_dp
      real(dp), parameter :: hpi    = 0.5_dp*pi_d
      real(dp), parameter :: rad    = pi_d/180.0_dp

      ! local variables
      real(dp) :: el,stl,elon,slon,as,bs,cs,ds,a,b,c,d, & 
           cdel,delt,sdel,caze,aze,azs,dif,cazs,eplat,stlat
      
      
      eplat = eplati
      stlat = stlati

      
      ! move points off the poles where the azimuths are undefined
      if(eplat ==  90.0_dp) eplat =  89.98_dp
      if(eplat == -90.0_dp) eplat = -89.98_dp
      if(stlat ==  90.0_dp) stlat =  89.98_dp
      if(stlat == -90.0_dp) stlat = -89.98_dp
      
      if(eplat == stlat .and. eplong == stlong) then
         delta = 0.0_dp
         azep  = 0.0_dp
         azst  = 0.0_dp
         return
      end if

      
      el  = atan(geoco*tan(eplat*rad))
      el  = hpi-el
      stl = atan(geoco*tan(stlat*rad))
      stl = hpi-stl

      elon = eplong*rad
      slon = stlong*rad

      as = cos(stl)
      bs = sin(stl)
      cs = cos(slon)
      ds = sin(slon)

      a = cos(el)
      b = sin(el)
      c = cos(elon)
      d = sin(elon)

      cdel = a*as+b*bs*(c*cs+d*ds)
      if(abs(cdel) > 1.) cdel = sign(1.0_dp,cdel)
      delt  = acos(cdel)
      delta = delt

      sdel = sin(delt)
      caze = (as-a*cdel)/(sdel*b)
      if(abs(caze) > 1.0_dp) caze = sign(1.0_dp,caze)
      aze = acos(caze)

      if(bs > 0.0_dp) cazs = (a-as*cdel)/(bs*sdel)
      if(bs == 0.0_dp) cazs = sign(1.0_dp,cazs)
      if(abs(cazs) > 1.0_dp) cazs = sign(1.0_dp,cazs)
      azs = acos(cazs)
      dif = ds*c-cs*d
      if(dif < 0.0_dp) aze = twopi_d-aze
      azep = aze
      if(dif > 0.0_dp) azs = twopi_d-azs
      azst = azs
    

      return
    end subroutine delaz


    function bessel_q(l,z)
      !------------------------------------------------------------!
      ! Returns q = x j_{l+1}(x)/j_{l}(x) for z = x**2, where j_l  !
      ! is the spherical Bessel function of the first kind. Since  !
      ! q is an even function of x it depends only on z, which may !
      ! be complex. The recurrence relation for j_l gives          !
      !                                                            !
      !            q_l(z) = z/(2l+3 - q_{l+1}(z)),                 !
      !                                                            !
      ! and the resulting continued fraction is evaluated using    !
      ! the modified Lentz method (Thompson & Barnett, J. Comput.  !
      ! Phys., 64, 490-509, 1986). The function is used to form    !
      ! the regular solutions in a homogeneous sphere that provide !
      ! the starting values for the radial integrations.           !
      !------------------------------------------------------------!
      use module_kinds
      implicit none
      integer(i4b), intent(in) :: l
      complex(dpc), intent(in) :: z
      complex(dpc) :: bessel_q
      integer(i4b), parameter :: max_terms = 100000
      real(dp), parameter :: eps = 1.0e-14_dp
      real(dp), parameter :: small = 1.0e-300_dp
      integer(i4b) :: k
      complex(dpc) :: f,c,d,delta,b

      ! evaluate the denominator 2l+3 - z/(2l+5 - z/(2l+7 - ...))
      f = 2*l+3
      c = f
      d = 0.0_dp
      do k = 1,max_terms
         b = 2*(l+k)+3
         d = b-z*d
         if(abs(d) < small) d = small
         c = b-z/c
         if(abs(c) < small) c = small
         d = 1.0_dp/d
         delta = c*d
         f = f*delta
         if(abs(delta-1.0_dp) < eps) exit
      end do
      if(k > max_terms) error stop 'bessel_q: continued fraction did not converge'
      bessel_q = z/f

      return
    end function bessel_q


end module module_function
