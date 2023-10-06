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
      
      
      use nrtype
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
      sfl3   = sqrt(real(l*(l+1)))
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
      x2 = -x2/sqrt(real(l*(l+1)))
      do i = 3,mp1
         k   = i-1
         f1  = sqrt(real(l+i-1)*(l-i+2))
         f2  = sqrt(real(l+i-2)*(l-i+3))
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


    subroutine xlm_cal(l,m_in,theta,xlm_in)
      ! This routine calculates the functions X_{lm} defined by equation 
      ! B.58 of Dahlen & Tromp. Upward recursion in l is used in the 
      ! calculation the associated Legendre function.
      ! Inputs:
      !        l     = angular degree 
      !        m     = angular order
      !        theta = angular argument in radians
      ! Output:
      !        xlm   = an array of size l+1 containg
      !                X_{0m}, X_{1m}, X_{2m}, ...,X_{lm}
      !
      ! The routine is based upon that in Numerical Recipes. 
      ! It seems to work well at least up to l = 4000.

      use nrtype
      implicit none
      integer(i4b), intent(in) :: l,m_in
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: xlm_in
      integer(i4b) :: ll,i,m
      real(dp) :: x,pmm,pmmp1,pll,fac,fac2
      real(dp), dimension(l+1) :: xlm



      if(m_in >= 0) then
         m = m_in
      else
         m = -m_in
      end if

!      if(size(xlm) < l+1) stop 'xlm_cal: input array too small'
      
      xlm(1:min(l+1,m))=0.0_dp
      if(l >= m) then
         x=cos(theta)   
         pmm=1.0_dp
         if(m > 0) then
            fac=1.0_dp
            fac2=(1.0_dp-x)*(1.0_dp+x)
            fac2=sqrt(fac2)
            do i=1,m
               pmm=pmm*fac*fac2
               fac=fac+2.0_dp
            end do
            if(mod(m,2) == 1) pmm=-pmm
         end if
         xlm(m+1)=pmm
         if(l >= m+1) then
            pmmp1=(real(2*m+1))*x*pmm
            xlm(m+2)=pmmp1
            if(l >= m+2) then
               do ll=m+2,l
                  pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                  pmm=pmmp1
                  pmmp1=pll
                  xlm(ll+1)=pll
               end do
            end if
         end if
      end if


      do i = 1,l+1
         ll = i-1
         if(ll-m+1 > 0) then
            fac = gammln(1.0_dp*(ll-m+1))-gammln(1.0_dp*(ll+m+1))            
            fac = exp(0.5_dp*fac)
            fac = sqrt((2*ll+1)/fourpi_d)*fac
         else
            fac = 0.0_dp
         end if
         xlm(i) = fac*xlm(i)
      end do

      if(m /= m_in) then
         xlm=xlm*(-1)**m_in
      end if
      i = size(xlm_in)
      xlm_in(1:i) = ((-1)**m)*xlm(l+2-i:l+1)


      return
    end subroutine xlm_cal




    subroutine angles(lats,lons,latr,lonr,delta,azep,azst,corr)
      ! this routine computes the epicentral angle and azimuth 
      ! from a given source to a given receiver. the inputs are
      ! the longitude and latitude of the two locations in degress.
      ! note that the answer is returned in radians, and that the 
      ! azimuth is measures anti-clockwise from south.

      ! if the optional input 'corr' is set to .true. that an 
      ! ellipticity correction is made 

      use nrtype
      implicit none
      real(dp), intent(in) :: latr,lonr,lats,lons
      real(dp), intent(out) :: delta,azep,azst
      logical(lgt), intent(in), optional :: corr
      integer(i4b) :: i,nr
      real(dp) :: conv,sgn
      real(dp) :: cltr,clts,clnr,clns,sltr,slnr,slts,slns,saz,sazr
      real(dp), dimension(3) :: rs,rr,ts,ps,ks,kr,tr,pr
      
      ! ellipticity correction factor = 1 - 2 / 298.3
      real(dp), parameter :: flat = 0.9932953402614817_dp


      conv=pi_d/180.0_dp

      ! convert angles to radians
      cltr=(90.0_dp-latr)*conv
      clts=(90.0_dp-lats)*conv
      clnr=lonr*conv
      clns=lons*conv


      ! compute the required sines and cosines
      sltr=sin(cltr); slts=sin(clts); cltr=cos(cltr); clts=cos(clts)
      slnr=sin(clnr); slns=sin(clns); clnr=cos(clnr); clns=cos(clns)

      ! compute the required unit vectors
      rs(1)=slts*clns; rs(2)=slts*slns; rs(3)=clts
      rr(1)=sltr*clnr; rr(2)=sltr*slnr; rr(3)=cltr
      ts(1)=clts*clns; ts(2)=clts*slns; ts(3)=-slts
      tr(1)=cltr*clnr; tr(2)=cltr*slnr; tr(3)=-sltr
      ps(1)=-slns;     ps(2)=clns;      ps(3)=0;
      pr(1)=-slnr;     pr(2)=clnr;      pr(3)=0;


      ! compute cos(delta)
      delta=dot_product(rr,rs)      
      ! compute the tangent vector to the source-recevier
      ! great circle a the source
      if(abs(delta) < 1.0_dp) then
         ks=(rr-delta*rs)/sqrt(1-delta*delta)
      else
         ks=0.0_dp
      end if
      if(abs(delta) < 1.0_dp) then
         kr=(delta*rr-rs)/sqrt(1-delta*delta)
      else
         kr=0.0_dp
      end if
      !  compute the sine and consine of the azimuth
      azep=dot_product(ks,ts)
      saz=dot_product(ks,ps)
      azst=dot_product(kr,tr)
      sazr=dot_product(kr,pr)
      ! check the deltas
      if(delta > 1.0_dp) delta=1.0_dp
      ! check the azimuths
      if(abs(azep) > 1.0_dp) azep=sign(1.0_dp,azep)
      if(abs(azst) > 1.0_dp) azst=sign(1.0_dp,azst)
      ! compute delta
      delta=acos(delta)
      ! compute azimuth picking the right quadrant
      if(saz >= 0.0_dp) then
         azep=acos(azep)
      else
         azep=acos(azep)
         azep=twopi_d-azep
      end if
      if(sazr >= 0.0_dp) then
         azst=acos(azst)
      else
         azst=acos(azst)
         azst=twopi_d-azst
      end if


      return
    end subroutine angles



    
    subroutine delaz(eplati,eplong,stlati,stlong,delta,azep,azst)
      use nrtype
      implicit none
      ! in/out
      real(dp), intent(in) :: eplati,eplong,stlati,stlong
      real(dp), intent(out) :: delta,azep,azst
      
      ! parameters
      real(dp), parameter :: geoco  = 0.993277_dp
      real(dp), parameter :: hpi    = 0.5_dp*pi_d
      real(dp), parameter :: rad    = 0.0174533_dp
      real(dp), parameter :: reprad = 57.29578_dp
      
      ! local variables
      real(dp) :: el,stl,elon,slon,as,bs,cs,ds,a,b,c,d, & 
           cdel,delt,sdel,caze,aze,azs,dif,cazs,eplat,stlat
      
      
      eplat = eplati
      stlat = stlati

      
      if(eplat == -90.0_dp) eplat = -89.98_dp
      if(stlat == -90.0_dp) stlat = -89.98_dp  !station at S pole
      
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





    subroutine sbessj(l,x,jl)
      ! This function returns the value of j_{l}(x) for integer l and
      ! complex x. This routine is based upon Guy Masters routine sphbr,
      ! which uses the continued fraction algorithm described by  W.J. 
      ! Lentz in Applied Optics, pg 668-671, vol. 15, No. 3, 1976.
      use nrtype
      implicit none
      integer(i4b), intent(in) :: l
      complex(dpc), intent(in) :: x
      complex(dpc), dimension(:), intent(out) :: jl
      integer(i4b) :: lp1,n,lout
      real(dp) :: nu
      real(dp), parameter :: eps=1.e-14_dp
      real(dp), parameter :: tiny=1.e-10_dp,huge=1.e+10_dp
      complex(dpc) :: rx,ans,denom,numer,ratio,a,etan,etad, & 
           num,den
      complex(dpc), dimension(l+1) :: j
      lp1=l+1
      lout=min(l+1,size(jl))
      j(1)=sin(x)/x
      if(l > 0) then
         do n=2,lp1
            rx=2.0_dp*one/x
            nu=n-0.5_dp
            ans=nu*rx
            nu=nu+1.0_dp
            rx=-rx
            denom=nu*rx
            numer=denom+one/ans
            do
               ratio=numer/denom
               ans=ans*ratio
               if(abs(abs(ratio)-1.0_dp) > eps) then
                  nu=nu+1.0_dp
                  rx=-rx
                  a=nu*rx
                  denom=a+one/denom
                  numer=a+one/numer
               else
                  exit
               end if
            end do
            j(n)=j(n-1)/ans
         end do
      end if
      jl(1:lout)=j(lp1-lout+1:lp1)

      return
    end subroutine sbessj


    subroutine sbessj2(l,x,ratio,dratio)
      ! this routine computes the ratio of spherical 
      ! Bessel functions of the first kind using a 
      ! continued fraction algorithm 
      !
      ! Inputs: 
      ! l ==  lower order of the Bessel function ratio
      ! x == argument of the Bessel functions
      ! Output:
      ! ratio == j_{l+1}(x)/j_{l}(x)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: l
      complex(dpc), intent(in) :: x
      complex(dpc), intent(out) :: ratio
      complex(dpc), intent(out), optional :: dratio
      integer(i4b) :: i,ll,sign
      real(dp), parameter :: tiny = 1.e-3_dp
      complex(dpc) :: xx,tanx,term

      ratio = x/real(2*l+3)

      if(present(dratio)) then
         dratio = 1.0_dp/real(2*l+3)
      end if

 !     tanx = sin(x)/cos(x)
  !    xx = 1.0_dp/x
  !    ratio = xx-1.0_dp/tanx   
  !    do ll=1,l
  !       ratio = real(2*ll+1)*xx-1.0_dp/ratio
  !    end do
!      xx=1/x
!      ratio=2*(l+0.5)*xx;
!      sign=-1;
!      do i=2,l
!         term=sign*2*(l+i-0.5)*xx
!         ratio=ratio+1/term;
!         sign=-sign;
!      end do
!      ratio=1/ratio

      return
    end subroutine sbessj2






    subroutine bfs_dum(l,xsq,ratio,dratio)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: l
      complex(dpc), intent(in) :: xsq
      complex(dpc), intent(out) :: ratio
      complex(dpc), intent(out), optional :: dratio

      real(dp), parameter :: eps = 1.e-9_dp
      complex(dpc) :: x,fp

      call bfs(l,xsq,eps,fp)
      x = sqrt(abs(xsq))
      ratio = (l-fp)/x
      if(present(dratio)) then
         dratio = 1.0_dp-((l+2)/x+fp)*ratio
      end if
      
      return
    end subroutine bfs_dum



    subroutine bfs(l,xsq,eps,fp)
      !  this routine calculates spherical bessel function of the ist kind.
      !  fp is equivalent to (r*dj/dr)/j
      !  where r is radius and j is the sbf of order l and argument x=k*r
      !  the technique employs the continued fraction approach
      !  described in w. lentz's article in applied optics, vol.15, #3, 1976
      use nrtype
      implicit none
      ! input output
      integer(i4b), intent(in) :: l 
      complex(dpc), intent(in) :: xsq
      real(dp), intent(in) :: eps
      complex(dpc), intent(out) :: fp
      ! local variables
      integer(i4b) :: lp1,nm1
      complex(dpc) :: numer,nu,x,a3,f,a,b,c,d, & 
           rx,ratio,rj,denom


      ! positive argument uses continued fraction
      if(real(xsq) > 0.0_dp) then
        x     = sqrt(xsq)
        lp1   = l+1
        rx    = 2.0_dp/x
        nu    = lp1-0.5_dp
        rj    = nu*rx
        rx    = -rx
        denom =(nu+1.0_dp)*rx
        numer = denom+1.0_dp/rj
        rj    = rj*numer/denom
        nm1   = 1
        do 
           nm1   = nm1+1
           rx    = -rx
           a3    = (nu+nm1)*rx
           denom = a3+1.0_dp/denom
           numer = a3+1.0_dp/numer
           ratio = numer/denom
           rj    = rj*ratio
           if(abs(abs(ratio)-1.0_dp) <= eps) exit
        end do
        fp =rj*x-lp1
     else
        !  series solution
        f  = 1.0_dp
        fp = real(l)
        a  = 1.0_dp
        b  = real(l+l)+1.0_dp
        c  = 2.0_dp
        d  = real(l)+2.0_dp
        do
           a  = -a*xsq/(c*(b+c))
           f  = f+a
           fp = fp+a*d
           if(abs(a*d) < eps) exit
           c = c+2.0_dp
           d = d+2.0_dp
        end do
        fp = fp/f
     end if
     return
   end subroutine bfs



    
    function dfact(n)
      ! This function returns the value of (2n+1)!! for integer n
      use nrtype
      implicit none
      real(dp) :: dfact
      integer(i4b), intent(in) :: n
      dfact=factln(2*n+1)-n*log(2.0_dp)-factln(n)
      return
    end function dfact


    function factln(n)
      ! This function returns the natural logarithm of n! for integer n
      use nrtype
      implicit none
      real(dp) :: factln
      integer(i4b), intent(in) :: n
      real(dp) :: xx
      xx=real(n+1)
      factln=gammln(xx)
      return
    end function factln


    function gammln(xx)
      use nrtype; use nrutil, only : arth,assert
      implicit none
      real(dp), intent(in) :: xx
      real(dp) :: gammln
      ! Returns the natural logarithm of the gamma function for xx>0
      real(dp) :: tmp,x
      real(dp) :: stp=2.5066282746310005_dp
      real(dp), dimension(6) :: coef=(/76.18009172947146_dp, &
           -86.50532032941677_dp,24.01409824083091_dp, &
           -1.231739572450155_dp,0.1208650973866179e-2_dp,&
           -0.5395239384953e-5_dp/)
      call assert(xx > 0.0,'gammln arg')
      x=xx
      tmp=x+5.5_dp
      tmp=(x+0.5_dp)*log(tmp)-tmp
      gammln=tmp+log(stp*(1.000000000190015_dp+&
           sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
    end function gammln





  
end module module_function
