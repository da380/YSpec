module module_int
  !==========================================================!
  ! This module contains a range of routines used in the     !
  ! numerical integration of the differential equations      !
  ! govering the derformation of a spherically symmetric     !
  ! earth model.                                             !
  !----------------------------------------------------------!
  ! In solid regions the following relationships are used    !
  ! for the 3rd order minor vectors                          !
  ! (1)  m_{7}  = -m_{3}                                     !
  ! (2)  m_{13} =  m_{2}                                     !
  ! (3)  m_{12} = -m_{5}                                     !
  ! (4)  m_{18} = -m_{14}                                    !
  ! (5)  m_{19} =  m_{8}                                     !
  ! (6)  m_{16} = -m_{9}                                     !
  ! While in fliud regions the relation                      !
  ! (1)  m_{5}  = -m_{2}                                     !
  ! is used.                                                 !
  ! The orderings used for the reduced system of minor       !
  ! vectors is then:                                         !
  ! In solid regions,                                        !
  ! m(1)  => m_{1}                                           !
  ! m(2)  => m_{2}                                           ! 
  ! m(3)  => m_{3}                                           ! 
  ! m(4)  => m_{4}                                           ! 
  ! m(5)  => m_{5}                                           ! 
  ! m(6)  => m_{6}                                           ! 
  ! m(7)  => m_{8}                                           ! 
  ! m(8)  => m_{9}                                           ! 
  ! m(9)  => m_{10}                                          ! 
  ! m(10) => m_{11}                                          ! 
  ! m(11) => m_{14}                                          ! 
  ! m(12) => m_{15}                                          ! 
  ! m(13) => m_{17}                                          ! 
  ! m(14) => m_{20}                                          ! 
  ! In fluid regions,                                        !
  ! m(1)  => m_{1}                                           !
  ! m(2)  => m_{2}                                           !
  ! m(3)  => m_{3}                                           !
  ! m(4)  => m_{4}                                           !
  ! m(5)  => m_{6}                                           !
  !                                                          !
  ! For the (n-1)th order minor vectors, written as b, we    !
  ! make use of the relation                                 !
  ! (1) b_{3}+b_{8}+b_{12}=0                                 !
  ! to reduce the 15th order system of equations to a 14th   !
  ! order system. The numbering of the resulting system is:  !
  ! b(1)  => b_{1}                                           !
  ! b(2)  => b_{2}                                           !
  ! b(3)  => b_{3}                                           !
  ! b(4)  => b_{4}                                           !
  ! b(5)  => b_{5}                                           !
  ! b(6)  => b_{6}                                           !
  ! b(7)  => b_{7}                                           !
  ! b(8)  => b_{8}                                           !
  ! b(9)  => b_{9}                                           !
  ! b(10) => b_{10}                                          !
  ! b(11) => b_{11}                                          !
  ! b(12) => b_{13}                                          !
  ! b(13) => b_{14}                                          !
  ! b(14) => b_{15}                                          !
  !==========================================================!



  use nrtype
  use module_spline
  implicit none
  
  !===========================================================!
  !            parameters used in the module                  !
  !===========================================================!

  integer(i4b), parameter, private :: nw = 4        ! this parameter is used 
                                                    ! in determining the 
                                                    ! step  size in the numerical 
                                                    ! integrations

                                            
  real(dp), parameter, public :: tiny=1.0e-20_dp    ! these parameters are used in the 
  real(dp), parameter, public :: huge=1/tiny        ! scaling of the numerical integrations.




                                                    


  !================================================!
  !       saved variables in the module            !
  !================================================!

  integer(i4b), save, public :: il       ! current level of the model
  integer(i4b), save, public :: ll       ! equal to the spherical harmonic degree l
  integer(i4b), save, public :: llp1     ! equal to l+1

  real(dp), save, public :: zeta         ! equal to sqrt(l(l+1))
  real(dp), save, public :: zeta2        ! equal to l(l+1)
  real(dp), save, public :: zeta2m2      ! equal to l(l+1)-2
  real(dp), save, public :: szeta2m2     ! equal to sqrt(l(l+1)-2)
  real(dp), save, public :: nu           ! equal to sqrt((2l+1)/4*pi)
  real(dp), save, public :: dr           ! initial step size for the numerical integrations
  real(dp), save, public :: ats          ! parameter used to scale 
                                         ! attenuation in the model
                                         ! ats=0 gives no attenuation, 
                                         ! ats=1 gives full attenuation.



  complex(dpc), save, public :: omega    ! current value of angular frequency w
  complex(dpc), save, public :: omega2   ! equal to w*w
  complex(dpc), save, public :: omegai   ! equal to 1/w
  complex(dpc), save, public :: omegai2  ! equal to 1/(w*w)
  complex(dpc), save, public :: lno      ! equal to  (2/pi)*ln(i*w/w0), with w0 the 
                                         ! reference frequency of the earth model.
                                         ! This value is used in the computation of
                                         ! the anelastic moduli.



  contains




    !========================================================!
    !========================================================!
    !                Toroidal mode section                   !
    !========================================================!
    !========================================================!



    subroutine tor_amat(rr,t,k,s)
      !------------------------------------------------!
      ! This routine returns the elements of the       !
      ! toroidal mode coefficient matrix at a given    !
      ! radius in the earth model.                     !
      !------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      real(dp), intent(out) :: t
      complex(dpc), intent(out) :: k,s
      real(dp) :: rrho,rlcon,rncon,rmu,rqmu,ri
      complex(dpc) :: clcon,cncon,qll
      ! evaluate cubic splines
      rrho  = splint_dis(r,rho,rho_cs,rr,il)
      rlcon = splint_dis(r,lcon,lcon_cs,rr,il)
      rncon = splint_dis(r,ncon,ncon_cs,rr,il)
      rmu   = splint_dis(r,mu,mu_cs,rr,il)
      rqmu  = splint_dis(r,qmu,qmu_cs,rr,il)      
      ! compute parameters
      qll = 1.0_dp+rqmu*lno
      clcon = rlcon*qll
      cncon = rncon*qll
      ri    = 1.0_dp/rr
      ! compute matrix elements
      t = 2.0_dp*ri
      k = 1.0_dp/clcon
      s = -omega2*rrho+ri*ri*zeta2m2*cncon
      return
    end subroutine tor_amat



    subroutine tor_derivs(rr,y,dydx)
      !---------------------------------------------!
      ! This routine evaluates the derivatives for  !
      ! the toroidal mode system.                   !
      !---------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: y
      complex(dpc), dimension(:), intent(inout) :: dydx
      real(dp) :: t
      complex(dpc) :: k,s
      call tor_amat(rr,t,k,s)
      dydx(1) = t*y(1)+k*y(2)
      dydx(2) = s*y(1)-t*y(2)
      return
    end subroutine tor_derivs



    subroutine tor_steps
      !--------------------------------------------!
      ! This routine returns an initial step size  !
      ! for the integration of the toroidal mode   !
      ! system. This guess is not very important,  !
      ! as the step size used will be modified by  !
      ! the numerical integration routine.         !
      !--------------------------------------------!
      use nrtype; use module_model
      use nrutil
      implicit none
      integer(i4b) :: im
      im=iminloc(beta(noc+1:nsl))      
      im=im+noc
      if(real(omega) == 0.0_dp) then
         dr=1.0_dp/nw
      else
         dr=nw*real(omega)/beta(im)
         dr=1.0_dp/dr
      end if
      return
    end subroutine tor_steps
      


    subroutine tor_ystart(is,rr,y)
      !-------------------------------------------!
      ! This routine returns the starting values  !
      ! for the toroidal mode equations in terms  !
      ! of spherical Bessel functions             !
      !-------------------------------------------!
      use nrtype; use module_function
      use module_spline; use module_model
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(2), intent(out) :: y
      real(dp) :: rrho,rmu,rqmu
      complex(dpc) :: cmu,cbeta,arg,ym,ratio,kbeta,sarg,carg,qll

      ! evaluate cubic splines
      rrho = splint_dis(r,rho,rho_cs,rr,is)
      rmu  = splint_dis(r,mu,mu_cs,rr,is)
      rqmu = splint_dis(r,qmu,qmu_cs,rr,is)      

      ! compute parameters
      qll = 1.0_dp+rqmu*lno
      cmu   = rmu*qll
      cbeta = sqrt(cmu/rrho) 
      kbeta = omega/cbeta

      ! compute spherical Bessel functions
      arg = kbeta*rr
      call sbessj2(llp1,arg,ratio)

      ! compute solutions
      y(1) = rr
      y(2) = cmu*((ll-1)-arg*ratio) 
      y    = zeta*y


      return
    end subroutine tor_ystart


    !=============================================================!
    !=============================================================!
    !                    Radial mode section                      !
    !=============================================================!
    !=============================================================!
    


    subroutine rad_amat(rr,t,k,s)
      !------------------------------------------------!
      ! This routine returns the elements of the       !
      ! radial mode coefficient matrix at a given      !
      ! radius in the earth model.                     !
      !------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t,k,s
      real(dp) :: rrho,racon,rccon,rfcon,rmu, & 
           rqmu,ri,gg,rkappa,rqkappa,rncon,rxlam,rxa2
      complex(dpc) :: cacon,cccon,cfcon,gamma,cncon,qaa,qll,qff
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,il)
      racon   = splint_dis(r,acon,acon_cs,rr,il)
      rccon   = splint_dis(r,ccon,ccon_cs,rr,il)
      rfcon   = splint_dis(r,fcon,fcon_cs,rr,il)
      rncon   = splint_dis(r,ncon,ncon_cs,rr,il)
      gg      = splint_dis(r,grav,grav_cs,rr,il)
      rmu     = splint_dis(r,mu,mu_cs,rr,il)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,il)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxlam    = splint_dis(r,xlam,xlam_cs,rr,il)
      rxa2    = splint_dis(r,xa2,xa2_cs,rr,il)
      ! compute parameters
      qaa = 1.0_dp+rxa2*lno
      qff = 1.0_dp+rxlam*lno
      qll = 1.0_dp+rqmu*lno
      cacon = racon*qaa
      cccon = rccon*qaa
      cfcon = rfcon*qff
      cncon = rncon*qll
      cccon = 1.0_dp/cccon
      gamma = cacon-cncon-cfcon*cfcon*cccon      
      ri    = 1.0_dp/rr
      ! compute matrix elements
      t = (1.0_dp-2.0_dp*cfcon*cccon)*ri
      k = cccon
      s = -omega2*rrho+4.0_dp*(gamma-rrho*gg*rr)*ri*ri
      return
    end subroutine rad_amat



    subroutine rad_amat_ng(rr,t,k,s)
      !------------------------------------------------!
      ! This routine returns the elements of the       !
      ! radial mode coefficient matrix at a given      !
      ! radius in the earth model.                     !
      !------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t,k,s
      real(dp) :: rrho,racon,rccon,rfcon,rmu, & 
           rqmu,ri,rkappa,rqkappa,rncon,rxlam,rxa2
      complex(dpc) :: cacon,cccon,cfcon,gamma,cncon,qaa,qll,qff
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,il)
      racon   = splint_dis(r,acon,acon_cs,rr,il)
      rccon   = splint_dis(r,ccon,ccon_cs,rr,il)
      rfcon   = splint_dis(r,fcon,fcon_cs,rr,il)
      rncon   = splint_dis(r,ncon,ncon_cs,rr,il)
      rmu     = splint_dis(r,mu,mu_cs,rr,il)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,il)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxlam    = splint_dis(r,xlam,xlam_cs,rr,il)
      rxa2    = splint_dis(r,xa2,xa2_cs,rr,il)
      ! compute parameters
      qaa = 1.0_dp+rxa2*lno
      qff = 1.0_dp+rxlam*lno
      qll = 1.0_dp+rqmu*lno
      cacon = racon*qaa
      cccon = rccon*qaa
      cfcon = rfcon*qff
      cncon = rncon*qll
      cccon = 1.0_dp/cccon
      gamma = cacon-cncon-cfcon*cfcon*cccon      
      ri    = 1.0_dp/rr
      ! compute matrix elements
      t = (1.0_dp-2.0_dp*cfcon*cccon)*ri
      k = cccon
      s = -omega2*rrho+4.0_dp*gamma*ri*ri
      return
    end subroutine rad_amat_ng
    





    subroutine rad_derivs(rr,y,dydx)
      !---------------------------------------------!
      ! This routine evaluates the derivatives for  !
      ! the radial mode system.                     !
      !---------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: y
      complex(dpc), dimension(:), intent(inout) :: dydx
      complex(dpc) :: t,k,s
      call rad_amat(rr,t,k,s)
      dydx(1) = t*y(1)+k*y(2)
      dydx(2) = s*y(1)-t*y(2)
      return
    end subroutine rad_derivs




    subroutine rad_derivs_ng(rr,y,dydx)
      !---------------------------------------------!
      ! This routine evaluates the derivatives for  !
      ! the radial mode system.                     !
      !---------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: y
      complex(dpc), dimension(:), intent(inout) :: dydx
      complex(dpc) :: t,k,s
      call rad_amat(rr,t,k,s)
      dydx(1) = t*y(1)+k*y(2)
      dydx(2) = s*y(1)-t*y(2)
      return
    end subroutine rad_derivs_ng

    subroutine rad_steps
      !--------------------------------------------!
      ! This routine returns an initial step size  !
      ! for the integration of the radial mode     !
      ! system. This guess is not very important,  !
      ! as the step size used will be modified by  !
      ! the numerical integration routine.         !
      !--------------------------------------------!
      use nrtype; use module_model
      use nrutil
      implicit none
      integer(i4b) :: im
      im=iminloc(alpha)      
      if(real(omega) == 0.0_dp) then
         dr = 1.0_dp/nw
      else
         dr = nw*real(omega)/alpha(im)
         dr = 1.0_dp/dr
      end if
      return
    end subroutine rad_steps


    subroutine rad_ystart(is,rr,y)
      !--------------------------------------------!
      ! this routine returns the starting values   !
      ! for the radial mode equations in terms of  !
      ! spherical Bessel functions                 !
      !--------------------------------------------!
      use nrtype; use module_function
      use module_spline; use module_model
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(2), intent(out) :: y
      real(dp) :: rrho,rmu,rqmu,ri,gg,rkappa,rqkappa
      complex(dpc) :: gamma,cmu,ckappa,calpha,ratio,arg
      complex(dpc), dimension(2) :: jl
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,is)
      gg      = splint_dis(r,grav,grav_cs,rr,is)
      rmu     = splint_dis(r,mu,mu_cs,rr,is)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,is)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,is)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,is)      
      ! compute parameters
      cmu    = rmu*(1.0_dp+rqmu*lno)
      ckappa = rkappa*(1.0_dp+rqkappa*lno)
      calpha = (ckappa+4.0_dp*cmu/3.0_dp)
      gamma  = omega2+16.0_dp*pibigg*rrho/3.0_dp
      gamma  = rrho*gamma/calpha
      ! compute spherical Bessel function ratios
      arg=gamma*rr      
      call sbessj2(1,arg,ratio)      
      ! compute solution
      y(1)=rr*ratio
      y(2)=rr*calpha*gamma-4.0_dp*cmu*ratio
      return      
    end subroutine rad_ystart


    subroutine rad_ystart_ng(is,rr,y)
      !--------------------------------------------!
      ! this routine returns the starting values   !
      ! for the radial mode equations in terms of  !
      ! spherical Bessel functions                 !
      !--------------------------------------------!
      use nrtype; use module_function
      use module_spline; use module_model
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(2), intent(out) :: y
      real(dp) :: rrho,rmu,rqmu,ri,rkappa,rqkappa
      complex(dpc) :: gamma,cmu,ckappa,calpha,ratio,arg
      complex(dpc), dimension(2) :: jl
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,is)
      rmu     = splint_dis(r,mu,mu_cs,rr,is)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,is)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,is)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,is)      
      ! compute parameters
      cmu    = rmu*(1.0_dp+rqmu*lno)
      ckappa = rkappa*(1.0_dp+rqkappa*lno)
      calpha = (ckappa+4.0_dp*cmu/3.0_dp)
      gamma  = omega2
      gamma  = rrho*gamma/calpha
      ! compute spherical Bessel function ratios
      arg=gamma*rr      
      call sbessj2(1,arg,ratio)      
      ! compute solution
      y(1)=rr*ratio
      y(2)=rr*calpha*gamma-4.0_dp*cmu*ratio
      return      
    end subroutine rad_ystart_ng


 
    !=============================================================!
    !=============================================================!
    !                   Spheroidal mode section                   !
    !=============================================================!
    !=============================================================!


    subroutine sph_solid_amat(rr,t11,t12,t21,t22,t31,t33, & 
         k11,k22,k33,s11,s12,s13,s22,s23)
      !-------------------------------------------------------!
      ! This routine returns the components of the spheroidal !
      ! mode coefficient matrix in solid parts of the earth   !
      ! model.                                                !
      !-------------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t11,t12,t21,t22,t31,t33, & 
         k11,k22,k33,s11,s12,s13,s22,s23
      real(dp) :: rrho,racon,rccon,rfcon,rmu, & 
           rqmu,ri,gg,rkappa,rqkappa,rncon,rlcon,rxlam,rxa2
      complex(dpc) :: cacon,cccon,cfcon,gamma,cncon,clcon,qll,qaa,qff
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,il)
      racon   = splint_dis(r,acon,acon_cs,rr,il)
      rccon   = splint_dis(r,ccon,ccon_cs,rr,il)
      rfcon   = splint_dis(r,fcon,fcon_cs,rr,il)
      rncon   = splint_dis(r,ncon,ncon_cs,rr,il)
      rlcon   = splint_dis(r,lcon,lcon_cs,rr,il)
      gg      = splint_dis(r,grav,grav_cs,rr,il)
      rmu     = splint_dis(r,mu,mu_cs,rr,il)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,il)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxlam = splint_dis(r,xlam,xlam_cs,rr,il)      
      rxa2 = splint_dis(r,xa2,xa2_cs,rr,il)      
      ! compute parameters
      qaa = 1.0_dp+rxa2*lno
      qff = 1.0_dp+rxlam*lno
      qll = 1.0_dp+rqmu*lno
      cacon = racon*qaa
      cccon = rccon*qaa
      cfcon = rfcon*qff
      cncon = rncon*qll
      clcon = rlcon*qll
      cccon = 1.0_dp/cccon
      gamma = cacon-cncon-cfcon*cfcon*cccon
      ri=1.0_dp/rr      
      ! compute matrix elements
      t11 = (1.0_dp-2.0_dp*cfcon*cccon)*ri
      t12 = zeta*ri*cfcon*cccon
      t21 = -zeta*ri
      t22 = 2.0_dp*ri
      t31 = -4*pibigg*rrho
      t33 = -ll*ri
      k11 = cccon
      k22 = 1.0_dp/clcon
      k33 = 4.0_dp*pibigg
      s11 = -omega2*rrho+4.0_dp*(gamma-rrho*gg*rr)*ri*ri      
      s12 = zeta*(rrho*gg*rr-2.0_dp*gamma)*ri*ri
      s13 = -rrho*llp1*ri
      s22 = -omega2*rrho+(zeta2*(gamma+cncon)-2.0_dp*cncon)*ri*ri
      s23 = rrho*zeta*ri


      return
    end subroutine sph_solid_amat


    
    subroutine sph_solid_amat_ng(rr,t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22)
      !-------------------------------------------------------!
      ! This routine returns the components of the spheroidal !
      ! mode coefficient matrix in solid parts of the earth   !
      ! model. Neglects self-gravitation                      !
      !-------------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22
      real(dp) :: rrho,racon,rccon,rfcon,rmu, & 
           rqmu,ri,rkappa,rqkappa,rncon,rlcon,rxlam,rxa2
      complex(dpc) :: cacon,cccon,cfcon,gamma,cncon,clcon,qff,qll,qaa
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,il)
      racon   = splint_dis(r,acon,acon_cs,rr,il)
      rccon   = splint_dis(r,ccon,ccon_cs,rr,il)
      rfcon   = splint_dis(r,fcon,fcon_cs,rr,il)
      rncon   = splint_dis(r,ncon,ncon_cs,rr,il)
      rlcon   = splint_dis(r,lcon,lcon_cs,rr,il)
      rmu     = splint_dis(r,mu,mu_cs,rr,il)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,il)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxlam = splint_dis(r,xlam,xlam_cs,rr,il)      
      rxa2 = splint_dis(r,xa2,xa2_cs,rr,il)      
      ! compute parameters
      qaa = 1.0_dp+rxa2*lno
      qff = 1.0_dp+rxlam*lno
      qll = 1.0_dp+rqmu*lno
      cacon = racon*qaa
      cccon = rccon*qaa
      cfcon = rfcon*qff
      cncon = rncon*qll
      clcon = rlcon*qll
      cccon = 1.0_dp/cccon
      gamma = cacon-cncon-cfcon*cfcon*cccon
      ri    = 1.0_dp/rr      
      ! compute matrix elements
      t11 = (1.0_dp-2.0_dp*cfcon*cccon)*ri
      t12 = zeta*ri*cfcon*cccon
      t21 = -zeta*ri
      t22 = 2.0_dp*ri
      k11 = cccon
      k22 = 1.0_dp/clcon
      s11 = -omega2*rrho+4.0_dp*gamma*ri*ri      
      s12 = -2.0_dp*zeta*gamma*ri*ri
      s22 = -omega2*rrho+(zeta2*(gamma+cncon)-2.0_dp*cncon)*ri*ri
      return
    end subroutine sph_solid_amat_ng


    subroutine sph_solid_amat_cow(rr,t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22)
      !-------------------------------------------------------!
      ! This routine returns the components of the spheroidal !
      ! mode coefficient matrix in solid parts of the earth   !
      ! model.  Uses Cowling approximation                    !
      !-------------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22
      real(dp) :: rrho,racon,rccon,rfcon,rmu, & 
           rqmu,ri,gg,rkappa,rqkappa,rncon,rlcon,rxlam,rxa2
      complex(dpc) :: cacon,cccon,cfcon,gamma,cncon,clcon,qaa,qff,qll
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,il)
      racon   = splint_dis(r,acon,acon_cs,rr,il)
      rccon   = splint_dis(r,ccon,ccon_cs,rr,il)
      rfcon   = splint_dis(r,fcon,fcon_cs,rr,il)
      rncon   = splint_dis(r,ncon,ncon_cs,rr,il)
      rlcon   = splint_dis(r,lcon,lcon_cs,rr,il)
      gg      = splint_dis(r,grav,grav_cs,rr,il)
      rmu     = splint_dis(r,mu,mu_cs,rr,il)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,il)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxlam = splint_dis(r,xlam,xlam_cs,rr,il)      
      rxa2 = splint_dis(r,xa2,xa2_cs,rr,il)      
      ! compute parameters
      qaa = 1.0_dp+rxa2*lno
      qff = 1.0_dp+rxlam*lno
      qll = 1.0_dp+rqmu*lno
      cacon = racon*qaa
      cccon = rccon*qaa
      cfcon = rfcon*qff
      cncon = rncon*qll
      clcon = rlcon*qll
      cccon = 1.0_dp/cccon
      gamma = cacon-cncon-cfcon*cfcon*cccon
      ri    = 1.0_dp/rr      
      ! compute matrix elements
      t11 = (1.0_dp-2.0_dp*cfcon*cccon)*ri
      t12 = zeta*ri*cfcon*cccon
      t21 = -zeta*ri
      t22 = 2.0_dp*ri
      k11 = cccon
      k22 = 1.0_dp/clcon
      s11 = -omega2*rrho+4.0_dp*(gamma-rrho*gg*rr)*ri*ri      
      s12 = zeta*(rrho*gg*rr-2.0_dp*gamma)*ri*ri
      s22 = -omega2*rrho+(zeta2*(gamma+cncon)-2.0_dp*cncon)*ri*ri
      return
    end subroutine sph_solid_amat_cow


    
    subroutine sph_fluid_amat(rr,t11,t12,t21,t22,k11,k22, & 
        s11,s12,s22)
      !-------------------------------------------------------!
      ! This routine returns the components of the spheroidal !
      ! mode coefficient matrix in fluid parts of the earth   !
      ! model.                                                !
      !-------------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22
      real(dp) :: rrho,rccon, & 
           ri,gg,rkappa,rqkappa,rxa2
      complex(dpc) :: cccon,qaa
      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,il)
      rccon   = splint_dis(r,ccon,ccon_cs,rr,il)
      gg      = splint_dis(r,grav,grav_cs,rr,il)     
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxa2 = splint_dis(r,xa2,xa2_cs,rr,il)      
      ! compute parameters
      qaa = 1.0_dp+rxa2*lno
      cccon = rccon*qaa
      cccon = 1.0_dp/cccon
      ri    = 1.0_dp/rr      
      ! compute matrix elements
      t11 = (gg*zeta2*omegai2*ri-1.0_dp)*ri
      t12 = zeta2*omegai2*ri*ri
      t21 = -4.0_dp*pibigg*rrho
      t22 = -ll*ri
      k11 = cccon-zeta2*omegai2*ri*ri/rrho
      k22 = 4.0_dp*pibigg
      s11 = -omega2*rrho+rrho*gg*(gg*zeta2*ri*omegai2-4.0_dp)*ri
      s12 = rrho*(gg*zeta2*ri*omegai2-llp1)*ri
      s22 = rrho*zeta2*ri*ri*omegai2
      return
    end subroutine sph_fluid_amat




    subroutine sph_fluid_amat_ng(rr,t,k,s)
      !-------------------------------------------------------!
      ! This routine returns the components of the spheroidal !
      ! mode coefficient matrix in fluid parts of the earth   !
      ! model.                                                !
      !-------------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t,k,s
      real(dp) :: rrho,rccon, & 
           ri,gg,rkappa,rqkappa,rxa2
      complex(dpc) :: cccon,qaa
      rrho=splint_dis(r,rho,rho_cs,rr,il)
      rccon=splint_dis(r,ccon,ccon_cs,rr,il)
      rkappa=splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa=splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxa2=splint_dis(r,xa2,xa2_cs,rr,il)      
      qaa   = 1.0_dp+rxa2*lno
      cccon = rccon*qaa
      cccon = 1.0_dp/cccon
      ri    = 1.0_dp/rr      
      t = -ri
      k = cccon-zeta2*omegai2*ri*ri/rrho
      s = -omega2*rrho     
      return
    end subroutine sph_fluid_amat_ng



    subroutine sph_fluid_amat_cow(rr,t,k,s)
      !-------------------------------------------------------!
      ! This routine returns the components of the spheroidal !
      ! mode coefficient matrix in fluid parts of the earth   !
      ! model.                                                !
      !-------------------------------------------------------!
      use nrtype; use module_model
      use module_spline
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), intent(out) :: t,k,s
      real(dp) :: rrho,racon,rccon, & 
           ri,gg,rkappa,rqkappa,rxa2
      complex(dpc) :: cacon,cccon,qaa

      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,il)
      rccon   = splint_dis(r,ccon,ccon_cs,rr,il)
      gg      = splint_dis(r,grav,grav_cs,rr,il)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,il)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,il)      
      rxa2=splint_dis(r,xa2,xa2_cs,rr,il)      

      ! compute some parameters
      qaa   = 1.0_dp+rxa2*lno
      cccon = rccon*qaa
      cccon = 1.0_dp/cccon
      ri=1.0_dp/rr      

      ! compute matrix elements
      t = (gg*zeta2*omegai2*ri-1.0_dp)*ri
      k = cccon-zeta2*omegai2*ri*ri/rrho
      s = -omega2*rrho+rrho*gg*(gg*zeta2*ri*omegai2-4.0_dp)*ri

      return
    end subroutine sph_fluid_amat_cow
    

    subroutine sph_solid_derivs(rr,y,dydr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 3rd   !
      ! order minors vector equations for the spheroidal  !
      ! mode system in solid parts of the earth model.    !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: y
      complex(dpc), dimension(:), intent(inout) :: dydr
      complex(dpc) :: t11,t12,t21,t22,t31,t33, & 
         k11,k22,k33,s11,s12,s13,s22,s23

      call sph_solid_amat(rr,t11,t12,t21,t22,t31,t33, &            
         k11,k22,k33,s11,s12,s13,s22,s23)      

      dydr(1) = t11*y(1) + t12*y(2) + k11*y(4)
      dydr(2) = t21*y(1) + t22*y(2) + k22*y(5)
      dydr(3) = t31*y(1) + t33*y(3) + k33*y(6)
      dydr(4) = s11*y(1) + s12*y(2) + s13*y(3) - t11*y(4) - t21*y(5) - t31*y(6)
      dydr(5) = s12*y(1) + s22*y(2) + s23*y(3) - t12*y(4) - t22*y(5)
      dydr(6) = s13*y(1) + s23*y(2) - t33*y(6)

      return
    end subroutine sph_solid_derivs



    subroutine sph_minors_solid_derivs(rr,m,dmdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 3rd   !
      ! order minors vector equations for the spheroidal  !
      ! mode system in solid parts of the earth model.    !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: m
      complex(dpc), dimension(:), intent(inout) :: dmdr
      complex(dpc) :: t11,t12,t21,t22,t31,t33, & 
         k11,k22,k33,s11,s12,s13,s22,s23,a123,a126,a135,a234

      call sph_solid_amat(rr,t11,t12,t21,t22,t31,t33, &            
         k11,k22,k33,s11,s12,s13,s22,s23)      

      ! compute some parameters
      a123 =  t11+t22+t33
      a126 =  t11+t22-t33
      a135 =  t11-t22+t33
      a234 = -t11+t22+t33

      ! compute derivatives
      dmdr(1)  = k33*m(4)-k22*m(6)+k11*m(10)+m(1)*a123
      dmdr(2)  = -k22*m(7)+m(1)*s13-m(3)*t21 & 
                 +m(2)*t22-m(4)*t31
      dmdr(3)  = -k11*m(11)+m(1)*s23+m(3)*t11-m(2)*t12
      dmdr(4)  = k22*m(9)-k11*m(12)+m(4)*a126
      dmdr(5)  = -k33*m(8)-m(1)*s12+m(10)*t12 & 
                 -m(6)*t21+m(3)*t31+m(5)*t33
      dmdr(6)  = -k33*m(9)-k11*m(13)-m(1)*s22 & 
                 -2.0_dp*m(5)*t12+m(6)*a135
      dmdr(7)  = m(3)*s12+m(6)*s13-m(2)*s22 & 
                 -m(5)*s23+m(11)*t12-m(7)*t22+m(9)*t31
      dmdr(8)  = m(4)*s12-m(3)*s13-m(2)*s23 & 
                 +m(12)*t12-m(9)*t21-m(8)*t33
      dmdr(9)  = k11*m(14)+m(4)*s22-2.0_dp*m(3)*s23 & 
                 -2.0_dp*m(8)*t12-m(9)*a234
      dmdr(10) = -k33*m(12)+k22*m(13)+m(1)*s11 & 
                 +2.0_dp*m(5)*t21-2.0_dp*m(2)*t31+m(10)*a234
      dmdr(11) = -m(3)*s11+m(2)*s12-m(5)*s13 & 
                 -m(10)*s23-m(11)*t11+m(7)*t21-m(8)*t31
      dmdr(12) = -k22*m(14)-m(4)*s11+2.0_dp*m(2)*s13 & 
                 +2.0_dp*m(8)*t21-m(12)*a135
      dmdr(13) = k33*m(14)-m(6)*s11+2.0_dp*m(5)*s12 & 
                 +m(10)*s22+2.0_dp*m(7)*t31-m(13)*a126
      dmdr(14) = m(9)*s11-2.0_dp*m(8)*s12+2.0_dp*m(7)*s13 & 
                -m(12)*s22+2.0_dp*m(11)*s23-m(14)*a123



      return
    end subroutine sph_minors_solid_derivs


    subroutine sph_minors_solid_derivs_ng(rr,m,dmdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 3rd   !
      ! order minors vector equations for the spheroidal  !
      ! mode system in solid parts of the earth model.    !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: m
      complex(dpc), dimension(:), intent(inout) :: dmdr
      complex(dpc) :: t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22,a11,a22

      call sph_solid_amat_ng(rr,t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22)

      ! compute some parameters
      a11 = t11+t22
      a22 = t11-t22

      ! compute derivatives
      dmdr(1) = a11*m(1)+k22*m(3)-k11*m(4)
      dmdr(2) = s12*m(1)-t21*m(3)+t12*m(4)
      dmdr(3) = s22*m(1)-2.0_dp*t12*m(2)+a22*m(3)+k11*m(5)
      dmdr(4) = -s11*m(1)+2.0_dp*t21*m(2)-a22*m(4)-k22*m(5)
      dmdr(5) = -2.0_dp*s12*m(2)+s11*m(3)-s22*m(4)-a11*m(5)     

      return
    end subroutine sph_minors_solid_derivs_ng



    subroutine sph_minors_solid_derivs_cow(rr,m,dmdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 3rd   !
      ! order minors vector equations for the spheroidal  !
      ! mode system in solid parts of the earth model.    !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: m
      complex(dpc), dimension(:), intent(inout) :: dmdr
      complex(dpc) :: t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22,a11,a22

      call sph_solid_amat_cow(rr,t11,t12,t21,t22, & 
         k11,k22,s11,s12,s22)

      ! compute some parameters
      a11 = t11+t22
      a22 = t11-t22

      ! compute derivatives
      dmdr(1)=a11*m(1)+k22*m(3)-k11*m(4)
      dmdr(2)=s12*m(1)-t21*m(3)+t12*m(4)
      dmdr(3)=s22*m(1)-2.0_dp*t12*m(2)+a22*m(3)+k11*m(5)
      dmdr(4)=-s11*m(1)+2.0_dp*t21*m(2)-a22*m(4)-k22*m(5)
      dmdr(5)=-2.0_dp*s12*m(2)+s11*m(3)-s22*m(4)-a11*m(5)

      return
    end subroutine sph_minors_solid_derivs_cow




    subroutine sph_minors_fluid_derivs(rr,m,dmdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 2nd   !
      ! order minors vector equations for the spheroidal  !
      ! mode system in fluid parts of the earth model.    !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: m
      complex(dpc), dimension(:), intent(inout) :: dmdr
      complex(dpc) :: t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22
     

      call sph_fluid_amat(rr,t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22)

      dmdr(1) = k22*m(3)-k11*m(4)+m(1)*(t11+t22)
      dmdr(2) = m(1)*s12+m(4)*t12-m(3)*t21
      dmdr(3) = k11*m(5)+m(1)*s22-2.0_dp*m(2)*t12+m(3)*(t11-t22)
      dmdr(4) = -k22*m(5)-m(1)*s11+2.0_dp*m(2)*t21+m(4)*(-t11+t22)
      dmdr(5) = m(3)*s11-2.0_dp*m(2)*s12-m(4)*s22-m(5)*(t11+t22)

  
      return
    end subroutine sph_minors_fluid_derivs



    subroutine sph_minors_fluid_derivs_ng(rr,m,dmdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 2nd   !
      ! order minors vector equations for the spheroidal  !
      ! mode system in fluid parts of the earth model.    !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: m
      complex(dpc), dimension(:), intent(inout) :: dmdr
      complex(dpc) :: t,k,s

      call sph_fluid_amat_ng(rr,t,k,s)

      dmdr(1) = t*m(1)+k*m(2)
      dmdr(2) = s*m(1)-t*m(2)

      return
    end subroutine sph_minors_fluid_derivs_ng


    subroutine sph_minors_fluid_derivs_cow(rr,m,dmdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 2nd   !
      ! order minors vector equations for the spheroidal  !
      ! mode system in fluid parts of the earth model.    !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: m
      complex(dpc), dimension(:), intent(inout) :: dmdr
      complex(dpc) :: t,k,s

      call sph_fluid_amat_cow(rr,t,k,s)

      dmdr(1) = t*m(1)+k*m(2)
      dmdr(2) = s*m(1)-t*m(2)

      return
    end subroutine sph_minors_fluid_derivs_cow



    subroutine sph_bvec_solid_derivs(rr,b,dbdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 2nd   !
      ! order 'b-vector' equations in solid parts of the  !
      ! earth model.                                      !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: b
      complex(dpc), dimension(:), intent(inout) :: dbdr
      complex(dpc) :: t11,t12,t21,t22,t31,t33, & 
         k11,k22,k33,s11,s12,s13,s22,s23,a123,a126,a135,a234

      call sph_solid_amat(rr,t11,t12,t21,t22,t31,t33, & 
         k11,k22,k33,s11,s12,s13,s22,s23)

      dbdr(1)  = s11*b(7)+s12*(b(8)-b(3))+s13*b(9)-s22*b(4) & 
                 -s23*b(5)-(t11+t22)*b(1)+t31*b(6)
      dbdr(2)  = s11*b(10)+s12*b(11)-s13*(2.0_dp*b(3)+b(8)) & 
                 -s23*b(4)-t21*b(6)-(t11+t33)*b(2)
      dbdr(3)  = s12*b(12)+s13*b(13)+t12*b(4)-t21*b(7)-t31*b(10)
      dbdr(4)  = -k22*b(1)-s11*b(12)+s13*b(14)+t21*(b(3)-b(8)) & 
                 +(t22-t11)*b(4)-t31*b(11)
      dbdr(5)  = -k33*b(2)-s11*b(13)-s12*b(14)-t21*b(9) & 
                 +t31*(2.0_dp*b(3)+b(8))+(t33-t11)*b(5)
      dbdr(6)  = s12*b(10)-s13*b(7)+s22*b(11) & 
                -s23*(b(3)+2.0_dp*b(8))-t12*b(2)-(t22+t33)*b(6)
      dbdr(7)  = k11*b(1)+s22*b(12)+s23*b(13)-t12*(b(3)-b(8)) & 
                 +(t11-t22)*b(7)
      dbdr(8)  = -s12*b(12)+s23*b(14)-t12*b(4)+t21*b(7)
      dbdr(9)  = -k33*b(6)-s12*b(13)-s22*b(14)-t12*b(5) & 
                 -t22*b(9)+t31*b(7)+t33*b(9)
      dbdr(10) = k11*b(2)+s23*b(12)+t12*b(11)+(t11-t33)*b(10)
      dbdr(11) = k22*b(6)-s13*b(12)+t21*b(10)+(t22-t33)*b(11)
      dbdr(12) = -k11*b(4)+k22*b(7)+(t11+t22)*b(12)
      dbdr(13) = -k11*b(5)+k33*b(10)+t12*b(14)+(t11+t33)*b(13)
      dbdr(14) = -k22*b(9)+k33*b(11)+t21*b(13)-t31*b(12) & 
                 +(t22+t33)*b(14)

      return
    end subroutine sph_bvec_solid_derivs


    subroutine sph_bvec_solid_derivs_ng(rr,b,dbdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 1st   !
      ! order 'b-vector' equations in fluid parts of the  !
      ! earth model.                                      !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: b
      complex(dpc), dimension(:), intent(inout) :: dbdr
      complex(dpc) :: t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22

      call sph_solid_amat_ng(rr,t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22)

      dbdr(1) = -s11*b(3)-s12*b(4)-t11*b(1)-t21*b(2)
      dbdr(2) = -s12*b(3)-s22*b(4)-t12*b(1)-t22*b(2)
      dbdr(3) = -k11*b(1)+t11*b(3)+t12*b(4)
      dbdr(4) = -k22*b(2)+t21*b(3)+t22*b(4)

      return
    end subroutine sph_bvec_solid_derivs_ng



    subroutine sph_bvec_solid_derivs_cow(rr,b,dbdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 1st   !
      ! order 'b-vector' equations in fluid parts of the  !
      ! earth model.                                      !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: b
      complex(dpc), dimension(:), intent(inout) :: dbdr
      complex(dpc) :: t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22

      call sph_solid_amat_cow(rr,t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22)

      dbdr(1) = -s11*b(3)-s12*b(4)-t11*b(1)-t21*b(2)
      dbdr(2) = -s12*b(3)-s22*b(4)-t12*b(1)-t22*b(2)
      dbdr(3) = -k11*b(1)+t11*b(3)+t12*b(4)
      dbdr(4) = -k22*b(2)+t21*b(3)+t22*b(4)

      return
    end subroutine sph_bvec_solid_derivs_cow


    


    subroutine sph_bvec_fluid_derivs(rr,b,dbdr)
      !---------------------------------------------------!
      ! This routine returns the derivatives of the 1st   !
      ! order 'b-vector' equations in fluid parts of the  !
      ! earth model.                                      !
      !---------------------------------------------------!
      use nrtype
      implicit none
      real(dp), intent(in) :: rr
      complex(dpc), dimension(:), intent(in) :: b
      complex(dpc), dimension(:), intent(inout) :: dbdr
      complex(dpc) :: t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22

      call sph_fluid_amat(rr,t11,t12,t21,t22,k11,k22, & 
           s11,s12,s22)

      dbdr(1) = -s11*b(3)-s12*b(4)-t11*b(1)-t21*b(2)
      dbdr(2) = -s12*b(3)-s22*b(4)-t12*b(1)-t22*b(2)
      dbdr(3) = -k11*b(1)+t11*b(3)+t12*b(4)
      dbdr(4) = -k22*b(2)+t21*b(3)+t22*b(4)

      return
    end subroutine sph_bvec_fluid_derivs





    subroutine sph_steps
      !-----------------------------------------!
      ! This routine returns an initial guess   !
      ! at the step size used in the spheroidal !
      ! mode numerical integrations             !
      !-----------------------------------------!
      use nrtype; use module_model
      use nrutil
      implicit none
      integer(i4b) :: im


      im=iminloc(alpha)      
      if(real(omega) == 0.0_dp) then
         dr = 1.0_dp/nw
      else
         dr = nw*real(omega)/alpha(im)
         dr = 1.0_dp/dr
      end if

      return
    end subroutine sph_steps


    subroutine sph_steps2(i1)
      !--------------------------------------------!
      ! This routine returns an initial step size  !
      ! for the integration of the toroidal mode   !
      ! system. This guess is not very important,  !
      ! as the step size used will be modified by  !
      ! the numerical integration routine.         !
      !--------------------------------------------!
      use nrtype; use module_model
      use nrutil
      implicit none
      integer(i4b), intent(in) :: i1
      if(real(omega) == 0.0_dp) then
         dr=1.0_dp/nw
      else
         dr=nw*real(omega)/alpha(i1)
         dr=1.0_dp/dr
      end if
      return
    end subroutine sph_steps2



    subroutine sph_solid_start(is,rr,y)
      !----------------------------------------------------!
      ! This routine returns the starting values for the   !
      ! spheroidal mode equations in solid parts of the    !
      ! earth model in terms of spherical Bessel functions !
      !----------------------------------------------------!
      use nrtype; use module_model
      use module_spline; use module_function
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(6,3), intent(out) :: y
      integer(i4b) :: i
      real(dp) :: rrho,rmu,rqmu,ri,gg,rkappa,rqkappa, & 
           sign,rl,rlm
      complex(dpc) :: gamma,cmu,ckappa,calpha,ratio, & 
           arg,cbeta,tmp1,tmp2,tmp3,zt,xi,gamma2,args

      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,is)
      gg      = splint_dis(r,grav,grav_cs,rr,is)
      rmu     = splint_dis(r,mu,mu_cs,rr,is)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,is)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,is)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,is)      
      ! compute parameters
      cmu    = rmu*(1.0_dp+rqmu*lno)
      ckappa = rkappa*(1.0_dp+rqkappa*lno)
      calpha = (ckappa+1.3333333333_dp*cmu)/rrho     
      cbeta  = cmu/rrho
      sign   = 1.0_dp
      tmp1   = omega2/cbeta
      tmp2   = omega2+5.3333333_dp*pibigg*rrho
      tmp2   = tmp2/calpha
      tmp3   = 2.666666666_dp*pibigg*zeta*rrho
      tmp3   = tmp3*tmp3/(calpha*cbeta)

      ! loop over first two solutions
      do i=1,2
         gamma  = (tmp1-tmp2)**2+tmp3
         gamma  = 0.5_dp*(tmp1+tmp2)  & 
                 +0.5_dp*sign*sqrt(gamma)
         gamma2 = gamma
         gamma  =sqrt(gamma)
         zt     = 0.75_dp*cbeta*(gamma2-tmp1)/(pibigg*rrho)
         xi     = zt-real(llp1)

         ! compute spherical Bessel function ratio
         arg    = gamma*rr
         args   = gamma2*rr*rr
         call bfs_dum(ll,args,ratio)
         ! compute solutions
         y(1,i) = ll*xi-zt*gamma*rr*ratio
         y(2,i) = zeta*(xi+rr*gamma*ratio)
         y(3,i) = -4.0_dp*pibigg*rrho*rr*zt
         y(4,i) = -rrho*calpha*zt*rr*gamma2 & 
                  +2.0_dp*cmu*ll*(ll-1)*xi/rr & 
                  +2.0_dp*cmu*(2.0_dp*zt+zeta2)*gamma*ratio
         y(5,i)  = rr*gamma2+2.0_dp*(ll-1)*xi/rr & 
                   -2.0_dp*(zt+1.0_dp)*gamma*ratio
         y(5,i) = cmu*zeta*y(5,i)
         y(6,i) = -rrho*(zeta2+llp1*zt)      
         
         sign   = -sign
      end do      

      ! compute the third solution
      y(1,3) = ll
      y(2,3) = zeta
      y(3,3) = (omega2-4.0_dp*pibigg*rrho*ll/3.0_dp)*rr      
      y(4,3) = 2.0_dp*ll*(ll-1)*cmu/rr
      y(5,3) = 2.0_dp*zeta*(ll-1)*cmu/rr
      y(6,3) = ((2*ll+1)*omega2 & 
               -8.0_dp*pibigg*ll*(ll-1)*rrho/3.0_dp)
      y(6,3) = 0.25_dp*y(6,3)/pibigg


      return
    end subroutine sph_solid_start



      
    subroutine sph_solid_start_ng(is,rr,y)
      !----------------------------------------------------!
      ! This routine returns the starting values for the   !
      ! spheroidal mode equations in solid parts of the    !
      ! earth model in terms of spherical Bessel functions !
      ! These formula come from Takeuchi and Saito 1972    !
      ! equations (104) and (105)                          !
      !----------------------------------------------------!
      use nrtype; use module_model
      use module_spline; use module_function
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(4,2), intent(out) :: y
      real(dp) :: rrho,rmu,rqmu,ri,gg,rkappa,rqkappa, & 
           sign,rl,rlm
      complex(dpc) :: gamma,cmu,ckappa,calpha,ratio, & 
           arg,cbeta,tmp1,tmp2,tmp3,zt,xi,gamma2
      complex(dpc), dimension(2) :: jl

      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,is)
      rmu     = splint_dis(r,mu,mu_cs,rr,is)
      rqmu    = splint_dis(r,qmu,qmu_cs,rr,is)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,is)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,is)      
      cmu     = rmu*(1.0_dp+rqmu*lno)
      
      ! compute some parameters
      ckappa = rkappa*(1.0_dp+rqkappa*lno)
      calpha = (ckappa+1.3333333333_dp*cmu)/rrho     
      cbeta  = cmu/rrho
      
      ! first solution
      arg = omega*rr/sqrt(calpha)
      call sbessj2(llp1,arg,ratio)
      jl(1)  = one
      jl(2)  = ratio      
      y(1,1) = ll*jl(1)-arg*jl(2)
      y(2,1) = zeta*jl(1)
      y(3,1) = -rrho*calpha*arg*arg*jl(1) & 
               +2.0_dp*cmu*(ll*(ll-1)*jl(1) & 
               +2.0_dp*arg*jl(2))
      y(3,1) = y(3,1)/rr
      y(4,1) = 2.0_dp*zeta*cmu*((ll-1)*jl(1) & 
               -arg*jl(2))
      y(4,1) = y(4,1)/rr

      ! second solution
      arg = omega*rr/sqrt(cbeta)
      call sbessj2(llp1,arg,ratio)
      jl(1)  = one
      jl(2)  = ratio
      y(1,2) = -ll*llp1*jl(1)
      y(2,2) = -zeta*llp1*jl(1)+arg*jl(2)
      y(3,2) = -2.0_dp*cmu*(ll*(ll*ll-1)*jl(1) & 
               -ll*llp1*arg*jl(2))
      y(3,2) = y(3,2)/rr
      y(4,2) = zeta*cmu*(arg*arg*jl(1) & 
               -2.0_dp*(ll*ll-1)*jl(1) & 
               -2.0_dp*arg*jl(2))
      y(4,2) = y(4,2)/rr



      return
    end subroutine sph_solid_start_ng


      

    subroutine sph_fluid_start(is,rr,y)
      !----------------------------------------------------!
      ! This routine returns the starting values for the   !
      ! spheroidal mode equations in fluid parts of the    !
      ! earth model in terms of spherical Bessel functions !
      !----------------------------------------------------!
      use nrtype; use module_model
      use module_spline; use module_function
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(4,2), intent(out) :: y
      integer(i4b) :: i
      real(dp) :: rrho,ri,gg,rkappa,rqkappa, & 
           sign,rl,rlm
      complex(dpc) :: gamma,cmu,ckappa,calpha,ratio, & 
           arg,tmp1,tmp2,tmp3,zt,xi,gamma2,gamg,args

      ! evaluate cubic splines
      rrho    = splint_dis(r,rho,rho_cs,rr,is)
      gg      = splint_dis(r,grav,grav_cs,rr,is)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,is)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,is)      

      ! compute parameters
      ckappa = rkappa*(1.0_dp+rqkappa*lno)
      calpha = ckappa/rrho
      gamg   = 4.0_dp*pibigg*rrho/3.0_dp
      gamma  = (omega2+4.0_dp*gamg-zeta2*gamg*gamg*omegai2)/calpha
      gamma2 = gamma
      gamma  = sqrt(gamma)
      zt     = -omega2/gamg
      xi     =  zt-llp1
      


      ! compute Bessel functions
      arg = gamma*rr
      args = gamma2*rr*rr
      call bfs_dum(ll,args,ratio)


      ! compute first solution
      y(1,1) = ll*xi-zt*gamma*rr*ratio
      y(2,1) = -4.0_dp*pibigg*rrho*rr*zt
      y(3,1) = -rrho*calpha*zt*rr*gamma2                   
      y(4,1) =  -rrho*(zeta2+llp1*zt)     

      ! compute second solution
      y(1,2) = ll
      y(2,2) = (omega2-4.0_dp*pibigg*rrho*ll/3.0_dp)*rr      
      y(3,2) = 0.0_dp
      y(4,2) = ((2*ll+1)*omega2 & 
               -8.0_dp*pibigg*ll*(ll-1)*rrho/3.0_dp)
      y(4,2) = 0.25_dp*y(4,2)/pibigg



      return
    end subroutine sph_fluid_start

    subroutine sph_fluid_start_ng(is,rr,y)
      !----------------------------------------------------!
      ! This routine returns the starting values for the   !
      ! spheroidal mode equations in solid parts of the    !
      ! earth model in terms of spherical Bessel functions !
      ! These formula come from Takeuchi and Saito 1972    !
      ! equations (104) and (105)                          !
      !----------------------------------------------------!
      use nrtype; use module_model
      use module_spline; use module_function
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(2), intent(out) :: y
      real(dp) :: rrho,ri,gg,rkappa,rqkappa, & 
           sign,rl,rlm
      complex(dpc) :: gamma,cmu,ckappa,calpha,ratio, & 
           arg,tmp1

      rrho    = splint_dis(r,rho,rho_cs,rr,is)
      rkappa  = splint_dis(r,kappa,kappa_cs,rr,is)
      rqkappa = splint_dis(r,qkappa,qkappa_cs,rr,is)      
      ckappa  = rkappa*(1.0_dp+rqkappa*lno)
      calpha  = ckappa/rrho     

      
      ! compute Bessel functions
      arg   = omega*rr/sqrt(calpha)
      
      call sbessj2(llp1,arg,ratio)

      ! compute solution
      y(1) = ll-arg*ratio
      y(2) = -rrho*calpha*arg*arg 
      y(2) = y(2)/rr




      return
    end subroutine sph_fluid_start_ng



    subroutine sph_minor_solid_start(is,rr,m)
      !--------------------------------------------------------!
      ! This routine uses the starting values for the          !
      ! spheroidal mode equations provided by sph_solid_start  !
      ! to compute the starting values for the 3rd order minor !
      ! vector equations in solid parts of the earth model     !
      !--------------------------------------------------------!
      use nrtype
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(14), intent(out) :: m
      complex(dpc), dimension(6,3) :: y
      integer(i4b), dimension(3) :: ind

      call sph_solid_start(is,rr,y)

      m(1)  =  det3(y((/1,2,3/),1),y((/1,2,3/),2),y((/1,2,3/),3))
      m(2)  =  det3(y((/1,2,4/),1),y((/1,2,4/),2),y((/1,2,4/),3))
      m(3)  =  det3(y((/1,2,5/),1),y((/1,2,5/),2),y((/1,2,5/),3))
      m(4)  =  det3(y((/1,2,6/),1),y((/1,2,6/),2),y((/1,2,6/),3))
      m(5)  =  det3(y((/1,3,4/),1),y((/1,3,4/),2),y((/1,3,4/),3))
      m(6)  =  det3(y((/1,3,5/),1),y((/1,3,5/),2),y((/1,3,5/),3))
      m(7)  =  det3(y((/1,4,5/),1),y((/1,4,5/),2),y((/1,4,5/),3))
      m(8)  =  det3(y((/1,4,6/),1),y((/1,4,6/),2),y((/1,4,6/),3))
      m(9)  =  det3(y((/1,5,6/),1),y((/1,5,6/),2),y((/1,5,6/),3))
      m(10) =  det3(y((/2,3,4/),1),y((/2,3,4/),2),y((/2,3,4/),3))
      m(11) =  det3(y((/2,4,5/),1),y((/2,4,5/),2),y((/2,4,5/),3))
      m(12) =  det3(y((/2,4,6/),1),y((/2,4,6/),2),y((/2,4,6/),3))
      m(13) =  det3(y((/3,4,5/),1),y((/3,4,5/),2),y((/3,4,5/),3))
      m(14) =  det3(y((/4,5,6/),1),y((/4,5,6/),2),y((/4,5,6/),3))      


      return
    end subroutine sph_minor_solid_start



    subroutine sph_minor_solid_start_ng(is,rr,m)      
      use nrtype
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(5), intent(out) :: m
      complex(dpc), dimension(4,2) :: y
      
      call sph_solid_start_ng(is,rr,y)      

      m(1) = det2(y((/1,2/),1),y((/1,2/),2))
      m(2) = det2(y((/1,3/),1),y((/1,3/),2))
      m(3) = det2(y((/1,4/),1),y((/1,4/),2))
      m(4) = det2(y((/2,3/),1),y((/2,3/),2))
      m(5) = det2(y((/3,4/),1),y((/3,4/),2))

      return
    end subroutine sph_minor_solid_start_ng


    subroutine sph_minor_fluid_start(is,rr,m)
      !--------------------------------------------------------!
      ! This routine uses the starting values for the          !
      ! spheroidal mode equations provided by sph_fluid_start  !
      ! to compute the starting values for the 2nd order minor !
      ! vector equations in fluid parts of the earth model     !
      !--------------------------------------------------------!
      use nrtype
      implicit none
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rr
      complex(dpc), dimension(5), intent(out) :: m
      complex(dpc), dimension(4,2) :: y


      call sph_fluid_start(is,rr,y)

      
      m(1) = det2(y((/1,2/),1),y((/1,2/),2))
      m(2) = det2(y((/1,3/),1),y((/1,3/),2))
      m(3) = det2(y((/1,4/),1),y((/1,4/),2))
      m(4) = det2(y((/2,3/),1),y((/2,3/),2))
      m(5) = det2(y((/3,4/),1),y((/3,4/),2))

      return
    end subroutine sph_minor_fluid_start


    function det3(y1,y2,y3)
      use nrtype
      implicit none
      complex(dpc) :: det3
      complex(dpc), dimension(3), intent(in) :: y1,y2,y3

      det3 = y1(1)*det2(y2(2:3),y3(2:3)) & 
             -y2(1)*det2(y1(2:3),y3(2:3)) & 
             +y3(1)*det2(y1(2:3),y2(2:3))

      return
    end function det3
      

    function det2(y1,y2)
      use nrtype
      implicit none
      complex(dpc) :: det2
      complex(dpc), dimension(2), intent(in) :: y1,y2

      det2 = y1(1)*y2(2)-y1(2)*y2(1)

      return
    end function det2
    

    subroutine b_solid_start(m,s,b)
      !-------------------------------------------------------!
      ! Given the relevant 3rd order minor vector m and       !
      ! a source vector s, this routine returns the starting  !
      ! value of the 'b-vector' equations in a solid part of  !
      ! the earth model.                                      !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(14), intent(in)  :: m
      complex(dpc), dimension(6), intent(in)  :: s
      complex(dpc), dimension(14), intent(out) :: b

      b(1)  = m(14)*s(3)-m(7)*s(4)-m(11)*s(5)-m(13)*s(6)      
      b(2)  = -m(14)*s(2)-m(8)*s(4)-m(12)*s(5)+m(11)*s(6)     
      b(3)  = m(7)*s(2)+m(8)*s(3)+m(2)*s(5)+m(5)*s(6)
      b(4)  = m(11)*s(2)+m(12)*s(3)-m(2)*s(4)+m(10)*s(6)
      b(5)  = m(13)*s(2)-m(11)*s(3)-m(5)*s(4)-m(10)*s(5)
      b(6)  = m(14)*s(1)-m(9)*s(4)+m(8)*s(5)-m(7)*s(6)
      b(7)  = -m(7)*s(1)+m(9)*s(3)+m(3)*s(5)+m(6)*s(6)
      b(8)  = -m(11)*s(1)-m(8)*s(3)-m(3)*s(4)-m(5)*s(6)
      b(9)  = -m(13)*s(1)+m(7)*s(3)-m(6)*s(4)+m(5)*s(5)
      b(10) = -m(8)*s(1)-m(9)*s(2)+m(4)*s(5)-m(3)*s(6)
      b(11) = -m(12)*s(1)+m(8)*s(2)-m(4)*s(4)+m(2)*s(6)
      b(12) = m(2)*s(1)+m(3)*s(2)+m(4)*s(3)-m(1)*s(6)
      b(13) = m(5)*s(1)+m(6)*s(2)-m(3)*s(3)+m(1)*s(5)
      b(14) = m(10)*s(1)-m(5)*s(2)+m(2)*s(3)-m(1)*s(4)      

      return
    end subroutine b_solid_start



    subroutine b_solid_start_ngcow(m,s,b)
      !-------------------------------------------------------!
      ! Given the relevant 3rd order minor vector m and       !
      ! a source vector s, this routine returns the starting  !
      ! value of the 'b-vector' equations in a solid part of  !
      ! the earth model.                                      !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(5), intent(in)  :: m
      complex(dpc), dimension(4), intent(in)  :: s
      complex(dpc), dimension(4), intent(out) :: b

      b(1) = -m(5)*s(2)-m(2)*s(3)-m(4)*s(4)
      b(2) = m(5)*s(1)-m(3)*s(3)+m(2)*s(4)
      b(3) = m(2)*s(1)+m(3)*s(2)-m(1)*s(4)
      b(4) = m(4)*s(1)-m(2)*s(2)+m(1)*s(3)

      return
    end subroutine b_solid_start_ngcow




    subroutine sph_minors_solid_2_fluid(m1,m2)
      !-------------------------------------------------------!
      ! This routine transforms a 3rd order minor vector in   !
      ! a solid region of the earth model into a 2nd order    !
      ! minor vector in a fluid region of the earth model     !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(14), intent(in) :: m1
      complex(dpc), dimension(5), intent(out) :: m2

      m2(1) = -m1(6)
      m2(2) = -m1(7)
      m2(3) = m1(9)
      m2(4) = -m1(13)
      m2(5) = m1(14)

      return
    end subroutine sph_minors_solid_2_fluid


    subroutine sph_minors_solid_2_fluid_ngcow(m1,m2)
      !-------------------------------------------------------!
      ! This routine transforms a 3rd order minor vector in   !
      ! a solid region of the earth model into a 2nd order    !
      ! minor vector in a fluid region of the earth model     !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(5), intent(in) :: m1
      complex(dpc), dimension(2), intent(out) :: m2

      m2(1) = -m1(3)
      m2(2) = -m1(5)

      return
    end subroutine sph_minors_solid_2_fluid_ngcow



    subroutine sph_minors_fluid_2_solid(m1,m2)
      !-------------------------------------------------------!
      ! This routine transforms a 2nd order minor vector in   !
      ! a fluid region of the earth model into a 3rd  order   !
      ! minor vector in a solid region of the earth model     !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(5), intent(in) :: m1
      complex(dpc), dimension(14), intent(out) :: m2

      m2(1)  = m1(1)
      m2(2)  = m1(2)
      m2(4)  = m1(3)
      m2(10) = -m1(4)
      m2(12) = -m1(5)    
  
      m2((/3,5,6,7,8,9,11,13,14/)) = zero

      return
    end subroutine sph_minors_fluid_2_solid


    subroutine sph_minors_fluid_2_solid_ngcow(m1,m2)
      !-------------------------------------------------------!
      ! This routine transforms a 2nd order minor vector in   !
      ! a fluid region of the earth model into a 3rd  order   !
      ! minor vector in a solid region of the earth model     !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(2), intent(in) :: m1
      complex(dpc), dimension(5), intent(out) :: m2

      m2(1) = m1(1)
      m2(4) = -m1(2)

      m2((/2,3,5/)) = zero

      return
    end subroutine sph_minors_fluid_2_solid_ngcow




    subroutine sph_delta_solid(m1,m2,delta)
      !------------------------------------------------------!
      ! Given the relevent 3rd order minor vectors, this     !
      ! routine returns the value of the secular determinant !
      ! of the system.                                       !
      !------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(14), intent(in) :: m1,m2
      complex(dpc), intent(out) :: delta

      delta = -m1(14)*m2(1)+2.0_dp*m1(7)*m2(2)      & 
              +2.0_dp*m1(11)*m2(3)+m1(13)*m2(4)      & 
              +2.0_dp*m1(8)*m2(5)+m1(12)*m2(6)       & 
              -2.0_dp*m1(2)*m2(7)-2.0_dp*m1(5)*m2(8) & 
              -m1(10)*m2(9)+m1(9)*m2(10)             & 
              -2.0_dp*m1(3)*m2(11)-m1(6)*m2(12)      & 
              -m1(4)*m2(13)+m1(1)*m2(14)

      return
    end subroutine sph_delta_solid




    subroutine sph_delta_solid_ngcow(m1,m2,delta)
      !------------------------------------------------------!
      ! Given the relevent 3rd order minor vectors, this     !
      ! routine returns the value of the secular determinant !
      ! of the system.                                       !
      !------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(5), intent(in) :: m1,m2
      complex(dpc), intent(out) :: delta

      delta = m1(5)*m2(1)+2.0_dp*m1(2)*m2(2) & 
              +m1(4)*m2(3)+m1(3)*m2(4)+m1(1)*m2(5)                  

      return
    end subroutine sph_delta_solid_ngcow


    subroutine sph_delta_fluid(m1,m2,delta)
      !------------------------------------------------------!
      ! Given the relevent 2nd order minor vectors, this     !
      ! routine returns the value of the secular determinant !
      ! of the system.                                       !
      !------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(5), intent(in) :: m1,m2
      complex(dpc), intent(out) :: delta

      delta = m1(5)*m2(1)+2.0_dp*m1(2)*m2(2) & 
              +m1(4)*m2(3)+m1(3)*m2(4)+m1(1)*m2(5)            

      return
    end subroutine sph_delta_fluid



    subroutine sph_delta_fluid_ngcow(m1,m2,delta)
      !------------------------------------------------------!
      ! Given the relevent 2nd order minor vectors, this     !
      ! routine returns the value of the secular determinant !
      ! of the system.                                       !
      !------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(2), intent(in) :: m1,m2
      complex(dpc), intent(out) :: delta

      delta=m2(2)*m1(1)-m2(1)*m1(2)

      return
    end subroutine sph_delta_fluid_ngcow



    subroutine sph_ysol(m,b,y)
      !-------------------------------------------------------!
      ! Given the relevant minor vector and b-vector, this    !
      ! routine returns the y-vector solution to the system   !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(14), intent(in) :: m,b
      complex(dpc), dimension(6), intent(out) :: y
      

      y(1) = m(3)*b(3)+m(1)*b(6)+m(2)*b(7)+2.0_dp*m(3)*b(8)+m(4)*b(9) & 
             +m(5)*b(10)+m(6)*b(11)+m(7)*b(12)+m(8)*b(13)+m(9)*b(14)
      y(2) = -m(1)*b(2)-m(2)*(2.0_dp*b(3)+b(8))-m(3)*b(4)-m(4)*b(5) & 
             -m(5)*b(11)-m(8)*b(14)+m(10)*b(10)+m(11)*b(12) & 
             +m(12)*b(13)
      y(3) = m(1)*b(1)-m(2)*b(9)+m(3)*b(5)+m(5)*(b(8)-b(3)) & 
             -m(6)*b(4)+m(7)*b(14)-m(10)*b(7)-m(11)*b(13) & 
             +m(13)*b(12)
      y(4) = m(2)*b(1)+m(5)*b(2)-m(7)*b(4)-m(8)*b(5) & 
             +m(10)*b(6)-m(11)*(b(3)+2.0_dp*b(8))-m(12)*b(9) & 
             -m(13)*b(11)+m(14)*b(14)
      y(5) = m(3)*b(1)-m(5)*b(6)+m(6)*b(2)+m(7)*(2.0_dp*b(3) & 
             +b(8))+m(8)*b(9)-m(9)*b(5)+m(11)*b(7)+m(13)*b(10) & 
             -m(14)*b(13)
      y(6) = m(2)*b(6)-m(3)*b(2)+m(4)*b(1)+m(7)*b(11) & 
             +m(8)*(b(3)-b(8))+m(9)*b(4)-m(11)*b(10)  & 
             +m(12)*b(7)+m(14)*b(12)

      return
    end subroutine sph_ysol


    

    subroutine sph_ysol_ngcow(m,b,y)
      !-------------------------------------------------------!
      ! Given the relevant minor vector and b-vector, this    !
      ! routine returns the y-vector solution to the system   !
      !-------------------------------------------------------!
      use nrtype
      implicit none
      complex(dpc), dimension(5), intent(in) :: m
      complex(dpc), dimension(4) :: b
      complex(dpc), dimension(4), intent(out) :: y

      y(1) = m(1)*b(2)+m(2)*b(3)+m(3)*b(4)
      y(2) = -m(1)*b(1)+m(4)*b(3)-m(2)*b(4)
      y(3) = -m(2)*b(1)-m(4)*b(2)+m(5)*b(4)
      y(4) = -m(3)*b(1)+m(2)*b(2)-m(5)*b(3)      

      return
    end subroutine sph_ysol_ngcow

    


!=============================================================!
!               numerical integration routines                !
!=============================================================!


    subroutine yint(i1,i2,r1,r2,y,derivs,scale)
      !---------------------------------------------------------!
      ! This routine integrates the given system of equations   !
      ! through the earth model between the specified  limits.  !
      !---------------------------------------------------------!
      use nrtype; use module_model
      implicit none
      integer(i4b), intent(in) :: i1,i2
      real(dp), intent(in) :: r1,r2
      complex(dpc), dimension(:), intent(inout) :: y
      integer(i4b), intent(inout), optional  :: scale
      interface
         subroutine derivs(x,y,dydx)
           use nrtype
           implicit none
           real(dp), intent(in) :: x
           complex(dpc), dimension(:), intent(in) :: y
           complex(dpc), dimension(:), intent(inout) :: dydx
         end subroutine derivs
      end interface
      logical(lgt) :: up
      integer(i4b) :: i,j,num,n
      real(dp) :: rr,rr1,rr2,h_use,h_did,h_next,h_min,drr
      complex(dpc), dimension(size(y,1)) :: dydx,yscal

      ! determine direction of integration
      if(r2 > r1) then
         up = .true.
      else if(r2 < r1) then
         up = .false.
      else
         return
      end if

      num = size(y,1)



      !--------------------------------------!
      !         upwards integration          !
      !--------------------------------------!
      if(up) then
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
               call derivs(rr,y,dydx)
               call rkck(y,dydx,rr,drr,y,derivs)
               if(present(scale)) then
                  if(all(abs(y) < tiny)) then
                     scale=scale+1
                     y=y*huge
                  end if
                  if(any(abs(y) > huge)) then
                     scale=scale-1
                     y=y*tiny
                  end if
               end if
            end do
         end do
         

      
      else
         
         !-----------------------------------!
         !        downwards integration      !
         !-----------------------------------! 
         do  i=i1,i2,-1
         
            if(i == i1) then
               rr1=r1
            else
               rr1=r(i+1)
            end if
            if(i == i2) then
               rr2=r2
            else
               rr2=r(i)
            end if
            if(rr2 == rr1) cycle
            call set_layer(i)
            rr=rr1
            drr = steps(i)
            n = 2+floor(abs(rr2-rr1)/drr)
            drr = (rr2-rr1)/(n-1)
            do j = 1,n-1
               rr = rr1+(j-1)*drr
               call derivs(rr,y,dydx)
               call rkck(y,dydx,rr,drr,y,derivs)
               if(present(scale)) then
                  if(all(abs(y) < tiny)) then
                     scale=scale+1
                     y=y*huge
                  end if
                  if(any(abs(y) > huge)) then
                     scale=scale-1
                     y=y*tiny
                  end if
               end if
            end do   
         end do
         
      end if

      return
    end subroutine yint
    
    
    

    subroutine rk4(y,dydx,x,h,yout,derivs)
      use nrtype; use nrutil, only : assert_eq
      implicit none
      complex(dpc), dimension(:), intent(in) :: y,dydx
      real(dp), intent(in) :: x,h
      complex(dpc), dimension(:), intent(out) :: yout
      interface
         subroutine derivs(x,y,dydx)
           use nrtype
           implicit none
           real(dp), intent(in) :: x
           complex(dpc), dimension(:), intent(in) :: y
           complex(dpc), dimension(:), intent(out) :: dydx
         end subroutine derivs
      end interface
      ! Given values for the n variables y and their derivatives dydx known
      ! at x, use the fourth-order Runge-Kutta method to advance the solution
      ! over an intervalh and return the incremented variables as yout, which
      ! need not be a distinct array from y. y, dydx and yout are all of
      ! length n. The user supplies the subroutine derivs, which returns the
      ! derivatives dydx at x.
      ! The routine has been modified putting calculations into double 
      ! precision and allowing the dependent variables to be complex valued. 
      integer(i4b) :: ndum
      real(dp) :: h6,hh,xh
      complex(dpc), dimension(size(y)) :: dym,dyt,yt
      ndum=assert_eq(size(y),size(dydx),size(yout),'rk4')
      hh=h*0.5_dp
      h6=h/6.0_dp
      xh=x+hh
      yt=y+hh*dydx
      call derivs(xh,yt,dyt)
      yt=y+hh*dyt
      call derivs(xh,yt,dym)
      yt=y+h*dym
      dym=dyt+dym
      call derivs(x+h,yt,dyt)
      yout=y+h6*(dydx+dyt+2.0_dp*dym)
    return
  end subroutine rk4
  


  subroutine rkck(y,dydx,x,h,yout,derivs)
    use nrtype; use nrutil, only : assert_eq
    implicit none
    complex(dpc), dimension(:), intent(in) :: y,dydx
    real(dp), intent(in) :: x,h
    complex(dpc), dimension(:), intent(out) :: yout
    interface
       subroutine derivs(x,y,dydx)
         use nrtype
         implicit none
         real(dp), intent(in) :: x
         complex(dpc), dimension(:), intent(in) :: y
         complex(dpc), dimension(:), intent(inout) :: dydx
       end subroutine derivs
    end interface
    ! Given values for the n variables y and their derivatives dydx known
    ! at x, use the fifth order Cash-Karp Runge-Kutta method  to advance
    ! the solution over an interval h and return the incremented variables
    ! as yout. Also return as estimate of the local truncation error in
    ! yout using the embedded fourth order method. The user supplies the 
    ! subroutine derivs, which returns derivatives dydx at x.
    integer(i4b) :: ndum
    complex(dpc), dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    real(dp), parameter :: a2=0.2_dp, a3=0.3_dp, a4=0.6_dp, a5=1.0_dp, &
         a6=0.875_dp, b21=0.2_dp, b31=3.0_dp/40.0_dp, b32=9.0_dp/40.0_dp, &
         b41=0.3_dp, b42=-0.9_dp, b43=1.2_dp, b51=-11.0_dp/54.0_dp, &
         b52=2.5_dp, b53=-70.0_dp/27.0_dp, b54=35.0_dp/27.0_dp, & 
         b61=1631.0_dp/55296.0_dp, b62=175.0_dp/512.0_dp, & 
         b63=575.0_dp/13824.0_dp, b64=44275.0_dp/110592.0_dp, &
         b65=253.0_dp/4096.0_dp, c1=37.0_dp/378.0_dp, &
         c3=250.0_dp/621.0_dp, c4=125.0_dp/594.0_dp, &
         c6=512.0_dp/1771.0_dp, dc1=c1-2825.0_dp/27648.0_dp, &
         dc3=c3-18575.0_dp/48384.0_dp, dc4=c4-13525.0_dp/55296.0_dp, &
         dc5=-277.0_dp/14336.0_dp, dc6=c6-0.25_dp
    ndum=assert_eq(size(y),size(dydx),size(yout),'rkck')
    ytemp=y+b21*h*dydx
    call derivs(x+a2*h,ytemp,ak2)
    ytemp=y+h*(b31*dydx+b32*ak2)
    call derivs(x+a3*h,ytemp,ak3)
    ytemp=y+h*(b41*dydx+b42*ak2+b43*ak3)
    call derivs(x+a4*h,ytemp,ak4)
    ytemp=y+h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
    call derivs(x+a5*h,ytemp,ak5)
    ytemp=y+h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)
    call derivs(x+a6*h,ytemp,ak6)
    yout=y+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6)
    ! estimate error as difference between fourth and fifth order methods
    return
  end subroutine rkck


!--------------------------------------------------------!
!              utility and storage routines              !
!--------------------------------------------------------!
    
    
    subroutine set_l(l_in)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: l_in
      ll=l_in
      llp1=ll+1
      zeta2=real(ll*llp1)
      zeta=sqrt(zeta2)
      zeta2m2=zeta2-2.0_dp
      if(zeta2m2 > 0.0_dp) then
         szeta2m2=sqrt(zeta2m2)
      else
         szeta2m2=0.0_dp
      end if
      nu=real(2*ll+1)/fourpi_d
      nu=sqrt(nu)
      return
    end subroutine set_l
  
    subroutine set_omega(omega_in) 
      use nrtype; use module_model
      implicit none
      complex(dpc), intent(in) :: omega_in
      omega   = omega_in
      omega2  = omega*omega
      omegai  = one/omega
      omegai2 = omegai*omegai
      lno     = 2.0_dp*log(ii*omega/tref)/pi_d
      lno     = ats*lno
      return
    end subroutine set_omega

    subroutine set_ats(ats_in)
      use nrtype
      implicit none
      real(dp), intent(in) :: ats_in
      ats  = ats_in 
      return
    end subroutine set_ats
    
    subroutine set_layer(il_in)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: il_in
      il=il_in
      return
    end subroutine set_layer


    function steps(i1)
      use nrtype
      use module_model
      implicit none
      real(dp) :: steps
      integer(i4b), intent(in) :: i1
      real(dp) :: tmp

      if(i1 < nknot) then
         tmp = max(alpha(i1),alpha(i1+1))
      else 
         tmp = alpha(i1)
      end if      
      steps = nw*real(omega)/tmp
      steps = 1.0_dp/steps      

      return
    end function steps

    
  end module module_int
