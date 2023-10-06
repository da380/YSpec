module module_model

  !===============================================================!
  ! This module contains routines that read in a given spherical  !
  ! earth model file and process it for use in the rest of the    !
  ! program. This module is also used for the storage of the      !
  ! model  file during run time.                                  !
  !===============================================================!


  use nrtype
  implicit none


  !===============================================================!
  !           normalization parameters for the model              !
  !===============================================================!
  real(dp), parameter :: r_norm=6.371e+6_dp                       !
  real(dp), parameter :: rho_norm=5.515e+3_dp                     ! 
  real(dp), parameter :: big=1.e+40_dp                            !
  real(dp), save      :: vel_norm                                 !
  real(dp), save      :: acl_norm                                 !
  real(dp), save      :: fre_norm                                 !
  real(dp), save      :: t_norm                                   !
  real(dp), save      :: con_norm                                 !
  real(dp), save      :: grav_norm                                !
  real(dp), save      :: moment_norm                              !
  real(dp), save      :: pot_norm                                 !
  real(dp), save      :: pibigg                                   !
  real(dp), parameter :: eptt=5.0_dp                              !
  !---------------------------------------------------------------!


  !===============================================================!
  !            arrays holding the model properties                !
  !===============================================================!
  ! The parameters are:                                           !
  ! ifanis == logical variable, true if the model is transversely !
  !           isotropic.                                          !
  ! ocean  == logical variable, true if the model has an ocean    !
  ! fic    == logical variable, true if the model has a fluid     !
  !           inner core                                          !
  ! nlayer == the number of distinct layers in the model, that    !
  !           is, the number of regions in the model separated by !
  !           parameter discontinities.                           !
  ! nknot  == the total number of radial knots in the model.      !
  ! nic    == the index of the solid side of the ICB              !
  ! noc    == the index of the fluid side of the CMB              !
  ! nsl    == the index of the solid side of the seafloor: set    !
  !           to nknot if there is no ocean present               !
  ! layer_index == an array containing the botom and top indices  !
  !                of a given layer in the model. The dimensions  !
  !                of the array are nlayer-by-2, and the storage  !
  !                format is such that:                           !
  !                layer_index(n,1) == botom index of nth layer,  !
  !                layer_index(n,2) == top index of nth layer.    !
  ! tref   == the reference period for the attenuation in the     !
  !           model.                                              !
  ! r      == array of dimension nknot containing the position of !
  !           each radial knot in the model.                      !
  ! rho    == array of dimension nknot containing the density at  !
  !           each radial knot in the model.                      !
  ! rho_cs == array of dimension knot containing the cubic spline !
  !           coefficient of density at each radial knot in the   !
  !           model.                                              !
  ! acon   == array of dimension nknot containing the A modulus   !
  !           at each radial knot in the model.                   !
  ! acon_cs ==  array of dimension nknot containing the cubic     !
  !             spline coefficient of the A modulus at each       !
  !             radial knot in the  model.                        !
  ! ccon   == array of dimension nknot containing the C modulus   !
  !           at each radial knot in the model.                   !
  ! ccon_cs ==  array of dimension nknot containing the cubic     !
  !             spline coefficient of the C modulus at each       !
  !             radial knot in the  model.                        !
  ! fcon   == array of dimension nknot containing the F modulus   !
  !           at each radial knot in the model.                   !
  ! fcon_cs ==  array of dimension knot containing the cubic      !
  !             spline coefficient of the F modulus at each       !
  !             radial knot in the  model.                        !
  ! lcon   == array of dimension nknot containing the L modulus   !
  !           at each radial knot in the model.                   !
  ! lcon_cs ==  array of dimension knot containing the cubic      !
  !             spline coefficient of the L modulus at each       !
  !             radial knot in the  model.                        !
  ! ncon   == array of dimension nknot containing the N modulus   !
  !           at each radial knot in the model.                   !
  ! ncon_cs ==  array of dimension knot containing the cubic      !
  !             spline coefficient of the N modulus at each       !
  !             radial knot in the  model.                        !
  ! kappa   == array of dimension nknot containing the bulk       !
  !            modulus at each radial knot in the model.          !
  ! kappa_cs == array of dimension knot containing the cubic      !
  !             spline coefficient of the bulk modulus at each    !
  !             radial knot in the  model.                        !
  ! mu       == array of dimension nknot containing the shear     !
  !             modulus at each radial knot in the model.         !
  ! mu_cs    == array of dimension knot containing the cubic      !
  !             spline coefficient of the shear modulus at each   !
  !             radial knot in the  model.                        !
  ! qkappa    == array of dimension nknot containing q_{kappa} at !
  !              each radial knot in the model.                   !
  ! qkappa_cs == array of dimension nknot containing the cubic    !
  !              spline coefficients for q_{kappa} at each radial !
  !              knot in the model.                               !
  ! qmu       == array of dimension nknot containing q_{mu} at    !
  !              each radial knot in the model.                   !
  ! qmu_cs    == array of dimension nknot containing the cubic    !
  !              spline coefficients for q_{mu} at each radial    !
  !              knot in the model.                               !
  ! grav      == array of dimension nknot containing the gravita- !
  !              tional acceleration at each radial knot in the   !
  !              model.                                           !
  ! grav_cs   == array of dimension nknot containing the cubic    !
  !              spline coffieient of the gravitational accelera- !
  !              tion at each radial knot in the model.           !
  !---------------------------------------------------------------!
  logical(lgt), save                        :: ifanis             !
  logical(lgt), save                        :: ocean              !
  logical(lgt), save                        :: fic                !
  integer(i4b), save                        :: nlayer             !
  integer(i4b), save                        :: nknot              !
  integer(i4b), save                        :: nic                !
  integer(i4b), save                        :: noc                !
  integer(i4b), save                        :: nsl                !
  integer(i4b), dimension(:,:), &                                 !
       allocatable , save                   :: layer_index        !
  real(dp), save                            :: tref               !
  real(dp), dimension(:), allocatable, save :: r                  !
  real(dp), dimension(:), allocatable, save :: rho                !
  real(dp), dimension(:), allocatable, save :: rho_cs             !
  real(dp), dimension(:), allocatable, save :: acon               !
  real(dp), dimension(:), allocatable, save :: acon_cs            !
  real(dp), dimension(:), allocatable, save :: ccon               !
  real(dp), dimension(:), allocatable, save :: ccon_cs            !
  real(dp), dimension(:), allocatable, save :: fcon               !
  real(dp), dimension(:), allocatable, save :: fcon_cs            !
  real(dp), dimension(:), allocatable, save :: lcon               !
  real(dp), dimension(:), allocatable, save :: lcon_cs            !
  real(dp), dimension(:), allocatable, save :: ncon               !
  real(dp), dimension(:), allocatable, save :: ncon_cs            !
  real(dp), dimension(:), allocatable, save :: kappa              !
  real(dp), dimension(:), allocatable, save :: kappa_cs           !
  real(dp), dimension(:), allocatable, save :: mu                 !
  real(dp), dimension(:), allocatable, save :: mu_cs              !
  real(dp), dimension(:), allocatable, save :: lambda             !
  real(dp), dimension(:), allocatable, save :: qkappa             !
  real(dp), dimension(:), allocatable, save :: qkappa_cs          !
  real(dp), dimension(:), allocatable, save :: qmu                !
  real(dp), dimension(:), allocatable, save :: qmu_cs             !
  real(dp), dimension(:), allocatable, save :: grav               !
  real(dp), dimension(:), allocatable, save :: grav_cs            !
  real(dp), dimension(:), allocatable, save :: beta               !
  real(dp), dimension(:), allocatable, save :: alpha              !
  real(dp), dimension(:), allocatable, save :: rat                !
  real(dp), dimension(:), allocatable, save :: rat_cs             !
  real(dp), dimension(:), allocatable, save :: xlam               !
  real(dp), dimension(:), allocatable, save :: xlam_cs            !
  real(dp), dimension(:), allocatable, save :: xa2                !
  real(dp), dimension(:), allocatable, save :: xa2_cs             ! 
  !---------------------------------------------------------------!
  
  !===============================================================!
  ! arrays used in numerical quadrature:                          !
  !===============================================================!
  real(dp), dimension(:), allocatable, save,  &                   !
       private :: xint,intg,qint                                  !
  !---------------------------------------------------------------!



  contains


    subroutine read_model(io)                                     
      !===========================================================!
      ! given a unit number to an open model file, this           !
      ! routine reads in the model file and then puts the         !
      ! model parameters into the required form for the           !
      ! program.                                                  !
      !===========================================================!
      use nrtype
      use module_spline
      implicit none
      ! Inputs: 
      integer(i4b), intent(in) :: io
      ! local variables: 
      character(len=256) :: ctmp
      integer(i4b) :: i,j,ians
      real(dp) :: fint


      ! read in the header to the model file
      read(io,*) ctmp            
      read(io,*) ians,tref
      tref=twopi_d/tref
      tref=tref/fre_norm
      if(ians == 1) then
         ifanis = .true.
      else
         ifanis = .false.
      end if
      read(io,*) nknot,nic,noc
      
      ! allocate the arrays needed for the model
      call allocate_model(nknot)

      ! read in the model file
      if(ifanis) then
         do i=1,nknot
            read(io,*) r(i),rho(i),ccon(i),lcon(i), & 
                 qkappa(i),qmu(i),acon(i),ncon(i),fcon(i)
         end do
      else
         do i=1,nknot
            read(io,*) r(i),rho(i),ccon(i),lcon(i), & 
                 qkappa(i),qmu(i)
            acon(i)=ccon(i)
            ncon(i)=lcon(i)
            fcon(i)=1.0_dp
         end do
      end if



      ! check for an ocean
      if(lcon(nknot) == 0.0_dp) then
         ocean = .true. 
      else
         ocean = .false. 
      end if
      if(ocean) then
         do i=nknot-1,1,-1
            if(lcon(i) > 0.0_dp) then
               nsl=i
               exit
            end if
         end do
      else
         nsl=nknot
      end if

      ! check for a fluid inner core
      if(lcon(1) == 0.0_dp) then
         fic = .true.
      else
         fic = .false.
      end if


      ! normalize the parameters and convert velocities
      ! to elastic constants (see section 8.9 of D&T)
      r=r/r_norm
      rho=rho/rho_norm
      acon=acon/vel_norm
      ccon=ccon/vel_norm
      lcon=lcon/vel_norm
      ncon=ncon/vel_norm
      acon=rho*acon**2
      ccon=rho*ccon**2
      lcon=rho*lcon**2
      ncon=rho*ncon**2
      fcon=(acon-2.0_dp*lcon)*fcon



      ! compute equivalent kappa and mu for the model
      ! (again, see section 8.9 of D&T)
      kappa=ccon+4.0_dp*acon-4.0_dp*ncon+4.0_dp*fcon
      kappa=kappa/9.0_dp
      mu=ccon+acon+6.0_dp*lcon+5.0_dp*ncon-2.0_dp*fcon
      mu=mu/15.0_dp
      lambda = kappa-2.0_dp*mu/3.0_dp

      ! compute p and s wave speeds
      beta=sqrt(mu/rho)
      alpha=sqrt((kappa+1.3333333333_dp*mu)/rho)

      

      ! convert from big Q to little q
      where(qkappa > 0.0_dp) qkappa=1.0_dp/qkappa
      where(qmu > 0.0_dp) qmu=1.0_dp/qmu


      ! compute some attenuation parameters as in minos
      rat = 4.0_dp*mu/(3.0_dp*(lambda+2.0_dp*mu))
      xlam = (1.0_dp-rat)*qkappa-0.5_dp*rat*qmu
      xlam = xlam/(1.0_dp-1.5_dp*rat)
      xa2 = (1.0_dp-rat)*qkappa+rat*qmu

      
      ! determine the number of layers in the model
      nlayer=1
      do i=2,nknot
         if(r(i) == r(i-1)) nlayer=nlayer+1
      end do
      allocate(layer_index(nlayer,2))
      layer_index(1,1)=1
      layer_index(nlayer,2)=nknot
      j=0
      do i=2,nknot
         if(r(i) == r(i-1)) then
            j=j+1
            layer_index(j,2)=i-1
            layer_index(j+1,1)=i
         end if
      end do


      ! compute some of the cubic spline coefficients 
      do i=1,nlayer
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              rho(layer_index(i,1):layer_index(i,2)), & 
              big,big,rho_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              acon(layer_index(i,1):layer_index(i,2)), & 
              big,big,acon_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              ccon(layer_index(i,1):layer_index(i,2)), & 
              big,big,ccon_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              fcon(layer_index(i,1):layer_index(i,2)), & 
              big,big,fcon_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              lcon(layer_index(i,1):layer_index(i,2)), & 
              big,big,lcon_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              ncon(layer_index(i,1):layer_index(i,2)), & 
              big,big,ncon_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              kappa(layer_index(i,1):layer_index(i,2)), & 
              big,big,kappa_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              mu(layer_index(i,1):layer_index(i,2)), & 
              big,big,mu_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              qkappa(layer_index(i,1):layer_index(i,2)), & 
              big,big,qkappa_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              qmu(layer_index(i,1):layer_index(i,2)), & 
              big,big,qmu_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              rat(layer_index(i,1):layer_index(i,2)), & 
              big,big,rat_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              xlam(layer_index(i,1):layer_index(i,2)), & 
              big,big,xlam_cs(layer_index(i,1):layer_index(i,2)))
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              xa2(layer_index(i,1):layer_index(i,2)), & 
              big,big,xa2_cs(layer_index(i,1):layer_index(i,2)))
      end do


      ! compute gravity using the radial integral formula
      allocate(xint(nknot),intg(nknot),qint(nknot))
      xint=r
      intg=r**2*rho
      do i=1,nlayer
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              intg(layer_index(i,1):layer_index(i,2)), & 
              big,big,qint(layer_index(i,1):layer_index(i,2)))
      end do
      grav(1)=0.0_dp
      fint=0.0_dp
      do i=1,nknot-1
         call gauslv(r(i),r(i+1),i,fint)
         grav(i+1)=fint
      end do
      where(r > 0.0_dp) grav=(4.0_dp*grav)/r**2
      deallocate(xint,intg,qint)

      ! compute cubic spline coefficients for gravity
      do i=1,nlayer
         call spline(r(layer_index(i,1):layer_index(i,2)), & 
              grav(layer_index(i,1):layer_index(i,2)), & 
              big,big,grav_cs(layer_index(i,1):layer_index(i,2)))
      end do



      return
    end subroutine read_model


    subroutine allocate_model(nknot_in,isign)
      !===========================================================!
      ! given a value for nknot this routine allocates            !
      ! the various arrays of the model. If the optional          !
      ! variable isign is set equal to 1 then the routine         !
      ! deallocates these same arrays.                            !
      !===========================================================!
      use nrtype
      implicit none
      integer(i4b), intent(in) :: nknot_in
      integer(i4b), intent(in), optional :: isign
      logical(lgt) :: ltmp
      ltmp=present(isign)
      if(ltmp) then
         if(isign /= 1) ltmp=.false.
      end if
      if(.not.ltmp) then
         allocate(r(nknot_in),     & 
              rho(nknot_in),       & 
              rho_cs(nknot_in),    & 
              acon(nknot_in),      & 
              acon_cs(nknot_in),   & 
              ccon(nknot_in),      & 
              ccon_cs(nknot_in),   & 
              fcon(nknot_in),      & 
              fcon_cs(nknot_in),   & 
              lcon(nknot_in),      & 
              lcon_cs(nknot_in),   & 
              ncon(nknot_in),      & 
              ncon_cs(nknot_in),   & 
              kappa(nknot_in),     & 
              kappa_cs(nknot_in),  & 
              mu(nknot_in),        & 
              mu_cs(nknot_in),     & 
              qkappa(nknot_in),    & 
              qkappa_cs(nknot_in), &
              qmu(nknot_in),       & 
              qmu_cs(nknot_in),    & 
              rat(nknot_in),       & 
              rat_cs(nknot_in),    & 
              xlam(nknot_in),       & 
              xlam_cs(nknot_in),    & 
              xa2(nknot_in),       & 
              xa2_cs(nknot_in),    & 
              grav(nknot_in),      & 
              grav_cs(nknot_in),   &
              lambda(nknot_in)   ,   &
              beta(nknot_in)   ,   &
              alpha(nknot_in)) 
      else
         deallocate(r,       & 
              rho,           & 
              rho_cs,        &
              acon,          & 
              acon_cs,       & 
              ccon,          & 
              ccon_cs,       & 
              fcon,          & 
              fcon_cs,       & 
              lcon,          & 
              lcon_cs,       & 
              ncon,          & 
              ncon_cs,       & 
              kappa,         & 
              kappa_cs,      & 
              mu,            & 
              mu_cs,         & 
              qkappa,        & 
              qkappa_cs,     &
              qmu,           &
              qmu_cs,        &
              rat,           &
              rat_cs,        &
              xlam,           &
              xlam_cs,        &
              xa2,           &
              xa2_cs,        &
              grav,          & 
              grav_cs,       &
              lambda,          & 
              beta,          &
              alpha,         & 
              layer_index)
      end if
      return
    end subroutine allocate_model


    subroutine set_parameters
      !===========================================================!
      ! this routine computes some of the normalization           !
      ! parameters used in the main program - it should           !
      ! always be the first routine called at run time.           !
      !===========================================================!
      use nrtype
      implicit none
      vel_norm=r_norm*sqrt(pi_d*bigg*rho_norm)
      acl_norm=pi_d*bigg*rho_norm*r_norm
      fre_norm=vel_norm/r_norm
      t_norm=1/fre_norm
      con_norm=vel_norm**2*rho_norm
      grav_norm=pi_d*bigg
      moment_norm=r_norm**5*rho_norm**2*grav_norm
      pot_norm = acl_norm*r_norm
      pibigg=pi_d*bigg/grav_norm
      return
    end subroutine set_parameters


    function intgds(y,iq)
      use nrtype; use module_spline
      implicit none
      real(dp) :: intgds
      real(dp), intent(in) :: y
      integer(i4b), intent(in) :: iq
      intgds=splint_dis(xint,intg,qint,y,iq)
      return
    end function intgds



    subroutine gauslv(r1,r2,iq,fint)
      !============================================================!
      ! Numerical integration by Gauss-Legendre quadrature.        !
      ! This routine has been copied from that in minos.           !
      !============================================================!
      use nrtype
      implicit none
      real(dp), intent(in) :: r1,r2
      integer(i4b), intent(in) :: iq
      real(dp), intent(inout) :: fint
      integer(i4b) :: i,j
      real(dp),dimension(2), parameter :: &
           w=(/.478628670499366_dp,.236926885056189_dp/),&
           x=(/.538469310105683_dp,.906179845938664_dp/)
      real(dp) :: t1,y1,y2
      real(dp) :: vals,vals1,sum
      y1=0.5_dp*(r2+r1)
      y2=0.5_dp*(r2-r1)
      vals=intgds(y1,iq)
      sum=.568888888888889_dp*vals
      do i=1,2
         t1=x(i)*y2
         vals=intgds(y1+t1,iq)
         vals1=intgds(y1-t1,iq)
         sum=sum+w(i)*(vals+vals1)
      end do
      fint=fint+y2*sum
      return
    end subroutine gauslv

    function find_radius_index(rr,up)
      use nrtype
      implicit none
      integer(i4b) :: find_radius_index
      real(dp), intent(in) :: rr
      integer(i4b), intent(in), optional :: up
      integer(i4b) :: iup,i
      if(present(up)) then
         iup=up
      else
         iup=1
      end if
      if(iup /= 1 .and. iup /= -1) stop 'bad input to fri'
      if(rr < 0.0_dp .or. rr > 1.0_dp) then
         find_radius_index=0
         return
      end if
      do i=1,nknot-1
         if(iup == 1) then
            if(rr == r(i) .and. rr == r(i+1)) then
               find_radius_index=i+1
               exit
            end if
            if(rr >= r(i) .and. rr < r(i+1)) then
               find_radius_index=i
               exit
            end if
         else
            if(rr == r(i) .and. rr == r(i+1)) then
               find_radius_index=i
               exit
            end if
            if(rr > r(i) .and. rr <= r(i+1)) then
               find_radius_index=i
               exit
            end if
         end if
      end do
      return
    end function find_radius_index

    subroutine grid_count(drr,rs,rr,ngrid)
      use nrtype
      implicit none
      real(dp), intent(in) :: drr
      real(dp), intent(in) :: rs
      real(dp), intent(in) :: rr
      integer(i4b), intent(out) :: ngrid
      logical(lgt) :: ls,lr
      integer(i4b) :: nr,i,j,k
      real(dp) :: r1,r2,rt,drt
      
      ! work out how many knots will be in the output
      ngrid = 0
      do i = 1,nknot-1
         r1 = r(i)
         r2 = r(i+1)
         if(r1 == r2) cycle
         ls = (rs > r1 .and. rs <= r2)
         lr = (rr > r1 .and. rr <= r2)
         if(ls .and. lr) then
            if(rs < rr) then
               nr = floor((rs-r1)/drr)+2
               ngrid = ngrid+nr
               nr = floor((rr-rs)/drr)+2
               ngrid = ngrid+nr-1
               nr = floor((r2-rr)/drr)+2
               ngrid = ngrid+nr-1
            elseif(rr < rs) then
               nr = floor((rr-r1)/drr)+2
               ngrid = ngrid+nr
               nr = floor((rs-rr)/drr)+2
               ngrid = ngrid+nr-1
               nr = floor((r2-rs)/drr)+2
               ngrid = ngrid+nr-1
            else
               nr = floor((rs-r1)/drr)+2
               ngrid = ngrid+nr
               nr = floor((r2-rs)/drr)+2
               ngrid=ngrid+nr-1
           end if               
         elseif(ls) then
           nr = floor((rs-r1)/drr)+2
            ngrid = ngrid+nr
            nr = floor((r2-rs)/drr)+2
           ngrid=ngrid+nr-1
         elseif(lr) then
            nr = floor((rr-r1)/drr)+2
            ngrid = ngrid+nr
            nr = floor((r2-rr)/drr)+2
            ngrid=ngrid+nr-1
         else
            nr = floor((r2-r1)/drr)+2
            ngrid = ngrid+nr
         end if         
      end do


     return
    end subroutine grid_count

    subroutine make_grid(ngrid,drr,rs,rr,is,ir,igrid,rgrid)
      ! given a maximimum radial step size (dr), 
      ! a source radius (rs), and a receiver radius
      ! (rr), this routine returns an array rgrid
      ! of radial knots in the model with spacing
      ! less than or equal to dr, and such that
      ! rs and rr are included. An integer array
      ! igrid is also returned which specifies the
      ! layer in which each radial knot lies, while
      ! the integers is and ir specify which knots
      ! correspond to the source and reciever points
      use nrtype
      implicit none
      integer(i4b), intent(in) :: ngrid
      real(dp), intent(in) :: drr
      real(dp), intent(in) :: rs
      real(dp), intent(in) :: rr
      integer(i4b), intent(out) :: is
      integer(i4b), intent(out) :: ir
      integer(i4b), dimension(ngrid), intent(out) :: igrid
      real(dp), dimension(ngrid),  intent(out) :: rgrid

      logical(lgt) :: ls,lr
      integer(i4b) :: nr,i,j,k
      real(dp) :: r1,r2,rt,drt



      ! build up the output arrays     
      k = 0
      do i = 1,nknot-1
         r1 = r(i)
         r2 = r(i+1)
         if(r1 == r2) cycle
         ls = (rs > r1 .and. rs <= r2)
         lr = (rr > r1 .and. rr <= r2)
         if(ls .and. lr) then
            if(rs < rr) then
               nr = floor((rs-r1)/drr)+2
               drt = (rs-r1)/(nr-1)
               do j = 1,nr
                  k = k+1
                  rgrid(k) = r1+(j-1)*drt
                  igrid(k) = i
               end do
               is = k
               nr = floor((rr-rs)/drr)+2
               drt = (rr-rs)/(nr-1)
               do j=2,nr
                  k = k+1
                  rgrid(k) = rs+(j-1)*drt
                  igrid(k) = i
               end do
               ir = k
               nr = floor((r2-rr)/drr)+2
               drt = (r2-rr)/(nr-1)
               do j = 2,nr
                  k = k+1
                  rgrid(k) = rr+(j-1)*drt
                  igrid(k) = i
               end do
           elseif(rr < rs) then
               nr = floor((rr-r1)/drr)+2
               drt = (rr-r1)/(nr-1)
               do j = 1,nr
                  k = k+1
                  rgrid(k) = r1+(j-1)*drt
                  igrid(k) = i
               end do
               ir = k
               nr = floor((rs-rr)/drr)+2
              drt = (rs-rr)/(nr-1)
              do j=2,nr
                 k = k+1
                  rgrid(k) = rr+(j-1)*drt
                  igrid(k) = i
              end do
              is = k
              nr = floor((r2-rs)/drr)+2
               drt = (r2-rs)/(nr-1)
               do j = 2,nr
                  k = k+1
                  rgrid(k) = rs+(j-1)*drt
                  igrid(k) = i
              end do
            else
               nr = floor((rs-r1)/drr)+2
               drt = (rs-r1)/(nr-1)
               do j = 1,nr
                  k = k+1
                  rgrid(k) = r1+(j-1)*drt
                  igrid(k) = i
               end do
               is = k
               ir = k               
               nr = floor((r2-rs)/drr)+2
               drt = (r2-rs)/(nr-1)
               do j = 2,nr
                  k = k+1
                  rgrid(k) = rs+(j-1)*drt
                  igrid(k) = i
               end do
            end if               
         elseif(ls) then
            nr = floor((rs-r1)/drr)+2
            drt = (rs-r1)/(nr-1)
            do j = 1,nr
              k = k+1
               rgrid(k) = r1+(j-1)*drt
               igrid(k) = i
            end do
            is = k
            nr = floor((r2-rs)/drr)+2
            drt = (r2-rs)/(nr-1)
            do j = 2,nr
               k = k+1
               rgrid(k) = rs+(j-1)*drt
               igrid(k) = i
            end do
         elseif(lr) then
            nr = floor((rr-r1)/drr)+2
            drt = (rr-r1)/(nr-1)
            do j = 1,nr
               k = k+1
               rgrid(k) = r1+(j-1)*drt
               igrid(k) = i
            end do
            ir = k
            nr = floor((r2-rr)/drr)+2
            drt = (r2-rr)/(nr-1)
            do j = 2,nr
               k = k+1
               rgrid(k) = rr+(j-1)*drt
               igrid(k) = i
            end do
         else
            nr = floor((r2-r1)/drr)+2
            drt = (r2-r1)/(nr-1)
            do j = 1,nr
               k = k+1
               rgrid(k) = r1+(j-1)*drt
               igrid(k) = i
            end do
         end if         
      end do


      
      return
    end subroutine make_grid

end module module_model
