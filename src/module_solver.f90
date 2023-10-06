module module_solver

  use nrtype
  implicit none

  ! switch to turn on or off attenuation:
  !
  ! ats_switch = 1 -- attenuation turned on
  ! ats_switch = 0 -- attenuation turned off
  !
  real(dp), save :: ats_switch 



  contains



    !=======================================================!
    !             Toroidal mode solver routines             !
    !=======================================================!
   


    subroutine tor_solver(ll_in,w,rr,rs,ir,is,ss,yout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rs,rr
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout
      integer(i4b) :: i,i1,scale1,scale2,scale3,ns
      real(dp) :: rt,rlcon,rmu,rqmu
      complex(dpc) :: delta,b,clcon,s1,s2
      complex(dpc), dimension(2) :: y1,y2,y3
      

     
      ! set parameters for the integration
      call set_l(ll_in)
      call set_ats(ats_switch)
      call set_omega(w)
      call tor_steps

      
      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if
     

      
      

      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then
         ! integrate up to the source
         call tor_start_level(ll,real(omega),is,1,i1)
         if(i1 == is .and. i1 > noc+1) i1 = i1-1
         rt = r(i1)
         if(rt > rs) then 
            i1 = rs-1
            rt = r(i1)
         end if
         if(i1 == noc+1) then
            y1(1)=one; y1(2)=zero;
         else
            call tor_ystart(i1,rt,y1)
         end if
         call yint(i1,is,rt,rs,y1,tor_derivs,scale1)

         
         ! carry on integration up to receiver
         y3=y1
         if(rr /= rs) then
            call yint(is,ir,rs,rr,y3,tor_derivs,scale3)
         end if

         ! integrate down to the reciever
         y2(1)=one; y2(2)=zero         
         if(rr /= r(nsl)) then
            call yint(nsl-1,ir,r(nsl),rr,y2,tor_derivs,scale2)
         end if
         ! compute delta
         delta=y2(2)*y3(1)-y2(1)*y3(2)
         delta=one/delta 
         ! put the solution together
         do i = 1,ns
            b = y1(2)*ss(1,i)-y1(1)*ss(2,i)
            yout(:,i) = -delta*y2*b
         end do
         if(scale3 /= 0) then
            yout=yout*huge**scale3
         end if        
      else
         ! integrate up to the reciever
         call tor_start_level(ll,real(omega),ir,1,i1)
         if(i1 == ir .and. i1 > noc+1) i1 = i1-1
         rt = r(i1)
         if(rt > rs) then 
            i1 = rs-1
            rt = r(i1)
         end if
         if(i1 == noc+1) then
            y1(1)=one; y1(2)=zero
         else
            call tor_ystart(i1,rt,y1)
         end if
         call yint(i1,ir,rt,rr,y1,tor_derivs,scale1)

         ! integrate down to the source
         y2(1)=one; y2(2)=zero
         if(rs /= r(nsl)) then
            call yint(nsl-1,is,r(nsl),rs,y2,tor_derivs,scale2)
         end if
         ! carry on integration to the receiver
         y3=y2
         if(rr /= rs) then
            call yint(is,ir,rs,rr,y3,tor_derivs,scale3)
         end if
         ! compute delta
         delta=y3(2)*y1(1)-y3(1)*y1(2)
         delta=one/delta         
         ! put the solution together
         do i=1,ns
            b = y2(2)*ss(1,i)-y2(1)*ss(2,i)
            yout(:,i) = -delta*y1*b
         end do        
         if(scale3 /= 0) then
            yout=yout*huge**scale3
         end if   
      end if

      
      return
    end subroutine tor_solver

    
    subroutine source_vector_tor(l,w,is,rs,mm,svt)
      !---------------------------------------------------------!
      ! returns source discontinuity vector for a given moment  !
      ! tensor source. The moment tensor is input as:           !
      ! m(1) = M_{r,r}                                          !
      ! m(2) = M_{r,theta}                                      !
      ! m(3) = M_{r,phi}                                        !
      ! m(4) = M_{theta,theta}                                  !
      ! m(5) = M_{theta,phi}                                    !
      ! m(6) = M_{phi,phi}                                      !
      !---------------------------------------------------------!
      use nrtype      
      use module_model
      use module_spline
      use module_int
      implicit none

      ! input/output
      integer(i4b), intent(in) :: l
      complex(dpc), intent(in) :: w
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rs
      real(dp), dimension(6), intent(in) :: mm
      complex(dpc), dimension(5,2), intent(out) :: svt

      ! local variables
      real(dp) :: rlcon,rmu,rqmu
      complex(dpc) :: clcon,qll


      call set_l(l)
      call set_ats(ats_switch)
      call set_omega(w)


      ! compute some parameters
      rlcon = splint_dis_r(r,lcon,lcon_cs,rs,is)
      rmu   = splint_dis_r(r,mu,mu_cs,rs,is)
      rqmu  = splint_dis_r(r,qmu,qmu_cs,rs,is)      
      qll = 1.0_dp+rqmu*lno
      clcon = rlcon*qll

      ! initialize the source vector
      svt = 0.0_dp


      if(l == 0) then
      
      else if(l == 1) then
         
         svt(2,1) = 0.5_dp*rs*(-mm(3)+ii*mm(2))/clcon
         svt(4,1) = 0.5_dp*rs*( mm(3)+ii*mm(2))/clcon

      else

         svt(2,1) = 0.5_dp*rs*(-mm(3)+ii*mm(2))/clcon
         svt(4,1) = 0.5_dp*rs*( mm(3)+ii*mm(2))/clcon
         
         svt(1,2) = 0.25_dp*szeta2m2*(-ii*mm(4)+ii*mm(6)+2.0_dp*mm(5))
         svt(5,2) = 0.25_dp*szeta2m2*( ii*mm(4)-ii*mm(6)+2.0_dp*mm(5))

      end if


      svt = nu*svt/(ii*w*rs*rs)




      return
    end subroutine source_vector_tor
    

    




    !=======================================================!
    !              Radial mode solver routines              !
    !=======================================================!
    subroutine rad_solver(w,rr,rs,ir,is,ss,yout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rs,rr
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout
      integer(i4b) :: i,i1,scale1,scale2,scale3,ns
      real(dp) :: rt,rccon,rfcon,rmu,rqmu,rkappa,rqkappa
      complex(dpc) :: delta,b
      complex(dpc), dimension(2) :: y1,y2,y3



      ! set parameters for the integration
      call set_l(0)
      call set_omega(w)
      call set_ats(ats_switch)
      call rad_steps


      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if





      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then
         ! integrate up to the source

         i1 = 1
         rt = min(0.1_dp*dr,0.1_dp*r(2))
         call rad_ystart(i1,rt,y1)
         call yint(i1,is,rt,rs,y1,rad_derivs,scale1)
         ! integrate down to the reciever
         y2(1)=one; y2(2)=zero
         if(rr /= r(nknot)) then
            call yint(nknot-1,ir,r(nknot),rr,y2,rad_derivs,scale2)
         end if
         ! carry on integration to the source
         y3=y2
         if(rs /= rr) then
            call yint(ir,is,rr,rs,y3,rad_derivs,scale3)
         end if

         ! compute delta
         delta=y3(2)*y1(1)-y3(1)*y1(2)
         if(abs(delta) < 1.0e-100_dp) then
            delta=zero
         else
            delta=one/delta
         end if
         ! put the solution together
         do i = 1,ns
            b = y1(2)*ss(1,i)-y1(1)*ss(2,i)
            yout(:,i) = -delta*y2*b
         end do
         if(scale3 /= 0) then
            yout=yout*huge**scale3
         end if     
      else
         ! integrate up to the reciever
         i1 = 1
         rt = min(0.1_dp*dr,0.1_dp*r(2))
         call rad_ystart(i1,rt,y1)
         call yint(i1,ir,rt,rr,y1,rad_derivs,scale1)
         ! integrate down to the source
         y2(1)=one; y2(2)=zero
         if(rs /= r(nsl)) then
            call yint(nsl-1,is,r(nsl),rs,y2,rad_derivs,scale2)
         end if
         ! carry on integration to the receiver
         y3=y2
         if(rr /= rs) then
            call yint(is,ir,rs,rr,y3,rad_derivs,scale3)
         end if
         ! compute delta
         delta=y3(2)*y1(1)-y3(1)*y1(2)
         if(abs(delta) < 1.0e-100_dp) then
            delta=zero
         else
            delta=one/delta
         end if
         ! put the solution together
         do i=1,ns
            b = y2(2)*ss(1,i)-y2(1)*ss(2,i)
            yout(:,i) = -delta*y1*b
         end do        
         if(scale3 /= 0) then
            do i=1,scale3
               yout=yout*huge
            end do
         end if
      end if
      
      return
    end subroutine rad_solver



    subroutine rad_solver_ng(w,rr,rs,ir,is,ss,yout)
      ! solves radial equations with no-gravity
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rs,rr
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout
      integer(i4b) :: i,i1,scale1,scale2,scale3,ns
      real(dp) :: rt,rccon,rfcon,rmu,rqmu,rkappa,rqkappa
      complex(dpc) :: delta,b
      complex(dpc), dimension(2) :: y1,y2,y3



      ! set parameters for the integration
      call set_l(0)
      call set_omega(w)
      call set_ats(ats_switch)
      call rad_steps


      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if





      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then
         ! integrate up to the source

         i1 = 1
         rt = min(0.1_dp*dr,0.1_dp*r(2))
         call rad_ystart_ng(i1,rt,y1)
         call yint(i1,is,rt,rs,y1,rad_derivs_ng,scale1)
         ! integrate down to the reciever
         y2(1)=one; y2(2)=zero
         if(rr /= r(nknot)) then
            call yint(nknot-1,ir,r(nknot),rr,y2,rad_derivs_ng,scale2)
         end if
         ! carry on integration to the source
         y3=y2
         if(rs /= rr) then
            call yint(ir,is,rr,rs,y3,rad_derivs_ng,scale3)
         end if

         ! compute delta
         delta=y3(2)*y1(1)-y3(1)*y1(2)
         if(abs(delta) < 1.0e-100_dp) then
            delta=zero
         else
            delta=one/delta
         end if
         ! put the solution together
         do i = 1,ns
            b = y1(2)*ss(1,i)-y1(1)*ss(2,i)
            yout(:,i) = -delta*y2*b
         end do
         if(scale3 /= 0) then
            yout=yout*huge**scale3
         end if     
      else
         ! integrate up to the reciever
         i1 = 1
         rt = min(0.1_dp*dr,0.1_dp*r(2))
         call rad_ystart_ng(i1,rt,y1)
         call yint(i1,ir,rt,rr,y1,rad_derivs_ng,scale1)
         ! integrate down to the source
         y2(1)=one; y2(2)=zero
         if(rs /= r(nsl)) then
            call yint(nsl-1,is,r(nsl),rs,y2,rad_derivs_ng,scale2)
         end if
         ! carry on integration to the receiver
         y3=y2
         if(rr /= rs) then
            call yint(is,ir,rs,rr,y3,rad_derivs_ng,scale3)
         end if
         ! compute delta
         delta=y3(2)*y1(1)-y3(1)*y1(2)
         if(abs(delta) < 1.0e-100_dp) then
            delta=zero
         else
            delta=one/delta
         end if
         ! put the solution together
         do i=1,ns
            b = y2(2)*ss(1,i)-y2(1)*ss(2,i)
            yout(:,i) = -delta*y1*b
         end do        
         if(scale3 /= 0) then
            do i=1,scale3
               yout=yout*huge
            end do
         end if
      end if
      
      return
    end subroutine rad_solver_ng


    

    subroutine rad_solver_wg(w,rr,rs,ir,is,ss,yout,pout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start      
      implicit none
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rs,rr
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout
      complex(dpc), dimension(:,:), intent(out) :: pout
      integer(i4b) :: i,i1,scale1,scale2,scale3,ns,irr,iss, & 
           ngrid,ig,j,igt,k
      integer(i4b), dimension(:), allocatable :: igrid
      real(dp) :: rt,rccon,rfcon,rmu,rqmu,rkappa,rqkappa      
      real(dp), dimension(:), allocatable :: rgrid
      complex(dpc) :: delta,b,ctmp1,ctmp2
      complex(dpc), dimension(2) :: y1,y2,y3
      complex(dpc), dimension(:,:), allocatable :: y1s,y2s,ys
      complex(dpc), dimension(:), allocatable :: u,u_cs

      ! set parameters for the integration
      call set_l(0)
      call set_omega(w)
      call set_ats(ats_switch)
      call rad_steps


      
      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if

      ! build the grid for storing the solutions
      call grid_count(dr,rs,rr,ngrid)
      allocate(rgrid(ngrid),igrid(ngrid))
      call make_grid(ngrid,dr,rs,rr,iss,irr,igrid,rgrid)
      allocate(y1s(2,ngrid),y2s(2,ngrid),ys(2,ngrid))

      
      
      ! integrate up to the source
      i1 = 1
      rt = min(0.1_dp*dr,0.1_dp*rgrid(2))
      call rad_ystart(i1,rt,y1)
      call yint(i1,igrid(2),rt,rgrid(2),y1,rad_derivs)
      y1s(:,2) = y1(:)
      do ig = 2,iss-1
         call yint(igrid(ig),igrid(ig+1),rgrid(ig), & 
              rgrid(ig+1),y1,rad_derivs)
         y1s(:,ig+1) = y1(:)
      end do
      
      ! integrate down to the source
      y2(1) = 1.0_dp
      y2(2) = 0.0_dp
      y2s(:,ngrid) = y2
      do ig = ngrid,iss+1,-1
         call yint(igrid(ig),igrid(ig-1),rgrid(ig), & 
              rgrid(ig-1),y2,rad_derivs)
         y2s(:,ig-1) = y2
      end do
      
      ! calculate delta
      delta=y2(2)*y1(1)-y2(1)*y1(2)
      delta = 1.0_dp/delta
      
      ! loop over the different sources
      do i = 1,ns
         
         ! compute the solution below the source
         b = y2(2)*ss(1,i)-y2(1)*ss(2,i)
         do ig = 1,iss
            ys(:,ig) = -delta*y1s(:,ig)*b
         end do

         ! compute the solution above the source
         b = y1(2)*ss(1,i)-y1(1)*ss(2,i)
         do ig = iss+1,ngrid
            ys(:,ig) = -delta*y2s(:,ig)*b
         end do
         
         ! extract the solution at the receiver
         yout(:,i) = ys(:,irr)
         
         ! calculate r \partial_{r}\phi at the 
         ! receiver
         pout(2,i) = -4.0_dp*pibigg*yout(1,i)
         
         ! calculate phi using radial integral formula
         ! performing the integration using the 
         ! trapeziod rule -- probably should make this
         ! better in the future
         pout(1,i) = 0.0_dp
         do ig = 1,irr-1
            if(ig == 1) then
               ctmp1 = 0.0_dp
               ctmp2 = -4.0_dp*pibigg*ys(1,ig+1)/rgrid(ig+1)
            else
               ctmp1 = -4.0_dp*pibigg*ys(1,ig)/rgrid(ig)
               ctmp2 = -4.0_dp*pibigg*ys(1,ig+1)/rgrid(ig+1)
            end if
            pout(1,i) = pout(1,i) + 0.5_dp*(ctmp1+ctmp2) & 
                 *(rgrid(ig+1)-rgrid(ig))
         end do
         pout(1,i) = rr*pout(1,i)
         
         
      end do
      

      return
    end subroutine rad_solver_wg




    subroutine source_vector_rad(w,is,rs,mm,s0)
      !---------------------------------------------------------!
      ! returns source discontinuity vector for a given moment  !
      ! tensor source. The moment tensor is input as:           !
      ! m(1) = M_{r,r}                                          !
      ! m(2) = M_{r,theta}                                      !
      ! m(3) = M_{r,phi}                                        !
      ! m(4) = M_{theta,theta}                                  !
      ! m(5) = M_{theta,phi}                                    !
      ! m(6) = M_{phi,phi}                                      !
      !---------------------------------------------------------!
      use nrtype
      use module_model
      use module_spline
      use module_int
      implicit none
      ! input/output
      complex(dpc), intent(in) :: w
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rs
      real(dp), dimension(6), intent(in) :: mm
      complex(dpc), dimension(2), intent(out) :: s0
      ! local variables
      real(dp) :: rccon,rfcon,rkappa,rqkappa, & 
           rmu,rqmu,rxa2,rxlam
      complex(dpc) :: cccon,cfcon,qaa,qff



      call set_l(0)
      call set_ats(ats_switch)
      call set_omega(w)


      ! compute some parameters
      rccon    = splint_dis_r(r,ccon,ccon_cs,rs,is)
      rfcon    = splint_dis_r(r,fcon,fcon_cs,rs,is)
      rmu      = splint_dis_r(r,mu,mu_cs,rs,is)
      rkappa   = splint_dis_r(r,kappa,kappa_cs,rs,is)
      rqmu     = splint_dis_r(r,qmu,qmu_cs,rs,is)     
      rqkappa  = splint_dis_r(r,qkappa,qkappa_cs,rs,is)           
      rxa2  = splint_dis_r(r,xa2,xa2_cs,rs,is)           
      rxlam  = splint_dis_r(r,xlam,xlam_cs,rs,is)           

      
      qaa = 1.0_dp+rxa2*lno
      qff = 1.0_dp+rxlam*lno
 
      cccon = rccon*qaa
      cccon = 1.0_dp/cccon
      cfcon = rfcon*qff


      s0(1) = rs*mm(1)*cccon
      s0(2) = -mm(4)-mm(6)+2.0_dp*mm(1)*cfcon*cccon

      s0 = nu*s0/(ii*w*rs*rs)



      return
    end subroutine source_vector_rad
        

    


    !=======================================================!
    !            Spheroidal mode solver routines            !
    !=======================================================!




    subroutine sph_solver(ll_in,w,rr,rs,ir,is,ss,yout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rs,rr
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout      
      integer(i4b) :: i,i1,scale1,scale2,scale3,scale4,iq,ns
      real(dp) :: rt,rsi,rccon,rlcon,rfcon,rmu,rkappa, & 
           rqmu,rqkappa
      complex(dpc) :: delta
      complex(dpc), dimension(5) :: mf1,mf2
      complex(dpc), dimension(14) :: ms1,ms2,ms3,b



      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if

      ! set parameters for the integration
      call set_l(ll_in)
      call set_ats(ats_switch)
      call set_omega(w)
      call sph_steps




      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then
         call sph_start_level(ll,real(omega),is,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= is) then
            i1 = is-1
            rt = r(i1)
         end if




         if(i1 <= nic) then 
            ! start in the inner core            
            call sph_minor_solid_start(i1,rt,ms1)
            ! integrate from start to ICB
            call yint(i1,nic,rt,r(nic),ms1, & 
                 sph_minors_solid_derivs,scale1)
            ! transform from solid to fluid systems
            call sph_minors_solid_2_fluid(ms1,mf1)
            ! integrate from ICB to CMB
            call yint(nic+1,noc,r(nic),r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, & 
                 sph_minors_solid_derivs,scale1)
         else if (i1 > nic .and. i1 <=  noc) then    
            ! start in the outer core
            call sph_minor_fluid_start(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)              
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, &
            sph_minors_solid_derivs,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            ! start in the mantle
            call sph_minor_solid_start(i1,rt,ms1)   
            ! integrate from start to source
            call yint(i1,is,rt,rs,ms1, & 
                 sph_minors_solid_derivs,scale1)        
         end if
         ms3=ms1
         ! carry on the integration to the receiver
         if(rr > rs) then
            call yint(is,ir,rs,rr,ms3, & 
                 sph_minors_solid_derivs,scale3)        
         end if
         ! integrate down to the receiver
         if(ocean) then
            mf2(1)=one; mf2(2:5)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs,scale2)
            call sph_minors_fluid_2_solid(mf2,ms2)
            call yint(nsl,ir,r(nsl), & 
                 rr,ms2,sph_minors_solid_derivs,scale2)
         else
            ms2(1)=one; ms2(2:14)=zero
            call yint(nknot-1,ir,r(nknot), & 
                 rr,ms2,sph_minors_solid_derivs,scale2)
         end if
         ! calculate delta
         call sph_delta_solid(ms3,ms2,delta)
         ! integrate the b-vector from the source to the receiver
         do i=1,ns
            scale4=0         
            call b_solid_start(ms1,ss(:,i),b)
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs,scale4)   
            ! put the solution together
            call sph_ysol(ms2,b,yout(:,i))  
            ! unscale the solution
            yout(:,i)=-yout(:,i)/delta
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      else
         ! integrate up to the receiver 
         call sph_start_level(ll,real(omega),ir,1,i1)
         rt = r(i1)
         if(i1 >= ir) then
            i1 = ir-1
            rt = r(i1)
         end if
         if(i1 <= nic) then
            call sph_minor_solid_start(i1,rt,ms1)
            call yint(i1,nic,rt,r(nic),ms1, & 
                 sph_minors_solid_derivs,scale1)
            call sph_minors_solid_2_fluid(ms1,mf1)
            call yint(nic+1,noc,r(nic),r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)
            call sph_minors_fluid_2_solid(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs,scale1)
         else if (i1 > nic .and. i1 <=  noc) then            
            call sph_minor_fluid_start(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)    
            call sph_minors_fluid_2_solid(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            call sph_minor_solid_start(i1,rt,ms1)            
            call yint(i1,ir,rt,rr,ms1, & 
                 sph_minors_solid_derivs,scale1)                    
         end if
         ms3=ms1
         ! carry on the integration to the source
         if(rs > rr) then
            call yint(ir,is,rr,rs,ms3, & 
                 sph_minors_solid_derivs,scale3)        
         end if
         ! integrate down to the source
         if(ocean) then
            mf2(1)=one; mf2(2:5)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs,scale2)
            call sph_minors_fluid_2_solid(mf2,ms2)
            call yint(nsl,is,r(nsl), & 
                 rs,ms2,sph_minors_solid_derivs,scale2)
         else
            ms2(1)=one; ms2(2:14)=zero
            call yint(nknot,is,r(nknot), & 
                 rs,ms2,sph_minors_solid_derivs,scale2)
        end if
         ! form delta
         call sph_delta_solid(ms3,ms2,delta)
         ! integrate the b-vector from the source to 
         ! the receiver
         do i=1,ns
            scale4=0
            call b_solid_start(ms2,ss(:,i),b)
            b=-b/delta
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs,scale4)     
            ! put the solution together
            call sph_ysol(ms1,b,yout(:,i))
            ! unscale the solution
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      end if
      
      return
    end subroutine sph_solver



    subroutine sph_solver_ng(ll_in,w,rr,rs,ir,is,ss,yout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rr,rs
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout

      integer(i4b) :: i,i1,scale1,scale2,scale3,scale4,iq,ns
      real(dp) :: rt
      complex(dpc) :: delta
      complex(dpc), dimension(4) :: b
      complex(dpc), dimension(2) :: mf1,mf2
      complex(dpc), dimension(5) :: ms1,ms2,ms3


      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if

      ! set parameters for the integration
      call set_l(ll_in)
      call set_ats(ats_switch)
      call set_omega(w)
      call sph_steps

      
      


      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then

         call sph_start_level(ll,real(omega),is,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= is) then
            i1 = is-1
            rt = r(i1)
         end if


         if(i1 <= nic) then 
            ! start in the inner core            
            call sph_minor_solid_start_ng(i1,rt,ms1)
            ! integrate from start to ICB
            call yint(i1,nic,rt,r(nic),ms1, & 
                 sph_minors_solid_derivs_ng,scale1)
            ! transform from solid to fluid systems
            call sph_minors_solid_2_fluid_ngcow(ms1,mf1)
            ! integrate from ICB to CMB
            call yint(nic+1,noc,r(nic),r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)
         else if (i1 > nic .and. i1 <=  noc) then    
            ! start in the outer core
            call sph_fluid_start_ng(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)              
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, &
            sph_minors_solid_derivs_ng,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            ! start in the mantle
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            ! integrate from start to source
            call yint(i1,is,rt,rs,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)        
         end if
         ms3=ms1
         ! carry on the integration to the receiver
         if(rr > rs) then
            call yint(is,ir,rs,rr,ms3, & 
                 sph_minors_solid_derivs_ng,scale3)        
         end if
         ! integrate down to the receiver
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_ng,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,ir,r(nsl), & 
                 rr,ms2,sph_minors_solid_derivs_ng,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,ir,r(nknot), & 
                 rr,ms2,sph_minors_solid_derivs_ng,scale2)
         end if
         ! calculate delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to the receiver
         do i=1,ns
            scale4=0         
            call b_solid_start_ngcow(ms1,ss(:,i),b)
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_ng,scale4)   
            ! put the solution together
            call sph_ysol_ngcow(ms2,b,yout(:,i))  
            ! unscale the solution
            yout(:,i)=yout(:,i)/delta
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      else
         ! integrate up to the receiver 
         call sph_start_level(ll,real(omega),ir,1,i1)
         rt = r(i1)
         if(i1 >= ir) then
            i1 = ir-1
            rt = r(i1)
         end if
         if(i1 <= nic) then
            call sph_minor_solid_start_ng(i1,rt,ms1)
            call yint(i1,nic,rt,r(nic),ms1, & 
                 sph_minors_solid_derivs_ng,scale1)
            call sph_minors_solid_2_fluid_ngcow(ms1,mf1)
            call yint(nic+1,noc,r(nic),r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)
         else if (i1 > nic .and. i1 <=  noc) then            
            call sph_fluid_start_ng(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)    
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_ng,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            call yint(i1,ir,rt,rr,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)                    
         end if
         ms3=ms1
         ! carry on the integration to the source
         if(rs > rr) then
            call yint(ir,is,rr,rs,ms3, & 
                 sph_minors_solid_derivs_ng,scale3)        
         end if
         ! integrate down to the source
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_ng,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,is,r(nsl), & 
                 rs,ms2,sph_minors_solid_derivs_ng,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,is,r(nknot), & 
                 rs,ms2,sph_minors_solid_derivs_ng,scale2)
        end if
         ! form delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to 
         ! the receiver
         do i=1,ns
            scale4=0
            call b_solid_start_ngcow(ms2,ss(:,i),b)
            b=-b/delta
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_ng,scale4)     
            ! put the solution together
            call sph_ysol_ngcow(ms1,b,yout(:,i))
            ! unscale the solution
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      end if

      return
    end subroutine sph_solver_ng



    subroutine sph_solver_cow(ll_in,w,rr,rs,ir,is,ss,yout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rr,rs
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout

      integer(i4b) :: i,i1,scale1,scale2,scale3,scale4,iq,ns
      real(dp) :: rt
      complex(dpc) :: delta
      complex(dpc), dimension(4) :: b
      complex(dpc), dimension(2) :: mf1,mf2
      complex(dpc), dimension(5) :: ms1,ms2,ms3


      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if

      ! set parameters for the integration
      call set_l(ll_in)
      call set_ats(ats_switch)
      call set_omega(w)
      call sph_steps




      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then

         call sph_start_level(ll,real(omega),is,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= is) then
            i1 = is-1
            rt = r(i1)
         end if

         if(i1 <= nic) then 
            ! start in the inner core            
            call sph_minor_solid_start_ng(i1,rt,ms1)
            ! integrate from start to ICB
            call yint(i1,nic,rt,r(nic),ms1, & 
                 sph_minors_solid_derivs_cow,scale1)
            ! transform from solid to fluid systems
            call sph_minors_solid_2_fluid_ngcow(ms1,mf1)
            ! integrate from ICB to CMB
            call yint(nic+1,noc,r(nic),r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)
         else if (i1 > nic .and. i1 <=  noc) then    
            ! start in the outer core
            call sph_fluid_start_ng(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)              
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, &
            sph_minors_solid_derivs_cow,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            ! start in the mantle
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            ! integrate from start to source
            call yint(i1,is,rt,rs,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)        
         end if
         ms3=ms1
         ! carry on the integration to the receiver
         if(rr > rs) then
            call yint(is,ir,rs,rr,ms3, & 
                 sph_minors_solid_derivs_cow,scale3)        
         end if
         ! integrate down to the receiver
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_cow,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,ir,r(nsl), & 
                 rr,ms2,sph_minors_solid_derivs_cow,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,ir,r(nknot), & 
                 rr,ms2,sph_minors_solid_derivs_cow,scale2)
         end if
         ! calculate delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to the receiver
         do i=1,ns
            scale4=0         
            call b_solid_start_ngcow(ms1,ss(:,i),b)
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_cow,scale4)   
            ! put the solution together
            call sph_ysol_ngcow(ms2,b,yout(:,i))  
            ! unscale the solution
            yout(:,i)=yout(:,i)/delta
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      else
         ! integrate up to the receiver 
         call sph_start_level(ll,real(omega),ir,1,i1)
         rt = r(i1)
         if(i1 >= ir) then
            i1 = ir-1
            rt = r(i1)
         end if
         if(i1 <= nic) then
            call sph_minor_solid_start_ng(i1,rt,ms1)
            call yint(i1,nic,rt,r(nic),ms1, & 
                 sph_minors_solid_derivs_cow,scale1)
            call sph_minors_solid_2_fluid_ngcow(ms1,mf1)
            call yint(nic+1,noc,r(nic),r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)
         else if (i1 > nic .and. i1 <=  noc) then            
            call sph_fluid_start_ng(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)    
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_cow,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            call yint(i1,ir,rt,rr,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)                    
         end if
         ms3=ms1
         ! carry on the integration to the source
         if(rs > rr) then
            call yint(ir,is,rr,rs,ms3, & 
                 sph_minors_solid_derivs_cow,scale3)        
         end if
         ! integrate down to the source
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_cow,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,is,r(nsl), & 
                 rs,ms2,sph_minors_solid_derivs_cow,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,is,r(nknot), & 
                 rs,ms2,sph_minors_solid_derivs_cow,scale2)
        end if
         ! form delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to 
         ! the receiver
         do i=1,ns
            scale4=0
            call b_solid_start_ngcow(ms2,ss(:,i),b)
            b=-b/delta
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_cow,scale4)     
            ! put the solution together
            call sph_ysol_ngcow(ms1,b,yout(:,i))
            ! unscale the solution
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      end if

      return
    end subroutine sph_solver_cow







    subroutine sph_solver_fic(ll_in,w,rr,rs,ir,is,ss,yout)
      ! solves the spheriodal system on the assumption that the 
      ! inner core of the earth model is fluid
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rs,rr
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout      
      integer(i4b) :: i,i1,scale1,scale2,scale3,scale4,iq,ns
      real(dp) :: rt,rsi,rccon,rlcon,rfcon,rmu,rkappa, & 
           rqmu,rqkappa
      complex(dpc) :: delta
      complex(dpc), dimension(5) :: mf1,mf2
      complex(dpc), dimension(14) :: ms1,ms2,ms3,b



      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if

      ! set parameters for the integration
      call set_l(ll_in)
      call set_ats(ats_switch)
      call set_omega(w)
      call sph_steps




      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then

         call sph_start_level(ll,real(omega),is,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= is) then
            i1 = is-1
            rt = r(i1)
         end if
 
         if(i1 <= nic) then 
            ! start in the inner core
            call sph_minor_fluid_start(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)   
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, &
            sph_minors_solid_derivs,scale1) 
         else if (i1 > nic .and. i1 <=  noc) then    
            ! start in the outer core
            call sph_minor_fluid_start(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)              
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, &
            sph_minors_solid_derivs,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            ! start in the mantle
            call sph_minor_solid_start(i1,rt,ms1)            
            ! integrate from start to source
            call yint(i1,is,rt,rs,ms1, & 
                 sph_minors_solid_derivs,scale1)        
         end if
         ms3=ms1
         ! carry on the integration to the receiver
         if(rr > rs) then
            call yint(is,ir,rs,rr,ms3, & 
                 sph_minors_solid_derivs,scale3)        
         end if
         ! integrate down to the receiver
         if(ocean) then
            mf2(1)=one; mf2(2:5)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs,scale2)
            call sph_minors_fluid_2_solid(mf2,ms2)
            call yint(nsl,ir,r(nsl), & 
                 rr,ms2,sph_minors_solid_derivs,scale2)
         else
            ms2(1)=one; ms2(2:14)=zero
            call yint(nknot-1,ir,r(nknot), & 
                 rr,ms2,sph_minors_solid_derivs,scale2)
         end if
         ! calculate delta
         call sph_delta_solid(ms3,ms2,delta)
         ! integrate the b-vector from the source to the receiver
         do i=1,ns
            scale4=0         
            call b_solid_start(ms1,ss(:,i),b)
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs,scale4)   
            ! put the solution together
            call sph_ysol(ms2,b,yout(:,i))  
            ! unscale the solution
            yout(:,i)=-yout(:,i)/delta
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      else
         ! integrate up to the receiver 
         call sph_start_level(ll,real(omega),ir,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= ir) then
            i1 = ir-1
            rt = r(i1)
         end if
         if(i1 <= nic) then
            call sph_minor_solid_start(i1,rt,ms1)
            call yint(i1,nic,rt,r(nic),ms1, & 
                 sph_minors_solid_derivs,scale1)
            call sph_minors_solid_2_fluid(ms1,mf1)
            call yint(nic+1,noc,r(nic),r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)
            call sph_minors_fluid_2_solid(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs,scale1)
         else if (i1 > nic .and. i1 <=  noc) then            
            call sph_minor_fluid_start(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs,scale1)    
            call sph_minors_fluid_2_solid(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            call sph_minor_solid_start(i1,rt,ms1)            
            call yint(i1,ir,rt,rr,ms1, & 
                 sph_minors_solid_derivs,scale1)                    
         end if
         ms3=ms1
         ! carry on the integration to the source
         if(rs > rr) then
            call yint(ir,is,rr,rs,ms3, & 
                 sph_minors_solid_derivs,scale3)        
         end if
         ! integrate down to the source
         if(ocean) then
            mf2(1)=one; mf2(2:5)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs,scale2)
            call sph_minors_fluid_2_solid(mf2,ms2)
            call yint(nsl,is,r(nsl), & 
                 rs,ms2,sph_minors_solid_derivs,scale2)
         else
            ms2(1)=one; ms2(2:14)=zero
            call yint(nknot,is,r(nknot), & 
                 rs,ms2,sph_minors_solid_derivs,scale2)
        end if
         ! form delta
         call sph_delta_solid(ms3,ms2,delta)
         ! integrate the b-vector from the source to 
         ! the receiver
         do i=1,ns
            scale4=0
            call b_solid_start(ms2,ss(:,i),b)
            b=-b/delta
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs,scale4)     
            ! put the solution together
            call sph_ysol(ms1,b,yout(:,i))
            ! unscale the solution
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      end if
      
      return
    end subroutine sph_solver_fic



    subroutine sph_solver_ng_fic(ll_in,w,rr,rs,ir,is,ss,yout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rr,rs
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout

      integer(i4b) :: i,i1,scale1,scale2,scale3,scale4,iq,ns
      real(dp) :: rt
      complex(dpc) :: delta
      complex(dpc), dimension(4) :: b
      complex(dpc), dimension(2) :: mf1,mf2
      complex(dpc), dimension(5) :: ms1,ms2,ms3


      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if

      ! set parameters for the integration
      call set_l(ll_in)
      call set_ats(ats_switch)
      call set_omega(w)
      call sph_steps

      
      


      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then

         call sph_start_level(ll,real(omega),is,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= is) then
            i1 = is-1
            rt = r(i1)
         end if


         if(i1 <= nic) then 
            ! start in the inner core            
            call sph_fluid_start_ng(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)
         else if (i1 > nic .and. i1 <=  noc) then    
            ! start in the outer core
            call sph_fluid_start_ng(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)              
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, &
            sph_minors_solid_derivs_ng,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            ! start in the mantle
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            ! integrate from start to source
            call yint(i1,is,rt,rs,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)        
         end if
         ms3=ms1
         ! carry on the integration to the receiver
         if(rr > rs) then
            call yint(is,ir,rs,rr,ms3, & 
                 sph_minors_solid_derivs_ng,scale3)        
         end if
         ! integrate down to the receiver
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_ng,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,ir,r(nsl), & 
                 rr,ms2,sph_minors_solid_derivs_ng,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,ir,r(nknot), & 
                 rr,ms2,sph_minors_solid_derivs_ng,scale2)
         end if
         ! calculate delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to the receiver
         do i=1,ns
            scale4=0         
            call b_solid_start_ngcow(ms1,ss(:,i),b)
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_ng,scale4)   
            ! put the solution together
            call sph_ysol_ngcow(ms2,b,yout(:,i))  
            ! unscale the solution
            yout(:,i)=yout(:,i)/delta
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      else
         ! integrate up to the receiver 
         call sph_start_level(ll,real(omega),ir,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= ir) then
            i1 = ir-1
            rt = r(i1)
         end if
         if(i1 <= nic) then
            call sph_fluid_start_ng(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)
         else if (i1 > nic .and. i1 <=  noc) then            
            call sph_fluid_start_ng(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_ng,scale1)    
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_ng,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            call yint(i1,ir,rt,rr,ms1, & 
                 sph_minors_solid_derivs_ng,scale1)                    
         end if
         ms3=ms1
         ! carry on the integration to the source
         if(rs > rr) then
            call yint(ir,is,rr,rs,ms3, & 
                 sph_minors_solid_derivs_ng,scale3)        
         end if
         ! integrate down to the source
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_ng,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,is,r(nsl), & 
                 rs,ms2,sph_minors_solid_derivs_ng,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,is,r(nknot), & 
                 rs,ms2,sph_minors_solid_derivs_ng,scale2)
        end if
         ! form delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to 
         ! the receiver
         do i=1,ns
            scale4=0
            call b_solid_start_ngcow(ms2,ss(:,i),b)
            b=-b/delta
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_ng,scale4)     
            ! put the solution together
            call sph_ysol_ngcow(ms1,b,yout(:,i))
            ! unscale the solution
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      end if

      return
    end subroutine sph_solver_ng_fic





    subroutine sph_solver_cow_fic(ll_in,w,rr,rs,ir,is,ss,yout)
      use nrtype; use module_model
      use module_int; use module_spline
      use module_start
      implicit none
      integer(i4b), intent(in) :: ll_in
      complex(dpc), intent(in) :: w
      real(dp), intent(in) :: rr,rs
      integer(i4b), intent(in) :: ir,is
      complex(dpc), dimension(:,:), intent(in) :: ss
      complex(dpc), dimension(:,:), intent(out) :: yout

      integer(i4b) :: i,i1,scale1,scale2,scale3,scale4,iq,ns
      real(dp) :: rt
      complex(dpc) :: delta
      complex(dpc), dimension(4) :: b
      complex(dpc), dimension(2) :: mf1,mf2
      complex(dpc), dimension(5) :: ms1,ms2,ms3


      ! number of sources
      ns = size(ss,2)

      if(ns /= size(yout,2)) then
         stop 'yout too small'
      end if

      ! set parameters for the integration
      call set_l(ll_in)
      call set_ats(ats_switch)
      call set_omega(w)
      call sph_steps




      scale1=0; scale2=0; scale3=0
      if(rr >= rs) then
         call sph_start_level(ll,real(omega),is,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= is) then
            i1 = is-1
            rt = r(i1)
         end if
         if(i1 <= nic) then 
            ! start in the inner core            
            call sph_fluid_start_ng(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)
         else if (i1 > nic .and. i1 <=  noc) then    
            ! start in the outer core
            call sph_fluid_start_ng(i1,rt,mf1)
            ! integrate from start to CMB
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)              
            ! transform from fluid to solid systems
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            ! integrate from CMB to the source
            call yint(noc+1,is,r(noc+1),rs,ms1, &
            sph_minors_solid_derivs_cow,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            ! start in the mantle
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            ! integrate from start to source
            call yint(i1,is,rt,rs,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)        
         end if
         ms3=ms1
         ! carry on the integration to the receiver
         if(rr > rs) then
            call yint(is,ir,rs,rr,ms3, & 
                 sph_minors_solid_derivs_cow,scale3)        
         end if
         ! integrate down to the receiver
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_cow,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,ir,r(nsl), & 
                 rr,ms2,sph_minors_solid_derivs_cow,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,ir,r(nknot), & 
                 rr,ms2,sph_minors_solid_derivs_cow,scale2)
         end if
         ! calculate delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to the receiver
         do i=1,ns
            scale4=0         
            call b_solid_start_ngcow(ms1,ss(:,i),b)
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_cow,scale4)   
            ! put the solution together
            call sph_ysol_ngcow(ms2,b,yout(:,i))  
            ! unscale the solution
            yout(:,i)=yout(:,i)/delta
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      else
         ! integrate up to the receiver 
         call sph_start_level(ll,real(omega),ir,1,i1)
         rt = r(i1)
         if(i1 == 1) rt = min(0.1_dp*dr,0.1_dp*r(2))
         if(i1 >= ir) then
            i1 = ir-1
            rt = r(i1)
         end if
         if(i1 <= nic) then
            call sph_fluid_start_ng(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)
         else if (i1 > nic .and. i1 <=  noc) then            
            call sph_fluid_start_ng(i1,rt,mf1)
            call yint(i1,noc,rt,r(noc),mf1, & 
                 sph_minors_fluid_derivs_cow,scale1)    
            call sph_minors_fluid_2_solid_ngcow(mf1,ms1)
            call yint(noc+1,ir,r(noc+1),rr,ms1, & 
                 sph_minors_solid_derivs_cow,scale1) 
         else if(i1 > noc .and. i1 <= nsl) then
            call sph_minor_solid_start_ng(i1,rt,ms1)            
            call yint(i1,ir,rt,rr,ms1, & 
                 sph_minors_solid_derivs_cow,scale1)                    
         end if
         ms3=ms1
         ! carry on the integration to the source
         if(rs > rr) then
            call yint(ir,is,rr,rs,ms3, & 
                 sph_minors_solid_derivs_cow,scale3)        
         end if
         ! integrate down to the source
         if(ocean) then
            mf2(1)=one; mf2(2)=zero
            call yint(nknot-1,nsl+1,r(nknot), & 
                 r(nsl+1),mf2,sph_minors_fluid_derivs_cow,scale2)
            call sph_minors_fluid_2_solid_ngcow(mf2,ms2)
            call yint(nsl,is,r(nsl), & 
                 rs,ms2,sph_minors_solid_derivs_cow,scale2)
         else
            ms2(1)=one; ms2(2:5)=zero
            call yint(nknot,is,r(nknot), & 
                 rs,ms2,sph_minors_solid_derivs_cow,scale2)
        end if
         ! form delta
         call sph_delta_solid_ngcow(ms3,ms2,delta)
         ! integrate the b-vector from the source to 
         ! the receiver
         do i=1,ns
            scale4=0
            call b_solid_start_ngcow(ms2,ss(:,i),b)
            b=-b/delta
            call yint(is,ir,rs,rr,b, & 
                 sph_bvec_solid_derivs_cow,scale4)     
            ! put the solution together
            call sph_ysol_ngcow(ms1,b,yout(:,i))
            ! unscale the solution
            yout(:,i)=yout(:,i)*huge**(scale3-scale4)
         end do
      end if

      return
    end subroutine sph_solver_cow_fic





    
    subroutine source_vector_sph(l,w,is,rs,mm,svs)
      !---------------------------------------------------------!
      ! returns source discontinuity vector for a given moment  !
      ! tensor source. The moment tensor is input as:           !
      ! m(1) = M_{r,r}                                          !
      ! m(2) = M_{r,theta}                                      !
      ! m(3) = M_{r,phi}                                        !
      ! m(4) = M_{theta,theta}                                  !
      ! m(5) = M_{theta,phi}                                    !
      ! m(6) = M_{phi,phi}                                      !
      !---------------------------------------------------------!
      use nrtype
      use module_model
      use module_spline
      use module_int
      implicit none
      ! input/output
      integer(i4b), intent(in) :: l
      complex(dpc), intent(in) :: w
      integer(i4b), intent(in) :: is
      real(dp), intent(in) :: rs
      real(dp), dimension(6), intent(in) :: mm
      complex(dpc), dimension(5,6), intent(out) :: svs
      ! local variables
      real(dp) :: rlcon,rccon,rfcon,rmu,rqmu,rkappa, & 
           rqkappa,rxlam,rxa2,qaa,qff,qll
      complex(dpc) :: clcon,cccon,cfcon 


      ! set some paramters
      call set_l(l)
      call set_ats(ats_switch)
      call set_omega(w)


      ! compute some parameters
      rlcon    = splint_dis_r(r,lcon,lcon_cs,rs,is)
      rccon    = splint_dis_r(r,ccon,ccon_cs,rs,is)
      rfcon    = splint_dis_r(r,fcon,fcon_cs,rs,is)
      rmu      = splint_dis_r(r,mu,mu_cs,rs,is)
      rkappa   = splint_dis_r(r,kappa,kappa_cs,rs,is)
      rqmu     = splint_dis_r(r,qmu,qmu_cs,rs,is)     
      rqkappa  = splint_dis_r(r,qkappa,qkappa_cs,rs,is)           
      rxlam    = splint_dis_r(r,xlam,xlam_cs,rs,is)           
      rxa2    = splint_dis_r(r,xa2,xa2_cs,rs,is)           

      qaa = 1.0_dp+rxa2*lno
      qff = 1.0_dp+rxlam*lno
      qll = 1.0_dp+rqmu*lno

      cccon = rccon*qaa
      cccon = 1.0_dp/cccon
      cfcon = rfcon*qff
      clcon = rlcon*qll


      ! initialize the source vector
      svs = 0.0_dp


      if(l == 1) then

         ! m = -1 term
         svs(2,2) = 0.5_dp*rs*(mm(2)+ii*mm(3))/clcon

         ! m = 0 term
         svs(3,1) = rs*mm(1)*cccon 
         svs(3,4) = -mm(4)-mm(6)+2.0_dp*mm(1)*cfcon*cccon
         svs(3,5) = 0.5_dp*zeta*(mm(4)+mm(6))-zeta*mm(1)*cfcon*cccon
         
         ! m = 1 term
         svs(4,2) = 0.5_dp*rs*(-mm(2)+ii*mm(3))/clcon

      else 


         ! m = -2 term
         svs(1,5) = 0.25_dp*szeta2m2*(-mm(4)+mm(6)-2.0_dp*ii*mm(5))

         ! m = -1 term
         svs(2,2) = 0.5_dp*rs*(mm(2)+ii*mm(3))/clcon

         ! m = 0 term
         svs(3,1) = rs*mm(1)*cccon
         svs(3,4) = -mm(4)-mm(6)+2.0_dp*mm(1)*cfcon*cccon
         svs(3,5) = 0.5_dp*zeta*(mm(4)+mm(6))-zeta*mm(1)*cfcon*cccon

         ! m = 1 term
         svs(4,2) = 0.5_dp*rs*(-mm(2)+ii*mm(3))/clcon

         ! m = 2 term
         svs(5,5) = 0.25_dp*szeta2m2*(-mm(4)+mm(6)+2.0_dp*ii*mm(5))

      end if

      
      svs = nu*svs/(ii*w*rs*rs)

 
      return
    end subroutine source_vector_sph



    !===========================================!
    !             utility rountines             !
    !===========================================!

    subroutine set_ats_switch(ats_switch_in)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: ats_switch_in      
      ats_switch = real(ats_switch_in)      
      return
    end subroutine set_ats_switch

end module module_solver
