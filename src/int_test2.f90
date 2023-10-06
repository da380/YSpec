program int_test2


  use nrtype                                                       
  use nrutil                                                       
  use module_util                                                  
  use module_model                                                 
  use module_solver                                                
  use module_fourier

  implicit none


  !================================================================!
  ! input/output files used:                                       !
  !================================================================!
  integer(i4b), parameter :: io1=7,io2=8,io3=9                     !
  ! model file name:                                               !
  character(len=256) :: model_file                                 !
  ! parameter file                                                 !
  character(len=256) :: param_file                                 !
  !----------------------------------------------------------------!

  logical(lgt) :: ltmp

  integer(i4b) :: ios,l,nw,iw,i1,i2,nt,mex,qex, & 
       is,ir

  real(dp) :: f1,f2,w1,w2,dw,ep,rs,rr,dt,tout,df, & 
       wmin,wmax,zeta,wr
  real(dp), dimension(6) :: mm

  complex(dpc) :: w
  complex(dpc), dimension(2) :: sr0
  complex(dpc), dimension(5,2) :: svt                              
  complex(dpc), dimension(5,6) :: svs                              
  complex(dpc), dimension(2,2) :: s2tmp,ytor,yrad,prad             
  complex(dpc), dimension(6,4) :: s6tmp,ysph,ysph2,ysph3

  !--------------------------------------------------!
  !      set up the normalization parameters         !
  !--------------------------------------------------!
  call set_parameters



  ! get the model file
  call get_string(' model file = ',model_file)
  inquire(file=trim(model_file),exist = ltmp)
  if(.not.ltmp) stop ' can''t find model file!'
  open(io2,file=model_file,form='formatted', & 
       action='read',iostat=ios)
  if(ios /= 0) stop ' problem opening model file'
  call  read_model(io2)
  close(io2)


  ! get value of l
  call get_integer(' l = ',l)
  
  ! get frequency range
  call get_float(' f1 (mHz) = ',f1)
  call get_float(' f2 (mHz) = ',f2)

  ! get the time step 
  call get_float(' dt (sec) = ',dt)


  ! get the time length 
  call get_float(' tout (min) = ',tout)

  ! set the moment tensor
  mm = 0.0_dp
  mm(1) = 1.0_dp
  mm(2) = 2.0_dp
  mm(3) = 3.0_dp

  ! get source depth
  call get_float(' source depth (km) = ',rs)
  
  ! get receiver depth 
  call get_float(' receiver depth (km) = ',rr)
  

  ! normalize the inputs
  f1 = f1*t_norm/1000.0_dp
  f2 = f2*t_norm/1000.0_dp
  dt = dt/t_norm
  rs=1.0_dp-rs*1000.0_dp/r_norm
  rr=1.0_dp-rr*1000.0_dp/r_norm
  tout=tout*60.0_dp/t_norm

  ! set parameters for frequency spacing
  mex = 5
  qex = 1

  ! get the frequency spacing
  call fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
  

  ! set some more parameters
  dw = twopi_d*df
  if(i1 == i2) stop 'too few frequency steps'  
  wmin=(i1-1)*dw
  wmax=(i2-1)*dw
  nw=i2-i1+1
  if(nw == 1) stop 'to few frequencies'

  
  !--------------------------------------------------!
  ! find layer numbers for the source and receivers  !
  !--------------------------------------------------!
  is=find_radius_index(rs,-1)
  if(is == 0)  stop ' source not in model!'
  if(is > nsl) stop ' source in the ocean!'
  ir=find_radius_index(rr,-1)
  if(ir == 0)  stop ' receiver not in model!'
  if(ir > nsl) stop ' receiver in the ocean!'



  ! set up dummy source vectors
  s2tmp(:,:) = 0.0_dp
  s2tmp(1,1) = 1.0_dp
  s2tmp(2,2) = 1.0_dp

  s6tmp(:,:) = 0.0_dp
  s6tmp(1,1)   = 1.0_dp
  s6tmp(2,2)   = 1.0_dp
  s6tmp(4,3)   = 1.0_dp
  s6tmp(5,4)   = 1.0_dp


  ! open the output file
  open(io1,file='int_test2.out1')
  open(io2,file='int_test2.out2')
  open(io3,file='int_test2.out3')
  

  ! start loop over omega
  wloop: do iw = 1,nw       

     print *, ' frequencies to go  = ',nw-iw+1
 
     w=wmin+(iw-1)*dw-ii*ep

                
     ! get the spheroidal source vectors
     call source_vector_sph(l,w,is,rs,mm,svs)
     
     ! solve the spheroidal equations
     
     call sph_solver_minors(l,w,rr,rs,ir,is,s6tmp,ysph)
     call sph_solver_simple(l,w,rr,rs,ir,is,s6tmp,ysph2)
     call sph_solver_symp(l,w,rr,rs,ir,is,s6tmp,ysph3)

     
     wr = real(w)*fre_norm*500.0_dp/pi_d
     write(io1,'(17e25.10e2)') wr,real(ysph(1,1)),imag(ysph(1,1)), & 
                                  real(ysph(2,1)),imag(ysph(2,1)), & 
                                  real(ysph(1,2)),imag(ysph(1,2)), & 
                                  real(ysph(2,2)),imag(ysph(2,2)), &
                                  real(ysph(1,3)),imag(ysph(1,3)), & 
                                  real(ysph(2,3)),imag(ysph(2,3)), &
                                  real(ysph(1,4)),imag(ysph(1,4)), & 
                                  real(ysph(2,4)),imag(ysph(2,4))
     
     write(io2,'(17e25.10e2)') wr,real(ysph2(1,1)),imag(ysph2(1,1)), & 
                                  real(ysph2(2,1)),imag(ysph2(2,1)), & 
                                  real(ysph2(1,2)),imag(ysph2(1,2)), & 
                                  real(ysph2(2,2)),imag(ysph2(2,2)), &
                                  real(ysph2(1,3)),imag(ysph2(1,3)), & 
                                  real(ysph2(2,3)),imag(ysph2(2,3)), &
                                  real(ysph2(1,4)),imag(ysph2(1,4)), & 
                                  real(ysph2(2,4)),imag(ysph2(2,4))

                    
     write(io3,'(17e25.10e2)') wr,real(ysph3(1,1)),imag(ysph3(1,1)), & 
                                  real(ysph3(2,1)),imag(ysph3(2,1)), & 
                                  real(ysph3(1,2)),imag(ysph3(1,2)), & 
                                  real(ysph3(2,2)),imag(ysph3(2,2)), &
                                  real(ysph3(1,3)),imag(ysph3(1,3)), & 
                                  real(ysph3(2,3)),imag(ysph3(2,3)), &
                                  real(ysph3(1,4)),imag(ysph3(1,4)), & 
                                  real(ysph3(2,4)),imag(ysph3(2,4))
        
     end do wloop

     close(io1)
     close(io2)
     close(io3)


contains 




  subroutine sph_solver_minors(ll_in,w,rr,rs,ir,is,ss,yout)
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
    call set_ats(0.0_dp)
    call set_omega(w)
    call sph_steps



    
    scale1=0; scale2=0; scale3=0
    call sph_start_level(ll,real(omega),is,1,i1)
    if(i1 <= noc) i1 = noc+1 ! force start in the mantle
    i1 = noc+1
    rt = r(i1)    
    if(i1 >= is) then
       i1 = is-1
       rt = r(i1)
    end if

    
    ! start in the mantle
    if(i1 == noc+1) then
       ms1 = 0.0_dp
       ms1(1) = 1.0_dp
    else
       call sph_minor_solid_start(i1,rt,ms1)   
    end if
    ! integrate from start to source
    call yint(i1,is,rt,rs,ms1, & 
         sph_minors_solid_derivs,scale1)        
    ms3=ms1
    ! carry on the integration to the receiver
    if(rr > rs) then
       call yint(is,ir,rs,rr,ms3, & 
            sph_minors_solid_derivs,scale3)        
    end if
    ! integrate down to the receiver
    ms2 = 0.0_dp
    ms2(1) = 1.0_dp
    call yint(nsl-1,ir,r(nsl), & 
         rr,ms2,sph_minors_solid_derivs,scale2)
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


    
    return
  end subroutine sph_solver_minors



  subroutine sph_solver_simple(ll_in,w,rr,rs,ir,is,ss,yout)
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
    integer(i4b) :: i,i1,scale3,scale4,iq,ns,k
    integer(i4b), dimension(3) :: scale1,scale2
    integer(i4b), dimension(6) :: indx

    real(dp) :: rt,rsi,rccon,rlcon,rfcon,rmu,rkappa, & 
         rqmu,rqkappa,d
    complex(dpc) :: delta
    
    complex(dpc), dimension(6) :: y3
    complex(dpc), dimension(6,3) :: y1,y2
    complex(dpc), dimension(6,6) :: yy

    


    ! number of sources
    ns = size(ss,2)
    
    if(ns /= size(yout,2)) then
       stop 'yout too small'
    end if
    
    ! set parameters for the integration
    call set_l(ll_in)
    call set_ats(0.0_dp)
    call set_omega(w)
    call sph_steps



    
    scale1=0; scale2=0; scale3=0
    call sph_start_level(ll,real(omega),is,1,i1)
    if(i1 <= noc) i1 = noc+1 ! force start in the mantle
    i1 = noc+1
    rt = r(i1)    
    if(i1 >= is) then
       i1 = is-1
       rt = r(i1)
    end if


    ! start in the mantle
    if(i1 == noc+1) then
       y1 = 0.0_dp
       do k = 1,3
          y1(k,k) = 1.0_dp
       end do
    else
       call sph_solid_start(i1,rt,y1)   
    end if    

    ! integrate three solutions upto source depth 
    do k = 1,3
       call yint(i1,is,rt,rs,y1(:,k), & 
            sph_solid_derivs,scale1(k))            
    end do


    ! integrate three solutions down to source depth
    y2 = 0.0_dp
    do k = 1,3
       y2(k,k) = 1.0_dp 
       call yint(nsl-1,is,r(nsl),rs,y2(:,k), & 
            sph_solid_derivs,scale2(k)) 
    end do


    ! put together the Y matrix
    yy(:,1:3) = y1(:,:)
    yy(:,4:6) = y2(:,:)
    

    ! lu decompose the Y matrix
    call ludcmp(yy,indx,d)

    delta = yy(1,1)
    do k = 2,6
       delta = delta*yy(k,k)
    end do
    delta = d*delta


    do k = 1,ns

       
       ! solve the linear system at the source
       y3 = ss(:,k)
       call lubksb(yy,indx,y3)   

       ! get the solution at the source
       y3 = matmul(y2,y3(4:6))

       ! propagate solution upto the receiver
       call  yint(is,ir,rs,rr,y3, & 
            sph_solid_derivs,scale3)


       ! unscale the solution 
       yout(:,k) = y3*huge**scale3
    

    end do



    
    return
  end subroutine sph_solver_simple



  subroutine sph_solver_symp(ll_in,w,rr,rs,ir,is,ss,yout)
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
    integer(i4b) :: i,i1,scale3,scale4,iq,ns,k
    integer(i4b), dimension(3) :: scale1,scale2
    integer(i4b), dimension(6) :: indx

    real(dp) :: rt,rsi,rccon,rlcon,rfcon,rmu,rkappa, & 
         rqmu,rqkappa,d
    complex(dpc) :: delta

    complex(dpc), dimension(3) :: st
    complex(dpc), dimension(6) :: y3
    complex(dpc), dimension(6,3) :: y1,y2,ys1,ys2
    complex(dpc), dimension(3,3) :: y1ty2,adj

    


    ! number of sources
    ns = size(ss,2)
    
    if(ns /= size(yout,2)) then
       stop 'yout too small'
    end if
    
    ! set parameters for the integration
    call set_l(ll_in)
    call set_ats(0.0_dp)
    call set_omega(w)
    call sph_steps



    
    scale1=0; scale2=0; scale3=0
    call sph_start_level(ll,real(omega),is,1,i1)
    if(i1 <= noc) i1 = noc+1 ! force start in the mantle
    i1 = noc+1
    rt = r(i1)    
    if(i1 >= is) then
       i1 = is-1
       rt = r(i1)
    end if


    ! start in the mantle
    if(i1 == noc+1) then
       y1 = 0.0_dp
       do k = 1,3
          y1(k,k) = 1.0_dp
       end do
    else
       call sph_solid_start(i1,rt,y1)   
    end if    

    ! integrate three solutions upto source depth 
    do k = 1,3
       call yint(i1,is,rt,rs,y1(:,k), & 
            sph_solid_derivs,scale1(k))            
    end do


    ! integrate three solutions down to source depth
    y2 = 0.0_dp
    do k = 1,3
       y2(k,k) = 1.0_dp 
       call yint(nsl-1,is,r(nsl),rs,y2(:,k), & 
            sph_solid_derivs,scale2(k)) 
    end do


    ! form the sigma matrices
    ys1(1:3,:) =  y1(4:6,:)
    ys1(4:6,:) = -y1(1:3,:)
    ys2(1:3,:) =  y2(4:6,:)
    ys2(4:6,:) = -y2(1:3,:)

    ! form the matrix \tilde{Y}^{(1)T}Y^{(2)}
    y1ty2 = matmul(transpose((ys1)),y2)
    
    
    ! compute the determinant of the system
    delta = det3(y1ty2(:,1),y1ty2(:,2),y1ty2(:,3))


    ! compute the elements of the adjoint matrix
     adj(1,1)=y1ty2(2,2)*y1ty2(3,3)-y1ty2(2,3)*y1ty2(3,2)
     adj(1,2)=y1ty2(3,2)*y1ty2(1,3)-y1ty2(1,2)*y1ty2(3,3)
     adj(1,3)=y1ty2(1,2)*y1ty2(2,3)-y1ty2(1,3)*y1ty2(2,2)
     adj(2,1)=y1ty2(2,3)*y1ty2(3,1)-y1ty2(2,1)*y1ty2(3,3)
     adj(2,2)=y1ty2(1,1)*y1ty2(3,3)-y1ty2(1,3)*y1ty2(3,1)
     adj(2,3)=y1ty2(1,3)*y1ty2(2,1)-y1ty2(1,1)*y1ty2(2,3)
     adj(3,1)=y1ty2(2,1)*y1ty2(3,2)-y1ty2(2,2)*y1ty2(3,1)
     adj(3,2)=y1ty2(1,2)*y1ty2(3,1)-y1ty2(1,1)*y1ty2(3,2)
     adj(3,3)=y1ty2(1,1)*y1ty2(2,2)-y1ty2(1,2)*y1ty2(2,1)


     
     ! put the solutions together
     do k = 1,ns
     
        st = matmul(transpose((ys1)),ss(:,k))
        st = matmul(adj,st)
        y3 = -matmul(y2,st)/delta


        ! integrate the solution to the receiver
        call  yint(is,ir,rs,rr,y3,sph_solid_derivs,scale3)            

        yout(:,k) = y3*huge**(scale3)


     end do



    
    return
  end subroutine sph_solver_symp










  
    subroutine ludcmp(a,indx,d)
      !---------------------------------------------------------------!
      ! this routine peforms an LU decomposition on the input matrix  !
      ! a using partial pivoting                                      !
      !---------------------------------------------------------------!
      use nrtype; use nrutil, only : assert_eq,imaxloc,nrerror, &     !
           outerprod,swap                                             !
      !---------------------------------------------------------------!
      implicit none                              
      complex(dpc), dimension(:,:), intent(inout) :: a                
      integer(I4B), dimension(:), intent(out) :: indx                 
      real(dp), intent(out) :: d                                      
      real(dp), dimension(size(a,1)) :: vv                            
      real(dP), parameter :: TINY=1.0e-20_sp                          
      integer(I4B) :: j,n,imax                                        
      n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')            
      d=1.0_sp                                                        
      vv=maxval(abs(a),dim=2)                                         
      if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')   
      vv=1.0_sp/vv                                                    
      do j=1,n                                                        
         imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))                   
         if (j /= imax) then                                         
            call swap(a(imax,:),a(j,:))                              
            d=-d                                                     
            vv(imax)=vv(j)                                           
         end if                                                      
         indx(j)=imax                                                
         if (abs(a(j,j)) == 0.0_sp) a(j,j)=TINY                      
         a(j+1:n,j)=a(j+1:n,j)/a(j,j)                                
         a(j+1:n,j+1:n)=a(j+1:n,j+1:n)- &                            
              outerprod(a(j+1:n,j),a(j,j+1:n))                       
      end do                                                         
    end subroutine ludcmp





    subroutine lubksb(a,indx,b)
      !---------------------------------------------------------------!
      ! this routine solves the linear system a*x=b for x where the   !
      ! matrix a has been LU decomposed into upper and lower triang-  !
      ! ular matrices using the routine ludcmp                        !
      !---------------------------------------------------------------!
      use nrtype; use nrutil, only : assert_eq                        !
      !---------------------------------------------------------------!
      implicit none    
      complex(dpc), dimension(:,:), intent(in) :: a                   
      integer(I4B), dimension(:), intent(in) :: indx                  
      complex(dpc), dimension(:), intent(inout) :: b                  
      integer(I4B) :: i,n,jj,ll                                       
      complex(dpc) :: summ                                            
      n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')            
      jj=0                                                            
      do i=1,n                                                        
         ll=indx(i)                                                   
         summ=b(ll)                                                   
         b(ll)=b(i)                                                   
         if (jj /= 0) then                                            
            summ=summ-dot_product(a(i,jj:i-1),b(jj:i-1))              
         else if (summ /= 0.0) then                                   
            jj=i                                                      
         end if                                                       
         b(i)=summ                                                    
      end do                                                          
      do i=n,1,-1                                                     
         b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)        
      end do                                                          
    end subroutine lubksb                                            



end program int_test2
