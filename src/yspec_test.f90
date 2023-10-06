program yspec_test


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
  integer(i4b), parameter :: io1=7,io2=8                           !
  ! model file name:                                               !
  character(len=256) :: model_file                                 !
  ! parameter file                                                 !
  character(len=256) :: param_file                                 !
  !----------------------------------------------------------------!


  integer(i4b) :: ios,l,nw,iw,i1,i2,nt,mex,qex, & 
       is,ir

  real(dp) :: f1,f2,w1,w2,dw,ep,rs,rr,dt,tout,df, & 
       wmin,wmax,zeta
  real(dp), dimension(6) :: mm

  complex(dpc) :: w
  complex(dpc), dimension(2) :: sr0
  complex(dpc), dimension(5,2) :: svt                              
  complex(dpc), dimension(5,6) :: svs                              
  complex(dpc), dimension(2,2) :: s2tmp,ytor,yrad,prad             
  complex(dpc), dimension(6,4) :: s6tmp,ysph                       

  !--------------------------------------------------!
  !      set up the normalization parameters         !
  !--------------------------------------------------!
  call set_parameters



  ! get the model file
  call get_string(' model file = ',model_file)
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




  ! set the attenuation to be on (arg == 1) or off (arg == 0)
  call set_ats_switch(1) 


  zeta  = sqrt(real(l*(l+1)))


  open(97,file='yspec_test.rad')
  open(98,file='yspec_test.tor')
  open(99,file='yspec_test.sph')

     ! start loop over omega
     wloop: do iw = 1,nw       

        print *, ' frequencies to go  = ',nw-iw+1
 
        w=wmin+(iw-1)*dw-ii*ep

        yrad = 0.0_dp
        ytor = 0.0_dp
        ysph = 0.0_dp

        if(l == 0) then

           !-----------------------------------------!
           !               radial modes              !
           !-----------------------------------------!
           
           !get the radial source vector
           call source_vector_rad(w,is,rs,mm,sr0)

           ! solve the radial equations
           call rad_solver_wg(w,rr,rs,ir,is,s2tmp,yrad,prad)

        else               

           !--------------------------------!
           !        toroidal modes          !
           !--------------------------------!

           if(l > 1) then
           
              ! get the source vectors
              call source_vector_tor(l,w,is,rs,mm,svt)
              
              ! solve the toroidal equations
              call tor_solver(l,w,rr,rs,ir,is,s2tmp,ytor)
              

           end if

           !------------------------------!
           !        spheroidal modes      !
           !------------------------------!
           
           ! get the spheroidal source vectors
           call source_vector_sph(l,w,is,rs,mm,svs)

           ! solve the spheroidal equations
           if(fic) then
              ! solver with fluid inner core
              call sph_solver_fic(l,w,rr,rs,ir,is,s6tmp,ysph)
           else
              ! solver with solid inner core
              call sph_solver(l,w,rr,rs,ir,is,s6tmp,ysph)
           end if
           
        end if
        
        write(97,*) real(w)*fre_norm*500.0_dp/pi, & 
                            real(yrad(1,1)),imag(yrad(1,1)), & 
                            real(yrad(1,2)),imag(yrad(1,2))
        
        write(98,*) real(w)*fre_norm*500.0_dp/pi, & 
                            real(ytor(1,1)),imag(ytor(1,1)), & 
                            real(ytor(1,2)),imag(ytor(1,2))

        write(99,*) real(w)*fre_norm*500.0_dp/pi, & 
                            real(ysph(1,1)),imag(ysph(1,1)), & 
                            real(ysph(2,1)),imag(ysph(2,1)), &
                            real(ysph(1,2)),imag(ysph(1,2)), & 
                            real(ysph(2,2)),imag(ysph(2,2)), &
                            real(ysph(1,3)),imag(ysph(1,3)), & 
                            real(ysph(2,3)),imag(ysph(2,3)), &
                            real(ysph(1,4)),imag(ysph(1,4)), & 
                            real(ysph(2,4)),imag(ysph(2,4))

        
     end do wloop

     close(97)
     close(98)
     close(99)


end program yspec_test
