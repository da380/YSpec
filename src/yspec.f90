program yspec
  !================================================================!
  ! This program computes the time domain solution to the          !
  ! seismic wave equation in a spherically symmetric earth model   !
  ! using the direct radial integration method described by        !
  ! Al-Attar & Woodhouse (2008).                                   !
  !================================================================!

  !================================================================!
  ! Modules used:                                                  !
  !================================================================!
  use nrtype                                                       !
  use nrutil                                                       !
  use module_util                                                  !
  use module_model                                                 !
  use module_solver                                                !
  use module_function                                              !
  use module_fourier                                               !
  !----------------------------------------------------------------!
  implicit none                                                    !
  !----------------------------------------------------------------!

  !================================================================!
  ! input/output files used:                                       !
  !================================================================!
  integer(i4b), parameter :: io1=7,io2=8,io3=9                     !
  ! model file name:                                               !
  character(len=256) :: model_file                                 !
  ! parameter file                                                 !
  character(len=256) :: param_file                                 !
  !----------------------------------------------------------------!

  
  !================================================================!
  ! variables local to the main program:                           !
  !================================================================!  
  character(len=10) :: string                                      !
  character(len=256) :: pref_out                                   !
  character(len=256) :: seis_out                                   !
  character(len=256) :: phi_coef_pref ='phi_coef.out'              !
  character(len=256) :: phi_coef_out                               !
  logical(lgt) :: ltmp                                             !  
  integer(i4b) :: i,j,k,l,lmin,lmax, &                             !
       nw,ios,i1,i2,is,ir,nt,nr,nout,mex, &                        !
       qex,m,im,isv,msign,ats_in,grav_switch, &                    !
       out_switch,cor_switch,phi_coef_switch                       !
  real(dp) :: wmin,wmax,dw,rr,rs,ep,lons,lats, &                   !
       zeta,wt,dt,df,t,f11,f12,f21,f22,f,tout,depth_rec, &         !
       depth_source,sqzm2,sqzm6,sqz,f1,f2,x,xp,xc                  !
  real(dp), dimension(6) :: mm                                     !
  real(dp), dimension(:), allocatable :: latr,lonr,delta, &        !
       azst,csc,azep                                               !
  real(dp), dimension(:,:), allocatable :: xl0,xlp1,xlp2,xlp3      !
  real(dp), dimension(:,:,:), allocatable :: xa,xpa,xca            !
  complex(dpc) :: w,tmp1,tmp2,ei                                   !
  complex(dpc), dimension(5) :: ww,vv,uu,pp                        !
  complex(dpc), dimension(5,2) :: svt                              !
  complex(dpc), dimension(5,6) :: svs                              !
  complex(dpc), dimension(:), allocatable :: s0,sm1,sp1,sm2,sp2    !
  complex(dpc), dimension(:), allocatable ::  eia                  !
  complex(dpc), dimension(2) :: stm2,stm1,stp1,stp2,sr0,ctmp3      !
  complex(dpc), dimension(4) :: ctmp1,ctmp2                        !  
  complex(dpc), dimension(6) :: ssm2,ssm1,ss0,ssp1,ssp2            !
  complex(dpc), dimension(2,2) :: s2tmp,ytor,yrad,prad             !
  complex(dpc), dimension(4,4) :: s4tmp,ysphng                     !
  complex(dpc), dimension(6,4) :: s6tmp,ysph                       !
  complex(dpc), dimension(:,:), allocatable :: ur,ut,up            !
  !----------------------------------------------------------------!

  character(len=:), allocatable :: ofile
  character(len=:), allocatable :: mfile


  !--------------------------------------------------!
  !      set up the normalization parameters         !
  !--------------------------------------------------!
  call set_parameters


  !--------------------------------------------------!
  !       deal with the command line arguments       !
  !--------------------------------------------------!
  
  if(.not.found_command_argument("-ofile",ofile)) then
     stop "-ofile is required input"
  end if

  if(.not.found_command_argument("-mfile",mfile)) then
     stop "-mfile is required input"
  else
     inquire(file=trim(mfile), exist=ltmp)
     if(.not.ltmp) stop "model file does not exist"
  end if

  

  


stop
  
  !--------------------------------------------------!
  !          open and read parameter file            !
  !--------------------------------------------------!
  call get_command_argument(1,param_file)
  inquire(file=trim(param_file),exist=ltmp)
  if(.not.ltmp) stop 'can''t find parameter file'

  open(io1,file=param_file,form='formatted', & 
       action='read',iostat=ios)
  if(ios /= 0) stop 'problem opening parameter file'

  ! get the prefix for the output files
  call read_string(io1,pref_out,4)

  ! read the model file
  call read_string(io1,model_file,2)
  open(io2,file=model_file,form='formatted', & 
       action='read',iostat=ios)
  if(ios /= 0) stop ' problem opening model file'
  call read_model(io2);
  close(io2)


  ! read computational parameters
  call read_int(io1,ats_in,2)
  call read_int(io1,grav_switch,2)
  call read_int(io1,out_switch,2)
  call read_int(io1,cor_switch,2)
  call read_int(io1,lmin,2)
  call read_int(io1,lmax,2)
  call read_float(io1,f1,2)
  call read_float(io1,f2,2)
  call read_float(io1,tout,2)
  call read_float(io1,dt,2)
  call read_float(io1,f11,2)
  call read_float(io1,f12,2)
  call read_float(io1,f21,2)
  call read_float(io1,f22,2)


  

  !  read source parameters
  call read_float(io1,rs,2)
  depth_source = rs
  call read_float(io1,lats,2)
  call read_float(io1,lons,2)
  do i=1,6
     call read_float(io1,mm(i),2)
  end do

  ! read receiver parameters
  call read_float(io1,rr,2)
  depth_rec = rr
  call read_int(io1,nr,2)
  allocate(latr(nr),lonr(nr))
  do i=1,2
     read(io1,*)
  end do
  do i=1,nr     
     read(io1,*) latr(i),lonr(i)
  end do
  
  close(io1)

  !--------------------------------------------------!
  !     compute epicentral angle and azimuths        !
  !--------------------------------------------------!
  allocate(delta(nr),azst(nr),azep(nr))
  if(lats == 90.0_dp) lats = 89.98_dp
  do i = 1,nr
     call delaz(latr(i),lonr(i),lats,lons,delta(i),azep(i),azst(i)) 
  end do
  azep = twopi_d-azep
  azst = pi_d-azst

  !--------------------------------------------------!
  !    compute the required spherical harmonics      !
  !--------------------------------------------------!
  allocate(xa(lmax+1,3,nr),xpa(lmax+1,3,nr),xca(lmax+1,3,nr))

  ! calculate all the legendre polynomials
  do i = 1,nr     
     call legendre(delta(i),0,0,xa(1,1:1,i), & 
                                xpa(1,1:1,i), & 
                                xca(1,1:1,i))                     
     call legendre(delta(i),1,1,xa(2,1:2,i), & 
                                xpa(2,1:2,i), & 
                                xca(2,1:2,i))                     
     do l = 2,lmax        
        call legendre(delta(i),l,2,xa(l+1,:,i), & 
             xpa(l+1,:,i), & 
             xca(l+1,:,i))                     
     end do
  end do
  
  ! get the complex exponentials
  allocate(eia(nr))
  eia(:)  = exp(ii*azst(:))

  !--------------------------------------------------!
  !                normalise the inputs              !
  !--------------------------------------------------!
  f1 = f1*t_norm/1000.0_dp
  f2 = f2*t_norm/1000.0_dp
  dt = dt/t_norm
  rs = r(nknot)-rs*1000.0_dp/r_norm
  rr = r(nknot)-rr*1000.0_dp/r_norm
  tout=tout*60.0_dp/t_norm
  mm = mm/moment_norm
  f11=f11/1000.0_dp*t_norm
  f12=f12/1000.0_dp*t_norm
  f21=f21/1000.0_dp*t_norm
  f22=f22/1000.0_dp*t_norm

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

  !----------------------------------------------------!
  !                   solve the equations              !
  !----------------------------------------------------!

  ! initialize displacement vector components
  allocate(ur(nw,nr),ut(nw,nr),up(nw,nr))
  ur = 0.0_dp
  ut = 0.0_dp
  up = 0.0_dp
  


  ! set up dummy source vectors
  s2tmp(:,:) = 0.0_dp
  s2tmp(1,1) = 1.0_dp
  s2tmp(2,2) = 1.0_dp

  s6tmp(:,:) = 0.0_dp
  s6tmp(1,1)   = 1.0_dp
  s6tmp(2,2)   = 1.0_dp
  s6tmp(4,3)   = 1.0_dp
  s6tmp(5,4)   = 1.0_dp

  s4tmp(:,:) = 0.0_dp
  do i = 1,4
     s4tmp(i,i) = 1.0_dp
  end do


  ! set the attenuation 
  call set_ats_switch(ats_in) 


  ! start loop over l
  l = lmin-1
  do i=1,lmax-lmin+1

     l = l+1
     zeta  = sqrt(real(l*(l+1)))
     write(6,'(" calculations for l = ",i5,", to go = ",i5)') l,lmax-l+1
     
     ! start loop over omega
     wloop: do j=1,nw       
 
        w=wmin+(j-1)*dw-ii*ep

        if(l == 0) then

           !-----------------------------------------!
           !               radial modes              !
           !-----------------------------------------!
           
           !get the radial source vector
           call source_vector_rad(w,is,rs,mm,sr0)

           ! solve the radial equations
           if(grav_switch == 0) then
              call rad_solver_ng(w,rr,rs,ir,is,s2tmp,yrad)
              prad = 0.0_dp
           else
              call rad_solver_wg(w,rr,rs,ir,is,s2tmp,yrad,prad)
           end if

           ! compute the uu coefficents
           uu(3) = (sr0(1)*yrad(1,1) + sr0(2)*yrad(1,2))/rr 
           

           if(out_switch == 1) then
              uu(3) = ii*w*uu(3)              
           elseif(out_switch == 2) then
              uu(3) = -w*w*uu(3)              
           end if

           ! sum the spherical harmonic series
           do k = 1,nr
              ur(j,k) = ur(j,k) + uu(3)*xa(1,1,k)              
           end do

        else               

           !--------------------------------!
           !        toroidal modes          !
           !--------------------------------!

           if(l > 1) then
           
              ! get the source vectors
              call source_vector_tor(l,w,is,rs,mm,svt)
              
              ! solve the toroidal equations
              call tor_solver(l,w,rr,rs,ir,is,s2tmp,ytor)
              
              if(out_switch == 1) then
                 ytor = ii*w*ytor
              elseif(out_switch == 2) then
                 ytor = -w*w*ytor
              end if
              
              ! put together the W coefficients
              ww(:) = (svt(:,1)*ytor(1,1)+svt(:,2)*ytor(1,2))/(zeta*rr)
              
              do k = 1,nr             
                 do m = -2,2
                    if(m == 0) cycle
                    if(m >= 0) then
                       x  =  xa(l+1,m+1,k)
                       xp = xpa(l+1,m+1,k)
                       xc = xca(l+1,m+1,k)
                    else
                       msign = (-1)**m
                       x  =  msign*xa(l+1,abs(m)+1,k)
                       xp = msign*xpa(l+1,abs(m)+1,k)
                       xc = msign*xca(l+1,abs(m)+1,k)
                    end if
                    ei = eia(k)**m                 
                    ut(j,k) = ut(j,k) + ii*m*ww(m+3)*x*ei              
                    up(j,k) = up(j,k) - ww(m+3)*xp*ei
                 end do
              end do
              
           end if

           !------------------------------!
           !        spheroidal modes      !
           !------------------------------!
           
           ! get the spheroidal source vectors
           call source_vector_sph(l,w,is,rs,mm,svs)

           ! solve the spheroidal equations
           if(fic) then
              ! solver with fluid inner core
              if(grav_switch == 0) then
                 call sph_solver_ng_fic(l,w,rr,rs,ir,is,s4tmp,ysphng)
              elseif(grav_switch == 1) then
                 call sph_solver_cow_fic(l,w,rr,rs,ir,is,s4tmp,ysphng)
              elseif(grav_switch == 2) then
                 call sph_solver_fic(l,w,rr,rs,ir,is,s6tmp,ysph)
              end if
           else
              ! solver with solid inner core
              if(grav_switch == 0) then
                 call sph_solver_ng(l,w,rr,rs,ir,is,s4tmp,ysphng)
              elseif(grav_switch == 1) then
                 call sph_solver_cow(l,w,rr,rs,ir,is,s4tmp,ysphng)
              elseif(grav_switch == 2) then
                 call sph_solver(l,w,rr,rs,ir,is,s6tmp,ysph)
              end if
           end if
           
           if(grav_switch /= 2) then
              ysph(1,:) = ysphng(1,:)
              ysph(2,:) = ysphng(2,:)
              ysph(3,:) = 0.0_dp
              ysph(4,:) = ysphng(3,:)
              ysph(5,:) = ysphng(4,:)
              ysph(6,:) = 0.0_dp
           end if
           

           if(out_switch == 1) then
              ysph(1,:) = ii*w*ysph(1,:)
              ysph(2,:) = ii*w*ysph(2,:)
              ysph(3,:) = ii*w*ysph(3,:)
           elseif(out_switch == 2) then              
              if(cor_switch == 0) then
                 ysph(1,:) = -w*w*ysph(1,:)
                 ysph(2,:) = -w*w*ysph(2,:)
                 ysph(3,:) = -w*w*ysph(3,:)
              else
                 ! make potential and tilt corrections
                 ctmp1 = -((l+1)*ysph(3,:)+2.0_dp*grav(ir)*ysph(1,:))/rr
                 ctmp2 = zeta*(ysph(3,:)+grav(ir)*ysph(1,:))/rr
                 ysph(1,:) = -w*w*ysph(1,:)+ctmp1
                 ysph(2,:) = -w*w*ysph(2,:)+ctmp2
              end if
           end if

           ! compute the U and V coefficients
           uu(:) = (ysph(1,1)*svs(:,1) + ysph(1,2)*svs(:,2) & 
                   +ysph(1,3)*svs(:,4) + ysph(1,4)*svs(:,5))/rr

           vv(:) = (ysph(2,1)*svs(:,1) + ysph(2,2)*svs(:,2) & 
                   +ysph(2,3)*svs(:,4) + ysph(2,4)*svs(:,5))/(zeta*rr)

           do k = 1,nr             
              do m = -2,2
                 if(m >= 0) then
                    x  =  xa(l+1,m+1,k)
                    xp = xpa(l+1,m+1,k)
                    xc = xca(l+1,m+1,k)
                 else
                    msign = (-1)**m
                    x  =  msign*xa(l+1,abs(m)+1,k)
                    xp = msign*xpa(l+1,abs(m)+1,k)
                    xc = msign*xca(l+1,abs(m)+1,k)
                 end if
                 ei = eia(k)**m                            
                 ur(j,k)  = ur(j,k)  + uu(m+3)*x*ei
                 ut(j,k)  = ut(j,k)  + vv(m+3)*xp*ei              
                 up(j,k)  = up(j,k)  + ii*m*vv(m+3)*xc*ei                     
              end do
           end do           

        end if
     end do wloop
     ! end loop over frequency
     
<<<<<<< HEAD
     ! inverse transform phi coefficients
     if(phi_coef_switch == 1) then

        ! perform the inverse fourier transformation
        call ifft_phi(nt,i1,i2,df,f11,f12,f12,f22, & 
             dt,phi_coef,phi_coef_time)

        ! undo the exponential decay
        do k=1,nt
           t = (k-1)*dt
           phi_coef_time(k,:)  = phi_coef_time(k,:)*exp(ep*t)
        end do

        ! dimensionalize the potential
        phi_coef_time = phi_coef_time*pot_norm

        ! undo the exponential decay
        do k=1,nt
           t = (k-1)*dt
           if(t > tout) exit
           write(io3,*) t*t_norm,real(phi_coef_time(k,:)), & 
                                 imag(phi_coef_time(k,:))
           
        end do

        
        close(io3)


     end if


=======
>>>>>>> main
  end do
  ! end loop over l


  !------------------------------------------!
  !   perform inverse Fourier transform      !
  !------------------------------------------!
  call ifft(nt,i1,i2,df,f11,f12,f21,f22,dt,ur,ut,up)
  nt = size(ur,1)
  
  !--------------------------------------------!
  ! undo the eponential decay and unnormalize  !
  ! the time-series                            !
  !--------------------------------------------!
  do i=1,nt
     t = (i-1)*dt
        ur(i,:)  = ur(i,:)*exp(ep*t)
        ut(i,:)  = ut(i,:)*exp(ep*t)
        up(i,:)  = up(i,:)*exp(ep*t)        
  end do


  ! dimensionalize the outputs
  if(out_switch == 0) then
     ur = ur*r_norm
     ut = ut*r_norm
     up = up*r_norm     
  elseif(out_switch == 1) then
     ur = ur*vel_norm
     ut = ut*vel_norm
     up = up*vel_norm     
  elseif(out_switch == 2) then
     ur = ur*acl_norm
     ut = ut*acl_norm
     up = up*acl_norm     
  end if

  !----------------------------------------------------------!
  ! rotate acceleration vector into geographic Up, North,    !
  ! East, co-ordinate system                                 !
  !----------------------------------------------------------!
  do i = 1,nt
     do  j = 1,nr
        tmp1 = ut(i,j) 
        tmp2 = up(i,j) 
        ut(i,j) = -tmp1*cos(azep(j))+tmp2*sin(azep(j))
        up(i,j) =  tmp1*sin(azep(j))+tmp2*cos(azep(j))
     end do
  end do

  !-------------------------------------------!
  !      write out seismograms to file        !
  !-------------------------------------------!

  nout = nt
  do i = 1,nt
     t = (i-1)*dt
     if(t > tout) then
        nout = i
        exit
     end if
  end do

  do k = 1,nr

     call string_cat_int(trim(pref_out)//'.',k,seis_out)
     open(io1,file=trim(seis_out),form='formatted')

     do i=1,nout
        t = (i-1)*dt
        t = t*t_norm
        write(io1,100) t,real(ur(i,k)),real(ut(i,k)), & 
             real(up(i,k))
     end do
     close(io1)

  end do
  
100 format(5e15.6)
  
  !-----------------------------------------!
  !         deallocate the model            !
  !-----------------------------------------!
  call allocate_model(nknot,1)

 
end program yspec
