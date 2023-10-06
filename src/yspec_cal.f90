program yspec_cal
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
  !----------------------------------------------------------------!
  implicit none                                                    !
  !----------------------------------------------------------------!

  !================================================================!
  ! input/output files used:                                       !
  !================================================================!
  integer(i4b), parameter :: io1=7,io2=8                           !
  ! model file name:                                               !
  character(len=256) :: model_file                                 !
  ! parameter file                                                 !
  character(len=256) :: param_file                                 !
  character(len=256) :: pre_pref = 'yspec_pre.'
  character(len=256) :: cal_pref = 'yspec_cal.'
  !----------------------------------------------------------------!

    
  !================================================================!
  ! variables local to the main program:                           !
  !================================================================!  
  character(len=10) :: string,string1
  character(len=256) :: pref_out                                   
  logical(lgt) :: ltmp                                             
  integer(i4b) :: i,j,k,l,ierr,lmin,lmax, &                        
       nw,ios,i1,i2,ians,is,ir,nt,nr,nout,sout,mex, &              
       qex,m,im,isv,msign,ijob,ats_in,grav_switch, & 
       out_switch,cor_switch
  real(dp) :: wmin,wmax,dw,wl,rr,rs,ep,lons,lats, &                !
       zeta,wt,dt,df,t,f11,f12,f21,f22,f,tout,depth_rec, &         !
       depth_source,sqzm2,sqzm6,sqz,f1,f2,mpar,x,xp,xc             !
  real(dp), dimension(6) :: mm                                     !
  real(dp), dimension(:), allocatable :: latr,lonr,delta, &        !
       azst,csc,azep                                               !
  real(dp), dimension(:,:), allocatable :: xl0,xlp1,xlp2,xlp3      !
  real(dp), dimension(:,:,:), allocatable :: xa,xpa,xca            !
  complex(dpc) :: w,tmp1,tmp2,ei                                   !
  complex(dpc), dimension(5) :: ww,vv,uu                           !
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

  !--------------------------------------------------!
  !      set up the normalization parameters         !
  !--------------------------------------------------!
  call set_parameters
  
  !-----------------------!
  ! get the job id number !
  !-----------------------!
  call get_command_argument(1,string)
  read(string,*) ijob


  !-----------------------------------------------! 
  ! open and read the parameter file for the job  !
  !-----------------------------------------------!
  write(string,'(i10)') ijob
  j = floor(log10(real(ijob)))
  open(io1,file=trim(pre_pref)//string(10-j:10), & 
          form='unformatted')

  read(io1) pref_out
  read(io1) model_file
  read(io1) lmin,lmax
  read(io1) f1,f2
  read(io1) tout
  read(io1) dt
  read(io1) f11,f12,f21,f22
  read(io1) rs
  read(io1) lats,lons
  read(io1) mm
  read(io1) rr
  read(io1) nr
  allocate(latr(nr),lonr(nr))
  read(io1) latr,lonr
  read(io1) ep
  read(io1) df
  read(io1) i1,i2
  read(io1) nt
  read(io1) ats_in
  read(io1) grav_switch
  read(io1) out_switch
  read(io1) cor_switch
  close(io1)

  ! read the model file  
  open(io2,file=model_file,form='formatted', & 
       action='read',iostat=ios)
  if(ios /= 0) stop ' problem opening model file'
  call read_model(io2);
  close(io2)


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
           else
              call rad_solver(w,rr,rs,ir,is,s2tmp,yrad)
           end if


           ! compute the uu coefficent
           uu(3) = 0.0_dp
           uu(1) = 0.0_dp
           do isv = 1,2
              uu(3) = uu(3) + sr0(isv)*yrad(1,isv)/rr
           end do


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
           elseif(out_switch == 2) then              
              if(cor_switch == 0) then
                 ysph(1,:) = -w*w*ysph(1,:)
                 ysph(2,:) = -w*w*ysph(2,:)
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
                 ur(j,k) = ur(j,k) + uu(m+3)*x*ei
                 ut(j,k) = ut(j,k) + vv(m+3)*xp*ei              
                 up(j,k) = up(j,k) + ii*m*vv(m+3)*xc*ei    
              end do
           end do

       end if
     end do wloop
  end do



  ! loop over the different receivers
  do i = 1,nr

     ! open the binary file for the spectra
     write(string,'(i10)') ijob
     write(string1,'(i10)') i
     j = floor(log10(real(ijob)))
     k = floor(log10(real(i)))
     open(io1,file=trim(cal_pref)//string(10-j:10)// & 
          '.'//string1(10-k:10),form='unformatted')
     write(io1) ur(1:i2-i1+1,i),ut(1:i2-i1+1,i),up(1:i2-i1+1,i)
     close(io1)
     
  end do



  close(99)





end program yspec_cal
