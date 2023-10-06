program yspec_pro
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
  integer(i4b), parameter :: io1=7,io2=8                           !
  ! model file name:                                               !
  character(len=256) :: model_file                                 !
  ! parameter file                                                 !
  character(len=256) :: param_file                                 !
  character(len=256) :: pre_pref = 'yspec_pre.'                    !
  character(len=256) :: cal_pref = 'yspec_cal.'                    !
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  
  !================================================================!
  ! variables local to the main program:                           !
  !================================================================!  
  character(len=10) :: string,string1                              !
  character(len=256) :: pref_out                                   !
  character(len=256) :: seis_out
  logical(lgt) :: ltmp                                             !  
  integer(i4b) :: i,j,k,l,ierr,lmin,lmax, &                        !
       nw,ios,i1,i2,ians,is,ir,nt,nr,nout,sout,mex,qex, & 
       njob,ijob,i11,i22,ats_in,grav_switch, & 
       out_switch,cor_switch
  real(dp) :: wmin,wmax,dw,wl,rr,rs,ep,lons,lats, &                !
       zeta,wt,dt,df,t,f11,f12,f21,f22,f,tout,depth_rec, &         !
       depth_source,sqzm2,sqzm6,sqz,f1,f2,mpar,ft1,ft2
  real(dp), dimension(6) :: mm                                     !
  real(dp), dimension(:), allocatable :: latr,lonr,delta, &        !
       azst,csc,azep                                               !
  real(dp), dimension(:,:), allocatable :: xl0,xlp1,xlp2,xlp3      !
  complex(dpc) :: w,tmp1,tmp2                                !
  complex(dpc), dimension(:), allocatable :: s0,sm1,sp1,sm2,sp2    !
  complex(dpc), dimension(:), allocatable :: eiaz,ieiaz,dtylm2, &  !
       dtylm1,dtyl0,dtylp1,dtylp2,dpylm2,dpylm1,dpylp1,dpylp2      !
  complex(dpc), dimension(:,:), allocatable ::  ylm3,ylm2,ylm1, &  !
       yl0,ylp1,ylp2,ylp3                                          !
  complex(dpc), dimension(2) :: stm2,stm1,stp1,stp2,sr0,ctmp3      !
  complex(dpc), dimension(4) :: ctmp1,ctmp2                        !  
  complex(dpc), dimension(6) :: ssm2,ssm1,ss0,ssp1,ssp2            !
  complex(dpc), dimension(2,2) :: s2tmp,ytor,yrad,prad             !
  complex(dpc), dimension(6,4) :: s6tmp,ysph                       !
  complex(dpc), dimension(:,:), allocatable :: ur,ut,up            !
  !----------------------------------------------------------------!
  

  !--------------------------------------------------!
  !      set up the normalization parameters         !
  !--------------------------------------------------!
  call set_parameters
  

  !------------------------!
  ! get the number of jobs !
  !------------------------!
  call get_command_argument(1,string)
  read(string,*) njob


  ! read in the parameters from the first job
  ijob = 1
  write(string,'(i10)') ijob
  j = floor(log10(real(ijob)))
  open(io1,file=trim(pre_pref)//string(10-j:10), & 
          form='unformatted')

  read(io1) pref_out
  read(io1) model_file
  read(io1) lmin,lmax
  read(io1) ft1,ft2
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
  read(io1) i11,i22
  read(io1) nt
  read(io1) ats_in
  read(io1) grav_switch
  read(io1) out_switch
  read(io1) cor_switch
  close(io1)
  

  ! get the minimun frequency
  f1 = ft1
  i1 = i11
  
  if(njob == 1) then
     f2 = ft2
     i2 = i22
  end if



  ! allocate the data arrays
  allocate(ur(nt,1),ut(nt,1),up(nt,1))




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


  ! write out a file with epicentral angles 
  open(io1,file=trim(pref_out)//'.angles', & 
       form='formatted')
  do i = 1,nr
     write(io1,*) delta(i)*180.0_dp/pi_d
  end do
  close(io1)
  


  ! start the loop over the receivers
  do ir = 1,nr

     ur = 0.0_dp
     up = 0.0_dp
     ut = 0.0_dp


     ! loop over the different job files, reading in
     ! the spectra
     do ijob = 1,njob

          write(string,'(i10)') ijob
          j = floor(log10(real(ijob)))
          open(io1,file=trim(pre_pref)//string(10-j:10), & 
               form='unformatted')          
          read(io1) pref_out
          read(io1) model_file
          read(io1) lmin,lmax
          read(io1) ft1,ft2
          read(io1) tout
          read(io1) dt
          read(io1) f11,f12,f21,f22
          read(io1) rs
          read(io1) lats,lons
          read(io1) mm
          read(io1) rr
          read(io1) nr
          read(io1) latr,lonr
          read(io1) ep
          read(io1) df
          read(io1) i11,i22
          read(io1) nt
          close(io1)



          ! get the maximum frequency
          if(ijob == njob .and. ir == 1) then
             f2 = ft2
             i2 = i22
          end if

          write(string1,'(i10)') ir
          k = floor(log10(real(ir)))          
          open(io1,file=trim(cal_pref)//string(10-j:10) & 
               //'.'//string1(10-k:10),form='unformatted')    
          read(io1) ur(i11-i1+1:i22-i1+1,1), & 
                    ut(i11-i1+1:i22-i1+1,1), & 
                    up(i11-i1+1:i22-i1+1,1)    
          close(io1)

       end do
       

       !------------------------------------------!
       !   perform inverse Fourier transform      !
       !------------------------------------------!


       call ifft(nt,i1,i2,df,f11,f12,f21,f22,dt,ur,ut,up)
       nt = size(ur,1)



  
       !--------------------------------------------!
       ! undo the eponential decay and unnormalize  !
       ! the accelerations                          !
       !--------------------------------------------!
       do i=1,nt
          t = (i-1)*dt
          ur(i,1) = ur(i,1)*exp(ep*t)
          ut(i,1) = ut(i,1)*exp(ep*t)
          up(i,1) = up(i,1)*exp(ep*t)
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
          tmp1 = ut(i,1) 
          tmp2 = up(i,1) 
          ut(i,1) = -tmp1*cos(azep(ir))+tmp2*sin(azep(ir))
          up(i,1) =  tmp1*sin(azep(ir))+tmp2*cos(azep(ir))
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


         call string_cat_int(trim(pref_out)//'.',ir,seis_out)
         open(io1,file=trim(seis_out),form='formatted')

         do i=1,nout
            t = (i-1)*dt
            t = t*t_norm
            write(io1,100) t,real(ur(i,1)),real(ut(i,1)), & 
                    real(up(i,1)) 
         end do
         close(io1)
               
  
100      format(4e15.6)



       


      end do
  

  
  !-----------------------------------------!
  !         deallocate the model            !
  !-----------------------------------------!
  call allocate_model(nknot,1)




end program yspec_pro
