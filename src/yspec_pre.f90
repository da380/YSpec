program yspec_pre
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
  !----------------------------------------------------------------!

  
  !================================================================!
  ! variables local to the main program:                           !
  !================================================================!  
  character(len=10) :: string                                      !
  character(len=256) :: pref_out                                   !
  logical(lgt) :: ltmp                                             !  
  integer(i4b) :: i,j,k,l,ierr,lmin,lmax, &                        !
       nw,ios,i1,i2,ians,is,ir,nt,nr,nout,sout,mex,qex, &          !
       njob,ncal,i11,i22,count,ijob,ats_in,grav_switch, & 
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

  integer(i4b), dimension(:), allocatable :: ncalj
  integer(i4b), dimension(:,:), allocatable :: indj
  

  !--------------------------------------------------!
  !      set up the normalization parameters         !
  !--------------------------------------------------!
  call set_parameters
  


  !--------------------------------------------------!
  !          open and read parameter file            !
  !--------------------------------------------------!
  call get_command_argument(1,param_file)
  inquire(file=trim(param_file),exist=ltmp)
  if(.not.ltmp) stop 'can''t find parameter file'

  call get_command_argument(2,string)
  read(string,*) njob

  open(io1,file=param_file,form='formatted', & 
       action='read',iostat=ios)
  if(ios /= 0) stop 'problem opening parameter file'

  

  ! get the prefix for the output files
  call read_string(io1,pref_out,4)


  ! read the model file
  call read_string(io1,model_file,2)
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
  !                normalise the inputs              !
  !--------------------------------------------------!
  f1 = f1*t_norm/1000.0_dp
  f2 = f2*t_norm/1000.0_dp
  dt = dt/t_norm
  rs=1.0_dp-rs*1000.0_dp/r_norm
  rr=1.0_dp-rr*1000.0_dp/r_norm
  tout=tout*60.0_dp/t_norm


  mm = mm/moment_norm

  f11=f11/1000.0_dp*t_norm
  f12=f12/1000.0_dp*t_norm
  f21=f21/1000.0_dp*t_norm
  f22=f22/1000.0_dp*t_norm

  mex = 5
  qex = 1

  call fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)


  ! divide the frequencies between the processors
  allocate(indj(njob,2))


  call jobs(i1,i2,njob,indj)

  do i = 1,njob
     print *, i1,indj(i,:),i2,indj(i,2)-indj(i,1)+1
  end do


  do i = 1,njob

       
     ! open the binary job file     
     write(string,'(i10)') i
     j = floor(log10(real(i)))
     open(io1,file=trim(pre_pref)//string(10-j:10), & 
          form='unformatted')


     i11 = indj(i,1)
     i22 = indj(i,2)

     ft1 = (i11-1)*df
     ft2 = (i22-1)*df


     ! write data to the job file
     write(io1) pref_out
     write(io1) model_file
     write(io1) lmin,lmax
     write(io1) ft1,ft2
     write(io1) tout
     write(io1) dt
     write(io1) f11,f12,f21,f22
     write(io1) rs
     write(io1) lats,lons
     write(io1) mm
     write(io1) rr
     write(io1) nr
     write(io1) latr,lonr
     write(io1) ep
     write(io1) df
     write(io1) i11,i22
     write(io1) nt
     write(io1) ats_in
     write(io1) grav_switch
     write(io1) out_switch
     write(io1) cor_switch

     close(io1)



  end do




  contains








  subroutine jobs(i1,i2,njob,indj)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: i1,i2
    integer(i4b), intent(in) :: njob
    integer(i4b), dimension(njob,2) :: indj


    integer(i4b) :: count,ijob,ncal,mn,mx,imn,imx
    integer(i4b), dimension(njob) :: ncalj
    real(dp) :: r,f,fac


    fac = 0.9_dp
    ncal = i2-i1+1

    ncalj = 2
    count = 2*njob
    ijob = 1
    do      
       if(count == ncal) exit
       call random_number(r)
       f = fac*real(ijob)/real(njob+1)
       if(r >= f) then
          ncalj(ijob) = ncalj(ijob)+1
          count = count+1
       end if
       if(ijob < njob) then
          ijob = ijob+1
       else
          ijob = 1
       end if
    end do

    ! check for null-jobs
    do 
       
       imn = iminloc(ncalj)
       mn = ncalj(imn)

       if(mn > 1) exit

       imx = imaxloc(ncalj)
       mx = ncalj(imx)

       if(mx == 3) stop ' too many jobs'

       ncalj(imn) = 2
       ncalj(imx) = ncalj(imx)-2

    end do


    
    indj(1,1) = i1
    do ijob = 1,njob
       indj(ijob,2) = indj(ijob,1)+ncalj(ijob)-1
       if(ijob < njob) then
          indj(ijob+1,1) = indj(ijob,2)+1
       end if
    end do


    
    
    
    
    return
  end subroutine jobs





end program yspec_pre
