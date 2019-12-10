!
!   others.f90
!   MarsInversion
!
!   Created by fuji on 04/12/2019.
!   Copyright 2019 nfuji. All rights reserved.
!
subroutine pinput

    use parameters
    implicit none
    character(200) :: argv
    integer :: argc
    character(200) :: tmpfile,metafile
    character(200) :: obsfile
    character(200) :: dummy
    !integer, external :: getpid
    integer :: iloop,it
    character(200) :: commandline

    call getarg(1,argv)
    metafile=argv
    call getarg(2,argv)
    workingDir=argv
    call getarg(3,argv)
    resultDir=argv
    call getarg(4,argv)
    inversionName=argv
  
    argc=iargc()
  
    if(argc.ne.4) then
        print *, "you need <metafile> <working directory> <result directory> <inversion name>"
        print *, "cheers"
        stop
    endif

  

  
    !write(tmpfile,"(Z5.5)") getpid()
    !tmpfile='tmpfileMarsInversion'//tmpfile
    tmpfile='tmpfileMarsInversion'
    open(unit=5,file=metafile,status='unknown')
    open(unit=1,file=tmpfile,status='unknown')
100 continue
    read(5,110) dummy
110 format(a200)
    if(dummy(1:1).eq.'#') goto 100
    if(dummy(1:3).eq.'end') goto 120
    write(1,110) dummy
    goto 100
120 continue
    close(1)
    close(5)

    open(unit=1,file=tmpfile,status='unknown')
    calculMode=0
    nmt=6

    read(1,110) dummy ! This dummy string determines i) test mode, ii) Alice normal, or iii) versionSGT mode

    if(dummy(1:10).eq.'versionSGT') then
        calculMode=2
        nmt=6
        read(1,110) SGTinfo
        SGTinfo = trim(SGTinfo)
        read(1,110) parentDir
        parentDir = trim(parentDir)
        read(1,110) eventName
        eventName = trim(eventName)
        read(1,110) stationName
        stationName = trim(stationName)
        read(1,*) stla, stlo
        read(1,*) dt
        read(1,*) tlen
        np=int(tlen/dt)

        read(1,*) npButterworth
        read(1,*) fmin
        read(1,*) fmax
        read(1,*) start, end ! onllly available for RSGT version
        read(1,*) ntwin
        allocate(twin(1:4,1:ntwin))
        allocate(itwin(1:4,1:ntwin))
        do iloop=1,ntwin
            read(1,*) twin(1,iloop),twin(2,iloop),twin(3,iloop),twin(4,iloop)
        enddo
        itwin=int(twin/dt)
        read(1,*) tlenData
        npData = int(tlenData/dt)

        allocate(obsRaw(1:npData,1:3))
        allocate(obsFilt(1:npData,1:3))



        do iloop=1,3
            read(1,110) obsfile
            obsfile=trim(workingDir)//"/"//trim(obsfile)
            open(unit=10,file=obsfile,status='unknown')
            do it=1,npData
                read(10,*) obsRaw(it,iloop)
            enddo
            close(10)
        enddo



        read(1,*) ntwinObs

        ! control the multiple moving windows
        read(1,110) dummy
        if(dummy(1:5).eq.'fixed') then
            NmovingWindowDimension=1
        elseif(dummy(1:).eq.'independent') then
            NmovingWindowDimension=ntwinObs
        endif





        allocate(twinObs(1:4,1:ntwinObs))
        allocate(itwinObs(1:4,1:ntwinObs))
        do iloop=1,ntwinObs
            read(1,*) twinObs(1,iloop),twinObs(2,iloop),twinObs(3,iloop),twinObs(4,iloop)
        enddo
        itwinObs=int(twinObs/dt)



        read(1,*) movingWindowStep
        ntStep=int(movingWindowStep/dt)


        ! Moving windows

        allocate(fMovingWindowStart(1:NmovingWindowDimension))
        allocate(fMovingWindowEnd(1:NmovingWindowDimension))

        allocate(iMovingWindowStart(1:NmovingWindowDimension))
        allocate(iMovingWindowEnd(1:NmovingWindowDimension))

        do iloop=1,NmovingWindowDimension
            read(1,*) fMovingWindowStart(iloop), fMovingWindowEnd(iloop)
        enddo
        iMovingWindowStart=int(fMovingWindowStart/dt)
        iMovingWindowEnd=int(fMovingWindowEnd/dt)
        
        do iloop=1,NmovingWindowDimension
            iMovingWindowStart(iloop)=itwinObs(1,iloop)+iMovingWindowStart(iloop)
            iMovingWindowEnd(iloop)=itwinObs(1,iloop)+iMovingWindowEnd(iloop)
            ! check whether syn and obs are available for these indices
            if(iMovingWindowEnd(iloop)<1) then
                print *, "no sufficient data points in syn data for the window", iloop
                stop
            endif
            if(itwinObs(1,iloop)<1) then
                print *, "no sufficient data points in obs data for the window", iloop
                stop
            endif
            if(iMovingWindowStart(iloop)>np) then
                print *, "no sufficient data points in syn data for the window", iloop
                stop
            endif
            if(itwinObs(4,iloop)>npData) then
                print *, "no sufficient data points in obs data for the window", iloop
                stop
            endif
        enddo






        commandline = 'mkdir -p '//trim(parentDir)
        call system(commandline)
        call pinputDSM(DSMconfFile,PoutputDir,psvmodel,modelname,tlenFull,rmin_,rmax_,rdelta_, &
                r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch, &
                synnswitch,SGTinfo)
        call readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
     


        call readpsvmodel(psvmodel,tmpfile)
        INFO_TSGT = trim(parentDir)//"/INFO_TSGT.TXT"
        INFO_RSGT = trim(parentDir)//"/INFO_RSGT.TXT"
        rsampletxt = trim(parentDir)//"/rsample.txt"
        modelcard = trim(parentDir)//"/"//trim(modelname)//".card"

        !synnfile = trim(parentDir)//"/"//trim(stationName)//"."//trim(eventName)//"."//trim(compo)//"s.dat"

     
        r_n = int((rmax_-rmin_)/rdelta_)+1
        allocate(r_(1:r_n))
        do iloop=1,r_n
            r_(iloop) = rmin_ + dble(iloop-1)*rdelta_
        enddo

        theta_n = int((thetamax-thetamin)/thetadelta)+1
        allocate(thetaD(1:theta_n))
        do iloop=1,theta_n
            thetaD(iloop) = thetamin + dble(iloop-1)*thetadelta
        enddo

        phi_n=int((phimax-phimin)/phidelta)+1
        allocate(phiD(1:phi_n))
        do iloop=1,phi_n
            phiD(iloop) = phimin + dble(iloop-1)*phidelta
        enddo
        
     
        ! lsmoothfinder for FFT
        np0=imax
        call lsmoothfinder(tlenFull,np0,samplingHz,lsmooth)
        iloop=1
        do while (iloop<lsmooth)
            iloop = iloop*2
        enddo
        lsmooth = iloop
        iloop = 0
        np1 = 1
        do while (np1<np0)
            np1 = np1*2
        enddo
        np1 = np1*lsmooth

        ! redefinition of samplingHz
        samplingHz = dble(2*np1)/tlenFull
        dtn = 1.d0/samplingHz
        iWindowStart = int(start*samplingHz)
        iWindowEnd   = int(end*samplingHz)



        ! allocate SGTs, synthetics in frequency
        allocate(omega(imin:imax))
        do iloop = imin, imax
            omega(iloop) = 2.d0*pi*dble(iloop)/tlenFull
        enddo
     
        nConfiguration=r_n*theta_n*phi_n

    elseif((dummy(1:6).eq.'normal').or.(dummy(1:4).eq.'test')) then
        if(dummy(1:4).eq.'test') calculMode=1
            read(1,*) dt
            read(1,*) tlen
            np=int(tlen/dt)
            
     
            read(1,*) npButterworth
            read(1,*) fmin
            read(1,*) fmax
            read(1,*) ntwin
            allocate(twin(1:4,1:ntwin))
            allocate(itwin(1:4,1:ntwin))
            do iloop=1,ntwin
                read(1,*) twin(1,iloop),twin(2,iloop),twin(3,iloop),twin(4,iloop)
            enddo
            itwin=int(twin/dt)
            read(1,*) tlenData
            npData = int(tlenData/dt)
     
            allocate(obsRaw(1:npData,1:3))
            allocate(obsFilt(1:npData,1:3))
     
     
     
            do iloop=1,3
                read(1,110) obsfile
                obsfile=trim(workingDir)//"/"//trim(obsfile)
                open(unit=10,file=obsfile,status='unknown')
                do it=1,npData
                    read(10,*) obsRaw(it,iloop)
                enddo
                close(10)
            enddo


            read(1,*) ntwinObs

            ! control the multiple moving windows
            read(1,110) dummy
            if(dummy(1:5).eq.'fixed') then
                NmovingWindowDimension=1
            elseif(dummy(1:).eq.'independent') then
                NmovingWindowDimension=ntwinObs
            endif


            


            allocate(twinObs(1:4,1:ntwinObs))
            allocate(itwinObs(1:4,1:ntwinObs))
            do iloop=1,ntwinObs
                read(1,*) twinObs(1,iloop),twinObs(2,iloop),twinObs(3,iloop),twinObs(4,iloop)
            enddo
            itwinObs=int(twinObs/dt)

            read(1,*) movingWindowStep
            ntStep=int(movingWindowStep/dt)


            ! Moving windows

            allocate(fMovingWindowStart(1:NmovingWindowDimension))
            allocate(fMovingWindowEnd(1:NmovingWindowDimension))

            allocate(iMovingWindowStart(1:NmovingWindowDimension))
            allocate(iMovingWindowEnd(1:NmovingWindowDimension))

            do iloop=1,NmovingWindowDimension
                read(1,*) fMovingWindowStart(iloop), fMovingWindowEnd(iloop)
            enddo
            iMovingWindowStart=int(fMovingWindowStart/dt)
            iMovingWindowEnd=int(fMovingWindowEnd/dt)
            
            do iloop=1,NmovingWindowDimension
                iMovingWindowStart(iloop)=itwinObs(1,iloop)+iMovingWindowStart(iloop)
                iMovingWindowEnd(iloop)=itwinObs(1,iloop)+iMovingWindowEnd(iloop)
                ! check whether syn and obs are available for these indices
                if(iMovingWindowEnd(iloop)<1) then
                    print *, "no sufficient data points in syn data for the window", iloop
                    stop
                endif
                if(itwinObs(1,iloop)<1) then
                    print *, "no sufficient data points in obs data for the window", iloop
                    stop
                endif
                if(iMovingWindowStart(iloop)>np) then
                    print *, "no sufficient data points in syn data for the window", iloop
                    stop
                endif
                if(itwinObs(4,iloop)>npData) then
                    print *, "no sufficient data points in obs data for the window", iloop
                    stop
                endif
            enddo


    
            read(1,*) nConfiguration
            allocate(filenames(nConfiguration))
            do iloop=1,nConfiguration
                read(1,110) filenames(iloop)
                filenames(iloop)=trim(workingDir)//"/"//trim(filenames(iloop))
            enddo
        endif
     
     
    endif
  
    close(1)

end subroutine pinput


subroutine pinputDSM(DSMconfFile,outputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch,SGTinfo)
  implicit none
  !character(120), parameter :: tmpfile='tmpworkingfile_for_SGTcalcul'
  character(120) :: dummy,outputDir,psvmodel,modelname,DSMconfFile,SGTinfo
  real(kind(0d0)) :: tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta
  real(kind(0d0)) :: thetamin,thetamax,thetadelta
  integer :: imin,imax,rsgtswitch,tsgtswitch,synnswitch
  integer, external :: getpid
  character(120) :: tmpfile

  !write(tmpfile, "(Z5.5)") getpid()
  !tmpfile='tmpworkingfile_for_SGTcalcul'//tmpfile
  tmpfile='tmpworkingfile_for_SGTcalcul'

  open(unit=2, file=SGTinfo)
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)

  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) DSMconfFile
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
  outputDir=trim(outputDir)
  psvmodel=trim(psvmodel)
  modelname=trim(modelname)
  read(1,*) tlen
  read(1,*) rmin_,rmax_,rdelta_
  read(1,*) r0min
  r0max=r0min
  r0delta=20.d0
  read(1,*) thetamin,thetamax,thetadelta
  read(1,*) imin,imax
  read(1,*) rsgtswitch,tsgtswitch,synnswitch
  close(1,status='delete')

end subroutine pinputDSM

subroutine readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
  implicit none
  !character(120), parameter :: tmpfile='tmpworkingfile_for_DSMconf'
  character(120) :: dummy,DSMconfFile
  real(kind(0d0)) :: re,ratc,ratl,omegai
  integer  :: maxlmax
  integer, external :: getpid
  character(120) :: tmpfile

  !write(tmpfile,"(Z5.5)") getpid()
  !tmpfile='tmpworkingfile_for_DSMconf'//tmpfile
  tmpfile='tmpworkingfile_for_DSMconf'

  open(unit=2, file=DSMconfFile, status='old',action='read',position='rewind')
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)


  open(unit=1,file=tmpfile,status='unknown')
  read(1,*) re
  read(1,*) ratc
  read(1,*) ratl
  read(1,*) omegai
  read(1,*) maxlmax
  close(1,status='delete')


end subroutine readDSMconf


subroutine readpsvmodel(psvmodel,tmpfile)
  implicit none
  character(120) :: psvmodel, tmpfile, dummy
  open(unit=2, file=psvmodel, status='old',action='read',position='rewind')
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)
end subroutine readpsvmodel


subroutine translat(geodetic,geocentric)

  implicit none
  real(kind(0d0)),parameter ::  flattening = 1.d0 / 298.25d0
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0
  real(kind(0d0)) :: geocentric, geodetic
  integer :: flag
  flag = 0
  if(geodetic .gt. 90.d0) then
     geodetic = 1.8d2 - geodetic
     flag = 1
  endif

  geodetic = geodetic / 1.8d2 * pi
  geocentric = datan((1.d0-flattening)*(1.d0-flattening)* dtan(geodetic) )
  geocentric = geocentric * 1.8d2 / pi

  if(flag .eq. 1) then
     geocentric = 1.8d2 - geocentric
  endif

  return
end subroutine translat


subroutine lsmoothfinder (tlen, np0, freq, lsmooth)

  implicit none
  real(kind(0d0)) :: tlen, freq
  integer :: np0, np, lsmooth, i
  np = 1
  do while (np<np0)
     np = np*2
  enddo
  lsmooth = int(0.5*tlen*freq/dble(np))
  i = 1

  do while (i<lsmooth)
     i = i*2
  enddo

  lsmooth = i

  return

end subroutine lsmoothfinder


subroutine inverse(aa,c,n)
  !============================================================
  ! Inverse matrix
  ! Method: Based on Doolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed
  ! during the calculation
  !===========================================================
  implicit none
  integer :: n
  real(kind(0d0)) :: a(n,n), c(n,n),aa(n,n)
  real(kind(0d0)) :: L(n,n), U(n,n), b(n), d(n), x(n)
  real(kind(0d0)) :: coeff
  integer :: i, j, k
  
  a=aa
  
  ! step 0: initialization for matrices L and U and b
  ! Fortran 90/95 aloows such operations on matrices
  L=0.0
  U=0.0
  b=0.0
  
  ! step 1: forward elimination
  do k=1, n-1
     do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
           a(i,j) = a(i,j)-coeff*a(k,j)
        enddo
     enddo
  enddo
  
  ! Step 2: prepare L and U matrices
  ! L matrix is a matrix of the elimination coefficient
  ! + the diagonal elements are 1.0
  do i=1,n
     L(i,i) = 1.0
  end do
  ! U matrix is the upper triangular part of A
  do j=1,n
     do i=1,j
        U(i,j) = a(i,j)
     enddo
  enddo
  
  ! Step 3: compute columns of the inverse matrix C
  do k=1,n
     b(k)=1.0
     d(1) = b(1)
     ! Step 3a: Solve Ld=b using the forward substitution
     do i=2,n
        d(i)=b(i)
        do j=1,i-1
           d(i) = d(i) - L(i,j)*d(j)
        enddo
     enddo
     ! Step 3b: Solve Ux=d using the back substitution
     x(n)=d(n)/U(n,n)
     do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
           x(i)=x(i)-U(i,j)*x(j)
        enddo
        x(i) = x(i)/u(i,i)
     enddo
     ! Step 3c: fill the solutions x(n) into column k of C
     do i=1,n
        c(i,k) = x(i)
     enddo
     b(k)=0.0
  enddo
end subroutine inverse

subroutine invbyCG(nd,ata,atd,eps,x)
  ! modified from invbyCG Jun. 2009 Nobuaki Fuji
  !
  !                       Nov. 2019 Nobuaki Fuji

  implicit none
  integer :: nd, ii
  real(kind(0d0)) :: ata(1:nd,1:nd)
  real(kind(0d0)) :: x(1:nd),r(1:nd),w(1:nd),z(1:nd),x0(1:nd),atd(1:nd)
  real(kind(0d0)) :: eps
  real(8) :: a, b, residual,initres,threshold,paap

  x0 = 0.d0
    
  !r = atd - matmul(ata,x0)
  r = atd
  w = -r
  z = matmul(ata,w)
  a = dot_product(r,w) / dot_product(w,z)
  x = x0 +a*w
  b = 0

  initres=dot_product(atd,atd)
  threshold=eps*initres
  

  
  do ii=1,nd

     r = r - a*z
     residual=dot_product(r,r)
   
     if(residual.le.threshold) exit
 
     b = dot_product(r,z)/dot_product (w,z)
     w = -r + b*w
     z = matmul(ata,w)
     
     !for pAAp calculation
     paap = dot_product(w,z)
     !end for pAAp calculation
     
     a = dot_product(r,w)/dot_product(w,z)
     x = x+a*w
     
  enddo


  

end subroutine invbyCG












subroutine bwfilt(x,y,dt,n,irek,norder,f1,f2)
  ! recursive filtering of data with butterworth filter
  ! x: input array
  ! y: output array
  ! dt: time increment
  ! n: number of data points
  ! irek=0: forward filtering only
  ! irek=1: forward and backward filtering
  ! norder: order of butterworth filter
  ! norder=0: only filtering, no determination of coefficients
  ! norder<0: no starplots of transfer function and impulse response
  ! f1: low cutoff frequency (Hz)
  ! f1=0: low pass filter
  ! f2: high cutoff frequency (Hz)
  ! f2>0.5/dt: high pass filter
  implicit none
  real(kind(0d0)), dimension(1) :: x,y
  real(kind(0d0)), dimension(10) :: a,b1,b2
  real(kind(0d0)) :: dt,f1,f2
  integer :: iunit,npoles,norder,irek,n,lx
  
  iunit = 3
  if (norder.ne.0) then
    npoles = iabs(norder)
    !determination of filter coefficients
    call bpcoeff(f1,f2,npoles,dt,a,b1,b2)
    if (norder.ge.0) then
      !plot of transfer function and impuulse response
      lx = 100
      !filtering
    endif
  endif
  if (n.ne.0) then
    call rekurs(x,y,n,a,b1,b2,npoles,irek)
  endif

  return

end subroutine bwfilt

subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
  ! performs recursive filtering of data in array x of length ndat
  ! filtered output in y
  ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
  ! npoles is the number of poles
  ! iflag=0: forward filtering only
  ! iflag.ne.0: forward and backward filtering
  implicit none
  real(kind(0d0)), dimension(10) :: z,z1,z2,a,b1,b2
  real(kind(0d0)) :: x1,x2
  integer :: ndat,npoles,iflag,n,i
  real(kind(0d0)), dimension(ndat) :: x,y
  
  !forward
  x1 = 0.d0
  x2 = 0.d0
  do i = 1,npoles
    z1(i) = 0.d0
    z2(i) = 0.d0
  enddo
  do n = 1,ndat
    z(1) = a(1)*(x(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
    do i = 2,npoles
      z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
    enddo
    x2 = x1
    x1 = x(n)
    do i = 1,npoles
      z2(i) = z1(i)
      z1(i) = z(i)
    enddo
    y(n) = z(npoles)
  enddo
  if (iflag.eq.0) then
    return
  endif
  !backward
  x1 = 0.d0
  x2 = 0.d0
  do i = 1,npoles
    z1(i) = 0.d0
    z2(i) = 0.d0
  enddo
  do n = ndat,1,-1
    z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
    do i = 2,npoles
      z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
    enddo
    x2 = x1
    x1 = y(n)
    do i = 1,npoles
       z2(i) = z1(i)
       z1(i) = z(i)
    enddo
    y(n) = z(npoles)
  enddo

  return

end subroutine rekurs

subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
  !determines filtercoefficients for recursive bandpassfilter
  implicit none
  real(kind(0d0)), dimension(10) :: a,b1,b2
  complex(kind(0d0)), dimension(20) :: s
  complex(kind(0d0)) :: t1,t2,p
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: f1,f2,dt,d2,w0,w1,w2,ssum,sprod,fact1,fact2,fact3
  integer :: i,npol2,n,npoles
  
  if (npoles.gt.10) stop 'npoles greater than 10: STOP'
  d2 = 2.d0/dt
  w1 = d2*tan(2.d0*pi*f1/d2)
  w2 = d2*tan(2.d0*pi*f2/d2)
  w0 = 0.5*(w2-w1)
  i = 1
  npol2 = npoles/2+1
  do n = 1,npoles
    p = exp(dcmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
    t1 = p*dcmplx(w0,0.d0)
    t2 = sqrt(t1*t1-dcmplx(w1*w2,0.d0))
    s(i) = t1+t2
    s(i+1) = t1-t2
    i = i+2
  enddo
  do n = 1,npoles
    ssum = 2*real(s(n))
    sprod = dble(s(n)*conjg(s(n)))
    fact1 = d2*d2-d2*ssum+sprod
    fact2 = 2.d0*(sprod-d2*d2)
    fact3 = d2*d2+d2*ssum+sprod
    a(n) = 2.d0*d2*w0/fact1
    b1(n) = fact2/fact1
    b2(n) = fact3/fact1
  enddo

  return

end subroutine bpcoeff
