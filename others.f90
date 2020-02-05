!
!   others.f90
!   MarsInversion
!
!   Created by fuji on 04/12/2019.
!   Copyright 2019 nfuji. All rights reserved.
!
subroutine pinput

    use parameters
    use angles
    implicit none
    character(200) :: argv
    integer :: argc
    character(200) :: tmpfile,metafile
    character(200) :: obsfile
    character(200) :: dummy
    real(kind(0d0)) :: fdummy,fdummy2,phirq
    !integer, external :: getpid
    integer :: iloop,it,icheck,jloop
    !character(200) :: commandline
    character(200) :: paramName

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
        close(1)
        calculMode=2
        nmt=6

        paramName="SGTinfo"
        call searchForParams(tmpfile,paramName,SGTinfo,0)
        print *, "SGTinfo is ", SGTinfo

        paramName="parentDir"
        call searchForParams(tmpfile,paramName,parentDir,0)
        print *, "parentDir is ", parentDir

        paramName="eventName"
        call searchForParams(tmpfile,paramName,eventName,0)
        print *, "eventName is ", trim(eventName)

        paramName="stationName"
        call searchForParams(tmpfile,paramName,stationName,0)
        print *, "stationName is ", stationName

        paramName="stlalo"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy, *) stla, stlo
        print *, "stla is ", stla, " and stlo is", stlo

        rlat=stla
        rlon=stlo
    
        paramName="searchAreaDistance"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) gcarcmin,gcarcmax,dgcarc
        print *, "gcarc min, max and interval are: ", gcarcmin,gcarcmax,dgcarc


        ntheta = int((gcarcmax-gcarcmin)/dgcarc)+1
        allocate(gcarc(1:ntheta))
        allocate(ithetaD(1:ntheta))
        do iloop=1,ntheta
            gcarc(iloop) = gcarcmin + dble(iloop-1)*dgcarc
        enddo

        paramName="searchAreaAzimuth"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) azimuthmin,azimuthmax,dazimuth
        print *, "azimuth min, max and interval are: ", azimuthmin,azimuthmax,dazimuth

        nphi = int((azimuthmax-azimuthmin)/dazimuth)+1
        allocate(azimuth(1:nphi))
        do iloop=1,nphi
            azimuth(iloop) = azimuthmin + dble(iloop-1)*dazimuth
        enddo


        paramName="searchAreaRadius"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) radiusmin,radiusmax,dradius
        print *, "source radius min, max and interval are: ", radiusmin,radiusmax,dradius


        nr = int((radiusmax-radiusmin)/dradius)+1
        allocate(radius(1:nr))
        allocate(iradiusD(1:nr))
        do iloop=1,nr
            radius(iloop) = radiusmin + dble(iloop-1)*dradius
        enddo
    
    
        paramName="dt"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) dt
        print *, "dt is ", dt, "(s)"
        samplingHz=1.d0/dt

        !paramName="tlenDSM"
        !call searchForParams(tmpfile,paramName,dummy,1)
        !read(dummy,*) tlenDSM
        !print *, "tlen in DSM synthetic is ", tlenDSM, " (s)"
        !npDSM=int(tlenDSM/dt)
        !print *, "  thus npDSM is ", npDSM

        paramName="tlenData"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) tlenData
        npData = int(tlenData/dt)
        print *, "tlenData is ", tlenData, " (s) and thus npData is ", npData



        paramName="movingWindowStep"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) movingWindowStep
        print *, "movingWindowStep for each MT inversion is ", movingWindowStep, " (s)"
        ntStep=int(movingWindowStep/dt)
        print *, "  thus npStep is ", ntStep

        paramName="npButterworth"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) npButterworth

        paramName="fmin"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) fmin

        paramName="fmax"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) fmax

        print *, "Butterworth parameters (npButterworth, fmin, fmax) are:"
        print *, "  ", npButterworth, fmin, fmax


        paramName="effectiveSynWindow"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) start, end
        print *, "effective DSM synthetic window is from ", start, " (s) to ", end, " (s)"

       
        paramName="numberofSynWindows"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) ntwin


        allocate(twin(1:4,1:ntwin))
        allocate(itwin(1:4,1:ntwin))


        do iloop=1,ntwin
            write(paramName,*) iloop
            paramName="synWindow"//trim(adjustl(paramName))
            
            call searchForParams(tmpfile,paramName,dummy,1)
            read(dummy,*) twin(1,iloop),twin(2,iloop),twin(3,iloop),twin(4,iloop)
            print *, iloop,"-th window in syn is characterised by:", &
                twin(1,iloop),twin(2,iloop),twin(3,iloop),twin(4,iloop)
        enddo

        itwin=int(twin/dt)




        allocate(obsRaw(1:npData,1:3))
        allocate(obsFilt(1:npData,1:3))

        paramName="obsZfile"
        call searchForParams(tmpfile,paramName,obsfile,0)
        obsfile=trim(workingDir)//"/"//trim(obsfile)
        print *, "Reading obsZfile: ", obsfile
        open(unit=10,file=obsfile,status='unknown')
        do it=1,npData
            read(10,*) obsRaw(it,1)
        enddo

        paramName="obsNfile"
        call searchForParams(tmpfile,paramName,obsfile,0)
        obsfile=trim(workingDir)//"/"//trim(obsfile)
        print *, "Reading obsNfile: ", obsfile
        open(unit=10,file=obsfile,status='unknown')
        do it=1,npData
            read(10,*) obsRaw(it,2)
        enddo

        paramName="obsEfile"
        call searchForParams(tmpfile,paramName,obsfile,0)
        obsfile=trim(workingDir)//"/"//trim(obsfile)
        print *, "Reading obsEfile: ", obsfile
        open(unit=10,file=obsfile,status='unknown')
        do it=1,npData
            read(10,*) obsRaw(it,3)
        enddo
        

        paramName="numberofObsWindows"
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) ntwinObs
        

        paramName="obsMovingWindowOption"
        call searchForParams(tmpfile,paramName,dummy,0)
        if(dummy(1:5).eq.'fixed') then
            NmovingWindowDimension=1
        elseif(dummy(1:).eq.'independent') then
            NmovingWindowDimension=ntwinObs
        endif

        print *, "the number of Obs Windows is ", ntwinObs, " and the dimension is ",  &
            NmovingWindowDimension, " since you chose ", trim(dummy), " option."



        allocate(twinObs(1:4,1:ntwinObs))
        allocate(itwinObs(1:4,1:ntwinObs))


        do iloop=1,ntwinObs
            write(paramName,*) iloop
            paramName="obsWindow"//trim(adjustl(paramName))
            
            call searchForParams(tmpfile,paramName,dummy,1)
            read(dummy,*) twinObs(1,iloop),twinObs(2,iloop),twinObs(3,iloop),twinObs(4,iloop)
            print *, iloop,"-th window in obs is characterised by:", &
                twinObs(1,iloop),twinObs(2,iloop),twinObs(3,iloop),twinObs(4,iloop)
        enddo
        itwinObs=int(twinObs/dt)


        
        !commandline = 'mkdir -p '//trim(parentDir)
        !call system(commandline)
        call pinputDSM(DSMconfFile,PoutputDir,psvmodel,&
                modelname,tlenFull,rmin_,rmax_,rdelta_, &
                r0min,r0max,r0delta,thetamin,thetamax,thetadelta, &
                imin,imax,rsgtswitch,tsgtswitch, &
                synnswitch,SGTinfo)
        call readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
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

        
        call makingIndependentWindow

        if(tlenFull<end) then
            print *, "DSM tlen is ", tlenFull, "and you asked to take window up to ", end
            print *, "it is impossible, sorry."
            stop
        endif

        tmpfile='tmpReadPSVmodel'
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

        print *, "r, theta are discretised in ", r_n, theta_n, "pieces in DSM pre-calculation"


        ! check r(:) and r_(:) because we don't interpolate for depth
        iradiusD=0
        do iloop=1,nr
            icheck=0
            do jloop=1,r_n
                if(dabs(r_(jloop)-radius(iloop))< 1.d-5*radius(iloop)) then
                    icheck=1
                    iradiusD(iloop)=jloop
                endif
            enddo
            if(icheck.eq.0) then
                print *, radius(iloop), "is not in the catalogue, sorry"
                stop
            endif
        enddo
        
        ! check gcarc(:) and thetaD(:) because we don't interpolate for depth
        ithetaD=0
        do iloop=1,ntheta
            icheck=0
            do jloop=1,theta_n
                if(dabs(thetaD(jloop)-gcarc(iloop))<1.d-5*gcarc(iloop)) then
                    icheck=1
                    ithetaD(iloop)=jloop
                endif
            enddo
            if(icheck.eq.0) then
                print *, gcarc(iloop), "is not in the catalogue, sorry"
                stop
            endif
        enddo
     
        allocate(latgeo(1:nphi,1:ntheta),longeo(1:nphi,1:ntheta))
        allocate(crq(1:nphi,1:ntheta),crq2(1:nphi,1:ntheta))
        allocate(srq(1:nphi,1:ntheta),srq2(1:nphi,1:ntheta))

        latgeo=0.d0
        longeo=0.d0
    
        do iloop=1,ntheta
            do jloop=1,nphi
                call geoCoordinates(stla*degree2radian,stlo*degree2radian, &
                    latgeo(jloop,iloop),longeo(jloop,iloop),azimuth(jloop)*degree2radian, &
                    fdummy,gcarc(iloop)*degree2radian) !!! geoCooridnates uses radian
            enddo
        enddo

        latgeo=latgeo*radian2degree
        longeo=longeo*radian2degree

        print *, "InSight: ", stla, stlo

        do iloop=1,ntheta
            do jloop=1,nphi
                print *, gcarc(iloop),azimuth(jloop), latgeo(jloop,iloop),longeo(jloop,iloop)
                call azimth(0,latgeo(jloop,iloop),longeo(jloop,iloop),stla,stlo,fdummy2,phirq,fdummy)
                !print *, " inverse: ",fdummy2,phirq,fdummy


                phirq=pi-phirq*degree2radian

                if(phirq.lt.0.d0) phirq=phirq+2.d0*pi
                if(phirq.gt.(2.d0*pi)) phirq=phirq-2.d0*pi
                
                crq(jloop,iloop)=dcos(phirq)
                srq(jloop,iloop)=dsin(phirq)
                crq2(jloop,iloop)=dcos(2.d0*phirq)
                srq2(jloop,iloop)=dsin(2.d0*phirq)
                
            enddo
        enddo


        



       



        ! allocate SGTs, synthetics in frequency
        allocate(omega(imin:imax))
        do iloop = imin, imax
            omega(iloop) = 2.d0*pi*dble(iloop)/tlenFull
        enddo
     
        !nConfiguration=r_n*theta_n*phi_n

        nConfiguration=nr*ntheta*nphi
        print *, "source locations to be tested (in r, theta, phi, total): "
        print *, "      ", nr,ntheta,nphi,nConfiguration
    elseif((dummy(1:6).eq.'normal').or.(dummy(1:4).eq.'test')) then


        print *, "no more support for non RSGT methods. Sorry."
        stop
     
    endif
  
  

end subroutine pinput




subroutine azimth(ellips,slat,slon,rlat,rlon,delta,azim,bazim)

!   This routine uses Euler angles to find the geocentric distance,
!   azimuth, and back azimuth for a source-reciever pair.
!
!   Input
!
!     slat  - source geographic latitude in decimal degrees
!     slon  - source longitude in decimal degrees
!     rlat  - reciever geographic latitude in decimal degrees
!     rlon  - reciever longitude in decimal degrees
!
!   Output
!
!     delta - geocentric source-reciever distance in decimal degrees of arc
!     azim  - geocentric azimuth from the source to the reciever
!     bazim - geocentric back azimuth from the reciever to the source
!
!   The distance calculated here delta is always between 0 and 180 degrees.
!   Accordingly, the azimuth and back azimuth are defined for the minor
!   arc between (slat,slon) and (rlat,rlon).
!
!   if ellips = 0 then geocentric = geocentric
!     because in NF version it is already taken into account so ellips should be 0 always


    implicit none
    real(kind(0d0)) :: dtor,e,slatra,slat,w,s,scolat,rlatra,rlat,rcolat,slonra,rlon,c2,s2,c1,s1,slatrc,x0,y0,z0,x1,y1,z1,x2,y2,z2,slon,rlonra,delta,azim,bazim,pi
    real(kind(0d0)), parameter :: flt = 298.25d0
    integer :: ellips


    dtor=4.d0*datan(1.d0)/180.d0
    pi=4.d0*datan(1.d0)

    if(ellips.ne.0) then
        e=1.d0/flt
    else
        e=0.d0
    endif


    !   Convert to geocentric coordinates and from latitude to colatitude.

    slatra=dtor*slat
    w=dsin(slatra)
    s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(slatra)
    !scolat=1.5707963d0-slatra+s
    scolat=pi*5.d-1-slatra+s
    rlatra=dtor*rlat
    w=dsin(rlatra)
    s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(rlatra)
    !rcolat=1.5707963d0-rlatra+s
    rcolat=pi*5.d-1-rlatra+s

    slonra=slon*dtor
    rlonra=rlon*dtor
    c2=dcos(scolat)
    s2=dsin(scolat)
    c1=dcos(slonra)
    s1=dsin(slonra)
    slatrc=dsin(rcolat)




    !  Find the azimuth and distance by rotating the source to the north pole.

    x0=slatrc*dcos(rlonra)
    y0=slatrc*dsin(rlonra)
    z0=dcos(rcolat)
    x1=c1*x0+s1*y0

    z0=dcos(rcolat)
    x1=c1*x0+s1*y0
    y1=-s1*x0+c1*y0
    z1=z0
    x2=c2*x1-s2*z1
    y2=y1
    z2=c2*z1+s2*x1
    call angles(x2,y2,z2,delta,azim)
    azim=180.d0-azim

    !  Find the back azimuth by rotating the receiver to the north pole.

    c2=dcos(rcolat)
    s2=dsin(rcolat)
    c1=dcos(rlonra)
    s1=dsin(rlonra)
    slatrc=dsin(scolat)
    x0=slatrc*dcos(slonra)
    y0=slatrc*dsin(slonra)
    z0=dcos(scolat)
    x1=c1*x0+s1*y0
    y1=-s1*x0+c1*y0
    z1=z0
    x2=c2*x1-s2*z1
    y2=y1
    z2=c2*z1+s2*x1
    call angles(x2,y2,z2,delta,bazim)
    bazim=180.d0-bazim

    return
end subroutine azimth

subroutine angles(x,y,z,theta,phi)

    !   Finds the angles theta and phi of a spherical polar coordinate
    !   system from the cartesion coordinates x, y, and z.

    implicit none
    real(kind(0d0)) :: pi,rtod,arg1,x,y,theta,phi,z
    ! real(kind(0d0)), parameter :: eps = 1.d-14
    real(kind(0d0)), parameter :: eps = 0.d0


    pi=4.d0*datan(1.d0)

    rtod=180.d0/pi
    arg1=dsqrt(x*x+y*y)
    theta=datan2(arg1,z)
    if(dabs(x).le.eps.and.dabs(y).le.eps) then
        phi=0.d0
    else
        phi=datan2(y,x)
    endif
    phi=phi*rtod
    theta=theta*rtod

    return
end subroutine angles


subroutine geoCoordinates(glat1,glon1,glat2,glon2,faz,baz,s)
    !
    ! *** solution of the geodetic direct problem after t.vincenty
    ! *** modified rainsford's method with helmert's elliptical terms
    ! *** effective in any azimuth and at any distance short of antipodal
    !
    ! *** a is the semi-major axis of the reference ellipsoid
    ! *** f is the flattening of the reference ellipsoid
    ! *** latitudes and longitudes in radians positive north and east
    ! *** azimuths in radians clockwise from north
    ! *** geodesic distance s assumed in units of semi-major axis a
    !
    ! *** programmed for cdc-6600 by lcdr l.pfeifer ngs rockville md 20feb75
    ! *** modified for system 360 by john g gergen ngs rockville md 750608
    !
    ! Here everything is in radian!
    implicit real*8 (a-h,o-z)
    !common/const/pi,rad
    !common/elipsoid/a,f
    a=1
    pi=4.d0*datan(1.d0)
    f=0
    data eps/0.5d-13/
    r=1.-f
    tu=r*dsin(glat1)/dcos(glat1)
    sf=dsin(faz)
    cf=dcos(faz)
    baz=0.
    if(abs(cf)>eps) baz=datan2(tu,cf)*2.
    cu=1./dsqrt(tu*tu+1.)
    su=tu*cu
    sa=cu*sf
    c2a=-sa*sa+1.
    x=dsqrt((1./r/r-1.)*c2a+1.)+1.
    x=(x-2.)/x
    c=1.-x
    c=(x*x/4.+1)/c
    d=(0.375d0*x*x-1.)*x
    tu=s/r/a/c
    y=tu
100 sy=dsin(y)
    cy=dcos(y)
    cz=dcos(baz+y)
    e=cz*cz*2.-1.
    c=y
    x=e*cy
    y=e+e-1.
    y=(((sy*sy*4.-3.)*y*cz*d/6.+x)*d/4.-cz)*sy*d+tu
    if(dabs(y-c).gt.eps)go to 100
    baz=cu*cy*cf-su*sy
    c=r*dsqrt(sa*sa+baz*baz)
    d=su*cy+cu*sy*cf
    glat2=datan2(d,c)
    c=cu*cy-su*sy*cf
    x=datan2(sy*sf,c)
    c=((-3.*c2a+4.)*f+4.)*c2a*f/16.
    d=((e*cy*c+cz)*sy*c+y)*sa
    glon2=glon1+x-(1.-c)*d*f
    baz=datan2(sa,baz)+pi
    return
end subroutine




subroutine searchForParams(filename,ParamName,textParam,paramisText)
    implicit none
    character(200) :: filename,textParam,text_line
    integer :: paramLength,textLength,paramisText,io
    character(200) :: ParamName
    integer :: iFind, jtemp, iCut
    filename=trim(filename)
    ParamName=trim(ParamName)
    paramLength=len_trim(ParamName)
    !print *, paramLength, ParamName
    iFind=0
    iCut=0
    open(20,file=filename,status='unknown')
    do while(iFind.eq.0)
        read(20,'(a)',IOSTAT=io) text_line
        if(io>0) then
            print *, "oh, no"
            print *, trim(ParamName), " is not found."
            stop
        endif
        textLength=len_trim(text_line)
        !print *, text_line(1:textLength)
        if(text_line(1:paramLength).eq.ParamName(1:paramLength)) then
            do jtemp = 1,textLength
                if(text_line(jtemp:jtemp).eq.'=') iCut = jtemp
            enddo
            !print *, iCut, textLength
            iFind=1
            !print *,text_line(iCut+1:textLength)
            textParam=text_line(iCut+1:textLength)
            if(paramisText.eq.0) then
                textParam=trim(textParam)
                textParam=adjustl(textParam)
            endif
            !print *, textParam
        endif
    enddo
    close(20)
    if(iFind.eq.0) then
        print *, ParamName, "is not found."
        stop
    endif

end subroutine


subroutine makingIndependentWindow
    use parameters
    implicit none
    integer :: iloop,jloop,tmpinteger
    integer, allocatable :: indexInWindow(:)
    character(200) ::  dummy, paramName, tmpfile

    ! Moving windows

    allocate(fMovingWindowStart(1:NmovingWindowDimension))
    allocate(fMovingWindowEnd(1:NmovingWindowDimension))

    allocate(iMovingWindowStart(1:NmovingWindowDimension))
    allocate(iMovingWindowEnd(1:NmovingWindowDimension))

    allocate(indexInWindow(1:NmovingWindowDimension))

    allocate(totalNumberInWindowDimension(1:NmovingWindowDimension))

    tmpfile='tmpfileMarsInversion'

    do iloop=1,NmovingWindowDimension
        write(paramName,*) iloop
        paramName="movingWindowRange"//trim(adjustl(paramName))
        
        call searchForParams(tmpfile,paramName,dummy,1)
        read(dummy,*) fMovingWindowStart(iloop), fMovingWindowEnd(iloop)
        print *, iloop,"-th moving range in obs is characterised by:", &
            fMovingWindowStart(iloop), fMovingWindowEnd(iloop)
    enddo

    iMovingWindowStart=int(fMovingWindowStart/dt)
    iMovingWindowEnd=int(fMovingWindowEnd/dt)
    

    nTimeCombination = 1
    do iloop=1,NmovingWindowDimension
        iMovingWindowStart(iloop)=itwinObs(1,iloop)+iMovingWindowStart(iloop)
        iMovingWindowEnd(iloop)=itwinObs(1,iloop)+iMovingWindowEnd(iloop)
        ! check whether syn and obs are available for these indices
        if(iMovingWindowEnd(iloop)<iWindowStart) then
            print *, "no sufficient data points in syn data for the window", iloop
            stop
        endif
        !if(itwinObs(1,iloop)<1) then
        !    print *, "no sufficient data points in obs data for the window", iloop
        !    stop
        !endif
        if(iMovingWindowStart(iloop)>iWindowEnd) then
            print *, "no sufficient data points in syn data for the window", iloop
            stop
        endif
        !if(itwinObs(4,iloop)+iMovingWindowEnd(iloop)>npData) then
        !    print *, "no sufficient data points in obs data for the window", iloop
        !    stop
        !endif
        totalNumberInWindowDimension(iloop)=(iMovingWindowEnd(iloop)-iMovingWindowStart(iloop))/ntStep+1
        nTimeCombination = nTimeCombination*totalNumberInWindowDimension(iloop)
        !print *, iloop, totalNumberInWindowDimension(iloop)
    enddo

    allocate(iEachWindowStart(1:nTimeCombination,1:NmovingWindowDimension))
    allocate(iEachWindowEnd(1:nTimeCombination,1:NmovingWindowDimension))

    do jloop=1,nTimeCombination
        tmpinteger=jloop-1
        do iloop=1,NmovingWindowDimension
            indexInWindow(iloop)=mod(tmpinteger,totalNumberInWindowDimension(iloop))+1

            tmpinteger=tmpinteger/totalNumberInWindowDimension(iloop)
            iEachWindowStart(jloop,iloop)=iMovingWindowStart(iloop)+ntStep*(indexInWindow(iloop)-1)
            iEachWindowEnd(jloop,iloop)=iMovingWindowEnd(iloop)+ntStep*(indexInWindow(iloop)-1)
            
        enddo
        !print *, jloop,indexInWindow(:),iEachWindowStart(jloop,:),iEachWindowEnd(jloop,:)
    enddo


end subroutine


subroutine pinputDSM(DSMconfFile,outputDir,psvmodel, &
    modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta, &
    thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch,SGTinfo)
  implicit none
  !character(120), parameter :: tmpfile='tmpworkingfile_for_SGTcalcul'
  character(200) :: dummy,outputDir,psvmodel,modelname,DSMconfFile,SGTinfo
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
  modelname=adjustl(modelname)
  outputDir=adjustl(outputDir)
  psvmodel=adjustl(psvmodel)
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
