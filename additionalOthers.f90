!
!   additionalOthers.f90
!   MarsMTInversion
!
!   Created by fuji on 26/03/2020.
!   Copyright 2020 nfuji. All rights reserved.
!


subroutine preprocessing

    use parameters
    use tmpSGTs
    use angles
    use mainparameters
    ! making taper for synthetics
    taperDSM=0.d0
    do iWindow=1,ntwin
        do it=itwin(1,iWindow),itwin(4,iWindow)
            if((it.gt.itwin(1,iWindow)).and.(it.lt.itwin(2,iWindow))) then
                xfwin=dsin(0.5d0*pi*dble(it-itwin(1,iWindow))/dble(itwin(2,iWindow)-itwin(1,iWindow)))
                taperDSM(it)=xfwin*xfwin
            elseif((it.ge.itwin(2,iWindow)).and.(it.le.itwin(3,iWindow))) then
                taperDSM(it)=1.d0
            elseif((it.gt.itwin(3,iWindow)).and.(it.lt.itwin(4,iWindow))) then
                xfwin=dsin(0.5d0*pi*dble(it-itwin(4,iWindow))/dble(itwin(3,iWindow)-itwin(4,iWindow)))
                taperDSM(it)=xfwin*xfwin
            endif
        enddo
    enddo


    ! For the observed data we taper the whole signal of a length of npData
    
    ! A priori taper for stabilising bwfilt
    taperOBS=1.d0
    do it=0,npData/20
        xfwin=dsin(0.5d0*pi*dble(it-1)/dble(npData/20-1))
        taperOBS(it)=xfwin*xfwin
    enddo
    do it=npData-npData/20,npData
        xfwin=dsin(0.5d0*pi*dble(it-npData)/dble(-npData/20))
        taperOBS(it)=xfwin*xfwin
    enddo
    
    ! For the observed data we filter the whole signal of a length of npData

    do icomp=1,3
        obsRaw(0:npData,icomp)=obsRaw(0:npData,icomp)*taperOBS(0:npData)
       call bwfilt(obsRaw(0:npData,icomp),obsFilt(0:npData,icomp),dt,npData+1,1,npButterworth,fmin,fmax)
    enddo
! making taper for observed data
   taperOBS=0.d0
   do iWindow=1,ntwinObs
       do it=itwinObs(1,iWindow),itwinObs(4,iWindow)
           if((it.gt.itwinObs(1,iWindow)).and.(it.lt.itwinObs(2,iWindow))) then
               xfwin=dsin(0.5d0*pi*dble(it-itwinObs(1,iWindow))/dble(itwinObs(2,iWindow)-itwinObs(1,iWindow)))
               taperOBS(it)=xfwin*xfwin
           elseif((it.ge.itwinObs(2,iWindow)).and.(it.le.itwinObs(3,iWindow))) then
               taperOBS(it)=1.d0
           elseif((it.gt.itwinObs(3,iWindow)).and.(it.lt.itwinObs(4,iWindow))) then
               xfwin=dsin(0.5d0*pi*dble(it-itwinObs(4,iWindow))/dble(itwinObs(3,iWindow)-itwinObs(4,iWindow)))
               taperOBS(it)=xfwin*xfwin
           endif
       enddo
   enddo

    ! we taper the filtered OBS

    ! For the observed data we filter the whole signal of a length of npData

    do icomp=1,3
        obsFiltTapered(0:npData,icomp)=obsFilt(0:npData,icomp)*taperOBS(0:npData)
    enddo
   

    obsArray=obsFiltTapered
    obsRawArray=obsFilt

    mtInverted=0.d0

end subroutine preprocessing


subroutine obsFiltWriting
    use parameters
    use tmpSGTs
    use angles
    use mainparameters
    ! Raw data and filtered data are written as fort.11-13 for references
       
    open(21,file="obsTrated.dat",status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)
    open(22,file="obsTratedFiltered.dat",status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)
    do it=0,npData
        tmpfloat(1)=sngl(dt*dble(it))
        tmpfloat(2)=sngl(obsRaw(it,1))
        tmpfloat(3)=sngl(obsRaw(it,2))
        tmpfloat(4)=sngl(obsRaw(it,3))
        write(21,rec=it+1) tmpfloat(1:4)
        tmpfloat(1)=sngl(dt*dble(it))
        tmpfloat(2)=sngl(obsFilt(it,1))
        tmpfloat(3)=sngl(obsFilt(it,2))
        tmpfloat(4)=sngl(obsFilt(it,3))
        write(22,rec=it+1) tmpfloat(1:4)
    enddo
    close(21)
    close(22)
end subroutine obsFiltWriting


