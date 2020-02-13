!
!   sgt.f90
!   MarsInversion
!
!   Created by fuji on 04/12/2019.
!   Copyright 2019 nfuji. All rights reserved.
!
subroutine rdsgtomega(rx,num_sgt,num_psv,ipsvorsh)

  use parameters
  use tmpSGTs
  implicit none
  real(kind(0d0)) :: rx ! in order to read the catalogue
  integer :: i,j,num_sgt,ipsvorsh,num_psv
  complex(kind(0e0)) :: sgtsngl(1:num_sgt,1:theta_n)
  complex(kind(0d0)) :: sgtdouble(1:num_sgt,1:theta_n)
  character(200) :: coutfile
  integer :: itheta
  ! about ipsvorsh
  ! PSV synn =  2; SH synn =  1
  ! PSV rsgt = 20; SH rsgt = 10
  ! PSV tsgt =200; SH tsgt =100

  if((ipsvorsh.eq.2).or.(ipsvorsh.eq.20).or.(ipsvorsh.eq.200)) then
     if(num_sgt.ne.num_sgt) then
        print *, "the number of SGTs is not propre."
        stop
     endif
  endif
  

  do i = imin,imax
     sgtsngl = cmplx(0.e0)
     sgtdouble = dcmplx(0.d0)
     if(i.ne.0) then
        !if(ipsvorsh.eq.2) then
        !   write(coutfile, '(I7,".",I7,".SYNN_PSV")') int(rx*1.d3),i

           ! just for depth600 because of BUGG in SGTpsv!!!!
           !write(coutfile, '(I7,".",I7,".SYNN_PSV")') 0,i
           !!!!!! ATTENTION !!!! PLEASE THIS IS NOT PROPRE!!!!!
           
        !elseif(ipsvorsh.eq.1) then
        !   write(coutfile, '(I7,".",I7,".SYNN_SH")') int(rx*1.d3),i
        !else
        if(ipsvorsh.eq.20) then
           write(coutfile, '(I7,".",I7,".RSGT_PSV")') int(rx*1.d3),i
        elseif(ipsvorsh.eq.10) then
           write(coutfile, '(I7,".",I7,".RSGT_SH")') int(rx*1.d3),i
        !elseif(ipsvorsh.eq.200) then
        !   write(coutfile, '(I7,".",I7,".",I7,".TSGT_PSV")') int(rx*1.d3),int(ry*1.d3),i
        !elseif(ipsvorsh.eq.100) then
        !   write(coutfile, '(I7,".",I7,".",I7,".TSGT_SH")') int(rx*1.d3),int(ry*1.d3),i
        endif

        
        !if((ipsvorsh.eq.200).or.(ipsvorsh.eq.100)) then ! TSGT
        !   do j = 1,29
        !         if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
        !   enddo
        !   coutfile = trim(modelname)//"."//coutfile
        !   coutfile = trim(PoutputDir)//"/TSGT/"//coutfile
        !else
        do j = 1,21
            if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
        enddo

        coutfile = trim(modelname)//"."//coutfile

        coutfile = trim(PoutputDir)//"/RSGT/"//coutfile ! RSGT & SYNN
        !endif

        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*num_sgt*kind(0e0)*theta_n)
        read(1,rec=1)sgtsngl(1:num_sgt,1:theta_n)
        
    
        sgtdouble(1:num_sgt,1:theta_n) = sgtsngl(1:num_sgt,1:theta_n)
        
        !open(1,file=coutfile,status='old',form='unformatted',action='read')
        !read(1) sgtsngl(1:num_sgt,1:theta_n)
        !sgtdouble(1:num_sgt,1:theta_n)=sgtsngl(1:num_sgt,1:theta_n)


        ! h9 in SH is 0 !
        if((ipsvorsh.eq.1).or.(ipsvorsh.eq.10)) then
            !sgtdouble(4,1:theta_n) = cmplx(0.d0)
           !sgtdouble(3:5,1:theta_n) = sgtdouble(3:5,1:theta_n)*5.d-1
           !sgtdouble(3:5,1:theta_n) = -sgtdouble(3:5,1:theta_n)
           
        endif
        
        if(ipsvorsh.eq.1) then
           sgtdouble(3:5,1:theta_n) = sgtdouble(3:5,1:theta_n)*5.d-1
           !sgtdouble(1:2,1:theta_n) = cmplx(0.d0)
        endif

        if(ipsvorsh.eq.2) then
           sgtdouble(2,1:theta_n) = sgtdouble(2,1:theta_n)*5.d-1
           sgtdouble(8:10,1:theta_n) = sgtdouble(8:10,1:theta_n)*5.d-1
        endif



        close(1)
            
        do itheta = 1,theta_n
           !redtime = -thetaD(itheta)*c_red_reci
           if(ipsvorsh.ge.10) then
              sgtdouble(1:num_sgt,itheta) = sgtdouble(1:num_sgt,itheta) !&
            !       *exp(cmplx(0.d0,-2.d0*pi*dble(i)/tlen*redtime))*1.d-21
           else
              sgtdouble(1:num_sgt,itheta) = sgtdouble(1:num_sgt,itheta) !&
            !       *exp(cmplx(0.d0,-2.d0*pi*dble(i)/tlen*redtime))*1.d-21
           endif
                 
        enddo
        !rsgtomega(1:1,i,1:theta_n) = rsgtomega(1,i,1:theta_n) + rsgtdouble(1,1:theta_n)
        !if(ipsvorsh.eq.2) then
        !   synnomega(1:num_psv,i,1:theta_n)=synnomega(1:num_psv,i,1:theta_n)+sgtdouble(1:num_psv,1:theta_n)
        !else
        if(ipsvorsh.eq.20) then
           rsgtomega(1:num_psv,i,1:theta_n)=rsgtomega(1:num_psv,i,1:theta_n)+sgtdouble(1:num_psv,1:theta_n)
        !elseif(ipsvorsh.eq.200) then
        !   tsgtomega(1:num_psv,i,1:theta_n)=tsgtomega(1:num_psv,i,1:theta_n)+sgtdouble(1:num_psv,1:theta_n)
        !elseif(ipsvorsh.eq.1) then
        !   synnomega(6:10,i,1:theta_n)=synnomega(6:10,i,1:theta_n)+sgtdouble(1:5,1:theta_n)
        elseif(ipsvorsh.eq.10) then
          rsgtomega(6:10,i,1:theta_n)=rsgtomega(6:10,i,1:theta_n)+sgtdouble(1:5,1:theta_n)
        !elseif(ipsvorsh.eq.100) then
        !   tsgtomega(7,i,1:theta_n) = tsgtomega(7,i,1:theta_n) + sgtdouble(1,1:theta_n)
        !   tsgtomega(12:20,i,1:theta_n)= tsgtomega(12:20,i,1:theta_n)+ sgtdouble(2:10,1:theta_n)
        endif
     endif

  enddo
 
  return
end subroutine rdsgtomega

subroutine tensorFFT_double(n,imin,imax,np1,ccvec,rvec,omegai,tlen,iWindowStart,iWindowEnd)
  ! this subroutine particularly calculates the FFT of the given tensor and make a double tensor
  
  implicit none
  integer :: iWindowStart,iWindowEnd
  integer :: i,j,n,imin,imax,np1,n1,m1
  complex(kind(0d0)) :: ccvec(1:n,imin:imax)
  complex(kind(0d0)) :: cvec(0:2*np1-1,1:n),eachcvec(0:2*np1-1)
  real(kind(0d0)) :: rvec(iWindowStart:iWindowEnd,1:n)
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: omegai, tlen,samplingHz
  real(kind(0d0)) :: omegai_log
  cvec = dcmplx(0.d0)
  cvec(imin:imax,1:n)=transpose(ccvec(1:n,imin:imax))

  rvec=0.d0

    !print *, "this is FFT", np1, imin,imax, omegai,np1-(np1-imax),np1-(np1-imin),dble(imax)/dble(2*np1),dble(2*np1)/dble(imax)
    samplingHz = dble(2*np1)/tlen
  
    omegai_log=-dlog(omegai)/tlen
  do j = 1,n
     do i = np1-imax, np1-imin
        n1 = np1 +i
        m1 = np1 -i
        !cvec(n1,j) = conjg(cvec(m1,j))
        !cvec(m1,j)=dconjg(cvec(n1,j))
        cvec(n1,j)=dconjg(cvec(m1,j))
     enddo
  enddo
  
  
  
  
  do j = 1,n
     eachcvec=dcmplx(0.d0)
     eachcvec(0:2*np1-1)=cvec(0:2*np1-1,j)

    
   

     call cdft(4*np1,dcos(pi/dble(2*np1)),dsin(pi/dble(2*np1)), eachcvec(0:2*np1-1))
     do i = iWindowStart, iWindowEnd
        rvec(i,j) = dble(eachcvec(i))*dexp(omegai_log*dble(i)/samplingHz)/tlen*1.d3
        !rvec(i,j) = dble(eachcvec(i))*dble(exp(omegai*tlen*dble(i)/dble(np1)))/tlen*1.d3
        !rvec(i,j)= dble(exp(omegai*dble(i)/samplingHz))/tlen*1.d3
      
     enddo
  enddo
  
  
  return
end subroutine tensorFFT_double


subroutine rsgt2h3time_authentic(ip,ith)
    use parameters
    use angles
    use tmpSGTs
    !use kernels
    implicit none
    integer ip,ith
    !   This subroutine calculates the 18 independent elements of the 3rd-order
    !   SSGT from the 10 azimuth independent coefficients and the sines and
    !   cosines of the azimuths.

! during RSGT calculation I did some misakes: phi_RQ =pi instead of 0



    !   Vertical component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.
        tmparray(iWindowStart:iWindowEnd,1,1)=rsgtTime(iWindowStart:iWindowEnd,1)
        tmparray(iWindowStart:iWindowEnd,1,2)=rsgtTime(iWindowStart:iWindowEnd,3)*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,4)
        tmparray(iWindowStart:iWindowEnd,1,3)=-rsgtTime(iWindowStart:iWindowEnd,3)*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,4)
        tmparray(iWindowStart:iWindowEnd,1,4)=-rsgtTime(iWindowStart:iWindowEnd,2)*crq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,1,5)=-rsgtTime(iWindowStart:iWindowEnd,2)*srq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,1,6)=rsgtTime(iWindowStart:iWindowEnd,3)*srq2(ip,ith)

     !   Radial component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.

        tmparray(iWindowStart:iWindowEnd,2,1)=rsgtTime(iWindowStart:iWindowEnd,5)
        tmparray(iWindowStart:iWindowEnd,2,2)=-(rsgtTime(iWindowStart:iWindowEnd,7))*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,6)
        tmparray(iWindowStart:iWindowEnd,2,3)=(rsgtTime(iWindowStart:iWindowEnd,7))*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,6)
        tmparray(iWindowStart:iWindowEnd,2,4)=(rsgtTime(iWindowStart:iWindowEnd,9))*crq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,2,5)=(rsgtTime(iWindowStart:iWindowEnd,9))*srq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,2,6)=-(rsgtTime(iWindowStart:iWindowEnd,7))*srq2(ip,ith)

     !   Transverse component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.
 
        tmparray(iWindowStart:iWindowEnd,3,1)=0.d0
        tmparray(iWindowStart:iWindowEnd,3,2)=(rsgtTime(iWindowStart:iWindowEnd,8))*srq2(ip,ith)
        tmparray(iWindowStart:iWindowEnd,3,3)=-tmparray(iWindowStart:iWindowEnd,3,3)
        tmparray(iWindowStart:iWindowEnd,3,4)=-(rsgtTime(iWindowStart:iWindowEnd,10))*srq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,3,5)=(rsgtTime(iWindowStart:iWindowEnd,10))*crq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,3,6)=-(rsgtTime(iWindowStart:iWindowEnd,8))*crq2(ip,ith)

  return
end subroutine rsgt2h3time_authentic

subroutine rsgt2h3time_adhoc(ip,ith)
    use parameters
    use angles
    use tmpSGTs
    !use kernels
    implicit none
    integer ip,ith
    !   This subroutine calculates the 18 independent elements of the 3rd-order
    !   SSGT from the 10 azimuth independent coefficients and the sines and
    !   cosines of the azimuths.



    !   Vertical component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp. Perfectly modified
        tmparray(iWindowStart:iWindowEnd,1,1)=rsgtTime(iWindowStart:iWindowEnd,1) ! valid
        tmparray(iWindowStart:iWindowEnd,1,2)=rsgtTime(iWindowStart:iWindowEnd,3)*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,4)
        tmparray(iWindowStart:iWindowEnd,1,3)=-rsgtTime(iWindowStart:iWindowEnd,3)*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,4)
        tmparray(iWindowStart:iWindowEnd,1,4)=2.d0*rsgtTime(iWindowStart:iWindowEnd,2)*crq(ip,ith) !modified
        tmparray(iWindowStart:iWindowEnd,1,5)=2.d0*rsgtTime(iWindowStart:iWindowEnd,2)*srq(ip,ith)
            !modified
        tmparray(iWindowStart:iWindowEnd,1,6)=2.d0*rsgtTime(iWindowStart:iWindowEnd,3)*srq2(ip,ith) ! modified

     !   Radial component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.

        tmparray(iWindowStart:iWindowEnd,2,1)=rsgtTime(iWindowStart:iWindowEnd,5)
        tmparray(iWindowStart:iWindowEnd,2,2)=-(rsgtTime(iWindowStart:iWindowEnd,7))*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,6)
        tmparray(iWindowStart:iWindowEnd,2,3)=(rsgtTime(iWindowStart:iWindowEnd,7))*crq2(ip,ith)-rsgtTime(iWindowStart:iWindowEnd,6)
        tmparray(iWindowStart:iWindowEnd,2,4)=(rsgtTime(iWindowStart:iWindowEnd,9))*crq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,2,5)=(rsgtTime(iWindowStart:iWindowEnd,9))*srq(ip,ith)
        tmparray(iWindowStart:iWindowEnd,2,6)=-(rsgtTime(iWindowStart:iWindowEnd,7))*srq2(ip,ith)

     !   Transverse component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp. Perfectly modified
 
        tmparray(iWindowStart:iWindowEnd,3,1)=0.d0
        tmparray(iWindowStart:iWindowEnd,3,2)=(rsgtTime(iWindowStart:iWindowEnd,8))*srq2(ip,ith)
        tmparray(iWindowStart:iWindowEnd,3,3)=-tmparray(iWindowStart:iWindowEnd,3,2)
        tmparray(iWindowStart:iWindowEnd,3,4)=2.d0*(rsgtTime(iWindowStart:iWindowEnd,10))*srq(ip,ith) !modified
        tmparray(iWindowStart:iWindowEnd,3,5)=-2.d0*(rsgtTime(iWindowStart:iWindowEnd,10))*crq(ip,ith) ! modified
        tmparray(iWindowStart:iWindowEnd,3,6)=-2.d0*(rsgtTime(iWindowStart:iWindowEnd,8))*crq2(ip,ith) ! modified

  return
end subroutine rsgt2h3time_adhoc



subroutine vectorFFT_double(imin,imax,np1,ccvec,rvec,omegai,tlen,iWindowStart,iWindowEnd)
  ! this subroutine particularly calculates the FFT of the given tensor and make a double tensor
  
  implicit none
  integer :: iWindowStart,iWindowEnd
  integer :: i,imin,imax,np1,n1,m1
  complex(kind(0d0)) :: ccvec(imin:imax)
  complex(kind(0d0)) :: cvec(0:2*np1-1)
  real(kind(0d0)) :: rvec(iWindowStart:iWindowEnd)
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: omegai, tlen,samplingHz
  
  cvec = dcmplx(0.d0)
  cvec(imin:imax)=ccvec(imin:imax)

  samplingHz = dble(2*np1)/tlen
  
  do i = imin, np1-1
     n1 = np1 +i
     m1 = np1 -i
     cvec(n1) = dconjg(cvec(m1))
  enddo
  

  call cdft(4*np1,dcos(pi/(2*np1)),dsin(pi/(2*np1)), cvec(0:2*np1-1))
  do i = iWindowStart, iWindowEnd
     rvec(i) = dble(dble(cvec(i))*dble(dexp(omegai*dble(i)/samplingHz))/tlen*1.d3)
  enddo

  
  
  return
end subroutine vectorFFT_double

subroutine cdft(n, wr, wi, c)
      
  integer :: n, i, j, k, l, m
  real(kind(0d0)) :: wr, wi, a(0 : n - 1), wmr, wmi, wkr, wki
  real(kind(0d0)) ::wdr, wdi, ss, xr, xi
  complex(kind(0d0)) :: c(0:n/2-1)
  
  do i = 0, n/2-1
     a(2*i) = dble(c(i))
     a(2*i+1) = imag(c(i))
  enddo
  

  wmr = wr
  wmi = wi
  m = n
  do while (m .gt. 4)
     l = m / 2
     wkr = 1
     wki = 0
     wdr = 1 - 2 * wmi * wmi
     wdi = 2 * wmi * wmr
     ss = 2 * wdi
     wmr = wdr
     wmi = wdi
     do j = 0, n - m, m
        i = j + l
        xr = a(j) - a(i)
        xi = a(j + 1) - a(i + 1)
        a(j) = a(j) + a(i)
        a(j + 1) = a(j + 1) + a(i + 1)
        a(i) = xr
        a(i + 1) = xi
        xr = a(j + 2) - a(i + 2)
        xi = a(j + 3) - a(i + 3)
        a(j + 2) = a(j + 2) + a(i + 2)
        a(j + 3) = a(j + 3) + a(i + 3)
        a(i + 2) = wdr * xr - wdi * xi
        a(i + 3) = wdr * xi + wdi * xr
     enddo
     do k = 4, l - 4, 4
        wkr = wkr - ss * wdi
        wki = wki + ss * wdr
        wdr = wdr - ss * wki
        wdi = wdi + ss * wkr
        do j = k, n - m + k, m
           i = j + l
           xr = a(j) - a(i)
           xi = a(j + 1) - a(i + 1)
           a(j) = a(j) + a(i)
           a(j + 1) = a(j + 1) + a(i + 1)
           a(i) = wkr * xr - wki * xi
           a(i + 1) = wkr * xi + wki * xr
           xr = a(j + 2) - a(i + 2)
           xi = a(j + 3) - a(i + 3)
           a(j + 2) = a(j + 2) + a(i + 2)
           a(j + 3) = a(j + 3) + a(i + 3)
           a(i + 2) = wdr * xr - wdi * xi
           a(i + 3) = wdr * xi + wdi * xr
        enddo
     enddo
     m = l
  enddo
  if (m .gt. 2) then
     do j = 0, n - 4, 4
        xr = a(j) - a(j + 2)
        xi = a(j + 1) - a(j + 3)
        a(j) = a(j) + a(j + 2)
        a(j + 1) = a(j + 1) + a(j + 3)
        a(j + 2) = xr
        a(j + 3) = xi
     enddo
  endif
  if (n .gt. 4) call bitrv2(n, a)

  
  do i = 0, n/2-1
     c(i) = dcmplx(a(2*i), a(2*i+1))
  enddo
  
  
end subroutine cdft



subroutine bitrv2(n, a)
  integer :: n, j, j1, k, k1, l, m, m2, n2
  real(kind(0d0)) :: a(0 : n - 1), xr, xi
  
  m = n / 4
  m2 = 2 * m
  n2 = n - 2
  k = 0
  do j = 0, m2 - 4, 4
     if (j .lt. k) then
        xr = a(j)
        xi = a(j + 1)
        a(j) = a(k)
        a(j + 1) = a(k + 1)
        a(k) = xr
        a(k + 1) = xi
     else if (j .gt. k) then
        j1 = n2 - j
        k1 = n2 - k
        xr = a(j1)
        xi = a(j1 + 1)
        a(j1) = a(k1)
        a(j1 + 1) = a(k1 + 1)
        a(k1) = xr
        a(k1 + 1) = xi
     endif
     k1 = m2 + k
     xr = a(j + 2)
     xi = a(j + 3)
     a(j + 2) = a(k1)
     a(j + 3) = a(k1 + 1)
     a(k1) = xr
     a(k1 + 1) = xi
     l = m
     do while (k .ge. l)
        k = k - l
        l = l / 2
     enddo
     k = k + l
  enddo
end subroutine bitrv2




