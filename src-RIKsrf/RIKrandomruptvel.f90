! Random rupture velocity k^-2 generator for
! Ruiz Integral K-squared (RIK) source model
! Coded by Frantisek Gallovic (2014)
!-------------------------------------------------------------------------------
! COORDINATE SYSTEM on the fault:
! The origin is in the left top corner of the fault while the strike direction is to the right,
! the x axes is positive to the strike direction and
! the y axes is positive in the down-dip direction.

    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535,vars=0.25d0  !maximum variations in percents
    REAL*8,PARAMETER:: UpperK=1.d0
    COMPLEX*16,ALLOCATABLE:: speq1(:),AC(:,:),speqd(:)
    REAL*8,ALLOCATABLE:: A(:,:),AA(:,:),C(:,:),D(:,:),D1(:,:)
    INTEGER i,j,k,NXX,NYY,M,N,FNN,FNM,NX,NY
    REAL*8 dkx,dky,L,W,kx,ky,epicx,epicy,topdepth,hypodepth,dip,vsbeta
    REAL*8 dxout,dyout,dd
    REAL*8 cx,cy,ms,dum,NyqL,NyqW,NLW,krad,corner,KCx,KCy,ran2
    INTEGER ndepth1,idum,pdfOption
    integer*4 :: nnx, nnz, iostat, ipoint
    integer*4, external :: Time_2d
    real*4 xg, zg, eps_init
    real*4, allocatable:: hs(:), t_rupt(:), depth1(:),vs1(:),rho1(:)
    REAL*8, allocatable:: xintrpl(:),tintrpl(:),yintrpl(:),xspln(:),yspln(:)
    INTEGER nlayers
    REAL, ALLOCATABLE:: hvr(:),vr(:)   ! layer-top location (from bottom), rupture velocities inside layers


    OPEN(101,FILE='ruptimegen.txt')
    OPEN(102,FILE='ruptvelgen.txt')
!    OPEN(105,FILE='specx.txt')
!    OPEN(106,FILE='specy.txt')
    OPEN(130,FILE='RIKstf.in')
!    OPEN(120,FILE='ruptimex.txt')
!    OPEN(121,FILE='ruptimey.txt')

!Startup settings

    write(*,*)'Variations in rupt. vel.: ',vars

    read(130,*)
    read(130,*)L,W
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)NX,NY
    read(130,*)
    read(130,*)epicx,epicy
    read(130,*)
  	read(130,*)hypodepth,dip
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)pdfOption
    if(pdfOption>1)read(130,*)
    read(130,*)
    read(130,*)idum,idum
    write(*,*)idum
    read(130,*)
    read(130,*)
    read(130,*)
    read(130,*)nlayers
    if(nlayers>0)then
      allocate(hvr(nlayers),vr(nlayers))
      do i=1,nlayers
        read(101,*)hvr(i),vr(i)
      enddo
    else
      read(130,*)vsbeta
    endif

!Schvalne dd, protoze se diskretizuje na ridsi siti a pak se interpoluje
    M=2**int(log(dble(NX))/log(2.d0)+1.d0)
    N=2**int(log(dble(NY))/log(2.d0)+1.d0)
    dxout=L/real(NX)
    dyout=W/real(NY)
    dd=min(L/real(NX),W/real(NY))
11  if(dd*M<L)then
      M=M*2
      goto 11
    endif
12  if(dd*N<W)then
      N=N*2
      goto 12
    endif
    FNN=N/2+1;FNM=M/2+1

    ALLOCATE(speq1(N),AC(M/2,N),speqd(N),A(M,N),AA(M,N),C(M,N),D(NX,NY),D1(M,NY))
    ALLOCATE(hs(M*N),t_rupt(M*N))
    ALLOCATE(xintrpl(M),xspln(M),yintrpl(N),tintrpl(N),yspln(N))

    dkx=1./dd/real(M);dky=1./dd/real(N)
    KCx=UpperK/L        !Corner wave-number for along-strike direction
    KCy=UpperK/W        !Corner wave-number for down-dip direction
    NyqL=dkx;NyqW=dky
    NLW=sqrt(NyqL**2+NyqW**2)

!Preparing white noise spectrum

      do i=1,M
        do j=1,N
          A(i,j)=ran2(idum)-.5d0
        enddo
      enddo
      CALL rlft3(A,speq1,M,N,1,1)
!      speq1=exp(cmplx(0.,atan2(imag(speq1),real(speq1))))
      do i=1,M/2
        AC(i,:)=cmplx(A(2*i-1,:),A(2*i,:))
      enddo

!Adding k^-2 by using the white noise spectrum

      do j=1,N
        if(j<=N/2+1)then
          ky=dky*real(j-1)
        else
          ky=-dky*real(N-j+1)
        endif
        do i=1,M/2+1
          kx=dkx*real(i-1)
          krad=sqrt(kx**2+ky**2)
          if(i<M/2+1.)then
            if(krad>=NLW)then
              AC(i,j)=AC(i,j)/sqrt(1.+((kx/KCx)**2+(ky/KCy)**2)**1)
            endif
          elseif(krad>=NLW)then
            speq1(j)=speq1(j)/sqrt(1.+((kx/KCx)**2+(ky/KCy)**2)**1)
          endif
        enddo
      enddo

!Back Fourier transform

      CALL rlft3(AC,speq1,M,N,1,-1)
      do i=1,M/2
        A(2*i-1,:)=real(AC(i,:))/M/N*2.
        A(2*i,:)=imag(AC(i,:))/M/N*2.
      enddo

!Adding the mean rupture velocity
      dum=sqrt(sum(A(:,:)**2)/M/N)
      A=A/dum*vars
      if(nlayers>0)then
        write(*,*)'Not yet coded' !- ty rychlosti se musi nacist z hvr a vr
        stop
      else
        write(*,*)'Obtaining rupture velocities from crustal.dat'
        open(10,FILE='crustal.dat')
        read(10,*)
        read(10,*)
        read(10,*)ndepth1
        allocate(depth1(ndepth1),vs1(ndepth1),rho1(ndepth1))
        read(10,*)
        read(10,*)
        do i=1,ndepth1
          read(10,*)depth1(i),dum,vs1(i),rho1(i)
        enddo
        close(10)
        do j=1,N
          dum=(dd*(j-1)+dd/2.)
          dum=(hypodepth+sin(dip/180.d0*PI)*(epicy-dum))
          if(dum>depth1(ndepth1))then
	          A(:,j)=vsbeta*vs1(ndepth1)*(1.+A(:,j))
          elseif(dum<=0.d0)then
            A(:,j)=vsbeta*vs1(1)*(1.+A(:,j))
	        else
	          do k=1,ndepth1
	            if(dum<depth1(k))exit
            enddo
            A(:,j)=vsbeta*vs1(k-1)*(1.+A(:,j))
          endif
        enddo
      endif

!Writing 2D rupture velocity distribution
    do j=1,int(W/dd)
      do i=1,int(L/dd)
        write(102,'(3E13.6)') real(i-1)*dd+.5*dd,real(j-1)*dd+.5*dd,A(i,j)
      enddo
      write(102,*)
    enddo
    
!Translating rupture velocities to rupture times

    do j = 1,N
      do i = 1,M
        ipoint = i + (j-1) * M
        hs(ipoint) = dd/A(i,j)
      enddo
    enddo
    eps_init = 0.
    nnx=M
    nnz=N
    xg=epicx/dd+1.
    zg=epicy/dd+1.

    iostat = Time_2d(hs, t_rupt, nnx, nnz, xg, zg, eps_init, 0);write(*,*)iostat

!PREINTERPOLOVANI DO VYSTUPNI DISKRETIZACE:
    do j=1,N
      yintrpl(j)=dd*(j-1)+dd/2.d0
    enddo
    do i=1,M
      do j=1,N
        ipoint = i + (j-1) * M
        tintrpl(j)=t_rupt(ipoint)
        A(i,j)=t_rupt(ipoint)
      enddo
      CALL spline(yintrpl,tintrpl,N,1.d30,1.d30,yspln)
      do j=1,NY
        dum=(j-1)*dyout+dyout/2.d0
        call splint(yintrpl,tintrpl,yspln,N,dum,D1(i,j))
      enddo
    enddo
    do i=1,M
      xintrpl(i)=dd*(i-1)+dd/2.d0  !Schvalne dy, protoze se diskretizuje na ridsi siti a pak se interpoluje
    enddo
    do j=1,NY
      CALL spline(xintrpl,D1(:,j),M,1.d30,1.d30,xspln)
      do i=1,NX
        dum=(i-1)*dxout+dxout/2.d0
        call splint(xintrpl,D1(:,j),xspln,M,dum,D(i,j))
      enddo
    enddo
    do j=1,NY
      do i=1,NX
        write(101,'(3E13.6)') real(i-1)*dxout+.5*dxout,real(j-1)*dyout+.5*dyout,D(i,j)
      enddo
      write(101,*)
    enddo

!Forward Fourier transform

    CALL rlft3(A,speq1,M,N,1,1)
    do i=1,M/2
      AC(i,:)=cmplx(A(2*i-1,:),A(2*i,:))
    enddo
    AC=AC/real(M*N/2);speq1=speq1/real(M*N/2)

!Writing amplitude spectrum along y:

!    do i=1,N/2+1
!      write(106,*)(i-1)*dky,abs(AC(1,i))
!    enddo

!Writing amplitude spectrum along x:

!    do i=1,M/2
!      write(105,*)(i-1)*dkx,abs(AC(i,1))
!    enddo
!      write(105,*)(M/2)*dkx,abs(speq1(1))

    END


!Subroutines from Numerical Receipes

      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX*16 data(nn1/2,nn2,nn3),speq(nn2,nn3)
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX*16 c1,c2,h1,h2,w
      c1=dcmplx(0.5d0,0.0d0)
      c2=dcmplx(0.0d0,-0.5d0*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=dcmplx(dble(wr),dble(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END




      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      DOUBLE PRECISION data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=dble(wr)*data(k2)-dble(wi)*data(k2+1)
                tempi=dble(wr)*data(k2+1)+dble(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=50000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      END

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=3.d-16,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END 