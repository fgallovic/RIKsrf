! Evaluates synthetic seismograms combining
! Ruiz Integral K-squared (RIK) source model
! (i.e. modified model by Ruiz et al., 2011)
! with DWN Green's functions from NEZsor.dat.
! Coded by Frantisek Gallovic (2015)
!-------------------------------------------------------------------------------

    MODULE crustaldat
    INTEGER ndepth
    REAL,ALLOCATABLE,DIMENSION(:):: depth,vp,vs,rho
    REAL hypodepth,dip
    END MODULE

    PROGRAM RIKSEISDWN
    USE crustaldat
    IMPLICIT NONE
    REAL,PARAMETER:: PI=3.1415926535d0
    REAL,PARAMETER:: rotateto=155.*(PI/180.)
    INTEGER nfmax,np,npRIK,NL,NW,NR
    INTEGER NLgf,NWgf
    REAL*8 T,artifDT,leng,widt,epicL,epicW,fc1,fc2
    REAL*8 elem,dL,dW,dt,df
    COMPLEX*16,DIMENSION(:,:),ALLOCATABLE:: gfN,gfE,gfZ
    COMPLEX*16,DIMENSION(:),ALLOCATABLE:: cseis
    COMPLEX*16,DIMENSION(:,:),ALLOCATABLE:: cirN,cirE,cirZ,sr,seisN,seisE,seisZ
    real, allocatable :: x1a(:),x2a(:)
    real, allocatable:: y2aN(:,:),ry2aN(:,:),dyry2aN(:),dyy2aN(:),rirN(:),iirN(:)
    real, allocatable:: y2aE(:,:),ry2aE(:,:),dyry2aE(:),dyy2aE(:),rirE(:),iirE(:)
    real, allocatable:: y2aZ(:,:),ry2aZ(:,:),dyry2aZ(:),dyy2aZ(:),rirZ(:),iirZ(:)
    real grre,grim,ll,ww
    REAL*4,DIMENSION(:),ALLOCATABLE:: fltr4
    REAL*8,DIMENSION(:,:),ALLOCATABLE:: mu
    INTEGER,ALLOCATABLE,DIMENSION(:):: fcsta
    real*4 dumarr(6),dtr4,f1r4,f2r4
    real*8 dum
    INTEGER i,j,jj,k,jw,jl,mm,dumi

    write(*,*)'Reading parameters...'

    open(10,file='input.dat',action='read')
    read(10,*)
    read(10,*) nfmax
    read(10,*)
    read(10,*) T
    read(10,*)
    read(10,*) artifDT
    read(10,*)
    read(10,*) NR
    read(10,*)
    read(10,*) NLgf,NWgf
    read(10,*)
    read(10,*) 
    read(10,*)
    read(10,*) dip,dip
    read(10,*)
    read(10,*) hypodepth
    read(10,*)
    read(10,*) leng,widt
    read(10,*)
    read(10,*) epicL,epicW
    read(10,*)
    read(10,*) np,npRIK
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*) dumi   !number of frequency bands
    if(dumi>1)stop 'Use a single frequency band!'
    read(10,*) fc1,fc2
    close(10)
    
    open(10,file='RIKsrf.in',action='read')
    do i=1,7
      read(10,*)
    enddo
    read(10,*) jj
    if(jj.ne.1)stop 'Only regular grid spacing is implemented now!'
    read(10,*) NL,NW
    if(NL==NLgf.and.NW==NWgf)write(*,*)'Note: Interpolation is NOT performed as it is not required.'
    dt=T/dble(np)
    df=1.d0/T
    dL=leng/dble(NL)
    dW=widt/dble(NW)
    elem=dL*dW

    allocate(mu(NL,NW))
    CALL read_crustaldat()
    do i=1,NW
      dum=(hypodepth+(epicW-widt/dble(NW)*(dble(i)-.5))*sin(dip/180.d0*PI))/1.d3
      if(dum>depth(ndepth))then
        mu(:,i)=rho(ndepth)*vs(ndepth)**2*1.d9
      else
        do j=1,ndepth
          if(dum<depth(j))exit
        enddo
        mu(:,i)=rho(j-1)*vs(j-1)**2*1.d9
      endif
    enddo

    write(*,*)'Reading slip rates...'
    ALLOCATE(sr(np,NL*NW))
    sr=0.d0
    open(201,FILE='sr.dat')
    i=0
    do jw=1,NW
      do jl=1,NL
        i=i+1
        do j=1,npRIK  !41
          read(201,*)sr(j,i)
        enddo
        sr(:,i)=sr(:,i)*mu(jl,jw)*elem
        read(201,*)
        read(201,*)
      enddo
    enddo
    close(201)

!$OMP parallel do private(i) DEFAULT(SHARED) SCHEDULE(DYNAMIC,1)
    do i=1,NL*NW
      call four1(sr(:,i),np,-1)
    enddo
!$omp end parallel do
    sr=sr*dt

    allocate(cirN(min(nfmax,np),NLgf*NWgf),cirE(min(nfmax,np),NLgf*NWgf),cirZ(min(nfmax,np),NLgf*NWgf))
    allocate(gfN(np,NLgf*NWgf),gfE(np,NLgf*NWgf),gfZ(np,NLgf*NWgf))
    allocate(seisN(np,NR),seisE(np,NR),seisZ(np,NR))
    seisN=0.d0;seisE=0.d0;seisZ=0.d0
    
    open(20,form='unformatted',file='NEZsor.dat')
    write(*,*)'Calculating seismograms...'    
!    write(*,*)'POZOR, PROHAZUJI GF HORIZONTALNE!'
    do jj=1,NR
      write(*,*)'Station ',jj
      j=0
      read(20) dumi
      do jw=1,NWgf
        do jl=1,NLgf
          j=j+1
!        do jl=NLgf,1,-1     !Tyto dva radky prohodi GF horizontalne!
!          j=(jw-1)*NLgf+jl
          read(20) dumi
          do k=1,nfmax
            read(20) dumarr
            if(k>np)cycle
            cirN(k,j)=cmplx(dumarr(1),dumarr(4))
            cirE(k,j)=cmplx(dumarr(2),dumarr(5))
            cirZ(k,j)=cmplx(dumarr(3),dumarr(6))
          enddo
        enddo
      enddo
      dtr4=dt
      f1r4=fc1
      f2r4=fc2
      do k=1,3
!$OMP parallel private(i,mm,cseis,fltr4) DEFAULT(SHARED)
        allocate(cseis(np),fltr4(np))
!$OMP do SCHEDULE(DYNAMIC,1)
        do i=1,NLgf*NWgf
          cseis=0.d0
          SELECT CASE (k)
          CASE(1)
            cseis(1:min(nfmax,np))=cirN(1:min(nfmax,np),i)
          CASE(2)
            cseis(1:min(nfmax,np))=cirE(1:min(nfmax,np),i)
          CASE(3)
            cseis(1:min(nfmax,np))=cirZ(1:min(nfmax,np),i)
          END SELECT
          cseis(np/2+2:np)=conjg(cseis(np/2:2:-1))
          cseis(np/2+1)=real(cseis(np/2+1))
          call four1(cseis,np,1)
          fltr4=real(cseis)*df
          do mm=1,int(artifDT/dt)  ! REMOVING DWN ARTIFACT FROM THE SEISMOGRAM BEGINNING
            fltr4(mm)=fltr4(mm)*(cos((dt*dble(mm-1)-artifDT)/artifDT*PI)/2.d0+.5d0);
          enddo
          if(f1r4>0.)then
            CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'BP', f1r4, f2r4, dtr4, 1, np)
          else
            CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'LP', f1r4, f2r4, dtr4, 1, np)
          endif

!          do mm=2,np   !time integration to get displacements
!            fltr4(mm)=fltr4(mm)+fltr4(mm-1)
!          enddo
!          fltr4=fltr4*dt
          
          cseis=fltr4*dt
          call four1(cseis,np,-1)
          SELECT CASE (k)
          CASE(1)
            gfN(1:np,i)=cseis(1:np)
          CASE(2)
            gfE(1:np,i)=cseis(1:np)
          CASE(3)
            gfZ(1:np,i)=cseis(1:np)
          END SELECT
        enddo
!$omp end do
        deallocate(cseis,fltr4)
!$omp end parallel
      enddo

! Convolution of SR with GF and integration along the fault
      if(NL==NLgf.and.NW==NWgf)then    !no interpolation needed
!$OMP parallel do private(i) DEFAULT(SHARED) SCHEDULE(DYNAMIC,1)
        do i=1,np
          seisN(i,jj)=sum(sr(i,:)*gfN(i,:))
          seisE(i,jj)=sum(sr(i,:)*gfE(i,:))
          seisZ(i,jj)=sum(sr(i,:)*gfZ(i,:))
        enddo
!$omp end parallel do
      else
!        write(*,*)'Interpolating GFs...'
        allocate(x1a(NLgf),x2a(NWgf))
        do jl=1,NLgf
          x1a(jl)=leng/real(NLgf)*(real(jl)-0.5)
        enddo
        do jw=1,NWgf
          x2a(jw)=widt/real(NWgf)*(real(jw)-0.5)
        enddo
!$OMP parallel  private(i,jw,jl,ww,ll,grre,grim,dyry2aN,dyy2aN,ry2aN,y2aN,rirN,iirN,dyry2aE,dyy2aE,ry2aE,y2aE,rirE,iirE,dyry2aZ,dyy2aZ,ry2aZ,y2aZ,rirZ,iirZ) DEFAULT(SHARED)
        allocate(dyry2aN(max(NLgf,NWgf)),dyy2aN(max(NLgf,NWgf)),y2aN(NWgf,NL),ry2aN(NWgf,NL),rirN(NLgf),iirN(NLgf))
        allocate(dyry2aE(max(NLgf,NWgf)),dyy2aE(max(NLgf,NWgf)),y2aE(NWgf,NL),ry2aE(NWgf,NL),rirE(NLgf),iirE(NLgf))
        allocate(dyry2aZ(max(NLgf,NWgf)),dyy2aZ(max(NLgf,NWgf)),y2aZ(NWgf,NL),ry2aZ(NWgf,NL),rirZ(NLgf),iirZ(NLgf))
!$OMP do SCHEDULE(DYNAMIC)
        do i=1,nfmax
          do jw=1,NWgf
            rirN(1:NLgf)=real(gfN(i,(jw-1)*NLgf+1:jw*NLgf))
            iirN(1:NLgf)=imag(gfN(i,(jw-1)*NLgf+1:jw*NLgf))
            rirE(1:NLgf)=real(gfE(i,(jw-1)*NLgf+1:jw*NLgf))
            iirE(1:NLgf)=imag(gfE(i,(jw-1)*NLgf+1:jw*NLgf))
            rirZ(1:NLgf)=real(gfZ(i,(jw-1)*NLgf+1:jw*NLgf))
            iirZ(1:NLgf)=imag(gfZ(i,(jw-1)*NLgf+1:jw*NLgf))
            call spline(x1a(1:NLgf),rirN(1:NLgf),NLgf,1.e30,1.e30,dyry2aN(1:NLgf))
            call spline(x1a(1:NLgf),iirN(1:NLgf),NLgf,1.e30,1.e30,dyy2aN(1:NLgf))
            call spline(x1a(1:NLgf),rirE(1:NLgf),NLgf,1.e30,1.e30,dyry2aE(1:NLgf))
            call spline(x1a(1:NLgf),iirE(1:NLgf),NLgf,1.e30,1.e30,dyy2aE(1:NLgf))
            call spline(x1a(1:NLgf),rirZ(1:NLgf),NLgf,1.e30,1.e30,dyry2aZ(1:NLgf))
            call spline(x1a(1:NLgf),iirZ(1:NLgf),NLgf,1.e30,1.e30,dyy2aZ(1:NLgf))
	        do jl=1,NL
	          ll=float(jl-1)*dL+dL/2.
	          call splint(x1a(1:NLgf),rirN(1:NLgf),dyry2aN(1:NLgf),NLgf,ll,ry2aN(jw,jl))
              call splint(x1a(1:NLgf),iirN(1:NLgf),dyy2aN(1:NLgf),NLgf,ll,y2aN(jw,jl))
	          call splint(x1a(1:NLgf),rirE(1:NLgf),dyry2aE(1:NLgf),NLgf,ll,ry2aE(jw,jl))
              call splint(x1a(1:NLgf),iirE(1:NLgf),dyy2aE(1:NLgf),NLgf,ll,y2aE(jw,jl))
	          call splint(x1a(1:NLgf),rirZ(1:NLgf),dyry2aZ(1:NLgf),NLgf,ll,ry2aZ(jw,jl))
              call splint(x1a(1:NLgf),iirZ(1:NLgf),dyy2aZ(1:NLgf),NLgf,ll,y2aZ(jw,jl))
            enddo
	      enddo
	      do jl=1,NL
            call spline(x2a(1:NWgf),ry2aN(1:NWgf,jl),NWgf,1.e30,1.e30,dyry2aN(1:NWgf))
            call spline(x2a(1:NWgf),y2aN(1:NWgf,jl),NWgf,1.e30,1.e30,dyy2aN(1:NWgf))
            call spline(x2a(1:NWgf),ry2aE(1:NWgf,jl),NWgf,1.e30,1.e30,dyry2aE(1:NWgf))
            call spline(x2a(1:NWgf),y2aE(1:NWgf,jl),NWgf,1.e30,1.e30,dyy2aE(1:NWgf))
            call spline(x2a(1:NWgf),ry2aZ(1:NWgf,jl),NWgf,1.e30,1.e30,dyry2aZ(1:NWgf))
            call spline(x2a(1:NWgf),y2aZ(1:NWgf,jl),NWgf,1.e30,1.e30,dyy2aZ(1:NWgf))
            do jw=1,NW
              ww=float(jw-1)*dW+dW/2.
              call splint(x2a(1:NWgf),ry2aN(1:NWgf,jl),dyry2aN(1:NWgf),NWgf,ww,grre)
              call splint(x2a(1:NWgf),y2aN(1:NWgf,jl),dyy2aN(1:NWgf),NWgf,ww,grim)
              seisN(i,jj)=seisN(i,jj)+cmplx(grre,grim)*sr(i,(jw-1)*NL+jl)
              call splint(x2a(1:NWgf),ry2aE(1:NWgf,jl),dyry2aE(1:NWgf),NWgf,ww,grre)
              call splint(x2a(1:NWgf),y2aE(1:NWgf,jl),dyy2aE(1:NWgf),NWgf,ww,grim)
              seisE(i,jj)=seisE(i,jj)+cmplx(grre,grim)*sr(i,(jw-1)*NL+jl)
              call splint(x2a(1:NWgf),ry2aZ(1:NWgf,jl),dyry2aZ(1:NWgf),NWgf,ww,grre)
              call splint(x2a(1:NWgf),y2aZ(1:NWgf,jl),dyy2aZ(1:NWgf),NWgf,ww,grim)
              seisZ(i,jj)=seisZ(i,jj)+cmplx(grre,grim)*sr(i,(jw-1)*NL+jl)
            enddo
          enddo
        enddo
!$omp end do
        deallocate(dyry2aN,dyy2aN,y2aN,ry2aN,rirN,iirN)
        deallocate(dyry2aE,dyy2aE,y2aE,ry2aE,rirE,iirE)
        deallocate(dyry2aZ,dyy2aZ,y2aZ,ry2aZ,rirZ,iirZ)
!$omp end parallel
        seisN(np/2+2:np,jj)=conjg(seisN(np/2:2:-1,jj))
        seisE(np/2+2:np,jj)=conjg(seisE(np/2:2:-1,jj))
        seisZ(np/2+2:np,jj)=conjg(seisZ(np/2:2:-1,jj))
        deallocate(x1a,x2a)
      endif
    
      call four1(seisN(:,jj),np,1)
      call four1(seisE(:,jj),np,1)
      call four1(seisZ(:,jj),np,1)
      
    enddo
    seisN=seisN*df
    seisE=seisE*df
    seisZ=seisZ*df

    close(20)
    deallocate(cirN,cirE,cirZ)

    write(*,*)'Saving seismograms...'        
    open(201,FILE='svseisNrik.dat')
    open(202,FILE='svseisErik.dat')
    open(203,FILE='svseisZrik.dat')
    do i=1,np
      write(201,'(1000E13.5)')dt*(i-1),dble(seisN(i,1:NR))*cos(rotateto)+dble(seisE(i,1:NR))*sin(rotateto)
      write(202,'(1000E13.5)')dt*(i-1),dble(seisE(i,1:NR))*cos(rotateto)-dble(seisN(i,1:NR))*sin(rotateto)
      write(203,'(1000E13.5)')dt*(i-1),dble(seisZ(i,1:NR))
    enddo
    
    END

    
    ! Using crustal.dat
    SUBROUTINE read_crustaldat()
    USE crustaldat
    IMPLICIT NONE
    INTEGER i
    if(allocated(depth))return
    open(10,FILE='crustal.dat',ACTION='READ',STATUS='OLD')
    write(*,*)'  (Using mu values from file crustal.dat)'
    read(10,*)
    read(10,*)
    read(10,*)ndepth
    allocate(depth(ndepth),vp(ndepth),vs(ndepth),rho(ndepth))
    read(10,*)
    read(10,*)
    do i=1,ndepth
      read(10,*)depth(i),vp(i),vs(i),rho(i)
    enddo
    close(10)
    END


!Numerical recipes
    
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=5000)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
    
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
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
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END

      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      DOUBLE PRECISION data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=dble(wr)*data(j)-dble(wi)*data(j+1)
            tempi=dble(wr)*data(j+1)+dble(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END

