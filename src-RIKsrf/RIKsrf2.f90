! Generator of slip rate functions for the
! Ruiz Integral K-squared (RIK) source model
! (i.e. modified model by Ruiz et al., 2011, see Gallovic, 2016)
! with random rupture front
! Coded by Frantisek Gallovic (2016)
! Version 2.0
!-------------------------------------------------------------------------------

    MODULE ruptveloc
    INTEGER nlayers
    REAL, ALLOCATABLE:: hvr(:),vr(:)   ! layer-top location (from bottom), rupture velocities inside layers
    REAL hl,hw                         ! location of hypocenter
    END MODULE
    
    MODULE crustaldat
    INTEGER ndepth
    REAL,ALLOCATABLE,DIMENSION(:):: depth,vp,vs,rho
    REAL hypodepth,dip
    END MODULE
  
    PROGRAM KKstf
    USE ruptveloc
    USE crustaldat
    IMPLICIT NONE
    REAL,PARAMETER:: PI=3.1415926535
    INTEGER NL,NW
    REAL mufix
!PDF for subsource position
    INTEGER,PARAMETER:: pdfNL=200,pdfNW=100
    REAL,ALLOCATABLE:: pdf2D(:,:),cpdf2D(:,:)
    REAL pdfDL,pdfDW,pdfGaussL,pdfGaussW,pdfGaussS
    REAL,ALLOCATABLE:: ruptimegen(:)
    CHARACTER*256 filename
    INTEGER ml(2),pdfOption,fileNL,fileNW
!SUBSOURCE PARAMETERS:
    INTEGER SRoption,SUBmax,SUBmin
    REAL,ALLOCATABLE,DIMENSION(:):: SUBposL,SUBposW,SUBsize,SUBslip,SUBmoment,SUBrisetime,SUBmvr,SUBruptime,SUBnuclL,SUBnuclW
    INTEGER,ALLOCATABLE,DIMENSION(:):: SUBno
    INTEGER SUBtot,idum1,idum2
!SOURCE PARAMETERS:
    REAL LF,WF,L,W,smL,smW,LWratio,M0,L0,vrsubfact,aparam,dt
    REAL,ALLOCATABLE,DIMENSION(:):: SRl,SRw,SRelem,SRmu,SRslip,SRmoment,SRstressdrop,SR,STF
    INTEGER NSR,NT
!others
    REAL ran2,dum,dumL,dumW,dumphi,dumr,totmoment,time,meanVR,ruptime,ruptimeSR,momentcorr
    INTEGER i,j,k,m,hits
    
    open(101,FILE='RIKsrf.in')
    read(101,*)
    read(101,*)LF, WF
    read(101,*)
    read(101,*)L, W, smL,smW
    if(smL+L>LF.or.smW+W>WF)then
        write(*,*)'Strong-motion area exceeds the fault size!'
        stop
    endif
    read(101,*)
    read(101,*)M0
    read(101,*)
    read(101,*)SRoption
    if(SRoption==1)then
      read(101,*)NL,NW,mufix
      NSR=NL*NW
      if(mufix==0.)CALL read_crustaldat()
    elseif(SRoption==2)then
      read(101,*)NL,NW,NSR
    else
      write(*,*)'Wrong SRoption!'
      stop
    endif
    read(101,*)
    read(101,*)hl,hw
    read(101,*)
    read(101,*)hypodepth,dip
    read(101,*)
    read(101,*)L0,vrsubfact,aparam
    read(101,*)
    read(101,*)SUBMIN,SUBMAX
    read(101,*)
    read(101,*)pdfOption
    if(pdfOption==2)read(101,*)pdfGaussL,pdfGaussW,pdfGaussS
    if(pdfOption==3)read(101,*)fileNL,fileNW,filename
    read(101,*)
    read(101,*)idum1,idum2
    read(101,*)
    read(101,*)dt,NT
    read(101,*)
    read(101,*)nlayers
    if(nlayers>0)then
      allocate(hvr(nlayers),vr(nlayers))
      do i=1,nlayers
        read(101,*)hvr(i),vr(i)
      enddo
    else
      call read_crustaldat()
      read(101,*)dum   !constant vr/vs
      nlayers=ndepth
      allocate(hvr(nlayers),vr(nlayers))
      do i=1,nlayers
        hvr(i)=(hypodepth-depth(nlayers-i+1))/sin(dip/180.d0*PI)+hw
        vr(i)=vs(nlayers-i+1)*dum
      enddo
    endif
    do i=2,nlayers
      if(vr(i)>vr(i-1))then
        write(*,*)'ERROR! VR should decrease upwards!'
      endif
    enddo
    if(hvr(nlayers)<WF)then
      write(*,*)'ERROR! Definition of VR does not cover the whole fault!'
      stop
    endif
    close(101)

    open(232,FILE='nucleationpoint.dat')
    write(232,*)hl,hw
    close(232)
    LWratio=L/W
    L0=L0*W

!Preparing PDF for distribution of subsources
    write(*,*)'Preparing PDF for subsource distribution...'
    ALLOCATE(pdf2D(pdfNL,pdfNW),cpdf2D(pdfNL,pdfNW))
    CALL fillpdf(pdfNL,pdfNW,LF,WF,L,W,smL,smW,pdf2D,cpdf2D,pdfOption,filename,fileNL,fileNW,pdfGaussL,pdfGaussW,pdfGaussS)
    pdfDL=LF/real(pdfNL)
    pdfDW=WF/real(pdfNW)
    write(*,*)'... done.'

!Calculate and read random variations of rupture time
    CALL randomruptvel(LF,WF,NL,NW,hl,hw,idum2)
    ALLOCATE(ruptimegen(NW*NL))
    open(232,FILE='ruptimegen.txt')
    k=0
    do j=1,NW
      do i=1,NL
        k=k+1
        read(232,*)dum,dum,ruptimegen(k)
      enddo
      read(232,*)
    enddo

!Reading location of the slip rate points
    ALLOCATE(SRl(NSR),SRw(NSR),SRelem(NSR),SRmu(NSR),SRslip(NSR),SRmoment(NSR),SRstressdrop(NSR))
    if(SRoption==1)then
      SRelem=LF*WF/real(NL*NW)
      k=0
      do j=1,NW
        do i=1,NL
          k=k+1
          SRl(k)=LF/real(NL)*(real(i)-.5)
          SRw(k)=WF/real(NW)*(real(j)-.5)
        enddo
      enddo
      if(mufix>0.)then
        SRmu=mufix
      else
        do k=1,NSR
          dum=(hypodepth+(hw-SRw(k))*sin(dip/180.d0*PI))
          if(dum>depth(ndepth))then
            SRmu(k)=rho(ndepth)*vs(ndepth)**2*1.d9
          else
            do j=1,ndepth
              if(dum<depth(j))exit
            enddo
            SRmu(k)=rho(j-1)*vs(j-1)**2*1.d9
          endif  
        enddo
      endif
    else
      open(101,FILE='srloc.dat')
      do i=1,NSR
        read(101,*)SRl(i),SRw(i),SRelem(i),SRmu(i)
      enddo
    endif
    SRelem=SRelem*1.e6
    
!Static subsource parameters
    write(*,*)'Preparing distribution and parameters of subsources...'
    ALLOCATE(SUBno(SUBmax))
    SUBno=0
    do i=SUBmin,SUBmax
      SUBno(i)=int(float(2*i-1)*LWratio)
    enddo
    SUBtot=sum(SUBno)
    ALLOCATE(SUBposL(SUBtot),SUBposW(SUBtot),SUBsize(SUBtot),SUBslip(SUBtot),SUBmoment(SUBtot),SUBrisetime(SUBtot))
    ALLOCATE(SUBnuclL(SUBtot),SUBnuclW(SUBtot),SUBmvr(SUBtot),SUBruptime(SUBtot))
    k=0
    SRslip=0.
    SRmoment=0.
    SRstressdrop=0.
    SUBslip=0.
    SUBmoment=0.
    do i=SUBmin,SUBmax
      do j=1,SUBno(i)
        k=k+1
        SUBsize(k)=W/float(i)/2.   !radius
       
! Locating the subsource according to the PDF
        do
!          if(k==1.and.pdfOption.ne.1)then   !Put the largest subsource at the position of the PDF maximum
!            ml=maxloc(pdf2D(:,:))
!          else
            ml=minloc(abs(cpdf2D(:,:)-ran2(idum1)))
!          endif
          SUBposL(k)=(real(ml(1))-.5)*pdfDL
          
          if(SUBmin==1.and.i==1)then
            SUBposW(k)=W/2.
            if(SUBposL(k)-SUBsize(k)>=0..and.SUBposL(k)+SUBsize(k)<=LF)exit
          else
            SUBposW(k)=(real(ml(2))-.5)*pdfDW
            if(SUBposL(k)-SUBsize(k)>=0..and.SUBposW(k)-SUBsize(k)>=0..and.SUBposL(k)+SUBsize(k)<=LF.and.SUBposW(k)+SUBsize(k)<=WF)exit
          endif
        enddo

        do m=1,NSR
          dum=SUBsize(k)**2-(SRl(m)-SUBposL(k))**2-(SRw(m)-SUBposW(k))**2
          if(dum>0.)then
            SUBmoment(k)=SUBmoment(k)+sqrt(dum)*SRelem(m)*SRmu(m)
            SUBslip(k)=SUBslip(k)+sqrt(dum)
            SRmoment(m)=SRmoment(m)+sqrt(dum)*SRelem(m)*SRmu(m)
            SRslip(m)=SRslip(m)+sqrt(dum)
            SRstressdrop(m)=SRstressdrop(m)-SRmu(m)/24.*7.*PI/1000. !(oprava na to, ze se pracuje se skluzem v km)
          endif
        enddo
      enddo
    enddo
    totmoment=sum(SRmoment(:))
    momentcorr=M0/totmoment
    write(*,*)totmoment*momentcorr
    SUBslip=SUBslip*momentcorr
    SUBmoment=SUBmoment*momentcorr
    SRslip=SRslip*momentcorr
    SRmoment=SRmoment*momentcorr
    SRstressdrop=SRstressdrop*momentcorr

    open(201,FILE='slipdistribution.dat')
    do i=1,NSR
      write(201,'(10E13.5)')SRl(i),SRw(i),SRslip(i),SRmoment(i),ruptimegen(i),SRstressdrop(i)
    enddo
    close(201)
    open(201,FILE='slipdistribution.gnuplot.dat')
    do i=1,NW
      write(201,'(1000E13.5)')SRslip((i-1)*NL+1:i*NL)
    enddo
    write(201,*);write(201,*)
    do i=1,NW
      write(201,'(1000E13.5)')SRstressdrop((i-1)*NL+1:i*NL)
    enddo
    close(201)

!Slip rates on subsources
    do k=1,SUBtot
      SUBmvr(k)=meanVR(SUBposL(k),SUBposW(k),SUBsize(k))
      if(2.*SUBsize(k)>=L0)then    !subsources start from the hypocenter
        SUBrisetime(k)=aparam*L0/SUBmvr(k)
      else                         !subsources start from a random point
        dumphi=ran2(idum1)*2.*pi;dumr=sqrt(ran2(idum1))*SUBsize(k)
        dumL=dumr*cos(dumphi);dumW=dumr*sin(dumphi)
        SUBnuclL(k)=dumL+SUBposL(k)
        SUBnuclW(k)=dumW+SUBposW(k)
         SUBruptime(k)=ruptimegen(int(SUBnuclW(k)/W*float(NW-1))*NL+int(SUBnuclL(k)/L*float(NL-1))+1)
        SUBmvr(k)=SUBmvr(k)*vrsubfact
        SUBrisetime(k)=aparam*2.*SUBsize(k)/SUBmvr(k)
      endif
    enddo
    write(*,*)'... done.'

    open(201,FILE='subsources.dat')
    do k=1,SUBtot
      write(201,'(10E13.5)')SUBposL(k),SUBposW(k),SUBsize(k),SUBmoment(k),SUBruptime(k),SUBrisetime(SUBtot),SUBmvr(SUBtot)
    enddo
    close(201)
    
!Evaluating slip rates
    write(*,*)'Preparing and saving slip rates...'
    allocate(sr(NT),stf(NT))
    open(201,FILE='sr.dat')
    totmoment=0.
    stf=0.
    do i=1,NSR
      sr=0.
      ruptimeSR=ruptimegen(i)    !Comment to go back to the version without rupt. vel. perturbations
!$OMP parallel do private(j,k,time,ruptime,dum) DEFAULT(SHARED)
      do j=1,NT
        time=dt*(j-1)
        do k=1,SUBtot
          dum=SUBsize(k)**2-(SRl(i)-SUBposL(k))**2-(SRw(i)-SUBposW(k))**2
          if(dum>0.)then
            if(2.*SUBsize(k)>=L0)then    !subsources start from the hypocenter
              ruptime=ruptimeSR
            else
              ruptime=SUBruptime(k)+sqrt((SRl(i)-SUBnuclL(k))**2+(SRw(i)-SUBnuclW(k))**2)/SUBmvr(k)
!ruptime=ruptimeSR    !Warning, uncomment if you want all subsources to start from the hypocenter
            endif
            if(time>ruptime.and.time<ruptime+SUBrisetime(k)*5.)then
              sr(j)=sr(j)+sqrt(dum)*momentcorr*(time-ruptime)*exp(-(time-ruptime)/SUBrisetime(k)*PI)/(SUBrisetime(k)/PI)**2
            endif
          endif
        enddo
      enddo
!$OMP end parallel do
      stf(:)=stf(:)+sr(:)*SRelem(i)*SRmu(i)
      totmoment=totmoment+sum(sr(:))*SRelem(i)*SRmu(i)*dt
      do j=1,NT
!        write(201,*)dt*(j-1),sr(j)
        write(201,*)sr(j)
      enddo
      write(201,*)
      write(201,*)
    enddo
    close(201)
    write(*,*)totmoment
    open(201,FILE='stf.dat')
    do i=1,NT
      write(201,*)dt*(i-1),stf(i)
    enddo
    write(*,*)'... done.'
    
    END PROGRAM
    
    
    SUBROUTINE fillpdf(pdfNL,pdfNW,LF,WF,L,W,smL,smW,pdf2D,cpdf2D,pdfOption,filename,fileNL,fileNW,pdfGaussL,pdfGaussW,pdfGaussS)  ! creates pdf and cumulative pdf for subsource distribution
    IMPLICIT NONE
    INTEGER pdfNL,pdfNW,pdfOption
    REAL pdf2D(pdfNL,pdfNW),cpdf2D(pdfNL,pdfNW)
    REAL LF,WF,L,W,smL,smW,pdfGaussL,pdfGaussW,pdfGaussS
    CHARACTER*256 filename
    INTEGER i,j,k,fileNL,fileNW
    REAL cumul,pdfDL,pdfDW,slipDL,slipDW
    REAL,ALLOCATABLE:: slip(:,:)
    INTEGER pifrom,pito,pjfrom,pjto
    open(229,FILE='strongmotionarea.dat')
    write(229,*)smL,smW;write(229,*)smL+L,smW;write(229,*)smL+L,smW+W;write(229,*)smL,smW+W;write(229,*)smL,smW
    close(229)
    pdfDL=LF/real(pdfNL)
    pdfDW=WF/real(pdfNW)
    pifrom=int(smL/pdfDL)+1
    pito=int((smL+L)/pdfDL+0.999)
    pjfrom=int(smW/pdfDW)+1
    pjto=int((smW+W)/pdfDW+0.999)
    pdf2D(:,:)=0.
    SELECT CASE(pdfOption)
    CASE(1)
      write(*,*)'Uniform spatial PDF for subsources'
      pdf2D(pifrom:pito,pjfrom:pjto)=1.
    CASE(2)
      write(*,*)'Gaussian PDF for subsources'
      do j=pjfrom,pjto
        do i=pifrom,pito
          pdf2D(i,j)=exp(-.5*((((real(i)-.5)*pdfDL-pdfGaussL)**2+((real(j)-.5)*pdfDL-pdfGaussW)**2)/pdfGaussS**2)**2)
        enddo
      enddo
    CASE(3)
      write(*,*)'Reading spatial PDF from',trim(filename)
      allocate(slip(fileNL,fileNW))
      slipDL=LF/real(fileNL)
      slipDW=WF/real(fileNW)
      open(329,FILE=trim(filename))
      do j=1,fileNW
        read(329,*)(slip(i,j),i=1,fileNL)
      enddo
      close(329)
      do j=pjfrom,pjto
        do i=pifrom,pito
          pdf2D(i,j)=slip(int(pdfDL/slipDL*(float(i)-0.5))+1,int(pdfDW/slipDW*(float(j)-0.5))+1)
        enddo
      enddo
      deallocate(slip)
    CASE DEFAULT
      write(*,*)'Wrong pdfOption!'
      stop
    END SELECT
!normalize and calculate cumulative distribution
    k=0
    cumul=0
    do j=1,pdfNW
      do i=1,pdfNL
        k=k+1
        cumul=cumul+pdf2D(i,j)
        cpdf2D(i,j)=cumul
      enddo
    enddo
    pdf2D=pdf2D/cumul
    cpdf2D=cpdf2D/cumul
    END SUBROUTINE
    
    
    FUNCTION meanVR(x,y,r)   !calculate mean rupture velocity (just slowness mean over subsource depth extent)
    USE ruptveloc
    IMPLICIT NONE
    REAL meanVR,x,y,r
    INTEGER itop,ibottom,i
    itop=1
    ibottom=1
    do i=1,nlayers
      if(hvr(i)<y+r) itop=i+1
      if(hvr(i)<y-r) ibottom=i+1
    enddo
    if(itop==ibottom)then    ! stf point in the layer with the hypocenter
      meanVR=VR(itop)
    else
      meanVR=(y+r-hvr(itop-1))/vr(itop)+(hvr(ibottom)-y+r)/vr(ibottom)
      do i=ibottom+1,itop-1
        meanVR=meanVR+(hvr(i)-hvr(i-1))/vr(i)
      enddo
      meanVR=2.*r/meanVR
    endif
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
    

    SUBROUTINE randomruptvel(L,W,NX,NY,epicx,epicy,idum)
! Random rupture velocity k^-2 generator for
! Ruiz Integral K-squared (RIK) source model
! Coded by Frantisek Gallovic (2014)
!-------------------------------------------------------------------------------
! COORDINATE SYSTEM on the fault:
! The origin is in the left top corner of the fault while the strike direction is to the right,
! the x axes is positive to the strike direction and
! the y axes is positive in the down-dip direction.
    USE ruptveloc
    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535,vars=0.25d0  !maximum variations in percents
    REAL*8,PARAMETER:: UpperK=1.d0
    COMPLEX*16,ALLOCATABLE:: speq1(:),AC(:,:),speqd(:)
    REAL*8,ALLOCATABLE:: A(:,:),AA(:,:),C(:,:),D(:,:),D1(:,:)
    INTEGER i,j,k,NXX,NYY,M,N,FNN,FNM,NX,NY
    REAL L,W,epicx,epicy
    REAL*8 dkx,dky,kx,ky
    REAL*8 dxout,dyout,dd
    REAL*8 cx,cy,ms,dum,NyqL,NyqW,NLW,krad,corner,KCx,KCy,ran2
    INTEGER ndepth1,idum,pdfOption
    integer*4 :: nnx, nnz, iostat, ipoint
    integer*4, external :: Time_2d
    real*4 xg, zg, eps_init
    real*4, allocatable:: hs(:), t_rupt(:)
    REAL*8, allocatable:: xintrpl(:),tintrpl(:),yintrpl(:),xspln(:),yspln(:)

    write(*,*)'Preparing rupture time variations with rupture velocity sigma ',vars

!Schvalne dd, protoze se diskretizuje na ridsi siti a pak se interpoluje
    M=2**int(log(dble(NX))/log(2.d0)+1.d0)
    N=2**int(log(dble(NY))/log(2.d0)+1.d0)
    dxout=L/dble(NX)
    dyout=W/dble(NY)
    dd=min(L/dble(NX),W/dble(NY))
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
      do j=1,N
        dum=(dd*(j-1)+dd/2.)
        if(dum>hvr(nlayers))then
          A(:,j)=vr(nlayers)*(1.+A(:,j))
        else
          do k=1,nlayers
            if(dum<hvr(k))then
              A(:,j)=vr(k)*(1.+A(:,j))
              exit
            endif
          enddo
        endif
      enddo

!Writing 2D rupture velocity distribution
    OPEN(102,FILE='ruptvelgen.txt')
    do j=1,int(W/dd)
      do i=1,int(L/dd)
        write(102,'(3E14.6)') real(i-1)*dd+.5*dd,real(j-1)*dd+.5*dd,A(i,j)
      enddo
      write(102,*)
    enddo
    close(102)

!Translating rupture velocities to rupture times
write(333,'(E13.5)')A
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

    iostat = Time_2d(hs, t_rupt, nnx, nnz, xg, zg, eps_init, 0)

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
    OPEN(101,FILE='ruptimegen.txt')
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
!    OPEN(106,FILE='specy.txt')
!    do i=1,N/2+1
!      write(106,*)(i-1)*dky,abs(AC(1,i))
!    enddo

!Writing amplitude spectrum along x:
!    OPEN(105,FILE='specx.txt')
!    do i=1,M/2
!      write(105,*)(i-1)*dkx,abs(AC(i,1))
!    enddo
!      write(105,*)(M/2)*dkx,abs(speq1(1))

    END    

    
! Numerical recipes    

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
      DOUBLE PRECISION AM,EPS,RNMX
      REAL ran2
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
    END FUNCTION
