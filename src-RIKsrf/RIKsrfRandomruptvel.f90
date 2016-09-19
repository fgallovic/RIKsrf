! Generator of slip rate functions for the
! Ruiz Integral K-squared (RIK) source model
! (i.e. modified model by Ruiz et al., 2011)
! with random rupture front
! Coded by Frantisek Gallovic (2014)
!-------------------------------------------------------------------------------

    MODULE ruptveloc
    INTEGER nlayers
    REAL, ALLOCATABLE:: hvr(:),vr(:)   ! layer-top location (from bottom), rupture velocities inside layers
    REAL hl,hw                         ! location of hypocenter
    INTEGER ihw                        ! layer number with hypocenter
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
    INTEGER SUBtot,idum
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
    read(101,*)idum
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
    if(hvr(1)<hw)then
      ihw=2
    else
      ihw=1
    endif
    do i=2,nlayers
      if(vr(i)>vr(i-1))then
        write(*,*)'ERROR! VR must decrease upwards!'
        stop
      endif
      if(hvr(i)<hw) ihw=i+1
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

!Read random variations of rupture time
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
            ml=minloc(abs(cpdf2D(:,:)-ran2(idum)))
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
!        do
!          dumL=ran2(idum)*SUBsize(k)-SUBsize(k)
!          dumW=ran2(idum)*SUBsize(k)-SUBsize(k)
!          if(dumL**2+dumW**2<=SUBsize(k)**2)exit
!        enddo
        dumphi=ran2(idum)*2.*pi;dumr=sqrt(ran2(idum))*SUBsize(k)
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
    
    
! Numerical recipes    
       
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
