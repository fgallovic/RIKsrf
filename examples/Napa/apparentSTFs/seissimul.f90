! Evaluates apparent source functions for the
! Ruiz Integral K-squared (RIK) source model
! (i.e. modified model by Ruiz et al., 2011)
! Coded by Frantisek Gallovic (2014)
! (strike is zero in the present implementation)
!-------------------------------------------------------------------------------

    MODULE crustaldat
    INTEGER ndepth
    REAL,ALLOCATABLE,DIMENSION(:):: depth,vp,vs,rho
    REAL hypodepth,dip
    END MODULE

    PROGRAM seissimul
    USE crustaldat
    IMPLICIT NONE
!TO GENERATE SLIP RATE POINTS - will not be needed later
    REAL,PARAMETER:: VSt=3.5
    REAL mufix
    INTEGER NL,NW
!PI
    REAL,PARAMETER:: PI=3.1415926535d0
!SOURCE PARAMETERS:
    INTEGER SRoption
    REAL L,W,M0,hl,hw,dt
    REAL,ALLOCATABLE,DIMENSION(:):: SRl,SRw,SRelem,SRmu
    REAL,ALLOCATABLE,DIMENSION(:,:):: sr
    INTEGER NSR,NT,idum
!STATIONS
    INTEGER,PARAMETER:: NR=10,NT2=4096  !Napa
    REAL staL(NR),staY(NR)
    REAL,ALLOCATABLE,DIMENSION(:,:):: seis
!others
    REAL dum
    INTEGER i,j,k,tshift
    
    open(101,FILE='RIKsrf.in')
    read(101,*)
    read(101,*)L, W
    read(101,*)
    read(101,*)
    read(101,*)
    read(101,*)M0
    read(101,*)
    read(101,*)SRoption
    if(SRoption==1)then
      read(101,*)NL,NW,mufix
      NSR=NL*NW
      if(mufix==0.)CALL read_crustaldat()
    elseif(SRoption==2)then
      read(101,*)NSR
    else
      write(*,*)'Wrong SRoption!'
      stop
    endif
    read(101,*)
    read(101,*)hl,hw
    read(101,*)
    read(101,*)hypodepth,dip
    read(101,*)
    read(101,*)
    read(101,*)
    read(101,*)
    read(101,*)
    read(101,*)i
    if(i>1)read(101,*)
    read(101,*)
    read(101,*)
    read(101,*)
    read(101,*)dt,NT
    close(101)
            
!Reading location of the slip rate points
    ALLOCATE(SRl(NSR),SRw(NSR),SRelem(NSR),SRmu(NSR))
    if(SRoption==1)then
      SRelem=L*W/real(NL*NW)
      k=0
      do j=1,NW
        do i=1,NL
          k=k+1
          SRl(k)=L/real(NL)*(real(i)-.5)
          SRw(k)=W/real(NW)*(real(j)-.5)
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

    ALLOCATE(sr(NT,NSR))
    open(201,FILE='sr.dat')
    do i=1,NSR
      do j=1,NT
        read(201,*)sr(j,i)
      enddo
      read(201,*)
      read(201,*)
    enddo
    close(201)
 
! Stations along the rupture direction and perpendicularly to the fault    
!    staL(1)=L*2.; staY(1)=0.
!    staL(2)=-L;  staY(2)=0.
!    staL(3)=L/2.;staY(3)=L
    
    open(201,FILE='stationsRIKseissimul.dat')
    do i=1,NR
      read(201,*)staL(i),staY(i)
    enddo
    
    ALLOCATE(seis(NT2,NR))
    seis=0.
!$OMP parallel do private(j,k,tshift) DEFAULT(SHARED)
    do j=1,NR
      do k=1,NSR
        tshift=int(sqrt((SRl(k)-hl-staL(j))**2+((hw-SRw(k))*cos(dip/180.*PI)-staY(j))**2+(hypodepth+(hw-SRw(k))*sin(dip/180.*PI))**2)/VSt/dt)
        seis(tshift+1:tshift+NT-1,j)=seis(tshift+1:tshift+NT-1,j)+sr(1:NT,k)*SRmu(k)*SRelem(k)
      enddo
    enddo
!$OMP end parallel do

    open(201,FILE='seismograms.dat')
    do i=1,NT2
      write(201,'(1000E13.5)')dt*(i-1),seis(i,:)
    enddo
    close(201)
    
    write(*,*)sum(seis(:,1))*dt  !Should give seismic moment (true?)
    
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
