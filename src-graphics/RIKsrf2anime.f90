    IMPLICIT NONE

    REAL*8, PARAMETER:: sizex=0.48d0,sizey=0.24d0
    REAL*8,ALLOCATABLE,DIMENSION(:):: L,W
    INTEGER,ALLOCATABLE,DIMENSION(:):: NL,NW
    INTEGER NT,NS
    REAL*8 dum,maxik,minik
    REAL*8,ALLOCATABLE:: srf(:,:,:)
    REAL*8 DL,DW,DT
    INTEGER i,j,k,kk,m,tfrom,tto,NSeg

    NSeg=1
    allocate(NL(NSeg),NW(NSeg),L(NSeg),W(NSeg))
    open(10,FILE='RIKsrf.in')
    read(10,*)
    read(10,*)(L(kk),W(kk),kk=1,NSeg)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)kk
    if(kk.ne.1)then
      write(*,*)'Only regular grid plotting is implemented'
      stop
    endif
    read(10,*)(NL(kk),NW(kk),kk=1,NSeg)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)kk
    if(kk.ne.1)read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)DT,NT
    close(10)

    NS=int(DT*NT)  !Number of slides equals the slip rate time window
    open(101,FILE='sr.dat');write(*,*)'Converting sr.dat'

    open(202,FILE='RIKsrf2anime.dat')
    open(201,FILE='RIKsrf2anime.gp')
    do kk=1,NSeg
      allocate(srf(NT,NL(kk),NW(kk)))
      DL=L(kk)/dble(NL(kk))
      DW=W(kk)/dble(NW(kk))

      do k=1,NW(kk)
        do j=1,NL(kk)
          do i=1,NT
            read(101,*)srf(i,j,k)
          enddo
        enddo
      enddo

      maxik=0.d0;minik=0.d0
      do i=1,NS
        do j=1,NW(kk)
          tfrom=int(dble((i-1)*NT)/dble(NS)+1)
          tto=tfrom+1
          write(202,'(1000E13.5)')(sum(srf(tfrom:tto,k,j))/dble(tto-tfrom+1),k=1,NL(kk)),sum(srf(tfrom:tto,NL(kk),j))/dble(tto-tfrom+1)
        enddo
        write(202,'(1000E13.5)')(sum(srf(tfrom:tto,k,NW(kk)))/dble(tto-tfrom+1),k=1,NL(kk)),sum(srf(tfrom:tto,NL(kk),NW(kk)))/dble(tto-tfrom+1)
        write(202,*);write(202,*)
        dum=maxval(sum(srf(tfrom:tto,1:NL(kk),1:NW(kk)),dim=1))/dble(tto-tfrom+1)
        if(dum>maxik)maxik=dum
        dum=minval(sum(srf(tfrom:tto,1:NL(kk),1:NW(kk)),dim=1))/dble(tto-tfrom+1)
        if(dum<minik)minik=dum
      enddo

      write(201,*)'set term postscript portrait color enh'
      write(201,'(A24,I1,A3)')'set output "RIKsrf2anime',kk,'.ps"'
      write(201,*)'set multiplot'
      write(201,*)'set size ',sizex,',',sizey
      write(201,*)'set size ratio -1'
      write(201,*)'set pm3d map corners2color c1'
      write(201,'(80A)')'set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )'
      write(201,*)'DW=',DW
      write(201,*)'DL=',DL
      write(201,*)' W=',W(kk)
      write(201,*)'set cbrange [',minik,':',maxik,']'
      write(201,*)'unset colorbox'
      write(201,*)'set xrange [0:',L(kk),']'
      write(201,*)'set yrange [0:',W(kk),']'
      write(201,*)'set xtics out scale 0.5 offset 0,.6'
      write(201,*)'set ytics out scale 0.5 offset 0.6,0'
      write(201,*)'set label 1 "" front at graph .98,graph .88 right'

      do i=1,NS
        write(201,*)'set format x ""'
        write(201,*)'set format y ""'
        write(201,*)'unset xlabel'
        write(201,*)'unset ylabel'
        if(i==NS/2)then
          write(201,*)'set format x "%g"'
          write(201,*)'set format y "%g"'
          write(201,*)'set xlabel "Along strike (km)" offset 0,1'
          write(201,*)'set ylabel "Up-dip (km)"'
        endif
        if(i==NS)then
          write(201,*)'set colorbox horizontal user origin 0.41,',real(1.-dble(mod(i-1,NS/2)+1)*sizey/1.7+.04),' size 0.3,0.01'
          write(201,*)'set cblabel "Slip rate (m/s)" offset 0,.5'
!          write(201,*)'set cbtics 1'
        endif
        write(201,'(A13,I2,A1,I2,A2)')'set label 1 "',i-1,'s"'
        write(201,*)'set origin ',sizex*((i-1)/(NS/2))/1.5,',',1.-dble(mod(i-1,NS/2)+1)*sizey/1.7
        write(201,'(A58,I5,A20)')'splot "RIKsrf2anime.dat" matrix u ($1*DL):($2*DW):3 index ',(kk-1)*NS+i-1,' notitle w pm3d,\'
        write(201,*)'"nucleationpoint.dat" u 1:2:(0) index ',kk-1,' notitle w p pt 3 lc 3 ps 1.'
      enddo
      write(201,*)'unset multiplot'
      deallocate(srf)
    enddo
    
    close(101)


    END
