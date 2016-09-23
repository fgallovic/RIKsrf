Vs=3.5     #S-wave velocity in km/s
M0=1.6e18  #M0 in Mm
Brune(x,M0,sd)=M0/(1.+(x/(4.9e6*Vs*(sd/M0/1.e7)**(1./3.)))**2)  #stress drop in bars

set term postscript portrait enhanced color solid 18
set output 'stfspec.ps'
set multiplot
set size 1,.45
set format y '%4.1f'
#set colors classic    #To go back to classic colors in new Gnuplot versions

set origin 0,.5
set key at screen 0.5,0.05
set xlabel 'Time (s)'
set ylabel 'Moment rate (*10^{18} Nms^{-1})'
plot [0.:8.] 'stf.dat' u ($1+0.0):($2/1.e18) notitle w l lt 1 lw 2

set origin 0,.05
set key top right vertical inside
set logscale xy
set xlabel 'Frequency (Hz)'
set ylabel 'Source spectrum (*10^{18} Nm)'
plot [0.05:5] Brune(x,M0,10)/1.e18 title 'Brune spectrum, {/Symbol Ds}=1MPa' w l lt 9 lw 2,\
Brune(x,M0,50)/1.e18 title 'Brune spectrum, {/Symbol Ds}=5MPa' w l lt 9 lw 4,\
Brune(x,M0,100)/1.e18 title 'Brune spectrum, {/Symbol Ds}=10MPa' w l lt 9 lw 2,\
'stfspecD.dat' u 1:($2/1.e18) notitle w l lt 1 lw 2
