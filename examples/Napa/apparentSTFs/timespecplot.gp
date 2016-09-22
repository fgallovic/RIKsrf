set term postscript landscape color solid 12
set output 'timespecplot.ps'
set multiplot
set size .5,.5

set xlabel 'Time (s)'
set ylabel 'Moment rate (Nm/s)'

set origin 0.,0.5
set xrange [0:10]
#set yrange [3e15:4e17]
set title 'Source time function'
plot 'stf.dat' u 1:2 notitle w l lt 9 lw 2

set origin 0.,0.0
set xrange [0:40]
#set yrange [1e15:1e18]
set title 'Apparent time functions'
plot 'seismograms.dat' u 1:2 t 'Right'w l lt 1 lw 2,'' u 1:3 t 'Left' w l lt 2 lw 2,'' u 1:4 t 'Perpendicular' w l lt 3 lw 2


set logscale xy
set xlabel 'Frequency (Hz)'
set ylabel 'Amplitude spectrum (m/s)'
set xrange [0.05:10]

set origin 0.5,0.5
set yrange [3e16:4e18]
set title 'Source acc. spectrum'
plot 'stfspecA.dat' u 1:2 notitle w l lt 9 lw 2

set origin 0.5,0.0
set yrange [1e16:1e19]
set title 'Apparent source acc. spectra'
plot 'seismogramsspecA.dat' u 1:2 t 'Right'w l lt 1 lw 2,'' u 1:3 t 'Left' w l lt 2 lw 2,'' u 1:4 t 'Perpendicular' w l lt 3 lw 2

unset multiplot
