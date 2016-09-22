set term postscript color solid 20
set output 'sr.ps'
set xlabel 'Time [s]'
set ylabel 'Slip rate [m/s]'
set xrange [0:12]

plot 'sr.dat' index 38200 t '38200' w l lw 3,\
'' index 45200 t '45200' w l lw 3,\
'' index 47300 t '47300' w l lw 3,\
'' index 46300 t '46300' w l lw 3,\
'' index 43100 t '43100' w l lw 3,\
'' index 43550 t '43550' w l lt 7 lw 3
