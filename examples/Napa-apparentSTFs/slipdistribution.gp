set term postscript
set output 'slipdistribution.ps'
set size ratio -1
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'
set cblabel 'Slip (m)'
set xrange [0:32]
set yrange [0:20]
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )   #nazvy barev zobrazi show colornames
plot 'slipdistribution.dat' u 1:2:3 notitle w p lc palette pt 6,'strongmotionarea.dat' u 1:2 notitle w l lt 9 lw 5,'nucleationpoint.dat' notitle w p pt 3 ps 4 lt 1
