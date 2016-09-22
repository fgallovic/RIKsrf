ifort -openmp -O3 -oRIKsrfRandomruptvel RIKsrfRandomruptvel.f90
gcc -C -Wall -DNO_IEEE_INFINITY -O3 -c Time_2d.c
ifort -C -oRIKrandomruptvel RIKrandomruptvel.f90 Time_2d.o

ifort -openmp -oseissimul seissimul.f90

./RIKstf
./seissimul
python stfspec.py
python seismogramsspec.py
gnuplot specplot.gp
gnuplot timespecplot.gp
