ifort -openmp -O3 -oRIKsrfRandomruptvel RIKsrfRandomruptvel.f90
gcc -C -Wall -DNO_IEEE_INFINITY -O3 -c Time_2d.c
ifort -C -oRIKrandomruptvel RIKrandomruptvel.f90 Time_2d.o

./RIKstf

gnuplot slipdistribution.gp
