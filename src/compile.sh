ifort -oRIKsrfRandomruptvel -openmp -O3 RIKsrfRandomruptvel.f90

gcc -C -Wall -DNO_IEEE_INFINITY -O3 -c Time_2d.c
ifort -C -oRIKrandomruptvel RIKrandomruptvel.f90 Time_2d.o

ifort -oRIKseisdwn -openmp RIKseisdwn.f90 nr.for filters.for
