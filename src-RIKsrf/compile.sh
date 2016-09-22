gcc  -c -Wall -DNO_IEEE_INFINITY Time_2d.c
ifort -openmp -oRIKsrf2 RIKsrf2.f90 Time_2d.o
