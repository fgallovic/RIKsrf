#./RIKsrf2

ifort -openmp -oseissimul seissimul.f90
./seissimul
python stfspec.py
python seismogramsspec.py
gnuplot specplot.gp
gnuplot timespecplot.gp
