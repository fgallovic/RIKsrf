ifort -O3 -autodouble -ocnv_nez cnv_nez.for
ifort -O3 -autodouble -ogr_nez gr_nez.for
ifort -O3 -oprepare prepare.f90
ifort -O3 -oresort resort.f90

./prepare
rm -fr dat
mkdir dat
