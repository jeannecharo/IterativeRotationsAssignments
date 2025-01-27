gfortran -O3 -march=native -fdefault-real-8 -fbounds-check -Wall -JObj -IObj -o compare_catalog.x compare_catalog.f90 libira.a -llapack
