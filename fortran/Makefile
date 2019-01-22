FC=gfortran
#FCFLAGS=-I/usr/include/ -lFPL
FCFLAGS=-I/usr/include/ -lnetcdf -lnetcdff

build:
	  $(FC) -c tools.F90 gen_potveg_CESM.F90 $(FCFLAGS)
		$(FC) -o gen_potveg_CESM.exe tools.o gen_potveg_CESM.o $(FCFLAGS)
		#$(FC) -c -o gen_potveg_CESM.exe gen_potveg_CESM.F90 $(FCFLAGS)
