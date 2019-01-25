FC=gfortran
#FCFLAGS=-I/usr/include/ -lFPL
FCFLAGS=-I/usr/include/ -lnetcdf -lnetcdff -ffree-line-length-none -ffixed-line-length-none

all: tools gen_potveg remap_rochedo

tools:
	$(FC) -c tools.F90 $(FCFLAGS)

gen_potveg: tools
	$(FC) -c gen_potveg_CESM.F90 $(FCFLAGS)
	$(FC) -o gen_potveg_CESM.exe tools.o gen_potveg_CESM.o $(FCFLAGS)

remap_rochedo: tools
	$(FC) -c remap_rochedo.F90 $(FCFLAGS)
	$(FC) -o remap_rochedo.exe tools.o remap_rochedo.o $(FCFLAGS)
	
clean:
	rm -f *.o *.mod *.exe

.PHONY: all clean
