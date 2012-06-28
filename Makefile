all: main
main: main.f90
	gfortran main.f90 -o main
depend: 
	perl ./f90makedepend main.f90 > .depend
clean:
	rm main *~ *.o *.mod

include .depend
