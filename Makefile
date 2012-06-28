all: main
main: main.f90
	gfortran main.f90 -o main
clean:
	rm main *~ *.o *.mod
