all: main
main: main.f90
	gfortran main.f90 -o main
depend: 
	perl ./f90makedepend $(SOURCES) > .depend
clean:
	rm $(TARGET) *~ *.o *.mod

include .depend
