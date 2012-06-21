SOURCES = vector.f90 proceed.f90 settings.f90 force.f90 main.f90 property.f90 physconst.f90 berendsen.f90 mpi_replica.f90
OBJECTS = ${patsubst %.f90, %.o, $(SOURCES)}
TARGET  = main
FC = mpif90 -m64
all: $(TARGET)
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@
main: $(OBJECTS)
	$(FC) $(OBJECTS) -o $@
depend: 
	perl ./f90makedepend $(SOURCES) > .depend
clean:
	rm $(TARGET) *~ *.o *.mod

include .depend
